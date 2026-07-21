# diag_sono_smoothing.R
# Diagnostic: compare different sono smoothing/interpolation approaches so we
# can choose the most sensible one before point-by-point comparison with
# encoder-predicted strain.
#
# Background: the DS3 sonomicrometer updates its analog output at its internal
# processing rate (~247 Hz for this session, stored in
# daq_sono_internal_sample_rate_hertz). The NI AI samples it at 1000 Hz,
# producing a staircase waveform (~4 AI samples per DS3 update). Any
# point-by-point comparison of sono strain with continuously-computed encoder
# strain requires smoothing/interpolating sono to remove the staircase first.
#
# Four representative conditions from trial 18 (passive, passive, all freqs;
# same session/specimen so geometry is constant across panels):
#   Panel A: small amplitude (~4.7 deg) + slow (1 Hz)
#   Panel B: small amplitude (~4.7 deg) + fast (7 Hz)
#   Panel C: large amplitude (~18.8 deg) + slow (1 Hz)
#   Panel D: large amplitude (~18.8 deg) + fast (7 Hz)
#
# Smoothing approaches compared (each overlaid in each panel):
#   1. Raw sono             -- shows the staircase plainly
#   2. Rolling mean 4 samp  -- one DS3 update period (1000/247 ~ 4 samples)
#   3. Butterworth LP 40 Hz -- below DS3 Nyquist; preserves muscle motion
#   4. Butterworth LP 120 Hz-- half DS3 internal rate (theoretical Nyquist)
#   5. Lin interp           -- detect DS3 step transitions, interpolate between
#
# Encoder-predicted strain is overlaid in all panels as the reference.
#
# Run with:  Rscript R/diag_sono_smoothing.R
# Outputs -> figs_diagnostic/sono_smoothing_comparison.png

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(cli); library(rhdf5); library(signal)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")

# Defaults come from paths_config.R (single source of truth) -- see that
# file if the OneDrive folder layout ever moves again.
SOURCE_DIR  <- raw_source_dir(BASS16_RAW_SUBFOLDER)
OUT_DIR     <- FIGS_DIAGNOSTIC_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

src("00_load_bender_flat.R")
src("01_calibrate.R")
src("muscle_geometry.R")
src("plot_strain_validation.R")      # attach_measured_strain
src("plot_angle_sono_validation.R")  # attach_sono_strain / .read_sono_right_mm_aligned

SAMPLE_RATE_HZ  <- 1000.0
DS3_RATE_HZ     <- 247.0   # daq_sono_internal_sample_rate_hertz from file attrs
BODY_WIDTH_MM   <- 42.0
MUSCLE_DEPTH_MM <- 0.0     # unset default -> resolve_muscle_depth_mm falls back to 1 mm
DCLAMP_MM       <- 41.0
LIDX_POS_MOTOR  <- -1.0
LIDX_RIGHT      <-  1.0

TRIAL_18 <- file.path(SOURCE_DIR, "2026-07-14_bass16_bender_18_dynamic.h5")

# =============================================================================
# Smoothing functions
# =============================================================================

#' Rolling mean over n samples (causal: averages current + n-1 prior samples).
.rollmean <- function(x, n) {
  if (n <= 1L) return(x)
  stats::filter(x, rep(1/n, n), sides = 1L) |> as.numeric()
}

#' Zero-phase Butterworth LP filter via filtfilt (requires signal package).
.butter_lp <- function(x, cutoff_hz, sample_rate_hz, order = 4L) {
  nyq <- sample_rate_hz / 2.0
  if (cutoff_hz >= nyq) return(x)
  ok <- is.finite(x)
  if (sum(ok) < 20L) return(x)
  filt <- signal::butter(order, cutoff_hz / nyq, type = "low")
  out  <- x
  out[ok] <- signal::filtfilt(filt, x[ok])
  out
}

#' Linear interpolation at the DS3's native update rate.
#' Strategy: subsample the raw staircase to one value per DS3-update period
#' (every round(sample_rate/ds3_rate) samples), then linearly interpolate
#' those anchor points back to the original 1000 Hz grid. This avoids the
#' staircase without requiring threshold-based step detection, which is
#' unreliable when ADC noise within a hold period has similar amplitude to
#' small DS3 steps.
.spline_interp <- function(x, t, ds3_rate_hz = DS3_RATE_HZ, sample_rate_hz = SAMPLE_RATE_HZ) {
  ok <- is.finite(x)
  if (sum(ok) < 10L) return(x)
  x_ok <- x[ok]; t_ok <- t[ok]
  step_n <- max(1L, as.integer(round(sample_rate_hz / ds3_rate_hz)))  # ~4 samples
  # pick one sample per DS3 update period as anchor
  anchor_idx <- seq(1L, length(x_ok), by = step_n)
  anchor_t   <- t_ok[anchor_idx]
  anchor_v   <- x_ok[anchor_idx]
  keep <- !duplicated(anchor_t)
  anchor_t <- anchor_t[keep]; anchor_v <- anchor_v[keep]
  if (length(anchor_t) < 2L) return(x)
  out <- x
  out[ok] <- approx(anchor_t, anchor_v, xout = t_ok, rule = 2L)$y
  out
}

# =============================================================================
# Load trial 18 and attach all strain signals
# =============================================================================
cli::cli_h1("Loading trial 18 (freq sweep / passive: 1-7 Hz, 4 amplitudes)")
td18 <- load_bender_flat(TRIAL_18, do_filter = FALSE, loadtorques = "x")
info18_h5 <- tryCatch({
  h5 <- H5Fopen(TRIAL_18, "H5F_ACC_RDONLY")
  on.exit(try(H5Fclose(h5), silent = TRUE), add = TRUE)
  ci   <- as.integer(h5read(h5, "/timeseries/cycle_index"))
  freq <- as.numeric(h5read(h5, "/metadata/index_cycle_frequency_hertz"))
  amp  <- as.numeric(h5read(h5, "/metadata/index_cycle_motor_amplitude_degree"))
  list(ci = ci, freq = freq, amp = amp)
}, error = function(e) NULL)

if (!is.null(info18_h5)) {
  td18$cycle_idx <- info18_h5$ci
  N <- nrow(td18); C <- length(info18_h5$freq)
  in_range <- !is.na(td18$cycle_idx) & td18$cycle_idx >= 1L & td18$cycle_idx <= C
  td18$freq_cycle <- NA_real_; td18$amp_cycle <- NA_real_
  if (any(in_range)) {
    idx <- td18$cycle_idx[in_range]
    td18$freq_cycle[in_range] <- info18_h5$freq[idx]
    td18$amp_cycle[in_range]  <- info18_h5$amp[idx]
  }
}

# Attach predicted strain (encoder-based) via existing helpers
td18 <- attach_predicted_strain(
  td18, local_body_width_mm = BODY_WIDTH_MM, measured_muscle_depth_mm = MUSCLE_DEPTH_MM,
  active_mask = rep(FALSE, nrow(td18))
)
td18 <- attach_measured_strain(td18)     # strain_measured_pct from enc.deg
td18 <- attach_sono_strain(td18, TRIAL_18, LIDX_RIGHT, LIDX_POS_MOTOR)

# strain_pred_encoder_right_pct = FORCE_SIGN_RIGHT * strain_measured_pct
# (already computed by attach_sono_strain -- used as the reference "predicted" column)

# =============================================================================
# Compute smoothed sono variants
# =============================================================================
# Work on the raw sono in mm (sono_right_mm), then convert to strain at end.
# We smooth BEFORE converting to strain so the V->mm calibration nonlinearity
# doesn't interact with the smoothing kernel.

cli::cli_h1("Computing smoothed sono variants")

# Read raw sono mm directly (attach_sono_strain already stored it in td18)
sono_raw_mm <- td18$sono_right_mm
t_vec       <- td18$t.s

# Reference length (L0) from near-zero commanded angle (same as attach_sono_strain uses)
L0 <- .sono_reference_length_mm(td18$angle.deg, sono_raw_mm)
cli::cli_alert_info("L0 (sono reference length) = {round(L0, 2)} mm")

# Strain conversion (applied to each smoothed version)
.mm_to_strain_pct <- function(mm_vec, L0) {
  -(mm_vec - L0) / L0 * 100.0
}

n_rollmean <- as.integer(round(SAMPLE_RATE_HZ / DS3_RATE_HZ))  # ~ 4 samples
cli::cli_alert_info("Rolling mean window = {n_rollmean} samples (~1 DS3 update period)")

td18$sono_raw_pct       <- .mm_to_strain_pct(sono_raw_mm, L0)
td18$sono_roll4_pct     <- .mm_to_strain_pct(.rollmean(sono_raw_mm, n_rollmean), L0)
td18$sono_lp40_pct      <- .mm_to_strain_pct(.butter_lp(sono_raw_mm, 40.0,  SAMPLE_RATE_HZ), L0)
td18$sono_lp120_pct     <- .mm_to_strain_pct(.butter_lp(sono_raw_mm, 120.0, SAMPLE_RATE_HZ), L0)
td18$sono_spline_pct    <- .mm_to_strain_pct(.spline_interp(sono_raw_mm, t_vec), L0)

cli::cli_alert_success("Smoothing complete")

# =============================================================================
# Panel selection: pick 2 cycles from the middle of each freq/amp block
# =============================================================================

# Target conditions: small/large amplitude x slow/fast frequency
# From survey: amp range ~ 4.7 (small) / 18.8 (large); freq: 1 Hz (slow) / 7 Hz (fast)
PANELS <- list(
  list(label = "A: small amp, slow\n(1 Hz, ~4.7 deg)", freq = 1.0,  amp_approx = 4.7),
  list(label = "B: small amp, fast\n(7 Hz, ~4.7 deg)", freq = 7.0,  amp_approx = 4.7),
  list(label = "C: large amp, slow\n(1 Hz, ~18.8 deg)",freq = 1.0,  amp_approx = 18.8),
  list(label = "D: large amp, fast\n(7 Hz, ~18.8 deg)",freq = 7.0,  amp_approx = 18.8)
)

SIGNAL_LEVELS <- c(
  "Predicted (encoder)",
  "Raw sono",
  "Roll-mean 4 samp",
  "LP 40 Hz",
  "LP 120 Hz",
  "Lin interp (DS3 steps)"
)
SIGNAL_COLORS <- c(
  "Predicted (encoder)"    = "black",
  "Raw sono"               = "#cccccc",
  "Roll-mean 4 samp"       = "#f97316",
  "LP 40 Hz"               = "#16a34a",
  "LP 120 Hz"              = "#2563eb",
  "Lin interp (DS3 steps)" = "#9333ea"
)
SIGNAL_LTYPE <- c(
  "Predicted (encoder)"    = "solid",
  "Raw sono"               = "solid",
  "Roll-mean 4 samp"       = "solid",
  "LP 40 Hz"               = "solid",
  "LP 120 Hz"              = "solid",
  "Lin interp (DS3 steps)" = "solid"
)

# Build tidy long tibble for ggplot, one block per panel
panel_dfs <- lapply(PANELS, function(p) {
  freq_target <- p$freq
  amp_target  <- p$amp_approx

  # find matching cycles (allow some tolerance on amp)
  mask <- !is.na(td18$freq_cycle) & !is.na(td18$amp_cycle) &
          abs(td18$freq_cycle - freq_target) < 0.01 &
          abs(td18$amp_cycle  - amp_target)  < 1.0   # ±1 deg tolerance

  sub <- td18[mask, , drop = FALSE]
  if (nrow(sub) < 20L) {
    cli::cli_alert_warning("Panel '{p$label}': too few samples, skip")
    return(NULL)
  }

  # 2 cycles from middle of block
  period_s <- 1.0 / freq_target
  t_mid    <- median(sub$t.s, na.rm = TRUE)
  t_lo     <- t_mid - period_s
  t_hi     <- t_mid + period_s
  slice    <- sub[sub$t.s >= t_lo & sub$t.s <= t_hi, , drop = FALSE]

  if (nrow(slice) < 20L) {
    # if too few, take the full block
    slice <- sub
    t_lo  <- min(sub$t.s, na.rm = TRUE)
  }

  # Normalize time to start at 0 within the slice for cleaner axis
  t_norm <- slice$t.s - min(slice$t.s, na.rm = TRUE)

  tibble(
    panel    = p$label,
    t        = t_norm,
    `Predicted (encoder)`    = slice$strain_pred_encoder_right_pct,
    `Raw sono`               = slice$sono_raw_pct,
    `Roll-mean 4 samp`       = slice$sono_roll4_pct,
    `LP 40 Hz`               = slice$sono_lp40_pct,
    `LP 120 Hz`              = slice$sono_lp120_pct,
    `Lin interp (DS3 steps)` = slice$sono_spline_pct
  )
})

panel_df <- dplyr::bind_rows(Filter(Negate(is.null), panel_dfs))

if (nrow(panel_df) == 0L) {
  cli::cli_abort("No panel data assembled -- check freq/amp matching")
}

# Pivot to long format for ggplot
panel_long <- panel_df |>
  tidyr::pivot_longer(
    cols      = all_of(SIGNAL_LEVELS),
    names_to  = "signal",
    values_to = "strain_pct"
  ) |>
  dplyr::filter(is.finite(.data$strain_pct)) |>
  dplyr::mutate(
    signal = factor(.data$signal, levels = SIGNAL_LEVELS),
    panel  = factor(.data$panel,  levels = vapply(PANELS, `[[`, "", "label"))
  )

# =============================================================================
# Plot
# =============================================================================
cli::cli_h1("Building plot")

# Line widths: predicted heavier, raw thinnest, smoothed versions medium
lwd_map <- c(
  "Predicted (encoder)" = 1.1,
  "Raw sono"            = 0.5,
  "Roll-mean 4 samp"    = 0.85,
  "LP 40 Hz"            = 0.85,
  "LP 120 Hz"           = 0.85,
  "Spline interp"       = 0.85
)
# Alpha: raw sono transparent so smoothed versions are visible on top
alpha_map <- c(
  "Predicted (encoder)" = 1.0,
  "Raw sono"            = 0.45,
  "Roll-mean 4 samp"    = 1.0,
  "LP 40 Hz"            = 1.0,
  "LP 120 Hz"           = 1.0,
  "Spline interp"       = 1.0
)

# We'll plot layer by layer so raw sono goes UNDER everything else.
df_raw    <- dplyr::filter(panel_long, .data$signal == "Raw sono")
df_pred   <- dplyr::filter(panel_long, .data$signal == "Predicted (encoder)")
df_smooth <- dplyr::filter(panel_long, !(.data$signal %in% c("Raw sono","Predicted (encoder)")))

p <- ggplot(panel_long, aes(x = .data$t, y = .data$strain_pct)) +
  # raw sono first (bottom layer, greyed out)
  geom_line(data = df_raw,
            aes(color = .data$signal),
            linewidth = lwd_map["Raw sono"], alpha = alpha_map["Raw sono"]) +
  # smoothed variants
  geom_line(data = df_smooth,
            aes(color = .data$signal, group = .data$signal),
            linewidth = lwd_map["LP 40 Hz"], alpha = 1.0) +
  # encoder predicted on top
  geom_line(data = df_pred,
            aes(color = .data$signal),
            linewidth = lwd_map["Predicted (encoder)"], alpha = 1.0, linetype = "dashed") +
  scale_color_manual(values = SIGNAL_COLORS, name = NULL,
                     breaks = SIGNAL_LEVELS) +
  facet_wrap(~panel, scales = "free", ncol = 2L) +
  labs(
    title    = "Sono smoothing comparison vs encoder-predicted strain (trial 18, passive)",
    subtitle = paste0(
      "Predicted (encoder) = curvature * r_m, right-folded (dashed black)  |  ",
      "Raw sono = DS3 staircase at 1000 Hz AI  |  DS3 internal rate ~247 Hz\n",
      "Roll-mean 4 samp = 1 DS3-update-period average  |  LP = zero-phase Butterworth  |  ",
      "Lin interp = linear between detected DS3 step transitions"
    ),
    x = "Time within shown window (s)",
    y = "Right-muscle strain (%, shortening positive)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.key.width = unit(1.6, "cm"),
    strip.text       = element_text(face = "bold", size = 10)
  ) +
  guides(color = guide_legend(nrow = 1L, override.aes = list(linewidth = 1.2)))

fout <- file.path(OUT_DIR, "sono_smoothing_comparison.png")
ggplot2::ggsave(fout, p, width = 13, height = 8, dpi = 150)
cli::cli_alert_success("Saved {fout}")

# =============================================================================
# Supplemental: correlation table for each panel x smoothing approach
# =============================================================================
cli::cli_h1("Correlation vs encoder-predicted, per panel x smoothing approach")

smooth_cols <- c(
  "Raw sono"               = "sono_raw_pct",
  "Roll-mean 4 samp"       = "sono_roll4_pct",
  "LP 40 Hz"               = "sono_lp40_pct",
  "LP 120 Hz"              = "sono_lp120_pct",
  "Lin interp (DS3 steps)" = "sono_spline_pct"
)

cor_rows <- list()
for (p_def in PANELS) {
  freq_target <- p_def$freq
  amp_target  <- p_def$amp_approx
  mask <- !is.na(td18$freq_cycle) & !is.na(td18$amp_cycle) &
          abs(td18$freq_cycle - freq_target) < 0.01 &
          abs(td18$amp_cycle  - amp_target)  < 1.0
  sub <- td18[mask, , drop = FALSE]
  pred <- sub$strain_pred_encoder_right_pct
  for (label in names(smooth_cols)) {
    meas <- sub[[smooth_cols[label]]]
    ok   <- is.finite(pred) & is.finite(meas)
    r    <- if (sum(ok) > 5) suppressWarnings(cor(pred[ok], meas[ok])) else NA_real_
    rmse <- if (sum(ok) > 5) sqrt(mean((pred[ok] - meas[ok])^2)) else NA_real_
    cor_rows[[length(cor_rows)+1]] <- tibble(
      panel  = p_def$label,
      signal = label,
      r      = round(r, 4),
      rmse   = round(rmse, 4),
      n      = sum(ok)
    )
  }
}

cor_tbl <- dplyr::bind_rows(cor_rows)
cat("\n=== r and RMSE(%) vs encoder-predicted, per panel x smoothing ===\n")
print(cor_tbl, n = Inf)

cli::cli_alert_success("diag_sono_smoothing.R complete -- output: {fout}")
