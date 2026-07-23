# diag_sono_amp_ratio_destaircased.R
# Diagnostic (not a pipeline deliverable): re-runs diag_sono_amp_ratio_
# multitrial.R's per-cycle amplitude-ratio check with the sono channel
# DE-STAIRCASED first, to test whether the amplitude-dependent excess found
# there (measured p2p >> predicted p2p at small cycle amplitudes, converging
# toward a ~1.15-1.3x floor at large amplitudes -- see that script and
# analysis_muscle_force_vector_log.md) is a real biological/mechanical
# effect or a signal-conditioning artifact of the sono channel itself.
#
# The concern (PI-raised, verified against file metadata): the DS3
# sonomicrometer updates internally at daq_sono_internal_sample_rate_hertz
# (~247 Hz, confirmed identical across bass16/17/18), but the NI AI clock
# samples ai6 at daq_ai_sample_rate_hertz (1000 Hz) -- i.e. every DS3 value
# is HELD across ~4 consecutive AI samples, producing a staircase. A raw
# per-cycle max()-min() over that staircase is an EXTREME-VALUE statistic
# on top of whatever ADC/electrical noise rides on each held plateau; if
# that noise has a roughly constant absolute magnitude regardless of true
# cycle amplitude, it would inflate the small-amplitude cycles' p2p far
# more (proportionally) than the large-amplitude ones -- which is exactly
# the amplitude-dependent pattern already observed. diag_sono_amp_ratio_
# multitrial.R's p2p was computed from the RAW staircase (attach_sono_strain
# / .read_sono_right_mm_aligned() applies calibration only, no smoothing),
# so that pattern is fully consistent with a signal-conditioning artifact
# and has NOT yet been ruled out.
#
# Three sono conditions per trial (smoothing functions ported from
# diag_sono_smoothing.R, applied to sono_right_mm BEFORE strain conversion
# so the V->mm calibration nonlinearity doesn't interact with the kernel):
#   1. Raw          -- staircase, as used by diag_sono_amp_ratio_multitrial.R
#   2. Roll-mean 4  -- averages across one DS3 update period (~4 AI samples)
#   3. Lin-interp   -- subsample to 1 pt/DS3-period, linearly interpolate
#                       back to 1000 Hz (removes the staircase edges without
#                       a frequency-domain filter's transient/phase effects)
#
# If de-staircasing collapses (or greatly shrinks) the amplitude-dependent
# floor, the excess was substantially a sono signal-conditioning artifact.
# If the floor persists near-unchanged, the raw-staircase hypothesis is
# NOT the explanation and the amplitude-dependent excess needs a different
# cause.
#
# Run with:  Rscript R/diag_sono_amp_ratio_destaircased.R
# Outputs -> figures: 02_processed/figs_diagnostic/ (FIGS_DIAGNOSTIC_DIR)
#            data:    02_processed/data_processed/

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2); library(cli); library(rhdf5)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")
src("00_load_bender_flat.R")
src("01_calibrate.R")
src("muscle_geometry.R")
src("plot_strain_validation.R")     # attach_measured_strain
src("plot_angle_sono_validation.R") # .read_sono_right_mm_aligned, .sono_reference_length_mm
src("plot_force_vs_time.R")         # .detect_stim_events (used by .stim_window_state_label)

OUT_DIR      <- FIGS_DIAGNOSTIC_DIR
DATA_OUT_DIR <- file.path(.crittergripper_root(), "02_processed", "data_processed")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SAMPLE_RATE_HZ <- 1000.0
MIN_P2P_PRED_PCT <- 0.5
AMP_BIN_BREAKS <- c(0, 2, 5, 10, 20, 40)
AMP_BIN_LABELS <- c("0-2%", "2-5%", "5-10%", "10-20%", "20-40%")

SPECIMEN_SUBFOLDERS <- c(
  bass16 = BASS16_RAW_SUBFOLDER,
  bass17 = BASS17_RAW_SUBFOLDER,
  bass18 = BASS18_RAW_SUBFOLDER
)

# =============================================================================
# Smoothing functions (ported from diag_sono_smoothing.R)
# =============================================================================

.rollmean <- function(x, n) {
  if (n <= 1L) return(x)
  out <- stats::filter(x, rep(1 / n, n), sides = 1L) |> as.numeric()
  # fill the leading NA-from-filter samples with the first valid value so
  # per-cycle max/min near a cycle's start sample isn't spuriously NA'd out
  first_ok <- which(is.finite(out))[1L]
  if (is.finite(first_ok) && first_ok > 1L) out[seq_len(first_ok - 1L)] <- out[first_ok]
  out
}

#' Subsample to 1 anchor per DS3 update period, linearly interpolate back to
#' the original (1000 Hz) grid -- avoids the staircase without a frequency-
#' domain filter's ringing/transient behavior at step edges.
.lin_interp_destaircase <- function(x, t, ds3_rate_hz, sample_rate_hz) {
  ok <- is.finite(x) & is.finite(t)
  if (sum(ok) < 10L) return(x)
  x_ok <- x[ok]; t_ok <- t[ok]
  step_n <- max(1L, as.integer(round(sample_rate_hz / ds3_rate_hz)))
  anchor_idx <- seq(1L, length(x_ok), by = step_n)
  anchor_t <- t_ok[anchor_idx]; anchor_v <- x_ok[anchor_idx]
  keep <- !duplicated(anchor_t)
  anchor_t <- anchor_t[keep]; anchor_v <- anchor_v[keep]
  if (length(anchor_t) < 2L) return(x)
  out <- x
  out[ok] <- approx(anchor_t, anchor_v, xout = t_ok, rule = 2L)$y
  out
}

# =============================================================================
# Helpers (geometry / cycle-freq reads, same pattern as the other diag_sono_*)
# =============================================================================

.dynamic_trial_files <- function(subfolder) {
  d <- raw_source_dir(subfolder)
  sort(list.files(d, pattern = "_dynamic\\.h5$", full.names = TRUE))
}

.read_geom_and_lidx <- function(h5path) {
  h5 <- H5Fopen(h5path, "H5F_ACC_RDONLY")
  on.exit(try(H5Fclose(h5), silent = TRUE), add = TRUE)
  m_attrs <- tryCatch(h5readAttributes(h5, "/metadata"), error = function(e) list())
  dbl1 <- function(key, default = NA_real_) {
    v <- suppressWarnings(as.numeric(m_attrs[[key]][1L]))
    if (length(v) == 0L || is.na(v)) default else v
  }
  list(
    width_mm       = dbl1("measurement_specimen_local_body_width_millimeter"),
    depth_mm       = dbl1("measurement_target_muscle_depth_millimeter"),
    lidx_pos_motor = dbl1("daq_specimen_lateral_index_on_positive_motor_side"),
    lidx_right     = dbl1("daq_specimen_side_index_right"),
    ds3_rate_hz    = dbl1("daq_sono_internal_sample_rate_hertz", default = 247.0)
  )
}

.read_cycle_freq <- function(h5path) {
  h5 <- H5Fopen(h5path, "H5F_ACC_RDONLY")
  on.exit(try(H5Fclose(h5), silent = TRUE), add = TRUE)
  ci_raw  <- tryCatch(h5read(h5, "/timeseries/cycle_index"), error = function(e) NULL)
  freq_ds <- tryCatch(h5read(h5, "/metadata/index_cycle_frequency_hertz"), error = function(e) NULL)
  list(ci = as.integer(ci_raw), freq_by_cycle = as.numeric(freq_ds))
}

.p2p_per_cycle <- function(x, ci_vec, cycles) {
  vapply(cycles, function(c) {
    v <- x[ci_vec == c]
    if (sum(is.finite(v)) < 4L) return(NA_real_)
    max(v, na.rm = TRUE) - min(v, na.rm = TRUE)
  }, numeric(1L))
}

.dominant_stim_state <- function(state_chr) {
  state_chr <- state_chr[!is.na(state_chr)]
  if (length(state_chr) == 0L) return(NA_character_)
  tab <- table(state_chr)
  names(tab)[which.max(tab)]
}

#' Process one dynamic trial -> per-cycle tibble with THREE p2p_meas columns
#' (raw / rollmean4 / lininterp), all sharing the same p2p_pred (encoder).
.process_trial <- function(specimen, h5path) {
  trial_label <- tools::file_path_sans_ext(basename(h5path))

  td <- tryCatch(load_bender_flat(h5path, do_filter = FALSE, loadtorques = "x"),
                 error = function(e) NULL)
  if (is.null(td) || nrow(td) == 0L) return(tibble())

  geom <- .read_geom_and_lidx(h5path)
  if (!is.finite(geom$width_mm) || !is.finite(geom$lidx_right) || !is.finite(geom$lidx_pos_motor)) return(tibble())

  info <- .read_cycle_freq(h5path)
  if (is.null(info$ci) || length(info$ci) != nrow(td)) return(tibble())
  td$cycle_idx <- info$ci
  td$freq_cycle <- dplyr::if_else(
    td$cycle_idx >= 1L & td$cycle_idx <= length(info$freq_by_cycle),
    info$freq_by_cycle[pmax(1L, td$cycle_idx)], NA_real_
  )

  td <- attach_predicted_strain(
    td, local_body_width_mm = geom$width_mm, measured_muscle_depth_mm = geom$depth_mm,
    active_mask = rep(FALSE, nrow(td))
  )
  td <- attach_measured_strain(td)

  # attach_sono_strain() normally builds strain_pred_encoder_right_pct as a
  # side effect of attaching sono strain -- since this script builds its own
  # (multi-condition) sono strain columns instead of calling that helper, the
  # encoder-predicted column has to be built explicitly here, identically.
  force_sign_right <- geom$lidx_right * geom$lidx_pos_motor
  td$strain_pred_encoder_right_pct <- force_sign_right * td$strain_measured_pct

  # Raw sono mm + reference length (same L0 for all smoothing variants, as
  # in diag_sono_smoothing.R, so the DC choice doesn't confound the
  # comparison -- L0 cancels out of p2p anyway, but keep parity).
  sono_raw_mm <- tryCatch(.read_sono_right_mm_aligned(h5path), error = function(e) NULL)
  if (is.null(sono_raw_mm) || length(sono_raw_mm) != nrow(td)) return(tibble())
  L0 <- .sono_reference_length_mm(td$angle.deg, sono_raw_mm)
  if (!is.finite(L0) || L0 <= 0) return(tibble())

  n_roll <- max(1L, as.integer(round(SAMPLE_RATE_HZ / geom$ds3_rate_hz)))

  mm_raw       <- sono_raw_mm
  mm_rollmean4 <- .rollmean(sono_raw_mm, n_roll)
  mm_lininterp <- .lin_interp_destaircase(sono_raw_mm, td$t.s, geom$ds3_rate_hz, SAMPLE_RATE_HZ)

  # NOTE: no force_sign_right fold here -- the sono side is a direct
  # physical length-change measurement (see attach_sono_strain(), the
  # vetted convention in plot_angle_sono_validation.R). This function is
  # only ever used for peak-to-peak amplitude downstream, which is
  # sign-invariant, so this was a latent (harmless-here) inconsistency,
  # not a bug that changed any output of this script.
  .mm_to_strain_right <- function(mm) -(mm - L0) / L0 * 100.0

  td$strain_sono_raw_pct       <- .mm_to_strain_right(mm_raw)
  td$strain_sono_rollmean4_pct <- .mm_to_strain_right(mm_rollmean4)
  td$strain_sono_lininterp_pct <- .mm_to_strain_right(mm_lininterp)

  stim_state_fac <- tryCatch(.stim_window_state_label(td), error = function(e) NULL)
  if (is.null(stim_state_fac)) return(tibble())
  td$stim_state <- as.character(stim_state_fac)

  cyc_list <- sort(unique(td$cycle_idx[td$cycle_idx >= 1L]))
  if (length(cyc_list) == 0L) return(tibble())

  p2p_pred      <- .p2p_per_cycle(td$strain_pred_encoder_right_pct, td$cycle_idx, cyc_list)
  p2p_raw       <- .p2p_per_cycle(td$strain_sono_raw_pct,            td$cycle_idx, cyc_list)
  p2p_rollmean4 <- .p2p_per_cycle(td$strain_sono_rollmean4_pct,      td$cycle_idx, cyc_list)
  p2p_lininterp <- .p2p_per_cycle(td$strain_sono_lininterp_pct,      td$cycle_idx, cyc_list)
  freq_by_c <- vapply(cyc_list, function(c) {
    v <- td$freq_cycle[td$cycle_idx == c]; v <- v[is.finite(v)]
    if (length(v) == 0L) NA_real_ else v[1L]
  }, numeric(1L))
  dom_state <- vapply(cyc_list, function(c) .dominant_stim_state(td$stim_state[td$cycle_idx == c]), character(1L))

  cli::cli_alert_success("{trial_label}: {length(cyc_list)} cycles")

  tibble(
    specimen = specimen, trial = trial_label, cycle = cyc_list, freq_hz = freq_by_c,
    stim_state = dom_state, p2p_pred_pct = p2p_pred,
    p2p_meas_raw_pct = p2p_raw, p2p_meas_rollmean4_pct = p2p_rollmean4, p2p_meas_lininterp_pct = p2p_lininterp
  )
}

# =============================================================================
# Batch
# =============================================================================
cli::cli_h1("Loading all dynamic trials (bass16/17/18), 3 sono conditions each")

all_files <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  tibble(specimen = specimen, h5path = .dynamic_trial_files(subfolder))
})
results <- purrr::pmap_dfr(all_files, function(specimen, h5path) {
  tryCatch(.process_trial(specimen, h5path),
           error = function(e) { cli::cli_alert_danger("{basename(h5path)}: {conditionMessage(e)}"); tibble() })
})
if (nrow(results) == 0L) cli::cli_abort("No usable cycles extracted.")

results <- results |>
  dplyr::mutate(
    amp_bin = cut(.data$p2p_pred_pct, breaks = AMP_BIN_BREAKS, labels = AMP_BIN_LABELS),
    ok_amp  = is.finite(.data$p2p_pred_pct) & .data$p2p_pred_pct >= MIN_P2P_PRED_PCT,
    stim_state = factor(.data$stim_state, levels = c("no stim", "left stim", "right stim"))
  )
write.csv(results, file.path(DATA_OUT_DIR, "amp_ratio_destaircased_all_dynamic_trials.csv"), row.names = FALSE)
cli::cli_alert_success("Saved {nrow(results)} cycles -> data_processed/amp_ratio_destaircased_all_dynamic_trials.csv")

# =============================================================================
# Summary: amp_ratio by amplitude bin, for each sono condition
# =============================================================================
long <- results |>
  dplyr::filter(.data$ok_amp) |>
  tidyr::pivot_longer(
    cols = c(p2p_meas_raw_pct, p2p_meas_rollmean4_pct, p2p_meas_lininterp_pct),
    names_to = "condition", values_to = "p2p_meas_pct"
  ) |>
  dplyr::mutate(
    condition = dplyr::recode(.data$condition,
      p2p_meas_raw_pct = "Raw (staircase)",
      p2p_meas_rollmean4_pct = "Roll-mean 4 samp",
      p2p_meas_lininterp_pct = "Lin-interp (de-staircased)"),
    condition = factor(.data$condition, levels = c("Raw (staircase)", "Roll-mean 4 samp", "Lin-interp (de-staircased)")),
    amp_ratio = .data$p2p_meas_pct / .data$p2p_pred_pct
  )

cli::cli_h1("Amplitude ratio by bin x condition")
summary_by_bin <- long |>
  dplyr::group_by(.data$condition, .data$amp_bin) |>
  dplyr::summarise(n = dplyr::n(), amp_ratio_mean = mean(.data$amp_ratio, na.rm = TRUE),
                    amp_ratio_median = median(.data$amp_ratio, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(.data$amp_bin, .data$condition)
print(summary_by_bin, n = 100)
write.csv(summary_by_bin, file.path(DATA_OUT_DIR, "amp_ratio_destaircased_summary_by_bin.csv"), row.names = FALSE)

cli::cli_h1("Linear fit (p2p_meas ~ p2p_pred) per condition, no-stim cycles only")
fit_tbl <- long |>
  dplyr::filter(.data$stim_state == "no stim") |>
  dplyr::group_by(.data$condition) |>
  dplyr::summarise(
    intercept = coef(lm(p2p_meas_pct ~ p2p_pred_pct))[1],
    slope     = coef(lm(p2p_meas_pct ~ p2p_pred_pct))[2],
    r2        = summary(lm(p2p_meas_pct ~ p2p_pred_pct))$r.squared,
    .groups = "drop"
  )
print(fit_tbl)

cli::cli_h1("Amplitude ratio by stim state (matched amplitude band 5-10%), per condition")
matched_tbl <- long |>
  dplyr::filter(.data$p2p_pred_pct >= 5, .data$p2p_pred_pct <= 10) |>
  dplyr::group_by(.data$condition, .data$stim_state) |>
  dplyr::summarise(n = dplyr::n(), amp_ratio_mean = mean(.data$amp_ratio, na.rm = TRUE), .groups = "drop")
print(matched_tbl, n = 100)

# =============================================================================
# Plot: amp_ratio vs predicted p2p, one panel per condition
# =============================================================================
p <- ggplot(long, aes(x = .data$p2p_pred_pct, y = .data$amp_ratio, color = .data$stim_state)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_point(shape = 1, size = 0.6, alpha = 0.25, stroke = 0.3) +
  scale_color_manual(values = c("no stim" = "grey60", "left stim" = "#1d4ed8", "right stim" = "#b91c1c")) +
  facet_wrap(~condition, ncol = 3) +
  coord_cartesian(ylim = c(0, 10)) +
  labs(title = "Amplitude ratio (measured p2p / predicted p2p) vs predicted amplitude, before vs after de-staircasing sono",
       subtitle = "All dynamic trials, bass16/17/18 pooled. If de-staircasing collapses the low-amplitude blowup, the excess was a sono signal-conditioning artifact.",
       x = "Predicted p2p (%, encoder)", y = "Amplitude ratio (measured/predicted), clipped at 10 for readability",
       color = "Stim state") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
fout <- file.path(OUT_DIR, "amp_ratio_destaircased_comparison.png")
ggplot2::ggsave(fout, p, width = 12, height = 4.8, dpi = 150)
cli::cli_alert_success("Saved {fout}")
cli::cli_alert_success("diag_sono_amp_ratio_destaircased.R complete")
