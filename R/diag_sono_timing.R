# diag_sono_timing.R
# Diagnostic (not a pipeline deliverable): determines WHY the dynamic panel of
# the sono validation identity plot shows a wide loop-shaped scatter (r~0.95)
# compared to isometric/isovelocity/frequency_sweep (which are much tighter).
#
# Hypothesis under test: the loop is a TIMING artifact from a fixed hardware
# lag between the DS3's internal segment-length computation and its analog
# output -- not a real strain-delivery error.  A fixed-ms lag produces a
# growing phase error (% cycle) at higher frequencies, which widens the
# instantaneous point-by-point scatter into a loop while leaving cycle-averaged
# amplitudes correct.
#
# Five diagnostic outputs:
#  1. Time-series overlay (predicted encoder strain vs sono strain, 2-3 cycles)
#     -- visual: same amplitude + phase-shifted vs different amplitude.
#  2. Cross-correlation lag per frequency block (ms and % cycle) across 0.5-5 Hz
#     -- tests fixed-ms (electronic latency) vs fixed-% (something else).
#  3. Amplitude ratio (measured p2p / predicted p2p per cycle) -- tests whether
#     amplitude is correct independent of the phase offset.
#  4. Isometric contrast -- near-zero expected lag (quasi-static).
#  5. 3-way identity refit: raw / lag-corrected / amplitude-only r and RMSE
#     -- if lag-correction collapses the loop, deviation was timing, not delivery.
#
# Key pipeline citations (read-only -- nothing below edits any existing file):
#  - Sono time base: same 1 kHz AI finite task as encoder; no resampling in
#    acquisition (bender_functions.py:2596-2607) or calibration
#    (R/01_calibrate.R:.cal_sono/.cal_sono_interp lines 135-188).
#  - DS3 internal rate (~247 Hz, daq_sono_internal_sample_rate_hertz attr) is
#    logged as metadata only; bender_functions.py:363 stores it in
#    self.sono_internal_rate but the field is never read again -- confirmed by
#    repo-wide grep.  Any DS3 hardware latency is therefore invisible to the
#    code and can only be measured empirically (this script).
#  - No alignment step exists: row i of encoder strain and row i of sono strain
#    are assumed simultaneous in R/plot_angle_sono_validation.R:attach_sono_strain
#    (lines 179-194) and in run_fv_fl_power_pipeline.R:414-433 (dynamic branch).
#
# Run with:  Rscript R/diag_sono_timing.R
# Outputs -> .diag_tmp/  (untracked, consistent with other diag_*.R scripts)

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2); library(cli); library(rhdf5)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")

# Raw-source default comes from paths_config.R (single source of truth) --
# see that file if the OneDrive folder layout ever moves again. OUT_DIR is
# local scratch, untracked, consistent with other diag_*.R scripts.
SOURCE_DIR <- raw_source_dir(BASS16_RAW_SUBFOLDER)
OUT_DIR    <- ".diag_tmp"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

src("00_load_bender_flat.R")
src("01_calibrate.R")
src("muscle_geometry.R")
src("plot_strain_validation.R")      # attach_measured_strain
src("plot_angle_sono_validation.R")  # attach_sono_strain / .read_sono_right_mm_aligned

SAMPLE_RATE_HZ <- 1000.0  # confirmed from metadata daq_ai_sample_rate_hertz

# Geometry (bass16, same specimen across all trials inspected)
BODY_WIDTH_MM   <- 42.0
MUSCLE_DEPTH_MM <- 0.0   # unset default -> resolve_muscle_depth_mm uses 1 mm fallback
DCLAMP_MM       <- 41.0  # measurement_clamp_separation_millimeter
LIDX_POS_MOTOR  <- -1.0  # daq_specimen_lateral_index_on_positive_motor_side
LIDX_LEFT       <- -1.0  # daq_specimen_side_index_left
LIDX_RIGHT      <-  1.0  # daq_specimen_side_index_right

FORCE_SIGN_RIGHT <- LIDX_RIGHT * LIDX_POS_MOTOR  # = -1; right shortens on neg angle

TRIAL_17 <- file.path(SOURCE_DIR, "2026-07-14_bass16_bender_17_dynamic.h5")
TRIAL_19 <- file.path(SOURCE_DIR, "2026-07-14_bass16_bender_19_dynamic.h5")
TRIAL_09 <- file.path(SOURCE_DIR, "2026-07-14_bass16_bender_09_isometric.h5")

# =============================================================================
# Helpers
# =============================================================================

#' Read cycle_index and frequency_hertz per sample (single_finite only).
.read_cycle_freq <- function(h5path) {
  h5 <- H5Fopen(h5path, "H5F_ACC_RDONLY")
  on.exit(try(H5Fclose(h5), silent = TRUE), add = TRUE)
  ci_raw  <- tryCatch(h5read(h5, "/timeseries/cycle_index"), error = function(e) NULL)
  freq_ds <- tryCatch(h5read(h5, "/metadata/index_cycle_frequency_hertz"), error = function(e) NULL)
  list(ci = as.integer(ci_raw), freq_by_cycle = as.numeric(freq_ds))
}

#' Attach predicted and sono strain to a single_finite td.
#' Returns td with strain_pred_encoder_right_pct and strain_sono_pct columns.
.attach_strains_dynamic <- function(td, h5path) {
  td <- attach_predicted_strain(
    td,
    local_body_width_mm      = BODY_WIDTH_MM,
    measured_muscle_depth_mm = MUSCLE_DEPTH_MM,
    active_mask              = rep(FALSE, nrow(td))   # geometry check, stim irrelevant
  )
  td <- attach_measured_strain(td)   # adds strain_measured_pct from enc.deg
  td <- attach_sono_strain(td, h5path, LIDX_RIGHT, LIDX_POS_MOTOR)
  td
}

#' Cross-correlate two finite numeric vectors; return lag at max |cor|.
#' Search window is +/- max_lag_samps; positive lag = y lags behind x.
.xcorr_lag <- function(x, y, max_lag_samps) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 20L) return(NA_integer_)
  x <- x[ok]; y <- y[ok]
  n  <- length(x)
  lags <- seq(-max_lag_samps, max_lag_samps)
  cors <- vapply(lags, function(lag) {
    if (lag >= 0) {
      xi <- x[seq_len(n - lag)]; yi <- y[(lag + 1):n]
    } else {
      xi <- x[(-lag + 1):n];     yi <- y[seq_len(n + lag)]
    }
    if (length(xi) < 10L) return(NA_real_)
    suppressWarnings(cor(xi, yi, use = "complete.obs"))
  }, numeric(1L))
  best_idx <- which.max(cors)
  if (length(best_idx) == 0L) return(NA_integer_)
  lag_at_max <- lags[best_idx]
  attr(lag_at_max, "cors") <- cors
  attr(lag_at_max, "lags") <- lags
  lag_at_max
}

#' Per-cycle peak-to-peak amplitude.
.p2p_per_cycle <- function(x, ci_vec, cycles) {
  vapply(cycles, function(c) {
    v <- x[ci_vec == c]
    if (sum(is.finite(v)) < 4L) return(NA_real_)
    max(v, na.rm = TRUE) - min(v, na.rm = TRUE)
  }, numeric(1L))
}

# =============================================================================
# Section 1 -- Load trial 17 and attach strains
# =============================================================================
cli::cli_h1("Loading trial 17 (dynamic, 1/2/3 Hz)")
td17 <- load_bender_flat(TRIAL_17, do_filter = FALSE, loadtorques = "x")
info17 <- .read_cycle_freq(TRIAL_17)
td17$cycle_idx <- info17$ci   # 1-based; -1 for pre-bending samples
td17$freq_cycle <- dplyr::if_else(
  td17$cycle_idx >= 1L & td17$cycle_idx <= length(info17$freq_by_cycle),
  info17$freq_by_cycle[pmax(1L, td17$cycle_idx)],
  NA_real_
)
cli::cli_alert_success("Loaded {nrow(td17)} samples; attaching strains...")
td17 <- .attach_strains_dynamic(td17, TRIAL_17)

# =============================================================================
# Section 2 -- Load trial 19 and attach strains (passive, 0.5/2.75/5 Hz)
# =============================================================================
cli::cli_h1("Loading trial 19 (dynamic passive, 0.5/2.75/5 Hz)")
td19 <- load_bender_flat(TRIAL_19, do_filter = FALSE, loadtorques = "x")
info19 <- .read_cycle_freq(TRIAL_19)
td19$cycle_idx <- info19$ci
td19$freq_cycle <- dplyr::if_else(
  td19$cycle_idx >= 1L & td19$cycle_idx <= length(info19$freq_by_cycle),
  info19$freq_by_cycle[pmax(1L, td19$cycle_idx)],
  NA_real_
)
td19 <- .attach_strains_dynamic(td19, TRIAL_19)

# =============================================================================
# Section 3 -- Load isometric trial 09 (segmented_finite; concat across steps)
# =============================================================================
cli::cli_h1("Loading trial 09 (isometric, segmented_finite)")
td09 <- load_bender_flat(TRIAL_09, do_filter = FALSE, loadtorques = "x")
# For isometric, attach_predicted_strain etc. need enc.deg / dclamp.m / muscle_strain_r_m
td09 <- attach_predicted_strain(
  td09,
  local_body_width_mm      = BODY_WIDTH_MM,
  measured_muscle_depth_mm = MUSCLE_DEPTH_MM
)
td09 <- attach_measured_strain(td09)
td09 <- attach_sono_strain(td09, TRIAL_09, LIDX_RIGHT, LIDX_POS_MOTOR)

# =============================================================================
# DIAGNOSTIC 1 -- Time-series overlay plots
# =============================================================================
cli::cli_h1("Diagnostic 1: time-series overlays")

# Helper: build 2-panel overlay plot for a slice
.plot_overlay <- function(df, title_str, freq_hz) {
  period_s <- 1.0 / freq_hz
  df <- df |> dplyr::filter(is.finite(.data$strain_pred_encoder_right_pct),
                             is.finite(.data$strain_sono_pct))
  if (nrow(df) < 10L) return(NULL)
  ggplot(df, aes(x = .data$t.s)) +
    geom_line(aes(y = .data$strain_pred_encoder_right_pct, color = "Predicted (encoder)"), linewidth = 0.7) +
    geom_line(aes(y = .data$strain_sono_pct, color = "Measured (sono)"), linewidth = 0.7, linetype = "dashed") +
    scale_color_manual(values = c("Predicted (encoder)" = "#1d4ed8", "Measured (sono)" = "#b91c1c"),
                       name = NULL) +
    labs(title = title_str,
         subtitle = sprintf("%.1f Hz | %d samples shown; predicted=encoder curvature*r, measured=sono right muscle",
                            freq_hz, nrow(df)),
         x = "Time (s)", y = "Right-muscle strain (%, shortening positive)") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
}

# Per frequency block in trial 17: show ~3 cycles centered on the middle of the block
for (fhz in c(1.0, 2.0, 3.0)) {
  period_s <- 1.0 / fhz
  sub <- td17 |> dplyr::filter(!is.na(.data$freq_cycle), abs(.data$freq_cycle - fhz) < 0.01)
  if (nrow(sub) < 20L) next
  # grab 3 cycles from the middle of the block
  t_mid  <- median(sub$t.s, na.rm = TRUE)
  t_lo   <- t_mid - 1.5 * period_s
  t_hi   <- t_mid + 1.5 * period_s
  slice  <- sub |> dplyr::filter(.data$t.s >= t_lo, .data$t.s <= t_hi)
  p <- .plot_overlay(slice, sprintf("Trial 17 -- %.0f Hz block, 3-cycle overlay", fhz), fhz)
  if (!is.null(p)) {
    fout <- file.path(OUT_DIR, sprintf("overlay_t17_%.0fHz.png", fhz))
    ggplot2::ggsave(fout, p, width = 9, height = 4, dpi = 150)
    cli::cli_alert_success("Saved {fout}")
  }
}

# Isometric: show one step (step_number 9 or first available active step)
if ("step_number" %in% names(td09)) {
  steps_available <- sort(unique(td09$step_number))
  for (sn in head(steps_available, 3L)) {
    sub09 <- td09 |> dplyr::filter(.data$step_number == sn,
                                   is.finite(.data$strain_pred_encoder_right_pct),
                                   is.finite(.data$strain_sono_pct))
    if (nrow(sub09) < 20L) next
    p09 <- ggplot(sub09, aes(x = .data$t.s)) +
      geom_line(aes(y = .data$strain_pred_encoder_right_pct, color = "Predicted (encoder)"), linewidth = 0.7) +
      geom_line(aes(y = .data$strain_sono_pct, color = "Measured (sono)"), linewidth = 0.7, linetype = "dashed") +
      scale_color_manual(values = c("Predicted (encoder)" = "#1d4ed8", "Measured (sono)" = "#b91c1c"),
                         name = NULL) +
      labs(title = sprintf("Trial 09 isometric -- step %d overlay", sn),
           subtitle = "Quasi-static hold: expect near-zero lag, tight tracking",
           x = "Time (s)", y = "Right-muscle strain (%, shortening positive)") +
      theme_bw(base_size = 11) + theme(legend.position = "bottom")
    ggplot2::ggsave(file.path(OUT_DIR, sprintf("overlay_t09_isometric_step%02d.png", sn)),
                   p09, width = 9, height = 4, dpi = 150)
    break  # just the first active one
  }
  cli::cli_alert_success("Saved isometric overlay")
}

# =============================================================================
# DIAGNOSTIC 2 & 3 -- Cross-correlation lag and amplitude ratio per freq block
# =============================================================================
cli::cli_h1("Diagnostic 2/3: lag estimation and amplitude ratio")

lag_rows <- list()

.estimate_block_lag <- function(td, h5path, trial_label, info) {
  rows <- list()
  unique_freqs <- sort(unique(info$freq_by_cycle[is.finite(info$freq_by_cycle)]))

  for (fhz in unique_freqs) {
    period_s    <- 1.0 / fhz
    max_lag_s   <- period_s * 0.5          # search +/- half period
    max_lag_smp <- as.integer(round(max_lag_s * SAMPLE_RATE_HZ))
    max_lag_smp <- max(max_lag_smp, 10L)   # at least 10 samples

    # gather all samples in this frequency block
    sub <- td |> dplyr::filter(!is.na(.data$freq_cycle), abs(.data$freq_cycle - fhz) < 0.01,
                               is.finite(.data$strain_pred_encoder_right_pct),
                               is.finite(.data$strain_sono_pct))
    if (nrow(sub) < 40L) {
      cli::cli_alert_warning("{trial_label} {fhz} Hz: too few samples ({nrow(sub)}), skip")
      next
    }

    pred <- sub$strain_pred_encoder_right_pct
    meas <- sub$strain_sono_pct

    # Cross-correlation: positive lag = meas lags pred (meas is delayed)
    lag_smp <- .xcorr_lag(pred, meas, max_lag_smp)
    if (is.na(lag_smp)) {
      cli::cli_alert_warning("{trial_label} {fhz} Hz: xcorr returned NA, skip")
      next
    }
    lag_ms     <- lag_smp / SAMPLE_RATE_HZ * 1000.0   # ms
    lag_pct_cyc <- lag_ms * fhz / 1000.0 * 100.0      # % of cycle

    # Amplitude ratio: p2p measured / p2p predicted, per cycle
    cyc_list <- sort(unique(sub$cycle_idx[sub$cycle_idx >= 1L]))
    p2p_pred <- .p2p_per_cycle(sub$strain_pred_encoder_right_pct, sub$cycle_idx, cyc_list)
    p2p_meas <- .p2p_per_cycle(sub$strain_sono_pct,               sub$cycle_idx, cyc_list)
    ok_amp   <- is.finite(p2p_pred) & is.finite(p2p_meas) & p2p_pred > 0
    amp_ratio <- if (sum(ok_amp) > 0) mean(p2p_meas[ok_amp] / p2p_pred[ok_amp], na.rm = TRUE) else NA_real_
    amp_ratio_sd <- if (sum(ok_amp) > 1) sd(p2p_meas[ok_amp] / p2p_pred[ok_amp], na.rm = TRUE) else NA_real_

    cli::cli_alert_info(
      "{trial_label} {fhz} Hz: lag={round(lag_ms,1)} ms ({round(lag_pct_cyc,1)}% cycle), amp_ratio={round(amp_ratio,3)} (sd={round(amp_ratio_sd,3)}), n_cyc={sum(ok_amp)}"
    )

    rows[[length(rows)+1]] <- tibble(
      trial = trial_label, freq_hz = fhz,
      lag_samples = lag_smp, lag_ms = lag_ms, lag_pct_cycle = lag_pct_cyc,
      amp_ratio_mean = amp_ratio, amp_ratio_sd = amp_ratio_sd,
      n_cycles = sum(ok_amp), n_samples = nrow(sub)
    )
  }
  if (length(rows) > 0) dplyr::bind_rows(rows) else tibble()
}

lag_t17 <- .estimate_block_lag(td17, TRIAL_17, "T17_dynamic", info17)
lag_t19 <- .estimate_block_lag(td19, TRIAL_19, "T19_dynamic_passive", info19)

# Isometric lag estimate: treat each step as a block
if ("step_number" %in% names(td09)) {
  iso_lag_rows <- list()
  for (sn in sort(unique(td09$step_number))) {
    sub <- td09 |>
      dplyr::filter(.data$step_number == sn,
                    is.finite(.data$strain_pred_encoder_right_pct),
                    is.finite(.data$strain_sono_pct))
    if (nrow(sub) < 40L) next
    lag_smp <- .xcorr_lag(sub$strain_pred_encoder_right_pct, sub$strain_sono_pct, 200L)
    if (is.na(lag_smp)) next
    lag_ms  <- lag_smp / SAMPLE_RATE_HZ * 1000.0
    p2p_p   <- max(sub$strain_pred_encoder_right_pct, na.rm = TRUE) - min(sub$strain_pred_encoder_right_pct, na.rm = TRUE)
    p2p_m   <- max(sub$strain_sono_pct, na.rm = TRUE)               - min(sub$strain_sono_pct, na.rm = TRUE)
    amp_r   <- if (is.finite(p2p_p) && p2p_p > 0) p2p_m / p2p_p else NA_real_
    iso_lag_rows[[length(iso_lag_rows)+1]] <- tibble(
      step = sn, lag_samples = lag_smp, lag_ms = lag_ms, amp_ratio = amp_r
    )
  }
  iso_lags <- if (length(iso_lag_rows) > 0) dplyr::bind_rows(iso_lag_rows) else tibble()
  if (nrow(iso_lags) > 0) {
    cli::cli_h2("Isometric lag summary (should be near zero):")
    print(iso_lags)
  }
}

lag_all <- dplyr::bind_rows(lag_t17, lag_t19)
if (nrow(lag_all) > 0) {
  cli::cli_h2("Combined lag table (dynamic trials):")
  print(lag_all)
}

# Lag vs frequency plot
if (nrow(lag_all) > 0) {
  p_lag <- ggplot(lag_all, aes(x = .data$freq_hz)) +
    geom_point(aes(y = .data$lag_ms, shape = .data$trial, color = .data$trial), size = 3) +
    geom_hline(yintercept = mean(lag_all$lag_ms, na.rm = TRUE),
               linetype = "dashed", color = "grey40") +
    labs(title = "Phase lag (ms) between encoder-predicted and sono strain vs. frequency",
         subtitle = "If lag is CONSTANT in ms across frequencies -> fixed hardware delay (DS3 latency), not delivery error\nIf lag grows proportionally -> something else",
         x = "Frequency (Hz)", y = "Lag (ms, positive = sono lags encoder)",
         shape = "Trial", color = "Trial") +
    theme_bw(base_size = 11)
  ggplot2::ggsave(file.path(OUT_DIR, "lag_vs_freq.png"), p_lag, width = 7, height = 4.5, dpi = 150)

  p_lag_pct <- ggplot(lag_all, aes(x = .data$freq_hz)) +
    geom_point(aes(y = .data$lag_pct_cycle, shape = .data$trial, color = .data$trial), size = 3) +
    labs(title = "Phase lag as % of cycle vs. frequency",
         subtitle = "If % cycle GROWS with frequency (but ms is constant) -> it is a fixed-delay timing artifact",
         x = "Frequency (Hz)", y = "Lag (% of cycle)",
         shape = "Trial", color = "Trial") +
    theme_bw(base_size = 11)
  ggplot2::ggsave(file.path(OUT_DIR, "lag_pct_vs_freq.png"), p_lag_pct, width = 7, height = 4.5, dpi = 150)
  cli::cli_alert_success("Saved lag_vs_freq.png and lag_pct_vs_freq.png")
}

# =============================================================================
# DIAGNOSTIC 5 -- 3-way identity fit: raw / lag-corrected / amplitude-only
# =============================================================================
cli::cli_h1("Diagnostic 5: 3-way identity refit (raw / lag-corrected / amplitude-only)")

.three_way_fit <- function(td, info, trial_label) {
  sub <- td |>
    dplyr::filter(!is.na(.data$freq_cycle),
                  is.finite(.data$strain_pred_encoder_right_pct),
                  is.finite(.data$strain_sono_pct))
  if (nrow(sub) < 20L) return(list())

  pred <- sub$strain_pred_encoder_right_pct
  meas <- sub$strain_sono_pct

  # 1) Raw
  r_raw   <- suppressWarnings(cor(pred, meas, use = "complete.obs"))
  rmse_raw <- sqrt(mean((pred - meas)^2, na.rm = TRUE))

  # Estimate a single representative lag: weighted average across frequency blocks
  best_lag_ms <- if (nrow(lag_all) > 0) {
    # Use median lag from whichever trial this is, fall back to all-trial median
    trial_prefix <- sub("_.*", "", trial_label)  # e.g. "T17" from "T17_dynamic"
    sub_lags <- lag_all |> dplyr::filter(grepl(trial_prefix, .data$trial, fixed = TRUE))
    if (nrow(sub_lags) == 0) median(lag_all$lag_ms, na.rm = TRUE) else median(sub_lags$lag_ms, na.rm = TRUE)
  } else NA_real_

  # 2) Lag-corrected
  r_lag <- NA_real_; rmse_lag <- NA_real_
  if (is.finite(best_lag_ms)) {
    lag_smp_cor <- as.integer(round(best_lag_ms / 1000.0 * SAMPLE_RATE_HZ))
    n <- nrow(sub)
    if (lag_smp_cor > 0 && lag_smp_cor < n - 10L) {
      pred_c <- pred[seq_len(n - lag_smp_cor)]
      meas_c <- meas[(lag_smp_cor + 1L):n]
      r_lag   <- suppressWarnings(cor(pred_c, meas_c, use = "complete.obs"))
      rmse_lag <- sqrt(mean((pred_c - meas_c)^2, na.rm = TRUE))
    } else if (lag_smp_cor < 0 && -lag_smp_cor < n - 10L) {
      shift <- -lag_smp_cor
      pred_c <- pred[(shift + 1L):n]
      meas_c <- meas[seq_len(n - shift)]
      r_lag   <- suppressWarnings(cor(pred_c, meas_c, use = "complete.obs"))
      rmse_lag <- sqrt(mean((pred_c - meas_c)^2, na.rm = TRUE))
    }
  }

  # 3) Amplitude-only: per-cycle p2p (one data point per cycle)
  cyc_list <- sort(unique(sub$cycle_idx[sub$cycle_idx >= 1L]))
  p2p_pred <- .p2p_per_cycle(pred, sub$cycle_idx, cyc_list)
  p2p_meas <- .p2p_per_cycle(meas, sub$cycle_idx, cyc_list)
  ok_amp <- is.finite(p2p_pred) & is.finite(p2p_meas) & p2p_pred > 0
  r_amp    <- suppressWarnings(cor(p2p_pred[ok_amp], p2p_meas[ok_amp], use = "complete.obs"))
  rmse_amp <- sqrt(mean((p2p_pred[ok_amp] - p2p_meas[ok_amp])^2, na.rm = TRUE))

  cli::cli_h2("3-way fit: {trial_label}")
  cli::cli_alert_info("  Raw:              r={round(r_raw,4)}, RMSE={round(rmse_raw,4)}%")
  cli::cli_alert_info("  Lag-corrected ({round(best_lag_ms,1)} ms shift): r={round(r_lag,4)}, RMSE={round(rmse_lag,4)}%")
  cli::cli_alert_info("  Amplitude-only ({sum(ok_amp)} cycles):  r={round(r_amp,4)}, RMSE={round(rmse_amp,4)}%")

  list(
    trial = trial_label,
    r_raw = r_raw, rmse_raw = rmse_raw,
    lag_applied_ms = best_lag_ms,
    r_lag = r_lag, rmse_lag = rmse_lag,
    r_amp = r_amp, rmse_amp = rmse_amp,
    n_points_raw = nrow(sub), n_cycles_amp = sum(ok_amp),
    sub = sub, p2p_pred = p2p_pred[ok_amp], p2p_meas = p2p_meas[ok_amp],
    best_lag_ms = best_lag_ms
  )
}

fit17 <- .three_way_fit(td17, info17, "T17_dynamic")
fit19 <- .three_way_fit(td19, info19, "T19_dynamic_passive")

# Scatter comparison plot: 3 panels side by side for trial 17
if (length(fit17) > 0 && !is.null(fit17$sub)) {
  sub17 <- fit17$sub
  lag_smp_use <- as.integer(round(fit17$best_lag_ms / 1000.0 * SAMPLE_RATE_HZ))
  n17 <- nrow(sub17)

  # raw scatter
  df_raw <- tibble(
    x = sub17$strain_pred_encoder_right_pct,
    y = sub17$strain_sono_pct,
    panel = "Raw"
  )
  # lag-corrected scatter
  df_lag <- if (!is.na(lag_smp_use) && abs(lag_smp_use) < n17 - 10L) {
    if (lag_smp_use >= 0) {
      tibble(x = sub17$strain_pred_encoder_right_pct[seq_len(n17 - lag_smp_use)],
             y = sub17$strain_sono_pct[(lag_smp_use+1L):n17],
             panel = "Lag-corrected")
    } else {
      shift <- -lag_smp_use
      tibble(x = sub17$strain_pred_encoder_right_pct[(shift+1L):n17],
             y = sub17$strain_sono_pct[seq_len(n17-shift)],
             panel = "Lag-corrected")
    }
  } else tibble(x = numeric(0), y = numeric(0), panel = character(0))

  # amplitude-only scatter
  df_amp <- tibble(
    x = fit17$p2p_pred,
    y = fit17$p2p_meas,
    panel = "Amplitude-only\n(per-cycle p2p)"
  )

  df_all <- dplyr::bind_rows(df_raw, df_lag, df_amp) |>
    dplyr::mutate(panel = factor(.data$panel,
                                 levels = c("Raw","Lag-corrected","Amplitude-only\n(per-cycle p2p)")))

  # annotation labels
  ann <- tibble::tribble(
    ~panel, ~label,
    "Raw",    sprintf("r=%.3f\nRMSE=%.2f%%\nn=%d pts", fit17$r_raw, fit17$rmse_raw, n17),
    "Lag-corrected", sprintf("r=%.3f\nRMSE=%.2f%%\nlag=%.0f ms", fit17$r_lag, fit17$rmse_lag, fit17$best_lag_ms),
    "Amplitude-only\n(per-cycle p2p)", sprintf("r=%.3f\nRMSE=%.2f%%\nn=%d cycles", fit17$r_amp, fit17$rmse_amp, fit17$n_cycles_amp)
  ) |>
    dplyr::mutate(panel = factor(.data$panel,
                                 levels = c("Raw","Lag-corrected","Amplitude-only\n(per-cycle p2p)")))

  p_3way <- ggplot(df_all, aes(x = .data$x, y = .data$y)) +
    geom_point(shape = 1, size = 0.8, alpha = 0.25, stroke = 0.3, color = "#1d4ed8") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.6) +
    geom_text(data = ann, aes(label = .data$label),
              x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = 2.9,
              inherit.aes = FALSE) +
    facet_wrap(~panel, scales = "free") +
    labs(title = "Trial 17 (dynamic, 1/2/3 Hz): 3-way sono vs encoder-predicted strain identity fit",
         subtitle = "Raw: all point-by-point samples | Lag-corrected: samples shifted by estimated hardware lag | Amplitude-only: 1 point per cycle (peak-to-peak measured vs predicted)",
         x = "Predicted strain (%, encoder right-folded)", y = "Measured strain (%, sono)") +
    theme_bw(base_size = 11)

  ggplot2::ggsave(file.path(OUT_DIR, "step8_detail.png"), p_3way, width = 13, height = 5.5, dpi = 150)
  cli::cli_alert_success("Saved step8_detail.png (3-way fit comparison)")
}

# Also do the same for trial 19
if (length(fit19) > 0 && !is.null(fit19$sub)) {
  sub19 <- fit19$sub
  lag_smp_use19 <- as.integer(round(fit19$best_lag_ms / 1000.0 * SAMPLE_RATE_HZ))
  n19 <- nrow(sub19)

  df_raw19 <- tibble(x=sub19$strain_pred_encoder_right_pct, y=sub19$strain_sono_pct, panel="Raw")
  df_amp19 <- tibble(x=fit19$p2p_pred, y=fit19$p2p_meas, panel="Amplitude-only\n(per-cycle p2p)")

  df_lag19 <- if (!is.na(lag_smp_use19) && abs(lag_smp_use19) < n19 - 10L) {
    if (lag_smp_use19 >= 0) {
      tibble(x=sub19$strain_pred_encoder_right_pct[seq_len(n19-lag_smp_use19)],
             y=sub19$strain_sono_pct[(lag_smp_use19+1L):n19], panel="Lag-corrected")
    } else {
      shift19 <- -lag_smp_use19
      tibble(x=sub19$strain_pred_encoder_right_pct[(shift19+1L):n19],
             y=sub19$strain_sono_pct[seq_len(n19-shift19)], panel="Lag-corrected")
    }
  } else tibble(x=numeric(0), y=numeric(0), panel=character(0))

  df_all19 <- dplyr::bind_rows(df_raw19, df_lag19, df_amp19) |>
    dplyr::mutate(panel = factor(.data$panel,
                                 levels = c("Raw","Lag-corrected","Amplitude-only\n(per-cycle p2p)")))
  ann19 <- tibble::tribble(
    ~panel, ~label,
    "Raw", sprintf("r=%.3f\nRMSE=%.2f%%", fit19$r_raw, fit19$rmse_raw),
    "Lag-corrected", sprintf("r=%.3f\nRMSE=%.2f%%\nlag=%.0f ms", fit19$r_lag, fit19$rmse_lag, fit19$best_lag_ms),
    "Amplitude-only\n(per-cycle p2p)", sprintf("r=%.3f\nRMSE=%.2f%%", fit19$r_amp, fit19$rmse_amp)
  ) |>
    dplyr::mutate(panel = factor(.data$panel,
                                 levels = c("Raw","Lag-corrected","Amplitude-only\n(per-cycle p2p)")))

  p_3way19 <- ggplot(df_all19, aes(x=.data$x, y=.data$y)) +
    geom_point(shape=1, size=0.8, alpha=0.25, stroke=0.3, color="#1d4ed8") +
    geom_abline(slope=1, intercept=0, linetype="dashed", linewidth=0.6) +
    geom_text(data=ann19, aes(label=.data$label), x=-Inf, y=Inf,
              hjust=-0.05, vjust=1.3, size=2.9, inherit.aes=FALSE) +
    facet_wrap(~panel, scales="free") +
    labs(title="Trial 19 (dynamic passive, 0.5/2.75/5 Hz): 3-way identity fit",
         x="Predicted strain (%, encoder right-folded)", y="Measured strain (%, sono)") +
    theme_bw(base_size=11)
  ggplot2::ggsave(file.path(OUT_DIR, "zoom_region.png"), p_3way19, width=13, height=5.5, dpi=150)
  cli::cli_alert_success("Saved zoom_region.png (trial 19 3-way fit)")
}

# =============================================================================
# Summary table
# =============================================================================
cli::cli_h1("Summary table")
if (nrow(lag_all) > 0) {
  cat("\n=== LAG + AMPLITUDE RATIO TABLE (dynamic trials) ===\n")
  print(lag_all |> dplyr::select(trial, freq_hz, lag_ms, lag_pct_cycle, amp_ratio_mean, amp_ratio_sd, n_cycles))
}
cat("\n=== 3-WAY FIT SUMMARY ===\n")
for (fit in list(fit17, fit19)) {
  if (length(fit) == 0) next
  cat(sprintf("\n%s:\n  Raw:              r=%.4f, RMSE=%.4f%%\n  Lag-corrected (%+.0f ms): r=%.4f, RMSE=%.4f%%\n  Amplitude-only:   r=%.4f, RMSE=%.4f%% (%d cycles)\n",
      fit$trial, fit$r_raw, fit$rmse_raw, fit$best_lag_ms, fit$r_lag, fit$rmse_lag,
      fit$r_amp, fit$rmse_amp, fit$n_cycles_amp))
}
cli::cli_alert_success("diag_sono_timing.R complete -- outputs in {OUT_DIR}/")
