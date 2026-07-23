# diag_sono_phase_lag_multitrial.R
# Diagnostic (not a pipeline deliverable): quantifies the phase/timing lag
# between encoder-predicted strain and sono-measured strain across EVERY
# dynamic trial in bass16/17/18 -- generalizing diag_sono_timing.R's
# cross-correlation approach (previously scoped to 2 bass16 trials and never
# run to completion/logged) to the full corpus, and adding a stim-state
# breakdown diag_sono_timing.R didn't have.
#
# Motivated by diag_sono_vibration_check.R's time-series overlay: predicted
# and sono traces are visibly out of sync (not just noisy) at 1-3 Hz, a
# separate phenomenon from the phase-BLIND per-cycle amplitude gain already
# quantified (peak-to-peak doesn't care about timing). This script measures
# that lag directly instead of eyeballing it.
#
# Per-cycle cross-correlation (search window = +/- half that cycle's period,
# same as diag_sono_timing.R), using the DE-STAIRCASED (lin-interp) sono
# strain -- the raw staircase's quantization can blur the cross-correlation
# peak and bias the lag estimate.
#
# Three questions:
#  1. Fixed-ms (hardware/DS3 latency) or fixed-%-of-cycle (something else)?
#     -- regress lag_ms against freq_hz; also inspect lag_pct_cycle vs freq.
#  2. Does STIMULATION add extra lag (excitation-contraction coupling delay,
#     a real biological signature) on top of whatever baseline/passive lag
#     exists? -- compare dominant-stim-state cycles' lag to no-stim cycles',
#     matched by frequency.
#  3. How much of the RAW POINTWISE fit (not the phase-blind p2p metric)
#     does lag-correction recover? -- shift sono by the estimated lag and
#     recompute pointwise r/RMSE vs. the uncorrected fit.
#
# Run with:  Rscript R/diag_sono_phase_lag_multitrial.R
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

SPECIMEN_SUBFOLDERS <- c(
  bass16 = BASS16_RAW_SUBFOLDER,
  bass17 = BASS17_RAW_SUBFOLDER,
  bass18 = BASS18_RAW_SUBFOLDER
)

# =============================================================================
# Shared helpers (same pattern as the other diag_sono_* scripts)
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

.dominant_stim_state <- function(state_chr) {
  state_chr <- state_chr[!is.na(state_chr)]
  if (length(state_chr) == 0L) return(NA_character_)
  tab <- table(state_chr)
  names(tab)[which.max(tab)]
}

#' Cross-correlate two finite numeric vectors; return lag at max cor.
#' Search window +/- max_lag_samps; positive lag = y lags behind x.
.xcorr_lag <- function(x, y, max_lag_samps) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 20L) return(list(lag = NA_integer_, peak_cor = NA_real_))
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  # max_lag_samps is derived from the cycle's NOMINAL period; if this
  # particular cycle lost samples to NA-filtering (short high-freq cycles
  # are most exposed), the nominal half-period can exceed the ACTUAL
  # finite-sample count n, driving n-lag negative. Cap it so at least a
  # third of the segment always overlaps.
  max_lag_samps <- min(max_lag_samps, max(1L, floor(n / 3)))
  lags <- seq(-max_lag_samps, max_lag_samps)
  cors <- vapply(lags, function(lag) {
    if (lag >= 0) {
      xi <- x[seq_len(n - lag)]; yi <- y[(lag + 1):n]
    } else {
      xi <- x[(-lag + 1):n]; yi <- y[seq_len(n + lag)]
    }
    if (length(xi) < 10L) return(NA_real_)
    suppressWarnings(cor(xi, yi, use = "complete.obs"))
  }, numeric(1L))
  best_idx <- which.max(cors)
  if (length(best_idx) == 0L) return(list(lag = NA_integer_, peak_cor = NA_real_))
  list(lag = lags[best_idx], peak_cor = cors[best_idx])
}

#' Process one dynamic trial -> per-cycle lag/amplitude tibble.
.process_trial <- function(specimen, h5path) {
  trial_label <- tools::file_path_sans_ext(basename(h5path))
  td <- tryCatch(load_bender_flat(h5path, do_filter = FALSE, loadtorques = "x"), error = function(e) NULL)
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

  td <- attach_predicted_strain(td, local_body_width_mm = geom$width_mm, measured_muscle_depth_mm = geom$depth_mm,
                                 active_mask = rep(FALSE, nrow(td)))
  td <- attach_measured_strain(td)
  force_sign_right <- geom$lidx_right * geom$lidx_pos_motor
  td$strain_pred_encoder_right_pct <- force_sign_right * td$strain_measured_pct

  sono_raw_mm <- tryCatch(.read_sono_right_mm_aligned(h5path), error = function(e) NULL)
  if (is.null(sono_raw_mm) || length(sono_raw_mm) != nrow(td)) return(tibble())
  L0 <- .sono_reference_length_mm(td$angle.deg, sono_raw_mm)
  if (!is.finite(L0) || L0 <= 0) return(tibble())
  mm_lininterp <- .lin_interp_destaircase(sono_raw_mm, td$t.s, geom$ds3_rate_hz, SAMPLE_RATE_HZ)
  # NOTE: force_sign_right folds ONLY the predicted (encoder) side, matching
  # attach_sono_strain()'s vetted convention (plot_angle_sono_validation.R).
  # The sono side is a direct physical length-change measurement and must
  # NOT also get force_sign_right -- doing so silently sign-flips the trace
  # whenever force_sign_right == -1, which mimics a ~half-cycle phase lag.
  td$strain_sono_lininterp_pct <- -(mm_lininterp - L0) / L0 * 100.0

  stim_state_fac <- tryCatch(.stim_window_state_label(td), error = function(e) NULL)
  td$stim_state <- if (is.null(stim_state_fac)) NA_character_ else as.character(stim_state_fac)

  cyc_list <- sort(unique(td$cycle_idx[td$cycle_idx >= 1L]))
  if (length(cyc_list) == 0L) return(tibble())

  rows <- vector("list", length(cyc_list))
  for (i in seq_along(cyc_list)) {
    c <- cyc_list[i]
    sub <- td[td$cycle_idx == c, , drop = FALSE]
    pred <- sub$strain_pred_encoder_right_pct
    meas <- sub$strain_sono_lininterp_pct
    fhz  <- sub$freq_cycle[is.finite(sub$freq_cycle)][1L]
    if (length(fhz) == 0L || !is.finite(fhz) || fhz <= 0) next
    p2p_pred <- if (sum(is.finite(pred)) >= 4L) max(pred, na.rm = TRUE) - min(pred, na.rm = TRUE) else NA_real_
    if (!is.finite(p2p_pred) || p2p_pred < MIN_P2P_PRED_PCT) next

    period_s   <- 1.0 / fhz
    max_lag_smp <- max(10L, as.integer(round(period_s * 0.5 * SAMPLE_RATE_HZ)))
    xc <- .xcorr_lag(pred, meas, max_lag_smp)
    if (is.na(xc$lag)) next
    lag_ms      <- xc$lag / SAMPLE_RATE_HZ * 1000.0
    lag_pct_cyc <- lag_ms * fhz / 1000.0 * 100.0

    rows[[i]] <- tibble(
      specimen = specimen, trial = trial_label, cycle = c, freq_hz = fhz,
      stim_state = .dominant_stim_state(sub$stim_state),
      p2p_pred_pct = p2p_pred, lag_samples = xc$lag, lag_ms = lag_ms,
      lag_pct_cycle = lag_pct_cyc, peak_xcorr = xc$peak_cor
    )
  }
  out <- dplyr::bind_rows(rows)
  cli::cli_alert_success("{trial_label}: {nrow(out)} cycles with lag estimate")
  out
}

# =============================================================================
# Batch
# =============================================================================
cli::cli_h1("Cross-correlation lag, per cycle, all dynamic trials bass16/17/18")

all_files <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  tibble(specimen = specimen, h5path = .dynamic_trial_files(subfolder))
})
lag_all <- purrr::pmap_dfr(all_files, function(specimen, h5path) {
  tryCatch(.process_trial(specimen, h5path),
           error = function(e) { cli::cli_alert_danger("{basename(h5path)}: {conditionMessage(e)}"); tibble() })
})
lag_all <- lag_all |> dplyr::mutate(stim_state = factor(.data$stim_state, levels = c("no stim", "left stim", "right stim")))
write.csv(lag_all, file.path(DATA_OUT_DIR, "sono_phase_lag_percycle.csv"), row.names = FALSE)
cli::cli_alert_success("Saved {nrow(lag_all)} cycles -> data_processed/sono_phase_lag_percycle.csv")

# =============================================================================
# Q1: fixed-ms or fixed-% of cycle?
# =============================================================================
cli::cli_h1("Q1: is lag fixed-ms (hardware latency) or fixed-% (something else)?")
by_freq <- lag_all |>
  dplyr::group_by(.data$freq_hz) |>
  dplyr::summarise(n = dplyr::n(), lag_ms_mean = mean(.data$lag_ms, na.rm = TRUE),
                    lag_ms_median = median(.data$lag_ms, na.rm = TRUE),
                    lag_pct_cycle_mean = mean(.data$lag_pct_cycle, na.rm = TRUE),
                    peak_xcorr_mean = mean(.data$peak_xcorr, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(.data$freq_hz)
print(by_freq, n = 100)
write.csv(by_freq, file.path(DATA_OUT_DIR, "sono_phase_lag_by_freq.csv"), row.names = FALSE)

fit_lag_freq <- lm(lag_ms ~ freq_hz, lag_all)
cli::cli_alert_info("lm(lag_ms ~ freq_hz): slope={round(coef(fit_lag_freq)[2],3)} ms per Hz, intercept={round(coef(fit_lag_freq)[1],3)} ms, R2={round(summary(fit_lag_freq)$r.squared,4)}")
cat("(slope ~0 -> lag_ms roughly CONSTANT across freq -> fixed-ms/hardware-latency-like.\n")
cat(" slope far from 0 -> lag_ms itself scales with freq -> NOT a simple fixed-ms delay.)\n\n")

# =============================================================================
# Q2: does stim add extra lag vs. no-stim, matched by frequency?
# =============================================================================
cli::cli_h1("Q2: does dominant stim state change the lag, at matched frequency?")
by_stim_freq <- lag_all |>
  dplyr::filter(!is.na(.data$stim_state)) |>
  dplyr::group_by(.data$freq_hz, .data$stim_state) |>
  dplyr::summarise(n = dplyr::n(), lag_ms_mean = mean(.data$lag_ms, na.rm = TRUE),
                    lag_ms_median = median(.data$lag_ms, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(.data$freq_hz, .data$stim_state)
print(by_stim_freq, n = 100)
write.csv(by_stim_freq, file.path(DATA_OUT_DIR, "sono_phase_lag_by_freq_stimstate.csv"), row.names = FALSE)

overall_by_stim <- lag_all |>
  dplyr::filter(!is.na(.data$stim_state)) |>
  dplyr::group_by(.data$stim_state) |>
  dplyr::summarise(n = dplyr::n(), lag_ms_mean = mean(.data$lag_ms, na.rm = TRUE),
                    lag_ms_median = median(.data$lag_ms, na.rm = TRUE),
                    lag_pct_cycle_mean = mean(.data$lag_pct_cycle, na.rm = TRUE), .groups = "drop")
cli::cli_h2("Pooled (all frequencies) lag by stim state:")
print(overall_by_stim)

# =============================================================================
# Q3: how much of the RAW POINTWISE fit does lag-correction recover?
# =============================================================================
cli::cli_h1("Q3: pointwise r/RMSE, raw vs. lag-corrected (using per-trial median lag)")

trial_lag <- lag_all |>
  dplyr::group_by(.data$specimen, .data$trial) |>
  dplyr::summarise(median_lag_ms = median(.data$lag_ms, na.rm = TRUE), .groups = "drop")

.recompute_pointwise_fit <- function(specimen, h5path, lag_ms_to_apply) {
  geom <- .read_geom_and_lidx(h5path)
  td <- tryCatch(load_bender_flat(h5path, do_filter = FALSE, loadtorques = "x"), error = function(e) NULL)
  if (is.null(td)) return(tibble())
  info <- .read_cycle_freq(h5path)
  if (is.null(info$ci) || length(info$ci) != nrow(td)) return(tibble())
  td$cycle_idx <- info$ci
  td <- attach_predicted_strain(td, local_body_width_mm = geom$width_mm, measured_muscle_depth_mm = geom$depth_mm, active_mask = rep(FALSE, nrow(td)))
  td <- attach_measured_strain(td)
  force_sign_right <- geom$lidx_right * geom$lidx_pos_motor
  td$strain_pred_encoder_right_pct <- force_sign_right * td$strain_measured_pct
  sono_raw_mm <- tryCatch(.read_sono_right_mm_aligned(h5path), error = function(e) NULL)
  if (is.null(sono_raw_mm) || length(sono_raw_mm) != nrow(td)) return(tibble())
  L0 <- .sono_reference_length_mm(td$angle.deg, sono_raw_mm)
  mm_lin <- .lin_interp_destaircase(sono_raw_mm, td$t.s, geom$ds3_rate_hz, SAMPLE_RATE_HZ)
  # See note in .process_trial(): force_sign_right belongs only on the
  # predicted side, not the sono side.
  td$strain_sono_lininterp_pct <- -(mm_lin - L0) / L0 * 100.0

  sub <- td |> dplyr::filter(.data$cycle_idx >= 1L, is.finite(.data$strain_pred_encoder_right_pct), is.finite(.data$strain_sono_lininterp_pct))
  if (nrow(sub) < 20L) return(tibble())
  pred <- sub$strain_pred_encoder_right_pct; meas <- sub$strain_sono_lininterp_pct
  r_raw <- suppressWarnings(cor(pred, meas)); rmse_raw <- sqrt(mean((pred - meas)^2))

  lag_smp <- as.integer(round(lag_ms_to_apply / 1000.0 * SAMPLE_RATE_HZ))
  n <- length(pred)
  r_lag <- NA_real_; rmse_lag <- NA_real_
  if (is.finite(lag_smp) && lag_smp > 0 && lag_smp < n - 10L) {
    pred_c <- pred[seq_len(n - lag_smp)]; meas_c <- meas[(lag_smp + 1L):n]
    r_lag <- suppressWarnings(cor(pred_c, meas_c)); rmse_lag <- sqrt(mean((pred_c - meas_c)^2))
  } else if (is.finite(lag_smp) && lag_smp < 0 && -lag_smp < n - 10L) {
    shift <- -lag_smp
    pred_c <- pred[(shift + 1L):n]; meas_c <- meas[seq_len(n - shift)]
    r_lag <- suppressWarnings(cor(pred_c, meas_c)); rmse_lag <- sqrt(mean((pred_c - meas_c)^2))
  }
  tibble(specimen = specimen, r_raw = r_raw, rmse_raw = rmse_raw, lag_ms_applied = lag_ms_to_apply, r_lag = r_lag, rmse_lag = rmse_lag)
}

pointwise_rows <- list()
for (i in seq_len(nrow(all_files))) {
  specimen <- all_files$specimen[i]; h5path <- all_files$h5path[i]
  trial_label <- tools::file_path_sans_ext(basename(h5path))
  lag_row <- trial_lag |> dplyr::filter(.data$specimen == !!specimen, .data$trial == !!trial_label)
  if (nrow(lag_row) == 0L || !is.finite(lag_row$median_lag_ms[1])) next
  res <- tryCatch(.recompute_pointwise_fit(specimen, h5path, lag_row$median_lag_ms[1]),
                  error = function(e) tibble())
  if (nrow(res) > 0) { res$trial <- trial_label; pointwise_rows[[length(pointwise_rows) + 1]] <- res }
}
pointwise_tbl <- dplyr::bind_rows(pointwise_rows)
write.csv(pointwise_tbl, file.path(DATA_OUT_DIR, "sono_phase_lag_pointwise_fit_comparison.csv"), row.names = FALSE)
cli::cli_h2("Pointwise fit: raw vs. lag-corrected (mean across trials)")
cli::cli_alert_info("Mean r_raw={round(mean(pointwise_tbl$r_raw, na.rm=TRUE),4)}  Mean r_lag_corrected={round(mean(pointwise_tbl$r_lag, na.rm=TRUE),4)}")
cli::cli_alert_info("Mean rmse_raw={round(mean(pointwise_tbl$rmse_raw, na.rm=TRUE),3)}%  Mean rmse_lag_corrected={round(mean(pointwise_tbl$rmse_lag, na.rm=TRUE),3)}%")

# =============================================================================
# Plots
# =============================================================================
p_lag_freq <- ggplot(lag_all, aes(x = .data$freq_hz, y = .data$lag_ms)) +
  geom_jitter(aes(color = .data$stim_state), width = 0.08, alpha = 0.3, size = 0.9) +
  stat_summary(fun = median, geom = "point", size = 3, color = "black") +
  stat_summary(fun = median, geom = "line", color = "black", linewidth = 0.7) +
  scale_color_manual(values = c("no stim" = "grey60", "left stim" = "#1d4ed8", "right stim" = "#b91c1c")) +
  labs(title = "Cross-correlation lag (sono relative to encoder-predicted) vs. cycle frequency",
       subtitle = "Black = per-frequency median across all cycles/trials/specimens. If flat -> fixed-ms hardware-like delay.",
       x = "Cycle frequency (Hz)", y = "Lag (ms, positive = sono lags encoder)", color = "Stim state") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggplot2::ggsave(file.path(OUT_DIR, "sono_phase_lag_vs_freq.png"), p_lag_freq, width = 8.5, height = 5.5, dpi = 150)

p_lag_stim <- ggplot(dplyr::filter(lag_all, !is.na(.data$stim_state)), aes(x = .data$stim_state, y = .data$lag_ms, fill = .data$stim_state)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_boxplot(outlier.alpha = 0.2) +
  scale_fill_manual(values = c("no stim" = "grey80", "left stim" = "#1d4ed8", "right stim" = "#b91c1c")) +
  labs(title = "Lag by dominant stim state (pooled across frequencies)",
       subtitle = "If stim adds extra lag (excitation-contraction coupling delay), stim boxes should sit above no-stim",
       x = NULL, y = "Lag (ms)") +
  theme_bw(base_size = 11) + theme(legend.position = "none")
ggplot2::ggsave(file.path(OUT_DIR, "sono_phase_lag_by_stimstate.png"), p_lag_stim, width = 6.5, height = 5, dpi = 150)

cli::cli_alert_success("Saved sono_phase_lag_vs_freq.png and sono_phase_lag_by_stimstate.png")
cli::cli_alert_success("diag_sono_phase_lag_multitrial.R complete")
