# diag_sono_vibration_check.R
# Diagnostic (not a pipeline deliverable): tests the MECHANICAL VIBRATION
# hypothesis for the persistent multiplicative gain (~1.15x even after
# de-staircasing, see diag_sono_amp_ratio_destaircased.R) between predicted
# (encoder) and measured (sono) strain in dynamic trials -- i.e. that the
# sono crystal leads / DS3 housing / clamp assembly physically vibrate under
# the dynamic protocol's continuous sinusoidal acceleration reversals (which
# isometric, with no motor motion, and isovelocity, one long constant-
# velocity ramp, don't produce), adding real mechanical noise to the sono
# reading that isn't present in the encoder-derived prediction.
#
# Two checks:
#  1. Time-series overlay of a few 5-10%-amplitude-bin cycles (the bin that
#     holds the bulk of the dynamic data and the clearest departure from
#     1:1) across different frequencies -- visual: does the sono trace show
#     higher-frequency ripple riding on top of the true cycle shape that
#     the encoder-predicted trace does not?
#  2. Per-cycle regression: does the RESIDUAL (measured p2p beyond what the
#     baseline gain-only fit predicts) correlate with peak commanded
#     ANGULAR VELOCITY or ACCELERATION within that cycle, independent of
#     amplitude? Vibration excited by motion should scale with velocity/
#     acceleration, not just with cycle amplitude -- and since dynamic
#     trials span a RANGE of frequencies at similar amplitudes, velocity/
#     acceleration and amplitude are not perfectly collinear here.
#
# Sono is DE-STAIRCASED (lin-interp condition, see diag_sono_amp_ratio_
# destaircased.R) before computing p2p, to isolate the gain component from
# the already-explained staircase/ADC-noise floor.
#
# Run with:  Rscript R/diag_sono_vibration_check.R
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

.p2p_per_cycle <- function(x, ci_vec, cycles) {
  vapply(cycles, function(c) {
    v <- x[ci_vec == c]
    if (sum(is.finite(v)) < 4L) return(NA_real_)
    max(v, na.rm = TRUE) - min(v, na.rm = TRUE)
  }, numeric(1L))
}

#' Peak |commanded angular velocity| (deg/s) and peak |angular acceleration|
#' (deg/s^2, finite-difference of velocity at SAMPLE_RATE_HZ) within a cycle.
.peak_kinematics_per_cycle <- function(anglevel_degps, ci_vec, cycles, sample_rate_hz) {
  accel_degps2 <- c(0, diff(anglevel_degps)) * sample_rate_hz
  peak_vel   <- vapply(cycles, function(c) { v <- anglevel_degps[ci_vec == c]; if (sum(is.finite(v)) < 4L) return(NA_real_); max(abs(v), na.rm = TRUE) }, numeric(1L))
  peak_accel <- vapply(cycles, function(c) { v <- accel_degps2[ci_vec == c]; if (sum(is.finite(v)) < 4L) return(NA_real_); max(abs(v), na.rm = TRUE) }, numeric(1L))
  list(peak_vel_degps = peak_vel, peak_accel_degps2 = peak_accel)
}

.dominant_stim_state <- function(state_chr) {
  state_chr <- state_chr[!is.na(state_chr)]
  if (length(state_chr) == 0L) return(NA_character_)
  tab <- table(state_chr)
  names(tab)[which.max(tab)]
}

#' Full per-trial processing: returns per-cycle tibble (amp + kinematics) AND
#' (only for trials matching `keep_samples_for`) the full per-sample td for
#' the time-series overlay panels.
.process_trial <- function(specimen, h5path, keep_samples = FALSE) {
  trial_label <- tools::file_path_sans_ext(basename(h5path))
  td <- tryCatch(load_bender_flat(h5path, do_filter = FALSE, loadtorques = "x"), error = function(e) NULL)
  if (is.null(td) || nrow(td) == 0L) return(list(cycles = tibble(), td = NULL))

  geom <- .read_geom_and_lidx(h5path)
  if (!is.finite(geom$width_mm) || !is.finite(geom$lidx_right) || !is.finite(geom$lidx_pos_motor)) return(list(cycles = tibble(), td = NULL))

  info <- .read_cycle_freq(h5path)
  if (is.null(info$ci) || length(info$ci) != nrow(td)) return(list(cycles = tibble(), td = NULL))
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
  if (is.null(sono_raw_mm) || length(sono_raw_mm) != nrow(td)) return(list(cycles = tibble(), td = NULL))
  L0 <- .sono_reference_length_mm(td$angle.deg, sono_raw_mm)
  if (!is.finite(L0) || L0 <= 0) return(list(cycles = tibble(), td = NULL))

  mm_lininterp <- .lin_interp_destaircase(sono_raw_mm, td$t.s, geom$ds3_rate_hz, SAMPLE_RATE_HZ)
  # NOTE: force_sign_right folds ONLY the predicted (encoder) side, matching
  # attach_sono_strain()'s vetted convention (plot_angle_sono_validation.R).
  # The sono side is a direct physical length-change measurement and must
  # NOT also get force_sign_right -- doing so silently sign-flips the trace
  # whenever force_sign_right == -1, which mimics a large phase lag in the
  # time-series overlay below (Part 1) without being one. p2p (Part 2) was
  # unaffected since peak-to-peak is sign-invariant.
  td$strain_sono_raw_pct       <- -(sono_raw_mm - L0) / L0 * 100.0
  td$strain_sono_lininterp_pct <- -(mm_lininterp - L0) / L0 * 100.0

  stim_state_fac <- tryCatch(.stim_window_state_label(td), error = function(e) NULL)
  td$stim_state <- if (is.null(stim_state_fac)) NA_character_ else as.character(stim_state_fac)

  cyc_list <- sort(unique(td$cycle_idx[td$cycle_idx >= 1L]))
  if (length(cyc_list) == 0L) return(list(cycles = tibble(), td = if (keep_samples) td else NULL))

  p2p_pred <- .p2p_per_cycle(td$strain_pred_encoder_right_pct, td$cycle_idx, cyc_list)
  p2p_meas <- .p2p_per_cycle(td$strain_sono_lininterp_pct,     td$cycle_idx, cyc_list)
  kin <- .peak_kinematics_per_cycle(td$anglevel.degps, td$cycle_idx, cyc_list, SAMPLE_RATE_HZ)
  freq_by_c <- vapply(cyc_list, function(c) { v <- td$freq_cycle[td$cycle_idx == c]; v <- v[is.finite(v)]; if (length(v) == 0L) NA_real_ else v[1L] }, numeric(1L))
  dom_state <- vapply(cyc_list, function(c) .dominant_stim_state(td$stim_state[td$cycle_idx == c]), character(1L))

  cli::cli_alert_success("{trial_label}: {length(cyc_list)} cycles")

  cycles_tbl <- tibble(
    specimen = specimen, trial = trial_label, cycle = cyc_list, freq_hz = freq_by_c, stim_state = dom_state,
    p2p_pred_pct = p2p_pred, p2p_meas_pct = p2p_meas,
    peak_vel_degps = kin$peak_vel_degps, peak_accel_degps2 = kin$peak_accel_degps2
  )
  list(cycles = cycles_tbl, td = if (keep_samples) td else NULL)
}

# =============================================================================
# Part 1: time-series overlay for visual ripple check
# =============================================================================
cli::cli_h1("Part 1: time-series overlay (5-10% amplitude bin, several frequencies)")

TRIAL_17_BASS16 <- file.path(raw_source_dir(BASS16_RAW_SUBFOLDER), "2026-07-14_bass16_bender_17_dynamic.h5")
res17 <- .process_trial("bass16", TRIAL_17_BASS16, keep_samples = TRUE)
td17 <- res17$td
cyc17 <- res17$cycles |> dplyr::filter(.data$p2p_pred_pct >= 5, .data$p2p_pred_pct <= 10)

overlay_panels <- list()
if (!is.null(td17) && nrow(cyc17) > 0) {
  for (fhz in sort(unique(cyc17$freq_hz))) {
    cand <- cyc17 |> dplyr::filter(abs(.data$freq_hz - fhz) < 0.01)
    if (nrow(cand) == 0L) next
    pick_cycle <- cand$cycle[which.min(abs(cand$p2p_pred_pct - 7.5))]  # near bin center
    sl <- td17 |> dplyr::filter(.data$cycle_idx == pick_cycle,
                                 is.finite(.data$strain_pred_encoder_right_pct),
                                 is.finite(.data$strain_sono_raw_pct))
    if (nrow(sl) < 10L) next
    t0 <- min(sl$t.s, na.rm = TRUE)
    overlay_panels[[length(overlay_panels) + 1L]] <- tibble(
      freq_hz = fhz, cycle = pick_cycle, t = sl$t.s - t0,
      `Predicted (encoder)` = sl$strain_pred_encoder_right_pct,
      `Raw sono (staircase)` = sl$strain_sono_raw_pct,
      `Sono, de-staircased` = sl$strain_sono_lininterp_pct
    )
  }
}
overlay_df <- dplyr::bind_rows(overlay_panels)

if (nrow(overlay_df) > 0) {
  overlay_long <- overlay_df |>
    tidyr::pivot_longer(cols = c(`Predicted (encoder)`, `Raw sono (staircase)`, `Sono, de-staircased`),
                         names_to = "signal", values_to = "strain_pct") |>
    dplyr::mutate(signal = factor(signal, levels = c("Predicted (encoder)", "Raw sono (staircase)", "Sono, de-staircased")),
                  panel = sprintf("%.0f Hz (cycle %d)", freq_hz, cycle))

  p_overlay <- ggplot(overlay_long, aes(x = t, y = strain_pct, color = signal, linetype = signal, linewidth = signal)) +
    geom_line() +
    scale_color_manual(values = c("Predicted (encoder)" = "black", "Raw sono (staircase)" = "grey65", "Sono, de-staircased" = "#b91c1c")) +
    scale_linetype_manual(values = c("Predicted (encoder)" = "dashed", "Raw sono (staircase)" = "solid", "Sono, de-staircased" = "solid")) +
    scale_linewidth_manual(values = c("Predicted (encoder)" = 0.9, "Raw sono (staircase)" = 0.5, "Sono, de-staircased" = 0.9)) +
    facet_wrap(~panel, scales = "free", ncol = 1) +
    labs(title = "Bass16 trial 17 (dynamic): one 5-10%-amplitude-bin cycle per frequency block",
         subtitle = "Looking for high-frequency ripple on the sono trace that the encoder-predicted trace does not have (mechanical-vibration signature)",
         x = "Time within cycle (s)", y = "Right-muscle strain (%, shortening positive)", color = NULL, linetype = NULL, linewidth = NULL) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")

  fout1 <- file.path(OUT_DIR, "sono_vibration_timeseries_overlay.png")
  ggplot2::ggsave(fout1, p_overlay, width = 9, height = 3.2 * length(unique(overlay_long$panel)) + 1.5, dpi = 150)
  cli::cli_alert_success("Saved {fout1}")
} else {
  cli::cli_alert_warning("No overlay panels assembled -- skipping Part 1 plot")
}

# =============================================================================
# Part 2: does the residual scale with velocity/acceleration, not just amplitude?
# =============================================================================
cli::cli_h1("Part 2: batch over all dynamic trials, peak velocity/acceleration per cycle")

all_files <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  tibble(specimen = specimen, h5path = .dynamic_trial_files(subfolder))
})
all_cycles <- purrr::pmap_dfr(all_files, function(specimen, h5path) {
  tryCatch(.process_trial(specimen, h5path, keep_samples = FALSE)$cycles,
           error = function(e) { cli::cli_alert_danger("{basename(h5path)}: {conditionMessage(e)}"); tibble() })
})

d <- all_cycles |>
  dplyr::filter(is.finite(.data$p2p_pred_pct), .data$p2p_pred_pct >= MIN_P2P_PRED_PCT,
                is.finite(.data$p2p_meas_pct), is.finite(.data$peak_vel_degps), is.finite(.data$peak_accel_degps2))
write.csv(d, file.path(DATA_OUT_DIR, "sono_vibration_check_percycle.csv"), row.names = FALSE)
cli::cli_alert_success("Saved {nrow(d)} usable cycles -> data_processed/sono_vibration_check_percycle.csv")

cli::cli_h2("Baseline model: p2p_meas ~ p2p_pred")
fit_base <- lm(p2p_meas_pct ~ p2p_pred_pct, d)
d$resid_base <- residuals(fit_base)
print(summary(fit_base)$coefficients)
cat("R^2 =", summary(fit_base)$r.squared, "\n")

cli::cli_h2("cor(residual, peak_vel_degps) and cor(residual, peak_accel_degps2)")
cat("cor with peak velocity:     ", round(cor(d$resid_base, d$peak_vel_degps, use = "complete.obs"), 4), "\n")
cat("cor with peak acceleration: ", round(cor(d$resid_base, d$peak_accel_degps2, use = "complete.obs"), 4), "\n")

cli::cli_h2("Does adding velocity/acceleration improve on amplitude alone?")
fit_vel   <- lm(p2p_meas_pct ~ p2p_pred_pct + peak_vel_degps, d)
fit_accel <- lm(p2p_meas_pct ~ p2p_pred_pct + peak_accel_degps2, d)
fit_both  <- lm(p2p_meas_pct ~ p2p_pred_pct + peak_vel_degps + peak_accel_degps2, d)
cli::cli_alert_info("R^2 base={round(summary(fit_base)$r.squared,4)}  +vel={round(summary(fit_vel)$r.squared,4)}  +accel={round(summary(fit_accel)$r.squared,4)}  +both={round(summary(fit_both)$r.squared,4)}")
cat("\n-- +velocity model --\n"); print(summary(fit_vel)$coefficients)
cat("\n-- +acceleration model --\n"); print(summary(fit_accel)$coefficients)

cli::cli_h2("Per-frequency-block check: does amp_ratio grow with frequency AT FIXED amplitude bin?")
d$amp_bin <- cut(d$p2p_pred_pct, breaks = c(0, 2, 5, 10, 20, 40), labels = c("0-2%","2-5%","5-10%","10-20%","20-40%"))
d$amp_ratio <- d$p2p_meas_pct / d$p2p_pred_pct
freq_within_bin <- d |>
  dplyr::filter(!is.na(.data$amp_bin)) |>
  dplyr::group_by(.data$amp_bin, freq_hz = round(.data$freq_hz, 1)) |>
  dplyr::summarise(n = dplyr::n(), amp_ratio_mean = mean(.data$amp_ratio, na.rm = TRUE),
                    peak_accel_mean = mean(.data$peak_accel_degps2, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(.data$amp_bin, freq_hz)
print(freq_within_bin, n = 100)
write.csv(freq_within_bin, file.path(DATA_OUT_DIR, "sono_vibration_check_freq_within_ampbin.csv"), row.names = FALSE)

# =============================================================================
# Plot: residual vs peak acceleration / velocity
# =============================================================================
p_resid <- ggplot(d, aes(x = .data$peak_accel_degps2, y = .data$resid_base, color = .data$specimen)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(alpha = 0.35, size = 1.0) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
  labs(title = "Residual from amplitude-only fit vs. peak commanded angular acceleration",
       subtitle = "If mechanical vibration drives the residual gain, expect a positive trend here (residual grows with acceleration, not just amplitude)",
       x = "Peak |angular acceleration| within cycle (deg/s^2)", y = "Residual: measured p2p - amplitude-only fitted p2p (%)",
       color = "Specimen") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
fout2 <- file.path(OUT_DIR, "sono_vibration_residual_vs_acceleration.png")
ggplot2::ggsave(fout2, p_resid, width = 8, height = 5.5, dpi = 150)
cli::cli_alert_success("Saved {fout2}")

cli::cli_alert_success("diag_sono_vibration_check.R complete")
