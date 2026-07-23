# diag_sono_progressive_slip.R
# Diagnostic (not a pipeline deliverable): tests the PROGRESSIVE/CYCLING-
# DRIVEN SLIP hypothesis for the persistent multiplicative gain (~1.15-1.7x
# depending on amplitude bin) between predicted (encoder) and measured
# (sono) strain in dynamic trials.
#
# Already ruled out: STIM-driven slip (diag_sono_amp_ratio_multitrial.R --
# amp_ratio is the same for stim vs no-stim cycles at matched amplitude).
# This script asks a DIFFERENT question: does repeated motor OSCILLATION
# itself (independent of whether the muscle was stimulated) progressively
# loosen the clamp -- i.e. does the gain trend UPWARD:
#   (a) with cycle index WITHIN a single trial (slip accumulating over the
#       ~10-60 s of one dynamic run), or
#   (b) with trial ORDER within a specimen's session (slip accumulating
#       across the ~1-3 hour session as the specimen is cycled trial after
#       trial), chronological order approximated by the bender_NN trial
#       number parsed from the filename (interleaved isovelocity/isometric/
#       frequency-sweep trials still preserve chronological order in NN,
#       since NN increments every run regardless of protocol type)?
#
# Reuses diag_sono_vibration_check.R's .process_trial() pattern verbatim
# (de-staircased sono, force_sign_right folded ONLY on the predicted side --
# see that script's header for the sign-bug history) to get per-cycle
# p2p_pred_pct / p2p_meas_pct / freq_hz / stim_state / cycle / trial /
# specimen, then adds:
#   - trial_num (parsed from filename, chronological proxy)
#   - amp_ratio = p2p_meas_pct / p2p_pred_pct
#   - resid_base2: residual from lm(p2p_meas_pct ~ p2p_pred_pct + freq_hz),
#     i.e. controlling for BOTH amplitude and frequency (frequency is
#     included here, unlike vibration_check's resid_base, because
#     diag_sono_vibration_check.R already showed acceleration/velocity add
#     nothing once amplitude is in the model, but progressive slip is a
#     within-trial/within-session TIME effect that frequency composition
#     across trials could otherwise confound).
#
# Two tests (linear regression + correlation), each on resid_base2 AND on
# the raw amp_ratio for transparency:
#   (a) within-trial:    resid ~ cycle_idx,          per trial + pooled
#                         (pooled model includes trial as a fixed effect,
#                         i.e. lm(resid2 ~ cycle_idx + factor(trial)), so
#                         the cycle_idx coefficient is a within-trial slope
#                         net of any across-trial baseline differences)
#   (b) trial order:     trial-mean(resid) ~ trial_num, per specimen
#
# A positive, significant slope in either test supports progressive slip.
# A flat/null result in both rules it out (in addition to the already-ruled-
# out stim-driven-slip test).
#
# Run with:  Rscript R/diag_sono_progressive_slip.R
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

SAMPLE_RATE_HZ   <- 1000.0
MIN_P2P_PRED_PCT <- 0.5
MIN_CYCLES_PER_TRIAL_FOR_WITHIN_FIT <- 6L

SPECIMEN_SUBFOLDERS <- c(
  bass16 = BASS16_RAW_SUBFOLDER,
  bass17 = BASS17_RAW_SUBFOLDER,
  bass18 = BASS18_RAW_SUBFOLDER
)

# =============================================================================
# Shared helpers (identical pattern to diag_sono_vibration_check.R /
# diag_sono_amp_ratio_multitrial.R -- see those files for provenance notes)
# =============================================================================

.dynamic_trial_files <- function(subfolder) {
  d <- raw_source_dir(subfolder)
  sort(list.files(d, pattern = "_dynamic\\.h5$", full.names = TRUE))
}

#' Chronological trial number parsed from the bender_NN filename token.
#' Every acquisition run (dynamic, isovelocity, isometric, frequency_sweep)
#' increments NN, so this is a valid within-session chronological proxy even
#' though this script only analyzes the "_dynamic" subset.
.trial_num_from_label <- function(trial_label) {
  m <- regmatches(trial_label, regexpr("bender_(\\d+)_", trial_label))
  if (length(m) == 0L || !nzchar(m)) return(NA_integer_)
  as.integer(gsub("[^0-9]", "", m))
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

.dominant_stim_state <- function(state_chr) {
  state_chr <- state_chr[!is.na(state_chr)]
  if (length(state_chr) == 0L) return(NA_character_)
  tab <- table(state_chr)
  names(tab)[which.max(tab)]
}

#' Per-trial processing -- identical formula/sign-convention to
#' diag_sono_vibration_check.R's .process_trial() (de-staircased sono,
#' force_sign_right on predicted side ONLY). Returns one row per cycle.
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
  # See diag_sono_vibration_check.R header for the sign-bug history this
  # avoids re-introducing. p2p / amp_ratio are sign-invariant regardless.
  td$strain_sono_lininterp_pct <- -(mm_lininterp - L0) / L0 * 100.0

  stim_state_fac <- tryCatch(.stim_window_state_label(td), error = function(e) NULL)
  td$stim_state <- if (is.null(stim_state_fac)) NA_character_ else as.character(stim_state_fac)

  cyc_list <- sort(unique(td$cycle_idx[td$cycle_idx >= 1L]))
  if (length(cyc_list) == 0L) return(tibble())

  p2p_pred  <- .p2p_per_cycle(td$strain_pred_encoder_right_pct, td$cycle_idx, cyc_list)
  p2p_meas  <- .p2p_per_cycle(td$strain_sono_lininterp_pct,     td$cycle_idx, cyc_list)
  freq_by_c <- vapply(cyc_list, function(c) { v <- td$freq_cycle[td$cycle_idx == c]; v <- v[is.finite(v)]; if (length(v) == 0L) NA_real_ else v[1L] }, numeric(1L))
  dom_state <- vapply(cyc_list, function(c) .dominant_stim_state(td$stim_state[td$cycle_idx == c]), character(1L))

  cli::cli_alert_success("{trial_label}: {length(cyc_list)} cycles")

  tibble(
    specimen = specimen, trial = trial_label, trial_num = .trial_num_from_label(trial_label),
    cycle = cyc_list, freq_hz = freq_by_c, stim_state = dom_state,
    p2p_pred_pct = p2p_pred, p2p_meas_pct = p2p_meas
  )
}

# =============================================================================
# Batch over every dynamic trial, every specimen
# =============================================================================
cli::cli_h1("Batch: per-cycle amp_ratio + trial_num, all dynamic trials bass16/17/18")

all_files <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  tibble(specimen = specimen, h5path = .dynamic_trial_files(subfolder))
})
all_cycles <- purrr::pmap_dfr(all_files, function(specimen, h5path) {
  tryCatch(.process_trial(specimen, h5path),
           error = function(e) { cli::cli_alert_danger("{basename(h5path)}: {conditionMessage(e)}"); tibble() })
})

d <- all_cycles |>
  dplyr::filter(is.finite(.data$p2p_pred_pct), .data$p2p_pred_pct >= MIN_P2P_PRED_PCT,
                is.finite(.data$p2p_meas_pct), is.finite(.data$trial_num))
d$amp_ratio <- d$p2p_meas_pct / d$p2p_pred_pct

cli::cli_h2("Baseline model controlling for BOTH amplitude and frequency: p2p_meas ~ p2p_pred + freq_hz")
fit_base2 <- lm(p2p_meas_pct ~ p2p_pred_pct + freq_hz, d)
d$resid_base2 <- residuals(fit_base2)
print(summary(fit_base2)$coefficients)
cat("R^2 =", summary(fit_base2)$r.squared, "\n\n")

write.csv(d, file.path(DATA_OUT_DIR, "sono_progressive_slip_percycle.csv"), row.names = FALSE)
cli::cli_alert_success("Saved {nrow(d)} usable cycles -> data_processed/sono_progressive_slip_percycle.csv")

# =============================================================================
# Test (a): within-trial cycle-index trend
# =============================================================================
cli::cli_h1("Test (a): does the gain trend UPWARD with cycle index WITHIN a trial?")

within_trial_slopes <- d |>
  dplyr::group_by(.data$specimen, .data$trial) |>
  dplyr::filter(dplyr::n() >= MIN_CYCLES_PER_TRIAL_FOR_WITHIN_FIT) |>
  dplyr::summarise(
    n_cycles = dplyr::n(),
    cor_resid_cycle = suppressWarnings(cor(.data$resid_base2, .data$cycle, use = "complete.obs")),
    slope_resid_cycle = tryCatch(coef(lm(resid_base2 ~ cycle, dplyr::cur_data()))[["cycle"]], error = function(e) NA_real_),
    p_resid_cycle = tryCatch(summary(lm(resid_base2 ~ cycle, dplyr::cur_data()))$coefficients["cycle", "Pr(>|t|)"], error = function(e) NA_real_),
    cor_ratio_cycle = suppressWarnings(cor(.data$amp_ratio, .data$cycle, use = "complete.obs")),
    slope_ratio_cycle = tryCatch(coef(lm(amp_ratio ~ cycle, dplyr::cur_data()))[["cycle"]], error = function(e) NA_real_),
    .groups = "drop"
  )
print(within_trial_slopes, n = 100)
write.csv(within_trial_slopes, file.path(DATA_OUT_DIR, "sono_progressive_slip_within_trial_slopes.csv"), row.names = FALSE)

n_trials_fit <- nrow(within_trial_slopes)
n_positive_resid_slope <- sum(within_trial_slopes$slope_resid_cycle > 0, na.rm = TRUE)
cli::cli_h2("Summary across {n_trials_fit} trials with >= {MIN_CYCLES_PER_TRIAL_FOR_WITHIN_FIT} cycles")
cat("Trials with POSITIVE within-trial resid~cycle slope: ", n_positive_resid_slope, " / ", n_trials_fit,
    " (", round(100 * n_positive_resid_slope / n_trials_fit, 1), "%)\n", sep = "")
cat("Mean within-trial resid~cycle slope:   ", signif(mean(within_trial_slopes$slope_resid_cycle, na.rm = TRUE), 4), " %-pt per cycle\n", sep = "")
cat("Median within-trial resid~cycle slope: ", signif(median(within_trial_slopes$slope_resid_cycle, na.rm = TRUE), 4), " %-pt per cycle\n", sep = "")
t_slopes <- tryCatch(t.test(within_trial_slopes$slope_resid_cycle), error = function(e) NULL)
if (!is.null(t_slopes)) {
  cat("One-sample t-test (slope != 0 across trials): t = ", round(t_slopes$statistic, 3),
      ", p = ", signif(t_slopes$p.value, 4), ", 95% CI [", signif(t_slopes$conf.int[1], 4), ", ", signif(t_slopes$conf.int[2], 4), "]\n\n", sep = "")
}

cli::cli_h2("Pooled model: resid_base2 ~ cycle + factor(trial) (within-trial slope net of trial-level baseline offsets)")
d_pooled <- d |> dplyr::group_by(.data$trial) |> dplyr::filter(dplyr::n() >= MIN_CYCLES_PER_TRIAL_FOR_WITHIN_FIT) |> dplyr::ungroup()
fit_pooled_cycle <- lm(resid_base2 ~ cycle + factor(trial), d_pooled)
pooled_cycle_coef <- summary(fit_pooled_cycle)$coefficients["cycle", ]
cli::cli_alert_info("Pooled within-trial cycle coefficient: {round(pooled_cycle_coef[['Estimate']], 5)} %-pt/cycle, SE={round(pooled_cycle_coef[['Std. Error']], 5)}, p={signif(pooled_cycle_coef[['Pr(>|t|)']], 4)}")

# =============================================================================
# Test (b): trial-order (session-chronology) trend
# =============================================================================
cli::cli_h1("Test (b): does the gain trend UPWARD with trial ORDER within a specimen's session?")

trial_means <- d |>
  dplyr::group_by(.data$specimen, .data$trial, .data$trial_num) |>
  dplyr::summarise(n_cycles = dplyr::n(),
                    amp_ratio_mean = mean(.data$amp_ratio, na.rm = TRUE),
                    resid_base2_mean = mean(.data$resid_base2, na.rm = TRUE),
                    .groups = "drop") |>
  dplyr::arrange(.data$specimen, .data$trial_num)
write.csv(trial_means, file.path(DATA_OUT_DIR, "sono_progressive_slip_trial_order.csv"), row.names = FALSE)
print(trial_means, n = 100)

trial_order_fits <- trial_means |>
  dplyr::group_by(.data$specimen) |>
  dplyr::filter(dplyr::n() >= 4L) |>
  dplyr::summarise(
    n_trials = dplyr::n(),
    cor_resid_trialnum = suppressWarnings(cor(.data$resid_base2_mean, .data$trial_num, use = "complete.obs")),
    slope_resid_trialnum = tryCatch(coef(lm(resid_base2_mean ~ trial_num, dplyr::cur_data()))[["trial_num"]], error = function(e) NA_real_),
    p_resid_trialnum = tryCatch(summary(lm(resid_base2_mean ~ trial_num, dplyr::cur_data()))$coefficients["trial_num", "Pr(>|t|)"], error = function(e) NA_real_),
    cor_ratio_trialnum = suppressWarnings(cor(.data$amp_ratio_mean, .data$trial_num, use = "complete.obs")),
    slope_ratio_trialnum = tryCatch(coef(lm(amp_ratio_mean ~ trial_num, dplyr::cur_data()))[["trial_num"]], error = function(e) NA_real_),
    p_ratio_trialnum = tryCatch(summary(lm(amp_ratio_mean ~ trial_num, dplyr::cur_data()))$coefficients["trial_num", "Pr(>|t|)"], error = function(e) NA_real_),
    .groups = "drop"
  )
cli::cli_h2("Per-specimen trial-order regression (trial-mean resid/ratio vs. bender_NN trial number)")
print(trial_order_fits, n = 100)
write.csv(trial_order_fits, file.path(DATA_OUT_DIR, "sono_progressive_slip_trial_order_fits.csv"), row.names = FALSE)

cli::cli_h2("Pooled (all specimens) model: resid_base2_mean ~ trial_num + factor(specimen)")
fit_pooled_trial <- tryCatch(lm(resid_base2_mean ~ trial_num + factor(specimen), trial_means), error = function(e) NULL)
if (!is.null(fit_pooled_trial)) {
  pooled_trial_coef <- summary(fit_pooled_trial)$coefficients["trial_num", ]
  cli::cli_alert_info("Pooled trial_num coefficient: {round(pooled_trial_coef[['Estimate']], 5)} %-pt/trial, SE={round(pooled_trial_coef[['Std. Error']], 5)}, p={signif(pooled_trial_coef[['Pr(>|t|)']], 4)}")
}

# =============================================================================
# Plots
# =============================================================================
cli::cli_h1("Plots")

p_cycle <- ggplot(d_pooled, aes(x = .data$cycle, y = .data$resid_base2, color = .data$specimen)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(alpha = 0.25, size = 0.9) +
  geom_smooth(aes(group = .data$trial), method = "lm", se = FALSE, linewidth = 0.35, alpha = 0.5) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linewidth = 1.0) +
  facet_wrap(~specimen, scales = "free_x") +
  labs(title = "Residual (amplitude- and frequency-controlled) vs. within-trial cycle index",
       subtitle = "Thin lines = per-trial fits, thick black = pooled fit. A progressive-slip signature should show a consistent positive slope.",
       x = "Cycle index within trial", y = "Residual: measured p2p - (amplitude+freq) fitted p2p (%)", color = "Specimen") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
fout_cycle <- file.path(OUT_DIR, "sono_progressive_slip_vs_cycle_index.png")
ggplot2::ggsave(fout_cycle, p_cycle, width = 10, height = 5.5, dpi = 150)
cli::cli_alert_success("Saved {fout_cycle}")

p_trial <- ggplot(trial_means, aes(x = .data$trial_num, y = .data$resid_base2_mean, color = .data$specimen)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
  facet_wrap(~specimen, scales = "free_x") +
  labs(title = "Trial-mean residual (amplitude+freq-controlled) vs. chronological trial order",
       subtitle = "Trial number parsed from bender_NN filename token; a session-progressive-slip signature should show a positive slope across the session",
       x = "Trial number (bender_NN, chronological)", y = "Trial-mean residual (%)", color = "Specimen") +
  theme_bw(base_size = 11) + theme(legend.position = "none")
fout_trial <- file.path(OUT_DIR, "sono_progressive_slip_vs_trial_order.png")
ggplot2::ggsave(fout_trial, p_trial, width = 10, height = 4.5, dpi = 150)
cli::cli_alert_success("Saved {fout_trial}")

cli::cli_alert_success("diag_sono_progressive_slip.R complete -- outputs in {OUT_DIR}/ and {DATA_OUT_DIR}/")
