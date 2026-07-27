# summary_fl_fv_tension_pooled_regression_baseline.R
# PI directive (2026-07-25): "let's draw FL and FV and isometric tension
# plots (figs_summary) using the pooled model data."
#
# CAVEAT communicated to the PI before building (see chat + analysis_
# muscle_force_vector_log.md 2026-07-25 addenda): the pooled position+time
# regression passive-baseline method was validated on exactly ONE step
# (bass17, step 16) and there it DISAGREED with the simpler pre_static/
# interp_linear methods -- opposite sign, ~6x magnitude. Using more data to
# fit it does not prove it is MORE ACCURATE, only that it uses more data.
# These are EXPLORATORY figures, NOT a replacement for the existing
# static-baseline canonical FL/FV (run_fv_fl_power_pipeline.R,
# plot_summary_profiles.R) or tension (diag_precondition_tension_vs_offset_
# isometric.R / summary_coughlin2000_bass_comparison.R) figures -- every
# title/subtitle below says so explicitly, and the filenames are suffixed
# "_pooledRegressionBaseline" so they cannot be confused with those.
#
# NAMING COLLISION WARNING: "pooled" here means the PASSIVE-BASELINE
# calculation method (regress passive torque on elapsed time + operating_
# point, POOLING all of one trial's own steps' pre/post windows) -- it does
# NOT mean "pooled across specimens/individuals" the way superplot_fl_
# pooled.R uses "pooled" (that file pools 3 DIFFERENT V=0 data SOURCES
# across fish using the vector-force method with F/F0 normalization -- an
# unrelated, pre-existing convention, not reused or modified here).
#
# METHOD: fit_trial_pooled_passive_model() -- per TRIAL (own fish, own
# session, own operating_point range) linear regression: passive torque ~
# elapsed real session time (wall_clock_start-based, since t.s resets per
# step and cannot be compared across steps -- caught as a bug earlier the
# same day) + operating_point (quadratic term, since passive torque vs.
# bend angle is not expected to be linear through zero). Trained on EVERY
# step's own pre-/post-baseline window mean in that SAME trial (all sides,
# all recruitment types -- more data per trial than any single step's own
# 2-window estimate). Predicts the passive reference at each step's own
# active-window elapsed time + operating_point; used in place of
# passive_force_Nm_static everywhere below (isometric AND isovelocity, for
# a consistent method across both -- isovelocity's own better-motivated
# velocity_matched baseline is deliberately NOT used here, so all 3 outputs
# share one baseline method).
#
# Reuses the EXISTING, unmodified plot builders (build_summary_plot_
# isometric()/_isovelocity(), plot_summary_profiles.R) via the documented
# "pass a muscle_force_Nm-substituted step_summary plus a distinguishing
# title" pattern (build_summary_plot_isometric()'s own docstring, 2026-07-16)
# -- exactly how the existing "_baselineInterp" variant is produced.
#
# Run with:  Rscript R/summary_fl_fv_tension_pooled_regression_baseline.R
# Outputs -> 02_processed/figs_summary/:
#   FL_isometric_pooledRegressionBaseline.png
#   FV_isovelocity_pooledRegressionBaseline.png
#   isometricTension_pooledRegressionBaseline.png

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(ggplot2); library(cli)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")
src("00_load_bender_flat.R")
src("01_calibrate.R")
src("02_deconvolve.R")
src("muscle_geometry.R")
src("plot_force_vs_time.R")
src("fit_fv_fl.R")             # fit_force_velocity_curve() -- REQUIRED by the isovelocity Hill refit below
src("03_analyze.R")
src("plot_summary_profiles.R") # build_summary_plot_isometric()/_isovelocity(), UNMODIFIED
src("parse_trial_filename.R")
src("dynamic_trial_precondition.R")  # extract_bender_trial_num()

OUT_DIR <- FIGS_SUMMARY_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SPECIMEN_SUBFOLDERS <- c(bass16 = BASS16_RAW_SUBFOLDER, bass17 = BASS17_RAW_SUBFOLDER, bass18 = BASS18_RAW_SUBFOLDER)
COUGHLIN_TENSION_MEAN_KNM2 <- 180
COUGHLIN_TENSION_SD_KNM2   <- 33.6

.load_one <- function(f) {
  td <- load_bender_flat(f, do_filter = TRUE, loadtorques = "x")
  tau <- deconvolve_bender(f, hub_path = NULL, verbose = FALSE)
  N <- min(nrow(td), length(tau)); td <- td[seq_len(N), , drop = FALSE]
  td$torque_inertia_corrected_Nm <- tau[seq_len(N)]
  attr(td, "Filename") <- f
  td
}

# ==========================================================================
# Pooled position+time passive-regression model -- ONE trial's step_summary
# in, same step_summary out with muscle_force_Nm REPLACED (original moved to
# muscle_force_Nm_static_orig for QA) and pooled_fit_r2 attached.
# ==========================================================================
fit_trial_pooled_passive_model <- function(steps, deactivation_window_s = 0.5) {
  wc_posix <- suppressWarnings(as.POSIXct(steps$wall_clock_start, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC"))
  t0_session <- suppressWarnings(min(wc_posix, na.rm = TRUE))
  if (!is.finite(t0_session)) { steps$muscle_force_Nm_pooled <- NA_real_; steps$pooled_fit_r2 <- NA_real_; return(steps) }
  steps$.elapsed_s     <- as.numeric(difftime(wc_posix, t0_session, units = "secs"))
  steps$.t_pre_mid_s   <- (steps$t_pre_baseline_start_s + steps$t_pre_baseline_end_s) / 2
  steps$.t_post_mid_s  <- (steps$t_post_baseline_start_s + steps$t_post_baseline_end_s) / 2
  steps$.t_active_mid_s <- (steps$stim_t0_s + (steps$stim_t1_s + deactivation_window_s)) / 2

  train <- dplyr::bind_rows(
    tibble::tibble(elapsed_s = steps$.elapsed_s + steps$.t_pre_mid_s, operating_point = steps$operating_point,
                  torque = steps$passive_force_Nm_static),
    tibble::tibble(elapsed_s = steps$.elapsed_s + steps$.t_post_mid_s, operating_point = steps$operating_point,
                torque = steps$post_force_Nm_static)
  ) |> dplyr::filter(is.finite(.data$elapsed_s), is.finite(.data$operating_point), is.finite(.data$torque))

  if (nrow(train) < 6L || dplyr::n_distinct(train$operating_point) < 3L) {
    steps$muscle_force_Nm_pooled <- NA_real_; steps$pooled_fit_r2 <- NA_real_
    steps$.elapsed_s <- NULL; steps$.t_pre_mid_s <- NULL; steps$.t_post_mid_s <- NULL; steps$.t_active_mid_s <- NULL
    return(steps)
  }
  fit <- tryCatch(stats::lm(torque ~ elapsed_s + operating_point + I(operating_point^2), data = train),
                  error = function(e) NULL)
  if (is.null(fit)) {
    steps$muscle_force_Nm_pooled <- NA_real_; steps$pooled_fit_r2 <- NA_real_
  } else {
    passive_pooled <- as.numeric(predict(fit, newdata = data.frame(
      elapsed_s = steps$.elapsed_s + steps$.t_active_mid_s, operating_point = steps$operating_point
    )))
    steps$muscle_force_Nm_pooled <- steps$force_sign * (steps$active_force_Nm - passive_pooled)
    steps$pooled_fit_r2 <- summary(fit)$r.squared
  }
  steps$.elapsed_s <- NULL; steps$.t_pre_mid_s <- NULL; steps$.t_post_mid_s <- NULL; steps$.t_active_mid_s <- NULL
  steps
}

# ==========================================================================
# Collect isometric steps, all specimens, pooled-regression baseline per trial
# ==========================================================================
cli::cli_h1("Collecting isometric steps (pooled-regression baseline)")

.collect_isometric_pooled <- function(specimen, source_dir) {
  manifest <- parse_trial_directory(source_dir)
  files <- manifest$fullpath[manifest$protocol == "isometric"]
  purrr::map_dfr(files, function(f) {
    trial_id <- tools::file_path_sans_ext(basename(f))
    td <- tryCatch(.load_one(f), error = function(e) NULL)
    if (is.null(td)) return(tibble())
    built <- tryCatch(build_segmented_step_summary(td, f), error = function(e) NULL)
    if (is.null(built) || nrow(built$step_summary) == 0L || !is.finite(built$r_m) || built$r_m <= 0) return(tibble())
    steps <- fit_trial_pooled_passive_model(built$step_summary)
    cli::cli_alert_success("{trial_id}: pooled fit R^2={round(mean(steps$pooled_fit_r2, na.rm=TRUE),3)}, {sum(is.finite(steps$muscle_force_Nm_pooled))} steps ok")
    steps$specimen  <- specimen
    steps$trial_id  <- trial_id
    steps$trial_num <- extract_bender_trial_num(trial_id)
    steps$r_m       <- built$r_m
    steps
  })
}

iso_pooled <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  .collect_isometric_pooled(specimen, raw_source_dir(subfolder))
})
iso_pooled$precondition <- classify_session_precondition(iso_pooled$specimen, iso_pooled$trial_num)

# Swap in the pooled-regression force for the shared plot builder (documented
# pattern, build_summary_plot_isometric() docstring) -- keep the static-
# baseline original for reference/QA, never overwritten.
iso_for_plot <- iso_pooled
iso_for_plot$muscle_force_Nm_static_orig <- iso_for_plot$muscle_force_Nm
iso_for_plot$muscle_force_Nm <- iso_for_plot$muscle_force_Nm_pooled

p_fl <- build_summary_plot_isometric(
  iso_for_plot,
  title = "Isometric Force-Length -- POOLED-REGRESSION passive baseline (EXPLORATORY)"
)
p_fl <- p_fl + labs(caption = paste0(
  "EXPLORATORY: passive baseline = per-trial regression on elapsed session time + operating_point (quadratic), ",
  "not the production static-baseline method. Disagreed with simpler methods on the one step tested (see ",
  "analysis_muscle_force_vector_log.md, 2026-07-25) -- NOT a validated replacement for the canonical FL figure.\n",
  "Mean pooled-fit R^2 = ", round(mean(iso_pooled$pooled_fit_r2, na.rm = TRUE), 3),
  " across ", dplyr::n_distinct(iso_pooled$trial_id), " isometric trial(s)."
))
fout_fl <- file.path(OUT_DIR, "FL_isometric_pooledRegressionBaseline.png")
ggplot2::ggsave(fout_fl, p_fl, width = 9, height = 6.5, dpi = 150)
cli::cli_alert_success("Saved {fout_fl}")

# ==========================================================================
# Collect isovelocity steps, all specimens, pooled-regression baseline
# ==========================================================================
cli::cli_h1("Collecting isovelocity steps (pooled-regression baseline)")

.collect_isovelocity_pooled <- function(specimen, source_dir) {
  manifest <- parse_trial_directory(source_dir)
  files <- manifest$fullpath[manifest$protocol == "isovelocity"]
  purrr::map_dfr(files, function(f) {
    trial_id <- tools::file_path_sans_ext(basename(f))
    td <- tryCatch(.load_one(f), error = function(e) NULL)
    if (is.null(td)) return(tibble())
    res <- tryCatch(analyze_isovelocity(td, f), error = function(e) NULL)
    if (is.null(res) || nrow(res$step_summary) == 0L || !is.finite(res$r_m) || res$r_m <= 0) return(tibble())
    steps <- fit_trial_pooled_passive_model(res$step_summary)
    cli::cli_alert_success("{trial_id}: pooled fit R^2={round(mean(steps$pooled_fit_r2, na.rm=TRUE),3)}, {sum(is.finite(steps$muscle_force_Nm_pooled))} steps ok")
    steps$specimen  <- specimen
    steps$trial_id  <- trial_id
    steps$trial_num <- extract_bender_trial_num(trial_id)
    steps$r_m       <- res$r_m
    steps
  })
}

isovel_pooled <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  .collect_isovelocity_pooled(specimen, raw_source_dir(subfolder))
})
isovel_pooled$precondition <- classify_session_precondition(isovel_pooled$specimen, isovel_pooled$trial_num)

isovel_for_plot <- isovel_pooled
isovel_for_plot$muscle_force_Nm_static_orig <- isovel_for_plot$muscle_force_Nm
isovel_for_plot$muscle_force_Nm <- isovel_for_plot$muscle_force_Nm_pooled

# Refit the Hill curve on the pooled-regression force (mirrors analyze_
# isovelocity()'s own per-side construction) -- specific-tension conversion
# uses a REPRESENTATIVE geometry (mean r_m/dclamp_m across trials actually
# pooled here) purely for the Hill annotation's specific-tension label; the
# per-point/per-trial force values above are NOT rescaled by this.
fits_pooled <- list()
for (sd in c("left", "right")) {
  ss <- dplyr::filter(isovel_for_plot, .data$muscle_side == sd, is.finite(.data$muscle_force_Nm))
  f0_row <- dplyr::filter(ss, abs(.data$shortening_value) < 1e-6)
  f0_iso <- if (nrow(f0_row) > 0L) mean(f0_row$muscle_force_Nm, na.rm = TRUE) else NA_real_
  ss_conc <- dplyr::filter(ss, .data$contraction_mode %in% c("concentric", "isometric_zero"))
  fits_pooled[[sd]] <- fit_force_velocity_curve(ss_conc$shortening_strain_pct, ss_conc$muscle_force_Nm,
                                                side_label = sd, f0_isometric = f0_iso)
}

p_fv <- build_summary_plot_isovelocity(isovel_for_plot, fits_pooled)
# build_summary_plot_isovelocity() has no title= param (unlike its isometric
# sibling) -- override its hardcoded title here so this plot can't be
# mistaken for the canonical FV figure at a glance.
p_fv <- p_fv + labs(title = "Isovelocity Force-Velocity -- POOLED-REGRESSION passive baseline (EXPLORATORY)")
p_fv <- p_fv + labs(caption = paste0(
  "EXPLORATORY: passive baseline = per-trial regression on elapsed session time + operating_point (quadratic), ",
  "NOT isovelocity's own better-motivated velocity-matched baseline and NOT the production static-baseline method. ",
  "Disagreed with simpler methods on the one step tested (see analysis_muscle_force_vector_log.md, 2026-07-25) -- ",
  "NOT a validated replacement for the canonical FV figure.\nMean pooled-fit R^2 = ",
  round(mean(isovel_pooled$pooled_fit_r2, na.rm = TRUE), 3),
  " across ", dplyr::n_distinct(isovel_pooled$trial_id), " isovelocity trial(s)."
))
fout_fv <- file.path(OUT_DIR, "FV_isovelocity_pooledRegressionBaseline.png")
ggplot2::ggsave(fout_fv, p_fv, width = 9, height = 6.5, dpi = 150)
cli::cli_alert_success("Saved {fout_fv}")

# ==========================================================================
# Isometric tension summary (pooled-regression baseline), right side only
# (production convention), all preconditioning phases shown/colored.
# ==========================================================================
cli::cli_h1("Isometric tension summary (pooled-regression baseline)")

iso_pooled$force_N_pooled <- iso_pooled$muscle_force_Nm_pooled / iso_pooled$r_m
iso_pooled$specific_tension_kNm2_pooled <- (iso_pooled$force_N_pooled / MEASURED_RED_MUSCLE_CSA_CM2) * 10

tension_right <- dplyr::filter(iso_pooled, .data$muscle_side == "right", is.finite(.data$specific_tension_kNm2_pooled))

p_tension <- ggplot(tension_right, aes(x = abs(.data$shortening_strain_pct), y = .data$specific_tension_kNm2_pooled)) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = COUGHLIN_TENSION_MEAN_KNM2 - COUGHLIN_TENSION_SD_KNM2,
           ymax = COUGHLIN_TENSION_MEAN_KNM2 + COUGHLIN_TENSION_SD_KNM2,
           fill = "grey70", alpha = 0.3) +
  geom_hline(yintercept = COUGHLIN_TENSION_MEAN_KNM2, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
  geom_point(aes(color = .data$specimen, shape = .data$precondition), size = 2.6, alpha = 0.85) +
  labs(title = "Isometric specific tension -- POOLED-REGRESSION passive baseline (EXPLORATORY)",
       subtitle = "Right-muscle steps, all isometric trials, all preconditioning phases shown.\nGrey ribbon/dashed = Coughlin (2000), 180 +/- 33.6 kN/m^2.",
       x = "|Predicted muscle shortening strain| (%L0)", y = "Specific tension (kN/m^2)",
       caption = "EXPLORATORY -- see FL/FV caveats above; NOT a validated replacement for the canonical tension panel (summary_coughlin2000_bass_comparison.R).") +
  theme_bw(base_size = 11) + theme(legend.position = "right")
fout_tension <- file.path(OUT_DIR, "isometricTension_pooledRegressionBaseline.png")
ggplot2::ggsave(fout_tension, p_tension, width = 9, height = 6, dpi = 150)
cli::cli_alert_success("Saved {fout_tension}")

cli::cli_alert_success("summary_fl_fv_tension_pooled_regression_baseline.R complete -- outputs in {OUT_DIR}/")
