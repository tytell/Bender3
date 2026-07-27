# summary_fl_fv_tension_optimizedbaseline.R
# PI directive (2026-07-27): "Let's abandon legacy and move forward with
# best practice. Don't delete figs yet. Just append 'optimizedbaseline' [to]
# the summary and diagnostic figures."
#
# "Best practice" = `passive_force_Nm_optimizedbaseline` /
# `muscle_force_Nm_optimizedbaseline` (build_segmented_step_summary(),
# 03_analyze.R, added 2026-07-27): prefer a REAL matched no-stim measurement
# (velocity_matched, isovelocity only) when one exists; otherwise use
# baselineInterp (pre+post linear interpolation) rather than the old
# pre-only static baseline ("legacy"). See analysis_muscle_force_vector_
# log.md 2026-07-27 addenda for the full justification (legacy's
# direction-dependent bias: overestimates concentric |force|, underestimates
# eccentric |force|; baselineInterp/optimizedbaseline removed that bias in
# all 3 specimens without a validated better alternative -- the v2
# anchoredRate_dirTime model CONVERGED onto baselineInterp rather than
# beating it).
#
# UNLIKE summary_fl_fv_tension_pooled_regression_baseline.R (still kept,
# NOT deleted, still labeled EXPLORATORY -- that pooled-position+time
# regression method was never validated as more accurate), these figures
# use the NOW-ENDORSED default method. Old legacy-baseline canonical
# figures (run_fv_fl_power_pipeline.R / plot_summary_profiles.R /
# summary_coughlin2000_bass_comparison.R) are LEFT IN PLACE, not deleted or
# overwritten -- per PI directive, this script only ADDS new
# "_optimizedbaseline"-suffixed siblings alongside them.
#
# Reuses the EXISTING, unmodified plot builders (build_summary_plot_
# isometric()/_isovelocity(), plot_summary_profiles.R) via the documented
# "pass a muscle_force_Nm-substituted step_summary plus a distinguishing
# title" pattern -- same pattern as the existing "_baselineInterp" and
# "_pooledRegressionBaseline" variants.
#
# Run with:  Rscript R/summary_fl_fv_tension_optimizedbaseline.R
# Outputs -> 02_processed/figs_summary/:
#   FL_isometric_optimizedbaseline.png
#   FV_isovelocity_optimizedbaseline.png
#   isometricTension_optimizedbaseline.png

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
src("03_analyze.R")            # defines muscle_force_Nm_optimizedbaseline directly, no per-trial refit needed
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
# Collect isometric steps, all specimens
# ==========================================================================
cli::cli_h1("Collecting isometric steps (optimizedbaseline)")

.collect_isometric <- function(specimen, source_dir) {
  manifest <- parse_trial_directory(source_dir)
  files <- manifest$fullpath[manifest$protocol == "isometric"]
  purrr::map_dfr(files, function(f) {
    trial_id <- tools::file_path_sans_ext(basename(f))
    td <- tryCatch(.load_one(f), error = function(e) NULL)
    if (is.null(td)) return(tibble())
    built <- tryCatch(build_segmented_step_summary(td, f), error = function(e) NULL)
    if (is.null(built) || nrow(built$step_summary) == 0L || !is.finite(built$r_m) || built$r_m <= 0) return(tibble())
    steps <- built$step_summary
    cli::cli_alert_success("{trial_id}: {sum(is.finite(steps$muscle_force_Nm_optimizedbaseline))} steps ok")
    steps$specimen  <- specimen
    steps$trial_id  <- trial_id
    steps$trial_num <- extract_bender_trial_num(trial_id)
    steps$r_m       <- built$r_m
    steps
  })
}

iso_opt <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  .collect_isometric(specimen, raw_source_dir(subfolder))
})
iso_opt$precondition <- classify_session_precondition(iso_opt$specimen, iso_opt$trial_num)

# Swap in the optimizedbaseline force for the shared plot builder (documented
# pattern, build_summary_plot_isometric() docstring) -- keep the legacy
# static-baseline original for reference/QA, never overwritten.
iso_for_plot <- iso_opt
iso_for_plot$muscle_force_Nm_legacy_orig <- iso_for_plot$muscle_force_Nm
iso_for_plot$muscle_force_Nm <- iso_for_plot$muscle_force_Nm_optimizedbaseline

p_fl <- build_summary_plot_isometric(
  iso_for_plot,
  title = "Isometric Force-Length -- optimizedbaseline (baselineInterp; legacy ABANDONED 2026-07-27)"
)
p_fl <- p_fl + labs(caption = paste0(
  "optimizedbaseline = velocity_matched where a real no-stim match exists, else baselineInterp (pre+post linear ",
  "interpolation) -- legacy (pre-only static baseline) abandoned 2026-07-27 for a confirmed direction-dependent ",
  "bias (see analysis_muscle_force_vector_log.md). Legacy-baseline canonical figure NOT deleted; kept for reference."
))
fout_fl <- file.path(OUT_DIR, "FL_isometric_optimizedbaseline.png")
ggplot2::ggsave(fout_fl, p_fl, width = 9, height = 6.5, dpi = 150)
cli::cli_alert_success("Saved {fout_fl}")

# ==========================================================================
# Collect isovelocity steps, all specimens
# ==========================================================================
cli::cli_h1("Collecting isovelocity steps (optimizedbaseline)")

.collect_isovelocity <- function(specimen, source_dir) {
  manifest <- parse_trial_directory(source_dir)
  files <- manifest$fullpath[manifest$protocol == "isovelocity"]
  purrr::map_dfr(files, function(f) {
    trial_id <- tools::file_path_sans_ext(basename(f))
    td <- tryCatch(.load_one(f), error = function(e) NULL)
    if (is.null(td)) return(tibble())
    res <- tryCatch(analyze_isovelocity(td, f), error = function(e) NULL)
    if (is.null(res) || nrow(res$step_summary) == 0L || !is.finite(res$r_m) || res$r_m <= 0) return(tibble())
    steps <- res$step_summary
    cli::cli_alert_success("{trial_id}: {sum(is.finite(steps$muscle_force_Nm_optimizedbaseline))} steps ok")
    steps$specimen  <- specimen
    steps$trial_id  <- trial_id
    steps$trial_num <- extract_bender_trial_num(trial_id)
    steps$r_m       <- res$r_m
    steps
  })
}

isovel_opt <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  .collect_isovelocity(specimen, raw_source_dir(subfolder))
})
isovel_opt$precondition <- classify_session_precondition(isovel_opt$specimen, isovel_opt$trial_num)

isovel_for_plot <- isovel_opt
isovel_for_plot$muscle_force_Nm_legacy_orig <- isovel_for_plot$muscle_force_Nm
isovel_for_plot$muscle_force_Nm <- isovel_for_plot$muscle_force_Nm_optimizedbaseline

# Refit the Hill curve on the optimizedbaseline force (mirrors analyze_
# isovelocity()'s own per-side construction).
fits_opt <- list()
for (sd in c("left", "right")) {
  ss <- dplyr::filter(isovel_for_plot, .data$muscle_side == sd, is.finite(.data$muscle_force_Nm))
  f0_row <- dplyr::filter(ss, abs(.data$shortening_value) < 1e-6)
  f0_iso <- if (nrow(f0_row) > 0L) mean(f0_row$muscle_force_Nm, na.rm = TRUE) else NA_real_
  ss_conc <- dplyr::filter(ss, .data$contraction_mode %in% c("concentric", "isometric_zero"))
  fits_opt[[sd]] <- fit_force_velocity_curve(ss_conc$shortening_strain_pct, ss_conc$muscle_force_Nm,
                                              side_label = sd, f0_isometric = f0_iso)
}

p_fv <- build_summary_plot_isovelocity(isovel_for_plot, fits_opt)
# build_summary_plot_isovelocity() has no title= param (unlike its isometric
# sibling) -- override its hardcoded title here so this plot can't be
# mistaken for the legacy-baseline canonical FV figure at a glance.
p_fv <- p_fv + labs(title = "Isovelocity Force-Velocity -- optimizedbaseline (velocity_matched + baselineInterp fallback)")
p_fv <- p_fv + labs(caption = paste0(
  "optimizedbaseline = velocity_matched (unchanged from production where a real no-stim match exists), else ",
  "baselineInterp instead of the old static-baseline fallback -- legacy fallback abandoned 2026-07-27 (see ",
  "analysis_muscle_force_vector_log.md). Legacy-baseline canonical figure NOT deleted; kept for reference."
))
fout_fv <- file.path(OUT_DIR, "FV_isovelocity_optimizedbaseline.png")
ggplot2::ggsave(fout_fv, p_fv, width = 9, height = 6.5, dpi = 150)
cli::cli_alert_success("Saved {fout_fv}")

# ==========================================================================
# Isometric tension summary (optimizedbaseline), right side only
# (production convention), all preconditioning phases shown/colored.
# ==========================================================================
cli::cli_h1("Isometric tension summary (optimizedbaseline)")

iso_opt$force_N_optimizedbaseline <- iso_opt$muscle_force_Nm_optimizedbaseline / iso_opt$r_m
iso_opt$specific_tension_kNm2_optimizedbaseline <- (iso_opt$force_N_optimizedbaseline / MEASURED_RED_MUSCLE_CSA_CM2) * 10

tension_right <- dplyr::filter(iso_opt, .data$muscle_side == "right", is.finite(.data$specific_tension_kNm2_optimizedbaseline))

p_tension <- ggplot(tension_right, aes(x = abs(.data$shortening_strain_pct), y = .data$specific_tension_kNm2_optimizedbaseline)) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = COUGHLIN_TENSION_MEAN_KNM2 - COUGHLIN_TENSION_SD_KNM2,
           ymax = COUGHLIN_TENSION_MEAN_KNM2 + COUGHLIN_TENSION_SD_KNM2,
           fill = "grey70", alpha = 0.3) +
  geom_hline(yintercept = COUGHLIN_TENSION_MEAN_KNM2, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
  geom_point(aes(color = .data$specimen, shape = .data$precondition), size = 2.6, alpha = 0.85) +
  labs(title = "Isometric specific tension -- optimizedbaseline (legacy ABANDONED 2026-07-27)",
       subtitle = "Right-muscle steps, all isometric trials, all preconditioning phases shown.\nGrey ribbon/dashed = Coughlin (2000), 180 +/- 33.6 kN/m^2.",
       x = "|Predicted muscle shortening strain| (%L0)", y = "Specific tension (kN/m^2)",
       caption = "optimizedbaseline = baselineInterp (isometric has no velocity_matched option). Legacy-baseline canonical\ntension panel (summary_coughlin2000_bass_comparison.R) NOT deleted; kept for reference.") +
  theme_bw(base_size = 11) + theme(legend.position = "right")
fout_tension <- file.path(OUT_DIR, "isometricTension_optimizedbaseline.png")
ggplot2::ggsave(fout_tension, p_tension, width = 9, height = 6, dpi = 150)
cli::cli_alert_success("Saved {fout_tension}")

cli::cli_alert_success("summary_fl_fv_tension_optimizedbaseline.R complete -- outputs in {OUT_DIR}/")
