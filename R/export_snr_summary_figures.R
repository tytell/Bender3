# export_snr_summary_figures.R
# Bulk export for ONE specimen: copies/regenerates isometric + isovelocity
# figures into the shared figs_summary/ folder, prefixed with the
# specimen ID, restricted to activation-SNR-passing data (PI-directed,
# 2026-07-18). Scope is isometric/isovelocity ONLY -- dynamic/frequency_sweep
# and pooled validation plots (angleValid, strainValidCmd,
# strainValidSono{Enc,Cmd}) have no per-step activation_snr and are out of
# scope for this export.
#
# Three handling rules, by plot type:
#   1. TRIAL plots (per-fish folder is flat, 2026-07-21 -- matched by
#      filename shape "{bassID}_bender_{NN}_{isometric|isovelocity}*.png",
#      not by subfolder): copied AS-IS (single-trial view, unfiltered) only
#      if a MAJORITY of that trial's steps have activation_snr >= MFV_UHAT_SNR_MIN.
#   2. SNR-AWARE summary plots (u_hat vector FL/FV, FV L0 sono, u_hat
#      empirical-vs-geometric, force-vs-time vector): REGENERATED from the
#      SNR-passing subset only (steps/points below threshold are DROPPED,
#      not just alpha-flagged), via the same builder functions
#      run_fv_fl_power_pipeline.R uses -- suffixed "_snrPass" (see
#      FIGURES_README.md for the naming convention).
#   3. LEGACY (pre-vector, no activation_snr field) summary plots
#      (FL_isometric_legacy.png, FV_isovelocity_legacy.png,
#      fatigueCheck_*_legacy.png, strainValidMeas.png,
#      forceTime_isometric_legacy.png, forceTime_isovelocity_legacy.png):
#      copied AS-IS -- there is no SNR to filter by on the legacy
#      z-torque-only path.
#
# Run once per specimen (mirrors run_fv_fl_power_pipeline.R's env-var
# pattern; sourcing that script regenerates the normal per-fish figures too):
#   BENDER3_BASS_ID=bass17 \
#   BENDER3_SOURCE_DIR=/path/to/2026-07-15_bass17_bender \
#   BENDER3_OUTPUT_DIR=/path/to/figs_bass17 \
#   Rscript R/export_snr_summary_figures.R

suppressPackageStartupMessages({library(dplyr); library(fs); library(cli); library(ggplot2)})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
source(file.path(.pipeline_root, "paths_config.R"))

BASS_ID       <- Sys.getenv("BENDER3_BASS_ID", "bass16")
# Default comes from paths_config.R (single source of truth) -- see that
# file if the OneDrive folder layout ever moves again.
SUMMARY_DEST  <- FIGS_SUMMARY_DIR
fs::dir_create(SUMMARY_DEST, recurse = TRUE)

source(file.path(.pipeline_root, "run_fv_fl_power_pipeline.R"))  # regenerates figures; populates iso_vec/isv_vec/etc.

cli::cli_h1("Exporting SNR-passing figures for {BASS_ID} -> figs_summary/")

snr_filter <- function(df) {
  if (is.null(df) || nrow(df) == 0L) return(df)
  dplyr::filter(df, is.finite(.data$activation_snr), .data$activation_snr >= MFV_UHAT_SNR_MIN)
}

# ---- 1. trial_plots: copy AS-IS only if a MAJORITY of the trial's steps pass ----
.trial_pass_fraction <- function(trial_id) {
  ss <- if (trial_id %in% names(isometric_steps_vec_all))       isometric_steps_vec_all[[trial_id]]
        else if (trial_id %in% names(isovelocity_steps_vec_all)) isovelocity_steps_vec_all[[trial_id]]
        else NULL
  if (is.null(ss)) return(NA_real_)
  ss <- dplyr::filter(ss, .data$muscle_side %in% c("left", "right"), is.finite(.data$activation_snr))
  if (nrow(ss) == 0L) return(NA_real_)
  mean(ss$activation_snr >= MFV_UHAT_SNR_MIN)
}

# Per-fish folders are now FLAT (2026-07-21, no trial_plots/ subfolder), so a
# bare "*isometric*.png"/"*isovelocity*.png" glob would also catch summary
# plots that happen to contain those substrings (FL_isometric_legacy.png,
# forceTime_isovelocity_uhatBoth.png, etc.) -- folder scoping no longer
# disambiguates trial vs. summary, so the glob must anchor on the trial-plot
# filename SHAPE instead: "{BASS_ID}_bender_{NN}_{protocol}[...].png".
trial_pngs <- unique(c(
  fs::dir_ls(TRIAL_PLOT_DIR, glob = sprintf("%s_bender_*isometric*.png", BASS_ID)),
  fs::dir_ls(TRIAL_PLOT_DIR, glob = sprintf("%s_bender_*isovelocity*.png", BASS_ID))
))
n_trial_copied <- 0L
for (f in trial_pngs) {
  fname    <- basename(f)
  trial_id <- sub("_baselineInterp$", "", sub("\\.png$", "", fname))
  frac     <- .trial_pass_fraction(trial_id)
  if (is.na(frac)) {
    cli::cli_alert_warning("{fname}: no SNR data for trial_id '{trial_id}' -- skipped")
    next
  }
  if (frac > 0.5) {
    # trial_plots filenames already start with the specimen ID (e.g.
    # bass16_bender_09_isometric.png) -- don't double it.
    dest_name <- if (startsWith(fname, paste0(BASS_ID, "_"))) fname else paste0(BASS_ID, "_", fname)
    fs::file_copy(f, file.path(SUMMARY_DEST, dest_name), overwrite = TRUE)
    n_trial_copied <- n_trial_copied + 1L
    cli::cli_alert_success("{fname}: {round(100*frac,1)}% of steps pass -> copied")
  } else {
    cli::cli_alert_info("{fname}: only {round(100*frac,1)}% of steps pass -- omitted")
  }
}
cli::cli_alert_info("Trial plots copied: {n_trial_copied} / {length(trial_pngs)}")

# ---- 2. legacy (no activation_snr field): copy AS-IS ----
legacy_summary <- c(
  "FL_isometric_legacy.png", "FL_isometric_baselineInterp.png",
  "FV_isovelocity_legacy.png",
  "fatigueCheck_isometric_legacy.png", "fatigueCheck_isometric_baselineInterp.png",
  "fatigueCheck_isovelocity_legacy.png",
  "strainValidMeas.png",
  "forceTime_isometric_legacy.png", "forceTime_isovelocity_legacy.png"
)
n_legacy <- 0L
for (fname in legacy_summary) {
  src <- file.path(SUMMARY_PLOT_DIR, fname)
  if (!file.exists(src)) next
  fs::file_copy(src, file.path(SUMMARY_DEST, paste0(BASS_ID, "_", fname)), overwrite = TRUE)
  n_legacy <- n_legacy + 1L
}
cli::cli_alert_info("Legacy (non-SNR) summary plots copied as-is: {n_legacy} (no activation_snr field on this path)")

# ---- 3. SNR-aware summary plots: regenerate from the SNR-passing subset ----
if (exists("iso_vec") && nrow(iso_vec) > 0L) {
  iso_vec_pass <- snr_filter(iso_vec)
  p <- if (!is.null(iso_vec_pass) && nrow(iso_vec_pass) > 0L) build_summary_plot_FL_vector(iso_vec_pass) else NULL
  if (!is.null(p)) {
    outp <- file.path(SUMMARY_DEST, paste0(BASS_ID, "_FL_isometric_uhatBoth_snrPass.png"))
    ggplot2::ggsave(outp, p, width = 12, height = 6, dpi = 150)
    cli::cli_alert_success("{basename(outp)}: {nrow(iso_vec_pass)}/{nrow(iso_vec)} points pass")
  } else cli::cli_alert_warning("No SNR-passing isometric vector points for {BASS_ID}")
}

if (exists("isv_vec") && nrow(isv_vec) > 0L) {
  isv_vec_pass <- snr_filter(isv_vec)
  p <- if (!is.null(isv_vec_pass) && nrow(isv_vec_pass) > 0L) build_summary_plot_FV_vector(isv_vec_pass) else NULL
  if (!is.null(p)) {
    outp <- file.path(SUMMARY_DEST, paste0(BASS_ID, "_FV_isovelocity_uhatBoth_snrPass.png"))
    ggplot2::ggsave(outp, p, width = 12, height = 6, dpi = 150)
    cli::cli_alert_success("{basename(outp)}: {nrow(isv_vec_pass)}/{nrow(isv_vec)} points pass")
  } else cli::cli_alert_warning("No SNR-passing isovelocity vector points for {BASS_ID}")
}

if (exists("fv_l0_all") && nrow(fv_l0_all) > 0L) {
  fv_l0_pass <- snr_filter(fv_l0_all)
  p <- if (!is.null(fv_l0_pass) && nrow(fv_l0_pass) > 0L) build_summary_plot_FV_L0_vector(fv_l0_pass) else NULL
  if (!is.null(p)) {
    outp <- file.path(SUMMARY_DEST, paste0(BASS_ID, "_FVl0_isovelocity_uhatBoth_snrPass.png"))
    ggplot2::ggsave(outp, p, width = 9, height = 6, dpi = 150)
    cli::cli_alert_success("{basename(outp)}: {nrow(fv_l0_pass)}/{nrow(fv_l0_all)} points pass")
  } else cli::cli_alert_warning("No SNR-passing FV L0 (sono) points for {BASS_ID}")
}

if (exists("uhat_tbl_all") && length(uhat_tbl_all) > 0L) {
  uhat_all <- dplyr::bind_rows(uhat_tbl_all)
  uhat_pass <- snr_filter(uhat_all)
  p <- if (!is.null(uhat_pass) && nrow(uhat_pass) > 0L) build_uhat_comparison_plot(uhat_pass) else NULL
  if (!is.null(p)) {
    outp <- file.path(SUMMARY_DEST, paste0(BASS_ID, "_uhatCompare_empVsGeom_snrPass.png"))
    ggplot2::ggsave(outp, p, width = 8, height = 7, dpi = 150)
    cli::cli_alert_success("{basename(outp)}: {nrow(uhat_pass)}/{nrow(uhat_all)} steps pass")
  } else cli::cli_alert_warning("No SNR-passing u_hat comparison steps for {BASS_ID}")
}

for (fam in names(force_ts_vec_all)) {
  ts_list <- force_ts_vec_all[[fam]]
  if (length(ts_list) == 0L) next
  ts_df <- dplyr::bind_rows(ts_list)
  ts_df_pass <- snr_filter(ts_df)
  if (is.null(ts_df_pass) || nrow(ts_df_pass) == 0L) {
    cli::cli_alert_warning("No SNR-passing {fam} force-vs-time rows for {BASS_ID}")
    next
  }
  facet_var <- if (fam == "isovelocity") "contraction_mode" else NULL
  p <- build_force_vs_time_vector_plot(
    ts_df_pass, title = sprintf("Muscle force along u_hat vs time (%s, SNR-passing only)", fam),
    facet_var = facet_var)
  if (!is.null(p)) {
    outp <- file.path(SUMMARY_DEST, paste0(BASS_ID, sprintf("_forceTime_%s_uhatBoth_snrPass.png", fam)))
    ggplot2::ggsave(outp, p, width = if (!is.null(facet_var)) 12 else 10, height = 6, dpi = 150)
    cli::cli_alert_success("{basename(outp)}: {nrow(ts_df_pass)}/{nrow(ts_df)} rows pass")
  }
}

cli::cli_alert_success("Done: {BASS_ID} SNR-passing export complete -> {SUMMARY_DEST}")
