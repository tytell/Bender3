# export_snr_summary_figures.R
# Regenerates SNR-filtered ("_snrPass") figure variants for ONE specimen,
# saved INTO THAT SPECIMEN'S OWN FOLDER (figs_{bassID}/) -- restricted to
# activation-SNR-passing data (PI-directed, 2026-07-18). Scope is
# isometric/isovelocity ONLY -- dynamic/frequency_sweep and the rig/system
# validation plots (angleValid, strainValidCmd, strainValidSono{Enc,Cmd} --
# relocated to figs_diagnostic/ 2026-07-21) have no per-step activation_snr
# and are out of scope for this export.
#
# NAME IS NOW A PARTIAL MISNOMER (2026-07-21): this script used to ALSO
# copy trial plots and legacy (non-SNR) summary plots AS-IS into the shared
# figs_summary/ folder, prefixed with the specimen ID. PI feedback
# (2026-07-21): those copies were pure DUPLICATES of files that already
# exist, unprefixed, in figs_{bassID}/ -- "completely unnecessary and
# unhelpful ... they should appear in the bass## folders" (where they
# already did). Both copy mechanisms were REMOVED, not relocated -- nothing
# was lost, the originals were never touched. This script now does exactly
# ONE thing: regenerate the SNR-passing-only ("_snrPass") variant of each
# vector-force summary plot and save it ALONGSIDE the unfiltered version in
# the SAME per-fish folder (figs_{bassID}/) -- it no longer writes to
# figs_summary/ at all. A rename (e.g. to `regenerate_snr_pass_figures.R`)
# would now more accurately describe it -- not done here since it wasn't
# requested and touches call-site references; flag if you want that too.
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

BASS_ID <- Sys.getenv("BENDER3_BASS_ID", "bass16")

source(file.path(.pipeline_root, "run_fv_fl_power_pipeline.R"))  # regenerates figures; populates iso_vec/isv_vec/etc.
# SUMMARY_PLOT_DIR (set by run_fv_fl_power_pipeline.R) IS this specimen's
# own folder (figs_{bassID}/, per the 2026-07-21 flatten) -- every ggsave
# below writes there, NOT to figs_summary/.

cli::cli_h1("Regenerating SNR/magnitude-passing figures for {BASS_ID} -> {SUMMARY_PLOT_DIR}")

# REVISED 2026-07-22 (PI-directed "SNR-based confidence gating audit" --
# implements analysis_muscle_force_vector_log.md's 2026-07-22 proposal). Two
# filters now, not one:
#   snr_filter()  -- UNCHANGED, ratio-only. Kept for uhat_tbl_all specifically:
#                    that table compares DIRECTIONS (angles), which have no
#                    magnitude to test against, and u_hat source selection is
#                    intentionally ratio-only by design (a low-SNR point's
#                    DIRECTION really is less reliable regardless of whether
#                    its force magnitude is real -- see
#                    muscle_force_vector.R's .mfv_finalize_step() use_emp
#                    logic / the audit's item 1 finding).
#   tier_filter()  -- NEW. For every OTHER "_snrPass" family (force/time-
#                    series data, which DOES have a magnitude to test): keeps
#                    a row if AT LEAST ONE of the supplied force columns'
#                    confidence tier (mfv_confidence_tier(), muscle_force_
#                    vector.R) is "confident" OR "confidently_small" --
#                    i.e. no longer drops a genuinely-small-but-real point
#                    (SNR < 3 but |force| >= its own noise floor) identically
#                    to pure noise, which the pure ratio test used to do.
snr_filter <- function(df) {
  if (is.null(df) || nrow(df) == 0L) return(df)
  dplyr::filter(df, is.finite(.data$activation_snr), .data$activation_snr >= MFV_UHAT_SNR_MIN)
}

tier_filter <- function(df, force_cols) {
  if (is.null(df) || nrow(df) == 0L) return(df)
  keep <- rep(FALSE, nrow(df))
  for (fc in force_cols) {
    if (!fc %in% names(df)) next
    tier <- mfv_confidence_tier(df[[fc]], df$activation_snr, df$baseline_force_noise_N)
    keep <- keep | tier %in% c("confident", "confidently_small")
  }
  dplyr::filter(df, keep)
}

# ---- SNR/magnitude-aware summary plots: regenerate from the passing subset ----
# No BASS_ID prefix needed (2026-07-21, was prefixed when these went into
# the shared figs_summary/ -- now saved into this specimen's OWN folder,
# where the "_snrPass" suffix alone already distinguishes them from the
# unfiltered/alpha-flagged uhatBoth versions, e.g. FL_isometric_uhatBoth.png
# vs. FL_isometric_uhatBoth_snrPass.png -- no collision).
if (exists("iso_vec") && nrow(iso_vec) > 0L) {
  iso_vec_pass <- tier_filter(iso_vec, c("muscle_force_vector_N", "muscle_force_vector_geom_N"))
  p <- if (!is.null(iso_vec_pass) && nrow(iso_vec_pass) > 0L) build_summary_plot_FL_vector(iso_vec_pass) else NULL
  if (!is.null(p)) {
    outp <- file.path(SUMMARY_PLOT_DIR, "FL_isometric_uhatBoth_snrPass.png")
    ggplot2::ggsave(outp, p, width = 12, height = 6, dpi = 150)
    cli::cli_alert_success("{basename(outp)}: {nrow(iso_vec_pass)}/{nrow(iso_vec)} points pass")
  } else cli::cli_alert_warning("No confidence-tier-passing isometric vector points for {BASS_ID}")
}

if (exists("isv_vec") && nrow(isv_vec) > 0L) {
  isv_vec_pass <- tier_filter(isv_vec, c("muscle_force_vector_N", "muscle_force_vector_geom_N"))
  p <- if (!is.null(isv_vec_pass) && nrow(isv_vec_pass) > 0L) build_summary_plot_FV_vector(isv_vec_pass) else NULL
  if (!is.null(p)) {
    outp <- file.path(SUMMARY_PLOT_DIR, "FV_isovelocity_uhatBoth_snrPass.png")
    ggplot2::ggsave(outp, p, width = 12, height = 6, dpi = 150)
    cli::cli_alert_success("{basename(outp)}: {nrow(isv_vec_pass)}/{nrow(isv_vec)} points pass")
  } else cli::cli_alert_warning("No confidence-tier-passing isovelocity vector points for {BASS_ID}")
}

if (exists("fv_l0_all") && nrow(fv_l0_all) > 0L) {
  fv_l0_pass <- tier_filter(fv_l0_all, "force_at_L0_N")
  p <- if (!is.null(fv_l0_pass) && nrow(fv_l0_pass) > 0L) build_summary_plot_FV_L0_vector(fv_l0_pass) else NULL
  if (!is.null(p)) {
    outp <- file.path(SUMMARY_PLOT_DIR, "FVl0_isovelocity_uhatBoth_snrPass.png")
    ggplot2::ggsave(outp, p, width = 9, height = 6, dpi = 150)
    cli::cli_alert_success("{basename(outp)}: {nrow(fv_l0_pass)}/{nrow(fv_l0_all)} points pass")
  } else cli::cli_alert_warning("No confidence-tier-passing FV L0 (sono) points for {BASS_ID}")
}

if (exists("uhat_tbl_all") && length(uhat_tbl_all) > 0L) {
  uhat_all <- dplyr::bind_rows(uhat_tbl_all)
  # ratio-only by design -- see snr_filter()'s header comment (u_hat DIRECTION
  # confidence has no magnitude counterpart to test against).
  uhat_pass <- snr_filter(uhat_all)
  p <- if (!is.null(uhat_pass) && nrow(uhat_pass) > 0L) build_uhat_comparison_plot(uhat_pass) else NULL
  if (!is.null(p)) {
    outp <- file.path(SUMMARY_PLOT_DIR, "uhatCompare_empVsGeom_snrPass.png")
    ggplot2::ggsave(outp, p, width = 8, height = 7, dpi = 150)
    cli::cli_alert_success("{basename(outp)}: {nrow(uhat_pass)}/{nrow(uhat_all)} steps pass")
  } else cli::cli_alert_warning("No SNR-passing u_hat comparison steps for {BASS_ID}")
}

for (fam in names(force_ts_vec_all)) {
  ts_list <- force_ts_vec_all[[fam]]
  if (length(ts_list) == 0L) next
  ts_df <- dplyr::bind_rows(ts_list)
  ts_df_pass <- tier_filter(ts_df, "muscle_force_vector_N")
  if (is.null(ts_df_pass) || nrow(ts_df_pass) == 0L) {
    cli::cli_alert_warning("No confidence-tier-passing {fam} force-vs-time rows for {BASS_ID}")
    next
  }
  facet_var <- if (fam == "isovelocity") "contraction_mode" else NULL
  p <- build_force_vs_time_vector_plot(
    ts_df_pass, title = sprintf("Muscle force along u_hat vs time (%s, confident + confidently-small only)", fam),
    facet_var = facet_var)
  if (!is.null(p)) {
    outp <- file.path(SUMMARY_PLOT_DIR, sprintf("forceTime_%s_uhatBoth_snrPass.png", fam))
    ggplot2::ggsave(outp, p, width = if (!is.null(facet_var)) 12 else 10, height = 6, dpi = 150)
    cli::cli_alert_success("{basename(outp)}: {nrow(ts_df_pass)}/{nrow(ts_df)} rows pass")
  }
}

cli::cli_alert_success("Done: {BASS_ID} SNR-passing figures regenerated in {SUMMARY_PLOT_DIR}")
