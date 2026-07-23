# superplot_fl_pooled_snr_passing.R
# SNR/magnitude-filtered companion to superplot_fl_pooled.R: keeps ONLY the
# points/steps whose confidence tier is "confident" or "confidently_small"
# ("GOOD" -- see mfv_confidence_tier(), muscle_force_vector.R), omitting only
# "unstable_magnitude"/"unconfirmable" ("BAD" -- indistinguishable from
# baseline force noise on at least one of the two tests).
#
# REVISED 2026-07-22 (PI-directed "SNR-based confidence gating audit" --
# implements the proposal from analysis_muscle_force_vector_log.md's
# 2026-07-22 entry): this file used to drop everything below
# `activation_snr >= MFV_UHAT_SNR_MIN`, a RATIO-only test that cannot
# distinguish "the noise floor is elevated" from "the true force is
# genuinely small" -- a weak-but-real contraction (e.g. a fatigued muscle, a
# brief V=0 hold) got dropped identically to pure noise. Now ALSO keeps a
# point whose force independently clears its OWN baseline_force_noise_N even
# if its SNR ratio is below threshold ("confidently_small"). A row is kept if
# EITHER its empirical OR geometric force estimate qualifies (mirrors the
# existing single-filter design, which already shares one activation_snr
# across both method facets of the same row) -- see `.tier` below.
#
# activation_snr (and now the magnitude test) is computed at STEP
# granularity (one value per isometric step, or per isovelocity active ramp)
# and broadcast to every point that step/ramp contributes -- so filtering
# drops entire steps/ramps, not individual samples within a kept step/ramp.
#
# Reuses ALL extraction logic from superplot_fl_pooled.R by sourcing it (this
# also regenerates the existing UNFILTERED diagnostic figures, unchanged) --
# then re-bins/re-plots only the passing subset of the same `pooled`
# tibble into figs_summary/, via the shared .build_fl_superplot() helper
# (2026-07-21 refactor) -- RAW and NORMALIZED (F/F0, see that file's header)
# companions, same as the unfiltered diagnostic pair.
#
# Run with:  Rscript R/superplot_fl_pooled_snr_passing.R

.root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
source(file.path(.root, "superplot_fl_pooled.R"))  # builds `pooled` (now incl. activation_snr); also sources paths_config.R

# Default comes from paths_config.R (single source of truth) -- see that
# file if the OneDrive folder layout ever moves again.
SUMMARY_OUTPUT_DIR <- FIGS_SUMMARY_DIR
fs::dir_create(SUMMARY_OUTPUT_DIR, recurse = TRUE)

n_total <- nrow(pooled)
KEEP_TIERS <- c("confident", "confidently_small")
.tier_emp  <- mfv_confidence_tier(pooled$force_emp_N,  pooled$activation_snr, pooled$baseline_force_noise_N, SNR_MIN)
.tier_geom <- mfv_confidence_tier(pooled$force_geom_N, pooled$activation_snr, pooled$baseline_force_noise_N, SNR_MIN)
passing <- dplyr::filter(pooled, .tier_emp %in% KEEP_TIERS | .tier_geom %in% KEEP_TIERS)
n_pass  <- nrow(passing)
n_confidently_small <- sum((.tier_emp == "confidently_small" | .tier_geom == "confidently_small") &
                          (.tier_emp %in% KEEP_TIERS | .tier_geom %in% KEEP_TIERS))
cli::cli_h1("SNR/magnitude-passing pooled FL superplot (confidence tier %in% {paste(KEEP_TIERS, collapse=', ')})")
cli::cli_alert_info("Kept {n_pass} / {n_total} points ({round(100*n_pass/n_total,1)}%) across {dplyr::n_distinct(passing$fish)} fish, {dplyr::n_distinct(passing$trial_id)} trials ({n_confidently_small} of those kept via the magnitude floor -- SNR ratio alone would have dropped them, see analysis_muscle_force_vector_log.md 2026-07-22)")

if (n_pass == 0L) cli::cli_abort("No confidence-tier-passing points -- nothing to plot.")

n_iso_p <- dplyr::n_distinct(dplyr::filter(passing, protocol == "isometric")$trial_id)
n_isv_p <- dplyr::n_distinct(dplyr::filter(passing, protocol == "isovelocity")$trial_id)
n_dyn_p <- dplyr::n_distinct(dplyr::filter(passing, protocol == "dynamic")$trial_id)

# Reuses .build_fl_superplot() from superplot_fl_pooled.R (2026-07-21 refactor)
# -- see that file for the RAW-vs-NORMALIZED (F/F0) rationale, and for the
# V=0-only rule (isometric + isovelocity's own V=0 holds + dynamic L0
# bookends -- NOT the isovelocity sweep, removed 2026-07-21).
p2 <- .build_fl_superplot(
  passing, "force_emp_N", "force_geom_N",
  title = "Pooled Force-Length superplot -- CONFIDENT + CONFIDENTLY-SMALL ONLY (SNR x magnitude tier)",
  subtitle = sprintf(
    "Same V=0-only pooling as the unfiltered diagnostic superplot, but %s (%d/%d, %.1f%%) omitted (confidence tier unstable_magnitude/unconfirmable -- fails BOTH activation_snr>=%g AND |force|>=own noise floor); %d points kept ONLY via the magnitude floor (SNR ratio alone would have dropped them, e.g. a weak/fatigued-but-real contraction)\n%d isometric + %d isovelocity + %d dynamic trial(s) contribute at least one passing point | %g%% length bins | connect-the-dots, NO fit | RAW absolute force",
    "low-confidence points", n_total - n_pass, n_total, 100 * (n_total - n_pass) / n_total, SNR_MIN, n_confidently_small,
    n_iso_p, n_isv_p, n_dyn_p, STRAIN_BIN_PCT),
  y_lab = "Muscle force along u_hat (N, + = reinforces passive direction)")

n_norm_avail_p <- dplyr::n_distinct(dplyr::filter(passing, is.finite(.data$force_emp_norm))$trial_id)
p2_norm <- .build_fl_superplot(
  passing, "force_emp_norm", "force_geom_norm",
  title = "Pooled Force-Length superplot, NORMALIZED -- CONFIDENT + CONFIDENTLY-SMALL ONLY",
  subtitle = sprintf(
    "Same confidence-tier-passing pooling, each point divided by that trial+side's own L0 force (F/F0, itself gated by mfv_gate_f0()) | %d/%d trials contributed a usable F0 | %g%% length bins | connect-the-dots, NO fit\nPROTOTYPE (2026-07-21) -- compare against the RAW file before treating this as canonical",
    n_norm_avail_p, dplyr::n_distinct(passing$trial_id), STRAIN_BIN_PCT),
  y_lab = "Muscle force / trial's own L0 force (F/F0, dimensionless)")

outfile2 <- file.path(SUMMARY_OUTPUT_DIR, "FLsuperplot_isometric_isovelocity_pooled_snrPass.png")
ggplot2::ggsave(outfile2, p2, width = 13, height = 6.5, dpi = 150)
cli::cli_alert_success("Wrote {outfile2}")

outfile2_norm <- file.path(SUMMARY_OUTPUT_DIR, "FLsuperplot_isometric_isovelocity_pooled_normalized_snrPass.png")
ggplot2::ggsave(outfile2_norm, p2_norm, width = 13, height = 6.5, dpi = 150)
cli::cli_alert_success("Wrote {outfile2_norm}")
