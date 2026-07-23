# diag_snr_magnitude_audit.R
# READ-ONLY DIAGNOSTIC (2026-07-22, PI-requested audit): does the ratio-only
# activation_snr gate (MFV_UHAT_SNR_MIN = 3.0) conflate "measurement noise"
# with "genuinely small real force"? See analysis_muscle_force_vector_log.md's
# 2026-07-22 "SNR-based confidence gating audit" entry for the full writeup;
# this script only builds the evidence plots that entry references.
#
# WHY: activation_snr is a RATIO (||dF_force|| / baseline_force_noise_N).
# A low ratio can mean EITHER (a) the noise floor is elevated, OR (b) the
# true force is genuinely small -- the ratio alone cannot distinguish these.
# muscle_force_vector.R's mfv_gate_f0() already fixes this for the F0
# denominator specifically (requires BOTH snr >= snr_min AND |force| >= noise
# -- a magnitude floor). This script applies that SAME two-condition logic as
# a QUADRANT classification to every finalized step/hold/ramp across the full
# 3-fish corpus, so the conflation (or lack of it) is visible directly on the
# force-vs-noise-floor plane, not just inferred from the ratio.
#
# QUADRANTS (mirrors mfv_gate_f0's two conditions, cross-tabulated instead of
# collapsed to one pass/fail):
#   confident              : snr_pass  AND  magnitude_pass
#   confidently small      : !snr_pass AND  magnitude_pass  <- the conflation
#                             case this audit is about: ratio-only gating
#                             would call this "low confidence"/noise, but the
#                             force is independently above ITS OWN noise floor.
#   unstable magnitude      : snr_pass  AND  !magnitude_pass <- the case
#                             mfv_gate_f0 was built to catch (high activation
#                             on the ||dF|| VECTOR, near-zero on the PROJECTED
#                             scalar force_geom_N -- different numerators).
#   unconfirmable           : !snr_pass AND  !magnitude_pass <- genuinely
#                             indistinguishable from noise on both tests.
#
# CAVEAT (flagged, not glossed over): activation_snr's numerator is
# ||dF_force|| (the 3-axis active-minus-passive FORCE vector magnitude,
# .mfv_finalize_step()'s `emp$mag`), NOT force_geom_N (the PROJECTED scalar
# this script plots on the y-axis and that mfv_gate_f0/most downstream
# consumers actually use). The two can diverge (that divergence is exactly
# why mfv_gate_f0 exists) -- so the dashed "ratio = 3" reference line here is
# illustrative of the SNR gate's threshold, not a literal recomputation of
# activation_snr from this plot's own y-axis.
#
# Canon token: snrmagnitudeaudit -- see FIGURES_README.md.
# Run: Rscript R/diag_snr_magnitude_audit.R

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(tibble); library(cli)
})

.root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
for (f in c("paths_config.R", "00_load_bender_flat.R", "01_calibrate.R", "02_deconvolve.R",
            "muscle_geometry.R", "fit_fv_fl.R", "03_analyze.R",
            "parse_trial_filename.R", "plot_strain_validation.R",
            "plot_angle_sono_validation.R", "muscle_force_vector.R",
            "plot_force_vs_time.R", "extract_dynamic_l0_bookends.R")) {
  source(file.path(.root, f))
}

# Reuse .load_trial_td() / .dynamic_l0_trial_rows() from superplot_fl_pooled.R
# (same pattern superplot_fv_pooled.R / diag_isovelocity_hillcheck.R already
# use) -- cut at "assemble" so only the helper defs are pulled in, not that
# script's own top-level extraction/plotting run.
.fl_pooled_lines <- readLines(file.path(.root, "superplot_fl_pooled.R"))
.cut_at <- grep("^# -+ assemble -+", .fl_pooled_lines)[1]
if (is.na(.cut_at)) cli::cli_abort("diag_snr_magnitude_audit.R: could not find the assemble-section marker in superplot_fl_pooled.R -- has its structure changed?")
.fl_pooled_defs_file <- tempfile(fileext = ".R")
writeLines(.fl_pooled_lines[seq_len(.cut_at - 1L)], .fl_pooled_defs_file)
source(.fl_pooled_defs_file)
unlink(.fl_pooled_defs_file)

SPECIMEN_DIRS <- c(
  bass16 = raw_source_dir(BASS16_RAW_SUBFOLDER),
  bass17 = raw_source_dir(BASS17_RAW_SUBFOLDER),
  bass18 = raw_source_dir(BASS18_RAW_SUBFOLDER)
)
OUT <- FIGS_DIAGNOSTIC_DIR
fs::dir_create(OUT, recurse = TRUE)
snr_min <- MFV_UHAT_SNR_MIN

# =============================================================================
# Full-corpus extraction: EVERY finalized step/hold/ramp, ALL strains/
# velocities (unlike the pooled superplots, which restrict to V=0 for FL or
# add moving ramps only for FV) -- this script wants the whole force-vs-noise
# landscape, not a curve.
# =============================================================================

#' Isometric: every left/right step, all strains (not just L0).
.extract_isometric_all <- function(fish, manifest) {
  iso <- dplyr::filter(manifest, protocol == "isometric")
  out <- list()
  for (i in seq_len(nrow(iso))) {
    tid <- iso$trial_id[i]; fp <- iso$fullpath[i]
    res <- tryCatch({ td <- .load_trial_td(fp); r <- analyze_isometric(td, filename = fp); r$trial_id <- tid; r },
                    error = function(e) { cli::cli_warn("iso load {tid}: {conditionMessage(e)}"); NULL })
    if (is.null(res)) next
    vec <- tryCatch(attach_vector_muscle_force(res, fp, "isometric"),
                    error = function(e) { cli::cli_warn("iso vec {tid}: {conditionMessage(e)}"); NULL })
    if (is.null(vec)) next
    ss <- dplyr::filter(vec$step_summary, .data$muscle_side %in% c("left", "right"))
    if (nrow(ss) == 0L) next
    out[[tid]] <- ss |>
      dplyr::transmute(fish = fish, trial_id = tid, category = "isometric",
                       muscle_side = .data$muscle_side,
                       shortening_strain_pct = .data$shortening_strain_pct,
                       force_emp_N = .data$muscle_force_vector_N,
                       force_geom_N = .data$muscle_force_vector_geom_N,
                       activation_snr = .data$activation_snr,
                       baseline_force_noise_N = .data$baseline_force_noise_N)
  }
  if (length(out) == 0L) tibble::tibble() else dplyr::bind_rows(out)
}

#' Isovelocity: every left/right V=0 hold AND moving (angle-matched only,
#' matching the FV superplot's own exclusion of static_baseline_fallback rows
#' -- see ytorquesignexamples' 2026-07-22 correction for why) ramp, via the
#' SAME cross-trial batch the production FV path uses.
.extract_isovelocity_all <- function(fish, manifest) {
  isv <- dplyr::filter(manifest, protocol == "isovelocity")
  if (nrow(isv) == 0L) return(tibble::tibble())
  iso_inputs <- list()
  for (i in seq_len(nrow(isv))) {
    tid <- isv$trial_id[i]; fp <- isv$fullpath[i]
    res <- tryCatch({ td <- .load_trial_td(fp); r <- analyze_isovelocity(td, filename = fp); r$trial_id <- tid; r },
                    error = function(e) { cli::cli_warn("isv load {tid}: {conditionMessage(e)}"); NULL })
    if (!is.null(res)) iso_inputs[[tid]] <- list(trial_id = tid, filename = fp, res = res)
  }
  if (length(iso_inputs) == 0L) return(tibble::tibble())
  batch <- tryCatch(compute_isovelocity_vector_batch(iso_inputs),
                    error = function(e) { cli::cli_warn("{fish} isovelocity batch: {conditionMessage(e)}"); NULL })
  if (is.null(batch)) return(tibble::tibble())
  out <- list()
  for (tid in names(batch$step_summaries)) {
    ss <- batch$step_summaries[[tid]] |>
      dplyr::filter(.data$muscle_side %in% c("left", "right")) |>
      dplyr::mutate(category = dplyr::case_when(
        .data$contraction_mode == "isometric_zero" ~ "isovelocity_V0",
        .data$contraction_mode %in% c("concentric", "eccentric") &
          grepl("^angle_matched", .data$passive_source) ~ "isovelocity_moving",
        .default = NA_character_))
    ss <- dplyr::filter(ss, !is.na(.data$category))
    if (nrow(ss) == 0L) next
    out[[tid]] <- ss |>
      dplyr::transmute(fish = fish, trial_id = tid, category = .data$category,
                       muscle_side = .data$muscle_side,
                       shortening_strain_pct = .data$shortening_strain_pct,
                       force_emp_N = .data$muscle_force_vector_N,
                       force_geom_N = .data$muscle_force_vector_geom_N,
                       activation_snr = .data$activation_snr,
                       baseline_force_noise_N = .data$baseline_force_noise_N)
  }
  if (length(out) == 0L) tibble::tibble() else dplyr::bind_rows(out)
}

#' Dynamic L0 bookends -- reuses .dynamic_l0_trial_rows() (superplot_fl_pooled.R).
.extract_dynamic_all <- function(fish, manifest) {
  dyn <- dplyr::filter(manifest, protocol == "dynamic")
  out <- list()
  for (i in seq_len(nrow(dyn))) {
    tid <- dyn$trial_id[i]; fp <- dyn$fullpath[i]
    rows <- tryCatch(.dynamic_l0_trial_rows(fp, tid),
                     error = function(e) { cli::cli_warn("dyn load {tid}: {conditionMessage(e)}"); NULL })
    if (is.null(rows) || nrow(rows) == 0L) next
    out[[tid]] <- rows |>
      dplyr::transmute(fish = fish, trial_id = tid, category = "dynamic_bookend",
                       muscle_side = .data$muscle_side,
                       shortening_strain_pct = 0.0,
                       force_emp_N = .data$muscle_force_vector_N,
                       force_geom_N = .data$muscle_force_vector_geom_N,
                       activation_snr = .data$activation_snr,
                       baseline_force_noise_N = .data$baseline_force_noise_N)
  }
  if (length(out) == 0L) tibble::tibble() else dplyr::bind_rows(out)
}

# ---------------------------------------------------------------- assemble ---

cli::cli_h1("SNR/magnitude conflation audit: force vs noise floor, full corpus")
all_rows <- list()
for (fish in names(SPECIMEN_DIRS)) {
  cli::cli_h2(fish)
  manifest <- parse_trial_directory(SPECIMEN_DIRS[[fish]])
  iso <- .extract_isometric_all(fish, manifest)
  isv <- .extract_isovelocity_all(fish, manifest)
  dyn <- .extract_dynamic_all(fish, manifest)
  cli::cli_alert_info("{fish}: isometric = {nrow(iso)}; isovelocity (V0+moving) = {nrow(isv)}; dynamic bookend = {nrow(dyn)}")
  all_rows[[fish]] <- dplyr::bind_rows(iso, isv, dyn)
}
d <- dplyr::bind_rows(all_rows) |>
  dplyr::filter(is.finite(.data$force_geom_N), is.finite(.data$activation_snr),
               is.finite(.data$baseline_force_noise_N), .data$baseline_force_noise_N > 0)
if (nrow(d) == 0L) cli::cli_abort("No finite rows extracted -- nothing to audit.")

CAT_LEVELS <- c("isometric", "isovelocity_V0", "isovelocity_moving", "dynamic_bookend")
d$category <- factor(d$category, levels = CAT_LEVELS)

d <- d |>
  dplyr::mutate(
    snr_pass = .data$activation_snr >= snr_min,
    mag_pass = abs(.data$force_geom_N) >= .data$baseline_force_noise_N,
    quadrant = dplyr::case_when(
      .data$snr_pass  & .data$mag_pass  ~ "confident",
      !.data$snr_pass & .data$mag_pass  ~ "confidently small (SNR fails, magnitude passes)",
      .data$snr_pass  & !.data$mag_pass ~ "unstable magnitude (SNR passes, magnitude fails)",
      TRUE                              ~ "unconfirmable (both fail)"))

QUAD_LEVELS <- c("confident", "confidently small (SNR fails, magnitude passes)",
                 "unstable magnitude (SNR passes, magnitude fails)", "unconfirmable (both fail)")
QUAD_PAL <- c("confident" = "#1b9e77",
             "confidently small (SNR fails, magnitude passes)" = "#7570b3",
             "unstable magnitude (SNR passes, magnitude fails)" = "#e7298a",
             "unconfirmable (both fail)" = "#999999")
d$quadrant <- factor(d$quadrant, levels = QUAD_LEVELS)

# =============================================================================
# Plot 1: force vs noise floor directly, log-log, quadrant-colored
# =============================================================================
ref_lines <- tibble::tibble(
  x = 10^seq(log10(min(d$baseline_force_noise_N)), log10(max(d$baseline_force_noise_N)), length.out = 50))
p1 <- ggplot(d, aes(baseline_force_noise_N, abs(force_geom_N))) +
  geom_abline(slope = 1, intercept = 0, color = "grey35", linewidth = 0.5) +
  geom_abline(slope = snr_min, intercept = 0, color = "grey35", linetype = "dashed", linewidth = 0.5) +
  geom_point(aes(color = quadrant, shape = fish), size = 1.7, alpha = 0.75) +
  annotate("text", x = min(d$baseline_force_noise_N), y = min(d$baseline_force_noise_N),
           label = "ratio = 1\n(magnitude floor)", hjust = 0, vjust = 1.3, size = 2.7, color = "grey35") +
  annotate("text", x = min(d$baseline_force_noise_N), y = snr_min * min(d$baseline_force_noise_N),
           label = sprintf("ratio = %.0f\n(SNR gate, illustrative)", snr_min), hjust = 0, vjust = -0.3, size = 2.7, color = "grey35") +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = QUAD_PAL, name = "quadrant") +
  facet_wrap(~category, nrow = 1) +
  labs(title = "Force magnitude vs. its OWN noise floor -- does a low activation_snr ratio mean noise, or small-but-real force?",
      subtitle = sprintf(
        "Log-log; each point is one finalized step/hold/ramp (|muscle_force_vector_geom_N| vs baseline_force_noise_N), all 3 fish, full corpus (not just V=0/pooled subsets) | n=%d | solid line = magnitude floor (ratio=1); dashed = SNR gate threshold (ratio=%.0f, ILLUSTRATIVE -- activation_snr's real numerator is ||dF_force||, not this axis's force_geom_N, see script header caveat)\n'confidently small' (purple) = points a ratio-only gate would call low-confidence/noise, but whose force independently clears its own noise floor -- exactly the conflation this audit is checking for",
        nrow(d), snr_min),
      x = "baseline_force_noise_N (pre-stim ||force-channel|| SD, N)",
      y = "|muscle_force_vector_geom_N| (N)") +
  theme_bw(10) + theme(plot.subtitle = element_text(size = 7.3), legend.position = "bottom")
ggsave(file.path(OUT, "snrmagnitudeaudit_1_forceVsNoiseFloor.png"), p1, width = 14, height = 5, dpi = 150)

# =============================================================================
# Plot 2: quadrant counts by category -- how big is this, not just anecdotal?
# =============================================================================
counts <- d |> dplyr::count(category, quadrant)
p2 <- ggplot(counts, aes(category, n, fill = quadrant)) +
  geom_col(position = "stack") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 2.6, color = "white") +
  scale_fill_manual(values = QUAD_PAL, name = "quadrant") +
  coord_flip() +
  labs(title = "How many points fall in each SNR-x-magnitude quadrant, by category?",
      subtitle = "Same classification as snrmagnitudeaudit_1 -- quantifies the conflation across the whole corpus rather than one anecdote.",
      x = NULL, y = "n points") +
  theme_bw(11) + theme(legend.position = "bottom")
ggsave(file.path(OUT, "snrmagnitudeaudit_2_quadrantCounts.png"), p2, width = 9, height = 5, dpi = 150)

# =============================================================================
# Plot 3: isovelocity V=0 notch, re-examined -- the isovhillcheck.png example
# =============================================================================
v0 <- dplyr::filter(d, .data$category == "isovelocity_V0")
p3 <- ggplot(v0, aes(baseline_force_noise_N, abs(force_geom_N))) +
  geom_abline(slope = 1, intercept = 0, color = "grey35", linewidth = 0.5) +
  geom_abline(slope = snr_min, intercept = 0, color = "grey35", linetype = "dashed", linewidth = 0.5) +
  geom_point(aes(color = quadrant, shape = muscle_side), size = 2.6, alpha = 0.85) +
  scale_x_log10() + scale_y_log10() +
  scale_color_manual(values = QUAD_PAL, name = "quadrant") +
  facet_wrap(~fish, nrow = 1) +
  labs(title = "Isovelocity V=0 holds re-examined: is the isovhillcheck V=0 notch noise, or small-but-real force?",
      subtitle = "Same axes/lines as snrmagnitudeaudit_1, isovelocity_V0 only. bass18's SNR-passing V=0 reps (weak trial, ~0.25-0.62 N) vs SNR-failing V=0 reps of a STRONGER trial (~1.69/2.41 N, consistent with its eccentric plateau) should land in DIFFERENT quadrants if this is a magnitude story, not a pure-noise one -- see console table + decision log for the exact per-trial numbers.",
      x = "baseline_force_noise_N (N)", y = "|muscle_force_vector_geom_N| (N)") +
  theme_bw(10) + theme(plot.subtitle = element_text(size = 7.3), legend.position = "bottom")
ggsave(file.path(OUT, "snrmagnitudeaudit_3_isovelocityV0notch.png"), p3, width = 12, height = 5, dpi = 150)

cli::cli_alert_success("Wrote snrmagnitudeaudit_1/2/3 to {OUT}")

# ---- console tables for the decision log -------------------------------
cli::cli_h2("Quadrant counts by category")
print(as.data.frame(counts), row.names = FALSE)

cli::cli_h2("bass18 isovelocity V=0 holds, per trial+side (the isovhillcheck notch)")
v0_b18 <- v0 |>
  dplyr::filter(.data$fish == "bass18") |>
  dplyr::arrange(.data$trial_id, .data$muscle_side) |>
  dplyr::transmute(trial_id, muscle_side,
                   force_geom_N = round(force_geom_N, 4),
                   activation_snr = round(activation_snr, 2),
                   baseline_force_noise_N = round(baseline_force_noise_N, 4),
                   snr_pass, mag_pass, quadrant)
print(as.data.frame(v0_b18), row.names = FALSE)

cli::cli_h2("All-fish isovelocity V=0 quadrant breakdown")
print(as.data.frame(v0 |> dplyr::count(fish, quadrant)), row.names = FALSE)
