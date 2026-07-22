# diag_ytorque_sign_examples.R
# READ-ONLY DIAGNOSTIC (2026-07-22, PI-requested). Does NOT modify any
# pipeline analysis: sources existing functions and only reads data + writes
# PNGs. Two canon tokens, one script:
#
# CORRECTED same day (PI follow-up: isovelocity's negative y-torque "happens
# well after the stimulus, likely velocity-ramp related -- inertial noise?"):
# v1 of this script captured isovelocity's non-V0 rows via extract_
# isovelocity_zero_points()'s internal per-step loop, which uses a STATIC
# pre-stim baseline for EVERY step (including moving ones) before discarding
# everything except contraction_mode=="isometric_zero" -- comparing a
# fast-moving window against a motionless pre-stim window is guaranteed to
# show a huge motion-linked deviation, unrelated to the REAL isovelocity path
# (compute_isovelocity_vector_batch(), ANGLE-matched passive, what actually
# feeds FV_isovelocity_uhatBoth.png). Now calls compute_isovelocity_vector_
# batch() directly per fish for the moving-step rows -- see
# ytorqueinertialtiming.png (R/diag_ytorque_inertial_timing.R) for the
# follow-up investigation this correction fed into: even with CORRECT
# angle-matched passive, isovelocity's moving steps still show a residual
# sign that is much closer to a coin flip (52-60% negative) than isometric/
# V0 holds' consistent 90-97% positive -- and the residual's magnitude/timing
# is a kinematic (inertial-coupling) feature of the ramp itself, present
# identically in a completely unstimulated no-stim ramp of the same speed.
#
#   1) ytorquesignexamples -- 3 real POSITIVE and 3 real NEGATIVE
#      force_yTorque_N examples (raw Ty trace, active window shaded, the
#      ACTUAL passive/active means the pipeline used marked) -- visual
#      evidence for the empirical finding (analysis_muscle_force_vector_log.md,
#      2026-07-22) that force_yTorque_N is NOT side-mirrored and NOT an L0
#      artifact: it splits almost entirely on STATIC (isometric hold, or
#      isovelocity's own V=0 holds -- consistently positive) vs. ACTIVELY
#      MOVING (isovelocity's concentric/eccentric ramp, using the REAL
#      angle-matched passive -- much closer to a coin flip than static
#      holds' consistent sign; see ytorqueinertialtiming.png for the
#      mechanism). All 6 examples are chosen with real, defensible SNR
#      (>= 4) so the sign shown is a genuine feature of the trace, not
#      pure noise -- exact pooled percentages are computed fresh each run
#      and printed in the figure's own subtitle, not hardcoded here.
#
#   2) axissnrcomparison -- for every real finalized step (328, 3 fish, all
#      3 categories), SNR computed identically (|force_axis_N| / the SAME
#      force-domain noise floor, .mfv_baseline_force_noise()) for all 4
#      axes side by side (moment/force/yTorque/zTorque) -- shows which axes
#      are trustworthy relative to each other, not just individually. Same
#      underlying "noise" this pipeline already uses for activation_snr, just
#      applied to each axis's own numerator instead of only the empirical
#      u_hat force magnitude.
#
# Neither panel re-implements any pipeline physics: both come from
# INSTRUMENTING (not copying) attach_vector_muscle_force()/
# compute_isovelocity_vector_batch()'s real calls to .mfv_finalize_step(),
# capturing its actual inputs/outputs -- avoids exactly the kind of silent
# reimplementation drift that caused the forceextractionmethod cross-step-
# contamination bug (see FIGURES_README.md).
#
# Run: Rscript R/diag_ytorque_sign_examples.R

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(ggplot2)
  library(cli); library(patchwork); library(tidyr)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")
OUT_DIR <- FIGS_DIAGNOSTIC_DIR

# Pull in superplot_fl_pooled.R's extraction helpers (extract_isometric_points()
# etc.) WITHOUT its bottom-of-file assemble/plot section -- these already walk
# every real trial of all 3 fish through attach_vector_muscle_force()/
# compute_isovelocity_vector_batch(), which is exactly the real call path we
# want to instrument (not re-derive).
.sp_lines <- readLines(file.path(.pipeline_root, "superplot_fl_pooled.R"))
.cut_at <- grep("^# -+ assemble -+", .sp_lines)[1]
if (is.na(.cut_at)) cli::cli_abort("superplot_fl_pooled.R's assemble-section marker not found -- check for upstream edits")
.defs_tmp <- tempfile(fileext = ".R")
writeLines(.sp_lines[seq_len(.cut_at - 1)], .defs_tmp)
source(.defs_tmp)  # brings in its own full source() chain, incl. muscle_force_vector.R

# =============================================================================
# Instrument .mfv_finalize_step() (capture real inputs/outputs, change nothing)
# =============================================================================
.capture <- list()
.passive_src <- list()
.orig_finalize <- .mfv_finalize_step
.mfv_finalize_step <- function(act, pass, noise, category, s, td6, step_rows, arms, geom,
                               trial_id, snr_min, deactivation_window_s,
                               baseline_pad_s, relaxation_s, sono_ctx = NULL) {
  fin <- .orig_finalize(act, pass, noise, category, s, td6, step_rows, arms, geom,
                        trial_id, snr_min, deactivation_window_s, baseline_pad_s, relaxation_s, sono_ctx)
  r <- fin$row
  denom <- if (is.finite(noise) && noise > 0) noise else NA_real_
  .capture[[length(.capture) + 1]] <<- tibble::tibble(
    fish = .CURRENT_FISH, category = category, side = s$muscle_side, cmode = s$contraction_mode,
    fp = .CURRENT_FP, step_number = s$step_number,
    stim_t0_s = s$stim_t0_s, stim_t1_s = s$stim_t1_s, deactivation_window_s = deactivation_window_s,
    pass_Ty = pass$T[2L], act_Ty = act$T[2L], op = s$operating_point,
    muscle_force_vector_N = r$muscle_force_vector_N, force_force_channel_N = r$force_force_channel_N,
    force_yTorque_N = r$force_yTorque_N, force_zTorque_N = r$force_zTorque_N,
    noise = denom, passive_source = NA_character_,  # back-filled by .mfv_assign_row() wrapper below
    snr_moment = abs(r$muscle_force_vector_N) / denom, snr_force = abs(r$force_force_channel_N) / denom,
    snr_yT = abs(r$force_yTorque_N) / denom, snr_zT = abs(r$force_zTorque_N) / denom
  )
  fin
}
.orig_load_trial_td <- .load_trial_td
.load_trial_td <- function(fullpath) { .CURRENT_FP <<- fullpath; .orig_load_trial_td(fullpath) }
# .mfv_assign_row() is the only place passive_source is known FOR THE ROW JUST
# FINALIZED, but it runs AFTER .mfv_finalize_step() -- track it in a side
# channel .LAST_PASSIVE_SOURCE so the capture above (inside finalize_step, one
# call earlier) can be back-filled correctly since calls happen in lockstep,
# one finalize_step() then one assign_row() per step, never interleaved.
.orig_assign_row <- .mfv_assign_row
.mfv_assign_row <- function(out_cols, i, fin, passive_source, angle_ok) {
  if (length(.capture) > 0L) .capture[[length(.capture)]]$passive_source <<- passive_source
  .orig_assign_row(out_cols, i, fin, passive_source, angle_ok)
}

SPECIMEN_DIRS_DIAG <- list(bass16 = raw_source_dir(BASS16_RAW_SUBFOLDER),
                           bass17 = raw_source_dir(BASS17_RAW_SUBFOLDER),
                           bass18 = raw_source_dir(BASS18_RAW_SUBFOLDER))
.CURRENT_FISH <- NA_character_; .CURRENT_FP <- NA_character_
cli::cli_h1("Instrumented pass over all 3 fish (real pipeline calls, capture-only)")
for (fish in names(SPECIMEN_DIRS_DIAG)) {
  .CURRENT_FISH <- fish
  manifest <- parse_trial_directory(SPECIMEN_DIRS_DIAG[[fish]])
  invisible(extract_isometric_points(fish, manifest))
  invisible(extract_dynamic_l0_points(fish, manifest))
  # NOTE (2026-07-22 correction): extract_isovelocity_zero_points() is
  # DELIBERATELY excluded from this instrumented pass -- its per-step loop
  # calls attach_vector_muscle_force() (STATIC pre-stim baseline) on EVERY
  # step of the trial, including non-V0 (moving) ones, before filtering down
  # to just contraction_mode=="isometric_zero" in its RETURN. For a moving
  # ramp, that static baseline compares a fast-moving active window against
  # a MOTIONLESS pre-stim window -- guaranteed to show a huge motion-linked
  # deviation regardless of muscle activity, NOT what force_yTorque_N
  # actually reports for those steps in the real pipeline. The REAL
  # isovelocity path (what feeds FV_isovelocity_uhatBoth.png) is
  # compute_isovelocity_vector_batch(), which uses ANGLE-MATCHED passive
  # instead -- called separately below, once per fish (it batches all of
  # that fish's isovelocity trials together for cross-trial passive
  # borrowing, so it cannot be called per-trial inside this loop).
  isv <- dplyr::filter(manifest, protocol == "isovelocity")
  iso_inputs <- list()
  for (i in seq_len(nrow(isv))) {
    tid <- isv$trial_id[i]; fpi <- isv$fullpath[i]
    res_i <- tryCatch({ td_i <- .load_trial_td(fpi); r_i <- analyze_isovelocity(td_i, filename = fpi); r_i$trial_id <- tid; r_i },
                      error = function(e) NULL)
    if (!is.null(res_i)) iso_inputs[[tid]] <- list(trial_id = tid, filename = fpi, res = res_i)
  }
  if (length(iso_inputs) > 0L) invisible(compute_isovelocity_vector_batch(iso_inputs))
}
allr <- dplyr::bind_rows(.capture) |> dplyr::filter(is.finite(force_yTorque_N))
# For isovelocity MOVING steps (concentric/eccentric), ONLY keep angle_matched
# (/_cross_trial) rows -- the fallback static_baseline path has the same
# motion-vs-motionless problem described above. isovelocity's OWN V=0 holds
# (contraction_mode=="isometric_zero") are exempt from this filter: they are
# genuinely motionless, so compute_isovelocity_vector_batch() correctly
# ALWAYS uses static_baseline_fallback for them (there is no ramp to
# angle-match against a hold) -- that is the appropriate baseline, not a bug.
allr <- allr |> dplyr::filter(
  category != "isovelocity" | cmode == "isometric_zero" | grepl("^angle_matched", passive_source)
)
cli::cli_alert_success("Captured {nrow(allr)} real finalized steps (angle-matched-only for isovelocity) across {dplyr::n_distinct(allr$fish)} fish")

# =============================================================================
# 1) ytorquesignexamples -- 3 positive + 3 negative, real, high-SNR
# =============================================================================
#' Pick the single highest-snr_yT row matching a filter, excluding dynamic
#' bookend files (those need a different, runtime-detected baseline window --
#' see extract_dynamic_l0_bookends.R -- kept out here so all 6 examples share
#' one simple segmented-file plotting path).
.pick <- function(df, ...) {
  df |> dplyr::filter(!grepl("dynamic", fp), snr_yT >= 4, ...) |> dplyr::arrange(dplyr::desc(snr_yT)) |> head(1)
}
sel <- dplyr::bind_rows(
  .pick(allr, fish == "bass18", side == "right", cmode == "concentric",  category == "isometric",   force_yTorque_N > 0),
  .pick(allr, fish == "bass18", side == "left",  cmode == "concentric",  category == "isometric",   force_yTorque_N > 0),
  .pick(allr, fish == "bass16", cmode %in% c("concentric", "eccentric"), category == "isometric",   force_yTorque_N > 0),
  .pick(allr, fish == "bass16", side == "left",  cmode == "eccentric",   category == "isovelocity", force_yTorque_N < 0),
  .pick(allr, fish == "bass17", side == "right", cmode %in% c("concentric", "eccentric"), category == "isovelocity", force_yTorque_N < 0),
  .pick(allr, fish == "bass17", side == "left",  cmode %in% c("concentric", "eccentric"), category == "isovelocity", force_yTorque_N < 0)
)
sel$example_label <- c("Positive 1", "Positive 2", "Positive 3", "Negative 1", "Negative 2", "Negative 3")
if (nrow(sel) != 6L || any(!is.finite(sel$force_yTorque_N))) {
  cli::cli_abort("ytorquesignexamples: expected 6 valid example rows, got {nrow(sel)} -- check .pick() filters against current data")
}
cli::cli_alert_info("Selected examples:")
print(as.data.frame(sel[, c("example_label", "fish", "category", "side", "cmode", "op",
                            "force_yTorque_N", "snr_yT")]))

#' One example panel: raw Ty trace (light smoothing overlay for legibility
#' only), active window shaded, the ACTUAL pass_Ty/act_Ty means the real
#' pipeline call used drawn as horizontal lines (not recomputed -- read
#' directly from the instrumented capture above, so this is guaranteed
#' consistent with the real force_yTorque_N value shown).
.build_example_panel <- function(row) {
  # NOTE: .load_trial_td() (superplot_fl_pooled.R) only loads X-torque (for
  # inertial correction, feeding the isometric/isovelocity FL analysis path)
  # -- it has no ytorque column at all. force_yTorque_N needs the full 6-axis
  # load instead (same one mfv_load_six_axis()/attach_vector_muscle_force()
  # use), which needs no deconvolution (isometric/isovelocity holds are
  # motionless -- see muscle_force_vector.R's module header).
  td <- load_bender_flat(row$fp, do_filter = TRUE, loadtorques = c("x", "y", "z"), loadforces = TRUE)
  Ty <- .mfv_col(td, "ytorque")
  if (is.null(Ty)) cli::cli_abort("No ytorque column in {basename(row$fp)}")
  pad_s <- 0.3
  win_end <- row$stim_t1_s - row$stim_t0_s + row$deactivation_window_s
  # step_number filter is REQUIRED for segmented files -- t.s resets to 0 at
  # the start of EVERY step (verified 2026-07-21, see diag_force_extraction_
  # baseline.R's .raw_window() docstring) -- without it, a naive time-range
  # filter silently pools rows from every step in the file.
  keep <- td$t.s >= (row$stim_t0_s - pad_s) & td$t.s <= (row$stim_t1_s + row$deactivation_window_s + pad_s) &
    td$step_number == row$step_number
  raw <- tibble::tibble(t_rel = td$t.s[keep] - row$stim_t0_s, Ty_Nm = Ty[keep])
  raw$Ty_smooth_Nm <- .smooth_trace_display_only(raw$Ty_Nm)

  sign_lab <- if (row$force_yTorque_N > 0) "POSITIVE" else "NEGATIVE"
  sign_col <- if (row$force_yTorque_N > 0) "#1d4ed8" else "#b91c1c"
  wrap40 <- function(s) paste(strwrap(s, width = 40), collapse = "\n")
  ggplot(raw, aes(t_rel, Ty_Nm)) +
    annotate("rect", xmin = 0, xmax = win_end, ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.12) +
    geom_vline(xintercept = c(0, win_end), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_line(color = "grey75", linewidth = 0.3) +
    geom_line(aes(y = Ty_smooth_Nm), color = "grey15", linewidth = 0.7) +
    geom_hline(yintercept = row$pass_Ty, color = "#7c3aed", linewidth = 0.8, linetype = "dotted") +
    geom_hline(yintercept = row$act_Ty, color = sign_col, linewidth = 1.0) +
    labs(
      title = wrap40(sprintf("%s: %s, %s, %s, op=%.1f", row$example_label, row$fish, row$category, row$side, row$op)),
      subtitle = wrap40(sprintf("force_yTorque_N = %+.3f N (SNR = %.1f) | dotted = passive Ty used, %s solid = active Ty used",
                                row$force_yTorque_N, row$snr_yT, sign_lab)),
      x = "Time relative to stim onset (s)", y = "Raw y-torque (N*m)"
    ) +
    theme_bw(base_size = 10) + theme(plot.title = element_text(size = 9.5), plot.subtitle = element_text(size = 7.5))
}

panels <- lapply(seq_len(nrow(sel)), function(i) .build_example_panel(sel[i, ]))
p_examples <- (panels[[1]] | panels[[2]] | panels[[3]]) / (panels[[4]] | panels[[5]] | panels[[6]]) +
  patchwork::plot_annotation(
    title = "force_yTorque_N sign: 3 real positive examples vs. 3 real negative examples",
    subtitle = paste(strwrap(
      sprintf("All 6 have real SNR >= 4 (not coin-flip noise; isovelocity's angle-matched SNR runs lower than isometric's, see axissnrcomparison.png). Top row (positive) = static holds (isometric, or isovelocity's own V=0 rep). Bottom row (negative) = isovelocity's actively-moving concentric/eccentric ramp, using the REAL angle-matched passive (compute_isovelocity_vector_batch(), same as FV_isovelocity_uhatBoth.png -- NOT the static-baseline artifact v1 of this figure used). Pooled across %d real finalized steps: STATIC categories are %.0f%% positive, ACTIVELY-MOVING isovelocity is %.0f%% negative -- even with angle-matching, moving-ramp sign is much closer to a coin flip than static holds' consistent sign; see ytorqueinertialtiming.png for why (uncorrected inertial coupling from the ramp itself, present identically in an unstimulated ramp of the same speed).",
              nrow(allr),
              100 * mean(allr$force_yTorque_N[allr$category != "isovelocity" | allr$cmode == "isometric_zero"] > 0, na.rm = TRUE),
              100 * mean(allr$force_yTorque_N[allr$category == "isovelocity" & allr$cmode %in% c("concentric", "eccentric")] < 0, na.rm = TRUE)),
      width = 130), collapse = "\n")
  )
ggplot2::ggsave(file.path(OUT_DIR, "ytorquesignexamples.png"), p_examples, width = 15, height = 9, dpi = 150)
cli::cli_alert_success("Wrote ytorquesignexamples.png")

# =============================================================================
# 2) axissnrcomparison -- SNR across all 4 axes, same noise floor, all steps
# =============================================================================
snr_long <- allr |>
  dplyr::filter(dplyr::if_all(c(snr_moment, snr_force, snr_yT, snr_zT), is.finite)) |>
  tidyr::pivot_longer(c(snr_moment, snr_force, snr_yT, snr_zT), names_to = "axis", values_to = "snr") |>
  dplyr::mutate(axis = factor(axis,
                              levels = c("snr_moment", "snr_force", "snr_yT", "snr_zT"),
                              labels = c("moment\n(muscle_force_vector_N)", "force channel\n(dF . u_hat)",
                                        "yTorque only\n(dTy/d)", "zTorque only\n(dTz/r_m)")))

p_snr_box <- ggplot(snr_long, aes(axis, snr, fill = axis)) +
  geom_boxplot(outlier.alpha = 0.25, width = 0.6) +
  geom_hline(yintercept = MFV_UHAT_SNR_MIN, linetype = "dashed", color = "grey30") +
  annotate("text", x = 0.6, y = MFV_UHAT_SNR_MIN, label = sprintf("SNR = %.0f (empirical u_hat gate)", MFV_UHAT_SNR_MIN),
          hjust = 0, vjust = -0.4, size = 3, color = "grey30") +
  scale_y_log10() +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(title = "SNR by axis, same force-domain noise floor for all 4 (log scale)",
      subtitle = sprintf("All %d real finalized steps, 3 fish, isometric + isovelocity-V0 + dynamic bookends pooled", nrow(snr_long) / 4L),
      x = NULL, y = "SNR = |force_axis_N| / noise (log10)") +
  theme_bw(base_size = 11)

# Per-step paired comparison: does a step's ranking of axes hold up across
# categories, or does e.g. yTorque only look good on ISOMETRIC and fall apart
# on ISOVELOCITY specifically?
p_snr_cat <- ggplot(snr_long, aes(axis, snr, fill = category)) +
  geom_boxplot(outlier.alpha = 0.2, position = position_dodge(0.75), width = 0.65) +
  scale_y_log10() +
  scale_fill_manual(values = c(isometric = "#1b9e77", isovelocity = "#d95f02", dynamic = "#7570b3")) +
  labs(title = "Same comparison, split by trial category",
      x = NULL, y = "SNR (log10)", fill = "Category") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

p_snr <- p_snr_box / p_snr_cat +
  patchwork::plot_annotation(
    title = "Axis-to-axis SNR comparison (all use the SAME baseline noise floor)")
ggplot2::ggsave(file.path(OUT_DIR, "axissnrcomparison.png"), p_snr, width = 10, height = 11, dpi = 150)
cli::cli_alert_success("Wrote axissnrcomparison.png")

med <- snr_long |> dplyr::group_by(axis) |> dplyr::summarise(median_snr = median(snr), .groups = "drop")
med_lab <- gsub("\n", " ", as.character(med$axis))
cli::cli_alert_info("Median SNR by axis: {paste(sprintf('%s=%.1f', med_lab, med$median_snr), collapse=', ')}")
