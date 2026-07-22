# diag_force_extraction_methods_compare.R
# READ-ONLY DIAGNOSTIC (2026-07-22, PI-requested follow-up to
# bass16_forceextractionmethod.png / diag_force_extraction_baseline.R).
#
# PI observation: "the max in window is probably NOT the best way to
# calculate muscle force. If anything, it should be calculated from the
# smoothed black line (or something like it)." bass16_forceextractionmethod.
# png's dynamic-bookend panel makes the case clearly: the smoothed trace
# rises to a clear peak (~0.22-0.3 N*m) then FADES within the same active
# window (twitch/tetanic decline, or simply a short stim burst followed by
# relaxation while the window stays open) -- the CURRENT method (MEAN over
# the WHOLE window, 0.09) dilutes that peak with the later off-peak decay,
# while raw MAX (0.31) is one noisy raw sample, not a robust peak.
#
# This script compares FOUR candidate scalar-extraction methods on the SAME
# raw traces, then checks how much the choice matters in aggregate across
# many real steps:
#   A) MEAN-of-full-window   -- CURRENT pipeline method (.mfv_window_means())
#   B) MAX-of-raw-window     -- single noisiest-possible sample (old
#                               hypothetical alternative in
#                               bass16_forceextractionmethod.png)
#   C) MAX-of-SMOOTHED-trace -- peak of the SAME display-smoothing already
#                               drawn on every trace in this pipeline's
#                               diagnostics (.smooth_trace_display_only()) --
#                               directly what the PI pointed at ("the
#                               smoothed black line")
#   D) NARROW-WINDOW-MEAN-AT-SMOOTHED-PEAK -- mean of RAW samples in a
#                               short window (0.15s) centered on the
#                               smoothed trace's peak time -- a compromise:
#                               robust to single-sample noise (unlike B),
#                               not diluted by post-peak decay (unlike A),
#                               and doesn't just report the smoothed line's
#                               own (already-attenuated-amplitude) value
#                               directly (unlike C).
#
# OUTCOME (2026-07-22, same day): PI picked D. Production now uses D (with
# a duration guard added afterward -- see muscle_force_vector.R's
# .mfv_window_peak_means() / 03_analyze.R's .legacy_peak_window_mean() and
# analysis_muscle_force_vector_log.md for the dynamic-bookend bug that
# guard fixes). This script itself is a READ-ONLY snapshot of the
# comparison AS IT WAS RUN that day and was NOT regenerated afterward --
# its own dynamic-bookend panel below predates the duration guard, so it
# does not reflect the fixed behavior for that specific window; the
# isometric/isovelocity panels and the sensitivity aggregate (isometric
# only) are unaffected by the guard and remain valid evidence for the D
# decision.
#
# Two tokens, one script:
#   1) muscleforcemethodcompare -- per-representative-step visual comparison
#      (mirrors bass16_forceextractionmethod.png's 3-panel layout, same
#      steps, now with all 4 candidates marked instead of 2).
#   2) muscleforcemethodsensitivity -- aggregate comparison across every
#      real SNR-passing isometric step (3 fish; isometric chosen because its
#      sustained hold is the cleanest test of "does the window include
#      post-peak decay", without isovelocity's separate motion-artifact
#      confound) -- how much would switching methods change the reported
#      force, and would it change the FL curve's SHAPE (not just scale)?
#
# Run: Rscript R/diag_force_extraction_methods_compare.R

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(ggplot2)
  library(cli); library(patchwork)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")
OUT_DIR <- FIGS_DIAGNOSTIC_DIR
src("00_load_bender_flat.R")
src("01_calibrate.R")
src("02_deconvolve.R")
src("muscle_geometry.R")
src("fit_fv_fl.R")
src("03_analyze.R")
src("parse_trial_filename.R")
src("plot_force_vs_time.R")  # .smooth_trace_display_only(), .detect_stim_events()
src("muscle_force_vector.R")
src("extract_dynamic_l0_bookends.R")

DEACTIVATION_WINDOW_S <- MFV_DEACTIVATION_WINDOW_S
PEAK_WINDOW_S <- 0.15  # method D's narrow window width, centered on the smoothed peak
PLOT_PAD_S <- 0.3
wrap_sub <- function(s, w = 110) paste(strwrap(s, width = w), collapse = "\n")

load_one <- function(f) {
  td <- load_bender_flat(f, do_filter = TRUE, loadtorques = "x")
  tau <- deconvolve_bender(f, hub_path = NULL, verbose = FALSE)
  N <- min(nrow(td), length(tau)); td <- td[seq_len(N), , drop = FALSE]
  td$torque_inertia_corrected_Nm <- tau[seq_len(N)]
  attr(td, "Filename") <- f
  td
}

#' Compute all 4 candidate scalar extractions from a raw window.
#' @param search_end End of the window searched for the PEAK (methods C/D) --
#'   defaults to win_end, but callers should pass the STIM duration alone
#'   (excluding the post-stim deactivation tail) when available: searching
#'   the full deactivation-extended window lets the trace's RECOVERY-toward-
#'   baseline transient right at the window's trailing edge get mistaken for
#'   a genuine peak (found empirically 2026-07-22, isometric panel -- a flat
#'   sustained contraction's smoothed trace drifts closer to baseline at the
#'   very edge of the deactivation tail, which can be "farther from the
#'   full-window mean" than the true plateau if the plateau itself is
#'   consistently offset from zero). A (mean_full) always still averages
#'   the FULL window regardless -- only the peak SEARCH range narrows.
#' @return list(mean_full, max_raw, max_smooth, narrow_mean, t_peak_rel,
#'   n_samples).
.extraction_candidates <- function(t_rel, torque_raw, win_end, peak_window_s = PEAK_WINDOW_S,
                                   search_end = win_end) {
  active <- t_rel >= 0 & t_rel <= win_end
  search <- t_rel >= 0 & t_rel <= search_end
  torque_smooth <- .smooth_trace_display_only(torque_raw)
  mean_full <- mean(torque_raw[active], na.rm = TRUE)
  max_val <- max(torque_raw[active], na.rm = TRUE); min_val <- min(torque_raw[active], na.rm = TRUE)
  max_raw <- if (abs(max_val - mean_full) >= abs(min_val - mean_full)) max_val else min_val
  search_idx <- which(search)
  smooth_search <- torque_smooth[search]
  # match the sign convention of max_raw: look for whichever extreme (peak
  # or trough of the SMOOTHED trace, WITHIN search_end only) is farther from
  # the full-window mean, so a genuinely negative-going contraction is still
  # tracked correctly.
  smax <- max(smooth_search, na.rm = TRUE); smin <- min(smooth_search, na.rm = TRUE)
  use_max <- abs(smax - mean_full) >= abs(smin - mean_full)
  max_smooth <- if (use_max) smax else smin
  t_peak_idx <- search_idx[if (use_max) which.max(smooth_search) else which.min(smooth_search)]
  t_peak_rel <- t_rel[t_peak_idx]
  narrow <- t_rel >= (t_peak_rel - peak_window_s / 2) & t_rel <= (t_peak_rel + peak_window_s / 2)
  narrow_mean <- mean(torque_raw[narrow], na.rm = TRUE)
  list(mean_full = mean_full, max_raw = max_raw, max_smooth = max_smooth, narrow_mean = narrow_mean,
      t_peak_rel = t_peak_rel, n_samples = sum(active, na.rm = TRUE))
}

.raw_window <- function(td, t0, t1, anchor_t, step_number) {
  keep <- td$t.s >= t0 & td$t.s <= t1
  if (!is.null(step_number)) keep <- keep & td$step_number == step_number
  tibble::tibble(t_rel = td$t.s[keep] - anchor_t, torque_Nm = td$torque_inertia_corrected_Nm[keep])
}

manifest <- parse_trial_directory(raw_source_dir(BASS16_RAW_SUBFOLDER))
iso_f <- manifest$fullpath[manifest$protocol == "isometric"][1]
iso_td <- load_one(iso_f)
iso_ss <- analyze_isometric(iso_td, filename = iso_f)$step_summary
iso_step <- iso_ss[iso_ss$step_number == 4, ][1, ]

isov_f <- manifest$fullpath[manifest$protocol == "isovelocity"][1]
isov_td <- load_one(isov_f)
isov_ss <- analyze_isovelocity(isov_td, filename = isov_f)$step_summary
isov_step <- isov_ss[isov_ss$step_number == 3, ][1, ]

dyn_f <- manifest$fullpath[manifest$protocol == "dynamic"][1]
dyn_td <- load_one(dyn_f)
bookends <- detect_dynamic_l0_bookends(dyn_td)
if (nrow(bookends) == 0L) {
  for (f in manifest$fullpath[manifest$protocol == "dynamic"][2:5]) {
    td_try <- load_one(f); be_try <- detect_dynamic_l0_bookends(td_try)
    if (nrow(be_try) > 0L) { dyn_f <- f; dyn_td <- td_try; bookends <- be_try; break }
  }
}
have_dynamic <- nrow(bookends) > 0L
dyn_step <- if (have_dynamic) bookends[1, ] else NULL

# =============================================================================
# 1) muscleforcemethodcompare -- per-step visual comparison, 4 methods
# =============================================================================
.build_compare_panel <- function(td, s, title, step_number = NULL) {
  dur <- s$stim_t1_s - s$stim_t0_s
  win_end <- dur + DEACTIVATION_WINDOW_S
  raw <- .raw_window(td, s$stim_t0_s - PLOT_PAD_S, s$stim_t1_s + DEACTIVATION_WINDOW_S + PLOT_PAD_S, s$stim_t0_s, step_number)
  active <- raw$t_rel >= 0 & raw$t_rel <= win_end
  cand <- .extraction_candidates(raw$t_rel[active], raw$torque_Nm[active], win_end, search_end = dur)
  raw$torque_smooth_Nm <- .smooth_trace_display_only(raw$torque_Nm)

  ggplot(raw, aes(t_rel, torque_Nm)) +
    annotate("rect", xmin = 0, xmax = win_end, ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.10) +
    annotate("rect", xmin = cand$t_peak_rel - PEAK_WINDOW_S / 2, xmax = cand$t_peak_rel + PEAK_WINDOW_S / 2,
             ymin = -Inf, ymax = Inf, fill = "#059669", alpha = 0.18) +
    geom_vline(xintercept = c(0, win_end), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_line(color = "grey75", linewidth = 0.3) +
    geom_line(aes(y = torque_smooth_Nm), color = "grey15", linewidth = 0.7) +
    geom_hline(aes(yintercept = cand$mean_full, color = "A: MEAN full window (CURRENT)"), linewidth = 1.0) +
    geom_hline(aes(yintercept = cand$max_raw, color = "B: MAX raw sample"), linewidth = 0.7, linetype = "dotted") +
    geom_hline(aes(yintercept = cand$max_smooth, color = "C: peak of SMOOTHED trace"), linewidth = 0.9, linetype = "dashed") +
    geom_hline(aes(yintercept = cand$narrow_mean, color = "D: narrow-window mean AT smoothed peak"), linewidth = 0.9) +
    scale_color_manual(name = NULL, values = c(
      "A: MEAN full window (CURRENT)" = "#1d4ed8", "B: MAX raw sample" = "#b91c1c",
      "C: peak of SMOOTHED trace" = "#b45309", "D: narrow-window mean AT smoothed peak" = "#059669")) +
    labs(title = title,
        subtitle = wrap_sub(sprintf("n=%d samples | green shade = D's %.2fs window at t=%.2fs | A=%.4g B=%.4g C=%.4g D=%.4g",
                           cand$n_samples, PEAK_WINDOW_S, cand$t_peak_rel, cand$mean_full, cand$max_raw, cand$max_smooth, cand$narrow_mean)),
        x = "Time relative to stim onset (s)", y = "Inertia-corrected torque (N*m)") +
    theme_bw(base_size = 10.5) + theme(plot.subtitle = element_text(size = 8), legend.position = "bottom")
}

p_iso <- .build_compare_panel(iso_td, iso_step, sprintf("Isometric (angle %.1f deg)", iso_step$operating_point), iso_step$step_number)
p_isov <- .build_compare_panel(isov_td, isov_step, sprintf("Isovelocity (%.1f deg/s)", isov_step$operating_point), isov_step$step_number)
p_dyn <- if (have_dynamic) {
  fake_s <- tibble::tibble(stim_t0_s = dyn_step$stim_t0_s, stim_t1_s = dyn_step$stim_t1_s)
  .build_compare_panel(dyn_td, fake_s, sprintf("Dynamic L0 bookend (%s)", dyn_step$muscle_side), NULL)
} else {
  ggplot() + annotate("text", x = 0, y = 0, label = "No dynamic bookend detected") + theme_void()
}
p_compare <- p_iso / p_isov / p_dyn +
  patchwork::plot_annotation(
    title = "Muscle-force extraction method comparison: 4 candidates on the same raw trace",
    subtitle = wrap_sub("bass16, one representative step per trial type | A (current) dilutes any post-peak decay across the WHOLE window; B is one noisy raw sample; C is the smoothed line's OWN peak value (attenuated amplitude by construction, low-pass filtering always shaves peaks); D is a robust plateau mean centered on the smoothed peak's TIMING but computed on RAW samples -- candidate replacement for A"))
ggplot2::ggsave(file.path(OUT_DIR, "muscleforcemethodcompare.png"), p_compare, width = 10, height = 12, dpi = 150)
cli::cli_alert_success("Wrote muscleforcemethodcompare.png")


# =============================================================================
# 2) muscleforcemethodsensitivity -- aggregate across every real SNR-passing
#    isometric step, 3 fish. Does the choice of method change the reported
#    force enough to matter, and does it change SHAPE (not just scale)?
# =============================================================================
.step_all_methods <- function(td, s, step_number) {
  dur <- s$stim_t1_s - s$stim_t0_s
  win_end <- dur + DEACTIVATION_WINDOW_S
  raw <- .raw_window(td, s$stim_t0_s, s$stim_t1_s + DEACTIVATION_WINDOW_S, s$stim_t0_s, step_number)
  if (nrow(raw) < 10L) return(NULL)
  cand <- .extraction_candidates(raw$t_rel, raw$torque_Nm, win_end, search_end = dur)
  # baseline (pre-stim static mean, unaffected by this comparison -- same
  # subtraction for all 4 candidates) so the comparison is on a genuine
  # ACTIVE-minus-PASSIVE force, not raw torque.
  pre <- td$t.s >= s$t_pre_baseline_start_s & td$t.s <= s$t_pre_baseline_end_s &
    (if (!is.null(step_number)) td$step_number == step_number else TRUE)
  base <- mean(td$torque_inertia_corrected_Nm[pre], na.rm = TRUE)
  tibble::tibble(step_number = s$step_number, muscle_side = s$muscle_side,
                strain_pct = s$shortening_strain_pct,
                A_mean_full = cand$mean_full - base, B_max_raw = cand$max_raw - base,
                C_smooth_peak = cand$max_smooth - base, D_narrow_mean = cand$narrow_mean - base)
}

SPECIMEN_DIRS_SENS <- list(bass16 = raw_source_dir(BASS16_RAW_SUBFOLDER), bass17 = raw_source_dir(BASS17_RAW_SUBFOLDER),
                          bass18 = raw_source_dir(BASS18_RAW_SUBFOLDER))
sens_rows <- list()
for (fish in names(SPECIMEN_DIRS_SENS)) {
  m <- parse_trial_directory(SPECIMEN_DIRS_SENS[[fish]])
  iso_files <- dplyr::filter(m, protocol == "isometric")
  for (i in seq_len(nrow(iso_files))) {
    fp <- iso_files$fullpath[i]; tid <- iso_files$trial_id[i]
    res <- tryCatch({ td <- load_one(fp); r <- analyze_isometric(td, filename = fp); r$trial_id <- tid; r },
                    error = function(e) NULL)
    if (is.null(res)) next
    ss <- res$step_summary |> dplyr::filter(.data$muscle_side %in% c("left", "right"))
    for (j in seq_len(nrow(ss))) {
      row <- .step_all_methods(res$td, ss[j, ], ss$step_number[j])
      if (!is.null(row)) sens_rows[[length(sens_rows) + 1L]] <- dplyr::mutate(row, fish = fish, trial_id = tid)
    }
  }
}
sens <- dplyr::bind_rows(sens_rows) |> dplyr::filter(is.finite(A_mean_full), is.finite(C_smooth_peak))
cli::cli_alert_success("muscleforcemethodsensitivity: {nrow(sens)} real isometric steps across {dplyr::n_distinct(sens$fish)} fish")

# Ratio-to-A blows up for any step with a near-zero A (dividing by ~noise) --
# restrict the RATIO view (not the scatter, which shows every point) to the
# top 3/4 of steps by |A|, so a handful of barely-active steps don't swamp
# the boxplot with +-huge outliers that are a division artifact, not a real
# method disagreement.
sens_long_all <- sens |>
  tidyr::pivot_longer(c(B_max_raw, C_smooth_peak, D_narrow_mean), names_to = "method", values_to = "force_alt") |>
  dplyr::mutate(method = dplyr::recode(method,
    B_max_raw = "B: MAX raw sample", C_smooth_peak = "C: peak of smoothed trace", D_narrow_mean = "D: narrow-window mean at peak"))
a_floor <- quantile(abs(sens$A_mean_full), 0.25, na.rm = TRUE)
sens_long <- sens_long_all |>
  dplyr::filter(abs(A_mean_full) >= a_floor) |>
  dplyr::mutate(ratio_to_A = force_alt / A_mean_full)

p_scatter <- ggplot(sens_long_all, aes(A_mean_full, force_alt, color = fish)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(alpha = 0.55, size = 1.6) +
  scale_color_manual(values = c(bass16 = "#1b9e77", bass17 = "#d95f02", bass18 = "#7570b3")) +
  facet_wrap(~method, ncol = 3, scales = "free") +
  labs(title = "Does the extraction method change the reported active-minus-passive force? (all real SNR-eligible isometric steps)",
      subtitle = sprintf("n=%d steps, 3 fish | dashed line = y=x (methods agree) | x-axis is the CURRENT method (A, full-window mean) throughout", nrow(sens)),
      x = "A: MEAN full window (current method, N*m)", y = "Alternative method (N*m)") +
  theme_bw(11) + theme(legend.position = "bottom")

p_ratio <- ggplot(sens_long, aes(method, ratio_to_A, fill = fish)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
  scale_fill_manual(values = c(bass16 = "#1b9e77", bass17 = "#d95f02", bass18 = "#7570b3")) +
  labs(title = "Alternative method / current method (A), per real step",
      subtitle = sprintf("ratio = 1 means the two methods agree for that step; systematic departure from 1 = the choice changes SCALE, not just noise | bottom quartile by |A| excluded (near-zero A -> division-artifact outliers), n=%d/%d steps shown", nrow(sens_long) %/% 3L, nrow(sens)),
      x = NULL, y = "force_alternative / force_A (ratio)") +
  theme_bw(11) + theme(legend.position = "bottom")

p_sensitivity <- p_scatter / p_ratio +
  patchwork::plot_annotation(
    title = "Muscle-force extraction method sensitivity (aggregate, real isometric steps)")
ggplot2::ggsave(file.path(OUT_DIR, "muscleforcemethodsensitivity.png"), p_sensitivity, width = 11, height = 11, dpi = 150)
cli::cli_alert_success("Wrote muscleforcemethodsensitivity.png")

med_ratio <- sens_long |> dplyr::group_by(method) |> dplyr::summarise(median_ratio = median(ratio_to_A, na.rm = TRUE), .groups = "drop")
cli::cli_alert_info("Median ratio to current method (A): {paste(sprintf('%s=%.2f', med_ratio$method, med_ratio$median_ratio), collapse=', ')}")
