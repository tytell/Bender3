# diag_force_extraction_baseline.R
# READ-ONLY DIAGNOSTIC (2026-07-21, PI-requested). Does NOT modify any
# pipeline analysis: sources existing functions and only reads data + writes
# PNGs. Two canon tokens, one script (they share the same representative
# steps/loading code):
#
#   1) forceextractionmethod -- HOW a step's scalar active force gets read
#      off its raw time series. AS-BUILT (2026-07-21) this plot shows the
#      then-current method = MEAN over the active window
#      [stim_t0, stim_t1 + deactivation_window]. SUPERSEDED 2026-07-22 (PI
#      decision, following muscleforcemethodcompare.png/
#      muscleforcemethodsensitivity.png below): production now uses
#      "Method D" (narrow-window mean centered on the smoothed trace's own
#      peak, duration-guarded) -- .mfv_window_peak_means() in
#      muscle_force_vector.R, .legacy_peak_window_mean() in 03_analyze.R.
#      This plot itself was NOT regenerated with the new method (it is a
#      historical record of the comparison that motivated the switch, kept
#      as-is); NOT max-in-window, NOT last-N-samples -- it shows what the
#      OLD plain-mean choice cost by also marking where MAX-in-window would
#      fall, and reports n_samples actually averaged.
#   2) passivebaselinemethod -- HOW the passive reference subtracted from
#      that active window is calculated, and how the mechanism differs by
#      trial type:
#        - isometric/isovelocity (segmented_finite): a FIXED window read
#          directly from protocol metadata (index_step_t_pre_baseline_*),
#          set at acquisition time -- STATIC (own-step pre-stim window
#          only) vs INTERPOLATED (linear between the pre- AND post-stim
#          window means, 03_analyze.R's passive_force_Nm_interp) are two
#          alternative READS of that same fixed window pair.
#        - dynamic bookends (single_finite): NO protocol-provided baseline
#          window exists at all (dynamic files have no index_step_* table).
#          extract_dynamic_l0_bookends.R computes one at RUNTIME by
#          shrinking a candidate window until no sample in it carries stim
#          on either side -- a fundamentally different mechanism, and it
#          has no post-baseline / interpolated counterpart.
#
# One representative fish (bass16 -- clean signal, avoids bass17's known
# low-SNR confound) illustrates the MECHANISM; this is about the METHOD,
# not a per-fish result, so it does not need all three fish (matches
# diag_within_step_force_development.R's single-fish precedent).
#
# Run: Rscript R/diag_force_extraction_baseline.R

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(ggplot2)
  library(cli); library(patchwork)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")

SRC_DIR <- raw_source_dir(BASS16_RAW_SUBFOLDER)
OUT_DIR <- FIGS_DIAGNOSTIC_DIR

src("00_load_bender_flat.R")
src("01_calibrate.R")
src("02_deconvolve.R")
src("muscle_geometry.R")
src("fit_fv_fl.R")
src("03_analyze.R")
src("parse_trial_filename.R")
src("plot_force_vs_time.R")  # .detect_stim_events() -- extract_dynamic_l0_bookends.R depends on this
src("muscle_force_vector.R")
src("extract_dynamic_l0_bookends.R")

DEACTIVATION_WINDOW_S <- MFV_DEACTIVATION_WINDOW_S  # 0.5s -- matches the real pipeline default exactly
PLOT_PAD_S <- 0.3  # extra display-only padding around the windows of interest
wrap_sub <- function(s, w = 100) paste(strwrap(s, width = w), collapse = "\n")

load_one <- function(f) {
  td <- load_bender_flat(f, do_filter = TRUE, loadtorques = "x")
  tau <- deconvolve_bender(f, hub_path = NULL, verbose = FALSE)
  N <- min(nrow(td), length(tau)); td <- td[seq_len(N), , drop = FALSE]
  td$torque_inertia_corrected_Nm <- tau[seq_len(N)]
  attr(td, "Filename") <- f
  td
}

manifest <- parse_trial_directory(SRC_DIR)

# ------------------------------------------------------------ isometric ---
iso_f  <- manifest$fullpath[manifest$protocol == "isometric"][1]
iso_td <- load_one(iso_f)
iso_ss <- analyze_isometric(iso_td, filename = iso_f)$step_summary
# Step 4: largest-angle rep of block 1 (left, ~17 deg) -- skip step 1 = L0,
# less illustrative of a genuine active pull against a bend.
iso_step <- iso_ss[iso_ss$step_number == 4, ][1, ]

# ---------------------------------------------------------- isovelocity ---
isov_f  <- manifest$fullpath[manifest$protocol == "isovelocity"][1]
isov_td <- load_one(isov_f)
isov_ss <- analyze_isovelocity(isov_td, filename = isov_f)$step_summary
isov_step <- isov_ss[isov_ss$step_number == 3, ][1, ]  # a genuine nonzero-velocity rep

# --------------------------------------------------- dynamic L0 bookend ---
dyn_f  <- manifest$fullpath[manifest$protocol == "dynamic"][1]
dyn_td <- load_one(dyn_f)
bookends <- detect_dynamic_l0_bookends(dyn_td)
if (nrow(bookends) == 0L) {
  # Try a few more dynamic files before giving up -- not every dynamic file
  # necessarily has a clean detectable bookend burst.
  more <- manifest$fullpath[manifest$protocol == "dynamic"][2:5]
  for (f in more) {
    td_try <- load_one(f)
    be_try <- detect_dynamic_l0_bookends(td_try)
    if (nrow(be_try) > 0L) { dyn_f <- f; dyn_td <- td_try; bookends <- be_try; break }
  }
}
have_dynamic <- nrow(bookends) > 0L
if (have_dynamic) {
  dyn_step <- bookends[1, ]
  cli::cli_alert_success("Dynamic bookend found in {basename(dyn_f)}: {dyn_step$muscle_side} at t={round(dyn_step$stim_t0_s,2)}s")
} else {
  cli::cli_alert_warning("No dynamic L0 bookend found in the first 5 dynamic files -- dynamic panel will be a placeholder")
}

cli::cli_alert_info("Isometric step {iso_step$step_number} ({basename(iso_f)}): angle {round(iso_step$operating_point,1)} deg")
cli::cli_alert_info("Isovelocity step {isov_step$step_number} ({basename(isov_f)}): rate {round(isov_step$operating_point,1)} deg/s")


# =============================================================================
# Shared raw-window reader (torque, plain -- NOT baseline-subtracted, so both
# diagnostics can show exactly what's fed INTO the extraction/baseline math)
# =============================================================================

#' @param td Loaded trial tibble (torque_inertia_corrected_Nm present).
#' @param t0,t1 Time window to slice (s) -- typically pre_baseline start
#'   through post_baseline end (or stim_t1+relax if no post-baseline).
#' @param anchor_t Time subtracted to make t_rel (usually stim_t0_s).
#' @param step_number REQUIRED for segmented (isometric/isovelocity) `td`:
#'   `t.s` RESETS TO 0 AT THE START OF EVERY STEP in a segmented_finite file
#'   (confirmed empirically 2026-07-21 -- NOT a global file-wide clock), so
#'   filtering on t.s alone silently picks up matching-time-range rows from
#'   EVERY step in the file, not just the intended one -- this was the ACTUAL
#'   cause of "solid black blob" traces in bass16_forceextractionmethod.png /
#'   bass16_passivebaselinemethod.png's isometric+isovelocity panels (PI
#'   correctly doubted the original "dense raw samples" explanation -- the
#'   real bug was ~16x too many rows from cross-step contamination, not
#'   genuine within-step sample density). Pass NULL for dynamic-bookend `td`,
#'   which has a single continuous monotonic t.s and no step_number column at
#'   all (verified empirically) -- no contamination is possible there, which
#'   is exactly why that panel already looked clean.
.raw_window <- function(td, t0, t1, anchor_t, step_number) {
  keep <- td$t.s >= t0 & td$t.s <= t1
  if (!is.null(step_number)) {
    if (!"step_number" %in% names(td)) {
      cli::cli_abort(".raw_window: step_number requested but td has no step_number column")
    }
    keep <- keep & td$step_number == step_number
  }
  tibble::tibble(t_rel = td$t.s[keep] - anchor_t, torque_Nm = td$torque_inertia_corrected_Nm[keep])
}


# =============================================================================
# 1) forceextractionmethod -- MEAN vs MAX in the active window, n_samples
# =============================================================================

#' One panel: raw torque, active window shaded, MEAN (solid, current method)
#' and MAX (dashed, hypothetical alternative) marked, n_samples annotated.
.build_extraction_panel <- function(td, s, title, step_number = NULL) {
  dur <- s$stim_t1_s - s$stim_t0_s
  win_end <- dur + DEACTIVATION_WINDOW_S
  raw <- .raw_window(td, s$stim_t0_s - PLOT_PAD_S, s$stim_t1_s + DEACTIVATION_WINDOW_S + PLOT_PAD_S, s$stim_t0_s, step_number)
  active <- raw$t_rel >= 0 & raw$t_rel <= win_end
  n_samples <- sum(active, na.rm = TRUE)
  mean_val <- mean(raw$torque_Nm[active], na.rm = TRUE)
  max_val  <- max(raw$torque_Nm[active], na.rm = TRUE)
  min_val  <- min(raw$torque_Nm[active], na.rm = TRUE)
  # whichever of max/min is farther from the mean is the more relevant
  # "extreme" comparison point (this force can be signed either way).
  extreme_val <- if (abs(max_val - mean_val) >= abs(min_val - mean_val)) max_val else min_val
  extreme_lab <- if (identical(extreme_val, max_val)) "MAX in window" else "MIN in window"
  # MEAN/MAX/MIN annotations are always computed on RAW (unsmoothed) torque
  # above, matching exactly what the pipeline itself averages -- this display
  # smoothing (DISPLAY_SMOOTH_HZ low-pass, plot_force_vs_time.R) is ONLY to
  # make the trace legible; it never feeds the annotated numbers.
  raw$torque_smooth_Nm <- .smooth_trace_display_only(raw$torque_Nm)

  ggplot(raw, aes(t_rel, torque_Nm)) +
    annotate("rect", xmin = 0, xmax = win_end, ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.12) +
    geom_vline(xintercept = c(0, win_end), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_line(color = "grey75", linewidth = 0.3) +
    geom_line(aes(y = torque_smooth_Nm), color = "grey15", linewidth = 0.7) +
    geom_hline(yintercept = mean_val, color = "#1d4ed8", linewidth = 1.0) +
    geom_hline(yintercept = extreme_val, color = "#b91c1c", linewidth = 0.8, linetype = "dotted") +
    labs(title = title,
        subtitle = wrap_sub(sprintf("active window = %.2fs (n = %d samples) | pale = raw, dark = display-smoothed only | blue solid = MEAN of the ACTIVE window (%.4g, CURRENT method -- NOT a baseline; see bass16_passivebaselinemethod.png for the separate pre-stim baseline mean) | red dotted = %s (%.4g)",
                           win_end, n_samples, mean_val, extreme_lab, extreme_val)),
        x = "Time relative to stim onset (s)", y = "Inertia-corrected torque (N*m)") +
    theme_bw(base_size = 10.5) + theme(plot.subtitle = element_text(size = 8))
}

p_ext_iso  <- .build_extraction_panel(iso_td,  iso_step,  sprintf("Isometric (angle %.1f deg)", iso_step$operating_point), step_number = iso_step$step_number)
p_ext_isov <- .build_extraction_panel(isov_td, isov_step, sprintf("Isovelocity (%.1f deg/s)", isov_step$operating_point), step_number = isov_step$step_number)

p_ext_dyn <- if (have_dynamic) {
  # Bookends have no index_step_* stim_t1_s -- use the detected burst's own
  # [stim_t0_s, stim_t1_s] directly (same active-window formula, different
  # source of stim_t0/stim_t1). step_number = NULL: dyn_td is single_finite
  # (one continuous monotonic t.s, no step_number column at all -- see
  # .raw_window()'s docstring), so no cross-step contamination is possible.
  fake_s <- tibble::tibble(stim_t0_s = dyn_step$stim_t0_s, stim_t1_s = dyn_step$stim_t1_s)
  .build_extraction_panel(dyn_td, fake_s, sprintf("Dynamic L0 bookend (%s)", dyn_step$muscle_side), step_number = NULL)
} else {
  ggplot() + annotate("text", x = 0, y = 0, label = "No dynamic bookend detected in sampled files") +
    theme_void()
}

p_extraction <- p_ext_iso / p_ext_isov / p_ext_dyn +
  patchwork::plot_annotation(
    title = "Force extraction method: MEAN over the active window (current) vs. MAX/MIN (hypothetical alternative)",
    subtitle = wrap_sub("bass16, one representative step per trial type | orange = active window [stim_t0, stim_t1 + deactivation_window] | this is what feeds every muscle_force_vector_N / muscle_force_Nm scalar in the pipeline", w = 110))

ggplot2::ggsave(file.path(OUT_DIR, "bass16_forceextractionmethod.png"), p_extraction,
                width = 9, height = 11, dpi = 150)
cli::cli_alert_success("Wrote bass16_forceextractionmethod.png")


# =============================================================================
# 2) passivebaselinemethod -- STATIC vs INTERPOLATED (segmented) vs
#    runtime-shrunk (dynamic bookend) baseline window mechanics
# =============================================================================

#' Segmented (isometric/isovelocity) panel: full trace from pre-baseline
#' through post-baseline, both windows shaded, STATIC (flat, pre-only) and
#' INTERPOLATED (pre->post linear, evaluated at the active-window midpoint)
#' baseline levels both drawn for direct comparison.
.build_baseline_panel_segmented <- function(td, s, title, step_number) {
  dur <- s$stim_t1_s - s$stim_t0_s
  act_mid <- (s$stim_t0_s + (s$stim_t1_s + DEACTIVATION_WINDOW_S)) / 2
  t_lo <- min(s$t_pre_baseline_start_s, s$stim_t0_s) - PLOT_PAD_S
  t_hi <- max(s$t_post_baseline_end_s, s$stim_t1_s + DEACTIVATION_WINDOW_S, na.rm = TRUE) + PLOT_PAD_S
  raw <- .raw_window(td, t_lo, t_hi, s$stim_t0_s, step_number)

  pre_win  <- raw$t_rel >= (s$t_pre_baseline_start_s - s$stim_t0_s) & raw$t_rel <= (s$t_pre_baseline_end_s - s$stim_t0_s)
  post_ok  <- is.finite(s$t_post_baseline_start_s) && is.finite(s$t_post_baseline_end_s)
  pre_mean  <- mean(raw$torque_Nm[pre_win], na.rm = TRUE)
  post_mean <- if (post_ok) {
    post_win <- raw$t_rel >= (s$t_post_baseline_start_s - s$stim_t0_s) & raw$t_rel <= (s$t_post_baseline_end_s - s$stim_t0_s)
    mean(raw$torque_Nm[post_win], na.rm = TRUE)
  } else NA_real_
  t_pre_mid_rel  <- (s$t_pre_baseline_start_s + s$t_pre_baseline_end_s) / 2 - s$stim_t0_s
  t_post_mid_rel <- if (post_ok) (s$t_post_baseline_start_s + s$t_post_baseline_end_s) / 2 - s$stim_t0_s else NA_real_
  act_mid_rel <- act_mid - s$stim_t0_s
  interp_at_act <- if (is.finite(post_mean) && t_post_mid_rel != t_pre_mid_rel) {
    pre_mean + (post_mean - pre_mean) * (act_mid_rel - t_pre_mid_rel) / (t_post_mid_rel - t_pre_mid_rel)
  } else pre_mean
  raw$torque_smooth_Nm <- .smooth_trace_display_only(raw$torque_Nm)  # display only, see forceextractionmethod panel

  p <- ggplot(raw, aes(t_rel, torque_Nm)) +
    annotate("rect", xmin = s$t_pre_baseline_start_s - s$stim_t0_s, xmax = s$t_pre_baseline_end_s - s$stim_t0_s,
             ymin = -Inf, ymax = Inf, fill = "#1d4ed8", alpha = 0.12) +
    annotate("rect", xmin = 0, xmax = dur, ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.10) +
    geom_line(color = "grey75", linewidth = 0.3) +
    geom_line(aes(y = torque_smooth_Nm), color = "grey15", linewidth = 0.7) +
    geom_hline(yintercept = pre_mean, color = "#1d4ed8", linewidth = 0.9) +
    geom_segment(x = t_pre_mid_rel, xend = act_mid_rel, y = pre_mean, yend = interp_at_act,
                color = "#7c3aed", linewidth = 0.9, linetype = "dashed") +
    geom_point(x = act_mid_rel, y = interp_at_act, color = "#7c3aed", size = 2.2)
  if (post_ok) {
    p <- p +
      annotate("rect", xmin = s$t_post_baseline_start_s - s$stim_t0_s, xmax = s$t_post_baseline_end_s - s$stim_t0_s,
               ymin = -Inf, ymax = Inf, fill = "#7c3aed", alpha = 0.12) +
      geom_hline(yintercept = post_mean, color = "#7c3aed", linewidth = 0.5, linetype = "dotted") +
      geom_segment(x = act_mid_rel, xend = t_post_mid_rel, y = interp_at_act, yend = post_mean,
                  color = "#7c3aed", linewidth = 0.9, linetype = "dashed")
  }
  p + labs(title = title,
          subtitle = wrap_sub(sprintf(
            "pale = raw, dark = display-smoothed only | blue shade = PRE-baseline window (fixed, from protocol metadata) | purple shade = POST-baseline window | orange = stim | blue line = STATIC baseline (%.4g, pre-only) | purple dashed/dot = INTERPOLATED baseline at active-window midpoint (%.4g)%s",
            pre_mean, interp_at_act, if (!post_ok) " -- NO post-baseline available, interpolated == static here" else "")),
          x = "Time relative to stim onset (s)", y = "Inertia-corrected torque (N*m)") +
    theme_bw(base_size = 10.5) + theme(plot.subtitle = element_text(size = 7.5))
}

p_base_iso  <- .build_baseline_panel_segmented(iso_td,  iso_step,  sprintf("Isometric (angle %.1f deg) -- FIXED metadata window", iso_step$operating_point), step_number = iso_step$step_number)
p_base_isov <- .build_baseline_panel_segmented(isov_td, isov_step, sprintf("Isovelocity (%.1f deg/s) -- FIXED metadata window", isov_step$operating_point), step_number = isov_step$step_number)

p_base_dyn <- if (have_dynamic) {
  # step_number = NULL: dyn_td is single_finite, see p_ext_dyn's matching comment above.
  raw <- .raw_window(dyn_td, dyn_step$t_pre_baseline_start_s - PLOT_PAD_S,
                     dyn_step$stim_t1_s + DEACTIVATION_WINDOW_S + PLOT_PAD_S, dyn_step$stim_t0_s, step_number = NULL)
  base_mean <- mean(dyn_td$torque_inertia_corrected_Nm[
    dyn_td$t.s >= dyn_step$t_pre_baseline_start_s & dyn_td$t.s <= dyn_step$t_pre_baseline_end_s], na.rm = TRUE)
  raw$torque_smooth_Nm <- .smooth_trace_display_only(raw$torque_Nm)  # display only, see forceextractionmethod panel
  ggplot(raw, aes(t_rel, torque_Nm)) +
    annotate("rect", xmin = dyn_step$t_pre_baseline_start_s - dyn_step$stim_t0_s,
             xmax = dyn_step$t_pre_baseline_end_s - dyn_step$stim_t0_s,
             ymin = -Inf, ymax = Inf, fill = "#1d4ed8", alpha = 0.12) +
    annotate("rect", xmin = 0, xmax = (dyn_step$stim_t1_s - dyn_step$stim_t0_s),
             ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.10) +
    geom_line(color = "grey75", linewidth = 0.3) +
    geom_line(aes(y = torque_smooth_Nm), color = "grey15", linewidth = 0.7) +
    geom_hline(yintercept = base_mean, color = "#1d4ed8", linewidth = 0.9) +
    labs(title = sprintf("Dynamic L0 bookend (%s) -- RUNTIME-SHRUNK window, no fixed metadata", dyn_step$muscle_side),
        subtitle = wrap_sub(sprintf(
          "pale = raw, dark = display-smoothed only | blue shade = baseline window found by shrinking a candidate window until stim-free (extract_dynamic_l0_bookends.R) -- window width %.2fs (vs. isometric/isovelocity's FIXED metadata width above) | blue line = STATIC baseline (%.4g) -- NO post-baseline / interpolated counterpart exists for a dynamic bookend",
          dyn_step$t_pre_baseline_end_s - dyn_step$t_pre_baseline_start_s, base_mean)),
        x = "Time relative to stim onset (s)", y = "Inertia-corrected torque (N*m)") +
    theme_bw(base_size = 10.5) + theme(plot.subtitle = element_text(size = 7.5))
} else {
  ggplot() + annotate("text", x = 0, y = 0, label = "No dynamic bookend detected in sampled files") + theme_void()
}

p_baseline <- p_base_iso / p_base_isov / p_base_dyn +
  patchwork::plot_annotation(
    title = "Passive baseline method: pre-stim reference window + STATIC vs. INTERPOLATED, by trial type",
    subtitle = wrap_sub("bass16, one representative step per trial type | segmented (isometric/isovelocity) uses a FIXED protocol-defined window; dynamic bookends compute one at runtime with NO post-baseline", w = 110))

ggplot2::ggsave(file.path(OUT_DIR, "bass16_passivebaselinemethod.png"), p_baseline,
                width = 10, height = 12, dpi = 150)
cli::cli_alert_success("Wrote bass16_passivebaselinemethod.png")
