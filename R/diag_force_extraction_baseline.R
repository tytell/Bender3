# diag_force_extraction_baseline.R
# READ-ONLY DIAGNOSTIC (2026-07-21, PI-requested). Does NOT modify any
# pipeline analysis: sources existing functions and only reads data + writes
# PNGs. Two canon tokens, one script (they share the same representative
# steps/loading code):
#
#   1) forceextractionmethod -- HOW a step's scalar active force gets read
#      off its raw time series. AS-BUILT (2026-07-21) this plot showed the
#      then-current method = MEAN over the full active window. REBUILT
#      2026-07-22 to show what CURRENT production actually does: "Method D"
#      (narrow-window mean of RAW samples centered on the active window's
#      OWN smoothed-trace peak, peak search restricted to the stim duration,
#      duration-guarded) -- .mfv_window_peak_means() in muscle_force_vector.R,
#      .legacy_peak_window_mean() in 03_analyze.R (adopted PI decision
#      2026-07-22). The panel now SHADES the actual averaging window (green,
#      matching muscleforcemethodcompare.png's convention) so the samples
#      Method D averages are visible as a region, not just a horizontal
#      value; the old full-window MEAN is KEPT as a secondary (blue dashed)
#      comparison line, deliberately -- Method D's duration guard FALLS BACK
#      to that plain mean for short bursts (the ~0.05s dynamic L0 bookends),
#      and showing both is exactly what makes that fallback legible in the
#      dynamic panel. The superseded MAX/MIN "hypothetical alternative" line
#      was dropped (it was never production and is covered by
#      muscleforcemethodcompare.png).
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
# 1) forceextractionmethod -- Method D (CURRENT production), showing WHERE the
#    averaged samples are (shaded narrow window), plus the old full-window MEAN
#    as a secondary comparison line
# =============================================================================

# Method D's narrow averaging-window width (s) -- MUST equal the production
# constant so this figure shows exactly what the pipeline does, not an
# approximation (MFV_PEAK_WINDOW_S in muscle_force_vector.R).
PEAK_WINDOW_S <- MFV_PEAK_WINDOW_S

#' Faithful re-implementation of production Method D (.mfv_window_peak_means()
#' in muscle_force_vector.R) on a single RAW torque trace, INCLUDING the
#' 2026-07-22 duration guard: find the smoothed trace's peak within the STIM
#' duration only [0, dur] (excluding the deactivation tail), then average RAW
#' samples in a MFV_PEAK_WINDOW_S-wide window centered on that peak's time
#' (clipped to the full active window). Falls back to the plain full-window
#' mean (Method A) when the stim-duration search window is itself narrower
#' than MFV_PEAK_WINDOW_S -- exactly the case for the ~0.05s dynamic L0
#' bookend bursts, so this panel shows the REAL production behavior (Method D
#' degrades to Method A there) rather than pretending a narrow window fits.
.method_d_on_trace <- function(t_rel, torque_raw, dur, win_end, peak_window_s = PEAK_WINDOW_S) {
  active <- t_rel >= 0 & t_rel <= win_end
  search <- t_rel >= 0 & t_rel <= dur
  full_mean <- mean(torque_raw[active], na.rm = TRUE)
  n_samples <- sum(active, na.rm = TRUE)
  t_search <- t_rel[search]; v_search <- torque_raw[search]
  # duration guard (mirrors .mfv_window_peak_means): a search window shorter
  # than the narrow averaging width can't confine the average to itself.
  fell_back <- !(sum(is.finite(v_search)) >= 3L &&
                   is.finite(diff(range(t_search, na.rm = TRUE))) &&
                   diff(range(t_search, na.rm = TRUE)) >= peak_window_s)
  if (fell_back) {
    return(list(d_val = full_mean, mean_full = full_mean, t_peak_rel = NA_real_, peak_val = NA_real_,
                narrow_lo = NA_real_, narrow_hi = NA_real_, fell_back = TRUE, n_samples = n_samples))
  }
  v_smooth <- .smooth_trace_display_only(v_search)
  smax <- max(v_smooth, na.rm = TRUE); smin <- min(v_smooth, na.rm = TRUE)
  use_max <- abs(smax - full_mean) >= abs(smin - full_mean)
  peak_idx <- if (use_max) which.max(v_smooth) else which.min(v_smooth)
  t_peak <- t_search[peak_idx]
  peak_val <- v_smooth[peak_idx]  # value ON the smoothed trace at the detected peak (dot marker)
  narrow <- active & t_rel >= (t_peak - peak_window_s / 2) & t_rel <= (t_peak + peak_window_s / 2)
  d_val <- mean(torque_raw[narrow], na.rm = TRUE)
  list(d_val = d_val, mean_full = full_mean, t_peak_rel = t_peak, peak_val = peak_val,
       narrow_lo = t_peak - peak_window_s / 2, narrow_hi = t_peak + peak_window_s / 2,
       fell_back = FALSE, n_samples = n_samples)
}

#' One panel: raw torque with the active window shaded (orange), Method D's
#' narrow averaging window shaded (green) so the PI can SEE which samples are
#' averaged, the Method D value as the PRIMARY (green solid) line, and the old
#' full-window MEAN (Method A) as a SECONDARY (blue dashed) comparison line.
#' The full-window mean is KEPT (not dropped) because Method D's duration
#' guard FALLS BACK to it for short bursts (dynamic bookends), so showing both
#' is what makes the fallback legible -- see this panel's dynamic row.
.build_extraction_panel <- function(td, s, title, step_number = NULL) {
  dur <- s$stim_t1_s - s$stim_t0_s
  win_end <- dur + DEACTIVATION_WINDOW_S
  raw <- .raw_window(td, s$stim_t0_s - PLOT_PAD_S, s$stim_t1_s + DEACTIVATION_WINDOW_S + PLOT_PAD_S, s$stim_t0_s, step_number)
  # Method D and Method A are both computed on RAW (unsmoothed) torque,
  # matching exactly what the pipeline averages -- the display smoothing
  # (DISPLAY_SMOOTH_HZ low-pass) only makes the trace legible and only feeds
  # Method D's PEAK-TIME search, never the averaged value itself.
  m <- .method_d_on_trace(raw$t_rel, raw$torque_Nm, dur, win_end)
  raw$torque_smooth_Nm <- .smooth_trace_display_only(raw$torque_Nm)

  p <- ggplot(raw, aes(t_rel, torque_Nm)) +
    annotate("rect", xmin = 0, xmax = win_end, ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.12)
  if (!m$fell_back) {
    p <- p + annotate("rect", xmin = m$narrow_lo, xmax = m$narrow_hi, ymin = -Inf, ymax = Inf,
                      fill = "#059669", alpha = 0.22)
  }
  p <- p +
    geom_vline(xintercept = c(0, win_end), linetype = "dashed", color = "grey40", linewidth = 0.4) +
    geom_line(color = "grey75", linewidth = 0.3) +
    geom_line(aes(y = torque_smooth_Nm), color = "grey15", linewidth = 0.7) +
    geom_hline(aes(yintercept = m$mean_full, colour = "A: MEAN full window (old, secondary)"),
               linewidth = 0.8, linetype = "dashed") +
    geom_hline(aes(yintercept = m$d_val, colour = "D: narrow-window mean at smoothed peak (CURRENT)"),
               linewidth = 1.1)
  if (!m$fell_back) {
    # dot on the smoothed trace at the DETECTED peak -- this is what anchors
    # the green averaging window (peak SEARCH restricted to [0, stim duration]).
    p <- p + annotate("point", x = m$t_peak_rel, y = m$peak_val, colour = "#059669",
                      size = 2.6, shape = 18) +
      annotate("segment", x = m$t_peak_rel, xend = m$t_peak_rel, y = -Inf, yend = m$peak_val,
               colour = "#059669", linewidth = 0.4, linetype = "dotted")
  }
  p <- p +
    scale_colour_manual(name = NULL, values = c(
      "A: MEAN full window (old, secondary)" = "#1d4ed8",
      "D: narrow-window mean at smoothed peak (CURRENT)" = "#059669"))
  d_note <- if (m$fell_back) {
    sprintf("stim duration %.2fs < %.2fs narrow window -> DURATION GUARD fires, Method D FALLS BACK to the full-window mean (D == A = %.4g); no green window drawn because none is used here",
            dur, PEAK_WINDOW_S, m$d_val)
  } else {
    sprintf("green diamond = detected smoothed-trace peak (search restricted to [0, stim dur]) at t=%.2fs; green shade = Method D's %.2fs averaging window centered there (D = %.4g); blue dashed = old full-window MEAN (A = %.4g). Production runs this PER-CHANNEL on all 6 axes; this single torque trace is illustrative.",
            m$t_peak_rel, PEAK_WINDOW_S, m$d_val, m$mean_full)
  }
  p + labs(title = title,
        subtitle = wrap_sub(sprintf("active window = %.2fs (n = %d samples) | pale = raw, dark = display-smoothed only | %s | this active reading is NOT a baseline (see bass16_passivebaselinemethod.png)",
                           win_end, m$n_samples, d_note)),
        x = "Time relative to stim onset (s)", y = "Inertia-corrected torque (N*m)") +
    theme_bw(base_size = 10.5) +
    theme(plot.subtitle = element_text(size = 8), legend.position = "bottom")
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
    title = "Force extraction method (CURRENT = Method D): narrow-window mean at the smoothed peak",
    subtitle = wrap_sub("bass16, one representative step per trial type | orange = active window [stim_t0, stim_t1 + deactivation_window]; green = the samples Method D actually averages | this is what feeds every muscle_force_vector_N / muscle_force_Nm scalar in the pipeline (2026-07-22 switch, replacing the old full-window mean shown in blue dashed for comparison)", w = 110))

ggplot2::ggsave(file.path(OUT_DIR, "bass16_forceextractionmethod.png"), p_extraction,
                width = 9, height = 11, dpi = 150)
cli::cli_alert_success("Wrote bass16_forceextractionmethod.png")


# =============================================================================
# 2) passivebaselinemethod -- STATIC vs INTERPOLATED (segmented) vs
#    runtime-shrunk (dynamic bookend) baseline window mechanics
# =============================================================================

#' A small zoomed inset over one baseline window [lo_rel, hi_rel] (relative to
#' stim onset): the RAW samples that actually fall inside the window as points,
#' plus the flat level (STATIC or POST mean) that the pipeline fits to them, so
#' the PI can judge BY EYE whether that flat line is a good fit to the real
#' passive scatter -- impossible to see on the full multi-second panel, where
#' a ~0.2-0.4s window is a hairline. PI note: it is EXPECTED for the raw
#' passive data to look scattered/fragmented here; that is the point, not a
#' defect to smooth over.
.baseline_window_inset <- function(raw, lo_rel, hi_rel, level_val, level_col, label) {
  d <- raw[is.finite(raw$t_rel) & raw$t_rel >= lo_rel & raw$t_rel <= hi_rel, , drop = FALSE]
  ggplot(d, aes(.data$t_rel, .data$torque_Nm)) +
    geom_point(color = "grey30", size = 0.45, alpha = 0.7) +
    geom_hline(yintercept = level_val, color = level_col, linewidth = 0.9) +
    labs(title = label, x = NULL, y = NULL) +
    theme_bw(base_size = 7) +
    theme(plot.title = element_text(size = 7, face = "bold"),
          axis.text = element_text(size = 6),
          plot.background = element_rect(fill = "white", color = "grey55"),
          panel.grid.minor = element_blank())
}

#' Segmented (isometric/isovelocity) panel: full trace from pre-baseline
#' through post-baseline, both windows shaded, STATIC (flat, pre-only) and
#' INTERPOLATED (pre->post linear, evaluated at the active-window midpoint)
#' baseline levels both drawn for direct comparison.
#'
#' The raw samples INSIDE each baseline window are ALSO overlaid as bold,
#' distinct-colored points (pre = dark red, post = purple) drawn ON TOP of the
#' smoothed line, plus a zoomed inset per window -- so the fit of the flat
#' baseline to the real passive scatter is judgeable regardless of the
#' compressed x-scale (added 2026-07-22: the pale grey75 raw trace WAS already
#' plotted underneath, but a ~0.2-0.4s window on a multi-second axis made the
#' in-window scatter effectively invisible, masked by the smoothed line on top).
.build_baseline_panel_segmented <- function(td, s, title, step_number) {
  dur <- s$stim_t1_s - s$stim_t0_s
  act_mid <- (s$stim_t0_s + (s$stim_t1_s + DEACTIVATION_WINDOW_S)) / 2
  t_lo <- min(s$t_pre_baseline_start_s, s$stim_t0_s) - PLOT_PAD_S
  t_hi <- max(s$t_post_baseline_end_s, s$stim_t1_s + DEACTIVATION_WINDOW_S, na.rm = TRUE) + PLOT_PAD_S
  raw <- .raw_window(td, t_lo, t_hi, s$stim_t0_s, step_number)

  pre_lo_rel <- s$t_pre_baseline_start_s - s$stim_t0_s
  pre_hi_rel <- s$t_pre_baseline_end_s - s$stim_t0_s
  pre_win  <- raw$t_rel >= pre_lo_rel & raw$t_rel <= pre_hi_rel
  post_ok  <- is.finite(s$t_post_baseline_start_s) && is.finite(s$t_post_baseline_end_s)
  pre_mean  <- mean(raw$torque_Nm[pre_win], na.rm = TRUE)
  post_lo_rel <- post_hi_rel <- NA_real_
  post_win <- rep(FALSE, nrow(raw))
  post_mean <- if (post_ok) {
    post_lo_rel <- s$t_post_baseline_start_s - s$stim_t0_s
    post_hi_rel <- s$t_post_baseline_end_s - s$stim_t0_s
    post_win <- raw$t_rel >= post_lo_rel & raw$t_rel <= post_hi_rel
    mean(raw$torque_Nm[post_win], na.rm = TRUE)
  } else NA_real_
  t_pre_mid_rel  <- (s$t_pre_baseline_start_s + s$t_pre_baseline_end_s) / 2 - s$stim_t0_s
  t_post_mid_rel <- if (post_ok) (s$t_post_baseline_start_s + s$t_post_baseline_end_s) / 2 - s$stim_t0_s else NA_real_
  act_mid_rel <- act_mid - s$stim_t0_s
  interp_at_act <- if (is.finite(post_mean) && t_post_mid_rel != t_pre_mid_rel) {
    pre_mean + (post_mean - pre_mean) * (act_mid_rel - t_pre_mid_rel) / (t_post_mid_rel - t_pre_mid_rel)
  } else pre_mean
  raw$torque_smooth_Nm <- .smooth_trace_display_only(raw$torque_Nm)  # display only, see forceextractionmethod panel

  pre_pts  <- raw[pre_win, , drop = FALSE]
  post_pts <- raw[post_win, , drop = FALSE]

  p <- ggplot(raw, aes(t_rel, torque_Nm)) +
    annotate("rect", xmin = pre_lo_rel, xmax = pre_hi_rel,
             ymin = -Inf, ymax = Inf, fill = "#1d4ed8", alpha = 0.12) +
    annotate("rect", xmin = 0, xmax = dur, ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.10) +
    geom_line(color = "grey75", linewidth = 0.3) +
    geom_line(aes(y = torque_smooth_Nm), color = "grey15", linewidth = 0.7) +
    # bold raw samples INSIDE the pre-baseline window -- pop regardless of zoom
    geom_point(data = pre_pts, aes(t_rel, torque_Nm), color = "#dc2626", size = 0.6, alpha = 0.75) +
    geom_hline(yintercept = pre_mean, color = "#1d4ed8", linewidth = 0.9) +
    geom_segment(x = t_pre_mid_rel, xend = act_mid_rel, y = pre_mean, yend = interp_at_act,
                color = "#7c3aed", linewidth = 0.9, linetype = "dashed") +
    geom_point(x = act_mid_rel, y = interp_at_act, color = "#7c3aed", size = 2.2)
  if (post_ok) {
    p <- p +
      annotate("rect", xmin = post_lo_rel, xmax = post_hi_rel,
               ymin = -Inf, ymax = Inf, fill = "#7c3aed", alpha = 0.12) +
      geom_point(data = post_pts, aes(t_rel, torque_Nm), color = "#9333ea", size = 0.6, alpha = 0.75) +
      geom_hline(yintercept = post_mean, color = "#7c3aed", linewidth = 0.5, linetype = "dotted") +
      geom_segment(x = act_mid_rel, xend = t_post_mid_rel, y = interp_at_act, yend = post_mean,
                  color = "#7c3aed", linewidth = 0.9, linetype = "dashed")
  }
  p <- p + labs(title = title,
          subtitle = wrap_sub(sprintf(
            "pale = raw, dark = display-smoothed | RED points = raw samples inside PRE window, PURPLE points = raw inside POST window (insets zoom each) | blue line = STATIC baseline (%.4g, pre-only) | purple dashed/dot = INTERPOLATED baseline at active midpoint (%.4g)%s",
            pre_mean, interp_at_act, if (!post_ok) " -- NO post-baseline available, interpolated == static here" else "")),
          x = "Time relative to stim onset (s)", y = "Inertia-corrected torque (N*m)") +
    theme_bw(base_size = 10.5) + theme(plot.subtitle = element_text(size = 7.5))

  # Zoomed inset(s): raw scatter vs the fitted flat level, per baseline window.
  pre_inset <- .baseline_window_inset(raw, pre_lo_rel, pre_hi_rel, pre_mean, "#1d4ed8",
                                      "PRE window (raw vs STATIC)")
  p <- p + patchwork::inset_element(pre_inset, left = 0.02, bottom = 0.62, right = 0.30, top = 0.99,
                                    align_to = "panel")
  if (post_ok) {
    post_inset <- .baseline_window_inset(raw, post_lo_rel, post_hi_rel, post_mean, "#7c3aed",
                                         "POST window (raw vs mean)")
    p <- p + patchwork::inset_element(post_inset, left = 0.70, bottom = 0.62, right = 0.98, top = 0.99,
                                      align_to = "panel")
  }
  p
}

p_base_iso  <- .build_baseline_panel_segmented(iso_td,  iso_step,  sprintf("Isometric (angle %.1f deg) -- FIXED metadata window", iso_step$operating_point), step_number = iso_step$step_number)
p_base_isov <- .build_baseline_panel_segmented(isov_td, isov_step, sprintf("Isovelocity (%.1f deg/s) -- FIXED metadata window", isov_step$operating_point), step_number = isov_step$step_number)

p_base_dyn <- if (have_dynamic) {
  # step_number = NULL: dyn_td is single_finite, see p_ext_dyn's matching comment above.
  raw <- .raw_window(dyn_td, dyn_step$t_pre_baseline_start_s - PLOT_PAD_S,
                     dyn_step$stim_t1_s + DEACTIVATION_WINDOW_S + PLOT_PAD_S, dyn_step$stim_t0_s, step_number = NULL)
  pre_lo_rel <- dyn_step$t_pre_baseline_start_s - dyn_step$stim_t0_s
  pre_hi_rel <- dyn_step$t_pre_baseline_end_s - dyn_step$stim_t0_s
  base_mean <- mean(dyn_td$torque_inertia_corrected_Nm[
    dyn_td$t.s >= dyn_step$t_pre_baseline_start_s & dyn_td$t.s <= dyn_step$t_pre_baseline_end_s], na.rm = TRUE)
  raw$torque_smooth_Nm <- .smooth_trace_display_only(raw$torque_Nm)  # display only, see forceextractionmethod panel
  pre_pts <- raw[is.finite(raw$t_rel) & raw$t_rel >= pre_lo_rel & raw$t_rel <= pre_hi_rel, , drop = FALSE]
  # Explicit in-panel annotation of the REAL limitation (PI-directed, so it
  # doesn't read as a rendering gap): dynamic single_finite files have no
  # index_step_* table, so there is genuinely no post-stim baseline window to
  # draw -- NOT fabricated to force symmetry with the segmented panels above.
  y_anno <- max(raw$torque_smooth_Nm, na.rm = TRUE)
  no_post_lab <- "NO post-baseline window exists for a dynamic bookend:\nsingle_finite files have no index_step_* table, so there is\nnothing to interpolate to. This is a real limitation, not a\nrendering gap -- a synthetic post-window is NOT fabricated."
  p <- ggplot(raw, aes(t_rel, torque_Nm)) +
    annotate("rect", xmin = pre_lo_rel, xmax = pre_hi_rel,
             ymin = -Inf, ymax = Inf, fill = "#1d4ed8", alpha = 0.12) +
    annotate("rect", xmin = 0, xmax = (dyn_step$stim_t1_s - dyn_step$stim_t0_s),
             ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.10) +
    geom_line(color = "grey75", linewidth = 0.3) +
    geom_line(aes(y = torque_smooth_Nm), color = "grey15", linewidth = 0.7) +
    geom_point(data = pre_pts, aes(t_rel, torque_Nm), color = "#dc2626", size = 0.6, alpha = 0.75) +
    geom_hline(yintercept = base_mean, color = "#1d4ed8", linewidth = 0.9) +
    annotate("label", x = pre_hi_rel + (dyn_step$stim_t1_s - dyn_step$stim_t0_s) * 0.5, y = y_anno,
             label = no_post_lab, hjust = 0, vjust = 1, size = 2.6,
             color = "grey20", fill = "#fff7ed") +
    labs(title = sprintf("Dynamic L0 bookend (%s) -- RUNTIME-SHRUNK window, NO post-baseline (see annotation)", dyn_step$muscle_side),
        subtitle = wrap_sub(sprintf(
          "pale = raw, dark = display-smoothed | RED points = raw samples inside the baseline window (inset zooms) | blue shade = baseline window found by shrinking a candidate window until stim-free (extract_dynamic_l0_bookends.R), width %.2fs (vs. isometric/isovelocity's FIXED metadata width above) | blue line = STATIC baseline (%.4g)",
          dyn_step$t_pre_baseline_end_s - dyn_step$t_pre_baseline_start_s, base_mean)),
        x = "Time relative to stim onset (s)", y = "Inertia-corrected torque (N*m)") +
    theme_bw(base_size = 10.5) + theme(plot.subtitle = element_text(size = 7.5))
  dyn_inset <- .baseline_window_inset(raw, pre_lo_rel, pre_hi_rel, base_mean, "#1d4ed8",
                                      "baseline window (raw vs STATIC)")
  p + patchwork::inset_element(dyn_inset, left = 0.02, bottom = 0.62, right = 0.32, top = 0.99,
                               align_to = "panel")
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
