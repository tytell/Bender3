# plot_fatigue_timeline.R
# Diagnostic + individual-tier (dual-tagged): L0 isometric muscle force vs.
# REAL elapsed session time (minutes since the fish's first stimulated step),
# across ALL isometric trials of one specimen. Distinct from
# plot_fatigue_check.R (force vs. stimulation ORDER within a block, no real
# time axis) -- that plot is unchanged and kept in parallel, not replaced.
#
# WHY: analysis_muscle_force_vector_log.md's point-selection design needs a
# "Gate A" (per-step trustworthiness) fatigue flag based on whole-session
# decay, not just per-step SNR/drift. This plot is step 1 -- gathering real
# decay data across bass16/17/18 -- BEFORE any threshold is chosen. Do not
# treat any hard cutoff derived here as final without a follow-up PI review.
#
# L0 DEFINITION (corrected 2026-07-21, supersedes the peak-force proxy from
# earlier the same day): every isometric trial in this protocol commands
# operating_point == 0 deg (straight/neutral bend) repeatedly, as the FIRST
# step of every recruitment block (verified directly in the raw HDF5
# index_step_operating_point/index_step_recruitment for bass16/17/18 -- e.g.
# bass16 trial 09 hits angle 0 at steps 1, 5, 9, 13, alternating
# left/right_unilateral). These are genuine STIMULATED contractions AT L0,
# deliberately repeated early and late within a trial (and across trials
# through the session) specifically so they can serve as pre-/post-stim
# fatigue bookends -- exactly the "L0 isometric pre- and post-stim
# contractions" this plot is for. L0 is therefore NOT re-derived from
# whichever angle happened to produce peak force among SNR-passing steps
# (that proxy silently EXCLUDED most of the true L0 reps whenever the
# peak-force angle wasn't exactly 0, e.g. bass16 landed on 5.6 deg and
# dropped 3 of its 4 true L0 reps from the "near-L0" set). "Near-L0" = a
# tight epsilon around the commanded 0, not a %-of-range search band --
# there's nothing to search for, the protocol tells us the L0 angle exactly.
#
# REAL ELAPSED TIME: session start = the EARLIEST wall_clock_start among all
# actually-stimulated steps (muscle_side left/right) across the WHOLE
# isometric corpus for this fish (not just the L0 subset), so elapsed time
# reflects true session progress even if the first L0 step comes later.
# wall_clock_start is read once, at its single source point
# (03_analyze.R::.read_segmented_step_geometry), and flows through unchanged
# -- see that file's comment.
#
# SNR-failing points are ALPHA-FLAGGED, never dropped (per-fish/individual-
# tier policy, analysis_muscle_force_vector_log.md).
#
# MISSING DATA (added 2026-07-21): every code path below returns a REAL
# ggplot object, never NULL/no file. When a fish's data can't support the
# real plot (e.g. bass17: zero left/right isometric steps, no valid
# wall_clock_start, or -- pre-fix -- no SNR-passing anchor), a placeholder
# figure is built and saved with the specific reason printed ON the image,
# so the absence of a real fatigue-timeline plot is visible in
# figs_diagnostic/ itself, not just in a console log the PI may not see.

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2); library(cli)
})

#' Blank placeholder figure carrying a build-in explanation of why the real
#' fatigue-timeline plot couldn't be built for this fish.
.fatigue_timeline_placeholder <- function(reason, title) {
  ggplot(tibble::tibble(x = 0, y = 0)) +
    annotate("text", x = 0, y = 0, label = paste0("No fatigue timeline plotted.\n\n", reason),
             size = 4.2, hjust = 0.5, vjust = 0.5, lineheight = 1.15) +
    xlim(-1, 1) + ylim(-1, 1) +
    labs(title = title) +
    theme_void(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, size = 13),
          panel.border = element_rect(color = "grey70", fill = NA))
}

#' @param step_summary_vec Combined vector-force step_summary across ALL
#'   isometric trials of ONE specimen (bind_rows(isometric_steps_vec_all),
#'   as built by run_fv_fl_power_pipeline.R) -- must have trial_id,
#'   muscle_side, operating_point, muscle_force_vector_N, activation_snr,
#'   wall_clock_start.
#' @param l0_angle_epsilon_deg Floating-point tolerance around the commanded
#'   L0 angle (0 deg) -- NOT a search band. All isometric protocols in this
#'   dataset command operating_point == 0 exactly, so this only absorbs
#'   floating-point noise from the gear-ratio/encoder scaling chain, not
#'   genuine angle uncertainty. Default (0.5 deg) is far below the smallest
#'   nonzero step spacing seen across bass16/17/18 (>= ~2.1 deg).
#' @param snr_min Activation-SNR pass threshold (alpha-flag only, never drops).
#' @return list(plot = <ggplot, always non-NULL>, is_placeholder = <logical>,
#'   l0_angle_deg, n_near_l0, session_start) -- the scalars are reported so
#'   the PI can see exactly what was plotted, not just the picture.
build_fatigue_timeline_plot <- function(step_summary_vec,
                                        l0_angle_epsilon_deg = 0.5,
                                        snr_min = 3.0,
                                        title = "Fatigue timeline: L0 isometric force vs. real elapsed time") {
  df <- step_summary_vec |>
    dplyr::filter(.data$muscle_side %in% c("left", "right"))

  if (nrow(df) == 0L) {
    reason <- "No left/right isometric steps in this specimen's data."
    cli::cli_warn("build_fatigue_timeline_plot: {reason}")
    return(list(plot = .fatigue_timeline_placeholder(reason, title), is_placeholder = TRUE,
                l0_angle_deg = NA_real_, n_near_l0 = 0L, session_start = as.POSIXct(NA)))
  }

  # Session start from the FULL stimulated corpus (before narrowing to L0),
  # so elapsed time reflects true session progress.
  wc_all <- suppressWarnings(as.POSIXct(df$wall_clock_start, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC"))
  session_start <- suppressWarnings(min(wc_all, na.rm = TRUE))
  if (!is.finite(session_start)) {
    reason <- "No valid wall_clock_start timestamps on any isometric step."
    cli::cli_warn("build_fatigue_timeline_plot: {reason}")
    return(list(plot = .fatigue_timeline_placeholder(reason, title), is_placeholder = TRUE,
                l0_angle_deg = NA_real_, n_near_l0 = 0L, session_start = as.POSIXct(NA)))
  }

  df <- df |>
    dplyr::filter(is.finite(.data$operating_point), is.finite(.data$muscle_force_vector_N)) |>
    dplyr::mutate(
      wall_clock_start_posix = suppressWarnings(as.POSIXct(.data$wall_clock_start, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC")),
      elapsed_min = as.numeric(difftime(.data$wall_clock_start_posix, session_start, units = "mins")),
      snr_pass = is.finite(.data$activation_snr) & .data$activation_snr >= snr_min
    ) |>
    dplyr::filter(is.finite(.data$elapsed_min))

  if (nrow(df) == 0L) {
    reason <- "No isometric steps with a finite elapsed session time."
    cli::cli_warn("build_fatigue_timeline_plot: {reason}")
    return(list(plot = .fatigue_timeline_placeholder(reason, title), is_placeholder = TRUE,
                l0_angle_deg = NA_real_, n_near_l0 = 0L, session_start = session_start))
  }

  l0_angle_deg <- 0
  near_l0 <- dplyr::filter(df, abs(.data$operating_point - l0_angle_deg) <= l0_angle_epsilon_deg)
  if (nrow(near_l0) == 0L) {
    reason <- sprintf("No isometric steps commanded within %.1f deg of L0 (0 deg) -- unexpected for this protocol template, check the raw index_step_operating_point values.", l0_angle_epsilon_deg)
    cli::cli_warn("build_fatigue_timeline_plot: {reason}")
    return(list(plot = .fatigue_timeline_placeholder(reason, title), is_placeholder = TRUE,
                l0_angle_deg = l0_angle_deg, n_near_l0 = 0L, session_start = session_start))
  }

  p <- ggplot(near_l0, aes(x = .data$elapsed_min, y = .data$muscle_force_vector_N,
                           color = .data$trial_id, alpha = .data$snr_pass)) +
    geom_point(size = 2.4) +
    geom_smooth(aes(group = 1), method = "loess", formula = y ~ x, se = FALSE,
               color = "black", linewidth = 0.6, span = 0.9, na.rm = TRUE) +
    scale_alpha_manual(values = c(`TRUE` = 1.0, `FALSE` = 0.3), name = paste0("SNR >= ", snr_min),
                       labels = c(`TRUE` = "pass", `FALSE` = "fail (shown, not dropped)")) +
    labs(title = title,
        subtitle = sprintf(
          "L0 = 0 deg (commanded neutral bend, repeated pre-/post-stim per protocol block); n = %d contraction(s)",
          nrow(near_l0)),
        x = "Elapsed session time (min, since first stimulated step)",
        y = "Muscle force along u_hat (N)", color = "Trial") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")

  list(plot = p, is_placeholder = FALSE, l0_angle_deg = l0_angle_deg,
       n_near_l0 = nrow(near_l0), session_start = session_start)
}
