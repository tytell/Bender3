# plot_force_vs_time.R
# Muscle force vs. time, pooled across trials, one figure per protocol
# category (isometric, isovelocity, dynamic). Time is expressed RELATIVE to
# each stimulation event's own onset (t_rel = 0 at stim on) so events from
# different steps/cycles/trials -- which occur at different absolute times
# and have different durations -- align for pooling.
#
# The plotted x-range is deliberately TRUNCATED, not the full trial: it runs
# from a short pre-stim baseline pad through the LONGEST stimulation duration
# observed in that category plus a fixed post-stim relaxation window (shorter
# stim events simply run out of data before that point and drop out of the
# pooled mean, which is expected and correctly shows fewer contributors at
# late t_rel).

library(dplyr)
library(ggplot2)

# PI-selected cutoff (2026-07-15, from diag_torque_smoothing_timedomain.png's
# green trace, the 8 Hz low-pass option) -- applied pipeline-wide to BOTH
# the display smoother and the pre-subtraction noise filter, per PI request
# to lock in one signal-processing choice across all trials.
DISPLAY_SMOOTH_HZ <- 8  # display-only low-pass cutoff for individual trace lines
# Isovelocity's pointwise (active-minus-matched-no-stim) subtraction
# differences two INDEPENDENTLY noisy raw torque traces, roughly doubling
# noise variance relative to isometric's scalar-minus-constant subtraction.
# Unlike DISPLAY_SMOOTH_HZ (cosmetic, applied only to the plotted line),
# NOISE_FILTER_HZ is applied to the RAW torque traces themselves BEFORE
# subtraction -- it changes the computed muscle_force_Nm values (and
# therefore the mean line too), not just how they're drawn.
NOISE_FILTER_HZ   <- 8

BASELINE_PAD_S    <- 0.2  # shown before stim onset, to see force rise from rest
RELAXATION_WINDOW_S <- 1.0  # shown after stim offset, to see the decay tail
# (RELAXATION_WINDOW_S is a visualization choice for showing decay dynamics;
# it is intentionally longer than the 0.5 s DEACTIVATION_WINDOW_S used
# elsewhere in the pipeline to compute the scalar "active force" mean.)

#' Stimulation-burst detector for dynamic trials, from the per-sample `stim`
#' column ("0"/"L"/"R"). Dynamic trials have no per-step stim_t0_s/
#' stim_t1_s metadata table (unlike isometric/isovelocity's index_step_*
#' attrs) -- BUT the sample-level `stim` signal here is NOT a sustained
#' on/off window like the segmented protocols; it is a sparse pulse TRAIN
#' (individual ~1-2-sample trigger pulses roughly 13 ms apart within one
#' stimulation burst, confirmed by inspection of a bass16 dynamic file).
#' A plain run-length on/off detector would treat every single pulse as its
#' own near-zero-duration "event". Consecutive pulses are merged into one
#' burst whenever the gap to the previous pulse is <= merge_gap_s (default
#' 0.1 s, comfortably above the ~13 ms intra-burst pulse spacing and well
#' below the >0.5 s gaps observed between bursts).
#'
#' Bursts are detected SEPARATELY PER SIDE (L, R): within a single work-loop
#' cycle the LEFT and RIGHT muscles are stimulated in antiphase alternation
#' (confirmed by inspection -- the same nominal `cycle` value contains both
#' L and R pulses), so a side-agnostic detector would merge two DIFFERENT
#' muscles' pulses into one "event" and mix their (oppositely-signed) force
#' contributions in a single trace.
#' @return tibble(event_id, muscle_side, stim_t0_s, stim_t1_s) -- one row per
#'   same-side burst.
.detect_stim_events <- function(td, merge_gap_s = 0.1) {
  purrr::map_dfr(c("L", "R"), function(side) {
    mask <- as.character(td$stim) == side
    mask[is.na(mask)] <- FALSE
    t_on <- td$t.s[mask]
    if (length(t_on) == 0L) return(NULL)
    grp <- cumsum(c(TRUE, diff(t_on) > merge_gap_s))
    tibble::tibble(t = t_on, grp = grp) |>
      dplyr::group_by(.data$grp) |>
      dplyr::summarise(stim_t0_s = min(.data$t), stim_t1_s = max(.data$t), .groups = "drop") |>
      dplyr::transmute(muscle_side = side, stim_t0_s = .data$stim_t0_s, stim_t1_s = .data$stim_t1_s)
  }) |>
    dplyr::mutate(event_id = dplyr::row_number())
}

#' Continuous per-sample muscle-force time series for isometric/isovelocity,
#' windowed to [stim_t0 - baseline_pad, stim_t1 + relaxation] per step, one
#' unit per (trial, step, side). Unlike .attach_segmented_muscle_force()
#' (which broadcasts a single step-level SCALAR across the active window for
#' the per-trial compound plot), this recomputes force CONTINUOUSLY per
#' sample -- force_sign * (torque_inertia_corrected_Nm - passive_force_Nm) --
#' so the rise/plateau/decay time-course is visible, not a flat line.
#' @param td The analyze_isometric()/analyze_isovelocity() `td` output
#'   (already has torque_inertia_corrected_Nm, step_number, t.s).
#' @param step_summary The matching step_summary (has muscle_side, force_sign,
#'   passive_force_Nm, stim_t0_s, stim_t1_s per step).
build_segmented_force_timeseries <- function(td, step_summary, trial_id,
                                              baseline_pad_s = BASELINE_PAD_S,
                                              relaxation_s = RELAXATION_WINDOW_S) {
  ss <- dplyr::filter(step_summary, .data$muscle_side %in% c("left", "right"),
                      is.finite(.data$stim_t0_s), is.finite(.data$stim_t1_s),
                      is.finite(.data$force_sign), is.finite(.data$passive_force_Nm))
  if (nrow(ss) == 0L) return(tibble::tibble())

  purrr::map_dfr(seq_len(nrow(ss)), function(i) {
    s <- ss[i, ]
    win <- td$step_number == s$step_number &
      td$t.s >= (s$stim_t0_s - baseline_pad_s) & td$t.s <= (s$stim_t1_s + relaxation_s)
    if (!any(win, na.rm = TRUE)) return(NULL)
    tibble::tibble(
      trial_id         = trial_id,
      unit_id          = paste0(trial_id, "_step", s$step_number),
      muscle_side      = s$muscle_side,
      t_rel            = td$t.s[win] - s$stim_t0_s,
      muscle_force_Nm  = s$force_sign * (td$torque_inertia_corrected_Nm[win] - s$passive_force_Nm),
      stim_duration_s  = s$stim_t1_s - s$stim_t0_s
    )
  })
}

#' Zero-phase low-pass filter, applied to a full contiguous trace (not just a
#' truncated window, to avoid edge distortion at the window boundary).
.lowpass_filtfilt <- function(y, cutoff_hz, fs_hz = 1000, order = 4) {
  ok <- is.finite(y)
  if (sum(ok) < 3L * order || cutoff_hz >= fs_hz / 2) return(y)
  filt <- signal::butter(order, cutoff_hz / (fs_hz / 2), type = "low")
  y[ok] <- tryCatch(signal::filtfilt(filt, y[ok]), error = function(e) y[ok])
  y
}

#' Continuous per-sample muscle-force time series for isovelocity steps.
#'
#' A single SCALAR passive baseline (as used in build_segmented_force_timeseries())
#' only cancels the specimen's own position/velocity-dependent passive bending
#' torque IN THE MEAN over a matched window -- it does NOT cancel pointwise,
#' because that passive torque itself varies continuously through the ramp.
#' Subtracting a scalar from a continuously-varying raw signal just reproduces
#' the raw ramp's shape (large swings unrelated to stimulation). Instead, this
#' interpolates the MATCHED no-stim step's own continuous torque trace onto
#' the active step's t_rel grid (both steps run the same commanded
#' operating_point, so they trace comparable motion profiles) and subtracts
#' pointwise -- analogous to how calc_muscle_torque() phase-bins dynamic
#' cycles rather than using one scalar.
build_isovelocity_force_timeseries <- function(td, step_summary, trial_id,
                                                baseline_pad_s = BASELINE_PAD_S,
                                                relaxation_s = RELAXATION_WINDOW_S) {
  ss <- dplyr::filter(step_summary, .data$muscle_side %in% c("left", "right"),
                      is.finite(.data$stim_t0_s), is.finite(.data$stim_t1_s), is.finite(.data$force_sign))
  if (nrow(ss) == 0L) return(tibble::tibble())

  is_no_stim_step <- vapply(step_summary$step_number, function(sn) {
    rows <- td$step_number == sn
    !any(as.character(td$stim[rows]) != "0", na.rm = TRUE)
  }, logical(1L))
  no_stim_steps <- step_summary$step_number[is_no_stim_step]
  no_stim_op    <- step_summary$operating_point[is_no_stim_step]

  purrr::map_dfr(seq_len(nrow(ss)), function(i) {
    s <- ss[i, ]
    match_idx <- which(abs(no_stim_op - s$operating_point) < 1e-3)
    if (length(match_idx) == 0L) return(NULL)  # no velocity-matched baseline available -- skip rather than fall back to a misleading scalar-subtracted trace
    match_step <- no_stim_steps[match_idx[1L]]
    match_row  <- step_summary[step_summary$step_number == match_step, ][1L, ]
    if (!is.finite(match_row$stim_t0_s)) return(NULL)

    # Filter each step's FULL torque trace (not just the truncated window)
    # before windowing/subtracting -- avoids edge distortion at the window
    # boundary and denoises the pointwise active-minus-matched subtraction
    # at the source (see NOISE_FILTER_HZ docstring above).
    act_step_rows <- td$step_number == s$step_number
    act_step_t    <- td$t.s[act_step_rows]
    act_step_filt <- .lowpass_filtfilt(td$torque_inertia_corrected_Nm[act_step_rows], NOISE_FILTER_HZ)
    act_win <- act_step_t >= (s$stim_t0_s - baseline_pad_s) & act_step_t <= (s$stim_t1_s + relaxation_s)
    if (!any(act_win, na.rm = TRUE)) return(NULL)
    t_rel_act   <- act_step_t[act_win] - s$stim_t0_s
    act_torque  <- act_step_filt[act_win]

    pass_win    <- td$step_number == match_step
    pass_filt   <- .lowpass_filtfilt(td$torque_inertia_corrected_Nm[pass_win], NOISE_FILTER_HZ)
    t_rel_pass  <- td$t.s[pass_win] - match_row$stim_t0_s
    pass_torque <- pass_filt
    if (sum(is.finite(t_rel_pass) & is.finite(pass_torque)) < 2L) return(NULL)

    pass_interp <- tryCatch(
      stats::approx(t_rel_pass, pass_torque, xout = t_rel_act, rule = 2)$y,
      error = function(e) rep(NA_real_, length(t_rel_act))
    )

    tibble::tibble(
      trial_id         = trial_id,
      unit_id          = paste0(trial_id, "_step", s$step_number),
      muscle_side      = s$muscle_side,
      contraction_mode = s$contraction_mode,
      t_rel            = t_rel_act,
      muscle_force_Nm  = s$force_sign * (act_torque - pass_interp),
      stim_duration_s  = s$stim_t1_s - s$stim_t0_s
    )
  })
}

#' Continuous per-sample muscle-force time series for dynamic trials.
#' td$muscle_force_Nm is ALREADY continuous+side-corrected here
#' (calc_muscle_torque() with include_all_active_samples = TRUE, force_sign
#' applied per-sample from that sample's own `stim` value -- see
#' .attach_dynamic_muscle_force() in the pipeline).
#'
#' Each event's window is additionally masked to `stim %in% c("0", side)` --
#' since L/R bursts interleave within the same cycle (see .detect_stim_events()
#' docstring), a naive time-window slice would also pick up the OTHER side's
#' pulses, which carry a real but DIFFERENT muscle's force_sign baked in and
#' must not be plotted as if they belonged to this event/side.
#' First finite value of `col` among the burst's own pulse samples (stim ==
#' side, between stim_t0_s and stim_t1_s) -- the commanded freq/amp/duty/
#' phase are constant for the duration of one burst, so any matching sample
#' gives the right value; NA if the column is absent or entirely NA there.
.dynamic_event_param <- function(td, col, t0, t1, side) {
  if (!col %in% names(td)) return(NA_real_)
  rows <- td$t.s >= t0 & td$t.s <= t1 & as.character(td$stim) == side
  v <- td[[col]][rows]
  v <- v[is.finite(v)]
  if (length(v) == 0L) NA_real_ else v[1L]
}

#' Continuous per-sample muscle-force time series for dynamic trials.
#' td$muscle_force_Nm is ALREADY continuous+side-corrected here
#' (calc_muscle_torque() with include_all_active_samples = TRUE, force_sign
#' applied per-sample from that sample's own `stim` value -- see
#' .attach_dynamic_muscle_force() in the pipeline).
#'
#' Each event's window is additionally masked to `stim %in% c("0", side)` --
#' since L/R bursts interleave within the same cycle (see .detect_stim_events()
#' docstring), a naive time-window slice would also pick up the OTHER side's
#' pulses, which carry a real but DIFFERENT muscle's force_sign baked in and
#' must not be plotted as if they belonged to this event/side.
#'
#' Also attaches the commanded freq_hz/amp_deg/duty/phase for this burst (PI
#' request, 2026-07-15) so dynamic trials -- which mix multiple commanded
#' conditions -- can be broken out into separate plots per condition instead
#' of one pooled plot mixing e.g. different frequencies together.
build_dynamic_force_timeseries <- function(td, trial_id,
                                            baseline_pad_s = BASELINE_PAD_S,
                                            relaxation_s = RELAXATION_WINDOW_S) {
  events <- .detect_stim_events(td)
  if (nrow(events) == 0L || !"muscle_force_Nm" %in% names(td)) return(tibble::tibble())

  purrr::map_dfr(seq_len(nrow(events)), function(i) {
    e <- events[i, ]
    win <- td$t.s >= (e$stim_t0_s - baseline_pad_s) & td$t.s <= (e$stim_t1_s + relaxation_s) &
      as.character(td$stim) %in% c("0", e$muscle_side)
    if (!any(win, na.rm = TRUE) || !any(is.finite(td$muscle_force_Nm[win]))) return(NULL)
    # e$muscle_side is "L"/"R" (matches the raw `stim` code, needed above for
    # the window mask); the OUTPUT muscle_side is spelled out "left"/"right"
    # to match the convention used everywhere else (FL/FV summary plots,
    # SIDE_COLORS) -- otherwise scale_color_manual(values = SIDE_COLORS)
    # silently fails to match ("No shared levels") and falls back to a
    # default palette, making side-adjustment impossible to verify visually.
    side_full <- c(L = "left", R = "right")[[e$muscle_side]]
    tibble::tibble(
      trial_id         = trial_id,
      unit_id          = paste0(trial_id, "_ev", e$event_id, "_", e$muscle_side),
      muscle_side      = side_full,
      freq_hz          = .dynamic_event_param(td, "freq.Hz", e$stim_t0_s, e$stim_t1_s, e$muscle_side),
      amp_deg          = .dynamic_event_param(td, "amp.deg", e$stim_t0_s, e$stim_t1_s, e$muscle_side),
      duty             = .dynamic_event_param(td, "duty",    e$stim_t0_s, e$stim_t1_s, e$muscle_side),
      phase            = .dynamic_event_param(td, "phase",   e$stim_t0_s, e$stim_t1_s, e$muscle_side),
      t_rel            = td$t.s[win] - e$stim_t0_s,
      muscle_force_Nm  = td$muscle_force_Nm[win],
      stim_duration_s  = e$stim_t1_s - e$stim_t0_s
    )
  })
}

#' Zero-phase low-pass smoother for ONE unit's trace, DISPLAY-ONLY (never
#' applied to the values feeding step_summary/FL/FV fits). Per-sample
#' inertial-corrected torque is inherently noisy; a scalar-minus-constant
#' subtraction (isometric) doesn't amplify it, but the isovelocity
#' pointwise matched-trace subtraction differences two independently noisy
#' signals, roughly doubling the noise VARIANCE -- this cutoff exists purely
#' so the individual raw-data lines are legible; the bold mean line is
#' always computed from the raw (unsmoothed) pooled data, not this.
.smooth_trace_display_only <- function(y, cutoff_hz = DISPLAY_SMOOTH_HZ, fs_hz = 1000, order = 4) {
  .lowpass_filtfilt(y, cutoff_hz = cutoff_hz, fs_hz = fs_hz, order = order)
}

SIDE_COLORS <- c(left = "#1d4ed8", right = "#b91c1c")

#' Pooled force-vs-time plot for one protocol category: each stim
#' event/step drawn as a thin, semi-transparent line (grouped by unit_id,
#' DISPLAY-smoothed -- see .smooth_trace_display_only()), overlaid with a
#' bold binned mean line (computed from raw, unsmoothed data) across all units.
#' @param df Combined output of build_segmented_force_timeseries()/
#'   build_isovelocity_force_timeseries()/build_dynamic_force_timeseries()
#'   across all trials in the category.
#' @param facet_var Optional column name to facet_wrap() by (e.g.
#'   "contraction_mode" for isovelocity's concentric/eccentric split, or
#'   "freq_hz"/"amp_deg"/"duty"/"phase" for dynamic's commanded-condition
#'   breakdowns) -- the mean +/- SE line is recomputed WITHIN each facet.
#' @param color_var Column to color individual lines by; "muscle_side" uses
#'   the same left=blue/right=red convention as the FL/FV summary plots, so
#'   side-adjustment can be checked visually (left and right should track
#'   each other within a facet, not mirror).
build_force_vs_time_plot <- function(df, title, facet_var = NULL, color_var = "trial_id",
                                      mean_bin_width_s = 0.02) {
  df <- dplyr::filter(df, is.finite(.data$t_rel), is.finite(.data$muscle_force_Nm))
  if (!is.null(facet_var)) df <- dplyr::filter(df, !is.na(.data[[facet_var]]))
  dur_range <- range(df$stim_duration_s, na.rm = TRUE)

  mean_grp <- c(if (!is.null(facet_var)) facet_var, "t_bin")
  mean_line <- df |>
    dplyr::mutate(t_bin = round(.data$t_rel / mean_bin_width_s) * mean_bin_width_s) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(mean_grp))) |>
    dplyr::summarise(
      y_mean = mean(.data$muscle_force_Nm, na.rm = TRUE),
      y_se   = stats::sd(.data$muscle_force_Nm, na.rm = TRUE) / sqrt(dplyr::n()),
      n      = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(y_se = dplyr::if_else(is.finite(.data$y_se), .data$y_se, 0))

  df_display <- df |>
    dplyr::arrange(.data$unit_id, .data$t_rel) |>
    dplyr::group_by(.data$unit_id) |>
    dplyr::mutate(muscle_force_display_Nm = .smooth_trace_display_only(.data$muscle_force_Nm)) |>
    dplyr::ungroup()

  # Stimulation-duration bar along the x-axis: darker band [0, min duration]
  # spans where EVERY pooled unit is still being stimulated; lighter band
  # [min, max duration] covers the extra range where only the longer-duration
  # units are still on (durations vary per unit -- see subtitle for n/range).
  p <- ggplot(df_display, aes(x = .data$t_rel, y = .data$muscle_force_display_Nm)) +
    annotate("rect", xmin = 0, xmax = dur_range[2L], ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.07) +
    annotate("rect", xmin = 0, xmax = dur_range[1L], ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.10) +
    geom_line(aes(group = .data$unit_id, color = .data[[color_var]]), alpha = 0.35, linewidth = 0.4) +
    geom_ribbon(data = mean_line, aes(x = .data$t_bin, y = .data$y_mean,
                                      ymin = .data$y_mean - .data$y_se, ymax = .data$y_mean + .data$y_se),
                inherit.aes = FALSE, alpha = 0.25, fill = "black") +
    geom_line(data = mean_line, aes(x = .data$t_bin, y = .data$y_mean),
              inherit.aes = FALSE, color = "black", linewidth = 1.1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    labs(title = title,
         subtitle = sprintf(
           "n = %d stim event(s)/step(s) across %d trial(s); stim duration range %.2f-%.2f s (orange band = stim ON, darker where ALL units still on); x-axis: %.1fs pre-stim baseline through longest stim + %.1fs relaxation; individual lines %d Hz low-pass (display only); bold = binned mean +/- SE per %s(from raw data), dashed vline = stim onset",
           dplyr::n_distinct(df$unit_id), dplyr::n_distinct(df$trial_id),
           dur_range[1L], dur_range[2L], BASELINE_PAD_S, RELAXATION_WINDOW_S, DISPLAY_SMOOTH_HZ,
           if (!is.null(facet_var)) paste0(facet_var, " ") else ""),
         x = "Time relative to stimulus onset (s)", y = "Muscle force (N*m)", color = color_var) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")

  if (identical(color_var, "muscle_side")) p <- p + scale_color_manual(values = SIDE_COLORS, name = "Muscle side")
  if (!is.null(facet_var)) p <- p + facet_wrap(stats::as.formula(paste0("~", facet_var)), scales = "free_y", labeller = label_both)
  p
}
