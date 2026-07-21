# extract_dynamic_l0_bookends.R
# Detects the FOUR static L0 (commanded 0 deg) isometric stim bursts that
# bracket every "dynamic" (continuous sinusoidal length-oscillation) trial:
# left+right BEFORE the cycling starts, left+right AFTER it ends.
#
# WHY (PI-directed, 2026-07-21): plot_fatigue_timeline.R's L0-vs-elapsed-time
# fatigue plot was fed only by isometric-trial data, but dynamic trial files
# ALSO contain genuine stimulated L0 contractions -- they are just not
# indexed anywhere in index_step_* metadata (dynamic/single_finite trials
# have exactly one index_step_number, dim=1). Recovering them requires
# detecting the bursts directly from the raw per-sample stim/angle columns.
#
# DETECTION: reuses .detect_stim_events() (plot_force_vs_time.R) rather than
# writing a second gap-tolerant burst merger -- that function already solves
# exactly the "stim pulses at 75 Hz, ~13 ms apart within one burst" merging
# problem for THIS SAME kind of dynamic-file per-sample stim signal (it is
# what the main dynamic-trial force path uses to find its cycling bursts),
# so a naive rle()-on-raw-stim reimplementation here would just risk a
# second, subtly different merge rule for the same underlying signal.
#
# What .detect_stim_events() cannot do on its own is tell a real ISOMETRIC
# bookend burst apart from one of the many CYCLING bursts during active
# oscillation -- both are "a burst of stim pulses on one side". The
# discriminator is sd(angle_commanded_degree) over the burst window: ~0 for
# a bookend (motor stationary at 0 deg) vs. 0.4-2+ deg for a cycling burst
# (motor actively oscillating) -- verified directly on
# 2026-07-14_bass16_bender_01_dynamic.h5. angle_sd_threshold_deg sits well
# below the smallest real cycling amplitude seen there.
#
# BASELINE-OVERLAP GOTCHA: the two bursts within one bookend pair (e.g.
# PRE-left then PRE-right) are only ~1s apart
# (protocol_prepoststim_separation_second), so a full-length passive
# baseline window immediately before a burst can swallow the OTHER burst's
# own active stim. .dynlb_clean_baseline_window() starts from a short
# window and shrinks it (from the start, keeping the end pinned at the
# burst onset) until no sample in it carries ANY stim -- verified
# empirically per bookend, not assumed safe from the nominal 1s spacing
# alone (onset jitter is real, not guaranteed by the protocol JSON).

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(cli)
})

#' Empty 0-row result with the full output schema (so callers can
#' unconditionally rbind/iterate without a separate nrow==0 branch).
.dynlb_empty <- function() {
  tibble::tibble(
    muscle_side = character(0), stim_t0_s = numeric(0), stim_t1_s = numeric(0),
    t_pre_baseline_start_s = numeric(0), t_pre_baseline_end_s = numeric(0),
    operating_point = numeric(0), contraction_mode = character(0),
    angle_sd_deg = numeric(0)
  )
}

#' Shrink a candidate pre-burst baseline window (end pinned at the burst
#' onset `t0`, start moved later in `step_s` increments) until no sample in
#' [start, t0) carries stim on EITHER side, or the window falls below
#' `min_window_s` (in which case no clean window exists and the caller
#' should drop this bookend rather than trust a contaminated baseline).
.dynlb_clean_baseline_window <- function(td, t0, window_s, min_window_s, step_s = 0.02) {
  start <- t0 - window_s
  while ((t0 - start) >= min_window_s) {
    in_win <- td$t.s >= start & td$t.s < t0
    contaminated <- any(as.character(td$stim[in_win]) != "0", na.rm = TRUE)
    if (!contaminated) return(list(start = start, end = t0, ok = TRUE))
    start <- start + step_s
  }
  list(start = NA_real_, end = NA_real_, ok = FALSE)
}

#' Detect the L0 (commanded 0 deg) isometric bookend contractions embedded
#' in one loaded dynamic trial's `td` (from load_bender_flat() -- needs
#' `t.s`, `angle.deg`, `stim`).
#'
#' @param td Loaded dynamic-trial tibble.
#' @param merge_gap_s Passed through to .detect_stim_events() -- max gap (s)
#'   between same-side stim samples still counted as one burst.
#' @param angle_sd_threshold_deg Max sd(angle.deg) over a burst window for it
#'   to be classified isometric (bookend) rather than cyclic. 0.05 deg sits
#'   well below the smallest real cycling amplitude observed and well above
#'   floating-point/encoder noise on a truly stationary command.
#' @param baseline_window_s,baseline_min_window_s Candidate/minimum pre-burst
#'   baseline window length (s) -- see .dynlb_clean_baseline_window().
#' @return tibble(muscle_side, stim_t0_s, stim_t1_s, t_pre_baseline_start_s,
#'   t_pre_baseline_end_s, operating_point (= 0), contraction_mode
#'   (= "isometric_zero"), angle_sd_deg) -- 0 to 4 rows.
detect_dynamic_l0_bookends <- function(td,
                                       merge_gap_s = 0.1,
                                       angle_sd_threshold_deg = 0.05,
                                       baseline_window_s = 0.4,
                                       baseline_min_window_s = 0.15) {
  if (!all(c("t.s", "angle.deg", "stim") %in% names(td))) {
    cli::cli_warn("detect_dynamic_l0_bookends: td missing t.s/angle.deg/stim -- skipped")
    return(.dynlb_empty())
  }

  events <- tryCatch(.detect_stim_events(td, merge_gap_s = merge_gap_s),
                     error = function(e) tibble::tibble())
  if (nrow(events) == 0L) return(.dynlb_empty())

  events$angle_sd_deg <- vapply(seq_len(nrow(events)), function(i) {
    win <- td$t.s >= events$stim_t0_s[i] & td$t.s <= events$stim_t1_s[i]
    stats::sd(td$angle.deg[win], na.rm = TRUE)
  }, numeric(1L))

  bookends <- events |>
    dplyr::filter(is.finite(.data$angle_sd_deg), .data$angle_sd_deg < angle_sd_threshold_deg) |>
    dplyr::arrange(.data$stim_t0_s)
  if (nrow(bookends) == 0L) return(.dynlb_empty())

  rows <- lapply(seq_len(nrow(bookends)), function(i) {
    b <- bookends[i, ]
    win <- .dynlb_clean_baseline_window(td, b$stim_t0_s, baseline_window_s, baseline_min_window_s)
    if (!win$ok) {
      cli::cli_warn(
        "detect_dynamic_l0_bookends: no stim-free pre-burst baseline for {if (b$muscle_side == 'L') 'left' else 'right'} burst at t={round(b$stim_t0_s, 2)}s -- dropped"
      )
      return(NULL)
    }
    tibble::tibble(
      muscle_side            = if (b$muscle_side == "L") "left" else "right",
      stim_t0_s               = b$stim_t0_s,
      stim_t1_s               = b$stim_t1_s,
      t_pre_baseline_start_s  = win$start,
      t_pre_baseline_end_s    = win$end,
      operating_point         = 0.0,
      contraction_mode        = "isometric_zero",
      angle_sd_deg            = b$angle_sd_deg
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) return(.dynlb_empty())
  out <- dplyr::bind_rows(rows)

  if (nrow(out) > 4L) {
    cli::cli_warn(
      "detect_dynamic_l0_bookends: {nrow(out)} low-angle-sd stim bursts detected (protocol design expects <= 4 -- pre/post x left/right) -- check angle_sd_threshold_deg"
    )
  }
  out
}
