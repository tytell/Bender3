# plot_angle_sono_validation.R
# Three PI-requested validation figures, independent of (and complementary
# to) plot_strain_validation.R's commanded-vs-encoder strain checks:
#   1) angleValid.png            -- measured (encoder) vs. commanded angle
#   2) strainValidSonoEnc.png   -- measured (sonomicrometry) vs.
#      predicted (from ENCODER angle) strain
#   3) strainValidSonoCmd.png -- measured (sonomicrometry) vs.
#      predicted (from COMMANDED angle) strain
#
# All three use ALL samples (no active-stim restriction) -- these are purely
# mechanical/geometric checks (does the motor track its command? does
# curvature-derived strain track real segment-length change?), not anything
# about stimulation, so passive-only trials (frequency_sweep, passive-only
# dynamic) contribute too.
#
# Sonomicrometry scope: this rig config wires exactly ONE physical sono
# channel, ai6 = "sono_right" (confirmed via daq_ai_channel_map/
# daq_sono_channels in every bass16 file inspected) -- there is a
# calibration_sono_left_volt_millimeter_breakpoints value in /metadata, but
# it does not correspond to any real acquired channel. Figures 2/3 therefore
# validate the RIGHT muscle only, using ALL samples/trials (not just
# right-recruited steps) -- the curvature-to-strain geometric relationship at
# the right muscle's position holds regardless of which side is being
# electrically stimulated at any given moment.

library(dplyr)
library(ggplot2)

PROTOCOL_FAMILY_LEVELS <- c("dynamic", "isovelocity", "isometric", "frequency_sweep")

# Three per-sample stimulation states (rig only ever drives "0"/"L"/"R" --
# confirmed no bilateral "B" code across the bass16 corpus). "no stim" is
# distinct from -- and does not imply -- resting position: it just means
# neither electrode was firing at that instant, which still matters here
# because a mechanical/geometric relationship (angle tracking, sono-vs-
# curvature strain) could plausibly differ under active muscle tone.
STIM_STATE_LEVELS <- c("no stim", "left stim", "right stim")
STIM_STATE_COLORS <- c("no stim" = "grey80", "left stim" = "#1d4ed8", "right stim" = "#b91c1c")

#' Per-sample stim state, WINDOW-based rather than literal-instant. The raw
#' `stim` column is a sparse pulse train (individual ~1-2-sample triggers,
#' confirmed ~75 pulses/s -- matches the 75 pps Grass S88 rate, see
#' .detect_stim_events() docstring in plot_force_vs_time.R): coloring by the
#' literal per-sample code would mark >90% of an actively-stimulated window
#' as "no stim" (the gaps between pulses), making the color essentially
#' invisible. Instead, each contiguous burst (pulses merged if <=
#' merge_gap_s apart, same rule as .detect_stim_events()) holds its side's
#' label for its FULL [stim_t0_s, stim_t1_s] span -- i.e. "is this muscle
#' currently being commanded to contract", not "is a pulse literally firing
#' on this exact sample".
.stim_window_state_label <- function(td, merge_gap_s = 0.1) {
  state <- rep("no stim", nrow(td))
  events <- tryCatch(.detect_stim_events(td, merge_gap_s = merge_gap_s), error = function(e) tibble::tibble())
  if (nrow(events) > 0L) {
    for (i in seq_len(nrow(events))) {
      e <- events[i, ]
      in_win <- td$t.s >= e$stim_t0_s & td$t.s <= e$stim_t1_s
      state[in_win] <- if (e$muscle_side == "L") "left stim" else "right stim"
    }
  }
  factor(state, levels = STIM_STATE_LEVELS)
}

# =============================================================================
# Sono channel read + calibration (segmented_finite-aware)
# =============================================================================

#' 0-based -> 1-based R column index of the "sono_right" channel within the
#' 8-channel aidata array, from the daq_ai_channel_map attribute.
.sono_right_channel_col <- function(m_a) {
  chan_map_raw <- m_a("daq_ai_channel_map")
  if (is.null(chan_map_raw)) return(NA_integer_)
  chan_map <- tryCatch(jsonlite::fromJSON(as.character(chan_map_raw)), error = function(e) NULL)
  if (is.null(chan_map)) return(NA_integer_)
  chan_names <- as.character(unname(chan_map))
  idx <- which(grepl("sono_right", chan_names, fixed = TRUE))
  if (length(idx) == 0L) return(NA_integer_)
  as.integer(names(chan_map)[idx[1L]]) + 1L
}

#' Read + calibrate the sono_right channel, aligned row-for-row with
#' load_bender_flat()'s td (same sample count/order). 01_calibrate.R's
#' calibrate_bender()/.cal_sono() only reads the top-level
#' "/timeseries/aidata" path, which exists for single_finite protocols
#' (dynamic, frequency_sweep) but NOT segmented_finite ones (isometric,
#' isovelocity store aidata per step_NNN/ subgroup instead) -- this mirrors
#' 00_load_bender_flat.R's .bfl_read_segmented_finite() step-concatenation
#' order so the result lines up with td.
#' @return Numeric vector (mm), or NULL if sono is unavailable/uncalibrated
#'   for this file.
.read_sono_right_mm_aligned <- function(raw_path) {
  h5 <- rhdf5::H5Fopen(raw_path, "H5F_ACC_RDONLY")
  on.exit(try(rhdf5::H5Fclose(h5), silent = TRUE), add = TRUE)

  m_attrs <- tryCatch(rhdf5::h5readAttributes(h5, "/metadata"), error = function(e) list())
  m_a  <- function(key) m_attrs[[key]]
  m_ds <- function(key) {
    path <- paste0("/metadata/", key)
    if (rhdf5::H5Lexists(h5, path)) tryCatch(rhdf5::h5read(h5, path), error = function(e) NULL) else NULL
  }

  bp <- tryCatch(as.numeric(m_ds("calibration_sono_right_volt_millimeter_breakpoints")), error = function(e) NULL)
  if (is.null(bp) || length(bp) < 4L || any(!is.finite(bp))) return(NULL)

  col_idx <- .sono_right_channel_col(m_a)
  if (is.na(col_idx)) return(NULL)

  .read_aidata_col <- function(path) {
    if (!rhdf5::H5Lexists(h5, path)) return(NULL)
    m <- tryCatch(rhdf5::h5read(h5, path), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    if (!is.matrix(m)) m <- as.matrix(m)
    if (ncol(m) != 8L && nrow(m) == 8L) m <- t(m)
    if (ncol(m) < col_idx) return(NULL)
    m[, col_idx]
  }

  sampling_mode <- .bfl_chr(m_a("protocol_sampling_mode"), "single_finite")
  v_raw <- if (sampling_mode == "single_finite") {
    .read_aidata_col("/timeseries/aidata")
  } else {
    ts_info <- tryCatch(rhdf5::h5ls(h5, recursive = TRUE), error = function(e) NULL)
    if (is.null(ts_info)) return(NULL)
    step_names <- sort(ts_info$name[ts_info$group == "/timeseries" & grepl("^step_\\d+$", ts_info$name)])
    if (length(step_names) == 0L) return(NULL)
    do.call(c, lapply(step_names, function(sname) .read_aidata_col(paste0("/timeseries/", sname, "/aidata"))))
  }
  if (is.null(v_raw) || length(v_raw) == 0L) return(NULL)

  .cal_sono_interp(v_raw, bp)  # reuse 01_calibrate.R's private breakpoint interpolator
}

# =============================================================================
# Sono strain (right muscle) + right-folded predicted strain
# =============================================================================

SONO_ZERO_ANGLE_TOL_DEG     <- 2.0   # widened progressively if too few near-zero samples
SONO_ZERO_ANGLE_MAX_TOL_DEG <- 20.0
SONO_ZERO_ANGLE_MIN_N       <- 5L

#' Resting-length reference (mm): mean sono_right_mm over samples where the
#' COMMANDED angle is near zero (straight body, ~zero curvature -- the same
#' zero-strain reference point the curvature-geometry formula assumes).
#' Tolerance widens (up to SONO_ZERO_ANGLE_MAX_TOL_DEG) if too few samples
#' qualify; falls back to the single closest-to-zero sample if none do.
.sono_reference_length_mm <- function(angle_deg, sono_mm) {
  ok <- is.finite(angle_deg) & is.finite(sono_mm)
  if (!any(ok)) return(NA_real_)
  tol <- SONO_ZERO_ANGLE_TOL_DEG
  near0 <- ok & abs(angle_deg) <= tol
  while (sum(near0) < SONO_ZERO_ANGLE_MIN_N && tol < SONO_ZERO_ANGLE_MAX_TOL_DEG) {
    tol <- tol * 2
    near0 <- ok & abs(angle_deg) <= tol
  }
  if (sum(near0) == 0L) {
    idx <- which(ok)[which.min(abs(angle_deg[ok]))]
    return(sono_mm[idx])
  }
  mean(sono_mm[near0], na.rm = TRUE)
}

#' Attach sono-derived RIGHT-muscle strain (strain_sono_pct) and the
#' RIGHT-folded predicted strain from both commanded and encoder angle
#' (strain_pred_commanded_right_pct / strain_pred_encoder_right_pct) to a td
#' that already has strain_pct/strain_measured_pct (attach_predicted_strain()
#' + attach_measured_strain() must run first).
#'
#' Sign convention: sonomicrometry measures the RIGHT segment's physical
#' length directly, so muscle SHORTENING = length DECREASE. strain_sono_pct
#' negates the raw fractional length change so that positive = shortening,
#' matching the convention used throughout the rest of the pipeline
#' (shortening_strain_pct, muscle_force_Nm's force_sign fold). The predicted
#' columns get the ANALOGOUS fold for the right muscle specifically
#' (force_sign_right = lidx_right * lidx_pos_motor -- same algebra as
#' resolve_step_contraction()'s force_sign, but evaluated for "right"
#' unconditionally rather than for whichever side was recruited that step,
#' since the geometric relationship between commanded/encoder angle and the
#' right muscle's position is purely mechanical and holds at all times).
attach_sono_strain <- function(td, raw_path, lidx_right, lidx_pos_motor) {
  sono_mm <- tryCatch(.read_sono_right_mm_aligned(raw_path), error = function(e) NULL)
  if (is.null(sono_mm) || length(sono_mm) != nrow(td)) {
    td$sono_right_mm    <- NA_real_
    td$strain_sono_pct  <- NA_real_
  } else {
    td$sono_right_mm <- sono_mm
    L0 <- .sono_reference_length_mm(td$angle.deg, sono_mm)
    td$strain_sono_pct <- if (is.finite(L0) && L0 > 0) -(sono_mm - L0) / L0 * 100.0 else NA_real_
  }

  force_sign_right <- if (is.finite(lidx_right) && is.finite(lidx_pos_motor)) lidx_right * lidx_pos_motor else NA_real_
  td$strain_pred_commanded_right_pct <- force_sign_right * td$strain_pct
  td$strain_pred_encoder_right_pct   <- if ("strain_measured_pct" %in% names(td)) force_sign_right * td$strain_measured_pct else NA_real_
  td
}

# =============================================================================
# Plot builders
# =============================================================================

#' Shared r/RMSE/n annotation, one row per protocol_family.
.validation_r2_labels <- function(df, x_col, y_col) {
  df |>
    dplyr::group_by(.data$protocol_family) |>
    dplyr::summarise(
      r    = suppressWarnings(stats::cor(.data[[x_col]], .data[[y_col]], use = "complete.obs")),
      rmse = sqrt(mean((.data[[x_col]] - .data[[y_col]])^2, na.rm = TRUE)),
      n    = dplyr::n(), .groups = "drop"
    ) |>
    dplyr::mutate(label = sprintf("r=%.3f, RMSE=%.2f, n=%s", .data$r, .data$rmse, format(.data$n, big.mark = ",")))
}

#' Figure 1: measured (E6 encoder) vs. commanded motor angle, all protocol
#' families, ALL samples (active + passive). Points are color-coded by
#' which muscle (if any) was being stimulated at that instant (PI request,
#' 2026-07-16) -- a purely mechanical/geometric relationship should not
#' depend on stim state, so this is a visual check as much as a color key.
build_angle_validation_plot <- function(df) {
  df <- dplyr::filter(df, is.finite(.data$angle_commanded), is.finite(.data$angle_measured))
  df$protocol_family <- factor(df$protocol_family, levels = intersect(PROTOCOL_FAMILY_LEVELS, unique(df$protocol_family)))

  lims <- range(c(df$angle_commanded, df$angle_measured), na.rm = TRUE)
  ref_df <- tibble::tibble(x = lims, y = lims)
  r2_labels <- .validation_r2_labels(df, "angle_commanded", "angle_measured")

  # "no stim" points are drawn in their OWN low-alpha layer first (background),
  # stim points drawn on top afterward at higher alpha/stroke -- stim samples
  # are a small minority everywhere (raw stim is a sparse ~75 pps pulse train
  # even across a whole burst window, see .stim_window_state_label()
  # docstring) and would otherwise be invisible underneath the much denser
  # "no stim" point cloud.
  df_no  <- dplyr::filter(df, .data$stim_state == "no stim")
  df_stim <- dplyr::filter(df, .data$stim_state != "no stim")

  ggplot(df, aes(x = .data$angle_commanded, y = .data$angle_measured, color = .data$stim_state)) +
    geom_point(data = df_no, shape = 1, size = 1.0, alpha = 0.25, stroke = 0.3) +
    geom_point(data = df_stim, shape = 1, size = 1.2, alpha = 0.6, stroke = 0.4) +
    geom_line(data = ref_df, aes(x = .data$x, y = .data$y), inherit.aes = FALSE,
              linetype = "dashed", color = "black", linewidth = 0.6) +
    geom_text(data = r2_labels, aes(x = -Inf, y = Inf, label = .data$label),
              inherit.aes = FALSE, hjust = -0.05, vjust = 1.3, size = 3.2, color = "black") +
    facet_wrap(~protocol_family, scales = "free") +
    scale_color_manual(values = STIM_STATE_COLORS, name = "Stim state", drop = FALSE) +
    labs(title = "Measured (E6 encoder) vs. commanded motor angle",
         subtitle = "ALL samples (active + passive), pooled per protocol category; dashed = 1:1 reference; color = which muscle (if any) was being commanded to contract over that stim burst's full window",
         x = "Commanded angle (deg)", y = "Measured angle (deg, E6 encoder)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
}

#' Figures 2/3: measured (sonomicrometry, RIGHT muscle) vs. predicted
#' (curvature geometry, right-folded) strain, from either the encoder or
#' commanded angle. Points are color-coded by which muscle (if any) was
#' being stimulated at that instant (PI request, 2026-07-16) -- lets the PI
#' see whether scatter/bias differs when the (sono-instrumented) right
#' muscle itself is active vs. the left muscle vs. neither.
#' @param predicted_col "strain_pred_encoder_right_pct" or
#'   "strain_pred_commanded_right_pct".
#' @param predicted_label Short label for titles/axis text, e.g. "encoder" or
#'   "commanded".
build_sono_strain_validation_plot <- function(df, predicted_col, predicted_label) {
  df <- dplyr::filter(df, is.finite(.data[[predicted_col]]), is.finite(.data$strain_sono_pct))
  df$protocol_family <- factor(df$protocol_family, levels = intersect(PROTOCOL_FAMILY_LEVELS, unique(df$protocol_family)))

  lims <- range(c(df[[predicted_col]], df$strain_sono_pct), na.rm = TRUE)
  ref_df <- tibble::tibble(x = lims, y = lims)
  r2_labels <- .validation_r2_labels(df, predicted_col, "strain_sono_pct")

  # See build_angle_validation_plot() -- "no stim" drawn first/low-alpha,
  # stim points layered on top/higher-alpha so the minority stim samples
  # aren't hidden under the denser no-stim cloud.
  df_no   <- dplyr::filter(df, .data$stim_state == "no stim")
  df_stim <- dplyr::filter(df, .data$stim_state != "no stim")

  ggplot(df, aes(x = .data[[predicted_col]], y = .data$strain_sono_pct, color = .data$stim_state)) +
    geom_point(data = df_no, shape = 1, size = 1.0, alpha = 0.25, stroke = 0.3) +
    geom_point(data = df_stim, shape = 1, size = 1.2, alpha = 0.6, stroke = 0.4) +
    geom_line(data = ref_df, aes(x = .data$x, y = .data$y), inherit.aes = FALSE,
              linetype = "dashed", color = "black", linewidth = 0.6) +
    geom_text(data = r2_labels, aes(x = -Inf, y = Inf, label = .data$label),
              inherit.aes = FALSE, hjust = -0.05, vjust = 1.3, size = 3.2, color = "black") +
    facet_wrap(~protocol_family, scales = "free") +
    scale_color_manual(values = STIM_STATE_COLORS, name = "Stim state", drop = FALSE) +
    labs(title = sprintf("Measured (sonomicrometry, right muscle) vs. predicted (%s angle) strain", predicted_label),
         subtitle = "RIGHT muscle only (only wired sono channel); ALL samples, right-folded so positive = shortening; dashed = 1:1 reference; color = which muscle (if any) was being commanded to contract over that stim burst's full window",
         x = sprintf("Predicted strain (%%, from %s angle, right-folded)", predicted_label),
         y = "Measured strain (%, sonomicrometry, right muscle)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
}
