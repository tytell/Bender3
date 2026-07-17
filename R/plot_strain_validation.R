# plot_strain_validation.R
# Measured (E6 encoder) vs. predicted (curvature-geometry) strain, restricted
# to samples/steps where the muscle is ACTUALLY being stimulated -- validates
# the assumption (used throughout the FL/FV x-axis) that commanded motor
# position tracks the specimen's real motion. Discrepancy sources: motor
# tracking error, gearbox backlash, specimen slip in the clamp.
#
# "Predicted" strain = curvature from the COMMANDED angle (angle.deg) via
# compute_predicted_strain() (muscle_geometry.R) -- already attached to td as
# strain_pct by attach_predicted_strain().
# "Measured" strain = the SAME formula run on the E6 ENCODER angle (enc.deg,
# angle_measured_degree) instead of the commanded angle -- computed here.

library(dplyr)
library(ggplot2)

#' Attach strain_measured_pct to a td that already has strain_pct/
#' muscle_strain_r_m from attach_predicted_strain() (muscle_geometry.R).
#' Reuses the SAME r_m (radial distance) so "measured" and "predicted" only
#' differ in which angle signal (encoder vs. commanded) drove the curvature.
attach_measured_strain <- function(td) {
  req <- c("enc.deg", "dclamp.m", "muscle_strain_r_m")
  missing <- setdiff(req, names(td))
  if (length(missing) > 0L) {
    cli::cli_abort("attach_measured_strain: missing column(s) {paste(missing, collapse=', ')} -- run attach_predicted_strain() first")
  }
  curve_measured_invm <- td$enc.deg * pi / 180.0 / td$dclamp.m
  td$strain_measured_pct <- curve_measured_invm * td$muscle_strain_r_m * 100.0
  td
}

#' 3-panel (dynamic / isovelocity / isometric, pooled) measured-vs-predicted
#' strain scatter, restricted to actively-stimulated samples/steps only.
#' frequency_sweep is excluded -- it is passive-only, so "where the muscle is
#' being stimulated" is an empty set for that category.
#' @param df Combined tibble (trial_id, protocol_family, strain_pct,
#'   strain_measured_pct), pre-filtered by the caller to active-stim rows only.
build_measured_vs_predicted_strain_plot <- function(df) {
  df <- dplyr::filter(df, is.finite(.data$strain_pct), is.finite(.data$strain_measured_pct))
  df$protocol_family <- factor(df$protocol_family, levels = c("dynamic", "isovelocity", "isometric"))

  lims <- range(c(df$strain_pct, df$strain_measured_pct), na.rm = TRUE)
  ref_df <- tibble::tibble(x = lims, y = lims)

  r2_labels <- df |>
    dplyr::group_by(.data$protocol_family) |>
    dplyr::summarise(
      r = suppressWarnings(stats::cor(.data$strain_pct, .data$strain_measured_pct, use = "complete.obs")),
      rmse = sqrt(mean((.data$strain_pct - .data$strain_measured_pct)^2, na.rm = TRUE)),
      n = dplyr::n(), .groups = "drop"
    ) |>
    dplyr::mutate(label = sprintf("r=%.3f, RMSE=%.2f%%, n=%s", .data$r, .data$rmse, format(.data$n, big.mark = ",")))

  # No trial_id legend -- with a dozen-plus pooled trials the per-trial color
  # key added noise without conveying anything useful (PI feedback, 2026-07-15);
  # points are plotted in one fixed color and rely on alpha for density.
  ggplot(df, aes(x = .data$strain_pct, y = .data$strain_measured_pct)) +
    geom_point(shape = 1, size = 1.4, alpha = 0.35, stroke = 0.4, color = "#1d4ed8") +
    geom_line(data = ref_df, aes(x = .data$x, y = .data$y), inherit.aes = FALSE,
              linetype = "dashed", color = "black", linewidth = 0.6) +
    geom_text(data = r2_labels, aes(x = -Inf, y = Inf, label = .data$label),
              inherit.aes = FALSE, hjust = -0.05, vjust = 1.3, size = 3.2) +
    facet_wrap(~protocol_family, scales = "free") +
    labs(title = "Measured (E6 encoder) vs. predicted (commanded-angle) muscle strain",
         subtitle = "CONTINUOUS per-sample strain, actively-stimulated samples only, pooled per protocol category; dashed = 1:1 reference",
         x = "Predicted strain (%, from commanded angle)", y = "Measured strain (%, from E6 encoder)") +
    theme_bw(base_size = 12)
}

# =============================================================================
# Step-level (FL/FV x-axis) measured-vs-predicted strain
# =============================================================================
# The plot above validates the CONTINUOUS per-sample strain_pct (used in the
# per-trial compound-plot strain panel). It does NOT validate
# shortening_strain_pct -- the STEP-LEVEL, SIGN-FOLDED strain that is the
# actual x-axis of the FL/FV summary curves (R/03_analyze.R,
# build_segmented_step_summary()) -- because shortening_strain_pct comes from
# the step's NOMINAL COMMANDED TARGET (index_step_operating_point), broadcast
# as one constant per step, not from any continuous signal. This section adds
# the analogous MEASURED (encoder) counterpart at that same step level, so the
# FL/FV curve's own x-axis assumption can be checked directly.

#' Per-step MEASURED-strain analogue of shortening_strain_pct.
#' Averages strain_measured_pct (continuous, encoder-derived; see
#' attach_measured_strain() above) over each step's active window
#' (stim_t0_s to stim_t1_s + deactivation_window_s, matching the window used
#' for active_force_Nm in build_segmented_step_summary()), then applies the
#' IDENTICAL sign fold used to build shortening_strain_pct from operating_point
#' (concentric -> abs, eccentric -> -abs, isometric_zero -> 0) so the two
#' columns are directly comparable on the same signed shortening axis.
#' @param td Per-trial td (must have step_number, t.s, enc.deg, dclamp.m,
#'   muscle_strain_r_m -- i.e. already run through attach_predicted_strain()).
#' @param steps step_summary tibble (must have step_number, contraction_mode,
#'   stim_t0_s, stim_t1_s, shortening_strain_pct, muscle_side).
#' @return steps with one new column, strain_measured_step_pct.
attach_step_measured_strain <- function(td, steps, deactivation_window_s = 0.5) {
  td <- attach_measured_strain(td)
  steps$strain_measured_step_pct <- vapply(seq_len(nrow(steps)), function(i) {
    s <- steps[i, ]
    if (is.na(s$contraction_mode) || !is.finite(s$stim_t0_s) || !is.finite(s$stim_t1_s)) return(NA_real_)
    win <- td$step_number == s$step_number &
      td$t.s >= s$stim_t0_s & td$t.s <= (s$stim_t1_s + deactivation_window_s)
    raw <- mean(td$strain_measured_pct[win], na.rm = TRUE)
    if (!is.finite(raw)) return(NA_real_)
    dplyr::case_when(
      s$contraction_mode == "concentric"    ~ abs(raw),
      s$contraction_mode == "eccentric"     ~ -abs(raw),
      s$contraction_mode == "isometric_zero" ~ 0.0,
      .default = NA_real_
    )
  }, numeric(1L))
  steps
}

#' Step-level measured-vs-predicted strain scatter -- validates
#' shortening_strain_pct (the FL curve x-axis) against the actual encoder
#' motion achieved during each step.
#'
#' ISOMETRIC ONLY: for isometric steps, operating_point is a commanded ANGLE
#' (deg), so shortening_strain_pct is a genuine position-based strain (%),
#' directly comparable to strain_measured_step_pct (also %, from the encoder
#' angle). For ISOVELOCITY, operating_point is a commanded VELOCITY (deg/s)
#' -- shortening_strain_pct there is actually a STRAIN-RATE (%/s; see its
#' axis label in build_summary_plot_isovelocity(), plot_summary_profiles.R),
#' not a position. Plotting a %/s quantity against a % positional strain on
#' the same 1:1 axis is comparing incompatible units and was found to
#' produce a nonsensical panel (values up to +-300 with no relationship to
#' the 1:1 line) -- isovelocity is deliberately EXCLUDED here, not just
#' filtered out incidentally.
#' @param df Combined step_summary rows, ISOMETRIC rows only (must have
#'   protocol_family, shortening_strain_pct, strain_measured_step_pct,
#'   muscle_side).
build_step_strain_validation_plot <- function(df) {
  df <- dplyr::filter(df, .data$protocol_family == "isometric", .data$muscle_side %in% c("left", "right"),
                      is.finite(.data$shortening_strain_pct), is.finite(.data$strain_measured_step_pct))
  df$protocol_family <- factor(df$protocol_family, levels = "isometric")

  lims <- range(c(df$shortening_strain_pct, df$strain_measured_step_pct), na.rm = TRUE)
  ref_df <- tibble::tibble(x = lims, y = lims)

  r2_labels <- df |>
    dplyr::group_by(.data$protocol_family) |>
    dplyr::summarise(
      r = suppressWarnings(stats::cor(.data$shortening_strain_pct, .data$strain_measured_step_pct, use = "complete.obs")),
      rmse = sqrt(mean((.data$shortening_strain_pct - .data$strain_measured_step_pct)^2, na.rm = TRUE)),
      n = dplyr::n(), .groups = "drop"
    ) |>
    dplyr::mutate(label = sprintf("r=%.3f, RMSE=%.2f%%, n=%s", .data$r, .data$rmse, format(.data$n, big.mark = ",")))

  # No trial_id legend (same PI feedback as the continuous plot above);
  # colored by muscle_side instead, which IS informative here (checks
  # whether tracking fidelity differs by side) and has only 2 levels.
  ggplot(df, aes(x = .data$shortening_strain_pct, y = .data$strain_measured_step_pct)) +
    geom_point(aes(color = .data$muscle_side), shape = 1, size = 2.2, alpha = 0.7, stroke = 0.6) +
    geom_line(data = ref_df, aes(x = .data$x, y = .data$y), inherit.aes = FALSE,
              linetype = "dashed", color = "black", linewidth = 0.6) +
    geom_text(data = r2_labels, aes(x = -Inf, y = Inf, label = .data$label),
              inherit.aes = FALSE, hjust = -0.05, vjust = 1.3, size = 3.2) +
    facet_wrap(~protocol_family, scales = "free") +
    scale_color_manual(values = c(left = "#1d4ed8", right = "#b91c1c"), name = "Muscle side") +
    labs(title = "Step-level (isometric only): measured (E6 encoder) vs. predicted (commanded target) shortening strain",
         subtitle = "One point per step; this is the FL curve's OWN x-axis (shortening_strain_pct) checked against the ACTUAL encoder motion achieved during that step. Isovelocity excluded -- its shortening_strain_pct is a strain-RATE (%/s), not a position, so it isn't comparable to encoder strain on this axis. Dashed = 1:1 reference.",
         x = "Predicted shortening strain (%, from step's commanded operating_point)",
         y = "Measured shortening strain (%, from E6 encoder, same sign fold)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
}
