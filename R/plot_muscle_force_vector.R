# plot_muscle_force_vector.R
# Summary + validation plots for the 6-axis line-of-action (u_hat) muscle-force
# path (see muscle_force_vector.R). ADDITIVE: these are new "_vector"-suffixed
# figures saved alongside the existing single-axis (zTorque) FL/FV/force-vs-time
# outputs, which are left untouched.
#
# Connect-the-dots per the FL/FV default (PI direction, 2026-07-16): NO fitted
# vertex/Vmax is drawn unless explicitly requested. Reuses bin_mean_se() /
# .mean_line_by_side() (plot_summary_profiles.R) and the left=blue/right=red
# SIDE_COLORS convention.

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

MFV_SIDE_COLORS <- c(left = "#1d4ed8", right = "#b91c1c")

# MFV_CONFIDENCE_ALPHA (4-tier: confident/confidently_small/unstable_magnitude/
# unconfirmable) and mfv_confidence_tier() now live in muscle_force_vector.R
# (sourced before this file everywhere) -- REPLACES the 2-level ratio-only
# version previously defined here (PI-directed 2026-07-22, "SNR-based
# confidence gating audit"): a ratio-only "low confidence" bucket conflated a
# genuinely-small-but-real point (e.g. a weak/fatigued muscle, a brief V=0
# hold) with a point that is indistinguishable from pure noise. See
# analysis_muscle_force_vector_log.md for the full rationale/evidence.

# Method labels for the empirical-vs-reference u_hat "both ways" comparison
# (PI A5, 2026-07-18) -- ordered so empirical (primary) is the left/first panel.
MFV_METHOD_LEVELS <- c("empirical u_hat", "geometric/longitudinal u_hat")

#' Stack the two force columns (empirical + reference u_hat) into one long
#' tibble with a `method` factor, for side-by-side faceting. Carries
#' baseline_force_noise_N through (added 2026-07-22, SNR-magnitude audit) so
#' the per-point confidence tier can be computed downstream -- a ratio
#' (activation_snr) alone cannot tell "elevated noise" apart from "genuinely
#' small real force," see mfv_confidence_tier() (muscle_force_vector.R).
.mfv_long_methods <- function(step_summary_vec) {
  base_cols <- intersect(c("trial_id", "muscle_side", "contraction_mode", "shortening_strain_pct",
                           "activation_snr", "baseline_force_noise_N"),
                         names(step_summary_vec))
  emp <- step_summary_vec |>
    dplyr::transmute(dplyr::across(dplyr::all_of(base_cols)),
                     force_N = .data$muscle_force_vector_N, method = MFV_METHOD_LEVELS[1L])
  geo <- step_summary_vec |>
    dplyr::transmute(dplyr::across(dplyr::all_of(base_cols)),
                     force_N = .data$muscle_force_vector_geom_N, method = MFV_METHOD_LEVELS[2L])
  out <- dplyr::bind_rows(emp, geo)
  out$method <- factor(out$method, levels = MFV_METHOD_LEVELS)
  out
}

#' 4-tier confidence flag for a point (PI-directed, 2026-07-18, REVISED
#' 2026-07-22 to be magnitude-aware -- see mfv_confidence_tier()
#' (muscle_force_vector.R) and MFV_CONFIDENCE_ALPHA above). Thin wrapper kept
#' for this file's call sites; NA/non-finite force or noise (e.g. u_hat fell
#' back to longitudinal, no baseline_force_noise_N available) is treated as
#' "unconfirmable" by mfv_confidence_tier() itself, never dropped.
.mfv_confidence_factor <- function(force, activation_snr, baseline_force_noise_N, snr_min) {
  mfv_confidence_tier(force, activation_snr, baseline_force_noise_N, snr_min = snr_min)
}

#' Shared builder for the FL / FV vector summaries -- both are
#' force-vs-shortening(-rate), connect-the-mean per side (NO model fit), and
#' faceted by method (empirical vs geometric/longitudinal u_hat) so the two
#' force estimates sit side by side. Points are flagged with the 4-tier
#' confidence factor (REVISED 2026-07-22, see mfv_confidence_tier() /
#' MFV_CONFIDENCE_ALPHA, muscle_force_vector.R) -- a LOW-ratio flag alone,
#' as used before this revision, conflated "elevated noise" with "genuinely
#' small real force" (e.g. bass17 isometric, max SNR 2.36, might be a weak-
#' but-real signal rather than noise -- the ratio alone can't say which). All
#' 4 tiers stay visible (never a filter); only their alpha differs.
.build_vector_xy_summary <- function(step_summary_vec, x_col, x_lab, title, shape_col,
                                     snr_min = MFV_UHAT_SNR_MIN) {
  long <- .mfv_long_methods(step_summary_vec)
  pts <- dplyr::filter(long, .data$muscle_side %in% c("left", "right"),
                       is.finite(.data$force_N), is.finite(.data[[x_col]]))
  if (nrow(pts) == 0L) return(NULL)
  pts$confidence <- .mfv_confidence_factor(pts$force_N, pts$activation_snr,
                                           pts$baseline_force_noise_N, snr_min)
  tier_n <- table(pts$confidence)

  trend <- purrr::map_dfr(levels(pts$method), function(m) {
    sub <- dplyr::filter(pts, .data$method == m)
    t <- .mean_line_by_side(sub, x_col, "force_N")
    if (nrow(t) == 0L) return(NULL)
    t$method <- factor(m, levels = MFV_METHOD_LEVELS); t
  })

  ggplot(pts, aes(x = .data[[x_col]], y = .data$force_N, color = .data$muscle_side)) +
    geom_point(aes(shape = .data[[shape_col]], alpha = .data$confidence), size = 2.0) +
    { if (nrow(trend) > 0L)
        geom_ribbon(data = trend, aes(x = .data$x_mid, y = .data$y_mean,
                                      ymin = .data$y_mean - .data$y_se, ymax = .data$y_mean + .data$y_se,
                                      color = NULL, fill = .data$muscle_side, shape = NULL, alpha = NULL),
                    alpha = 0.2, inherit.aes = FALSE) } +
    { if (nrow(trend) > 0L)
        geom_line(data = trend, aes(x = .data$x_mid, y = .data$y_mean, color = .data$muscle_side,
                                    shape = NULL, alpha = NULL),
                  inherit.aes = FALSE, linewidth = 1.0) } +
    facet_wrap(~method) +
    scale_color_manual(values = MFV_SIDE_COLORS, name = "Muscle side") +
    scale_fill_manual(values = MFV_SIDE_COLORS, guide = "none") +
    scale_alpha_manual(values = MFV_CONFIDENCE_ALPHA, name = "Confidence tier (SNR x magnitude)", drop = FALSE) +
    labs(title = title,
         subtitle = sprintf("n = %d points across %d trial(s); confidence tier counts: confident=%d, confidently_small=%d, unstable_magnitude=%d, unconfirmable=%d (SNR>=%.1f AND/OR |force|>=own noise floor, see mfv_confidence_tier()); force = projected active moment / |r x u_hat|, computed BOTH ways; solid = mean per bin (NO model fit, includes ALL tiers), band = +/-SE",
                            dplyr::n_distinct(paste(pts$trial_id, pts$muscle_side, pts[[x_col]])),
                            dplyr::n_distinct(pts$trial_id),
                            tier_n[["confident"]], tier_n[["confidently_small"]],
                            tier_n[["unstable_magnitude"]], tier_n[["unconfirmable"]], snr_min),
         x = x_lab, y = "Muscle force along u_hat (N)", shape = shape_col) +
    theme_bw(base_size = 12)
}

#' Force-Length (vector), empirical vs geometric/longitudinal u_hat side-by-side.
build_summary_plot_FL_vector <- function(step_summary_vec,
                                         title = "Isometric summary: Force-Length [6-axis LOA vector, N]") {
  .build_vector_xy_summary(step_summary_vec, "shortening_strain_pct",
                           "Muscle shortening strain (%, folded by recruited side)",
                           title, shape_col = "trial_id")
}

#' Force-Velocity (vector), empirical vs geometric/longitudinal u_hat side-by-side.
build_summary_plot_FV_vector <- function(step_summary_vec,
                                         title = "Isovelocity summary: Force-Velocity [6-axis LOA vector, N]") {
  .build_vector_xy_summary(step_summary_vec, "shortening_strain_pct",
                           "Muscle shortening strain-rate (%/s, folded by recruited side; concentric > 0)",
                           title, shape_col = "contraction_mode")
}

#' Force-Velocity (vector), sono-confirmed at the L0 (resting-length) crossing.
#' RIGHT muscle only (sono is the sole wired right-side channel): for each
#' isovelocity ramp the muscle force is sampled at the instant its (40 Hz
#' low-pass) sono length crosses L0, so all FV points share a common length.
#' Connect-the-mean, no Hill fit unless explicitly requested.
#' Points are flagged with the 4-tier confidence factor the same way as
#' .build_vector_xy_summary() (REVISED 2026-07-22 -- see MFV_CONFIDENCE_ALPHA /
#' mfv_confidence_tier(), muscle_force_vector.R).
build_summary_plot_FV_L0_vector <- function(fv_l0_df,
                                            title = "Isovelocity FV [6-axis LOA vector, N] -- sono-confirmed at L0 crossing (right muscle)",
                                            snr_min = MFV_UHAT_SNR_MIN) {
  df <- dplyr::filter(fv_l0_df, is.finite(.data$force_at_L0_N), is.finite(.data$shortening_strain_pct),
                      .data$muscle_side == "right")
  if (nrow(df) == 0L) return(NULL)
  df$confidence <- .mfv_confidence_factor(df$force_at_L0_N, df$activation_snr,
                                          df$baseline_force_noise_N, snr_min)
  tier_n <- table(df$confidence)
  trend <- .mean_line_by_side(df, "shortening_strain_pct", "force_at_L0_N")

  ggplot(df, aes(x = .data$shortening_strain_pct, y = .data$force_at_L0_N)) +
    geom_point(aes(shape = .data$contraction_mode, alpha = .data$confidence),
              color = MFV_SIDE_COLORS[["right"]], size = 2.4) +
    { if (nrow(trend) > 0L)
        geom_ribbon(data = trend, aes(x = .data$x_mid, y = .data$y_mean,
                                      ymin = .data$y_mean - .data$y_se, ymax = .data$y_mean + .data$y_se),
                    inherit.aes = FALSE, alpha = 0.2, fill = MFV_SIDE_COLORS[["right"]]) } +
    { if (nrow(trend) > 0L)
        geom_line(data = trend, aes(x = .data$x_mid, y = .data$y_mean),
                  inherit.aes = FALSE, color = MFV_SIDE_COLORS[["right"]], linewidth = 1.0) } +
    scale_alpha_manual(values = MFV_CONFIDENCE_ALPHA, name = "Confidence tier (SNR x magnitude)", drop = FALSE) +
    labs(title = title,
         subtitle = sprintf("n = %d sono-confirmed L0 crossing(s) across %d trial(s); confidence tier counts: confident=%d, confidently_small=%d, unstable_magnitude=%d, unconfirmable=%d; force sampled where sono length = L0; solid = mean per bin (NO model fit, includes ALL tiers), band = +/-SE",
                            nrow(df), dplyr::n_distinct(df$trial_id),
                            tier_n[["confident"]], tier_n[["confidently_small"]],
                            tier_n[["unstable_magnitude"]], tier_n[["unconfirmable"]]),
         x = "Muscle shortening strain-rate (%/s; concentric > 0)",
         y = "Muscle force along u_hat at L0 (N)", shape = "Contraction mode") +
    theme_bw(base_size = 12)
}

#' u_hat validation: empirical (active-minus-passive force direction) vs
#' geometric (beam-tangent) line-of-action angle, per bend angle. Agreement
#' validates u_hat independently; divergence is diagnostic (PI A5).
build_uhat_comparison_plot <- function(uhat_tbl,
                                       title = "Line of action: empirical vs geometric") {
  df <- dplyr::filter(uhat_tbl, is.finite(.data$uhat_angle_emp_deg), is.finite(.data$uhat_angle_geom_deg))
  if (nrow(df) == 0L) return(NULL)
  lims <- range(c(df$uhat_angle_emp_deg, df$uhat_angle_geom_deg), na.rm = TRUE)
  ref <- tibble::tibble(x = lims, y = lims)
  r <- suppressWarnings(stats::cor(df$uhat_angle_geom_deg, df$uhat_angle_emp_deg, use = "complete.obs"))

  ggplot(df, aes(x = .data$uhat_angle_geom_deg, y = .data$uhat_angle_emp_deg,
                 color = .data$muscle_side, size = .data$activation_snr)) +
    geom_line(data = ref, aes(x = .data$x, y = .data$y), inherit.aes = FALSE,
              linetype = "dashed", color = "black") +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = MFV_SIDE_COLORS, name = "Muscle side") +
    scale_size_continuous(name = "Activation SNR") +
    labs(title = title,
         subtitle = sprintf("r = %.3f, n = %d step(s); dashed = 1:1; point size = activation SNR; empirical primary, geometric = beam-tangent theta/2 cross-check (approx)", r, nrow(df)),
         x = "Geometric u_hat angle in XY (deg, beam tangent ~theta/2)",
         y = "Empirical u_hat angle in XY (deg, dF direction)") +
    theme_bw(base_size = 12)
}

#' Muscle force (vector, N) vs time relative to stim onset, pooled per category.
#' Thin per-step lines + bold binned mean +/- SE. Same layout as
#' plot_force_vs_time.R's builder but y-axis in Newtons (vector LOA force).
build_force_vs_time_vector_plot <- function(ts_df, title, facet_var = NULL,
                                            mean_bin_width_s = 0.02) {
  df <- dplyr::filter(ts_df, is.finite(.data$t_rel), is.finite(.data$muscle_force_vector_N))
  if (!is.null(facet_var)) df <- dplyr::filter(df, !is.na(.data[[facet_var]]))
  if (nrow(df) == 0L) return(NULL)
  dur_range <- range(df$stim_duration_s, na.rm = TRUE)

  mean_grp <- c(if (!is.null(facet_var)) facet_var, "t_bin")
  mean_line <- df |>
    dplyr::mutate(t_bin = round(.data$t_rel / mean_bin_width_s) * mean_bin_width_s) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(mean_grp))) |>
    dplyr::summarise(y_mean = mean(.data$muscle_force_vector_N, na.rm = TRUE),
                     y_se = stats::sd(.data$muscle_force_vector_N, na.rm = TRUE) / sqrt(dplyr::n()),
                     .groups = "drop") |>
    dplyr::mutate(y_se = dplyr::if_else(is.finite(.data$y_se), .data$y_se, 0))

  p <- ggplot(df, aes(x = .data$t_rel, y = .data$muscle_force_vector_N)) +
    annotate("rect", xmin = 0, xmax = dur_range[2L], ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.07) +
    annotate("rect", xmin = 0, xmax = dur_range[1L], ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.10) +
    geom_line(aes(group = .data$unit_id, color = .data$muscle_side), alpha = 0.35, linewidth = 0.4) +
    geom_ribbon(data = mean_line, aes(x = .data$t_bin, y = .data$y_mean,
                                      ymin = .data$y_mean - .data$y_se, ymax = .data$y_mean + .data$y_se),
                inherit.aes = FALSE, alpha = 0.25, fill = "black") +
    geom_line(data = mean_line, aes(x = .data$t_bin, y = .data$y_mean),
              inherit.aes = FALSE, color = "black", linewidth = 1.1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
    scale_color_manual(values = MFV_SIDE_COLORS, name = "Muscle side") +
    labs(title = title,
         subtitle = sprintf("n = %d step(s) across %d trial(s); orange = stim ON; bold = binned mean +/- SE; force along u_hat",
                            dplyr::n_distinct(df$unit_id), dplyr::n_distinct(df$trial_id)),
         x = "Time relative to stimulus onset (s)", y = "Muscle force along u_hat (N)") +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")

  if (!is.null(facet_var)) p <- p + facet_wrap(stats::as.formula(paste0("~", facet_var)),
                                               scales = "free_y", labeller = ggplot2::label_both)
  p
}
