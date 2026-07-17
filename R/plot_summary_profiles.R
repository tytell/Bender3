# plot_summary_profiles.R
# Required "Summary Profiles" visualization (task spec): ONE master summary
# plot per protocol category, showing all individual raw data points,
# mean +/- SE ribbon, and extracted properties (FL curve / FV curve / power
# curve, as applicable to that category). Pools across ALL trial files in
# the category -- these are cross-trial summaries, not per-trial plots.
#
# Category -> property mapping (PI-confirmed, 2026-07-15; FL default trend
# updated 2026-07-16):
#   isometric        -> Force-Length: connect-the-mean line (mean +/- SE per
#                       binned strain level, per muscle side) -- NO model fit
#                       by default. A parabolic fit to this sparse, noisy
#                       data previously invented force "peaks" (fitted
#                       vertices) that lay outside the tested strain range
#                       entirely -- a fit artifact, not a measurement. See
#                       fit_fv_fl.R module header.
#   isovelocity      -> Force-Velocity: same connect-the-mean default, PLUS
#                       the Hill fit + peak power overlay (Hill was
#                       EXPLICITLY requested by name at project start --
#                       fit_rigor: hill -- so it remains the one exception;
#                       its drawn curve is clipped to the tested velocity
#                       range, never extrapolated).
#   dynamic          -> per-cycle work/power vs. frequency (work-loop power curve;
#                       mean +/- SE bin line, no model fit)
#   frequency_sweep  -> passive stiffness/damping vs. frequency (no active fit --
#                       protocol is passive-only; see analyze_frequency_sweep.R)

library(dplyr)
library(ggplot2)

#' Bin x into `n_bins` equal-width bins and compute mean +/- SE of y per bin
#' (for the required "mean +/- SE ribbon"). Returns one row per non-empty bin.
bin_mean_se <- function(x, y, n_bins = 12L) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 2L) return(tibble::tibble(x_mid = numeric(0), y_mean = numeric(0), y_se = numeric(0), n = integer(0)))
  if (min(x) == max(x)) {
    # Degenerate range (e.g. only one distinct x value pooled, such as a
    # single-frequency dynamic corpus) -- seq()/cut() would otherwise produce
    # duplicate breaks. Collapse to one bin spanning all points.
    bin <- rep(1L, length(x))
  } else {
    brks <- seq(min(x), max(x), length.out = n_bins + 1L)
    brks[length(brks)] <- brks[length(brks)] + 1e-9
    bin <- cut(x, breaks = brks, include.lowest = TRUE, labels = FALSE)
  }
  tibble::tibble(x = x, y = y, bin = bin) |>
    dplyr::group_by(.data$bin) |>
    dplyr::summarise(
      x_mid  = mean(.data$x, na.rm = TRUE),
      y_mean = mean(.data$y, na.rm = TRUE),
      y_se   = stats::sd(.data$y, na.rm = TRUE) / sqrt(dplyr::n()),
      n      = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::filter(is.finite(.data$y_se) | .data$n == 1L) |>
    dplyr::mutate(y_se = dplyr::if_else(is.finite(.data$y_se), .data$y_se, 0))
}

#' Generic (species-nonspecific) vertebrate skeletal-muscle active
#' force-length SHAPE, for visual reference only -- NOT a fitted model and
#' NOT calibrated to any real force/pressure units (per-side y-scale below
#' is an arbitrary shape-comparison rescale, not a stress normalization).
#' No largemouth bass (or any fish) red-muscle FL curve exists in the
#' literature (Carroll & Wainwright 2006, on this exact muscle, state "The
#' force-length properties of this muscle are not known"; the nearest fish
#' red-muscle analog is a single point -- Rome 1998, common carp, ~20%
#' shortening -> ~15% force loss -- not a full curve). This uses the
#' classical Gordon/Huxley/Julian (1966) sarcomere active-FL relationship in
#' its common Hill-model Gaussian parametrization (Zajac 1989 / Thelen 2003):
#'   f(l_norm) = exp(-((l_norm - 1) / gamma)^2),  l_norm = 1 + strain_pct/100
#' gamma = 0.45 is a typical generic vertebrate skeletal-muscle plateau width
#' -- an assumption, not a bass-specific measurement.
generic_skeletal_muscle_FL_shape <- function(strain_pct, gamma = 0.45) {
  l_norm <- 1 + strain_pct / 100.0
  exp(-((l_norm - 1) / gamma)^2)
}

#' Rescale the generic FL shape onto an observed muscle_force_Nm range purely
#' for SHAPE comparison (min(shape)->min(observed), max(shape)->max(observed))
#' -- explicitly not a stress/pressure calibration.
.overlay_generic_FL_curve <- function(observed_force, x_range, gamma = 0.45, n = 200L) {
  xs <- seq(x_range[1L], x_range[2L], length.out = n)
  shape <- generic_skeletal_muscle_FL_shape(xs, gamma = gamma)
  rng <- range(observed_force, na.rm = TRUE)
  y <- rng[1L] + shape * diff(rng)
  tibble::tibble(x = xs, y = y)
}

#' DEFAULT trend representation for sparse FL/FV/summary x-y data (PI
#' direction, 2026-07-16): bin x WITHIN each side group and take bin_mean_se()
#' (mean +/- SE per bin, see its docstring above) -- explicitly NOT a
#' polynomial/spline/model fit. Pools ACROSS trial replicates at each x
#' value/bin (the earlier version fit a separate polynomial per trial; trial
#' replicates are now aggregated together into one mean line per side, per
#' the "connect the means" rule). If the underlying relationship is
#' monotonic over the tested range, the line is honestly monotonic -- no
#' peak/vertex is imposed that the data don't contain.
#'
#' n_bins is capped to the number of distinct (rounded) x values actually
#' present so a handful of discrete step levels aren't over-split into
#' near-empty bins, and floored at 4 so a single degenerate bin doesn't
#' collapse informative structure.
#' @return tibble(x_mid, y_mean, y_se, n, <side_col>) stacked across sides,
#'   or an empty tibble if no side has enough finite points.
.mean_line_by_side <- function(pts, x_col, y_col, side_col = "muscle_side", n_bins = 10L) {
  sides <- unique(pts[[side_col]])
  purrr::map_dfr(sides, function(sd) {
    sub <- pts[pts[[side_col]] == sd, , drop = FALSE]
    nb <- min(n_bins, max(4L, dplyr::n_distinct(round(sub[[x_col]], 1))))
    r <- bin_mean_se(sub[[x_col]], sub[[y_col]], n_bins = nb)
    if (nrow(r) == 0L) return(tibble::tibble())
    r[[side_col]] <- sd
    r
  })
}

#' Mass-/area-specific-property suffix (PI-directed, 2026-07-16) -- appends
#' specific tension (N/cm^2) and mass-specific peak power (W/kg) when
#' add_specific_properties_to_fit() (muscle_geometry.R) has populated them on
#' a Hill FV fit; empty string ("") when geometry was unavailable, so the
#' base label above is unchanged in that case. FV-only: FL no longer fits a
#' model by default (see .mean_line_by_side()), so there is no FL fit object
#' to annotate.
.fv_specific_property_suffix <- function(fit) {
  if (is.null(fit$specific_tension_Ncm2) || !is.finite(fit$specific_tension_Ncm2)) return("")
  if (!is.null(fit$mass_specific_peak_power_Wkg) && is.finite(fit$mass_specific_peak_power_Wkg)) {
    sprintf(", %.3g N/cm^2 specific tension, %.3g W/kg peak power", fit$specific_tension_Ncm2, fit$mass_specific_peak_power_Wkg)
  } else {
    sprintf(", %.3g N/cm^2 specific tension", fit$specific_tension_Ncm2)
  }
}

#' Hill Force-Velocity fit annotation. The Hill model was EXPLICITLY
#' requested by name at project start (fit_rigor: hill), so per the
#' "no unrequested fit" rule (PI direction, 2026-07-16) it remains the one
#' model overlay drawn by default on a summary plot -- but Vmax and peak
#' power are, by the Hill model's own definition, extrapolated to F=0 (a
#' velocity essentially never directly tested), so they are flagged
#' EXTRAPOLATED here rather than presented as directly observed measurements.
#' The drawn curve itself is clipped to the tested velocity range (see
#' build_summary_plot_isovelocity()), never extended out to Vmax.
.fv_fit_annotation_label <- function(fit) {
  if (is.null(fit) || identical(fit$status, "failed")) {
    return(sprintf("%s fit: FAILED (%s)", toupper(fit$side %||% "?"), fit$reason %||% "unknown reason"))
  }
  sprintf("%s Hill fit OK: F0=%.4g N*m, Vmax=%.2f%%/s [EXTRAPOLATED, not observed], Ppeak=%.4g, n=%d%s",
          toupper(fit$side), fit$F0, fit$Vmax, fit$peak_power, fit$n,
          .fv_specific_property_suffix(fit))
}
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Wrap a long subtitle string onto multiple lines so it doesn't get clipped
#' at the plot's right edge (ggplot2 doesn't auto-wrap plot.subtitle). Used
#' for the isovelocity summary subtitle, which got long once it started
#' combining the descriptive point count with the Hill fit annotation.
.wrap_subtitle <- function(s, width = 78L) paste(strwrap(s, width = width), collapse = "\n")

#' Force-Length summary: pools step_summary rows across all isometric trials.
#' DEFAULT trend is connect-the-mean (PI direction, 2026-07-16): NO model fit
#' is drawn or annotated -- see .mean_line_by_side()/fit_fv_fl.R module
#' header. A parabolic fit to this sparse, noisy data previously invented
#' force "peaks" (fitted vertices) outside the tested strain range; the
#' connect-the-mean line shows the data's true (possibly monotonic) shape
#' honestly instead.
#' @param step_summary Combined tibble (all trials, rbind of analyze_isometric()$step_summary,
#'   with a trial_id column added by the caller), filtered to muscle_side %in% c("left","right").
#' @param title Plot title override (default unchanged from the original
#'   static-baseline caption) -- lets callers reuse this builder unmodified
#'   for the "_interpbaseline" variant (2026-07-16) by passing a
#'   muscle_force_Nm-substituted step_summary plus a distinguishing title.
build_summary_plot_isometric <- function(step_summary, title = "Isometric summary: Force-Length") {
  side_colors <- c(left = "#1d4ed8", right = "#b91c1c")
  pts <- dplyr::filter(step_summary, .data$muscle_side %in% c("left", "right"),
                       is.finite(.data$muscle_force_Nm), is.finite(.data$shortening_strain_pct))

  # Connect-the-mean trend, pooled ACROSS trial replicates per side (see
  # .mean_line_by_side() docstring) -- the DEFAULT and only trend shown.
  trend <- .mean_line_by_side(pts, "shortening_strain_pct", "muscle_force_Nm")

  p <- ggplot(pts, aes(x = .data$shortening_strain_pct, y = .data$muscle_force_Nm, color = .data$muscle_side)) +
    geom_point(aes(shape = .data$trial_id), size = 2.2, alpha = 0.6) +
    { if (nrow(trend) > 0L)
        geom_ribbon(data = trend, aes(x = .data$x_mid, y = .data$y_mean,
                                      ymin = .data$y_mean - .data$y_se, ymax = .data$y_mean + .data$y_se,
                                      color = NULL, fill = .data$muscle_side, shape = NULL),
                    alpha = 0.2, inherit.aes = FALSE)
    } +
    { if (nrow(trend) > 0L)
        geom_line(data = trend, aes(x = .data$x_mid, y = .data$y_mean, color = .data$muscle_side, shape = NULL),
                  inherit.aes = FALSE, linewidth = 1.0)
    } +
    scale_color_manual(values = side_colors, name = "Muscle side") +
    scale_fill_manual(values = side_colors, guide = "none") +
    labs(title = title,
         subtitle = sprintf(
           "n = %d points pooled across %d trial(s); solid line = mean per binned strain level (NO model fit), band = +/-SE",
           nrow(pts), dplyr::n_distinct(pts$trial_id)),
         x = "Muscle shortening strain (%, folded by recruited side)", y = "Muscle force (N*m)",
         shape = "Trial") +
    theme_bw(base_size = 12)

  # Generic (non-bass-specific) skeletal-muscle FL SHAPE overlay -- see
  # generic_skeletal_muscle_FL_shape() docstring: no bass/fish red-muscle FL
  # curve exists in the literature, so this is a shape reference only, y-axis
  # rescaled per side purely for visual comparison (not a stress calibration).
  ref_df <- purrr::map_dfr(c("left", "right"), function(sd) {
    obs <- pts$muscle_force_Nm[pts$muscle_side == sd]
    if (length(obs) < 2L) return(NULL)
    rng <- range(pts$shortening_strain_pct[pts$muscle_side == sd], na.rm = TRUE)
    .overlay_generic_FL_curve(obs, rng) |> dplyr::mutate(muscle_side = sd)
  })
  if (nrow(ref_df) > 0L) {
    p <- p + geom_line(data = ref_df, aes(x = .data$x, y = .data$y, color = .data$muscle_side, shape = NULL),
                        inherit.aes = FALSE, linetype = "dotted", linewidth = 0.9, alpha = 0.6) +
      labs(caption = "Dotted = generic vertebrate skeletal-muscle FL SHAPE (Gordon-Huxley-Julian/Thelen Gaussian, gamma=0.45), y-rescaled for shape comparison only -- no bass/fish red-muscle FL curve exists in the literature (Carroll & Wainwright 2006)")
  }
  p
}

#' Force-Velocity + peak power summary: pools step_summary across all
#' isovelocity trials. DEFAULT trend is connect-the-mean (same as FL, PI
#' direction 2026-07-16); ADDITIONALLY overlays the Hill fit (dashed, clipped
#' to the tested velocity range) because the Hill model was EXPLICITLY
#' requested by name at project start -- see .fv_fit_annotation_label()
#' docstring.
build_summary_plot_isovelocity <- function(step_summary, fits) {
  side_colors <- c(left = "#1d4ed8", right = "#b91c1c")
  pts <- dplyr::filter(step_summary, .data$muscle_side %in% c("left", "right"),
                       is.finite(.data$muscle_force_Nm), is.finite(.data$shortening_strain_pct))
  pts_conc <- dplyr::filter(pts, .data$contraction_mode %in% c("concentric", "isometric_zero"))

  # Connect-the-mean trend on the FULL (concentric + eccentric) data, pooled
  # across trial replicates per side -- the DEFAULT trend (see
  # .mean_line_by_side()). The Hill fit below is drawn/annotated separately
  # and ADDITIONALLY, restricted to the concentric limb per its own scope.
  trend <- .mean_line_by_side(pts, "shortening_strain_pct", "muscle_force_Nm")

  p <- ggplot(pts, aes(x = .data$shortening_strain_pct, y = .data$muscle_force_Nm, color = .data$muscle_side)) +
    geom_point(aes(shape = .data$contraction_mode), size = 2.2, alpha = 0.6) +
    { if (nrow(trend) > 0L)
        geom_ribbon(data = trend, aes(x = .data$x_mid, y = .data$y_mean,
                                      ymin = .data$y_mean - .data$y_se, ymax = .data$y_mean + .data$y_se,
                                      color = NULL, fill = .data$muscle_side, shape = NULL),
                    alpha = 0.2, inherit.aes = FALSE)
    } +
    { if (nrow(trend) > 0L)
        geom_line(data = trend, aes(x = .data$x_mid, y = .data$y_mean, color = .data$muscle_side, shape = NULL),
                  inherit.aes = FALSE, linewidth = 1.0)
    } +
    scale_color_manual(values = side_colors, name = "Muscle side") +
    scale_fill_manual(values = side_colors, guide = "none") +
    labs(title = "Isovelocity summary: Force-Velocity",
         subtitle = .wrap_subtitle(paste0(
           sprintf("n = %d points pooled across %d trial(s); solid = mean per bin (NO model fit), band = +/-SE  |  ",
                   nrow(pts), dplyr::n_distinct(pts$trial_id)),
           paste(sapply(fits, .fv_fit_annotation_label), collapse = "  |  "),
           "  |  dashed = Hill fit (explicitly requested by name), clipped to tested velocity range"
         )),
         x = "Muscle shortening strain-rate (%/s, folded by recruited side; concentric > 0)",
         y = "Muscle force (N*m)", shape = "Contraction mode") +
    theme_bw(base_size = 12)

  for (sd in c("left", "right")) {
    f <- fits[[sd]]
    if (!is.null(f) && identical(f$status, "ok")) {
      xs <- seq(1e-6, max(pts_conc$shortening_strain_pct[pts_conc$muscle_side == sd], na.rm = TRUE), length.out = 200)
      curve_df <- tibble::tibble(x = xs, y = f$predict(xs), muscle_side = sd)
      p <- p + geom_line(data = curve_df, aes(x = .data$x, y = .data$y, color = .data$muscle_side, shape = NULL),
                          inherit.aes = FALSE, linetype = "dashed", linewidth = 1.0)
    }
  }
  p
}

#' Dynamic summary: per-cycle work-loop power vs. cycling frequency, pooled
#' across all active (non-passive-only) dynamic trials.
#' @param cycle_summary Combined summarize_muscle_cycles() output across trials
#'   (trial_id column added by caller).
build_summary_plot_dynamic <- function(cycle_summary) {
  cs <- dplyr::filter(cycle_summary, is.finite(.data$freq.Hz), is.finite(.data$avg_power.W))
  ribbon <- bin_mean_se(cs$freq.Hz, cs$avg_power.W, n_bins = min(12L, max(3L, dplyr::n_distinct(cs$freq.Hz))))

  ggplot(cs, aes(x = .data$freq.Hz, y = .data$avg_power.W)) +
    geom_point(aes(color = .data$trial_id), size = 2.0, alpha = 0.65) +
    { if (nrow(ribbon) > 0L)
        geom_ribbon(data = ribbon, aes(x = .data$x_mid, y = .data$y_mean,
                                       ymin = .data$y_mean - .data$y_se, ymax = .data$y_mean + .data$y_se),
                    inherit.aes = FALSE, alpha = 0.25, fill = "#1d4ed8")
    } +
    { if (nrow(ribbon) > 0L)
        geom_line(data = ribbon, aes(x = .data$x_mid, y = .data$y_mean),
                  inherit.aes = FALSE, color = "#1d4ed8", linewidth = 1.0)
    } +
    labs(title = "Dynamic (work-loop) summary: mean power vs. cycling frequency",
         subtitle = sprintf("n = %d cycles pooled across %d trial(s); ribbon = mean +/- SE per frequency bin",
                            nrow(cs), dplyr::n_distinct(cs$trial_id)),
         x = "Cycling frequency (Hz)", y = "Mean muscle power (W)", color = "Trial") +
    theme_bw(base_size = 12)
}

#' Mass-specific version of build_summary_plot_dynamic() (PI-directed,
#' 2026-07-16, see compute_muscle_mass_and_csa()/muscle_geometry.R): plots
#' avg_power.Wkg (now populated once muscle_mass.kg is attached -- see
#' run_fv_fl_power_pipeline.R's dynamic branch) instead of raw avg_power.W,
#' with a shaded reference band for Coughlin et al. 1996 scup red-muscle
#' steady-state power (coughlin_steady_state_power_limits(), 03_analyze.R)
#' so weak/strong twitches are visible at a glance against a literature
#' benchmark. ADDITIONAL to build_summary_plot_dynamic() (which is
#' unchanged) -- saved as a separate "_massspecific" file.
build_summary_plot_dynamic_massspecific <- function(cycle_summary) {
  cs <- dplyr::filter(cycle_summary, is.finite(.data$freq.Hz), is.finite(.data$avg_power.Wkg))
  if (nrow(cs) == 0L) return(NULL)
  ribbon <- bin_mean_se(cs$freq.Hz, cs$avg_power.Wkg, n_bins = min(12L, max(3L, dplyr::n_distinct(cs$freq.Hz))))
  lims <- coughlin_steady_state_power_limits()
  xr <- range(cs$freq.Hz, na.rm = TRUE)

  ggplot(cs, aes(x = .data$freq.Hz, y = .data$avg_power.Wkg)) +
    annotate("rect", xmin = xr[1L] - 0.5, xmax = xr[2L] + 0.5, ymin = lims$lo, ymax = lims$hi,
             fill = "grey40", alpha = 0.15) +
    geom_hline(yintercept = lims$mean, linetype = "dashed", color = "grey40") +
    geom_point(aes(color = .data$trial_id), size = 2.0, alpha = 0.65) +
    { if (nrow(ribbon) > 0L)
        geom_ribbon(data = ribbon, aes(x = .data$x_mid, y = .data$y_mean,
                                       ymin = .data$y_mean - .data$y_se, ymax = .data$y_mean + .data$y_se),
                    inherit.aes = FALSE, alpha = 0.25, fill = "#1d4ed8")
    } +
    { if (nrow(ribbon) > 0L)
        geom_line(data = ribbon, aes(x = .data$x_mid, y = .data$y_mean),
                  inherit.aes = FALSE, color = "#1d4ed8", linewidth = 1.0)
    } +
    labs(title = "Dynamic (work-loop) summary: mass-specific power vs. cycling frequency",
         subtitle = sprintf("n = %d cycles pooled across %d trial(s); grey band = %s (%.0f-%.0f W/kg)",
                            nrow(cs), dplyr::n_distinct(cs$trial_id), lims$reference, lims$lo, lims$hi),
         x = "Cycling frequency (Hz)", y = "Mean muscle power (W/kg)", color = "Trial") +
    theme_bw(base_size = 12)
}

#' Frequency-sweep summary: passive stiffness (EI1) and damping (etaI1) vs.
#' instantaneous frequency, pooled across all frequency_sweep trials.
#' @param cycle_results Combined analyze_frequency_sweep() output across
#'   trials (trial_id column added by caller), status == "ok" rows only.
build_summary_plot_frequency_sweep <- function(cycle_results) {
  cr <- dplyr::filter(cycle_results, .data$status == "ok",
                      is.finite(.data$frequency_hz), is.finite(.data$EI1.Nm2))
  rib_k <- bin_mean_se(cr$frequency_hz, cr$EI1.Nm2, n_bins = 14L)
  rib_d <- bin_mean_se(cr$frequency_hz, cr$etaI1.Nm2s, n_bins = 14L)

  p_k <- ggplot(cr, aes(x = .data$frequency_hz, y = .data$EI1.Nm2)) +
    geom_point(aes(color = .data$trial_id), size = 1.6, alpha = 0.5) +
    geom_ribbon(data = rib_k, aes(x = .data$x_mid, y = .data$y_mean,
                                  ymin = .data$y_mean - .data$y_se, ymax = .data$y_mean + .data$y_se),
                inherit.aes = FALSE, alpha = 0.25, fill = "#166534") +
    geom_line(data = rib_k, aes(x = .data$x_mid, y = .data$y_mean), inherit.aes = FALSE,
              color = "#166534", linewidth = 1.0) +
    labs(y = "Passive stiffness EI (N*m^2)", x = NULL, color = "Trial") +
    theme_bw(base_size = 12)

  p_d <- ggplot(cr, aes(x = .data$frequency_hz, y = .data$etaI1.Nm2s)) +
    geom_point(aes(color = .data$trial_id), size = 1.6, alpha = 0.5) +
    geom_ribbon(data = rib_d, aes(x = .data$x_mid, y = .data$y_mean,
                                  ymin = .data$y_mean - .data$y_se, ymax = .data$y_mean + .data$y_se),
                inherit.aes = FALSE, alpha = 0.25, fill = "#b45309") +
    geom_line(data = rib_d, aes(x = .data$x_mid, y = .data$y_mean), inherit.aes = FALSE,
              color = "#b45309", linewidth = 1.0) +
    labs(y = "Passive damping etaI (N*m^2*s)", x = "Frequency (Hz)", color = "Trial") +
    theme_bw(base_size = 12)

  patchwork::wrap_plots(p_k, p_d, ncol = 1) +
    patchwork::plot_annotation(
      title = "Frequency-sweep summary: passive stiffness & damping vs. frequency",
      subtitle = sprintf("n = %d cycles (status=ok) pooled across %d trial(s) -- passive-only protocol, no active FV/FL/power fit",
                         nrow(cr), dplyr::n_distinct(cr$trial_id))
    )
}
