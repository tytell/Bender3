# fit_fv_fl.R
# Force-Velocity (Hill hyperbola) curve fit for segmented_finite
# (isovelocity) per-muscle-side step summaries.
#
# Force-Length (isometric) summary plots do NOT fit a model by default (PI
# direction, 2026-07-16): a quadratic/parabolic fit to sparse, noisy FL data
# was found to invent force "peaks" (fitted vertices) that lay outside the
# tested strain range entirely -- a fit artifact, not a measurement. The FL
# default is now "connect the mean" (see .mean_line_by_side()/
# build_summary_plot_isometric(), plot_summary_profiles.R): plot the mean
# (+/- SE) at each tested strain level with a straight connecting line and no
# imposed functional form. If a monotonic FL model fit is genuinely needed,
# it must be requested explicitly by name, and any drawn curve/CI must never
# extrapolate past the tested strain range.
#
# The Hill Force-Velocity fit below is kept as the default FV overlay
# because the Hill model was explicitly requested by name at project start
# (fit_rigor: hill) -- an EXPLICITLY-named fit is allowed to be drawn by
# default. Its drawn curve is still clipped to the tested velocity range
# (see build_summary_plot_isovelocity()); Vmax/peak power are, by the model's
# own definition, extrapolated to F=0 and are flagged as such wherever
# reported.
#
# Every fit function returns status = "ok" or "failed" explicitly -- callers
# (and the final pipeline report) must surface "failed" fits, never silently
# drop them or report success.

library(dplyr)
library(broom)
library(minpack.lm)

.fit_min_points <- 4L

#' Hill Force-Velocity fit on the CONCENTRIC (shortening, x > 0) limb only --
#' the classical Hill hyperbola is a shortening-quadrant model and is not fit
#' to eccentric (lengthening) points, which are reported as raw data alongside
#' but excluded from the nonlinear fit (standard practice; documented here as
#' the FV model-scope decision).
#'
#' F(V) = (F0 + a) * b / (V + b) - a
#' Vmax  = b * F0 / a        (V at F = 0)
#' Peak power = max_{0<V<Vmax} F(V) * V, found numerically.
#'
#' @param v Signed shortening velocity/strain-rate (positive = concentric).
#' @param force Muscle force at that step.
#' @param f0_isometric Optional externally-measured isometric force (force at
#'   v ~ 0) to anchor F0 when there are too few points to fit F0 freely.
fit_force_velocity_curve <- function(v, force, side_label = NA_character_, f0_isometric = NA_real_) {
  ok <- is.finite(v) & is.finite(force) & v >= 0
  v <- v[ok]; force <- force[ok]
  n_pts <- dplyr::n_distinct(v)

  if (n_pts < .fit_min_points) {
    return(list(status = "failed",
                reason = sprintf("only %d distinct concentric velocities (need >= %d)", n_pts, .fit_min_points),
                side = side_label, n = length(v)))
  }

  F0_fixed <- is.finite(f0_isometric) && f0_isometric > 0
  F0_start <- if (F0_fixed) f0_isometric else max(force, na.rm = TRUE)
  if (!is.finite(F0_start) || F0_start <= 0) {
    return(list(status = "failed", reason = "no positive isometric force anchor available",
                side = side_label, n = length(v)))
  }

  vmax_guess <- max(v, na.rm = TRUE) * 2.0
  if (!is.finite(vmax_guess) || vmax_guess <= 0) vmax_guess <- 1.0

  fit_res <- if (F0_fixed) {
    # F0 anchored to the externally-measured isometric force (v ~ 0); passed
    # into the model frame as a fixed column so only a/b are estimated --
    # keeps the fit well-determined with as few as 4 concentric points.
    tryCatch(
      minpack.lm::nlsLM(
        force ~ (F0v + a) * b / (v + b) - a,
        data = data.frame(v = v, force = force, F0v = F0_start),
        start = list(a = 0.25 * F0_start, b = 0.25 * vmax_guess),
        lower = c(a = 1e-9, b = 1e-9),
        upper = c(a = 1e6, b = 1e6),
        control = minpack.lm::nls.lm.control(maxiter = 500L)
      ),
      error = function(e) e
    )
  } else {
    tryCatch(
      minpack.lm::nlsLM(
        force ~ (F0 + a) * b / (v + b) - a,
        data = data.frame(v = v, force = force),
        start = list(F0 = F0_start, a = 0.25 * F0_start, b = 0.25 * vmax_guess),
        lower = c(F0 = 1e-9, a = 1e-9, b = 1e-9),
        upper = c(F0 = 1e6, a = 1e6, b = 1e6),
        control = minpack.lm::nls.lm.control(maxiter = 500L)
      ),
      error = function(e) e
    )
  }

  if (inherits(fit_res, "error") || is.null(fit_res)) {
    reason <- if (inherits(fit_res, "error")) conditionMessage(fit_res) else "fit returned NULL"
    return(list(status = "failed", reason = paste("Hill nlsLM did not converge:", reason),
                side = side_label, n = length(v)))
  }

  co <- stats::coef(fit_res)
  a  <- co[["a"]]
  b  <- co[["b"]]
  F0 <- if (F0_fixed) F0_start else co[["F0"]]

  if (!all(is.finite(c(a, b, F0))) || a <= 0 || b <= 0 || F0 <= 0) {
    return(list(status = "failed", reason = "non-finite or non-physical fitted Hill parameters",
                side = side_label, n = length(v)))
  }

  Vmax <- b * F0 / a
  hill_force <- function(vv) (F0 + a) * b / (vv + b) - a
  power_fn   <- function(vv) hill_force(vv) * vv
  opt <- tryCatch(
    stats::optimize(power_fn, interval = c(1e-6, Vmax * 0.999), maximum = TRUE),
    error = function(e) NULL
  )
  if (is.null(opt) || !is.finite(opt$objective)) {
    return(list(status = "failed", reason = "peak-power optimization failed / non-finite",
                side = side_label, n = length(v), F0 = F0, a = a, b = b, Vmax = Vmax))
  }

  glance <- tryCatch(broom::glance(fit_res), error = function(e) NULL)
  tidy   <- tryCatch(broom::tidy(fit_res), error = function(e) NULL)

  list(
    status = "ok", side = side_label, n = length(v),
    F0 = F0, a = a, b = b, Vmax = Vmax,
    peak_power = opt$objective, velocity_at_peak_power = opt$maximum,
    F0_fixed = F0_fixed, tidy = tidy, glance = glance,
    predict = hill_force
  )
}
