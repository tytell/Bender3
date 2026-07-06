# 03_analyze.R
# Muscle / mechanics analysis pipeline for Bender flat-schema HDF5 files.
#
# Ported from the validated scup corpus:
#   ~/Desktop/bender_projects/bender_scup_muscle/R/
#   ~/Desktop/bender_projects/ScupMechanics/R/
#
# single_finite (dynamic, sweep, constant-freq, combo): fully ported.
# segmented_finite (isometric, isovelocity, FL, FV): stub entry points — FL curve
#   fitting and FV power estimation are deferred to a follow-on session.
#
# Dependencies: dplyr, tibble, pracma, signal (optional LP filter), cli
#
# Typical pipeline for single_finite:
#   source("R/load_bender_flat.R")
#   source("R/03_analyze.R")
#   td  <- load_bender_flat("trial.h5")
#   td  <- set_cycle_types(td)
#   msc <- calc_muscle_torque(td)
#   cyc <- summarize_muscle_cycles(msc)

library(dplyr)
library(tibble)
library(pracma)
library(rlang)


# =============================================================================
# set_cycle_types
# Label whole cycles (act/pass) and half-cycles for muscle analysis.
# =============================================================================

#' Label bending cycles and half-cycles for muscle analysis.
set_cycle_types <- function(df) {
  if ("is_active_by_cycle" %in% names(df)) {
    df1 <- df |>
      dplyr::group_by(dplyr::across(dplyr::any_of(c("fishcode", "trial", "cycle")))) |>
      dplyr::mutate(
        cycletype = dplyr::if_else(
          dplyr::first(.data$is_active_by_cycle %in% TRUE),
          "act", "pass"
        )
      ) |>
      dplyr::ungroup()
  } else {
    df1 <- df |>
      dplyr::group_by(dplyr::across(dplyr::any_of(c("fishcode", "trial", "cycle")))) |>
      dplyr::mutate(cycletype = dplyr::if_else(any(.data$stim != "0"), "act", "pass")) |>
      dplyr::ungroup()
  }

  halfcycdata <- df1 |>
    dplyr::group_by(dplyr::across(dplyr::any_of(c("fishcode", "trial", "halfcycle")))) |>
    dplyr::summarize(
      halfcycletype = dplyr::if_else(any(.data$stim != "0"), "act", "pass"),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      halfcycleorder = dplyr::case_when(
        .data$halfcycletype == "act" & dplyr::lag(.data$halfcycletype) == "pass" ~ "first",
        .data$halfcycletype == "act" & dplyr::lead(.data$halfcycletype) == "pass" ~ "last",
        .data$halfcycletype == "pass" & dplyr::lag(.data$halfcycletype) == "act" ~ "first",
        .default = "middle"
      )
    )

  dplyr::left_join(
    df1, halfcycdata,
    by = intersect(c("fishcode", "trial", "halfcycle"), names(df1))
  )
}


# =============================================================================
# detect_trial_type / detect_muscle_group_keys
# =============================================================================

#' Classify trial structure for downstream grouping.
detect_trial_type <- function(data) {
  stopifnot(is.data.frame(data), nrow(data) > 0L)
  if (!"stimulus_type" %in% names(data)) {
    return(list(stimulus_type = NA_character_, layout = "unknown",
                group_keys = character(0), n_freq = NA_integer_,
                n_amp = NA_integer_, n_duty = NA_integer_, n_phase = NA_integer_))
  }
  stimulus_type <- unique(stats::na.omit(data$stimulus_type))[1L]
  n_freq  <- dplyr::n_distinct(data$freq.Hz, na.rm = TRUE)
  n_amp   <- if ("amp.deg" %in% names(data)) dplyr::n_distinct(data$amp.deg, na.rm = TRUE) else NA_integer_
  n_duty  <- dplyr::n_distinct(data$duty, na.rm = TRUE)
  n_phase <- dplyr::n_distinct(data$phase, na.rm = TRUE)
  layout <- dplyr::case_when(
    identical(stimulus_type, "constant_frequency") & n_freq <= 1L & n_duty <= 1L & n_phase <= 1L ~ "single_condition",
    identical(stimulus_type, "constant_frequency") ~ "constant_multi_param",
    identical(stimulus_type, "frequency_amplitude_combo") ~ "freq_amp_combo",
    identical(stimulus_type, "frequency_sweep") ~ "frequency_sweep",
    .default = "unknown"
  )
  group_keys <- detect_muscle_group_keys(data, stimulus_type = stimulus_type, layout = layout)
  list(stimulus_type = stimulus_type, layout = layout, group_keys = group_keys,
       n_freq = n_freq, n_amp = n_amp, n_duty = n_duty, n_phase = n_phase)
}

#' Grouping keys for passive torque templates by trial layout.
detect_muscle_group_keys <- function(data, stimulus_type = NULL, layout = NULL) {
  if (is.null(layout)) layout <- detect_trial_type(data)$layout
  keys <- if (layout == "freq_amp_combo") {
    curve_col <- if ("curvature.invm" %in% names(data)) "curvature.invm" else "curve.invm"
    c("freq.Hz", curve_col)
  } else if (layout == "constant_multi_param") {
    k <- character(0)
    if (dplyr::n_distinct(data$duty, na.rm = TRUE) > 1L) k <- c(k, "duty")
    if (dplyr::n_distinct(data$phase, na.rm = TRUE) > 1L) k <- c(k, "phase")
    k
  } else {
    character(0)
  }
  keys[keys %in% names(data)]
}


# =============================================================================
# Muscle mass
# =============================================================================

#' Estimate red-muscle CSA and mass from fish TL.
#' Coefficients from scupBender/data_muscle.Rmd (CSA regression on TL).
estimate_muscle_mass <- function(fishlength.m, dclamp.m) {
  muscle_csa.m2 <- (.001023 * fishlength.m^2) + (.0001249 * fishlength.m) + 9.043e-7
  muscle_mass.kg <- muscle_csa.m2 * dclamp.m * 1060.0
  list(muscle_csa.m2 = muscle_csa.m2, muscle_mass.kg = muscle_mass.kg)
}

#' Add muscle_csa.m2 and muscle_mass.kg columns.
add_muscle_mass <- function(df, length_col = "fishlength.m", clamp_col = "dclamp.m") {
  if (!all(c(length_col, clamp_col) %in% names(df))) {
    cli::cli_abort("Need {length_col} and {clamp_col} columns")
  }
  mm <- Map(estimate_muscle_mass,
            fishlength.m = df[[length_col]],
            dclamp.m     = df[[clamp_col]])
  df$muscle_csa.m2  <- vapply(mm, `[[`, numeric(1L), "muscle_csa.m2")
  df$muscle_mass.kg <- vapply(mm, `[[`, numeric(1L), "muscle_mass.kg")
  df
}


# =============================================================================
# Curvature helpers
# =============================================================================

#' Ensure nominal curvature.invm is present; backfill from amplitude or peak.
add_nominal_curvature_if_missing <- function(df) {
  if ("curvature.invm" %in% names(df) && any(is.finite(df$curvature.invm))) return(df)
  if (!("curve.invm" %in% names(df))) return(df)
  grp <- intersect(c("trial", "filename", "cycle"), names(df))
  has_amp <- "curveamp.invm" %in% names(df)
  df |>
    dplyr::group_by(dplyr::across(dplyr::any_of(grp))) |>
    dplyr::mutate(curvature.invm = abs(dplyr::coalesce(
      if (has_amp) .data$curveamp.invm else NA_real_,
      if ("amp.deg" %in% names(df) && "dclamp.m" %in% names(df))
        .data$amp.deg * pi / 180.0 / dplyr::first(.data$dclamp.m)
      else NA_real_,
      suppressWarnings(max(abs(.data$curve.invm), na.rm = TRUE))
    ))) |>
    dplyr::ungroup()
}


# =============================================================================
# Filter helpers
# =============================================================================

filter_cycles_if_set <- function(df, cycles_keep) {
  if (is.null(cycles_keep)) df else dplyr::filter(df, .data$cycle %in% cycles_keep)
}

#' Drop pass/act transition half-cycles (ramp up/down at block edges).
filter_transition_halfcycles <- function(df, require_act = TRUE, warn = TRUE) {
  if (!"halfcycleorder" %in% names(df)) {
    if (isTRUE(warn)) cli::cli_warn("halfcycleorder column missing; transition half-cycle filter skipped")
    return(df)
  }
  out <- dplyr::filter(df, .data$halfcycleorder == "middle")
  if (isTRUE(require_act) && "halfcycletype" %in% names(out)) {
    out <- dplyr::filter(out, .data$halfcycletype == "act")
  }
  out
}

#' Filter by nominal tailbeat frequency and curvature.
filter_poster_kinematics <- function(df, freq_target = 1, curve_target = 4,
                                     freq_tol = 0.15, curve_tol = 0.5) {
  df <- add_nominal_curvature_if_missing(df)
  out <- df |> dplyr::filter(abs(.data$freq.Hz - freq_target) <= freq_tol)
  if (!is.null(curve_target) && !is.na(curve_target)) {
    out <- out |> dplyr::filter(abs(.data$curvature.invm - curve_target) <= curve_tol)
  }
  out
}

#' Filter by activation parameters (phase and duty).
filter_poster_activation <- function(df, phase_target = -0.13, duty_target = 0.30, tol = 0.02) {
  df |>
    dplyr::filter(
      abs(.data$phase - phase_target) <= tol,
      abs(.data$duty  - duty_target)  <= tol
    )
}


# =============================================================================
# calc_work
# =============================================================================

#' Mechanical work via trapezoidal integration.
#' angle must be in radians.
calc_work <- function(angle, torque) {
  if (max(abs(angle), na.rm = TRUE) >= 2 * pi) {
    stop("angle range seems high -- check units (should be radians, not degrees)")
  }
  pracma::trapz(angle, torque)
}


# =============================================================================
# Stiffness / damping (ported from ScupMechanics/R/stiffness_damping.R)
# =============================================================================

.avg_mech <- function(x, torque) {
  pc <- prcomp(matrix(c(x, torque), ncol = 2L))
  as.numeric(pc$rotation[2L, 1L] / pc$rotation[1L, 1L])
}

calc_avg_stiffness  <- function(curve, torque)      .avg_mech(curve, torque)
calc_avg_damping    <- function(curverate, torque)   .avg_mech(curverate, torque)

.high_mech <- function(x, torque, sign = 1, qhigh = 0.99) {
  c1 <- sign * x; t1 <- sign * torque
  ishigh <- c1 >= quantile(c1, qhigh)
  median(t1[ishigh], na.rm = TRUE) / median(c1[ishigh], na.rm = TRUE)
}

calc_high_curve_stiffness    <- function(curve, torque, ...)     .high_mech(curve, torque, ...)
calc_high_curverate_damping  <- function(curverate, torque, ...) .high_mech(curverate, torque, ...)

.low_mech <- function(t, x, torque, sign = 1, qlow = 0.08) {
  c1 <- sign * x; t1 <- sign * torque
  d <- quantile(abs(c1), qlow)
  trange <- range(t); dur <- trange[2L] - trange[1L]
  islow <- (abs(c1) < d) & dplyr::between(t, trange[1L] + 0.25 * dur, trange[2L] - 0.25 * dur)
  if (!any(islow) || all(!is.finite(c1[islow]))) return(NA_real_)
  X <- matrix(c(c1[islow], rep_len(1, sum(islow))), ncol = 2L)
  as.numeric(lm.fit(X, t1[islow])$coefficients[1L])
}

calc_low_curve_stiffness   <- function(t, curve, torque, ...)    .low_mech(t, curve, torque, ...)
calc_low_curverate_damping <- function(t, curverate, torque, ...) .low_mech(t, curverate, torque, ...)


# =============================================================================
# get_mechanics_by_half_cycle (ported from ScupMechanics)
# =============================================================================

#' Bending mechanics (stiffness, damping, work) per half-cycle.
get_mechanics_by_half_cycle <- function(df) {
  empty_mech <- function(df) {
    df |>
      dplyr::slice_head(n = 1L) |>
      dplyr::mutate(sgn = NA_real_,
                    EI1.Nm2 = NA_real_, EIL.Nm2 = NA_real_, EIM.Nm2 = NA_real_,
                    etaI1.Nm2s = NA_real_, etaIL.Nm2s = NA_real_, etaIM.Nm2s = NA_real_,
                    work.Nm = NA_real_)
  }

  tq_col <- intersect(c("bodytorque.Nm", "xtorque.Nm", "muscle_torque.Nm"), names(df))[1L]
  if (is.na(tq_col) || all(is.na(df[[tq_col]]))) return(empty_mech(df))

  mech <- df |>
    dplyr::filter(!is.na(.data$halfcycle)) |>
    dplyr::group_by(.data$halfcycle) |>
    dplyr::summarize(
      sgn = sign(median(.data$angle.deg)),
      dplyr::across(
        -dplyr::any_of(c("angle.deg", "anglevel.degps", "curve.invm", "curverate.invms",
                         "xtorque.Nm", "ytorque.Nm", "ztorque.Nm", "bodytorque.Nm",
                         "xforce.N", "yforce.N", "zforce.N", "muscle_torque.Nm")),
        dplyr::first
      ),
      .groups = "drop"
    ) |>
    dplyr::mutate(EI1.Nm2 = NA_real_, EIL.Nm2 = NA_real_, EIM.Nm2 = NA_real_,
                  etaI1.Nm2s = NA_real_, etaIL.Nm2s = NA_real_, etaIM.Nm2s = NA_real_,
                  work.Nm = NA_real_)

  for (i in seq(1L, nrow(mech) - 1L)) {
    if (!all(c("halfcycleorder", "halfcycletype") %in% names(mech))) next
    if (all(mech$halfcycleorder[i:(i + 1L)] != "first") &&
        mech$halfcycletype[i] == mech$halfcycletype[i + 1L]) {
      dcyc <- df |>
        dplyr::filter(.data$halfcycle >= mech$halfcycle[i],
                      .data$halfcycle <  mech$halfcycle[i] + 1.0)
      if (all(is.finite(dcyc[[tq_col]])) && all(is.finite(dcyc$curve.invm)) &&
          all(is.finite(dcyc$curverate.invms))) {
        tq <- dcyc[[tq_col]]; cv <- dcyc$curve.invm; cr <- dcyc$curverate.invms
        mech$EI1.Nm2[i]    <- calc_avg_stiffness(cv, tq)
        mech$EIM.Nm2[i]    <- calc_high_curve_stiffness(cv, tq, sign = mech$sgn[i])
        mech$EIL.Nm2[i]    <- calc_low_curve_stiffness(dcyc$t.s, cv, tq, sign = mech$sgn[i])
        mech$etaI1.Nm2s[i] <- calc_avg_damping(cr, tq)
        mech$etaIM.Nm2s[i] <- calc_high_curverate_damping(cr, tq, sign = mech$sgn[i])
        mech$etaIL.Nm2s[i] <- calc_low_curverate_damping(dcyc$t.s, cr, tq, sign = mech$sgn[i])
        mech$work.Nm[i]    <- suppressWarnings(
          calc_work(dcyc$angle.deg * pi / 180.0, tq)
        )
      }
    }
  }
  mech
}


# =============================================================================
# interpolate_even_phase (ported from ScupMechanics)
# =============================================================================

#' Interpolate time series to evenly spaced normalized phase grid.
interpolate_even_phase <- function(df, tnorm, tsec, vars, npercyc = 64L,
                                   .halfcycle = halfcycle) {
  tnorm_q    <- rlang::enquo(tnorm)
  halfcyc_q  <- rlang::enquo(.halfcycle)

  beforedata <- df |> dplyr::ungroup() |> dplyr::filter(t.s <= 0)
  afterdata  <- df |> dplyr::ungroup() |> dplyr::filter(t.s > 0 & t.norm == 0)
  df         <- df |> dplyr::ungroup() |> dplyr::filter((t.s == 0) | (t.s > 0 & t.norm != 0))

  ncyc     <- floor(max(dplyr::pull(df, !!tnorm_q), na.rm = TRUE))
  tnormeven <- seq(0, ncyc, by = 1.0 / npercyc)
  halfcyc   <- floor(tnormeven * 2.0) / 2.0

  constcols <- df |>
    dplyr::group_by(!!halfcyc_q) |>
    dplyr::summarize(dplyr::across(dplyr::everything(), \(x) all(x == dplyr::first(x, na_rm = TRUE), na.rm = TRUE))) |>
    dplyr::ungroup() |>
    dplyr::summarize(dplyr::across(dplyr::everything(), all)) |>
    dplyr::select(where(isTRUE)) |>
    names()

  cycdata <- df |>
    dplyr::ungroup() |>
    dplyr::select(!!halfcyc_q, dplyr::any_of(constcols)) |>
    dplyr::group_by(!!halfcyc_q) |>
    dplyr::summarize(dplyr::across(dplyr::everything(), dplyr::first))

  interp1 <- function(x, v, xout) {
    good <- !is.na(v) & !is.na(x)
    span <- xout >= min(x, na.rm = TRUE) & xout <= max(x, na.rm = TRUE)
    vs   <- rep_len(NA_real_, length(xout))
    if (sum(good) >= 2L) {
      fn  <- splinefun(x[good], v[good])
      vs[span] <- fn(xout[span])
    }
    vs
  }

  tnorm_vec <- dplyr::pull(df, !!tnorm_q)
  t_out <- interp1(tnorm_vec, df[[rlang::as_name(rlang::enquo(tsec))]], tnormeven)
  var_names <- tidyselect::eval_select(rlang::enquo(vars), df) |> names()
  interp_df <- tibble::tibble("{rlang::quo_name(tnorm_q)}" := tnormeven)
  interp_df[[rlang::as_name(rlang::enquo(tsec))]] <- t_out
  for (v in var_names) {
    interp_df[[v]] <- interp1(tnorm_vec, df[[v]], tnormeven)
  }
  interp_df[[rlang::quo_name(halfcyc_q)]] <- halfcyc

  interp_df <- dplyr::left_join(interp_df, cycdata,
                                 by = rlang::quo_name(halfcyc_q))
  all_cols <- names(interp_df)
  dplyr::bind_rows(
    dplyr::select(beforedata, dplyr::any_of(all_cols)),
    interp_df,
    dplyr::select(afterdata,  dplyr::any_of(all_cols))
  )
}


# =============================================================================
# estimate_body_torque (ported from bender_scup_muscle)
# =============================================================================

#' Select and rename body torque channel from calibrated F/T tibble.
estimate_body_torque <- function(df, method = "xtorque") {
  if (!method %in% c("xtorque", "ztorque", "yforce", "none")) {
    stop("method must be one of: xtorque, ztorque, yforce, none")
  }
  if ("xtorque.Nm" %in% names(df)) df <- df |> dplyr::mutate(bodytorque_xt.Nm = -.data$xtorque.Nm)
  if ("ztorque.Nm" %in% names(df)) df <- df |> dplyr::mutate(bodytorque_zt.Nm = .data$ztorque.Nm)
  if ("yforce.N"   %in% names(df) && "dclamp.m" %in% names(df)) {
    df <- df |> dplyr::mutate(bodytorque_yf.Nm = .data$yforce.N * .data$dclamp.m)
  }
  if ("bodytorque.Nm" %in% names(df)) {
    cli::cli_alert_warning("Overwriting existing bodytorque.Nm")
    df <- df |> dplyr::select(-"bodytorque.Nm")
  }
  if (method == "xtorque" && "bodytorque_xt.Nm" %in% names(df)) {
    df <- df |> dplyr::rename(bodytorque.Nm = bodytorque_xt.Nm)
  } else if (method == "ztorque" && "bodytorque_zt.Nm" %in% names(df)) {
    df <- df |> dplyr::rename(bodytorque.Nm = bodytorque_zt.Nm)
  } else if (method == "yforce" && "bodytorque_yf.Nm" %in% names(df)) {
    df <- df |> dplyr::rename(bodytorque.Nm = bodytorque_yf.Nm)
  }
  df
}


# =============================================================================
# detect_passive_mode / calc_muscle_torque
# =============================================================================

detect_passive_mode <- function(data, cycles_keep = NULL) {
  bend <- dplyr::filter(data, .data$cycle > 0, .data$t.s > 0)
  if (!is.null(cycles_keep)) bend <- dplyr::filter(bend, .data$cycle %in% cycles_keep)
  if (nrow(bend) == 0L) return("none")
  if ("cycletype" %in% names(bend)) {
    if (any(bend$cycletype == "pass", na.rm = TRUE) && any(bend$cycletype == "act", na.rm = TRUE)) {
      return("halfcycle")
    }
  }
  if ("halfcycletype" %in% names(bend)) {
    if (any(bend$halfcycletype == "pass", na.rm = TRUE) && any(bend$halfcycletype == "act", na.rm = TRUE)) {
      return("halfcycle")
    }
  }
  if (any(as.character(bend$stim) == "0", na.rm = TRUE) &&
      any(as.character(bend$stim) != "0", na.rm = TRUE)) return("stim_off")
  "none"
}


#' Muscle torque via active minus passive subtraction.
calc_muscle_torque <- function(data, torque_col = NULL, sample_freq.Hz = NULL,
                               group_keys = NULL, cycles_keep = NULL,
                               passive_mode = "auto",
                               include_all_active_samples = FALSE) {
  if (is.null(torque_col)) {
    torque_col <- intersect(
      c("xtorque.Nm", "xtorque0.Nm", "bodytorque.Nm"),
      names(data)
    )[1L]
  }
  if (is.na(torque_col) || !torque_col %in% names(data)) {
    stop("Torque column not found; available: ",
         paste(intersect(names(data), c("xtorque.Nm", "xtorque0.Nm", "bodytorque.Nm")), collapse = ", "))
  }
  if (is.null(sample_freq.Hz)) {
    sample_freq.Hz <- attr(data, "SampleFrequency.Hz")
    if (is.null(sample_freq.Hz)) sample_freq.Hz <- 1000.0
  }

  data <- add_nominal_curvature_if_missing(data)
  if (is.null(group_keys)) group_keys <- detect_muscle_group_keys(data)
  if (passive_mode == "auto") passive_mode <- detect_passive_mode(data, cycles_keep)

  phase_bin_width <- 0.02
  df <- data |>
    dplyr::group_by(.data$cycle) |>
    dplyr::filter(.data$cycle > 0) |>
    dplyr::mutate(
      t.interval.s = 1.0 / sample_freq.Hz,
      t.perc       = t.interval.s * (0:(dplyr::n() - 1L)),
      phase_frac   = (dplyr::row_number() - 1L) / pmax(dplyr::n() - 1L, 1L),
      phase_bin    = cut(.data$phase_frac, breaks = seq(0, 1, by = phase_bin_width),
                         include.lowest = TRUE, labels = FALSE),
      pos.rad      = .data$angle.deg * (pi / 180.0),
      dist.rad     = .data$pos.rad - dplyr::lag(.data$pos.rad),
      .torque      = .data[[torque_col]]
    ) |>
    dplyr::ungroup()

  phase_col  <- if (passive_mode == "stim_off") "phase_bin" else "t.perc"
  pass_grp   <- unique(c(group_keys, intersect(c("phase", "duty"), names(df)), phase_col))
  act_grp    <- c(pass_grp, "stim", "cycle")

  if (passive_mode == "halfcycle") {
    pass_rows <- dplyr::filter(df, .data$cycletype == "pass")
    act_rows  <- dplyr::filter(df, .data$cycletype == "act")
  } else if (passive_mode == "stim_off") {
    pass_rows <- dplyr::filter(df, as.character(.data$stim) == "0")
    act_rows  <- df
  } else {
    pass_rows <- df[0L, ]; act_rows <- df[0L, ]
  }
  if (!is.null(cycles_keep)) act_rows <- filter_cycles_if_set(act_rows, cycles_keep)

  passtorque <- pass_rows |>
    dplyr::group_by(dplyr::across(dplyr::all_of(pass_grp))) |>
    dplyr::summarise(avg_pass_torque.Nm = mean(.data$.torque, na.rm = TRUE), .groups = "drop")

  if (isTRUE(include_all_active_samples)) {
    join_keys <- intersect(pass_grp, names(act_rows))
    return(act_rows |>
      dplyr::left_join(passtorque, by = join_keys) |>
      dplyr::mutate(
        muscle_torque.Nm = .data$.torque - .data$avg_pass_torque.Nm,
        pass_torque.Nm   = .data$avg_pass_torque.Nm,
        t.interval.s     = 1.0 / sample_freq.Hz,
        torque_col       = torque_col,
        passive_mode     = passive_mode
      ))
  }

  acttorque <- act_rows |>
    dplyr::group_by(dplyr::across(dplyr::all_of(act_grp))) |>
    dplyr::summarise(avg_act_torque.Nm = dplyr::first(.data$.torque), .groups = "drop")

  combo <- dplyr::left_join(acttorque, passtorque,
                             by = intersect(pass_grp, names(acttorque))) |>
    dplyr::mutate(
      muscle_torque.Nm = .data$avg_act_torque.Nm - .data$avg_pass_torque.Nm,
      pass_torque.Nm   = .data$avg_pass_torque.Nm,
      t.interval.s     = 1.0 / sample_freq.Hz,
      torque_col       = torque_col,
      passive_mode     = passive_mode
    )

  combo
}


# =============================================================================
# add_muscle_instantaneous / summarize_muscle_cycles
# =============================================================================

#' Instantaneous muscle power per sample. Normalizes by muscle_mass.kg when present.
add_muscle_instantaneous <- function(muscletorque, mass_col = "muscle_mass.kg") {
  grp <- intersect(c("trial", "cycle", "filename"), names(muscletorque))
  out <- muscletorque |>
    dplyr::group_by(dplyr::across(dplyr::any_of(grp))) |>
    dplyr::mutate(
      insta_power.W = .data$muscle_torque.Nm *
        ((.data$pos.rad - dplyr::lag(.data$pos.rad)) / .data$t.interval.s)
    ) |>
    dplyr::ungroup()
  if (mass_col %in% names(out) && any(is.finite(out[[mass_col]]))) {
    out <- out |>
      dplyr::mutate(insta_power.Wkg = .data$insta_power.W / .data[[mass_col]])
  } else {
    out$insta_power.Wkg <- NA_real_
  }
  out
}

#' Per-cycle muscle work (trapz) and mean power.
summarize_muscle_cycles <- function(muscletorque, cycles_keep = NULL,
                                    mass_normalize = TRUE, mass_col = "muscle_mass.kg") {
  df <- muscletorque
  if (!is.null(cycles_keep)) df <- filter_cycles_if_set(df, cycles_keep)
  df <- add_muscle_instantaneous(df, mass_col = mass_col)

  df |>
    dplyr::filter(!is.na(.data$muscle_torque.Nm), !is.na(.data$pos.rad)) |>
    dplyr::group_by(dplyr::across(dplyr::any_of(c(
      "filename", "fishcode", "trial", "cycle", "duty", "phase",
      "freq.Hz", "curvature.invm", "stim"
    )))) |>
    dplyr::summarise(
      work.J         = calc_work(.data$pos.rad, .data$muscle_torque.Nm),
      avg_power.W    = mean(.data$insta_power.W, na.rm = TRUE),
      peak_torque.Nm = max(abs(.data$muscle_torque.Nm), na.rm = TRUE),
      muscle_mass.kg = dplyr::first(.data[[mass_col]]),
      .groups        = "drop"
    ) |>
    dplyr::mutate(
      work.Jkg    = if (mass_normalize && is.finite(.data$muscle_mass.kg[1L]))
        .data$work.J / .data$muscle_mass.kg else NA_real_,
      avg_power.Wkg = if (mass_normalize && is.finite(.data$muscle_mass.kg[1L]))
        .data$avg_power.W / .data$muscle_mass.kg else NA_real_
    )
}


# =============================================================================
# Coughlin power QC (ported from filter_coughlin_power.R)
# =============================================================================

coughlin_steady_state_power_limits <- function() {
  list(mean = 133.0, sd = 19.0, lo = 133.0 - 19.0, hi = 133.0 + 19.0,
       reference = "Coughlin et al. 1996, scup red muscle, steady-state")
}

#' Flag steps exceeding the Coughlin steady-state upper power limit (W/kg).
flag_coughlin_power_violations <- function(muscle_df, active_only = TRUE) {
  limits <- coughlin_steady_state_power_limits()
  df <- add_nominal_curvature_if_missing(muscle_df) |>
    add_muscle_instantaneous()
  if (isTRUE(active_only) && "cycletype" %in% names(df)) {
    df <- dplyr::filter(df, .data$cycletype == "act")
  }
  grp_keys <- intersect(c("filename", "trial", "freq.Hz", "curvature.invm"), names(df))
  df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grp_keys))) |>
    dplyr::summarise(
      peak_power_Wkg = max(abs(.data$insta_power.Wkg), na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(exceeds_coughlin = .data$peak_power_Wkg > limits$hi)
}


# =============================================================================
# filter_cycle_quality (sinusoid / torque-curvature QC)
# =============================================================================

sinusoid_quality_defaults <- function() {
  list(max_abs_r_decoupled = 0.25, min_r_harmonic_decoupled = 0.55,
       min_r_harmonic = 0.40, max_roughness = 0.12, max_abs_r_rough = 0.30,
       max_abs_r_weak = 0.20, min_r_harmonic_weak = 0.35, min_n = 20L)
}

#' Quality metrics for one trial x frequency x curvature step.
score_step_sinusoid_quality <- function(df, min_n = 20L,
                                        torque_col = "muscle_torque.Nm",
                                        curvature_col = "curve.invm") {
  if (nrow(df) < min_n || !torque_col %in% names(df)) return(NULL)
  tq <- df[[torque_col]]
  if (all(!is.finite(tq))) return(NULL)
  cv  <- if (curvature_col %in% names(df)) df[[curvature_col]] else NA_real_
  ang <- if ("angle.deg" %in% names(df)) df[["angle.deg"]] else NA_real_
  r_tc <- if (all(is.finite(cv)))  suppressWarnings(stats::cor(tq, cv,  use = "complete.obs")) else NA_real_
  r_ta <- if (all(is.finite(ang))) suppressWarnings(stats::cor(tq, ang, use = "complete.obs")) else NA_real_
  tibble::tibble(n = nrow(df), r_torque_curve = r_tc, r_torque_angle = r_ta)
}


# =============================================================================
# segmented_finite stubs
# =============================================================================

#' (STUB) Force-length curve analysis for isometric steps.
#'
#' Parses step_manifest + index_step_* from the loaded tibble and returns a
#' per-step summary. FL curve fitting and Hill-equation fit are deferred to a
#' follow-on session (M3 approved scope boundary).
#'
#' @param td Tibble from load_bender_flat() for a segmented_finite isometric file.
#' @return Per-step tibble with step_number, target_angle_degree, and raw torque
#'   statistics (mean passive, mean active). Fit columns are NA until implemented.
analyze_isometric <- function(td) {
  stopifnot(is.data.frame(td))
  if (!"step_number" %in% names(td)) {
    cli::cli_abort("analyze_isometric: step_number column absent; load with load_bender_flat()")
  }

  tq_col <- intersect(c("xtorque.Nm", "bodytorque.Nm"), names(td))[1L]
  if (is.na(tq_col)) {
    cli::cli_abort("analyze_isometric: no torque column found (need xtorque.Nm or bodytorque.Nm)")
  }

  td |>
    dplyr::group_by(.data$step_number) |>
    dplyr::summarise(
      n_samples       = dplyr::n(),
      target_angle_deg = dplyr::first(dplyr::coalesce(
        .data[["index_step_target_angle_degree"]], NA_real_
      )),
      torque_mean_Nm  = mean(.data[[tq_col]], na.rm = TRUE),
      torque_sd_Nm    = stats::sd(.data[[tq_col]], na.rm = TRUE),
      torque_peak_Nm  = max(abs(.data[[tq_col]]), na.rm = TRUE),
      # FL fit deferred -- TODO implement Hill fit
      FL_Fo_Nm        = NA_real_,
      FL_Lo_deg       = NA_real_,
      .groups = "drop"
    )
}


#' (STUB) Force-velocity curve analysis for isovelocity steps.
#'
#' Returns a per-step tibble with velocity and torque statistics.
#' P0 estimation and Hill FV fit are deferred to a follow-on session.
#'
#' @param td Tibble from load_bender_flat() for a segmented_finite isovelocity file.
#' @return Per-step tibble with step_number, velocity_deg_s, and torque statistics.
analyze_isovelocity <- function(td) {
  stopifnot(is.data.frame(td))
  if (!"step_number" %in% names(td)) {
    cli::cli_abort("analyze_isovelocity: step_number column absent; load with load_bender_flat()")
  }

  tq_col <- intersect(c("xtorque.Nm", "bodytorque.Nm"), names(td))[1L]
  if (is.na(tq_col)) {
    cli::cli_abort("analyze_isovelocity: no torque column found (need xtorque.Nm or bodytorque.Nm)")
  }

  td |>
    dplyr::group_by(.data$step_number) |>
    dplyr::summarise(
      n_samples        = dplyr::n(),
      velocity_deg_s   = dplyr::first(dplyr::coalesce(
        .data[["index_step_velocity_degree_per_second"]], NA_real_
      )),
      torque_mean_Nm   = mean(.data[[tq_col]], na.rm = TRUE),
      torque_sd_Nm     = stats::sd(.data[[tq_col]], na.rm = TRUE),
      torque_peak_Nm   = max(abs(.data[[tq_col]]), na.rm = TRUE),
      # FV fit deferred -- TODO implement Hill FV + P0/Vmax estimates
      FV_P0_Nm         = NA_real_,
      FV_Vmax_deg_s    = NA_real_,
      FV_power_max_W   = NA_real_,
      .groups = "drop"
    )
}
