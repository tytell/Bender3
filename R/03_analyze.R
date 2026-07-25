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
#   source("R/00_load_bender_flat.R")
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
#'
#' Cycle-level act/pass is determined from OBSERVED stim (does `stim != "0"`
#' occur ANYWHERE within this physical `cycle`, i.e. `floor(t.norm)`), never
#' from the `is_active_by_cycle` design flag, even when that column is
#' present. FIXED 2026-07-24 (PI-directed dynamic-power audit): `cycle`
#' (`00_load_bender_flat.R`'s `floor(t.norm)`) and `cycle_index` (the raw
#' HDF5 `/metadata/index_cycle_*` design-table index that `is_active_by_cycle`
#' is actually attached against) are TWO INDEPENDENT counters with no shared
#' time origin -- collapsing is_active_by_cycle to a per-`cycle` label via
#' `first()` silently mislabelled the active window by ~5 cycles (~1.7 s) on
#' EVERY active dynamic trial audited across bass16/17/18 (real stim
#' delivery ran ~1.7 s LATER and ~1.7 s LONGER than the design flag claimed,
#' every time -- not a one-off). Observed stim is authoritative regardless
#' of any such design/delivery timing offset. See
#' analysis_muscle_force_vector_log.md, "dynamic-cycle act/pass
#' cycle/cycle_index misalignment" for the audit.
set_cycle_types <- function(df) {
  df1 <- df |>
    dplyr::group_by(dplyr::across(dplyr::any_of(c("fishcode", "trial", "cycle")))) |>
    dplyr::mutate(cycletype = dplyr::if_else(any(.data$stim != "0"), "act", "pass")) |>
    dplyr::ungroup()

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
  # NOTE: actual bass16 corpus `metadata/protocol_type` values are the schema-locked
  # snake_case strings "dynamic"/"frequency_sweep" (verified against real files,
  # 2026-07-15) -- NOT the earlier-assumed "constant_frequency"/
  # "frequency_amplitude_combo" Title/snake sub-type strings, which never occur on
  # disk. "dynamic" covers both single-condition and multi-parameter-combo trials;
  # sub-layout is derived from the realized n_freq/n_amp/n_duty/n_phase counts.
  layout <- dplyr::case_when(
    identical(stimulus_type, "frequency_sweep") ~ "frequency_sweep",
    identical(stimulus_type, "dynamic") &
      n_freq <= 1L & n_duty <= 1L & n_phase <= 1L & (is.na(n_amp) | n_amp <= 1L) ~ "single_condition",
    identical(stimulus_type, "dynamic") & n_freq > 1L & !is.na(n_amp) & n_amp > 1L ~ "freq_amp_combo",
    identical(stimulus_type, "dynamic") ~ "constant_multi_param",
    # legacy/alternate protocol_type spellings, kept for forward compatibility
    identical(stimulus_type, "constant_frequency") & n_freq <= 1L & n_duty <= 1L & n_phase <= 1L ~ "single_condition",
    identical(stimulus_type, "constant_frequency") ~ "constant_multi_param",
    identical(stimulus_type, "frequency_amplitude_combo") ~ "freq_amp_combo",
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
#' @param phase_match_all_rows Only used with include_all_active_samples =
#'   TRUE and passive_mode = "halfcycle" (dynamic trials). Whole-CYCLE
#'   act/pass gating (from index_cycle_active) marks a cycle "pass" the
#'   moment the protocol's DESIGNED active block ends -- but real muscle
#'   force does not vanish at that instant; it decays over the following
#'   relaxation period, which typically falls inside the very next
#'   "pass"-labeled cycle(s). With the default FALSE, those rows are
#'   dropped entirely (not just set NA -- absent from the returned tibble),
#'   truncating any relaxation-tail display at the block boundary -- most
#'   visible at higher cycle frequencies (e.g. 3 Hz => 0.33 s/cycle) where
#'   the fixed post-stim relaxation window spans several cycles and hits a
#'   "pass" cycle almost immediately. Setting TRUE returns EVERY row
#'   (both "act"- and "pass"-labeled cycles) subtracted against the SAME
#'   phase(t.perc)-matched passive baseline used for "act" rows -- for
#'   true rest cycles this correctly yields a near-zero residual; for the
#'   relaxation tail immediately after a block, it reveals the decaying
#'   real force instead of an artificial cutoff. Does not change which
#'   rows count as "pass" when building the passive baseline itself.
calc_muscle_torque <- function(data, torque_col = NULL, sample_freq.Hz = NULL,
                               group_keys = NULL, cycles_keep = NULL,
                               passive_mode = "auto",
                               include_all_active_samples = FALSE,
                               phase_match_all_rows = FALSE) {
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
    act_rows  <- if (isTRUE(phase_match_all_rows)) df else dplyr::filter(df, .data$cycletype == "act")
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
#'
#' ONE ROW PER PHYSICAL CYCLE (FIXED 2026-07-24, PI-directed: "it is NOT
#' accurate to quantify cycle power/work for only the duration of the right
#' stimulus ... you have to take the average difference over the full range
#' of active cycles, which includes the left-right pairs working together").
#' Previously grouped by `cycle` AND `stim` together, which SPLIT every
#' physical bending cycle into up to 2 fragments (an L-window row and an
#' R-window row), each integrated over only that side's own
#' `[stim_onset, stim_offset + relaxation_s]` slice -- never the whole
#' stroke both muscles drive together. `stim` is no longer a grouping key;
#' work/power now integrate over every non-NA sample of the whole cycle
#' (both the L- and R-attributed portions, phase-matched act-minus-passive
#' torque from calc_muscle_torque()). `sides_present` replaces the old
#' single-valued `stim` output column so which side(s) actually contributed
#' to a given cycle stays visible for QA.
summarize_muscle_cycles <- function(muscletorque, cycles_keep = NULL,
                                    mass_normalize = TRUE, mass_col = "muscle_mass.kg") {
  df <- muscletorque
  if (!is.null(cycles_keep)) df <- filter_cycles_if_set(df, cycles_keep)
  df <- add_muscle_instantaneous(df, mass_col = mass_col)

  df |>
    dplyr::filter(!is.na(.data$muscle_torque.Nm), !is.na(.data$pos.rad)) |>
    dplyr::group_by(dplyr::across(dplyr::any_of(c(
      "filename", "fishcode", "trial", "cycle", "duty", "phase",
      "freq.Hz", "curvature.invm"
    )))) |>
    dplyr::summarise(
      work.J         = calc_work(.data$pos.rad, .data$muscle_torque.Nm),
      avg_power.W    = mean(.data$insta_power.W, na.rm = TRUE),
      peak_power.W   = max(abs(.data$insta_power.W), na.rm = TRUE),
      peak_torque.Nm = max(abs(.data$muscle_torque.Nm), na.rm = TRUE),
      muscle_mass.kg = if (mass_col %in% names(df)) dplyr::first(.data[[mass_col]]) else NA_real_,
      sides_present  = if ("stim" %in% names(df))
        paste(sort(unique(stats::na.omit(as.character(.data$stim)))), collapse = "+")
        else NA_character_,
      .groups        = "drop"
    ) |>
    dplyr::mutate(
      work.Jkg    = if (mass_normalize && is.finite(.data$muscle_mass.kg[1L]))
        .data$work.J / .data$muscle_mass.kg else NA_real_,
      avg_power.Wkg = if (mass_normalize && is.finite(.data$muscle_mass.kg[1L]))
        .data$avg_power.W / .data$muscle_mass.kg else NA_real_,
      # Peak (not mean) instantaneous mass-specific power, for direct
      # comparison against Coughlin et al. 1996 scup red-muscle steady-state
      # power limits (coughlin_steady_state_power_limits(), below) -- "is
      # this twitch weak or strong" per PI direction, 2026-07-16.
      peak_power.Wkg = if (mass_normalize && is.finite(.data$muscle_mass.kg[1L]))
        .data$peak_power.W / .data$muscle_mass.kg else NA_real_,
      coughlin_limit_lo_Wkg = coughlin_steady_state_power_limits()$lo,
      coughlin_limit_hi_Wkg = coughlin_steady_state_power_limits()$hi,
      exceeds_coughlin_hi = .data$peak_power.Wkg > .data$coughlin_limit_hi_Wkg,
      below_coughlin_lo   = .data$peak_power.Wkg < .data$coughlin_limit_lo_Wkg
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
# segmented_finite: isometric (FL) / isovelocity (FV) analysis
#
# COMPLETED 2026-07-15 (was a stub through the M1 dispatcher milestone;
# BUILD_PLAN.md M2+ "Segmented analysis" scope). Requires
# R/muscle_geometry.R and R/fit_fv_fl.R to be sourced first.
#
# Windowing (PI-confirmed, 2026-07-15):
#   Force_passive = mean torque over the step's pre-stim baseline window
#     [index_step_t_pre_baseline_start_second, _end_second].
#   Force_active  = mean torque over [index_step_stim_t0_second,
#     index_step_stim_t1_second + 0.5s] -- the 0.5 s tail implements the
#     "continue passive subtraction for 0.5s after stimulus removal"
#     deactivation-window rule (twitch decay outlasts the electrical pulse).
#   Force_muscle  = Force_active - Force_passive, torque_inertia_corrected_Nm.
#
# ADDITIONAL interpolated-baseline method (PI-directed, 2026-07-16, isometric
# only -- see build_segmented_step_summary()'s passive_force_Nm_interp
# comment): the static pre-stim-only baseline above is contaminated by
# ongoing viscoelastic stress relaxation that hasn't settled by the time the
# active window is sampled, biasing Force_muscle by an amount that scales
# with |displacement| and can mimic a monotonic FL trend. muscle_force_Nm_interp
# (analyze_isometric()) linearly interpolates the passive reference between
# the pre- and post-baseline window means instead. This is KEPT SEPARATE from
# Force_muscle above (not a replacement) and is only surfaced in
# "_baselineInterp"-suffixed plots (run_fv_fl_power_pipeline.R;
# see FIGURES_README.md for the naming convention).
#
# NOTE (PI direction, 2026-07-16): isometric FL summary plots no longer fit a
# model (parabola) to either muscle_force_Nm or muscle_force_Nm_interp by
# default -- see fit_fv_fl.R module header / build_summary_plot_isometric()
# (plot_summary_profiles.R). analyze_isometric() therefore returns EMPTY
# fits/fits_interp lists; they are kept as fields (rather than removed) only
# so the shared isometric/isovelocity per-trial reporting code in
# run_fv_fl_power_pipeline.R doesn't need a separate code path.
#
# Muscle side / contraction mode: see resolve_step_contraction()
# (muscle_geometry.R). Unilateral steps fold into ONE curve per stimulated
# side (concentric + eccentric together); bilateral_simultaneous steps are
# reported separately and EXCLUDED from the per-side FL/FV fits (net torque
# can't be attributed to one muscle when both are stimulated at once).
#
# FL x-axis / FV x-axis use predicted MUSCLE strain / strain-rate (not raw
# joint angle/velocity), from the commanded operating_point run through the
# same curvature -> strain geometry as the rest of the pipeline (see
# compute_predicted_strain(), muscle_geometry.R): shortening_strain_pct =
# shortening_value * pi/180 / dclamp_m * r_m * 100.
# =============================================================================

#' Read the per-step index_step_* geometry/timing table + file-level geometry
#' attrs needed for segmented FL/FV analysis.
.read_segmented_step_geometry <- function(filename) {
  h5 <- rhdf5::H5Fopen(filename, "H5F_ACC_RDONLY")
  on.exit(try(rhdf5::H5Fclose(h5), silent = TRUE), add = TRUE)

  m_attrs <- tryCatch(rhdf5::h5readAttributes(h5, "/metadata"), error = function(e) list())
  m_a  <- function(key) m_attrs[[key]]
  m_ds <- function(key) {
    path <- paste0("/metadata/", key)
    if (rhdf5::H5Lexists(h5, path)) tryCatch(rhdf5::h5read(h5, path), error = function(e) NULL) else NULL
  }
  dbl1 <- function(v, default = NA_real_) {
    v <- suppressWarnings(as.numeric(v[1L]))
    if (length(v) == 0L || is.na(v)) default else v
  }

  step_number <- as.integer(m_ds("index_step_number"))
  if (is.null(step_number) || length(step_number) == 0L) {
    cli::cli_abort(".read_segmented_step_geometry: index_step_number absent in {basename(filename)}")
  }

  # block_number/direction are only present when use_block_stim was active for
  # the trial (schema SS4 "Block metadata") -- default to NA so callers (e.g.
  # the fatigue-check plot) can fall back to treating the whole trial as one
  # block rather than erroring on older/non-block protocols.
  block_number <- suppressWarnings(as.numeric(m_ds("index_step_block_number")))
  if (is.null(block_number) || length(block_number) != length(step_number)) {
    block_number <- rep(NA_real_, length(step_number))
  }
  block_direction <- m_ds("index_step_block_direction")
  if (is.null(block_direction) || length(block_direction) != length(step_number)) {
    block_direction <- rep(NA_character_, length(step_number))
  }

  # wall_clock_start (ISO-8601 string, per-step) -- single point where this
  # gets surfaced into the pipeline; flows through build_segmented_step_summary()
  # -> analyze_isometric()/analyze_isovelocity() -> attach_vector_muscle_force()
  # unchanged (plain left_join column), reaching plot_fatigue_timeline.R without
  # further plumbing. Added 2026-07-21 for the near-L0-vs-real-elapsed-time
  # fatigue timeline (see analysis_muscle_force_vector_log.md, Gate A).
  wall_clock_start_raw <- m_ds("index_step_wall_clock_start")
  if (is.null(wall_clock_start_raw) || length(wall_clock_start_raw) != length(step_number)) {
    wall_clock_start_raw <- rep(NA_character_, length(step_number))
  }

  steps <- tibble::tibble(
    step_number             = step_number,
    operating_point          = as.numeric(m_ds("index_step_operating_point")),
    operating_point_units    = as.character(m_ds("index_step_operating_point_units")),
    recruitment              = as.character(m_ds("index_step_recruitment")),
    wall_clock_start         = as.character(wall_clock_start_raw),
    stim_t0_s                = as.numeric(m_ds("index_step_stim_t0_second")),
    stim_t1_s                = as.numeric(m_ds("index_step_stim_t1_second")),
    t_pre_baseline_start_s   = as.numeric(m_ds("index_step_t_pre_baseline_start_second")),
    t_pre_baseline_end_s     = as.numeric(m_ds("index_step_t_pre_baseline_end_second")),
    t_post_baseline_start_s  = as.numeric(m_ds("index_step_t_post_baseline_start_second")),
    t_post_baseline_end_s    = as.numeric(m_ds("index_step_t_post_baseline_end_second")),
    block_number             = as.integer(block_number),
    block_direction          = as.character(block_direction)
  )

  local_body_width_mm  <- dbl1(m_a("measurement_specimen_local_body_width_millimeter"))
  local_body_height_mm <- dbl1(m_a("measurement_specimen_local_body_height_millimeter"))
  dclamp_mm            <- dbl1(m_a("measurement_clamp_separation_millimeter"))
  density_g_per_mm3    <- dbl1(m_a("measurement_specimen_density_gram_per_cubic_millimeter"))

  # Mass/CSA estimate (PI-directed, 2026-07-16) for mass-/area-specific
  # properties -- see compute_muscle_mass_and_csa() (muscle_geometry.R) for
  # the oval-CSA x test-section-length x 3%-fraction recipe and caveats.
  muscle <- compute_muscle_mass_and_csa(local_body_width_mm, local_body_height_mm,
                                        dclamp_mm, density_g_per_mm3)

  list(
    steps               = steps,
    local_body_width_mm = local_body_width_mm,
    local_body_height_mm = local_body_height_mm,
    muscle_depth_mm_raw  = dbl1(m_a("measurement_target_muscle_depth_millimeter")),
    dclamp_mm            = dclamp_mm,
    density_g_per_mm3    = density_g_per_mm3,
    muscle               = muscle,
    lidx_pos_motor       = dbl1(m_a("daq_specimen_lateral_index_on_positive_motor_side")),
    lidx_left            = dbl1(m_a("daq_specimen_side_index_left")),
    lidx_right           = dbl1(m_a("daq_specimen_side_index_right"))
  )
}


# ACTIVE-window scalar-force extraction width (s) for "Method D" (PI-
# directed, 2026-07-22, following review of
# muscleforcemethodcompare.png/muscleforcemethodsensitivity.png): same value
# as muscle_force_vector.R's MFV_PEAK_WINDOW_S, defined separately here
# because this file does not always co-source muscle_force_vector.R (some
# callers use only the legacy zTorque-only path). Keep the two in sync if
# either changes.
LEGACY_PEAK_WINDOW_S <- 0.15

#' "Method D" extraction for a step's ACTIVE window: mean of RAW
#' (unsmoothed) `v` samples in a window `peak_window_s` wide, centered on
#' the SMOOTHED trace's own peak (via .smooth_trace_display_only(),
#' plot_force_vs_time.R -- co-sourced with this file everywhere it is
#' called), with the peak SEARCH restricted to `search_mask` (the true stim
#' duration, EXCLUDING the post-stim deactivation tail -- searching into the
#' tail let a settling/recovery transient there get mistaken for the peak,
#' found empirically on muscleforcemethodcompare.png's isometric panel)
#' while the narrow averaging window itself is bounded by `avg_mask` (the
#' full active+deactivation window). Falls back to the plain mean over
#' `avg_mask` whenever there are too few finite samples to smooth/search
#' reliably, OR whenever the search window itself is SHORTER than
#' `peak_window_s` -- see muscle_force_vector.R's .mfv_window_peak_means()
#' for why (found empirically 2026-07-22 on the ~0.05s dynamic L0 bookends:
#' a smoothed trace still rising at a too-short search window's own edge
#' pins the "peak" to that edge, and the narrow window then balloons past
#' the search boundary into the deactivation tail, landing on transients
#' unrelated to the pulse itself). Replaces a plain full-window MEAN
#' (PI-directed, 2026-07-22): validated on 92 real isometric steps
#' (muscleforcemethodsensitivity.png) to track the plain mean almost
#' exactly for FLAT/sustained holds, while correctly capturing genuine
#' rise-then-decay TRANSIENTS (isovelocity's moving ramps) that a plain
#' mean dilutes.
.legacy_peak_window_mean <- function(v, t, avg_mask, search_mask, peak_window_s = LEGACY_PEAK_WINDOW_S) {
  v_avg <- v[avg_mask]
  full_mean <- mean(v_avg, na.rm = TRUE)
  if (sum(is.finite(v_avg)) < 3L) return(full_mean)
  v_search <- v[search_mask]; t_search <- t[search_mask]
  if (sum(is.finite(v_search)) < 3L) return(full_mean)
  if (diff(range(t_search, na.rm = TRUE)) < peak_window_s) return(full_mean)
  v_smooth <- .smooth_trace_display_only(v_search)
  smax <- max(v_smooth, na.rm = TRUE); smin <- min(v_smooth, na.rm = TRUE)
  if (!is.finite(smax) || !is.finite(smin)) return(full_mean)
  use_max <- abs(smax - full_mean) >= abs(smin - full_mean)
  t_peak <- t_search[if (use_max) which.max(v_smooth) else which.min(v_smooth)]
  narrow <- avg_mask & t >= (t_peak - peak_window_s / 2) & t <= (t_peak + peak_window_s / 2)
  if (sum(is.finite(v[narrow])) < 1L) return(full_mean)
  mean(v[narrow], na.rm = TRUE)
}

#' Deactivation-window active-sample mask for a segmented (isometric /
#' isovelocity) td: TRUE within [stim_t0, stim_t1 + 0.5s] of the sample's step.
.segmented_active_mask <- function(td, step_geom, deactivation_window_s = 0.5) {
  active <- rep(FALSE, nrow(td))
  steps  <- step_geom$steps
  for (i in seq_len(nrow(steps))) {
    if (!is.finite(steps$stim_t0_s[i]) || !is.finite(steps$stim_t1_s[i])) next
    in_step <- td$step_number == steps$step_number[i]
    in_win  <- td$t.s >= steps$stim_t0_s[i] & td$t.s <= (steps$stim_t1_s[i] + deactivation_window_s)
    active[in_step & in_win] <- TRUE
  }
  active
}


#' Velocity-matched passive torque for isovelocity steps (PI-directed fix,
#' 2026-07-15): a STATIC pre-motion baseline can't capture the
#' velocity-dependent viscoelastic drag present once the specimen is
#' actually moving, which biases active-minus-passive "muscle force" at
#' higher commanded speeds (concentric trends spuriously negative, eccentric
#' spuriously inflated positive). Instead, match each STIMULATED step to the
#' no-stim step in the same trial commanded to the SAME operating_point
#' (same speed AND bend direction) and use that step's own torque over the
#' analogous window as the passive reference.
#'
#' No-stim steps are detected from the ACTUAL delivered stim (td$stim != "0"
#' during the step's samples), not from index_step_recruitment -- this
#' corpus's bilateral_simultaneous "calibration" steps list a recruited side
#' but never actually fire the stimulator (stim_side stays "none"
#' throughout), so they are exactly the velocity-matched no-stim ramps this
#' fix needs, despite their recruitment label.
#'
#' @return Named numeric vector (by step_number) of matched passive torque;
#'   NA for steps with no same-operating_point no-stim match (caller should
#'   fall back to the static baseline in that case rather than dropping the
#'   step).
.isovelocity_velocity_matched_passive <- function(td, steps, deactivation_window_s) {
  is_no_stim_step <- vapply(steps$step_number, function(sn) {
    rows <- td$step_number == sn
    !any(as.character(td$stim[rows]) != "0", na.rm = TRUE)
  }, logical(1L))

  matched <- vapply(seq_len(nrow(steps)), function(i) {
    if (is_no_stim_step[i]) return(NA_real_)
    match_idx <- which(is_no_stim_step & abs(steps$operating_point - steps$operating_point[i]) < 1e-3)
    if (length(match_idx) == 0L) return(NA_real_)
    match_step <- steps$step_number[match_idx[1L]]
    rows <- td$step_number == match_step
    win <- td$t.s[rows] >= steps$stim_t0_s[i] & td$t.s[rows] <= (steps$stim_t1_s[i] + deactivation_window_s)
    if (!any(win, na.rm = TRUE)) return(NA_real_)
    mean(td$torque_inertia_corrected_Nm[rows][win], na.rm = TRUE)
  }, numeric(1L))

  stats::setNames(matched, steps$step_number)
}


#' Shared per-step force/strain/side summary for isometric + isovelocity.
#' Returns list(step_summary, td) -- td gains strain_active_pct/
#' strain_passive_pct/torque_inertia_corrected_Nm columns for the required
#' per-trial compound plot; step_summary is one row per step with muscle
#' force, side, contraction mode, and the FL/FV x-axis value.
#'
#' @param passive_mode "static_baseline" (default; mean torque over the
#'   step's own pre-stim baseline window -- correct when the specimen isn't
#'   moving, i.e. isometric) or "velocity_matched" (isovelocity: passive
#'   torque taken from a same-operating_point no-stim step instead; falls
#'   back to the static baseline for any step lacking a no-stim match).
build_segmented_step_summary <- function(td, filename = attr(td, "Filename"),
                                         torque_col = "torque_inertia_corrected_Nm",
                                         deactivation_window_s = 0.5,
                                         passive_mode = c("static_baseline", "velocity_matched")) {
  passive_mode <- match.arg(passive_mode)
  stopifnot(is.data.frame(td))
  if (!"step_number" %in% names(td)) {
    cli::cli_abort("build_segmented_step_summary: step_number column absent; load with load_bender_flat()")
  }
  if (!torque_col %in% names(td)) {
    cli::cli_abort("build_segmented_step_summary: {torque_col} column absent -- run deconvolve_bender() and attach it first")
  }
  geom <- .read_segmented_step_geometry(filename)

  active_mask <- .segmented_active_mask(td, geom, deactivation_window_s)
  td <- attach_predicted_strain(td, geom$local_body_width_mm, geom$muscle_depth_mm_raw,
                                active_mask = active_mask)
  td$torque_inertia_corrected_Nm <- td[[torque_col]]

  depth <- resolve_muscle_depth_mm(geom$muscle_depth_mm_raw)
  dclamp_m <- geom$dclamp_mm / 1000.0

  side_tbl <- resolve_step_contraction(
    geom$steps$recruitment, geom$steps$operating_point,
    geom$lidx_pos_motor, geom$lidx_left, geom$lidx_right
  )
  steps <- dplyr::bind_cols(geom$steps, side_tbl)

  # commanded-operating-point -> muscle-strain(-rate) x-axis, using the same
  # curvature geometry as the rest of the pipeline (r_m from strain formula).
  # d(curvature)/dt = d(angle)/dt * pi/180 / dclamp_m, so running this same
  # formula on operating_point gives STRAIN (%) for isometric (operating_point
  # is an angle) and STRAIN RATE (%/s) for isovelocity (operating_point is an
  # angular velocity) -- shared column name, units follow the caller's protocol.
  strain_res <- compute_predicted_strain(1.0, geom$local_body_width_mm, depth$depth_mm) # for r_m only
  r_m <- strain_res$r_m
  steps$shortening_strain_pct <- steps$shortening_value * pi / 180.0 / dclamp_m * r_m * 100.0

  force_by_step <- td |>
    dplyr::group_by(.data$step_number) |>
    dplyr::summarise(
      passive_force_Nm_static = {
        s <- geom$steps[geom$steps$step_number == dplyr::first(.data$step_number), ]
        win <- .data$t.s >= s$t_pre_baseline_start_s & .data$t.s <= s$t_pre_baseline_end_s
        mean(.data$torque_inertia_corrected_Nm[win], na.rm = TRUE)
      },
      # Post-stim baseline mean (2026-07-16 addition): NOT used by the
      # original static-baseline muscle_force_Nm below (kept unchanged), but
      # feeds passive_force_Nm_interp just below -- see that block's comment
      # for why the pre-stim-only baseline is biased.
      post_force_Nm_static = {
        s <- geom$steps[geom$steps$step_number == dplyr::first(.data$step_number), ]
        win <- .data$t.s >= s$t_post_baseline_start_s & .data$t.s <= s$t_post_baseline_end_s
        mean(.data$torque_inertia_corrected_Nm[win], na.rm = TRUE)
      },
      active_force_Nm = {
        s <- geom$steps[geom$steps$step_number == dplyr::first(.data$step_number), ]
        win        <- .data$t.s >= s$stim_t0_s & .data$t.s <= (s$stim_t1_s + deactivation_window_s)
        search_win <- .data$t.s >= s$stim_t0_s & .data$t.s <= s$stim_t1_s
        if (any(win, na.rm = TRUE) && any(.data$is_active_sample[win], na.rm = TRUE)) {
          .legacy_peak_window_mean(.data$torque_inertia_corrected_Nm, .data$t.s, win, search_win)
        } else {
          NA_real_
        }
      },
      n_samples = dplyr::n(),
      .groups = "drop"
    )

  steps <- dplyr::left_join(steps, force_by_step, by = "step_number")

  # Interpolated-baseline alternative (2026-07-16, PI-directed investigation
  # of the bass17/bass16 isometric FL curve resembling a monotonic "total
  # tension" shape instead of a bell): binning raw torque across a full
  # isometric hold shows it keeps decaying continuously from the pre-stim
  # baseline window all the way through the post-stim baseline window
  # (viscoelastic stress relaxation of the bent specimen, amplitude scaling
  # with |displacement|) -- there is no settled plateau. Because
  # passive_force_Nm_static is sampled EARLY in that decay and
  # active_force_Nm is sampled slightly LATER along the SAME decay,
  # active - passive_force_Nm_static is biased by an amount that grows with
  # |displacement| regardless of stimulation, which masquerades as a
  # monotonic force-vs-strain trend. This interpolates the passive
  # reference (evaluated at the active window's own mean time) linearly
  # between the pre- and post-baseline window means, cancelling a
  # (locally) linear relaxation trend instead of assuming the passive
  # torque is flat across the whole hold. Falls back to the static value
  # when no valid post-baseline mean exists (e.g. trailing step cut off).
  # Kept STRICTLY ADDITIONAL: passive_force_Nm/muscle_force_Nm (below) are
  # untouched, so all existing plots/fits reproduce exactly as before --
  # this only adds *_interp columns, surfaced in separate
  # "_baselineInterp"-suffixed outputs.
  steps <- steps |>
    dplyr::mutate(
      .t_pre_mid_s    = (.data$t_pre_baseline_start_s + .data$t_pre_baseline_end_s) / 2,
      .t_post_mid_s   = (.data$t_post_baseline_start_s + .data$t_post_baseline_end_s) / 2,
      .t_active_mid_s = (.data$stim_t0_s + (.data$stim_t1_s + deactivation_window_s)) / 2,
      passive_force_Nm_interp = dplyr::if_else(
        is.finite(.data$post_force_Nm_static) & is.finite(.data$.t_post_mid_s) &
          (.data$.t_post_mid_s != .data$.t_pre_mid_s),
        .data$passive_force_Nm_static +
          (.data$post_force_Nm_static - .data$passive_force_Nm_static) *
            (.data$.t_active_mid_s - .data$.t_pre_mid_s) / (.data$.t_post_mid_s - .data$.t_pre_mid_s),
        .data$passive_force_Nm_static
      )
    ) |>
    dplyr::select(-".t_pre_mid_s", -".t_post_mid_s", -".t_active_mid_s")

  if (passive_mode == "velocity_matched") {
    matched <- .isovelocity_velocity_matched_passive(td, steps, deactivation_window_s)
    steps$passive_force_Nm_matched <- unname(matched[as.character(steps$step_number)])
    # fall back to the static baseline for any step lacking a same-velocity
    # no-stim match, rather than silently dropping it from the FL/FV data.
    steps$passive_force_Nm <- dplyr::if_else(
      is.finite(steps$passive_force_Nm_matched),
      steps$passive_force_Nm_matched, steps$passive_force_Nm_static
    )
    steps$passive_force_source <- dplyr::if_else(
      is.finite(steps$passive_force_Nm_matched), "velocity_matched", "static_baseline_fallback"
    )
  } else {
    steps$passive_force_Nm <- steps$passive_force_Nm_static
    steps$passive_force_source <- "static_baseline"
  }

  steps <- steps |>
    dplyr::mutate(
      # force_sign (muscle_geometry.R::resolve_step_contraction) re-expresses
      # RAW lab-frame torque (which follows commanded-angle sign -- opposite
      # for left- vs right-bend regardless of which muscle is active) into
      # the recruited muscle's OWN contraction-positive frame. Without this,
      # left/right FL & FV curves are mirror images of each other purely from
      # the lab-frame sign flip, not from real physiological sidedness.
      muscle_force_Nm       = .data$force_sign * (.data$active_force_Nm - .data$passive_force_Nm),
      # Interpolated-baseline counterpart of muscle_force_Nm above (see
      # passive_force_Nm_interp comment) -- always computed from the static
      # pre/post windows regardless of passive_mode, so it is available for
      # isovelocity steps too even though only the isometric path currently
      # surfaces it in plots/fits.
      muscle_force_Nm_interp = .data$force_sign * (.data$active_force_Nm - .data$passive_force_Nm_interp),
      muscle_depth_mm_used  = depth$depth_mm,
      muscle_depth_assumed  = depth$assumed
    )

  list(step_summary = steps, td = td, r_m = r_m, dclamp_m = dclamp_m, muscle = geom$muscle,
       muscle_depth_mm_used = depth$depth_mm, muscle_depth_assumed = depth$assumed)
}


#' Force-Length curve analysis for isometric steps (completes the M1 stub).
#'
#' Windows Force_active/Force_passive per step (deactivation-window rule) and
#' folds unilateral steps by stimulated muscle side (concentric + eccentric
#' together); bilateral_simultaneous steps are reported but excluded from the
#' per-side grouping (see module header).
#'
#' Does NOT fit a model (parabola) to the resulting FL points by default (PI
#' direction, 2026-07-16 -- see fit_fv_fl.R module header): the default FL
#' visualization is a connect-the-mean line with no imposed functional form
#' (build_summary_plot_isometric(), plot_summary_profiles.R). `fits` /
#' `fits_interp` below are therefore always EMPTY lists -- kept as fields
#' only so callers (run_fv_fl_power_pipeline.R) share one code path with
#' analyze_isovelocity(), which still returns real Hill fits.
#'
#' @param td Tibble from load_bender_flat() with torque_inertia_corrected_Nm
#'   already attached (see deconvolve_bender()).
#' @param filename Raw H5 path (defaults to attr(td, "Filename")).
#' @return list(step_summary, td, fits = list(), fits_interp = list()).
analyze_isometric <- function(td, filename = attr(td, "Filename")) {
  built <- build_segmented_step_summary(td, filename)
  steps <- built$step_summary

  list(step_summary = steps, td = built$td, fits = list(), fits_interp = list(),
       r_m = built$r_m, dclamp_m = built$dclamp_m, muscle = built$muscle,
       muscle_depth_mm_used = built$muscle_depth_mm_used,
       muscle_depth_assumed = built$muscle_depth_assumed)
}


#' Force-Velocity curve analysis for isovelocity steps (completes the M1 stub).
#'
#' Same side-folding as analyze_isometric(), but Force_passive uses a
#' VELOCITY-MATCHED no-stim reference step (same commanded operating_point)
#' rather than a static pre-motion baseline -- a static baseline can't
#' capture velocity-dependent viscoelastic drag once the specimen is moving,
#' which otherwise biases muscle_force_Nm at higher commanded speeds
#' (PI-directed fix, 2026-07-15; see build_segmented_step_summary()'s
#' passive_mode doc and .isovelocity_velocity_matched_passive()). Fits a
#' Hill hyperbola per side on the CONCENTRIC (shortening) limb, anchored at
#' F0 = the step's own zero-velocity (isometric) muscle force when available.
#'
#' @param td Tibble from load_bender_flat() with torque_inertia_corrected_Nm attached.
#' @param filename Raw H5 path (defaults to attr(td, "Filename")).
#' @return list(step_summary, td, fits = list(left=, right=)).
analyze_isovelocity <- function(td, filename = attr(td, "Filename")) {
  built <- build_segmented_step_summary(td, filename, passive_mode = "velocity_matched")
  steps <- built$step_summary

  fits <- list()
  for (sd in c("left", "right")) {
    ss <- dplyr::filter(steps, .data$muscle_side == sd, is.finite(.data$muscle_force_Nm))
    f0_row <- dplyr::filter(ss, abs(.data$shortening_value) < 1e-6)
    f0_iso <- if (nrow(f0_row) > 0L) mean(f0_row$muscle_force_Nm, na.rm = TRUE) else NA_real_
    ss_conc <- dplyr::filter(ss, .data$contraction_mode %in% c("concentric", "isometric_zero"))
    fits[[sd]] <- fit_force_velocity_curve(
      ss_conc$shortening_strain_pct, ss_conc$muscle_force_Nm,
      side_label = sd, f0_isometric = f0_iso
    )
    # Mass-/area-specific properties (PI-directed, 2026-07-16) -- specific
    # tension (N/cm^2, from F0) AND mass-specific peak power (W/kg, from the
    # Hill fit's composite peak_power) -- see add_specific_properties_to_fit()
    # (muscle_geometry.R) for the torque/strain-rate -> physical-units
    # derivation. No-op if the fit failed or geometry is unavailable.
    fits[[sd]] <- add_specific_properties_to_fit(fits[[sd]], built$r_m, built$dclamp_m, built$muscle, kind = "FV")
  }

  list(step_summary = steps, td = built$td, fits = fits,
       r_m = built$r_m, dclamp_m = built$dclamp_m, muscle = built$muscle,
       muscle_depth_mm_used = built$muscle_depth_mm_used,
       muscle_depth_assumed = built$muscle_depth_assumed)
}
