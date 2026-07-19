# muscle_force_vector.R
# Muscle force along its line of action (u_hat), recovered from the FULL
# 6-axis ATI Mini40 wrench instead of a single bending-torque axis.
#
# WHY (PI-directed, 2026-07-18): the existing FL/FV force path uses only the
# lateral bending moment (zTorque, via torque_inertia_corrected_Nm), which is
# dominated by PASSIVE body bending stiffness -- the active muscle signal is
# buried. The muscle runs longitudinally (sensor X) at a dorsoventral offset
# from the sensor origin, so its tension shows up LARGELY as a moment about
# sensor Y (yTorque) and secondarily across the other force/torque channels.
# Recovering muscle force therefore needs the whole wrench, not one axis.
#
# COORDINATE FRAME (ATI Mini40 raw sensor frame; NO sensor->lab transform in
# code -- see 01_calibrate.R, muscle_force_vector_physics.md):
#   X = longitudinal (+ anterior), Y = mediolateral (+ right),
#   Z = dorsoventral (+ ventral). Bending is lateral (XY plane); bend axis = Z.
#
# MOMENT ARMS (see muscle_force_vector_physics.md):
#   r_m = local_body_width/2 - muscle_depth   (mediolateral arm; Fx -> Tz)
#   d   = clamp_offset_vertical               (dorsoventral arm;  Fx -> Ty)
# The muscle attaches at offset r = (0, r_m, d) from the sensor origin, so a
# tension F along u_hat produces the wrench [F*u_hat ; r x (F*u_hat)].
#
# METHOD (matches PI A2/A5, 2026-07-18):
#   u_hat  : EMPIRICAL primary -- unit direction of the active-minus-passive
#            FORCE vector (dFx,dFy,dFz). GEOMETRIC cross-check -- beam tangent
#            from the commanded bend angle. Divergence is diagnostic.
#   Force  : projected-moment / moment-arm. The effective moment arm for a
#            unit force along u_hat at attachment r is |r x u_hat| (the true
#            perpendicular distance from the sensor origin to the line of
#            action); the projected active moment is dTau . (r x u_hat)/|..|.
#            This FUSES the yTorque (arm d) and zTorque (arm r_m) estimates
#            into one number: for u_hat ~ X it reduces to
#            (dTy*d - dTz*r_m)/(d^2 + r_m^2) == Fx when the model holds.
#   A pure force-channel estimate (dF . u_hat) is reported alongside as an
#   independent cross-check.
#
# All active-minus-passive: isometric uses the step's own pre-stim baseline
# window; isovelocity uses an ANGLE-matched no-stim ramp (match by enc.deg,
# not time) so the 5-deg-vs-0-deg start offset washes out. When no same-trial
# angle-overlapping passive exists, the ramp is borrowed CROSS-TRIAL but ONLY
# from the SAME individual (same pipeline run = one specimen) at the same
# signed velocity (compute_isovelocity_vector_batch()).
#
# SIGN CONVENTION (PI-directed, 2026-07-18; REVISED same day per PI
# biomechanical correction): left and right muscle forces are standardized
# to ONE convention so the two sides are directly comparable -- see the
# standardization block in .mfv_finalize_step(). The reference is the RAW
# PASSIVE torque's own direction at the matched angle/window (NOT an
# arbitrary "always positive" rule, superseded): a step whose active
# deviation REINFORCES the passive direction (concentric-consistent) comes
# out positive and more extreme in that direction; a step whose deviation
# OPPOSES the passive direction (eccentric-consistent) comes out negative or
# small. This preserves the physical relative-to-passive relationship instead
# of erasing it. See `pass_moment_N` / `tension_relative_to_passive` in the
# step summary for the per-step reference value and classification.
#
# FV SONO CONFIRMATION (PI-directed, 2026-07-18): the DS3 sono length (RIGHT
# muscle only -- the sole wired channel) is low-pass filtered at 40 Hz and used
# to sample each isovelocity ramp's force AT the L0 (resting-length) crossing,
# so FV points are taken at a common length rather than as a whole-ramp mean.
#
# NOT inertia-corrected: the z-axis inertial deconvolution (02_deconvolve.R)
# applies only to the bending axis. Isometric holds are motionless (alpha ~ 0)
# so the raw calibrated, bias-subtracted 6-axis wrench needs no inertia term;
# isovelocity's angle-matched passive subtraction removes motion-correlated
# (viscoelastic + inertial) artifacts componentwise. This is FLAGGED, not
# silently assumed.

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(rhdf5); library(cli)
})

# Minimum activation SNR (||dF_force|| / pre-stim force noise) for a step's
# EMPIRICAL u_hat to be trusted; below this the geometric u_hat is used for the
# direction instead and the step is flagged. Visible + vetoable per PI A6.
MFV_UHAT_SNR_MIN <- 3.0

# Deactivation tail (s) appended to the active window, matching the rest of
# the pipeline's DEACTIVATION_WINDOW_S.
MFV_DEACTIVATION_WINDOW_S <- 0.5

# Sono (DS3) low-pass cutoff (Hz) for de-staircasing the segment-length signal
# before L0-crossing detection. 40 Hz sits below the DS3 internal update rate
# (~247 Hz) yet preserves muscle-length dynamics -- the value validated as best
# in diag_sono_smoothing.R (PI-directed, 2026-07-18).
MFV_SONO_LP_HZ <- 40.0

# Signed-velocity tolerance (deg/s) for matching an active isovelocity step to a
# no-stim ramp of the same commanded velocity (within- or cross-trial).
MFV_VELOCITY_MATCH_TOL <- 1.0


# =============================================================================
# Geometry
# =============================================================================

#' Read the file-level geometry needed for the vector-force moment arms.
#' @return list(width_mm, depth_mm_raw, dvert_mm, lidx_pos_motor, lidx_left,
#'   lidx_right).
.mfv_read_geom <- function(filename) {
  h5 <- rhdf5::H5Fopen(filename, "H5F_ACC_RDONLY")
  on.exit(try(rhdf5::H5Fclose(h5), silent = TRUE), add = TRUE)
  m_attrs <- tryCatch(rhdf5::h5readAttributes(h5, "/metadata"), error = function(e) list())
  dbl1 <- function(key, default = NA_real_) {
    v <- suppressWarnings(as.numeric(m_attrs[[key]][1L]))
    if (length(v) == 0L || is.na(v)) default else v
  }
  list(
    width_mm       = dbl1("measurement_specimen_local_body_width_millimeter"),
    depth_mm_raw   = dbl1("measurement_target_muscle_depth_millimeter"),
    dvert_mm       = dbl1("measurement_clamp_offset_vertical_millimeter"),
    lidx_pos_motor = dbl1("daq_specimen_lateral_index_on_positive_motor_side"),
    lidx_left      = dbl1("daq_specimen_side_index_left"),
    lidx_right     = dbl1("daq_specimen_side_index_right")
  )
}

#' Muscle moment arms (meters) from file geometry.
#' r_m  = half body width - muscle depth (mediolateral, muscle-specific arm;
#'        resolve_muscle_depth_mm() supplies the assumed depth when unset).
#' r_body = half body width (mediolateral, passive-body arm; reported only).
#' d    = clamp_offset_vertical (dorsoventral arm).
#' @return list(r_m_m, r_body_m, d_m, depth_mm, depth_assumed).
resolve_muscle_moment_arms <- function(width_mm, depth_mm_raw, dvert_mm) {
  depth <- resolve_muscle_depth_mm(depth_mm_raw)  # muscle_geometry.R
  half_width_mm <- width_mm / 2.0
  r_m_m    <- (half_width_mm - depth$depth_mm) / 1000.0
  r_body_m <- half_width_mm / 1000.0
  d_m      <- dvert_mm / 1000.0
  list(r_m_m = r_m_m, r_body_m = r_body_m, d_m = d_m,
       depth_mm = depth$depth_mm, depth_assumed = depth$assumed)
}

.mfv_cross <- function(a, b) {
  c(a[2L] * b[3L] - a[3L] * b[2L],
    a[3L] * b[1L] - a[1L] * b[3L],
    a[1L] * b[2L] - a[2L] * b[1L])
}


# =============================================================================
# Channel access (prefer filtered <base>.Nm/.N over raw <base>0.Nm/0.N)
# =============================================================================

.MFV_FORCE_BASES  <- c(x = "xforce",  y = "yforce",  z = "zforce")
.MFV_TORQUE_BASES <- c(x = "xtorque", y = "ytorque", z = "ztorque")

.mfv_col <- function(td, base) {
  filt <- if (grepl("force", base)) paste0(base, ".N") else paste0(base, ".Nm")
  raw  <- if (grepl("force", base)) paste0(base, "0.N") else paste0(base, "0.Nm")
  if (filt %in% names(td)) td[[filt]] else if (raw %in% names(td)) td[[raw]] else NULL
}

#' Are all six calibrated F/T channels present on td?
.mfv_has_six_axis <- function(td) {
  all(vapply(c(.MFV_FORCE_BASES, .MFV_TORQUE_BASES),
             function(b) !is.null(.mfv_col(td, b)), logical(1L)))
}

#' Load a trial's full 6-axis calibrated (bias-subtracted, LP-filtered) wrench,
#' row-aligned to an already-loaded reference td (same file, same sample order).
#' Kept SEPARATE from the pipeline's primary load_bender_flat(..., loadtorques
#' = "x") call so the existing single-axis force path is byte-for-byte
#' unaffected -- this is a second, isolated read used only by the vector path.
#' @return td6 truncated to nrow(ref_td); NULL if load/alignment fails.
mfv_load_six_axis <- function(filename, ref_td) {
  td6 <- tryCatch(
    load_bender_flat(filename, do_filter = TRUE,
                     loadtorques = c("x", "y", "z"), loadforces = TRUE),
    error = function(e) { cli::cli_warn("mfv_load_six_axis: {conditionMessage(e)}"); NULL }
  )
  if (is.null(td6) || nrow(td6) < nrow(ref_td) || !.mfv_has_six_axis(td6)) return(NULL)
  td6 <- td6[seq_len(nrow(ref_td)), , drop = FALSE]
  # sanity: time axis must line up with the reference td
  if ("t.s" %in% names(td6) && "t.s" %in% names(ref_td)) {
    dmax <- suppressWarnings(max(abs(td6$t.s - ref_td$t.s), na.rm = TRUE))
    if (!is.finite(dmax) || dmax > 1e-6) {
      cli::cli_warn("mfv_load_six_axis: time axis misaligned (max dt {signif(dmax,3)}s) -- vector path skipped")
      return(NULL)
    }
  }
  td6
}


# =============================================================================
# u_hat + force from a wrench delta
# =============================================================================

#' Empirical line of action from an active-minus-passive FORCE vector.
#' @return list(uhat (unit 3-vec), angle_xy_deg (atan2(dFy,dFx)), mag (||dF||)).
uhat_empirical <- function(dF) {
  mag <- sqrt(sum(dF^2))
  if (!is.finite(mag) || mag <= 0) {
    return(list(uhat = c(NA_real_, NA_real_, NA_real_), angle_xy_deg = NA_real_, mag = NA_real_))
  }
  uhat <- dF / mag
  list(uhat = uhat, angle_xy_deg = atan2(dF[2L], dF[1L]) * 180 / pi, mag = mag)
}

#' Geometric (beam-tangent) line-of-action cross-check from the commanded bend
#' angle. FIRST-ORDER model (FLAGGED as approximate): the test section bends
#' through the joint angle theta across the clamp span, and the muscle's mean
#' line of action rotates in the XY plane by ~theta/2 relative to the fixed
#' longitudinal axis (half the end-to-end tangent rotation of a uniformly-bent
#' beam). Purely in-plane (uz = 0) -- any real dorsoventral fiber tilt shows up
#' only in the empirical u_hat, which is exactly the cross-check's purpose.
#' @param bend_angle_deg Signed commanded joint angle at the step (deg).
#' @return list(uhat, angle_xy_deg).
uhat_geometric <- function(bend_angle_deg) {
  if (!is.finite(bend_angle_deg)) {
    return(list(uhat = c(NA_real_, NA_real_, NA_real_), angle_xy_deg = NA_real_))
  }
  phi <- (bend_angle_deg / 2.0) * pi / 180.0
  list(uhat = c(cos(phi), sin(phi), 0.0), angle_xy_deg = phi * 180 / pi)
}

#' Muscle force from a wrench delta, given a line of action and moment arms.
#'
#' Primary (F_moment_N): projected-moment / moment-arm. Effective arm for a
#'   unit force along u_hat at attachment r = (0, r_m, d) is a = |r x u_hat|;
#'   projected active moment is dTau . (r x u_hat)/a; force = that / a.
#' Cross-checks: F_force_N = dF . u_hat (force channels only);
#'   F_yT_N = dTy / d and F_zT_N = -dTz / r_m (single-arm longitudinal-force
#'   estimates, valid when u_hat ~ X).
#'
#' Sign: F_moment_N/F_force_N are in the RAW lab frame (they follow bend
#' direction). Callers fold sides by multiplying with force_sign
#' (resolve_step_contraction(), muscle_geometry.R), exactly as the existing
#' muscle_force_Nm path does.
#'
#' @param dF,dTau length-3 active-minus-passive force (N) / torque (N*m) deltas.
#' @return list(F_moment_N, F_force_N, F_yT_N, F_zT_N, eff_arm_m, rxu). `rxu`
#'   (the r x u_hat axis used for the moment projection) is returned so callers
#'   can project OTHER torque vectors (e.g. the raw passive torque, for the
#'   passive-relative sign convention in .mfv_finalize_step()) onto the exact
#'   same axis/arm without recomputing it.
solve_muscle_force_from_wrench <- function(dF, dTau, uhat, r_m_m, d_m) {
  r <- c(0.0, r_m_m, d_m)
  rxu <- .mfv_cross(r, uhat)
  eff_arm <- sqrt(sum(rxu^2))
  F_moment <- if (is.finite(eff_arm) && eff_arm > 0) sum(dTau * rxu) / (eff_arm^2) else NA_real_
  F_force  <- if (all(is.finite(uhat))) sum(dF * uhat) else NA_real_
  F_yT <- if (is.finite(d_m)   && d_m   > 0) dTau[2L] / d_m    else NA_real_
  F_zT <- if (is.finite(r_m_m) && r_m_m > 0) -dTau[3L] / r_m_m else NA_real_
  list(F_moment_N = F_moment, F_force_N = F_force, F_yT_N = F_yT, F_zT_N = F_zT,
       eff_arm_m = eff_arm, rxu = rxu)
}


# =============================================================================
# Per-step vector force (isometric + isovelocity)
# =============================================================================

#' Mean of each of the 6 channels over a time window on one step.
.mfv_window_means <- function(td6, step_rows, win_mask) {
  sel <- step_rows & win_mask
  get_mean <- function(base) {
    v <- .mfv_col(td6, base)
    if (is.null(v)) NA_real_ else mean(v[sel], na.rm = TRUE)
  }
  list(
    F = c(get_mean("xforce"),  get_mean("yforce"),  get_mean("zforce")),
    T = c(get_mean("xtorque"), get_mean("ytorque"), get_mean("ztorque"))
  )
}

#' Pre-stim force-noise SD (N) over a step's baseline window -- denominator for
#' the empirical-u_hat activation SNR.
.mfv_baseline_force_noise <- function(td6, step_rows, base_mask) {
  sel <- step_rows & base_mask
  if (sum(sel, na.rm = TRUE) < 3L) return(NA_real_)
  comp_sd <- vapply(c("xforce", "yforce", "zforce"), function(b) {
    v <- .mfv_col(td6, b); if (is.null(v)) NA_real_ else stats::sd(v[sel], na.rm = TRUE)
  }, numeric(1L))
  sqrt(sum(comp_sd^2, na.rm = TRUE))
}

#' Identify no-stim steps (delivered stim all "0") in a td6.
.mfv_no_stim_steps <- function(td6, step_numbers) {
  is_ns <- vapply(step_numbers, function(sn) {
    rows <- td6$step_number == sn
    !any(as.character(td6$stim[rows]) != "0", na.rm = TRUE)
  }, logical(1L))
  step_numbers[is_ns]
}

#' Extract one step's channel-vs-angle "ramp" (enc.deg + the 6 channels over
#' the step's rows) -- the reusable unit for angle-matched passive subtraction,
#' both within-trial and cross-trial.
.mfv_ramp_from_step <- function(td6, rows) {
  gc6 <- function(b) { v <- .mfv_col(td6, b); if (is.null(v)) rep(NA_real_, sum(rows)) else v[rows] }
  list(ang = td6$enc.deg[rows],
       F = list(gc6("xforce"),  gc6("yforce"),  gc6("zforce")),
       T = list(gc6("xtorque"), gc6("ytorque"), gc6("ztorque")))
}

#' Interpolate a passive ramp's 6 channels onto an active window's angles
#' (match by ANGLE, so a 5-deg vs 0-deg start offset washes out). Requires
#' >= 50% angle overlap, else returns NULL (interpolation would be pure
#' extrapolation). Works identically for within- and cross-trial ramps.
#' @return list(F=3, T=3) angle-matched passive means, or NULL.
.mfv_interp_ramp_onto <- function(ang_act, ramp) {
  if (length(ang_act) == 0L || all(!is.finite(ang_act))) return(NULL)
  ang_p <- ramp$ang
  if (sum(is.finite(ang_p)) < 2L) return(NULL)
  rng <- range(ang_p, na.rm = TRUE)
  if (mean(ang_act >= rng[1L] & ang_act <= rng[2L], na.rm = TRUE) < 0.5) return(NULL)
  interp1 <- function(vp) {
    ok <- is.finite(ang_p) & is.finite(vp)
    if (sum(ok) < 2L) return(NA_real_)
    agg <- stats::aggregate(vp[ok], by = list(a = ang_p[ok]), FUN = mean)
    if (nrow(agg) < 2L || sum(is.finite(agg$a) & is.finite(agg$x)) < 2L) return(NA_real_)
    yi <- tryCatch(stats::approx(agg$a, agg$x, xout = ang_act, rule = 2)$y,
                   error = function(e) rep(NA_real_, length(ang_act)))
    mean(yi, na.rm = TRUE)
  }
  am <- list(F = vapply(ramp$F, interp1, numeric(1L)), T = vapply(ramp$T, interp1, numeric(1L)))
  if (!all(is.finite(c(am$F, am$T)))) return(NULL)
  am
}


# =============================================================================
# Sono (DS3) segment length: 40 Hz low-pass + L0-crossing sampler
# =============================================================================

#' Zero-phase Butterworth low-pass (filtfilt); returns x unchanged if `signal`
#' is unavailable or too few finite samples.
.mfv_butter_lp <- function(x, cutoff_hz, sample_rate_hz, order = 4L) {
  nyq <- sample_rate_hz / 2.0
  if (!is.finite(cutoff_hz) || cutoff_hz >= nyq) return(x)
  ok <- is.finite(x)
  if (sum(ok) < 20L || !requireNamespace("signal", quietly = TRUE)) return(x)
  filt <- signal::butter(order, cutoff_hz / nyq, type = "low")
  out <- x; out[ok] <- signal::filtfilt(filt, x[ok]); out
}

#' Load the RIGHT-muscle DS3 sono length (mm), aligned to ref_td, SMOOTHED +
#' FILTERED with a 40 Hz low-pass (PI-directed, 2026-07-18: de-staircases the
#' ~247 Hz DS3 output sampled at 1 kHz while preserving muscle-length dynamics).
#' Sono is the ONLY wired channel for the right muscle, so this is right-only.
#' @return list(mm = filtered length, L0_mm = resting length, raw_mm) or NULL.
mfv_load_sono_lp40 <- function(filename, ref_td, cutoff_hz = MFV_SONO_LP_HZ) {
  mm <- tryCatch(.read_sono_right_mm_aligned(filename), error = function(e) NULL)
  if (is.null(mm) || length(mm) < nrow(ref_td)) return(NULL)
  mm <- mm[seq_len(nrow(ref_td))]
  dt <- suppressWarnings(stats::median(diff(ref_td$t.s), na.rm = TRUE))
  sr <- if (is.finite(dt) && dt > 0) 1.0 / dt else 1000.0
  mm_lp <- .mfv_butter_lp(mm, cutoff_hz, sr)
  L0 <- .sono_reference_length_mm(ref_td$angle.deg, mm_lp)
  list(mm = mm_lp, L0_mm = L0, raw_mm = mm)
}

#' Sample the continuous force at the sono L0 (resting-length) crossing within a
#' ramp. Muscle SHORTENING = length DECREASE, so the ramp passes through L0 as
#' the muscle crosses its resting length; we take the FIRST sign change of
#' (sono_len - L0) and linearly interpolate both time and force there.
#' @return list(force_at_L0_N, crossing_t_abs) or NULL if no crossing.
.mfv_fv_l0_crossing <- function(t_abs, f_std, sono_len, L0) {
  ok <- is.finite(t_abs) & is.finite(f_std) & is.finite(sono_len)
  if (sum(ok) < 3L || !is.finite(L0)) return(NULL)
  t_abs <- t_abs[ok]; f_std <- f_std[ok]; d <- sono_len[ok] - L0
  o <- order(t_abs); t_abs <- t_abs[o]; f_std <- f_std[o]; d <- d[o]
  sgn <- sign(d)
  idx <- which(sgn[-length(sgn)] * sgn[-1L] < 0)
  if (length(idx) == 0L) return(NULL)
  i <- idx[1L]
  frac <- d[i] / (d[i] - d[i + 1L])
  list(force_at_L0_N = f_std[i] + frac * (f_std[i + 1L] - f_std[i]),
       crossing_t_abs = t_abs[i] + frac * (t_abs[i + 1L] - t_abs[i]))
}


# =============================================================================
# Shared per-step force solve + SIGN STANDARDIZATION + continuous ts + FV L0
# =============================================================================

#' Solve one step's vector muscle force from pre-computed active/passive channel
#' means, standardize the sign, build the continuous force-vs-time trace, and
#' (when sono is supplied for a right-recruited ramp) sample the FV L0 point.
#' Shared by the isometric loop and the isovelocity cross-trial batch so the
#' force physics + sign convention live in ONE place.
#'
#' @param act,pass list(F=3, T=3) active / passive channel means (N, N*m).
#' @param sono_ctx optional list(mm, L0_mm) for FV L0 crossing (right only).
#' @return list(row = named vector of step columns, ts = force-vs-time tibble or
#'   NULL, fv_l0 = one-row tibble or NULL).
.mfv_finalize_step <- function(act, pass, noise, category, s, td6, step_rows, arms,
                               trial_id, snr_min, deactivation_window_s,
                               baseline_pad_s, relaxation_s, sono_ctx = NULL) {
  dF <- act$F - pass$F
  dT <- act$T - pass$T

  # -- u_hat: empirical primary, geometric/longitudinal cross-check+fallback --
  emp <- uhat_empirical(dF)
  snr <- if (is.finite(noise) && noise > 0) emp$mag / noise else NA_real_
  geo <- if (identical(category, "isometric")) uhat_geometric(s$operating_point)
         else list(uhat = c(NA_real_, NA_real_, NA_real_), angle_xy_deg = NA_real_)
  use_emp <- all(is.finite(emp$uhat)) && is.finite(snr) && snr >= snr_min
  if (use_emp) {
    uhat <- emp$uhat; uhat_source <- "empirical"
  } else if (all(is.finite(geo$uhat))) {
    uhat <- geo$uhat; uhat_source <- "geometric_fallback"
  } else {
    uhat <- c(1.0, 0.0, 0.0); uhat_source <- "longitudinal_fallback"
  }

  fr     <- solve_muscle_force_from_wrench(dF, dT, uhat, arms$r_m_m, arms$d_m)
  uhat_ref <- if (identical(category, "isometric") && all(is.finite(geo$uhat))) geo$uhat else c(1.0, 0.0, 0.0)
  fr_ref <- solve_muscle_force_from_wrench(dF, dT, uhat_ref, arms$r_m_m, arms$d_m)

  # ==========================================================================
  # SIGN-CONVENTION STANDARDIZATION (left muscle == right muscle) -- PI-directed
  # 2026-07-18, REVISED 2026-07-18 per PI biomechanical correction (superseding
  # the earlier "always positive" draft below the line). The legacy z-torque
  # path folds sides with force_sign = rec_lidx * lidx_pos_motor, correct for
  # the BEND-DIRECTION-dependent z (lateral) torque; this 6-axis vector force
  # is dominated by the y (dorsoventral-offset) torque instead, whose polarity
  # from muscle tension alone is the SAME for both muscles (both pull their
  # posterior attachment anteriorly) -- applying the bend-direction force_sign
  # here would MIRROR-FLIP one side, leaving left/right with opposite signs.
  #
  # PI's correction: concentric contraction is NOT "more extreme" in an
  # absolute sense and eccentric "less extreme" -- concentric is more extreme
  # IN THE DIRECTION THE PASSIVE BASELINE TORQUE ALREADY POINTS at that angle
  # (the muscle's own pull REINFORCES the passive bending moment there), while
  # eccentric deviates OPPOSITE that passive direction (the muscle's pull
  # OPPOSES/reduces it). The correct sign reference is therefore the PASSIVE
  # torque's own direction at the matched angle/window -- not an arbitrary
  # per-step "make this step's own value positive" rule.
  #
  # tension_sign is fixed to sign(pass_moment) -- the sign of the RAW passive
  # torque projected onto the SAME r x u_hat axis used for the active moment
  # (fr$rxu / fr$eff_arm_m). Multiplying by that fixed reference means a step
  # whose active deviation shares the passive torque's sign (reinforcing /
  # concentric-consistent) comes out POSITIVE, and a step whose deviation
  # opposes it (opposing / eccentric-consistent) comes out NEGATIVE or small --
  # exactly the relative-to-passive relationship described by the PI, instead
  # of erasing it. When the passive projected moment is too small to
  # distinguish from baseline force noise (e.g. an isometric hold at ~0 deg,
  # where the passive reference direction is itself ambiguous), fall back to
  # the previous always-positive-per-step convention rather than manufacture a
  # sign from noise. (`noise`, a FORCE-domain SD, is used as an approximate --
  # FLAGGED -- floor for this force-equivalent moment quantity; it is not a
  # torque-domain noise estimate.)
  # ==========================================================================
  pass_moment <- if (is.finite(fr$eff_arm_m) && fr$eff_arm_m > 0) {
    sum(pass$T * fr$rxu) / (fr$eff_arm_m^2)
  } else NA_real_
  pass_moment_ambiguous <- !is.finite(pass_moment) || !is.finite(noise) || abs(pass_moment) <= noise
  if (pass_moment_ambiguous) {
    tension_sign <- if (is.finite(fr$F_moment_N) && fr$F_moment_N < 0) -1 else 1
    tension_relative_to_passive <- "ambiguous_fallback"
  } else {
    tension_sign <- sign(pass_moment)
    tension_relative_to_passive <- if (is.finite(fr$F_moment_N) && sign(fr$F_moment_N) == sign(pass_moment))
      "reinforcing" else "opposing"
  }

  row <- list(
    step_number = s$step_number,
    dFx = dF[1L], dFy = dF[2L], dFz = dF[3L], dTx = dT[1L], dTy = dT[2L], dTz = dT[3L],
    uhat_x = uhat[1L], uhat_y = uhat[2L], uhat_z = uhat[3L],
    uhat_angle_emp_deg = emp$angle_xy_deg, uhat_angle_geom_deg = geo$angle_xy_deg,
    activation_snr = snr, eff_arm_m = fr$eff_arm_m,
    muscle_force_vector_N      = tension_sign * fr$F_moment_N,
    muscle_force_vector_geom_N = tension_sign * fr_ref$F_moment_N,
    force_force_channel_N      = tension_sign * fr$F_force_N,
    force_yTorque_N            = tension_sign * fr$F_yT_N,
    force_zTorque_N            = tension_sign * fr$F_zT_N,
    uhat_source = uhat_source,
    pass_moment_N = pass_moment,
    tension_relative_to_passive = tension_relative_to_passive
  )

  # -- continuous per-sample projected force (standardized), + FV L0 sample --
  r <- c(0.0, arms$r_m_m, arms$d_m); rxu <- .mfv_cross(r, uhat)
  eff_arm <- sqrt(sum(rxu^2))
  ts <- NULL; fv_l0 <- NULL
  if (is.finite(eff_arm) && eff_arm > 0) {
    disp_win <- step_rows & td6$t.s >= (s$stim_t0_s - baseline_pad_s) &
      td6$t.s <= (s$stim_t1_s + relaxation_s)
    if (any(disp_win)) {
      Tx <- .mfv_col(td6, "xtorque"); Ty <- .mfv_col(td6, "ytorque"); Tz <- .mfv_col(td6, "ztorque")
      f_t <- tension_sign * (((Tx[disp_win] - pass$T[1L]) * rxu[1L]) +
                             ((Ty[disp_win] - pass$T[2L]) * rxu[2L]) +
                             ((Tz[disp_win] - pass$T[3L]) * rxu[3L])) / (eff_arm^2)
      t_abs <- td6$t.s[disp_win]
      cmode_val <- if ("contraction_mode" %in% names(s)) s$contraction_mode else NA_character_
      ts <- tibble::tibble(
        trial_id = trial_id, unit_id = paste0("step", s$step_number),
        muscle_side = s$muscle_side, contraction_mode = cmode_val,
        t_rel = t_abs - s$stim_t0_s, muscle_force_vector_N = f_t,
        stim_duration_s = s$stim_t1_s - s$stim_t0_s,
        activation_snr = snr
      )
      # FV L0 crossing (right muscle only -- sono is right-only). Sample force
      # where the sono length crosses L0 within the active ramp window.
      if (!is.null(sono_ctx) && identical(s$muscle_side, "right") &&
          identical(category, "isovelocity")) {
        act_disp <- t_abs >= s$stim_t0_s & t_abs <= (s$stim_t1_s + deactivation_window_s)
        cr <- .mfv_fv_l0_crossing(t_abs[act_disp], f_t[act_disp],
                                  sono_ctx$mm[disp_win][act_disp], sono_ctx$L0_mm)
        if (!is.null(cr)) {
          fv_l0 <- tibble::tibble(
            trial_id = trial_id, step_number = s$step_number, muscle_side = s$muscle_side,
            velocity_operating_point = s$operating_point,
            shortening_strain_pct = s$shortening_strain_pct,
            contraction_mode = cmode_val,
            force_at_L0_N = cr$force_at_L0_N, L0_mm = sono_ctx$L0_mm,
            crossing_t_rel_s = cr$crossing_t_abs - s$stim_t0_s, sono_confirmed = TRUE,
            activation_snr = snr)
        }
      }
    }
  }
  list(row = row, ts = ts, fv_l0 = fv_l0)
}

#' Empty per-step output tibble with the vector columns (NA-filled).
.mfv_empty_out_cols <- function(step_numbers) {
  tibble::tibble(
    step_number = step_numbers,
    dFx = NA_real_, dFy = NA_real_, dFz = NA_real_,
    dTx = NA_real_, dTy = NA_real_, dTz = NA_real_,
    uhat_x = NA_real_, uhat_y = NA_real_, uhat_z = NA_real_,
    uhat_angle_emp_deg = NA_real_, uhat_angle_geom_deg = NA_real_,
    uhat_source = NA_character_, activation_snr = NA_real_, eff_arm_m = NA_real_,
    muscle_force_vector_N = NA_real_, muscle_force_vector_geom_N = NA_real_,
    force_force_channel_N = NA_real_, force_yTorque_N = NA_real_, force_zTorque_N = NA_real_,
    pass_moment_N = NA_real_, tension_relative_to_passive = NA_character_,
    passive_source = NA_character_, angle_matched_ok = NA
  )
}

#' Write one .mfv_finalize_step() row into out_cols at position i.
.mfv_assign_row <- function(out_cols, i, fin, passive_source, angle_ok) {
  rw <- fin$row
  for (nm in c("dFx","dFy","dFz","dTx","dTy","dTz","uhat_x","uhat_y","uhat_z",
               "uhat_angle_emp_deg","uhat_angle_geom_deg","activation_snr","eff_arm_m",
               "muscle_force_vector_N","muscle_force_vector_geom_N","force_force_channel_N",
               "force_yTorque_N","force_zTorque_N","pass_moment_N")) {
    out_cols[[nm]][i] <- as.numeric(rw[[nm]])
  }
  out_cols$uhat_source[i]    <- rw[["uhat_source"]]
  out_cols$tension_relative_to_passive[i] <- rw[["tension_relative_to_passive"]]
  out_cols$passive_source[i] <- passive_source
  out_cols$angle_matched_ok[i] <- angle_ok
  out_cols
}

#' Assemble the uhat comparison table from a vector step_summary.
.mfv_uhat_tbl <- function(step_summary) {
  step_summary |>
    dplyr::filter(.data$muscle_side %in% c("left", "right"),
                  is.finite(.data$uhat_angle_emp_deg) | is.finite(.data$uhat_angle_geom_deg)) |>
    dplyr::transmute(.data$step_number, .data$muscle_side, .data$contraction_mode,
                     bend_angle_deg = .data$operating_point,
                     .data$uhat_angle_emp_deg, .data$uhat_angle_geom_deg,
                     .data$activation_snr, .data$uhat_source)
}


# =============================================================================
# ISOMETRIC: per-trial vector muscle force (step's own pre-stim baseline)
# =============================================================================

#' Compute per-step vector muscle force for one ISOMETRIC trial (each step's
#' own pre-stim window is the passive baseline). Isovelocity is handled by
#' compute_isovelocity_vector_batch() so it can borrow passive ramps
#' cross-trial within the same individual.
#' @return list(step_summary, uhat_tbl, force_ts, arms, snr_min) or NULL.
attach_vector_muscle_force <- function(res, filename, category,
                                       deactivation_window_s = MFV_DEACTIVATION_WINDOW_S,
                                       baseline_pad_s = 0.2, relaxation_s = 1.0,
                                       snr_min = MFV_UHAT_SNR_MIN,
                                       td6_override = NULL) {
  # td6_override: diagnostic hook only (e.g. testing an alternate F/T
  # calibration matrix against the same trial). NULL preserves the normal
  # embedded-calibration load path unchanged.
  td6 <- if (!is.null(td6_override)) td6_override else mfv_load_six_axis(filename, res$td)
  if (is.null(td6)) {
    cli::cli_warn("attach_vector_muscle_force: 6-axis wrench unavailable for {basename(filename)} -- vector path skipped")
    return(NULL)
  }
  arms <- {
    geom <- .mfv_read_geom(filename)
    resolve_muscle_moment_arms(geom$width_mm, geom$depth_mm_raw, geom$dvert_mm)
  }
  ss <- res$step_summary
  trial_id <- if (is.null(res$trial_id)) NA_character_ else res$trial_id
  out_cols <- .mfv_empty_out_cols(ss$step_number)
  ts_list <- list()

  for (i in seq_len(nrow(ss))) {
    s <- ss[i, ]
    if (!s$muscle_side %in% c("left", "right")) next
    if (!is.finite(s$stim_t0_s) || !is.finite(s$stim_t1_s)) next
    step_rows <- td6$step_number == s$step_number
    if (!any(step_rows)) next
    act_win  <- td6$t.s >= s$stim_t0_s & td6$t.s <= (s$stim_t1_s + deactivation_window_s)
    base_win <- td6$t.s >= s$t_pre_baseline_start_s & td6$t.s <= s$t_pre_baseline_end_s

    act  <- .mfv_window_means(td6, step_rows, act_win)
    pass <- .mfv_window_means(td6, step_rows, base_win)
    noise <- .mfv_baseline_force_noise(td6, step_rows, base_win)

    fin <- .mfv_finalize_step(act, pass, noise, category, s, td6, step_rows, arms,
                              trial_id, snr_min, deactivation_window_s, baseline_pad_s, relaxation_s)
    out_cols <- .mfv_assign_row(out_cols, i, fin, "static_baseline", NA)
    if (!is.null(fin$ts)) ts_list[[length(ts_list) + 1L]] <- fin$ts
  }

  step_summary <- dplyr::left_join(ss, out_cols, by = "step_number")
  list(step_summary = step_summary, uhat_tbl = .mfv_uhat_tbl(step_summary),
       force_ts = if (length(ts_list) > 0L) dplyr::bind_rows(ts_list) else tibble::tibble(),
       fv_l0 = tibble::tibble(), arms = arms, snr_min = snr_min)
}


# =============================================================================
# ISOVELOCITY: cross-trial (same-individual) batch with angle-matched passive
# + sono-confirmed FV L0 sampling
# =============================================================================

#' Compute vector muscle force for ALL isovelocity trials of ONE specimen
#' together, so an active ramp with no same-trial angle-overlapping passive can
#' borrow an angle-matched no-stim ramp from ANOTHER trial of the SAME
#' individual at the same signed velocity (PI-directed, 2026-07-18). Also
#' samples the sono-confirmed FV L0 crossing (right muscle).
#'
#' @param iso_inputs list of list(trial_id, filename, res) -- one per
#'   isovelocity trial of the specimen (res from analyze_isovelocity()).
#' @return list(step_summaries (named by trial_id), uhat_tbl, force_ts, fv_l0,
#'   snr_min) or NULL if nothing usable.
compute_isovelocity_vector_batch <- function(iso_inputs,
                                             deactivation_window_s = MFV_DEACTIVATION_WINDOW_S,
                                             baseline_pad_s = 0.2, relaxation_s = 1.0,
                                             snr_min = MFV_UHAT_SNR_MIN,
                                             velocity_tol = MFV_VELOCITY_MATCH_TOL) {
  if (length(iso_inputs) == 0L) return(NULL)

  # -- load each trial's 6-axis wrench + sono once ------------------------
  trials <- list()
  for (inp in iso_inputs) {
    td6 <- mfv_load_six_axis(inp$filename, inp$res$td)
    if (is.null(td6)) {
      cli::cli_warn("isovelocity batch: 6-axis wrench unavailable for {inp$trial_id} -- skipped")
      next
    }
    geom <- .mfv_read_geom(inp$filename)
    arms <- resolve_muscle_moment_arms(geom$width_mm, geom$depth_mm_raw, geom$dvert_mm)
    sono <- tryCatch(mfv_load_sono_lp40(inp$filename, inp$res$td), error = function(e) NULL)
    trials[[inp$trial_id]] <- list(trial_id = inp$trial_id, td6 = td6, arms = arms,
                                   sono = sono, ss = inp$res$step_summary)
  }
  if (length(trials) == 0L) return(NULL)

  # -- build the SAME-INDIVIDUAL passive-ramp library (no-stim ramps) ------
  # Every entry comes from this specimen's own isovelocity trials; keyed by the
  # step's signed commanded velocity (operating_point) so borrowing only ever
  # matches a ramp traversing the same angle range in the same direction.
  passive_library <- list()
  for (tr in trials) {
    ns <- .mfv_no_stim_steps(tr$td6, unique(tr$td6$step_number))
    for (sn in ns) {
      op <- tr$ss$operating_point[match(sn, tr$ss$step_number)]
      if (!is.finite(op)) next
      passive_library[[length(passive_library) + 1L]] <- list(
        trial_id = tr$trial_id, operating_point = op,
        ramp = .mfv_ramp_from_step(tr$td6, tr$td6$step_number == sn))
    }
  }

  step_summaries <- list(); ts_list <- list(); fv_list <- list()

  for (tr in trials) {
    td6 <- tr$td6; ss <- tr$ss; arms <- tr$arms
    no_stim_steps <- .mfv_no_stim_steps(td6, unique(td6$step_number))
    out_cols <- .mfv_empty_out_cols(ss$step_number)

    for (i in seq_len(nrow(ss))) {
      s <- ss[i, ]
      if (!s$muscle_side %in% c("left", "right")) next
      if (!is.finite(s$stim_t0_s) || !is.finite(s$stim_t1_s)) next
      step_rows <- td6$step_number == s$step_number
      if (!any(step_rows)) next
      act_win  <- td6$t.s >= s$stim_t0_s & td6$t.s <= (s$stim_t1_s + deactivation_window_s)
      base_win <- td6$t.s >= s$t_pre_baseline_start_s & td6$t.s <= s$t_pre_baseline_end_s
      ang_act  <- td6$enc.deg[step_rows & act_win]

      act   <- .mfv_window_means(td6, step_rows, act_win)
      noise <- .mfv_baseline_force_noise(td6, step_rows, base_win)

      # passive: within-trial angle match -> cross-trial (same specimen) angle
      # match at same signed velocity -> static pre-stim baseline fallback.
      passive_source <- NA_character_; angle_ok <- NA; pass <- NULL
      within_ns <- no_stim_steps[abs(
        ss$operating_point[match(no_stim_steps, ss$step_number)] - s$operating_point) < velocity_tol]
      if (length(within_ns) > 0L) {
        pass <- .mfv_interp_ramp_onto(ang_act, .mfv_ramp_from_step(td6, td6$step_number == within_ns[1L]))
        if (!is.null(pass)) { passive_source <- "angle_matched"; angle_ok <- TRUE }
      }
      if (is.null(pass)) {
        for (lib in passive_library) {
          if (identical(lib$trial_id, tr$trial_id)) next  # already tried within-trial
          if (abs(lib$operating_point - s$operating_point) >= velocity_tol) next
          cand <- .mfv_interp_ramp_onto(ang_act, lib$ramp)
          if (!is.null(cand)) { pass <- cand; passive_source <- "angle_matched_cross_trial"; angle_ok <- TRUE; break }
        }
      }
      if (is.null(pass)) {
        pass <- .mfv_window_means(td6, step_rows, base_win)
        passive_source <- "static_baseline_fallback"; angle_ok <- FALSE
      }

      sono_ctx <- if (!is.null(tr$sono)) tr$sono else NULL
      fin <- .mfv_finalize_step(act, pass, noise, "isovelocity", s, td6, step_rows, arms,
                                tr$trial_id, snr_min, deactivation_window_s,
                                baseline_pad_s, relaxation_s, sono_ctx = sono_ctx)
      out_cols <- .mfv_assign_row(out_cols, i, fin, passive_source, angle_ok)
      if (!is.null(fin$ts))    ts_list[[length(ts_list) + 1L]] <- fin$ts
      if (!is.null(fin$fv_l0)) fv_list[[length(fv_list) + 1L]] <- fin$fv_l0
    }

    step_summaries[[tr$trial_id]] <- dplyr::left_join(ss, out_cols, by = "step_number") |>
      dplyr::mutate(trial_id = tr$trial_id)
  }

  uhat_tbl <- purrr::map_dfr(names(step_summaries), function(tid)
    .mfv_uhat_tbl(step_summaries[[tid]]) |> dplyr::mutate(trial_id = tid))

  list(step_summaries = step_summaries,
       uhat_tbl = uhat_tbl,
       force_ts = if (length(ts_list) > 0L) dplyr::bind_rows(ts_list) else tibble::tibble(),
       fv_l0    = if (length(fv_list) > 0L) dplyr::bind_rows(fv_list) else tibble::tibble(),
       snr_min = snr_min)
}


# =============================================================================
# Activation ranking (PI A6) + specific-tension sanity check (PI feature 5)
# =============================================================================

#' Rank isometric steps by activation magnitude (||active-minus-passive force
#' vector||) and SNR, so the operator can eyeball/veto which steps are "clearly
#' active" before curves/u_hat are trusted. Prints an ASCII table.
rank_isometric_steps_by_activation <- function(step_summary_vec, trial_id = NA_character_,
                                               snr_min = MFV_UHAT_SNR_MIN) {
  r <- step_summary_vec |>
    dplyr::filter(.data$muscle_side %in% c("left", "right")) |>
    dplyr::mutate(force_delta_mag_N = sqrt(.data$dFx^2 + .data$dFy^2 + .data$dFz^2)) |>
    dplyr::arrange(dplyr::desc(.data$activation_snr)) |>
    dplyr::transmute(.data$step_number, .data$muscle_side,
                     bend_angle_deg = round(.data$operating_point, 2),
                     force_delta_mag_N = signif(.data$force_delta_mag_N, 3),
                     activation_snr = round(.data$activation_snr, 2),
                     muscle_force_vector_N = signif(.data$muscle_force_vector_N, 3),
                     uhat_source = .data$uhat_source,
                     clearly_active = .data$activation_snr >= snr_min)
  cli::cli_h3("Activation ranking{if (!is.na(trial_id)) paste0(' [', trial_id, ']') else ''} (SNR >= {snr_min} = clearly active; eyeball/veto per PI A6)")
  print(as.data.frame(r), row.names = FALSE)
  invisible(r)
}

#' Specific-tension sanity check: convert the single most-activated step's
#' vector muscle force through the muscle CSA and report N/cm^2 against the
#' typical vertebrate skeletal-muscle range (~15-30 N/cm^2). REPORTED, never
#' tuned to hit it (PI feature 5).
muscle_force_specific_tension_check <- function(step_summary_vec, muscle,
                                                trial_id = NA_character_) {
  cand <- step_summary_vec |>
    dplyr::filter(.data$muscle_side %in% c("left", "right"),
                  is.finite(.data$muscle_force_vector_N), is.finite(.data$activation_snr)) |>
    dplyr::arrange(dplyr::desc(.data$activation_snr))
  if (nrow(cand) == 0L) {
    cli::cli_alert_warning("specific-tension check: no finite vector-force step available")
    return(invisible(NULL))
  }
  s <- cand[1L, ]
  csa_cm2 <- muscle$csa_muscle_cm2
  st <- if (is.finite(csa_cm2) && csa_cm2 > 0) abs(s$muscle_force_vector_N) / csa_cm2 else NA_real_
  cli::cli_alert_info(
    "Specific-tension sanity{if (!is.na(trial_id)) paste0(' [', trial_id, ']') else ''}: step {s$step_number} ({s$muscle_side}, SNR {round(s$activation_snr,2)}) F={signif(s$muscle_force_vector_N,3)} N / CSA {signif(csa_cm2,3)} cm^2 = {signif(st,3)} N/cm^2 (typical vertebrate skeletal muscle ~15-30 N/cm^2; REPORTED, not tuned)"
  )
  invisible(tibble::tibble(trial_id = trial_id, step_number = s$step_number,
                           muscle_side = s$muscle_side, force_N = s$muscle_force_vector_N,
                           csa_muscle_cm2 = csa_cm2, specific_tension_Ncm2 = st))
}
