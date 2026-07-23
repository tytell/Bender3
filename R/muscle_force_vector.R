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
# All active-minus-passive: isometric uses a RELAXATION-AWARE passive baseline
# (PI-directed 2026-07-22, "M2", .mfv_isometric_relaxation_passive()) -- the
# step's own quiescent samples fit the viscoelastic relaxation curve, evaluated
# inside the stim gap at each channel's active peak, replacing the stale static
# pre-stim window mean (see that helper's header + R/diag_isometric_passive_
# models.R). Isovelocity uses an ANGLE-matched no-stim ramp (match by enc.deg,
# not time) so the 5-deg-vs-0-deg start offset washes out. When no same-trial
# angle-overlapping passive exists, the ramp is borrowed CROSS-TRIAL but ONLY
# from the SAME individual (same pipeline run = one specimen) at the same
# signed velocity (compute_isovelocity_vector_batch()).
#
# SIGN CONVENTION (PI-directed, 2026-07-22 -- REPLACES the 2026-07-18
# "reinforcing/opposing relative to passive torque" rule below, which is now
# deleted, not just revised): empirical check across all 328 real finalized
# steps (3 fish, isometric + isovelocity-V0 + dynamic bookends) showed the
# reinforcing/opposing rule was MANUFACTURING negative values at some
# off-L0 angles -- the underlying wrench-derived quantities (F_moment_N,
# F_force_N, F_yT_N; muscle_force_vector_N/_geom_N, force_force_channel_N,
# force_yTorque_N) are NOT reliably side- or angle-mirrored to begin with
# (mostly positive with a real negative minority, concentrated in
# isovelocity -- not explained by side or by concentric/eccentric position),
# so there is no real ambiguity for a passive-relative rule to resolve, and
# no separate L0 special case is needed either: these five quantities are
# now reported in their RAW, AS-COMPUTED sign, with NO per-step sign
# decision at all.
#
# The ONE quantity that DOES need a fixed correction is force_zTorque_N: the
# same empirical check showed it cleanly MIRRORS by anatomical side (left
# always one sign, right always the other) across all 3 fish, for a reason
# unrelated to L0/concentric/eccentric -- exactly the SAME bend-direction
# mirroring the legacy zTorque-only path already corrects for via
# force_sign = rec_lidx * lidx_pos_motor (resolve_step_contraction(),
# muscle_geometry.R). force_zTorque_N applies that identical, FIXED
# per-side/per-trial correction (computed directly from `geom`'s
# lidx_left/lidx_right/lidx_pos_motor, not from any step-level decision) --
# see .mfv_finalize_step()'s zt_sign.
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

#' SNR + MAGNITUDE gated F0 (F/F0 denominator) from a set of L0 / V=0 rep rows.
#'
#' PI-directed 2026-07-22 (magnitude floor). The activation-SNR gate alone
#' (activation_snr >= snr_min) measures the force VECTOR's amplitude above the
#' baseline noise, but the F0 VALUE used as the F/F0 denominator is the
#' PROJECTED muscle_force_vector_N -- a rep can clear the SNR gate yet still
#' project to a near-zero, sign-unstable force (the bass17 isovelocity-V=0
#' case: activation_snr 3.3-5.5 but muscle_force_vector_geom_N ~ +-0.0005 N),
#' which detonates F/F0 (one point hit ~ -64). This ADDITIONALLY requires each
#' surviving rep's |force| to exceed its OWN baseline force-noise floor
#' (baseline_force_noise_N, the pre-stim ||force-channel|| SD), i.e. the L0
#' tension estimate must itself be above the noise it is measured against, not
#' merely be the projection of an above-noise vector. Reps failing EITHER gate
#' are dropped; F0 is the mean of survivors, or NA (normalized points become
#' NA, RAW points are untouched) when none qualify.
#' @param force,snr,noise equal-length rep vectors (the F0 quantity being
#'   normalized, its activation_snr, and its baseline_force_noise_N).
#' @return scalar gated F0 (mean of survivors) or NA_real_.
mfv_gate_f0 <- function(force, snr, noise, snr_min = MFV_UHAT_SNR_MIN) {
  keep <- is.finite(force) & is.finite(snr) & snr >= snr_min &
          is.finite(noise) & abs(force) >= noise
  if (!any(keep)) return(NA_real_)
  mean(force[keep], na.rm = TRUE)
}

# Ordered tier levels for mfv_confidence_tier() below -- the single shared
# vocabulary every SNR-gating site in the codebase should reuse (PI-directed
# 2026-07-22, "SNR-based confidence gating audit",
# analysis_muscle_force_vector_log.md). Never reorder without checking every
# consumer's factor level assumptions (alpha scales, filter %in% lists).
MFV_CONFIDENCE_TIERS <- c("confident", "confidently_small", "unstable_magnitude", "unconfirmable")

# Point transparency for the 4-tier confidence flag (REPLACES the old 2-level
# MFV_CONFIDENCE_ALPHA, which conflated "confidently small" with
# "unconfirmable" -- see mfv_confidence_tier()). Shared by every plot that
# alpha-flags points instead of dropping them (plot_muscle_force_vector.R,
# diag_isovelocity_hillcheck.R, plot_fatigue_timeline.R) so a real-but-small
# point reads visibly differently from a point that could be pure noise.
MFV_CONFIDENCE_ALPHA <- c(confident = 0.85, confidently_small = 0.55,
                          unstable_magnitude = 0.35, unconfirmable = 0.15)
names(MFV_CONFIDENCE_ALPHA) <- MFV_CONFIDENCE_TIERS

#' 4-level SNR x magnitude confidence tier (PI-directed 2026-07-22 audit --
#' generalizes mfv_gate_f0()'s two-condition test, which is ALREADY correct
#' for the F0 denominator, into a labeled factor every OTHER SNR-gating site
#' in the codebase can share). A ratio-only SNR gate (`snr >= snr_min`) cannot
#' distinguish "the noise floor is elevated" from "the true force is
#' genuinely small" -- both produce a low ratio. Crossing SNR-pass with an
#' ABSOLUTE magnitude-pass (the SAME `abs(force) >= k * noise` test
#' mfv_gate_f0 already uses) resolves the ambiguity:
#'   confident           snr_pass  & mag_pass  -- as trustworthy as before.
#'   confidently_small    !snr_pass & mag_pass  -- REAL, just quiet; a
#'                         ratio-only gate would misfile this as noise (the
#'                         conflation this function exists to fix).
#'   unstable_magnitude   snr_pass  & !mag_pass -- high activation on the
#'                         ||dF|| VECTOR, near-zero on the PROJECTED scalar
#'                         `force` -- the case mfv_gate_f0 was built to catch
#'                         for F0; this generalizes the same caution to any
#'                         OTHER site that reads a projected force.
#'   unconfirmable        !snr_pass & !mag_pass -- indistinguishable from
#'                         pure noise on EITHER test; keep the pipeline's
#'                         existing never-drop, alpha-flag-faintest policy.
#' Evidence: R/diag_snr_magnitude_audit.R + analysis_muscle_force_vector_log.md
#' (2026-07-22 entry) -- 42 real "confidently_small" points and 13
#' "unstable_magnitude" points exist across the 3-fish corpus outside the
#' F0-denominator context mfv_gate_f0 already protects.
#' @param force,snr,noise equal-length vectors (the force being judged, its
#'   activation_snr, and its baseline_force_noise_N).
#' @return factor, levels = MFV_CONFIDENCE_TIERS. Non-finite input on EITHER
#'   test -> "unconfirmable" (the most conservative label, never NA -- a
#'   point with missing SNR/noise data is not silently dropped from the
#'   classification, it is treated as not-yet-trustworthy).
mfv_confidence_tier <- function(force, snr, noise, snr_min = MFV_UHAT_SNR_MIN, k = 1.0) {
  snr_pass <- is.finite(snr) & snr >= snr_min
  mag_pass <- is.finite(force) & is.finite(noise) & abs(force) >= k * noise
  lvl <- dplyr::case_when(
    snr_pass  & mag_pass  ~ MFV_CONFIDENCE_TIERS[1L],
    !snr_pass & mag_pass  ~ MFV_CONFIDENCE_TIERS[2L],
    snr_pass  & !mag_pass ~ MFV_CONFIDENCE_TIERS[3L],
    .default = MFV_CONFIDENCE_TIERS[4L])
  factor(lvl, levels = MFV_CONFIDENCE_TIERS)
}

# Deactivation tail (s) appended to the active window, matching the rest of
# the pipeline's DEACTIVATION_WINDOW_S.
MFV_DEACTIVATION_WINDOW_S <- 0.5

# Sono (DS3) low-pass cutoff (Hz) for de-staircasing the segment-length signal
# before L0-crossing detection. 40 Hz sits below the DS3 internal update rate
# (~247 Hz) yet preserves muscle-length dynamics -- the value validated as best
# in diag_sono_smoothing.R (PI-directed, 2026-07-18).
MFV_SONO_LP_HZ <- 40.0

# ACTIVE-window scalar-force extraction (PI-directed, 2026-07-22, following
# review of muscleforcemethodcompare.png/muscleforcemethodsensitivity.png):
# narrow-window width (s) for "Method D" -- mean of RAW samples in a window
# this wide, centered on the ACTIVE window's own smoothed-trace peak time.
# See .mfv_window_peak_means().
MFV_PEAK_WINDOW_S <- 0.15

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
#' Sign: all four are in the RAW lab frame (they follow bend direction) --
#' see .mfv_finalize_step() for which of the four get a fixed correction
#' before being written to the step summary (as of 2026-07-22: only F_zT_N,
#' via the SAME force_sign = rec_lidx * lidx_pos_motor the legacy muscle_force_Nm
#' path already uses; F_moment_N/F_force_N/F_yT_N are reported raw).
#'
#' @param dF,dTau length-3 active-minus-passive force (N) / torque (N*m) deltas.
#' @return list(F_moment_N, F_force_N, F_yT_N, F_zT_N, eff_arm_m, rxu). `rxu`
#'   (the r x u_hat axis used for the moment projection) is returned so callers
#'   can project OTHER torque vectors (e.g. the raw passive torque, historically
#'   used for a passive-relative sign convention -- see git history/module
#'   header, superseded 2026-07-22) onto the exact same axis/arm without
#'   recomputing it.
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

#' Mean of each of the 6 channels over a time window on one step. Used for
#' the PASSIVE/baseline window (a steady reference has no "peak" to chase --
#' see .mfv_window_peak_means() for the ACTIVE window instead).
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

#' "Method D" extraction for the ACTIVE window (PI-directed, 2026-07-22,
#' replacing the plain full-window mean previously used here): for each of
#' the 6 channels independently, find that channel's own smoothed-trace peak
#' (via .smooth_trace_display_only(), plot_force_vs_time.R -- sourced before
#' this file in every real pipeline entry point), searched ONLY within
#' [stim_t0_s, stim_t1_s] (the true stim duration, EXCLUDING the post-stim
#' deactivation tail -- searching into the tail let a settling/recovery
#' transient there get mistaken for the peak, found empirically on the
#' isometric panel of muscleforcemethodcompare.png), then take the mean of
#' RAW (unsmoothed) samples in a window MFV_PEAK_WINDOW_S wide centered on
#' that peak's time, clipped to `win_mask` (the full active+deactivation
#' window actually available for this step). Validated 2026-07-22
#' (muscleforcemethodsensitivity.png, 92 real isometric steps, 3 fish):
#' tracks the plain full-window mean almost exactly for FLAT/sustained holds
#' (isometric), where the plain mean was already fine, while correctly
#' capturing genuine rise-then-decay TRANSIENTS (isovelocity's moving ramps)
#' that a plain full-window mean dilutes.
#' Falls back to the plain window mean whenever there are too few finite
#' samples to smooth/search reliably, OR whenever the search window itself
#' is SHORTER than `peak_window_s` (found empirically 2026-07-22 on the
#' dynamic L0 bookends -- their stim bursts are only ~0.05s, so a smoothed
#' trace still RISING at the search window's own edge gets its "peak" pinned
#' to that edge; the narrow window then balloons ~3x past the search
#' boundary into the deactivation tail, which can land on a large transient
#' unrelated to the tiny bookend pulse itself -- observed >2N spurious
#' outliers vs a real ~0.1-0.8N range on bass18's dynamic bookends). This
#' guard is about DURATION, not sample count: a short, fast burst needs the
#' plain mean; a long, sustained/ramping window is where Method D helps.
#' @param search_win_mask The narrower [stim_t0_s, stim_t1_s]-only mask used
#'   for the peak SEARCH; `win_mask` (full active+deactivation window) still
#'   bounds the narrow averaging window and is what the plain-mean fallback
#'   uses.
#' @param passive_pw Optional named list of full-length td6 POINTWISE passive
#'   vectors (keyed by .mfv_col() base, from .mfv_ramp_passive_pointwise()). When
#'   supplied for a channel, it is subtracted from that channel BEFORE peak
#'   finding + averaging, so the returned value is Method D of the
#'   ACTIVE-minus-passive DELTA (used by the isovelocity pointwise passive). NULL
#'   (default) preserves the raw-active behavior used everywhere else.
.mfv_window_peak_means <- function(td6, step_rows, win_mask, search_win_mask,
                                   peak_window_s = MFV_PEAK_WINDOW_S,
                                   passive_pw = NULL) {
  sel <- step_rows & win_mask
  search_sel <- step_rows & search_win_mask
  t <- td6$t.s
  get_peak_mean <- function(base) {
    v <- .mfv_col(td6, base)
    if (is.null(v)) return(NA_real_)
    if (!is.null(passive_pw) && !is.null(passive_pw[[base]])) v <- v - passive_pw[[base]]
    v_sel <- v[sel]
    full_mean <- mean(v_sel, na.rm = TRUE)
    if (sum(is.finite(v_sel)) < 3L) return(full_mean)
    v_search <- v[search_sel]; t_search <- t[search_sel]
    if (sum(is.finite(v_search)) < 3L) return(full_mean)
    if (diff(range(t_search, na.rm = TRUE)) < peak_window_s) return(full_mean)
    v_smooth <- .smooth_trace_display_only(v_search)
    smax <- max(v_smooth, na.rm = TRUE); smin <- min(v_smooth, na.rm = TRUE)
    if (!is.finite(smax) || !is.finite(smin)) return(full_mean)
    use_max <- abs(smax - full_mean) >= abs(smin - full_mean)
    t_peak <- t_search[if (use_max) which.max(v_smooth) else which.min(v_smooth)]
    narrow <- sel & t >= (t_peak - peak_window_s / 2) & t <= (t_peak + peak_window_s / 2)
    if (sum(is.finite(v[narrow])) < 1L) return(full_mean)
    mean(v[narrow], na.rm = TRUE)
  }
  list(
    F = c(get_peak_mean("xforce"),  get_peak_mean("yforce"),  get_peak_mean("zforce")),
    T = c(get_peak_mean("xtorque"), get_peak_mean("ytorque"), get_peak_mean("ztorque"))
  )
}

#' ISOMETRIC relaxation-aware passive baseline (PI-directed, 2026-07-22, "M2";
#' REPLACES the static pre-stim window mean previously used for the isometric
#' passive). Prototype settled in R/diag_isometric_passive_models.R: the static
#' pre-stim mean is sampled ~0.5-1 s BEFORE the active window on a still-relaxing
#' passive, so it is STALE -- the leftover viscoelastic creep scales with |bend|
#' and masquerades as the FL "arms" (a model-free zero-muscle control reproduced
#' 121-169% of the bass16/17 |bend| slope with NO activation). This instead fits
#' the passive's relaxation curve to the QUIESCENT (no-activation) samples --
#' the pre-stim window .. stim onset, plus after-deactivation .. post-baseline
#' end -- and evaluates it INSIDE the unobserved stim gap at EACH channel's own
#' active peak time (found exactly as .mfv_window_peak_means() finds it, so the
#' passive is subtracted at the same instant the active value is sampled). This
#' is the scalar-per-channel equivalent of pointwise passive subtraction; on the
#' validated fish it drops cor(F,|strain|) bass16 +0.19->+0.03, bass17
#' +0.93->+0.25, bass18 +0.57->+0.42 (bass18's monotonic rise is genuine, not
#' the artifact) WITHOUT inventing a spurious bell.
#'
#' Per channel it falls back to the plain static pre-stim mean (identical to the
#' old .mfv_window_means() passive) whenever the relaxation fit is unreliable
#' (fewer than 8 quiescent samples, search window shorter than the narrow peak
#' window, or the loess fails). CAVEAT (FLAGGED): the fit INTERPOLATES the
#' passive across the ~0.3 s stim gap it cannot observe, assuming the contraction
#' does not discontinuously perturb the passive; low-force fish sit below this
#' passive-drift floor and their FL SHAPE is not resolvable -- a magnitude/SNR
#' gate on that is a SEPARATE, still-pending decision (not made here).
#' `fits_T` (list of 3 fitted loess models, x/y/z torque order, NULL per
#' channel wherever that channel fell back to the static mean) is returned so
#' the continuous force-vs-time DISPLAY trace (.mfv_finalize_step()'s `ts`,
#' ported 2026-07-22) can subtract this SAME fitted relaxation curve
#' POINTWISE at every sample instead of the scalar peak-time value used for
#' the row -- one fit per channel, not refit twice.
#' @return list(F=3, T=3) passive channel values (N, N*m) at each channel's own
#'   active-peak time, used_relaxation flag, and fits_T (see above).
.mfv_isometric_relaxation_passive <- function(td6, step_rows, s, win_mask, search_win_mask,
                                              relaxation_s, base_mask,
                                              peak_window_s = MFV_PEAK_WINDOW_S) {
  t <- td6$t.s
  sel <- step_rows & win_mask
  search_sel <- step_rows & search_win_mask
  base_sel <- step_rows & base_mask
  quies <- step_rows & (
    (t >= s$t_pre_baseline_start_s & t < s$stim_t0_s) |
    (t >= (s$stim_t1_s + relaxation_s) & t <= s$t_post_baseline_end_s))
  tq <- t[quies]
  fit_relaxation <- function(v) {
    vq <- v[quies]; ok <- is.finite(tq) & is.finite(vq)
    if (sum(ok) < 8L) return(NULL)
    tqi <- tq[ok]; vqi <- vq[ok]
    idx <- seq(1L, length(tqi), by = max(1L, floor(length(tqi) / 400)))
    df <- data.frame(t = tqi[idx], g = vqi[idx])
    tryCatch(stats::loess(g ~ t, data = df, span = 0.6, degree = 2,
                          control = stats::loess.control(surface = "direct")),
             error = function(e) NULL)
  }
  predict_at <- function(fit, at) {
    if (is.null(fit)) return(NA_real_)
    as.numeric(tryCatch(stats::predict(fit, newdata = data.frame(t = at)),
                        error = function(e) NA_real_))[1L]
  }
  n_fit <- 0L
  fits_T <- list(x = NULL, y = NULL, z = NULL)
  get_pass <- function(base, torque_slot = NA_character_) {
    v <- .mfv_col(td6, base)
    if (is.null(v)) return(NA_real_)
    static_mean <- mean(v[base_sel], na.rm = TRUE)  # M0 fallback
    v_sel <- v[sel]
    if (sum(is.finite(v_sel)) < 3L) return(static_mean)
    v_search <- v[search_sel]; t_search <- t[search_sel]
    if (sum(is.finite(v_search)) < 3L) return(static_mean)
    if (diff(range(t_search, na.rm = TRUE)) < peak_window_s) return(static_mean)
    full_mean <- mean(v_sel, na.rm = TRUE)
    v_smooth <- .smooth_trace_display_only(v_search)
    smax <- max(v_smooth, na.rm = TRUE); smin <- min(v_smooth, na.rm = TRUE)
    if (!is.finite(smax) || !is.finite(smin)) return(static_mean)
    use_max <- abs(smax - full_mean) >= abs(smin - full_mean)
    t_peak <- t_search[if (use_max) which.max(v_smooth) else which.min(v_smooth)]
    fit <- fit_relaxation(v)
    pv <- predict_at(fit, t_peak)
    if (!is.finite(pv)) return(static_mean)
    n_fit <<- n_fit + 1L
    if (!is.na(torque_slot)) fits_T[[torque_slot]] <<- fit
    pv
  }
  list(
    F = c(get_pass("xforce"),  get_pass("yforce"),  get_pass("zforce")),
    T = c(get_pass("xtorque", "x"), get_pass("ytorque", "y"), get_pass("ztorque", "z")),
    used_relaxation = n_fit > 0L,
    fits_T = fits_T
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

#' POINTWISE angle-matched passive (PI-directed 2026-07-22; REPLACES the old
#' window-MEAN .mfv_interp_ramp_onto). Interpolate a velocity-matched no-stim
#' ramp's 6 channels onto EACH of a step's own samples BY ANGLE -- not collapsed
#' to a scalar window mean -- returning FULL-LENGTH td6 vectors (0 off-step),
#' keyed by the .mfv_col() base names ("xforce".."ztorque"). The caller then
#' subtracts these POINTWISE and runs Method D on the delta, so the passive is
#' removed at each sample's own angle (the isovelocity analog of the isometric
#' relaxation-fit-in-time). Because the active window SWEEPS through angle, the
#' passive varies by up to ~6 N across it; a single mean is a poor stand-in for
#' the passive at the Method-D peak's own angle (diag_isovelocity_passive_models).
#' Matches by ANGLE (a 5-deg vs 0-deg start offset washes out); requires >= 50%
#' of the ACTIVE-window angles to fall within the ramp's angle range, else NULL
#' (interpolation would be pure extrapolation) so the caller falls back. Works
#' identically for within- and cross-trial ramps (angle is comparable across the
#' same specimen's trials).
#' @return named list of 6 full-length td6 passive vectors, or NULL.
.mfv_ramp_passive_pointwise <- function(td6, step_rows, ang_act, ramp) {
  if (length(ang_act) == 0L || all(!is.finite(ang_act))) return(NULL)
  ang_p <- ramp$ang
  if (sum(is.finite(ang_p)) < 2L) return(NULL)
  rng <- range(ang_p, na.rm = TRUE)
  if (mean(ang_act >= rng[1L] & ang_act <= rng[2L], na.rm = TRUE) < 0.5) return(NULL)
  idx <- which(step_rows)
  ang_step <- td6$enc.deg[idx]
  n <- nrow(td6)
  bases <- c("xforce", "yforce", "zforce", "xtorque", "ytorque", "ztorque")
  chans <- c(ramp$F, ramp$T)  # same base order as `bases`
  out <- list()
  for (k in seq_along(bases)) {
    vp <- chans[[k]]
    ok <- is.finite(ang_p) & is.finite(vp)
    if (sum(ok) < 2L) return(NULL)
    agg <- stats::aggregate(vp[ok], by = list(a = ang_p[ok]), FUN = mean)
    if (nrow(agg) < 2L) return(NULL)
    yi <- tryCatch(stats::approx(agg$a, agg$x, xout = ang_step, rule = 2)$y,
                   error = function(e) NULL)
    if (is.null(yi) || any(!is.finite(yi))) return(NULL)
    full <- numeric(n); full[idx] <- yi
    out[[bases[k]]] <- full
  }
  out
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
#' @param geom File geometry from .mfv_read_geom(filename) -- needs
#'   lidx_left/lidx_right/lidx_pos_motor for force_zTorque_N's side
#'   correction (see zt_sign below). Same object attach_vector_muscle_force()/
#'   compute_isovelocity_vector_batch() already read for `arms`; passed
#'   through here rather than re-read.
#' @param sono_ctx optional list(mm, L0_mm) for FV L0 crossing (right only).
#' @param passive_curve_fits optional list(x, y, z) fitted loess relaxation
#'   models (from .mfv_isometric_relaxation_passive()$fits_T), used to
#'   subtract the passive POINTWISE (at each ts sample's own time) in the
#'   continuous force-vs-time trace below, instead of the single constant
#'   `pass$T[i]` used everywhere else. NULL (the default, and what
#'   isovelocity/dynamic callers pass implicitly by omission) preserves the
#'   original constant-passive ts behavior exactly. Per-channel NULL entries
#'   (e.g. a channel that fell back to the static mean) fall back to that
#'   channel's constant `pass$T[i]` too -- ts is never left with a channel-
#'   dropout gap.
#' @param passive_pw_T optional named list of full-length td6 POINTWISE passive
#'   torque vectors (keyed "xtorque"/"ytorque"/"ztorque", from
#'   .mfv_ramp_passive_pointwise()), used by the ISOVELOCITY path to subtract the
#'   angle-matched passive at each ts sample's own ANGLE. When supplied, the
#'   isovelocity caller also pre-subtracts it into `act` (so `pass` is zero and
#'   the scalar row is Method D of the delta), keeping the ts trace and the
#'   scalar row on the identical pointwise baseline. Takes precedence over
#'   `passive_curve_fits`; NULL (default) preserves prior behavior.
#' @return list(row = named vector of step columns, ts = force-vs-time tibble or
#'   NULL, fv_l0 = one-row tibble or NULL).
.mfv_finalize_step <- function(act, pass, noise, category, s, td6, step_rows, arms, geom,
                               trial_id, snr_min, deactivation_window_s,
                               baseline_pad_s, relaxation_s, sono_ctx = NULL,
                               passive_curve_fits = NULL, passive_pw_T = NULL) {
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
  # SIGN CONVENTION (PI-directed, 2026-07-22 -- see module header for the full
  # empirical rationale). Two DIFFERENT rules for two DIFFERENT quantities:
  #
  # 1) muscle_force_vector_N / muscle_force_vector_geom_N / force_force_channel_N
  #    / force_yTorque_N -- NO correction. Reported exactly as computed
  #    (fr$F_moment_N etc., RAW lab-frame sign). Confirmed empirically (328
  #    real steps, 3 fish, all 3 categories) NOT to be side- or L0-mirrored;
  #    the real negative minority (concentrated in isovelocity) is genuine
  #    signal, not a sign-convention artifact -- forcing it positive, or
  #    flipping it via a passive-torque comparison, would erase real data.
  #
  # 2) force_zTorque_N -- the ONE quantity confirmed to cleanly mirror by
  #    anatomical side (left one sign, right the other, consistently, across
  #    all 3 fish) -- gets the SAME fixed per-side correction the legacy
  #    zTorque-only path already uses: force_sign = rec_lidx * lidx_pos_motor
  #    (resolve_step_contraction(), muscle_geometry.R). Computed directly from
  #    `geom` (file-level attrs, present for isometric/isovelocity/dynamic
  #    alike) rather than depending on a step_summary column, so this works
  #    identically for dynamic L0 bookends' synthetic step rows too.
  # ==========================================================================
  rec_lidx <- if (identical(s$muscle_side, "left")) geom$lidx_left
             else if (identical(s$muscle_side, "right")) geom$lidx_right
             else NA_real_
  zt_sign <- rec_lidx * geom$lidx_pos_motor

  row <- list(
    step_number = s$step_number,
    dFx = dF[1L], dFy = dF[2L], dFz = dF[3L], dTx = dT[1L], dTy = dT[2L], dTz = dT[3L],
    uhat_x = uhat[1L], uhat_y = uhat[2L], uhat_z = uhat[3L],
    uhat_angle_emp_deg = emp$angle_xy_deg, uhat_angle_geom_deg = geo$angle_xy_deg,
    activation_snr = snr, eff_arm_m = fr$eff_arm_m,
    muscle_force_vector_N      = fr$F_moment_N,
    muscle_force_vector_geom_N = fr_ref$F_moment_N,
    force_force_channel_N      = fr$F_force_N,
    force_yTorque_N            = fr$F_yT_N,
    force_zTorque_N            = zt_sign * fr$F_zT_N,
    baseline_force_noise_N     = noise,
    uhat_source = uhat_source
  )

  # -- continuous per-sample projected force (standardized), + FV L0 sample --
  r <- c(0.0, arms$r_m_m, arms$d_m); rxu <- .mfv_cross(r, uhat)
  eff_arm <- sqrt(sum(rxu^2))
  ts <- NULL; fv_l0 <- NULL
  if (is.finite(eff_arm) && eff_arm > 0) {
    disp_win <- step_rows & td6$t.s >= (s$stim_t0_s - baseline_pad_s) &
      td6$t.s <= (s$stim_t1_s + relaxation_s)
    if (any(disp_win)) {
      # Matches muscle_force_vector_N's own RAW-sign convention above (no
      # per-step correction) -- this ts trace IS that same quantity's
      # per-sample time series.
      Tx <- .mfv_col(td6, "xtorque"); Ty <- .mfv_col(td6, "ytorque"); Tz <- .mfv_col(td6, "ztorque")
      t_abs <- td6$t.s[disp_win]
      # POINTWISE passive (PI-directed, 2026-07-22): subtract the moving passive
      # at EVERY sample so this display trace matches the corrected SCALAR
      # muscle_force_vector_N/_geom_N computed above. Two pointwise sources:
      #   ISOMETRIC -- per-channel relaxation loess fit (passive_curve_fits),
      #                evaluated at each sample's own TIME.
      #   ISOVELOCITY -- per-td6-row angle-matched ramp (passive_pw_T), indexed
      #                by disp_win (subtracted at each sample's own ANGLE).
      # Falls back to the constant pass$T[i] per channel (unchanged
      # dynamic/static-fallback behavior) when neither is supplied or a value is
      # non-finite -- ts is never left with a channel-dropout gap.
      passive_at <- function(slot, base, idx) {
        if (!is.null(passive_pw_T) && !is.null(passive_pw_T[[base]])) {
          pv <- passive_pw_T[[base]][disp_win]
          return(ifelse(is.finite(pv), pv, pass$T[idx]))
        }
        fit <- if (!is.null(passive_curve_fits)) passive_curve_fits[[slot]] else NULL
        if (is.null(fit)) return(rep(pass$T[idx], length(t_abs)))
        pv <- tryCatch(as.numeric(stats::predict(fit, newdata = data.frame(t = t_abs))),
                      error = function(e) rep(NA_real_, length(t_abs)))
        ifelse(is.finite(pv), pv, pass$T[idx])
      }
      passTx <- passive_at("x", "xtorque", 1L); passTy <- passive_at("y", "ytorque", 2L)
      passTz <- passive_at("z", "ztorque", 3L)
      f_t <- (((Tx[disp_win] - passTx) * rxu[1L]) +
             ((Ty[disp_win] - passTy) * rxu[2L]) +
             ((Tz[disp_win] - passTz) * rxu[3L])) / (eff_arm^2)
      cmode_val <- if ("contraction_mode" %in% names(s)) s$contraction_mode else NA_character_
      ts <- tibble::tibble(
        trial_id = trial_id, unit_id = paste0("step", s$step_number),
        muscle_side = s$muscle_side, contraction_mode = cmode_val,
        t_rel = t_abs - s$stim_t0_s, muscle_force_vector_N = f_t,
        stim_duration_s = s$stim_t1_s - s$stim_t0_s,
        activation_snr = snr, baseline_force_noise_N = noise
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
            activation_snr = snr, baseline_force_noise_N = noise)
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
    baseline_force_noise_N = NA_real_,
    passive_source = NA_character_, angle_matched_ok = NA
  )
}

#' Write one .mfv_finalize_step() row into out_cols at position i.
.mfv_assign_row <- function(out_cols, i, fin, passive_source, angle_ok) {
  rw <- fin$row
  for (nm in c("dFx","dFy","dFz","dTx","dTy","dTz","uhat_x","uhat_y","uhat_z",
               "uhat_angle_emp_deg","uhat_angle_geom_deg","activation_snr","eff_arm_m",
               "muscle_force_vector_N","muscle_force_vector_geom_N","force_force_channel_N",
               "force_yTorque_N","force_zTorque_N","baseline_force_noise_N")) {
    out_cols[[nm]][i] <- as.numeric(rw[[nm]])
  }
  out_cols$uhat_source[i]    <- rw[["uhat_source"]]
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
  geom <- .mfv_read_geom(filename)
  arms <- resolve_muscle_moment_arms(geom$width_mm, geom$depth_mm_raw, geom$dvert_mm)
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
    act_win    <- td6$t.s >= s$stim_t0_s & td6$t.s <= (s$stim_t1_s + deactivation_window_s)
    stim_win   <- td6$t.s >= s$stim_t0_s & td6$t.s <= s$stim_t1_s
    base_win <- td6$t.s >= s$t_pre_baseline_start_s & td6$t.s <= s$t_pre_baseline_end_s

    act  <- .mfv_window_peak_means(td6, step_rows, act_win, stim_win)
    # Passive = relaxation-aware baseline evaluated at each channel's active peak
    # (M2, PI-directed 2026-07-22), replacing the stale static pre-stim mean.
    # Fixes the SCALAR muscle_force_vector_N/_geom_N (what the FL/FV superplots
    # sample) AND, via pass$fits_T passed to .mfv_finalize_step() below (ported
    # 2026-07-22), the continuous force_ts DISPLAY trace, which now subtracts
    # this same fitted relaxation curve POINTWISE at every sample instead of a
    # single constant -- both use the identical M2 baseline, no display/scalar
    # mismatch.
    pass <- .mfv_isometric_relaxation_passive(td6, step_rows, s, act_win, stim_win,
                                              relaxation_s, base_win)
    noise <- .mfv_baseline_force_noise(td6, step_rows, base_win)

    fin <- .mfv_finalize_step(act, pass, noise, category, s, td6, step_rows, arms, geom,
                              trial_id, snr_min, deactivation_window_s, baseline_pad_s, relaxation_s,
                              passive_curve_fits = pass$fits_T)
    out_cols <- .mfv_assign_row(out_cols, i, fin, if (isTRUE(pass$used_relaxation)) "relaxation_fit" else "static_baseline", NA)
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
    trials[[inp$trial_id]] <- list(trial_id = inp$trial_id, td6 = td6, arms = arms, geom = geom,
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
      stim_win <- td6$t.s >= s$stim_t0_s & td6$t.s <= s$stim_t1_s
      base_win <- td6$t.s >= s$t_pre_baseline_start_s & td6$t.s <= s$t_pre_baseline_end_s
      ang_act  <- td6$enc.deg[step_rows & act_win]

      noise <- .mfv_baseline_force_noise(td6, step_rows, base_win)

      # passive: POINTWISE angle+velocity match (PI-directed 2026-07-22) --
      # within-trial -> cross-trial (same specimen) at the same SIGNED velocity
      # (operating_point, so same |v| AND direction = same angle range swept the
      # same way) -> static pre-stim baseline fallback. The matched no-stim ramp
      # is interpolated onto THIS step's samples BY ANGLE and subtracted
      # POINTWISE; Method D then runs on the ACTIVE-minus-passive delta. Replaces
      # the old window-MEAN collapse, which subtracted one scalar from the
      # Method-D peak while the passive swept up to ~6 N across the window (a
      # velocity-growing residual -- diag_isovelocity_passive_models). Velocity
      # matching itself is unchanged; only the mean-collapse became pointwise.
      passive_source <- NA_character_; angle_ok <- NA; passive_pw <- NULL
      within_ns <- no_stim_steps[abs(
        ss$operating_point[match(no_stim_steps, ss$step_number)] - s$operating_point) < velocity_tol]
      if (length(within_ns) > 0L) {
        passive_pw <- .mfv_ramp_passive_pointwise(td6, step_rows, ang_act,
                        .mfv_ramp_from_step(td6, td6$step_number == within_ns[1L]))
        if (!is.null(passive_pw)) { passive_source <- "angle_matched"; angle_ok <- TRUE }
      }
      if (is.null(passive_pw)) {
        for (lib in passive_library) {
          if (identical(lib$trial_id, tr$trial_id)) next  # already tried within-trial
          if (abs(lib$operating_point - s$operating_point) >= velocity_tol) next
          cand <- .mfv_ramp_passive_pointwise(td6, step_rows, ang_act, lib$ramp)
          if (!is.null(cand)) { passive_pw <- cand; passive_source <- "angle_matched_cross_trial"; angle_ok <- TRUE; break }
        }
      }
      if (is.null(passive_pw)) {
        # NEAREST same-sign velocity stim-off ramp (PI-directed 2026-07-22): the
        # experimental design provides stim-off ("passive") ramps only at a SUBSET
        # of the commanded velocities; steps at the remaining velocities have NO
        # exact-velocity ramp. Isovelocity ramps all sweep the SAME angle range
        # (same amplitude, differing speed), so the closest same-DIRECTION stim-off
        # ramp still traverses the matching angles and gives a position+velocity
        # baseline vastly better than the old static single-angle fallback (which
        # subtracted one pre-stim angle from a full swept ramp -> the large
        # wrong-signed FV overshoot). Same-sign guard preserves sweep direction.
        same_sign <- Filter(function(lib) is.finite(lib$operating_point) &&
                              sign(s$operating_point) != 0 &&
                              sign(lib$operating_point) == sign(s$operating_point),
                            passive_library)
        if (length(same_sign) > 0L) {
          ord <- order(abs(vapply(same_sign, function(l) l$operating_point, 0.0) -
                             s$operating_point))
          for (li in ord) {
            cand <- .mfv_ramp_passive_pointwise(td6, step_rows, ang_act, same_sign[[li]]$ramp)
            if (!is.null(cand)) {
              passive_pw <- cand; passive_source <- "angle_matched_nearest_v"; angle_ok <- TRUE; break
            }
          }
        }
      }
      if (!is.null(passive_pw)) {
        # delta = active - pointwise passive, THEN Method D (pass is zero: the
        # subtraction already happened per-sample inside .mfv_window_peak_means).
        act  <- .mfv_window_peak_means(td6, step_rows, act_win, stim_win, passive_pw = passive_pw)
        pass <- list(F = c(0.0, 0.0, 0.0), T = c(0.0, 0.0, 0.0))
        passive_pw_T <- passive_pw[c("xtorque", "ytorque", "ztorque")]
      } else {
        act  <- .mfv_window_peak_means(td6, step_rows, act_win, stim_win)
        pass <- .mfv_window_means(td6, step_rows, base_win)
        passive_source <- "static_baseline_fallback"; angle_ok <- FALSE
        passive_pw_T <- NULL
      }

      sono_ctx <- if (!is.null(tr$sono)) tr$sono else NULL
      fin <- .mfv_finalize_step(act, pass, noise, "isovelocity", s, td6, step_rows, arms, tr$geom,
                                tr$trial_id, snr_min, deactivation_window_s,
                                baseline_pad_s, relaxation_s, sono_ctx = sono_ctx,
                                passive_pw_T = passive_pw_T)
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
#'
#' `confidence_tier` (REPLACES the old ratio-only `clearly_active` boolean,
#' 2026-07-22 SNR-magnitude audit -- see mfv_confidence_tier()) shows the
#' human reviewer the magnitude distinction too: a step can be
#' "confidently_small" (SNR below threshold, but its projected force is
#' independently above its own noise floor -- a real, quiet contraction, not
#' obviously noise) rather than reading identically to "unconfirmable" the way
#' a bare SNR<3 boolean did.
rank_isometric_steps_by_activation <- function(step_summary_vec, trial_id = NA_character_,
                                               snr_min = MFV_UHAT_SNR_MIN) {
  r <- step_summary_vec |>
    dplyr::filter(.data$muscle_side %in% c("left", "right")) |>
    dplyr::mutate(force_delta_mag_N = sqrt(.data$dFx^2 + .data$dFy^2 + .data$dFz^2),
                  confidence_tier = mfv_confidence_tier(.data$muscle_force_vector_N,
                                                        .data$activation_snr,
                                                        .data$baseline_force_noise_N,
                                                        snr_min = snr_min)) |>
    dplyr::arrange(dplyr::desc(.data$activation_snr)) |>
    dplyr::transmute(.data$step_number, .data$muscle_side,
                     bend_angle_deg = round(.data$operating_point, 2),
                     force_delta_mag_N = signif(.data$force_delta_mag_N, 3),
                     activation_snr = round(.data$activation_snr, 2),
                     muscle_force_vector_N = signif(.data$muscle_force_vector_N, 3),
                     uhat_source = .data$uhat_source,
                     confidence_tier = .data$confidence_tier)
  cli::cli_h3("Activation ranking{if (!is.na(trial_id)) paste0(' [', trial_id, ']') else ''} (confidence tier = SNR x magnitude, see mfv_confidence_tier(); confident/confidently_small = magnitude-real; eyeball/veto per PI A6)")
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
