# muscle_geometry.R
# Muscle-side / contraction-mode resolution and predicted-strain geometry for
# segmented_finite (isometric, isovelocity) analysis.
#
# Scope (per PI decisions, 2026-07-15):
#  - Predicted strain uses HALF of local body WIDTH (lateral bending plane;
#    daq_bending_axis_specimen == "lateral" for this rig) minus an assumed
#    muscle depth, per PI direction: measurement_target_muscle_depth_millimeter
#    reads 0.0 in every bass16 file inspected (looks like an unset default,
#    not a real measurement) -- the PI directed using an ASSUMED 1 mm muscle
#    depth in that case instead of trusting the literal 0.
#  - Unilateral steps are folded by STIMULATED muscle side (left muscle -> one
#    curve, right muscle -> a separate curve), with contraction mode
#    (concentric/eccentric) resolved from bend direction relative to the
#    stimulated side: stim left + bend left = concentric (shortening);
#    stim left + bend right = eccentric (lengthening). Both modes contribute
#    to the SAME per-side curve (signed shortening axis), per PI direction.
#  - bilateral_simultaneous steps stimulate both sides at once; net torque
#    cannot be cleanly attributed to one muscle, so they are EXCLUDED from the
#    per-side FL/FV fits and reported separately as their own group.

library(dplyr)

# Assumed muscle depth (mm) beneath the body surface, used only when the
# exported measurement_target_muscle_depth_millimeter reads exactly 0 (the
# observed unset-default value in this corpus). This is a stated ASSUMPTION,
# not a measured quantity -- flagged wherever strain is reported.
DEFAULT_MUSCLE_DEPTH_MM <- 1.0

#' Resolve muscle depth to use for strain geometry.
#' Uses the file's own measurement_target_muscle_depth_millimeter when it is a
#' positive finite value; otherwise falls back to DEFAULT_MUSCLE_DEPTH_MM and
#' flags the fallback so callers can report the assumption explicitly.
resolve_muscle_depth_mm <- function(measured_depth_mm) {
  depth <- suppressWarnings(as.numeric(measured_depth_mm[1L]))
  if (length(depth) == 0L || !is.finite(depth) || depth <= 0) {
    list(depth_mm = DEFAULT_MUSCLE_DEPTH_MM, assumed = TRUE)
  } else {
    list(depth_mm = depth, assumed = FALSE)
  }
}

#' Predicted muscle strain (fraction) from instantaneous curvature.
#' strain = curvature (1/m) * r (m), r = half body width - muscle depth,
#' i.e. the radial distance from the neutral (mid-width) axis out to the
#' muscle's position beneath the lateral body surface.
#' @param curvature_invm Instantaneous curvature, 1/m (e.g. curve.invm).
#' @param local_body_width_mm measurement_specimen_local_body_width_millimeter.
#' @param muscle_depth_mm Resolved muscle depth (mm); see resolve_muscle_depth_mm().
#' @return list(strain = numeric vector (fraction), r_m = radial distance (m))
compute_predicted_strain <- function(curvature_invm, local_body_width_mm, muscle_depth_mm) {
  half_width_mm <- local_body_width_mm / 2.0
  r_m <- (half_width_mm - muscle_depth_mm) / 1000.0
  if (!is.finite(r_m) || r_m <= 0) {
    cli::cli_warn(
      "compute_predicted_strain: non-positive radial distance ({signif(r_m,4)} m) -- half body width smaller than assumed muscle depth; strain set to NA"
    )
    return(list(strain = rep(NA_real_, length(curvature_invm)), r_m = NA_real_))
  }
  list(strain = curvature_invm * r_m, r_m = r_m)
}

#' Resolve, per segmented step, which muscle (side) was stimulated and whether
#' the commanded bend direction shortens (concentric) or lengthens (eccentric)
#' that muscle -- generalized from the rig's own sign-convention attrs (works
#' for any daq_positive_motor_direction / lateral-index config, not just this
#' specimen's left-positive wiring).
#'
#' @param recruitment character vector, e.g. index_step_recruitment
#'   ("left_unilateral"/"right_unilateral"/"bilateral_simultaneous").
#' @param operating_point numeric vector, signed target angle (isometric) or
#'   velocity (isovelocity) for the step (index_step_operating_point).
#' @param lidx_pos_motor daq_specimen_lateral_index_on_positive_motor_side attr.
#' @param lidx_left daq_specimen_side_index_left attr.
#' @param lidx_right daq_specimen_side_index_right attr.
#' @return tibble(muscle_side, contraction_mode, shortening_value, force_sign)
#'   -- one row per step_number. shortening_value is the operating_point
#'   signed so that POSITIVE always means concentric/shortening for that
#'   muscle's own curve (standard Hill-curve convention), NA for
#'   bilateral/zero-magnitude steps.
#'   force_sign (+-1, NA for bilateral) re-expresses RAW lab-frame torque
#'   (which follows the commanded angle's sign, i.e. the SAME sign for both
#'   sides -- bending left and bending right necessarily read opposite raw
#'   torque regardless of which muscle is active) into the recruited muscle's
#'   OWN contraction-positive frame: force_sign = rec_lidx * lidx_pos_motor.
#'   Without this, left- and right-side FL/FV curves come out as mirror
#'   images of each other purely from the lab-frame sign flip between
#'   opposite bend directions -- not from any real physiological asymmetry.
#'   Callers must multiply their raw active-minus-passive torque by
#'   force_sign before folding left/right into one FL/FV axis.
resolve_step_contraction <- function(recruitment, operating_point,
                                      lidx_pos_motor, lidx_left, lidx_right) {
  recruitment <- as.character(recruitment)
  op          <- as.numeric(operating_point)

  side <- dplyr::case_when(
    grepl("^left_unilateral$",  recruitment) ~ "left",
    grepl("^right_unilateral$", recruitment) ~ "right",
    grepl("bilateral",          recruitment) ~ "bilateral",
    .default = NA_character_
  )
  rec_lidx <- dplyr::case_when(
    side == "left"  ~ lidx_left,
    side == "right" ~ lidx_right,
    .default = NA_real_
  )
  bend_lidx <- dplyr::if_else(op == 0, NA_real_, sign(op) * lidx_pos_motor)

  contraction_mode <- dplyr::case_when(
    side == "bilateral" ~ NA_character_,
    is.na(bend_lidx)     ~ "isometric_zero",
    bend_lidx == rec_lidx  ~ "concentric",
    bend_lidx == -rec_lidx ~ "eccentric",
    .default = NA_character_
  )
  shortening_value <- dplyr::case_when(
    side == "bilateral"                ~ NA_real_,
    contraction_mode == "isometric_zero" ~ 0.0,
    contraction_mode == "concentric"    ~ abs(op),
    contraction_mode == "eccentric"     ~ -abs(op),
    .default = NA_real_
  )

  # rec_lidx is NA for bilateral steps, so force_sign is NA there too --
  # matches shortening_value's bilateral-exclusion (folded per-side curves
  # never include bilateral steps; see analyze_isometric/analyze_isovelocity).
  force_sign <- rec_lidx * lidx_pos_motor

  tibble::tibble(
    muscle_side       = side,
    contraction_mode  = contraction_mode,
    shortening_value  = shortening_value,
    force_sign        = force_sign
  )
}


# =============================================================================
# Dynamic-trial (single_finite) contraction-phase classification
# =============================================================================

#' Map a dynamic-trial event-resolved recruited side ("L"/"R") to the
#' rig-geometry lidx_left/lidx_right constant for that side.
#'
#' TEMPORARY COMPENSATION (2026-07-25, PI-confirmed -- see
#' analysis_muscle_force_vector_log.md "2026-07-25 addendum -- root cause
#' found: make_stimuli() metadata/config-consistency bug"). For
#' single_finite/dynamic HDF5 files recorded before that Python-side bug is
#' fixed at the source (bender_functions.py::make_stimuli() -- follow-up
#' TODO, not yet done), the DAQ-recorded `stim_side` ("L"/"R") and the
#' independently-set rig-geometry attrs `daq_specimen_side_index_left`/
#' `_right` do NOT agree on which physical muscle is "left." Verified
#' directly against the raw commanded angle.deg trace, 3/3 specimens
#' (bass16/17/18): a `stim_side == "L"` event's own peak STRETCH sits at
#' the lidx_right extremum, not lidx_left -- backward from the intended
#' design. This function applies the SWAPPED mapping to compensate at the
#' point of use.
#'
#' SINGLE POINT OF TRUTH: every dynamic-trial caller that resolves a
#' recruited side to its rig-geometry lidx constant (force_sign display
#' correction in run_fv_fl_power_pipeline.R::.attach_dynamic_muscle_force(),
#' contraction-phase classification below) MUST go through this function --
#' do not re-inline the row_side -> lidx mapping elsewhere. The PI's planned
#' fix is to repair the raw HDF5 files directly (rather than leave this as a
#' permanent post-hoc compensation); once that lands, THIS is the one place
#' to update -- either restore the direct (unswapped) mapping, or gate it
#' behind a per-file acquisition-date/rig-software-version check if only
#' some files get repaired.
dynamic_recruited_side_to_lidx <- function(row_side, lidx_left, lidx_right) {
  dplyr::case_when(
    row_side == "L" ~ lidx_right,
    row_side == "R" ~ lidx_left,
    .default = NA_real_
  )
}

#' Per-sample concentric/eccentric contraction-phase classification for a
#' dynamic trial, from the COMMANDED bend direction (sign(dist.rad) --
#' dist.rad derives from `angle_commanded_degree` via 00_load_bender_flat.R,
#' so its sign is deterministic by construction, unlike raw torque) crossed
#' with the event-resolved recruited side (row_side, "L"/"R" -- the burst
#' window side from .detect_stim_events(), NOT the sparse per-sample stim
#' pulse column).
#'
#' Reuses the EXACT SAME convention resolve_step_contraction() (above)
#' already applies for isometric/isovelocity: concentric when
#' bend_lidx == rec_lidx (bend direction toward the recruited muscle's own
#' shortening extreme), eccentric when bend_lidx == -rec_lidx. NA when
#' velocity is zero/near-zero (direction undefined) or no side is resolved
#' for that sample.
#'
#' PI-approved 2026-07-25 port into production (see
#' analysis_muscle_force_vector_log.md) of the diagnostic prototyped in
#' R/diag_dynamic_power_contraction_phase.R -- consumed by
#' add_muscle_instantaneous()/summarize_muscle_cycles() (03_analyze.R) via
#' the `contraction` column a caller (currently only
#' run_fv_fl_power_pipeline.R::.attach_dynamic_muscle_force()) attaches to
#' the calc_muscle_torque() output before summarizing.
classify_dynamic_contraction <- function(row_side, dist.rad, lidx_pos_motor, lidx_left, lidx_right) {
  rec_lidx <- dynamic_recruited_side_to_lidx(row_side, lidx_left, lidx_right)
  bend_lidx <- dplyr::if_else(
    !is.finite(dist.rad) | dist.rad == 0, NA_real_,
    sign(dist.rad) * lidx_pos_motor
  )
  dplyr::case_when(
    !is.finite(rec_lidx) | !is.finite(bend_lidx) ~ NA_character_,
    bend_lidx == rec_lidx  ~ "concentric",
    bend_lidx == -rec_lidx ~ "eccentric",
    .default = NA_character_
  )
}


# =============================================================================
# Per-sample predicted strain, split into ACTIVE / PASSIVE columns
# =============================================================================

#' Attach predicted-strain columns to a loaded td tibble.
#' Adds: strain_pct (continuous, from curve.invm), is_active_sample (logical
#' mask), strain_active_pct (= strain_pct where active, else NA),
#' strain_passive_pct (= strain_pct where NOT active, else NA), and
#' muscle_depth_assumed_mm / muscle_depth_was_assumed (geometry provenance).
#'
#' @param td Tibble from load_bender_flat() (must have curve.invm).
#' @param local_body_width_mm measurement_specimen_local_body_width_millimeter.
#' @param measured_muscle_depth_mm measurement_target_muscle_depth_millimeter
#'   (raw, possibly the unset-default 0).
#' @param active_mask Logical vector, same length as nrow(td), TRUE where the
#'   sample counts as "active" for this protocol's windowing rule. If NULL,
#'   falls back to stim_type == "active" (best-effort default) when present.
# =============================================================================
# Muscle mass / cross-sectional-area estimate (PI-directed, 2026-07-16), for
# converting torque-based FL/FV/power results into mass- and area-specific
# properties (W/kg, N/cm^2) comparable to literature vertebrate red-muscle
# benchmarks -- to judge whether the recorded twitches are physiologically
# weak or strong, not just "positive" or "negative" in raw torque.
# =============================================================================

# Fraction of the total oval body-segment CSA/volume assumed to be the
# target (red) muscle -- a coarse literature-free ESTIMATE (not a
# dissected/weighed measurement), per PI direction.
DEFAULT_MUSCLE_VOLUME_FRACTION <- 0.03

# =============================================================================
# MEASURED red-muscle CSA (PI-directed, 2026-07-24): "I agree about the CSA
# [geometric estimate being the likely cause of the tension gap vs. Coughlin
# 2000]. Use the CSA in this table as an estimate for now" --
# 01_inputs/bass_csa_measurements/bass_csa_measurements.xlsx.
#
# Source: image-analysis CSA measurements (red / white / whole muscle,
# left/right) on FIVE cross-sectional "chunks" (each with an anterior and a
# posterior face) dissected from ONE reference specimen, "bass07" -- NOT one
# of bass16/17/18 used for the mechanics rig. The spreadsheet has NO chunk-
# length or bass07-total-length metadata, so there is no exact way to map a
# chunk to a %L position.
#
# WORKING ASSUMPTION (stated explicitly because it is unverified): with 5
# roughly evenly-spaced chunks spanning the body, chunks 2-4 (the middle
# three) are assumed to bracket the ~50%L clamp position the PI reported for
# bass16/17/18 (chunk 1 = most anterior, chunk 5 = most posterior, matching
# the whole-CSA taper: chunk1/2 ~31-37 cm^2 widest, chunk5 ~12-20 cm^2
# narrowest toward the tail).
#
# MEASURED_RED_MUSCLE_CSA_CM2 = mean of ALL red-left + red-right measurements
# (both faces) across chunks 2, 3 and 4 (n=11 of 17 total red measurements;
# chunk 1 and 5 excluded as likely outside the 50%L bracket) = 0.55 cm^2.
# Left/right are pooled (not averaged as a pair) to get a single "one side"
# CSA estimate, matching the single-side (right muscle) force this pipeline
# already divides by CSA (see diag_precondition_tension_vs_offset_isometric.R,
# summary_coughlin2000_bass_comparison.R). Range across those 11 values is
# large (0.257-1.50 cm^2) -- reflects a real anterior-posterior red-muscle
# gradient plus segmentation noise on thin tissue, not a single sharp value;
# treat MEASURED_RED_MUSCLE_CSA_CM2 as a rough "for now" point estimate, not
# a precise per-specimen measurement (bass16/17/18's OWN CSA was never
# imaged -- this substitutes a same-species reference value for all three).
MEASURED_RED_MUSCLE_CSA_CM2 <- 0.55

#' Estimate muscle mass (kg) and cross-sectional area (cm^2) from gross body
#' geometry: approximate the body cross-section at the test section as an
#' OVAL (local body width x height), multiply by the clamp-to-clamp
#' "test section" length for a volume, then take muscle_fraction (default
#' 3%) of the total CSA and of the total volume INDEPENDENTLY (i.e. muscle
#' CSA = 3% of total oval CSA; muscle mass = 3% of total oval volume x
#' specimen density) -- per PI direction, 2026-07-16.
#'
#' NOT a substitute for a dissected/weighed muscle sample -- this exists
#' purely to sanity-check whether the ALREADY-measured torque/power values
#' correspond to a physiologically weak or strong twitch once normalized by
#' a plausible muscle mass/area, against literature benchmarks (e.g.
#' Coughlin et al. 1996 scup red muscle ~114-152 W/kg steady-state power;
#' typical vertebrate skeletal-muscle specific tension ~15-30 N/cm^2).
#'
#' @param local_body_width_mm measurement_specimen_local_body_width_millimeter.
#' @param local_body_height_mm measurement_specimen_local_body_height_millimeter.
#' @param test_section_length_mm measurement_clamp_separation_millimeter (the
#'   clamp-to-clamp span that is actually bent/tested -- same dclamp_mm used
#'   elsewhere in this pipeline for the curvature -> strain geometry).
#' @param density_g_per_mm3 measurement_specimen_density_gram_per_cubic_millimeter
#'   (whole-body density from the file; used as a stand-in for muscle
#'   tissue density -- not a muscle-specific measurement).
#' @param muscle_fraction Fraction of total oval CSA/volume assumed to be the
#'   target muscle (default DEFAULT_MUSCLE_VOLUME_FRACTION, 3%).
#' @return list(csa_total_cm2, volume_total_cm3, csa_muscle_cm2,
#'   volume_muscle_cm3, muscle_mass_g, muscle_mass_kg, muscle_fraction) --
#'   all NA if any required geometry input is missing/non-positive.
compute_muscle_mass_and_csa <- function(local_body_width_mm, local_body_height_mm,
                                        test_section_length_mm, density_g_per_mm3,
                                        muscle_fraction = DEFAULT_MUSCLE_VOLUME_FRACTION) {
  w   <- suppressWarnings(as.numeric(local_body_width_mm[1L]))
  h   <- suppressWarnings(as.numeric(local_body_height_mm[1L]))
  l   <- suppressWarnings(as.numeric(test_section_length_mm[1L]))
  rho <- suppressWarnings(as.numeric(density_g_per_mm3[1L]))
  ok <- length(w) == 1L && is.finite(w) && w > 0 &&
        length(h) == 1L && is.finite(h) && h > 0 &&
        length(l) == 1L && is.finite(l) && l > 0 &&
        length(rho) == 1L && is.finite(rho) && rho > 0
  if (!ok) {
    return(list(csa_total_cm2 = NA_real_, volume_total_cm3 = NA_real_,
                csa_muscle_cm2 = NA_real_, volume_muscle_cm3 = NA_real_,
                muscle_mass_g = NA_real_, muscle_mass_kg = NA_real_,
                muscle_fraction = muscle_fraction))
  }

  csa_total_mm2     <- pi * (w / 2) * (h / 2)   # oval cross-section area
  volume_total_mm3  <- csa_total_mm2 * l        # test-section volume
  csa_muscle_mm2    <- muscle_fraction * csa_total_mm2
  volume_muscle_mm3 <- muscle_fraction * volume_total_mm3
  muscle_mass_g     <- volume_muscle_mm3 * rho  # rho is g/mm^3

  list(
    csa_total_cm2     = csa_total_mm2 / 100.0,
    volume_total_cm3  = volume_total_mm3 / 1000.0,
    csa_muscle_cm2    = csa_muscle_mm2 / 100.0,
    volume_muscle_cm3 = volume_muscle_mm3 / 1000.0,
    muscle_mass_g     = muscle_mass_g,
    muscle_mass_kg    = muscle_mass_g / 1000.0,
    muscle_fraction   = muscle_fraction
  )
}

#' Convert a Force-Length/Force-Velocity fit's torque-based outputs into
#' physical force/velocity/power and mass-/area-specific properties, using
#' the SAME muscle lever arm (r_m, meters) and test-section length
#' (dclamp_m, meters) already used to build that fit's own strain(-rate)
#' x-axis (compute_predicted_strain()/03_analyze.R shortening_strain_pct).
#'
#' Force0_N = F0_Nm / r_m (lever-arm torque -> force).
#' specific_tension_Ncm2 = Force0_N / csa_muscle_cm2.
#' FV fits only: the Hill fit's peak_power is F0(torque)*Vmax in a composite
#' Nm*(%/s) unit (NOT physical power, since %/s != rad/s) -- recovering real
#' power requires re-introducing the SAME dclamp_m/r_m scaling that produced
#' the %/s strain-rate axis in the first place (v_muscle = omega*r_m, and
#' strain_rate_pct_per_s = omega*r_m/dclamp_m*100, so
#' omega = Vmax_pct_per_s*dclamp_m/(100*r_m)); real muscle power at that
#' operating point is peak_power_composite * dclamp_m / (100 * r_m).
#' mass_specific_peak_power_Wkg = that / muscle_mass_kg.
#'
#' @param fit A model-fit result in the fit_fv_fl.R convention (list with a
#'   `status` field, "ok" or "failed") -- currently only
#'   fit_force_velocity_curve() output is passed in practice (FL no longer
#'   fits a model by default, PI direction 2026-07-16; see fit_fv_fl.R
#'   module header), but the `kind = "FL"` branch below is kept generic in
#'   case an FL model fit is explicitly requested again. Returned unchanged,
#'   just with extra fields added, if status != "ok" or any geometry input
#'   is non-finite.
#' @param r_m Muscle lever arm (meters) -- build_segmented_step_summary()'s r_m.
#' @param dclamp_m Test-section length (meters) -- geom$dclamp_mm/1000.
#' @param muscle Result of compute_muscle_mass_and_csa().
#' @param kind "FL" or "FV".
add_specific_properties_to_fit <- function(fit, r_m, dclamp_m, muscle, kind = c("FL", "FV")) {
  kind <- match.arg(kind)
  if (is.null(fit) || !identical(fit$status, "ok")) return(fit)
  ok_geom <- is.finite(r_m) && r_m > 0 && is.finite(muscle$csa_muscle_cm2) && muscle$csa_muscle_cm2 > 0
  if (!ok_geom) return(fit)

  fit$force0_N <- fit$F0 / r_m
  fit$specific_tension_Ncm2 <- fit$force0_N / muscle$csa_muscle_cm2
  fit$muscle_csa_cm2 <- muscle$csa_muscle_cm2

  if (kind == "FV" && is.finite(dclamp_m) && dclamp_m > 0 &&
      is.finite(muscle$muscle_mass_kg) && muscle$muscle_mass_kg > 0) {
    power_conv <- dclamp_m / (100.0 * r_m)   # Nm*(%/s) composite -> real W
    fit$peak_power_W <- fit$peak_power * power_conv
    fit$mass_specific_peak_power_Wkg <- fit$peak_power_W / muscle$muscle_mass_kg
    fit$muscle_mass_kg <- muscle$muscle_mass_kg
    fit$Vmax_m_per_s <- (fit$Vmax / 100.0) * dclamp_m
  }
  fit
}


attach_predicted_strain <- function(td, local_body_width_mm, measured_muscle_depth_mm,
                                     active_mask = NULL) {
  depth <- resolve_muscle_depth_mm(measured_muscle_depth_mm)
  strain_res <- if ("curve.invm" %in% names(td)) {
    compute_predicted_strain(td$curve.invm, local_body_width_mm, depth$depth_mm)
  } else {
    list(strain = rep(NA_real_, nrow(td)), r_m = NA_real_)
  }
  td$strain_pct <- strain_res$strain * 100.0

  if (is.null(active_mask)) {
    active_mask <- if ("stim_type" %in% names(td)) {
      tolower(trimws(as.character(td$stim_type))) == "active"
    } else {
      rep(FALSE, nrow(td))
    }
  }
  active_mask[is.na(active_mask)] <- FALSE

  td$is_active_sample     <- active_mask
  td$strain_active_pct    <- dplyr::if_else(active_mask, td$strain_pct, NA_real_)
  td$strain_passive_pct   <- dplyr::if_else(!active_mask, td$strain_pct, NA_real_)
  td$muscle_strain_r_m         <- strain_res$r_m
  td$muscle_depth_assumed_mm   <- depth$depth_mm
  td$muscle_depth_was_assumed  <- depth$assumed
  td
}
