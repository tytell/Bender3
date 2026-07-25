# diag_dynamic_power_contraction_phase.R
# PROTOTYPE / DIAGNOSTIC ONLY -- does NOT feed any production figure.
#
# Follow-up to the 2026-07-24 whole-cycle + L/R sign-consistency fixes (see
# analysis_muscle_force_vector_log.md). After those two fixes, dynamic
# muscle power came out CONSISTENTLY negative-leaning on both L and R
# (no more artificial L-vs-R sign split) -- but still uniformly negative,
# which is physiologically odd for an actively driven work loop. PI's
# critical methodological point on why:
#
#   "including stim-only components of the active cycles will inevitably
#   miss much of the work that is performed by the residual forces of the
#   muscle ... a passively lengthening muscle gets stimmed (eccentric
#   contraction), then begins to shorten (as it is still stimmed, concentric
#   contraction), then continues producing force as it shortens even after
#   the stimulus is removed."
#
# The PRODUCTION path (calc_muscle_torque() + add_muscle_instantaneous())
# trusts the RAW SIGNED torque*velocity product's sign, sample by sample,
# over a window that mixes concentric, eccentric, and post-stim residual-
# decay phases together -- and never resolves whether the raw calibration
# sign convention (torque axis vs. commanded-angle axis) even agrees with
# "positive = this muscle producing shortening force" in the first place.
#
# This script instead ports the STRUCTURAL approach from the legacy scup
# codebase (~/Desktop/bender_projects/scupBender/2025-08-13_bender_
# functions.R::calculate_muscle_time(), lines 888-940): classify each
# sample's contraction type from the COMMANDED motion (not the possibly-
# ambiguous raw torque sign), then FORCE muscle_torque/power/work's sign
# from that classification (magnitude-fold + impose sign), instead of
# trusting the raw product's sign at all.
#
# Bender3 adaptation (does NOT reuse scup's hardcoded phase-quarter
# lookup -- that assumed a fixed nominal sinusoid shape with no real
# per-sample commanded-angle data; Bender3 has the real thing):
#   - "sideward" (scup) -> derived here directly from
#     sign(dist.rad) (dist.rad = diff of `pos.rad`, and `pos.rad` comes from
#     `angle.deg`, which 00_load_bender_flat.R populates from the HDF5
#     `angle_commanded_degree` dataset -- i.e. this IS the commanded
#     waveform, not a noisy measured signal; its sign is fully
#     deterministic/trustworthy by construction, unlike raw torque).
#   - crossed with `side` (scup) -> here, msc$stim, the EVENT-resolved
#     "L"/"R" muscle recruited for this sample (same resolution
#     .attach_dynamic_muscle_force() already uses).
#   - the concentric/eccentric decision rule reuses the EXACT SAME
#     rig-geometry convention already audited for isometric/isovelocity
#     (resolve_step_contraction(), muscle_geometry.R): bend_lidx =
#     sign(bend direction) * lidx_pos_motor; rec_lidx = lidx_left/lidx_right
#     for the recruited side; concentric when bend_lidx == rec_lidx. No new
#     geometry constant invented -- just the existing force_sign logic
#     applied per-sample from commanded velocity instead of per-step from a
#     signed operating_point.
#   - Because classification depends on COMMANDED bend direction and
#     event-resolved side (not on whether stim is literally "on" right
#     now), post-stim relaxation-tail samples get classified the same way
#     as the rest of that same directional half-cycle -- naturally covering
#     the PI's "continues producing force as it shortens even after the
#     stimulus is removed" case, with NO separate windowing rule needed.
#
# Output: a comparison table + CSV, RAW (trusted-sign, current production
# convention) vs. CLASSIFIED (magnitude-folded + phase-imposed sign) --
# does NOT modify 03_analyze.R / run_fv_fl_power_pipeline.R. If the
# classified numbers look physically sensible (net positive mean power,
# no L-vs-R artifact), that is the signal to port this into production
# next, as its own separate, explicitly-approved change.
#
# Run with:  Rscript R/diag_dynamic_power_contraction_phase.R

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(cli); library(rhdf5); library(readr)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")
src("00_load_bender_flat.R")
src("01_calibrate.R")
src("02_deconvolve.R")
src("muscle_geometry.R")
src("03_analyze.R")
src("parse_trial_filename.R")
src("plot_force_vs_time.R")       # .detect_stim_events(), RELAXATION_WINDOW_S
src("dynamic_trial_precondition.R")  # extract_bender_trial_num()

DATA_OUT_DIR <- file.path(.crittergripper_root(), "02_processed", "data_processed")
dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SPECIMEN_SUBFOLDERS <- c(bass16 = BASS16_RAW_SUBFOLDER, bass17 = BASS17_RAW_SUBFOLDER, bass18 = BASS18_RAW_SUBFOLDER)

# =============================================================================
# Duplicated helpers (deliberate duplication pattern, see
# diag_precondition_power_check.R / run_fv_fl_power_pipeline.R -- no shared
# module boundary exists yet). Keep in sync with those two copies.
# =============================================================================
.read_file_level_geometry <- function(filename) {
  h5 <- rhdf5::H5Fopen(filename, "H5F_ACC_RDONLY")
  on.exit(try(rhdf5::H5Fclose(h5), silent = TRUE), add = TRUE)
  m_attrs <- tryCatch(rhdf5::h5readAttributes(h5, "/metadata"), error = function(e) list())
  dbl1 <- function(v, default = NA_real_) {
    v <- suppressWarnings(as.numeric(v[1L]))
    if (length(v) == 0L || is.na(v)) default else v
  }
  local_body_width_mm  <- dbl1(m_attrs[["measurement_specimen_local_body_width_millimeter"]])
  local_body_height_mm <- dbl1(m_attrs[["measurement_specimen_local_body_height_millimeter"]])
  dclamp_mm            <- dbl1(m_attrs[["measurement_clamp_separation_millimeter"]])
  density_g_per_mm3    <- dbl1(m_attrs[["measurement_specimen_density_gram_per_cubic_millimeter"]])
  muscle <- compute_muscle_mass_and_csa(local_body_width_mm, local_body_height_mm, dclamp_mm, density_g_per_mm3)
  list(muscle = muscle,
       lidx_pos_motor = dbl1(m_attrs[["daq_specimen_lateral_index_on_positive_motor_side"]]),
       lidx_left      = dbl1(m_attrs[["daq_specimen_side_index_left"]]),
       lidx_right     = dbl1(m_attrs[["daq_specimen_side_index_right"]]))
}

#' Same fixed (2026-07-24) version as diag_precondition_power_check.R:
#' msc$muscle_torque.Nm stays RAW; msc$stim is the event-resolved "L"/"R"
#' side (NOT the sparse per-sample stim pulse column), needed here as the
#' `side` half of the concentric/eccentric classification below.
.attach_dynamic_muscle_force <- function(td, torque_col, lidx_pos_motor, lidx_left, lidx_right,
                                          relaxation_s = RELAXATION_WINDOW_S) {
  td$.row_id <- seq_len(nrow(td))
  td2 <- set_cycle_types(td)
  msc <- tryCatch(
    calc_muscle_torque(td2, torque_col = torque_col, include_all_active_samples = TRUE, phase_match_all_rows = TRUE),
    error = function(e) { cli::cli_warn("calc_muscle_torque failed: {conditionMessage(e)}"); NULL }
  )
  if (is.null(msc) || nrow(msc) == 0L || !(".row_id" %in% names(msc))) return(NULL)
  # Capture the RAW per-sample stim pulse ("0"/"L"/"R", the sparse literal
  # trigger column) BEFORE it gets overwritten below by the event-resolved
  # row_side -- needed as a QA breakdown (literal-pulse vs. tail-only
  # samples within the same event window; see Check 3 below), distinct
  # from diag_precondition_power_check.R's copy of this helper.
  msc$stim_literal <- as.character(msc$stim)
  events <- tryCatch(.detect_stim_events(td), error = function(e) tibble::tibble())
  row_side <- rep(NA_character_, nrow(td))
  if (nrow(events) > 0L) {
    events <- dplyr::arrange(events, .data$stim_t0_s)
    for (i in seq_len(nrow(events))) {
      e <- events[i, ]
      in_win <- td$t.s >= e$stim_t0_s & td$t.s <= (e$stim_t1_s + relaxation_s)
      row_side[in_win] <- e$muscle_side
    }
  }
  msc$stim <- row_side[msc$.row_id]
  if ("cycletype" %in% names(msc)) msc <- dplyr::filter(msc, .data$cycletype == "act")
  msc
}

.load_one <- function(f) {
  td <- load_bender_flat(f, do_filter = TRUE, loadtorques = "x")
  tau <- deconvolve_bender(f, hub_path = NULL, verbose = FALSE)
  N <- min(nrow(td), length(tau)); td <- td[seq_len(N), , drop = FALSE]
  td$torque_inertia_corrected_Nm <- tau[seq_len(N)]
  attr(td, "Filename") <- f
  td
}

# =============================================================================
# Contraction classification (the new piece this script tests)
# =============================================================================

#' Per-sample concentric/eccentric classification from COMMANDED bend
#' direction (sign(dist.rad), deterministic -- angle.deg is
#' angle_commanded_degree) crossed with the event-resolved recruited side
#' (msc$stim, "L"/"R"). Reuses resolve_step_contraction()'s (muscle_
#' geometry.R) exact convention: bend_lidx = sign(bend dir) * lidx_pos_motor,
#' rec_lidx = lidx_left/lidx_right; concentric when bend_lidx == rec_lidx,
#' eccentric when bend_lidx == -rec_lidx. NA for zero/near-zero velocity
#' samples (direction undefined) or samples with no resolved side.
#' 2026-07-25 PI-confirmed fix: `stim_side` ("L"/"R", DAQ-recorded at
#' acquisition time) and `daq_specimen_side_index_left`/`_right` (separate,
#' independently-set rig-geometry metadata) do NOT agree on which physical
#' muscle is "left" for single_finite/dynamic files. Verified directly
#' against the raw commanded angle.deg trace, 3/3 specimens (bass16/17/18,
#' same lidx_pos_motor=-1/lidx_left=-1/lidx_right=+1 attrs each): "L"-labeled
#' stim pulses fire at mean angle.deg ~+2.0 to +2.4 deg (i.e. straddling the
#' POSITIVE/lidx_left-side extremum), "R"-labeled pulses at ~-1.8 to -2.1 deg
#' (the lidx_right-side extremum) -- exactly backward from the PI-confirmed
#' intended design ("phase=0 = stim centered on the RECRUITED muscle's own
#' peak STRETCH", analysis_muscle_force_vector_log.md 2026-07-25 addendum).
#' A stim_side=="L" event's own peak stretch is therefore at the
#' lidx_right-side extremum -- i.e. row_side=="L" resolves to lidx_right,
#' not lidx_left. Swapped below (was lidx_left/lidx_right, unswapped).
.classify_contraction <- function(stim, dist.rad, lidx_pos_motor, lidx_left, lidx_right) {
  rec_lidx <- dplyr::case_when(
    stim == "L" ~ lidx_right,
    stim == "R" ~ lidx_left,
    .default = NA_real_
  )
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

#' Collects BOTH per-cycle and per-sample tables in one pass (map_dfr can
#' only return one tibble shape per call, and we need two -- explicit loop
#' instead of trying to bind_rows two different row shapes together).
.collect_dynamic_all <- function(specimen, source_dir) {
  manifest <- parse_trial_directory(source_dir)
  files <- manifest$fullpath[manifest$protocol == "dynamic"]
  cyc_list <- vector("list", length(files))
  samp_list <- vector("list", length(files))
  for (i in seq_along(files)) {
    f <- files[[i]]
    trial_id <- tools::file_path_sans_ext(basename(f))
    td <- tryCatch(.load_one(f), error = function(e) NULL)
    if (is.null(td) || all(as.character(td$stim) == "0")) next
    geom_f <- .read_file_level_geometry(f)
    msc <- tryCatch(
      .attach_dynamic_muscle_force(td, "torque_inertia_corrected_Nm", geom_f$lidx_pos_motor, geom_f$lidx_left, geom_f$lidx_right),
      error = function(e) { cli::cli_alert_danger("{trial_id}: {conditionMessage(e)}"); NULL }
    )
    if (is.null(msc) || nrow(msc) == 0L) next

    msc$contraction <- .classify_contraction(msc$stim, msc$dist.rad,
                                              geom_f$lidx_pos_motor, geom_f$lidx_left, geom_f$lidx_right)
    msc$insta_power_raw.W    <- msc$muscle_torque.Nm * (msc$dist.rad / msc$t.interval.s)
    msc$work_increment_raw.J <- msc$muscle_torque.Nm * msc$dist.rad
    msc$insta_power_classified.W <- dplyr::case_when(
      msc$contraction == "concentric" ~ abs(msc$insta_power_raw.W),
      msc$contraction == "eccentric"  ~ -abs(msc$insta_power_raw.W),
      .default = NA_real_
    )
    msc$work_increment_classified.J <- dplyr::case_when(
      msc$contraction == "concentric" ~ abs(msc$work_increment_raw.J),
      msc$contraction == "eccentric"  ~ -abs(msc$work_increment_raw.J),
      .default = NA_real_
    )
    msc$muscle_mass.kg <- geom_f$muscle$muscle_mass_kg

    cyc <- msc |>
      dplyr::filter(!is.na(.data$muscle_torque.Nm), !is.na(.data$pos.rad)) |>
      dplyr::group_by(dplyr::across(dplyr::any_of(c(
        "filename", "fishcode", "trial", "cycle", "duty", "phase", "freq.Hz", "curvature.invm"
      )))) |>
      dplyr::summarise(
        n_samples        = dplyr::n(),
        n_concentric     = sum(.data$contraction == "concentric", na.rm = TRUE),
        n_eccentric      = sum(.data$contraction == "eccentric",  na.rm = TRUE),
        frac_concentric  = .data$n_concentric / pmax(.data$n_concentric + .data$n_eccentric, 1L),
        work_raw.J              = sum(.data$work_increment_raw.J, na.rm = TRUE),
        work_classified.J       = sum(.data$work_increment_classified.J, na.rm = TRUE),
        avg_power_raw.W         = mean(.data$insta_power_raw.W, na.rm = TRUE),
        avg_power_classified.W  = mean(.data$insta_power_classified.W, na.rm = TRUE),
        muscle_mass.kg   = dplyr::first(.data$muscle_mass.kg),
        sides_present    = paste(sort(unique(stats::na.omit(as.character(.data$stim)))), collapse = "+"),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        avg_power_raw.Wkg        = .data$avg_power_raw.W / .data$muscle_mass.kg,
        avg_power_classified.Wkg = .data$avg_power_classified.W / .data$muscle_mass.kg,
        work_raw.Jkg             = .data$work_raw.J / .data$muscle_mass.kg,
        work_classified.Jkg      = .data$work_classified.J / .data$muscle_mass.kg
      )
    if (nrow(cyc) > 0L) {
      cyc$specimen  <- specimen
      cyc$trial_id  <- trial_id
      cyc$trial_num <- extract_bender_trial_num(trial_id)
      cyc_list[[i]] <- cyc
      cli::cli_alert_success("{trial_id}: {nrow(cyc)} active cycles (mean frac_concentric {round(mean(cyc$frac_concentric, na.rm=TRUE), 2)})")
    }

    samp_list[[i]] <- msc |>
      dplyr::filter(!is.na(.data$muscle_torque.Nm), is.finite(.data$dist.rad), .data$stim %in% c("L", "R")) |>
      dplyr::transmute(
        specimen = specimen, trial_id = trial_id,
        stim = .data$stim, contraction = .data$contraction,
        is_literal_pulse = .data$stim_literal == .data$stim,  # TRUE = within the literal sparse
        # stim-pulse train itself; FALSE = event-window tail sample (relaxation
        # decay or pre-onset baseline) attributed to this side but with no
        # literal pulse AT this exact sample.
        torque_raw.Nm = .data$muscle_torque.Nm,
        velocity_raw.rads = .data$dist.rad / .data$t.interval.s,
        insta_power_raw.W = .data$insta_power_raw.W,
        insta_power_classified.W = .data$insta_power_classified.W
      )
  }
  list(cyc = dplyr::bind_rows(cyc_list), samples = dplyr::bind_rows(samp_list))
}

cli::cli_h1("Collecting dynamic trials, bass16/17/18, RAW vs. CLASSIFIED (contraction-phase) power/work")
collected <- purrr::imap(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  .collect_dynamic_all(specimen, raw_source_dir(subfolder))
})
dyn_cycles  <- dplyr::bind_rows(purrr::map(collected, "cyc"))
dyn_samples <- dplyr::bind_rows(purrr::map(collected, "samples"))

write.csv(dyn_cycles,  file.path(DATA_OUT_DIR, "diag_dynamic_power_contraction_phase_percycle.csv"),  row.names = FALSE)
write.csv(dyn_samples, file.path(DATA_OUT_DIR, "diag_dynamic_power_contraction_phase_persample.csv"), row.names = FALSE)

# =============================================================================
# Empirical check 1 (per-sample, pooled): does the L-vs-R sign flip-flop
# (the fix #2 bug) stay fixed under classification, and does classification
# resolve the "why uniformly negative" question?
# =============================================================================
cli::cli_h1("Check 1: per-sample pooled means by recruited side (RAW vs. CLASSIFIED)")
chk1 <- dyn_samples |>
  dplyr::group_by(.data$stim) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_torque_raw.Nm     = mean(.data$torque_raw.Nm, na.rm = TRUE),
    mean_velocity_raw.rads = mean(.data$velocity_raw.rads, na.rm = TRUE),
    mean_power_raw.W       = mean(.data$insta_power_raw.W, na.rm = TRUE),
    pct_positive_raw        = 100 * mean(.data$insta_power_raw.W > 0, na.rm = TRUE),
    mean_power_classified.W = mean(.data$insta_power_classified.W, na.rm = TRUE),
    pct_concentric          = 100 * mean(.data$contraction == "concentric", na.rm = TRUE),
    .groups = "drop"
  )
print(as.data.frame(chk1))
write.csv(chk1, file.path(DATA_OUT_DIR, "diag_dynamic_power_contraction_phase_check1_bysample.csv"), row.names = FALSE)

# =============================================================================
# Empirical check 2 (per-cycle): whole-cycle mean/net power+work, RAW vs.
# CLASSIFIED, pooled and by specimen -- the number that would actually feed
# the Coughlin comparison figure if this is ported to production.
# =============================================================================
cli::cli_h1("Check 2: per-cycle mean power (W) and work (J), RAW vs. CLASSIFIED -- pooled")
chk2_pooled <- dyn_cycles |>
  dplyr::summarise(
    n_cycles = dplyr::n(),
    mean_frac_concentric      = mean(.data$frac_concentric, na.rm = TRUE),
    mean_avg_power_raw.W      = mean(.data$avg_power_raw.W, na.rm = TRUE),
    pct_cycles_positive_raw   = 100 * mean(.data$avg_power_raw.W > 0, na.rm = TRUE),
    mean_avg_power_classified.W = mean(.data$avg_power_classified.W, na.rm = TRUE),
    pct_cycles_positive_classified = 100 * mean(.data$avg_power_classified.W > 0, na.rm = TRUE),
    mean_work_raw.J           = mean(.data$work_raw.J, na.rm = TRUE),
    mean_work_classified.J    = mean(.data$work_classified.J, na.rm = TRUE)
  )
print(as.data.frame(chk2_pooled))

cli::cli_h1("Check 2b: per-cycle mean power (W), RAW vs. CLASSIFIED -- by specimen")
chk2_specimen <- dyn_cycles |>
  dplyr::group_by(.data$specimen) |>
  dplyr::summarise(
    n_cycles = dplyr::n(),
    mean_frac_concentric        = mean(.data$frac_concentric, na.rm = TRUE),
    mean_avg_power_raw.W        = mean(.data$avg_power_raw.W, na.rm = TRUE),
    mean_avg_power_classified.W = mean(.data$avg_power_classified.W, na.rm = TRUE),
    mean_avg_power_raw.Wkg        = mean(.data$avg_power_raw.Wkg, na.rm = TRUE),
    mean_avg_power_classified.Wkg = mean(.data$avg_power_classified.Wkg, na.rm = TRUE),
    .groups = "drop"
  )
print(as.data.frame(chk2_specimen))
write.csv(dplyr::bind_rows(dplyr::mutate(chk2_pooled, specimen = "ALL"), chk2_specimen),
          file.path(DATA_OUT_DIR, "diag_dynamic_power_contraction_phase_check2_percycle.csv"), row.names = FALSE)


# =============================================================================
# Empirical check 3: is frac_concentric/mean power driven by the LITERAL
# stim-pulse portion of the window vs. the tail (relaxation decay/pre-onset
# baseline) portion? Properly isolates "this side's own window, sub-phased"
# instead of "literal pulse vs. everything else" (the PI's flagged caveat on
# today's earlier, flawed version of this comparison).
# =============================================================================
cli::cli_h1("Check 3: literal stim-pulse samples vs. tail-only samples, by side")
chk3 <- dyn_samples |>
  dplyr::group_by(.data$stim, .data$is_literal_pulse) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_velocity_raw.rads = mean(.data$velocity_raw.rads, na.rm = TRUE),
    pct_concentric          = 100 * mean(.data$contraction == "concentric", na.rm = TRUE),
    mean_power_raw.W        = mean(.data$insta_power_raw.W, na.rm = TRUE),
    mean_power_classified.W = mean(.data$insta_power_classified.W, na.rm = TRUE),
    .groups = "drop"
  )
print(as.data.frame(chk3))
write.csv(chk3, file.path(DATA_OUT_DIR, "diag_dynamic_power_contraction_phase_check3_pulseVsTail.csv"), row.names = FALSE)

cli::cli_alert_success("diag_dynamic_power_contraction_phase.R complete -- CSVs in {DATA_OUT_DIR}/")
