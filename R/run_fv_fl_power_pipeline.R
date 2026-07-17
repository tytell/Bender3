# run_fv_fl_power_pipeline.R
# Orchestrates the full FV / FL / power-output extraction pipeline for a
# bender specimen corpus: load -> deconvolve -> dispatch-by-protocol -> analyze
# -> per-trial compound plot -> per-category summary plot -> verification.
# Specimen-specific SOURCE_DIR / OUTPUT_DIR default to the bass16 corpus but
# can be overridden via BENDER3_SOURCE_DIR / BENDER3_OUTPUT_DIR env vars to
# run the same pipeline against another specimen's raw-data directory.
#
# Reuses (does not duplicate): R/00_load_bender_flat.R, R/01_calibrate.R,
# R/02_deconvolve.R, R/03_analyze.R (set_cycle_types/calc_muscle_torque/
# summarize_muscle_cycles/analyze_isometric/analyze_isovelocity),
# R/muscle_geometry.R, R/fit_fv_fl.R, R/analyze_frequency_sweep.R,
# R/parse_trial_filename.R, R/plot_trial_compound.R, R/plot_summary_profiles.R.
#
# Run with:  Rscript R/run_fv_fl_power_pipeline.R
# Or for another specimen:
#   BENDER3_SOURCE_DIR=/path/to/other_specimen_bender \
#   BENDER3_OUTPUT_DIR=/path/to/other_specimen_figures \
#   Rscript R/run_fv_fl_power_pipeline.R

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(stringr)
  library(fs); library(ggplot2); library(cli); library(rhdf5)
})

## SOURCE_DIR / OUTPUT_DIR can be overridden via env vars so the same pipeline
## can be pointed at a different specimen's raw-data / figures directories
## without editing this file (e.g. Sys.setenv(BENDER3_SOURCE_DIR=...) or
## `BENDER3_SOURCE_DIR=... BENDER3_OUTPUT_DIR=... Rscript run_fv_fl_power_pipeline.R`).
SOURCE_DIR <- Sys.getenv(
  "BENDER3_SOURCE_DIR",
  "/Users/yjimenez/Library/CloudStorage/OneDrive-ProvidenceCollege/01_JimenezLab/01_PermanentArchive/bender_crittergripper/2026-07-14_bass16_bender"
)
OUTPUT_DIR <- Sys.getenv(
  "BENDER3_OUTPUT_DIR",
  "/Users/yjimenez/Library/CloudStorage/OneDrive-ProvidenceCollege/01_JimenezLab/02_ResearchHub/proj_crittergripper/figures/bass16_figures"
)
TRIAL_PLOT_DIR   <- file.path(OUTPUT_DIR, "trial_plots")
SUMMARY_PLOT_DIR <- file.path(OUTPUT_DIR, "summary_plots")

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("00_load_bender_flat.R")
src("01_calibrate.R")
src("02_deconvolve.R")
src("muscle_geometry.R")
src("fit_fv_fl.R")
src("03_analyze.R")
src("analyze_frequency_sweep.R")
src("parse_trial_filename.R")
src("plot_trial_compound.R")
src("plot_summary_profiles.R")
src("plot_fatigue_check.R")
src("plot_strain_validation.R")
src("plot_angle_sono_validation.R")
src("plot_force_vs_time.R")

fs::dir_create(TRIAL_PLOT_DIR, recurse = TRUE)
fs::dir_create(SUMMARY_PLOT_DIR, recurse = TRUE)

DEACTIVATION_WINDOW_S <- 0.5


#' Minimal file-level geometry reader for single_finite (dynamic /
#' frequency_sweep) trials -- same two attrs .read_segmented_step_geometry()
#' reads for segmented trials, without the index_step_* timing table those
#' protocols don't have. Also reads the three rig-geometry attrs
#' (daq_specimen_lateral_index_on_positive_motor_side/side_index_left/right)
#' needed for the same force_sign side-correction
#' resolve_step_contraction()/muscle_geometry.R applies for isometric/
#' isovelocity -- these are FILE-level (not step-level) attrs, so they exist
#' in dynamic files too even though dynamic has no index_step_* table.
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
  # Mass/CSA estimate (PI-directed, 2026-07-16), same recipe as
  # .read_segmented_step_geometry() -- see compute_muscle_mass_and_csa()
  # (muscle_geometry.R). Feeds td$muscle_mass.kg for the dynamic-protocol
  # mass-specific power calc (summarize_muscle_cycles()).
  muscle <- compute_muscle_mass_and_csa(local_body_width_mm, local_body_height_mm,
                                        dclamp_mm, density_g_per_mm3)
  list(
    local_body_width_mm  = local_body_width_mm,
    local_body_height_mm = local_body_height_mm,
    muscle_depth_mm_raw  = dbl1(m_attrs[["measurement_target_muscle_depth_millimeter"]]),
    dclamp_mm            = dclamp_mm,
    density_g_per_mm3    = density_g_per_mm3,
    muscle               = muscle,
    lidx_pos_motor       = dbl1(m_attrs[["daq_specimen_lateral_index_on_positive_motor_side"]]),
    lidx_left            = dbl1(m_attrs[["daq_specimen_side_index_left"]]),
    lidx_right           = dbl1(m_attrs[["daq_specimen_side_index_right"]])
  )
}

#' Broadcast per-cycle calc_muscle_torque() output back onto a full-length
#' td$muscle_force_Nm column via an explicit row-id join (calc_muscle_torque
#' preserves row order/identity through its group_by/filter/mutate chain, so
#' this recovers exact per-sample alignment without inventing new logic).
#'
#' Applies the SAME side-correction as resolve_step_contraction()
#' (muscle_geometry.R) -- force_sign = rec_lidx * lidx_pos_motor -- BEFORE
#' returning msc, so both the continuous td$muscle_force_Nm AND
#' summarize_muscle_cycles() (called on this same msc downstream) see
#' side-corrected values. Without this, dynamic muscle_force_Nm follows raw
#' lab-frame torque, which flips sign with bend direction regardless of
#' which muscle was actually recruited -- exactly the mirroring bug already
#' fixed for isometric/isovelocity (see resolve_step_contraction() docstring).
#'
#' rec_lidx is resolved from EVENT-level burst side (.detect_stim_events(),
#' also used by build_dynamic_force_timeseries()), NOT the row's own
#' instantaneous `stim` value. The raw `stim` signal is a SPARSE pulse train
#' (individual ~1-2-sample triggers ~13 ms apart -- see .detect_stim_events()
#' docstring): using the per-sample value directly left ~95% of samples
#' WITHIN an active burst (the "0" gaps between pulses) -- and 100% of the
#' relaxation tail after the last pulse -- with rec_lidx = NA and therefore
#' muscle_force_Nm = NA, which is what produced the "sudden cutoff" right
#' after (or even mid-) stimulation instead of a real decay tail. Each row
#' within [event_stim_t0_s, event_stim_t1_s + relaxation_s] now inherits
#' that event's side, so the whole burst plus its relaxation window gets a
#' signed value.
#'
#' Uses calc_muscle_torque(..., phase_match_all_rows = TRUE): the whole-CYCLE
#' act/pass gate (index_cycle_active) also used to hard-truncate rows the
#' instant a DESIGNED active block ended, even though the relaxation tail
#' (and, at high cycle frequencies, most of the tail) falls inside the very
#' next "pass"-labeled cycle. phase_match_all_rows = TRUE instead returns
#' every row -- both act- and pass-labeled cycles -- subtracted against the
#' same phase-matched passive baseline (see calc_muscle_torque() docstring).
.attach_dynamic_muscle_force <- function(td, torque_col, lidx_pos_motor, lidx_left, lidx_right,
                                          relaxation_s = RELAXATION_WINDOW_S) {
  td$.row_id <- seq_len(nrow(td))
  td2 <- set_cycle_types(td)
  msc <- tryCatch(
    calc_muscle_torque(td2, torque_col = torque_col, include_all_active_samples = TRUE,
                       phase_match_all_rows = TRUE),
    error = function(e) { cli::cli_warn("calc_muscle_torque failed: {conditionMessage(e)}"); NULL }
  )
  td$muscle_force_Nm <- NA_real_
  if (!is.null(msc) && nrow(msc) > 0L && ".row_id" %in% names(msc)) {
    events <- tryCatch(.detect_stim_events(td), error = function(e) tibble::tibble())
    row_side <- rep(NA_character_, nrow(td))
    if (nrow(events) > 0L) {
      # Chronological order so a later event's window wins any (rare) overlap
      # with the tail of an immediately-preceding same-cycle opposite-side event.
      events <- dplyr::arrange(events, .data$stim_t0_s)
      for (i in seq_len(nrow(events))) {
        e <- events[i, ]
        in_win <- td$t.s >= e$stim_t0_s & td$t.s <= (e$stim_t1_s + relaxation_s)
        row_side[in_win] <- e$muscle_side
      }
    }
    rec_lidx <- dplyr::case_when(
      row_side == "L" ~ lidx_left,
      row_side == "R" ~ lidx_right,
      .default = NA_real_
    )
    force_sign <- rec_lidx * lidx_pos_motor
    if (all(!is.finite(force_sign))) {
      cli::cli_warn("Dynamic force_sign side-correction unavailable (missing rig-geometry attrs or no stim events detected) -- muscle_force_Nm left in RAW lab-frame sign convention")
    } else {
      msc$muscle_torque.Nm <- force_sign[msc$.row_id] * msc$muscle_torque.Nm
      msc$stim <- row_side[msc$.row_id]  # event-resolved side, for downstream (summarize_muscle_cycles groups by `stim`)
    }
    td$muscle_force_Nm[msc$.row_id] <- msc$muscle_torque.Nm
  }
  # msc (the return value below) stays restricted to DESIGNED-active cycles
  # (cycletype == "act"), same scope as before this fix -- summarize_muscle_cycles()'s
  # per-cycle work/power table is a property of one designed activation, not
  # of the relaxation tail. Only td$muscle_force_Nm (the continuous per-sample
  # column, above) benefits from the full phase_match_all_rows = TRUE scope,
  # since that is what feeds the compound-plot and force_vs_time relaxation
  # display.
  msc_active <- if (!is.null(msc) && "cycletype" %in% names(msc)) dplyr::filter(msc, .data$cycletype == "act") else msc
  list(td = td, msc = msc_active)
}

#' Broadcast per-step Force_muscle (from build_segmented_step_summary) back
#' onto the segmented td's continuous samples, restricted to that step's own
#' active window (stim_t0 .. stim_t1 + deactivation window) -- the only
#' interval in which "muscle force" is actually a measured quantity for
#' isometric/isovelocity steps (elsewhere it is left NA, matching the
#' passive-only-trial NA-panel convention in build_trial_compound_plot()).
#'
#' @param force_col Column in step_summary to broadcast (default
#'   "muscle_force_Nm", the original static-baseline method). Pass
#'   "muscle_force_Nm_interp" to build the interpolated-baseline variant
#'   (2026-07-16) for the "_interpbaseline"-suffixed compound plot instead --
#'   always written to td$muscle_force_Nm either way, since
#'   build_trial_compound_plot() only knows that one column name.
.attach_segmented_muscle_force <- function(td, step_summary, deactivation_window_s = DEACTIVATION_WINDOW_S,
                                           force_col = "muscle_force_Nm") {
  td$muscle_force_Nm <- NA_real_
  for (i in seq_len(nrow(step_summary))) {
    s <- step_summary[i, ]
    if (!is.finite(s[[force_col]])) next
    idx <- td$step_number == s$step_number &
      td$t.s >= s$stim_t0_s & td$t.s <= (s$stim_t1_s + deactivation_window_s)
    td$muscle_force_Nm[idx] <- s[[force_col]]
  }
  td
}


# =============================================================================
# 1. Discover + validate trial manifest
# =============================================================================

manifest <- parse_trial_directory(SOURCE_DIR)
n_h5_files <- nrow(manifest)
cli::cli_h1("Bender FV/FL/power pipeline: {n_h5_files} trial file(s) discovered")
print(manifest |> select(trial_id, protocol, trial_number))

manifest$category <- dplyr::case_when(
  manifest$protocol == "isometric"       ~ "isometric",
  manifest$protocol == "isovelocity"     ~ "isovelocity",
  manifest$protocol == "frequency_sweep" ~ "frequency_sweep",
  manifest$protocol == "dynamic"         ~ "dynamic",
  .default = "unknown"
)
if (any(manifest$category == "unknown")) {
  cli::cli_warn("Unrecognized protocol(s) in filename, will be skipped: {manifest$filename[manifest$category=='unknown']}")
}


# =============================================================================
# 2. Per-trial loop: load -> deconvolve -> analyze -> compound plot
# =============================================================================

trial_reports         <- list()
isometric_steps_all   <- list()
isovelocity_steps_all <- list()
dynamic_cycles_all    <- list()
freqsweep_cycles_all  <- list()
strain_check_all      <- list()
angle_check_all       <- list()
sono_strain_check_all <- list()
force_ts_all          <- list(isometric = list(), isovelocity = list(), dynamic = list())
n_plots_written       <- 0L

for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, ]
  tid <- row$trial_id
  cli::cli_h2("[{i}/{nrow(manifest)}] {tid} (category={row$category})")

  report <- list(trial_id = tid, protocol = row$protocol, category = row$category,
                 status = "ok", note = NA_character_, fit_status = NA_character_)

  result <- tryCatch({
    if (row$category == "unknown") stop("unrecognized protocol in filename; skipped")

    td <- load_bender_flat(row$fullpath, do_filter = TRUE, loadtorques = "x")
    if (is.null(td) || nrow(td) == 0L) stop("load_bender_flat returned no data")

    tau <- deconvolve_bender(row$fullpath, hub_path = NULL, verbose = FALSE)
    if (is.null(tau)) stop("deconvolve_bender failed (see warning above)")
    N <- min(nrow(td), length(tau))
    td <- td[seq_len(N), , drop = FALSE]
    td$torque_inertia_corrected_Nm <- tau[seq_len(N)]
    attr(td, "Filename") <- row$fullpath

    # Read once, shared across all 4 protocol categories below: lidx_right/
    # lidx_pos_motor for the sono-strain right-muscle sign fold (see
    # attach_sono_strain(), plot_angle_sono_validation.R), and the raw
    # commanded-vs-measured angle pair for the angle validation figure --
    # neither needs any category-specific processing.
    geom_common <- .read_file_level_geometry(row$fullpath)
    angle_check_common <- tibble::tibble(
      trial_id = tid, protocol_family = row$category,
      angle_commanded = td$angle.deg, angle_measured = td$enc.deg,
      stim_state = .stim_window_state_label(td)
    )

    if (row$category %in% c("isometric", "isovelocity")) {
      analyze_fn <- if (row$category == "isometric") analyze_isometric else analyze_isovelocity
      res <- analyze_fn(td, filename = row$fullpath)
      # Adds strain_measured_step_pct (encoder-based, same sign fold as
      # shortening_strain_pct) so the FL/FV curve's own commanded-target
      # x-axis can be validated against actual motion (strain_validation_measured.png).
      steps <- attach_step_measured_strain(res$td, res$step_summary) |> dplyr::mutate(trial_id = tid, .before = 1L)

      # Isometric no longer fits a model by default (PI direction,
      # 2026-07-16 -- see 03_analyze.R analyze_isometric() docstring /
      # fit_fv_fl.R module header): res$fits is always an empty list for
      # isometric, so report_status just states that plainly. Isovelocity
      # still fits the (explicitly-requested) Hill model, so its per-side
      # status is surfaced as before.
      if (length(res$fits) == 0L) {
        report$fit_status <- "no model fit (connect-the-mean default; see summary plot)"
      } else {
        fit_bits <- purrr::imap_chr(res$fits, function(f, side) {
          sprintf("%s=%s%s", side, f$status, if (identical(f$status, "failed")) paste0(" (", f$reason, ")") else "")
        })
        report$fit_status <- paste(fit_bits, collapse = "; ")
        if (any(purrr::map_chr(res$fits, "status") == "failed")) {
          cli::cli_alert_warning("{tid}: {row$category} Hill fit FAILED for at least one side -- {report$fit_status}")
        }
      }

      td_plot <- .attach_segmented_muscle_force(res$td, res$step_summary)
      td_plot$t_stitched <- stitch_step_time_local(td_plot, row$fullpath)

      plt <- build_trial_compound_plot(
        td_plot, time_col = "t_stitched",
        title = tid, subtitle = sprintf("%s | muscle depth %s%.2g mm | fit: %s",
                                        row$category,
                                        if (isTRUE(res$muscle_depth_assumed)) "ASSUMED " else "measured ",
                                        res$muscle_depth_mm_used, report$fit_status)
      )

      # ADDITIONAL interpolated-baseline compound plot (2026-07-16,
      # isometric only -- see 03_analyze.R module header / build_segmented_
      # step_summary()'s passive_force_Nm_interp comment). Original `plt`
      # above is untouched; this is a second, separately-saved figure
      # (run_fv_fl_power_pipeline.R section 2 save block appends the
      # "_interpbaseline" filename suffix), not a replacement.
      plt_interp <- NULL
      if (row$category == "isometric") {
        td_plot_interp <- .attach_segmented_muscle_force(res$td, res$step_summary, force_col = "muscle_force_Nm_interp")
        td_plot_interp$t_stitched <- stitch_step_time_local(td_plot_interp, row$fullpath)
        plt_interp <- build_trial_compound_plot(
          td_plot_interp, time_col = "t_stitched",
          title = paste0(tid, " [interpolated baseline]"),
          subtitle = sprintf("%s | muscle depth %s%.2g mm | no model fit (connect-the-mean default)",
                             row$category,
                             if (isTRUE(res$muscle_depth_assumed)) "ASSUMED " else "measured ",
                             res$muscle_depth_mm_used)
        )
      }

      # Measured (encoder) vs. predicted (commanded-angle) strain, restricted
      # to samples where the muscle was ACTUALLY stimulated (real delivered
      # stim_side/stim, not just window membership -- .segmented_active_mask
      # marks the whole nominal window TRUE even for no-stim calibration steps).
      res$td <- attach_measured_strain(res$td)
      strain_check <- tryCatch({
        res$td |>
          dplyr::filter(as.character(.data$stim) != "0",
                        is.finite(.data$strain_pct), is.finite(.data$strain_measured_pct)) |>
          dplyr::transmute(trial_id = tid, protocol_family = row$category,
                            strain_pct = .data$strain_pct, strain_measured_pct = .data$strain_measured_pct)
      }, error = function(e) tibble::tibble())

      # Sonomicrometry (RIGHT muscle only -- the only wired sono channel in
      # this rig config) is the sole INDEPENDENT length measurement in the
      # rig; angle_check/sono_check use ALL samples (no active-stim
      # restriction), since they validate purely mechanical/geometric
      # tracking, not anything about stimulation.
      sono_stim_state <- .stim_window_state_label(res$td)
      sono_check <- tryCatch({
        attach_sono_strain(res$td, row$fullpath, geom_common$lidx_right, geom_common$lidx_pos_motor) |>
          dplyr::transmute(trial_id = tid, protocol_family = row$category,
                            strain_pred_commanded_right_pct = .data$strain_pred_commanded_right_pct,
                            strain_pred_encoder_right_pct   = .data$strain_pred_encoder_right_pct,
                            strain_sono_pct                 = .data$strain_sono_pct,
                            stim_state                      = sono_stim_state)
      }, error = function(e) tibble::tibble())

      # Isometric: nothing moves, so a scalar passive baseline is valid
      # pointwise. Isovelocity: the specimen moves continuously, so the
      # scalar baseline used for the step-level FV fit only cancels passive
      # bending torque IN THE MEAN -- the continuous trace needs the
      # velocity-matched no-stim step's own trace subtracted pointwise
      # instead (see build_isovelocity_force_timeseries() docstring).
      force_ts <- tryCatch({
        if (row$category == "isovelocity") {
          build_isovelocity_force_timeseries(res$td, res$step_summary, trial_id = tid)
        } else {
          build_segmented_force_timeseries(res$td, res$step_summary, trial_id = tid)
        }
      }, error = function(e) tibble::tibble())

      list(category = row$category, plot = plt, plot_interp = plt_interp, steps = steps,
           strain_check = strain_check, force_ts = force_ts,
           angle_check = angle_check_common, sono_check = sono_check)

    } else if (row$category == "frequency_sweep") {
      geom <- .read_file_level_geometry(row$fullpath)
      td <- attach_predicted_strain(td, geom$local_body_width_mm, geom$muscle_depth_mm_raw,
                                    active_mask = rep(FALSE, nrow(td)))
      td <- attach_measured_strain(td)
      freqsweep_stim_state <- .stim_window_state_label(td)
      sono_check <- tryCatch({
        attach_sono_strain(td, row$fullpath, geom$lidx_right, geom$lidx_pos_motor) |>
          dplyr::transmute(trial_id = tid, protocol_family = row$category,
                            strain_pred_commanded_right_pct = .data$strain_pred_commanded_right_pct,
                            strain_pred_encoder_right_pct   = .data$strain_pred_encoder_right_pct,
                            strain_sono_pct                 = .data$strain_sono_pct,
                            stim_state                      = freqsweep_stim_state)
      }, error = function(e) tibble::tibble())
      td$muscle_force_Nm <- NA_real_

      cyc <- analyze_frequency_sweep(td, torque_col = "torque_inertia_corrected_Nm")
      cyc <- dplyr::mutate(cyc, trial_id = tid, .before = 1L)
      n_failed <- sum(cyc$status == "failed")
      report$fit_status <- sprintf("%d/%d cycles ok, %d failed", sum(cyc$status == "ok"), nrow(cyc), n_failed)
      if (n_failed > 0L) cli::cli_alert_warning("{tid}: {n_failed} chirp cycle(s) FAILED stiffness/damping estimation")

      plt <- build_trial_compound_plot(
        td, time_col = "t.s", title = tid,
        subtitle = sprintf("frequency_sweep (passive-only) | %s", report$fit_status),
        muscle_force_note = "frequency_sweep trial is passive-only (no stimulation) -- Force_muscle undefined"
      )
      list(category = row$category, plot = plt, cycles = cyc, angle_check = angle_check_common, sono_check = sono_check)

    } else { # dynamic
      geom <- .read_file_level_geometry(row$fullpath)
      is_passive_only <- all(as.character(td$stim) == "0")
      active_mask <- as.character(td$stim) != "0"
      td <- attach_predicted_strain(td, geom$local_body_width_mm, geom$muscle_depth_mm_raw,
                                    active_mask = active_mask)
      # attach_measured_strain()/attach_sono_strain() run UNCONDITIONALLY
      # (even for passive-only trials) -- angle/strain tracking and the
      # sono-vs-curvature-geometry relationship are purely mechanical and
      # hold whether or not the muscle is being stimulated.
      td <- attach_measured_strain(td)
      dynamic_stim_state <- .stim_window_state_label(td)
      sono_check <- tryCatch({
        attach_sono_strain(td, row$fullpath, geom$lidx_right, geom$lidx_pos_motor) |>
          dplyr::transmute(trial_id = tid, protocol_family = "dynamic",
                            strain_pred_commanded_right_pct = .data$strain_pred_commanded_right_pct,
                            strain_pred_encoder_right_pct   = .data$strain_pred_encoder_right_pct,
                            strain_sono_pct                 = .data$strain_sono_pct,
                            stim_state                      = dynamic_stim_state)
      }, error = function(e) tibble::tibble())

      strain_check <- tibble::tibble()
      force_ts <- tibble::tibble()
      if (is_passive_only) {
        td$muscle_force_Nm <- NA_real_
        report$note <- "passive-only dynamic trial (stim_side all 'none') -- no active/passive split, no work-loop cycles"
        cyc_summary <- tibble::tibble()
      } else {
        strain_check <- tryCatch({
          td |>
            dplyr::filter(.data$is_active_sample, is.finite(.data$strain_pct), is.finite(.data$strain_measured_pct)) |>
            dplyr::transmute(trial_id = tid, protocol_family = "dynamic",
                              strain_pct = .data$strain_pct, strain_measured_pct = .data$strain_measured_pct)
        }, error = function(e) tibble::tibble())
        bc <- .attach_dynamic_muscle_force(td, "torque_inertia_corrected_Nm",
                                           geom$lidx_pos_motor, geom$lidx_left, geom$lidx_right)
        td <- bc$td
        force_ts <- tryCatch(build_dynamic_force_timeseries(td, trial_id = tid), error = function(e) tibble::tibble())
        # muscle_mass.kg (PI-directed, 2026-07-16 -- see .read_file_level_geometry()'s
        # compute_muscle_mass_and_csa() call) feeds summarize_muscle_cycles()'s
        # ALREADY-EXISTING but previously-dormant mass-normalization (avg_power.Wkg,
        # work.Jkg, peak_power.Wkg, Coughlin comparison) -- broadcast as a constant
        # (same specimen/geometry for every cycle in this trial).
        if (!is.null(bc$msc) && nrow(bc$msc) > 0L) bc$msc$muscle_mass.kg <- geom$muscle$muscle_mass_kg
        cyc_summary <- if (!is.null(bc$msc) && nrow(bc$msc) > 0L) {
          summarize_muscle_cycles(bc$msc) |> dplyr::mutate(trial_id = tid, .before = 1L)
        } else {
          tibble::tibble()
        }
        if (nrow(cyc_summary) == 0L) {
          report$note <- "stimulated trial but 0 active cycles resolved by calc_muscle_torque"
          cli::cli_alert_warning("{tid}: {report$note}")
        }
        report$fit_status <- sprintf("%d active cycles summarized", nrow(cyc_summary))
      }

      plt <- build_trial_compound_plot(
        td, time_col = "t.s", title = tid,
        subtitle = sprintf("dynamic%s | %s", if (is_passive_only) " (passive-only)" else "", report$fit_status %||% ""),
        muscle_force_note = if (is_passive_only) "passive-only dynamic trial -- no active stimulation, Force_muscle undefined" else NULL
      )
      list(category = row$category, plot = plt, cycles = cyc_summary, strain_check = strain_check, force_ts = force_ts,
           angle_check = angle_check_common, sono_check = sono_check)
    }
  }, error = function(e) {
    report$status <<- "error"
    report$note   <<- conditionMessage(e)
    cli::cli_alert_danger("{tid}: FAILED -- {conditionMessage(e)}")
    NULL
  })

  if (!is.null(result)) {
    out_path <- file.path(TRIAL_PLOT_DIR, paste0(tid, ".png"))
    ggplot2::ggsave(out_path, result$plot, width = 9, height = 9, dpi = 150)
    n_plots_written <- n_plots_written + 1L

    # ADDITIONAL interpolated-baseline plot (isometric only, see plot_interp
    # construction above) -- saved alongside, NOT counted in n_plots_written
    # (the verification check below requires that count to equal the
    # trial-file count for the ORIGINAL method's plots).
    if (!is.null(result$plot_interp)) {
      out_path_interp <- file.path(TRIAL_PLOT_DIR, paste0(tid, "_interpbaseline.png"))
      ggplot2::ggsave(out_path_interp, result$plot_interp, width = 9, height = 9, dpi = 150)
    }

    if (result$category == "isometric")   isometric_steps_all[[tid]]   <- result$steps
    if (result$category == "isovelocity") isovelocity_steps_all[[tid]] <- result$steps
    if (result$category == "dynamic" && nrow(result$cycles) > 0L) dynamic_cycles_all[[tid]] <- result$cycles
    if (result$category == "frequency_sweep") freqsweep_cycles_all[[tid]] <- result$cycles
    if (!is.null(result$strain_check) && nrow(result$strain_check) > 0L) strain_check_all[[tid]] <- result$strain_check
    if (!is.null(result$angle_check) && nrow(result$angle_check) > 0L) angle_check_all[[tid]] <- result$angle_check
    if (!is.null(result$sono_check) && nrow(result$sono_check) > 0L) sono_strain_check_all[[tid]] <- result$sono_check
    if (!is.null(result$force_ts) && nrow(result$force_ts) > 0L && result$category %in% names(force_ts_all)) {
      force_ts_all[[result$category]][[tid]] <- result$force_ts
    }
  }

  trial_reports[[tid]] <- report
}

trial_report_tbl <- dplyr::bind_rows(trial_reports)


# =============================================================================
# 3. Verification: individual-plot count must equal trial-file count
# =============================================================================

cli::cli_h1("Verification: per-trial plot count")
n_ok <- sum(trial_report_tbl$status == "ok")
cli::cli_alert_info("{n_h5_files} trial files discovered; {n_ok} processed ok; {n_plots_written} compound plots written")
if (n_plots_written != n_ok) {
  cli::cli_abort("Plot-count mismatch: {n_plots_written} plots written but {n_ok} trials reported ok -- investigate before trusting outputs")
}
if (any(trial_report_tbl$status == "error")) {
  cli::cli_alert_danger("{sum(trial_report_tbl$status=='error')} trial(s) FAILED entirely (no plot produced):")
  print(dplyr::filter(trial_report_tbl, status == "error") |> dplyr::select(trial_id, note))
}


# =============================================================================
# 4. Per-category summary plots
# =============================================================================

cli::cli_h1("Building per-category summary plots")

# Specimen-level geometry (PI-directed, 2026-07-16) -- local body
# width/height, test-section (clamp-to-clamp) length, and whole-body
# density are specimen constants, identical across every trial file for
# this specimen, so reading them once from the FIRST discovered file is
# equivalent to reading them per-trial. Feeds add_specific_properties_to_fit()
# below for the pooled FL/FV summary fits (specific tension N/cm^2,
# mass-specific peak power W/kg -- see muscle_geometry.R).
specimen_geom <- .read_file_level_geometry(manifest$fullpath[1L])
# Muscle lever arm (r_m, meters) -- SAME formula compute_predicted_strain()
# uses (half body width minus resolved muscle depth), needed alongside
# specimen_geom$muscle to convert the pooled fits' torque (N*m) to force
# (N) via add_specific_properties_to_fit().
specimen_depth <- resolve_muscle_depth_mm(specimen_geom$muscle_depth_mm_raw)
specimen_r_m <- (specimen_geom$local_body_width_mm / 2 - specimen_depth$depth_mm) / 1000.0
specimen_dclamp_m <- specimen_geom$dclamp_mm / 1000.0
cli::cli_alert_info(
  "Specimen geometry: body {specimen_geom$local_body_width_mm}x{specimen_geom$local_body_height_mm} mm (WxH oval), test section {specimen_geom$dclamp_mm} mm, density {specimen_geom$density_g_per_mm3} g/mm^3 -> est. muscle mass {signif(specimen_geom$muscle$muscle_mass_g,3)} g, muscle CSA {signif(specimen_geom$muscle$csa_muscle_cm2,3)} cm^2 ({specimen_geom$muscle$muscle_fraction*100}% of total {signif(specimen_geom$muscle$volume_total_cm3,3)} cm^3 / {signif(specimen_geom$muscle$csa_total_cm2,3)} cm^2)"
)

if (length(isometric_steps_all) > 0L) {
  iso_steps <- dplyr::bind_rows(isometric_steps_all)
  # No model fit is computed for pooled FL by default (PI direction,
  # 2026-07-16 -- see fit_fv_fl.R module header / 03_analyze.R
  # analyze_isometric() docstring): the summary plot's connect-the-mean
  # line (build_summary_plot_isometric(), plot_summary_profiles.R) is the
  # sole trend representation, so there is nothing to fit/report here.
  cli::cli_alert_info(
    "Isometric pooled FL: n = {nrow(iso_steps)} points pooled across {dplyr::n_distinct(iso_steps$trial_id)} trial(s) -- no model fit (connect-the-mean default, see summary plot)"
  )
  p_iso <- build_summary_plot_isometric(iso_steps)
  ggplot2::ggsave(file.path(SUMMARY_PLOT_DIR, "summary_isometric_FL.png"), p_iso, width = 9, height = 6, dpi = 150)

  # ADDITIONAL interpolated-baseline FL summary (2026-07-16, see
  # 03_analyze.R module header / build_segmented_step_summary()'s
  # passive_force_Nm_interp comment) -- saved as a SEPARATE file, original
  # summary_isometric_FL.png above is untouched. Reuses
  # build_summary_plot_isometric() unmodified by substituting
  # muscle_force_Nm with muscle_force_Nm_interp before calling it. Same
  # no-model-fit default applies (see above).
  iso_steps_interp <- dplyr::mutate(iso_steps, muscle_force_Nm = .data$muscle_force_Nm_interp)
  p_iso_interp <- build_summary_plot_isometric(iso_steps_interp,
                                                title = "Isometric summary: Force-Length [interpolated baseline]")
  ggplot2::ggsave(file.path(SUMMARY_PLOT_DIR, "summary_isometric_FL_interpbaseline.png"), p_iso_interp,
                  width = 9, height = 6, dpi = 150)
}

if (length(isovelocity_steps_all) > 0L) {
  isv_steps <- dplyr::bind_rows(isovelocity_steps_all)
  isv_fits <- purrr::map(c(left = "left", right = "right"), function(sd) {
    ss <- dplyr::filter(isv_steps, muscle_side == sd, is.finite(muscle_force_Nm))
    f0_row <- dplyr::filter(ss, abs(shortening_value) < 1e-6)
    f0_iso <- if (nrow(f0_row) > 0L) mean(f0_row$muscle_force_Nm, na.rm = TRUE) else NA_real_
    ss_conc <- dplyr::filter(ss, contraction_mode %in% c("concentric", "isometric_zero"))
    f <- fit_force_velocity_curve(ss_conc$shortening_strain_pct, ss_conc$muscle_force_Nm, side_label = sd, f0_isometric = f0_iso)
    add_specific_properties_to_fit(f, specimen_r_m, specimen_dclamp_m, specimen_geom$muscle, kind = "FV")
  })
  for (sd in names(isv_fits)) {
    f <- isv_fits[[sd]]
    cli::cli_alert(if (identical(f$status, "ok"))
      "Isovelocity pooled FV fit [{sd}]: OK (Hill model fit; Vmax={round(f$Vmax,2)}%/s [EXTRAPOLATED, not observed], Ppeak={signif(f$peak_power,3)}{if (!is.null(f$mass_specific_peak_power_Wkg)) sprintf(' [%.3g N/cm^2 specific tension, %.3g W/kg peak -- Coughlin scup red muscle ~114-152 W/kg]', f$specific_tension_Ncm2, f$mass_specific_peak_power_Wkg) else ''})"
      else "Isovelocity pooled FV fit [{sd}]: FAILED -- {f$reason}")
  }
  p_isv <- build_summary_plot_isovelocity(isv_steps, isv_fits)
  # Wider than the other summary plots (width=11 vs. 9) -- its subtitle
  # combines the descriptive connect-the-mean note with the full per-side
  # Hill fit annotation, so it wraps to more lines (see .wrap_subtitle(),
  # plot_summary_profiles.R) and needs the extra width to stay legible.
  ggplot2::ggsave(file.path(SUMMARY_PLOT_DIR, "summary_isovelocity_FV.png"), p_isv, width = 11, height = 6.5, dpi = 150)
}

if (length(dynamic_cycles_all) > 0L) {
  dyn_cycles <- dplyr::bind_rows(dynamic_cycles_all)
  p_dyn <- build_summary_plot_dynamic(dyn_cycles)
  ggplot2::ggsave(file.path(SUMMARY_PLOT_DIR, "summary_dynamic_power.png"), p_dyn, width = 9, height = 6, dpi = 150)

  # ADDITIONAL mass-specific power summary (PI-directed, 2026-07-16, see
  # compute_muscle_mass_and_csa()/muscle_geometry.R + .attach_dynamic_muscle_force
  # muscle_mass.kg wiring above) -- original summary_dynamic_power.png (raw
  # Watts) is untouched; this is a separately-saved W/kg version with a
  # Coughlin et al. 1996 reference band for weak/strong-twitch context.
  if (any(is.finite(dyn_cycles$avg_power.Wkg))) {
    p_dyn_msp <- build_summary_plot_dynamic_massspecific(dyn_cycles)
    if (!is.null(p_dyn_msp)) {
      ggplot2::ggsave(file.path(SUMMARY_PLOT_DIR, "summary_dynamic_power_massspecific.png"), p_dyn_msp,
                      width = 9, height = 6, dpi = 150)
      lims <- coughlin_steady_state_power_limits()
      n_exceed <- sum(dyn_cycles$exceeds_coughlin_hi, na.rm = TRUE)
      n_below  <- sum(dyn_cycles$below_coughlin_lo, na.rm = TRUE)
      # NOTE: this compares PEAK instantaneous mass-specific power (Coughlin's
      # own metric) against the reference band -- NOT the plotted MEAN
      # (cycle-averaged) power above, which is typically much smaller/can be
      # near zero or negative (a work loop's torque and angular velocity are
      # in opposite phase for part of each cycle). Peak-vs-mean divergence is
      # expected for oscillatory power, not a contradiction between the two.
      cli::cli_alert_info(
        "Dynamic mass-specific PEAK instantaneous power (est. muscle mass {signif(specimen_geom$muscle$muscle_mass_g,3)} g): {n_exceed}/{nrow(dyn_cycles)} cycles ABOVE Coughlin scup red-muscle steady-state band ({lims$lo}-{lims$hi} W/kg), {n_below}/{nrow(dyn_cycles)} BELOW -- NOT the same as the plotted mean/cycle-averaged power"
      )
    }
  } else {
    cli::cli_alert_warning("summary_dynamic_power_massspecific.png skipped -- avg_power.Wkg all NA (muscle mass estimate unavailable, check specimen geometry attrs)")
  }
}

if (length(freqsweep_cycles_all) > 0L) {
  fs_cycles <- dplyr::bind_rows(freqsweep_cycles_all)
  p_fs <- build_summary_plot_frequency_sweep(fs_cycles)
  ggplot2::ggsave(file.path(SUMMARY_PLOT_DIR, "summary_frequency_sweep_stiffness_damping.png"), p_fs, width = 9, height = 8, dpi = 150)
}

# Diagnostic (not an originally-required deliverable): fatigue check --
# muscle force vs. stimulation order within block, isometric + isovelocity
# pooled, to test whether the monotonic FL/FV shapes could be a fatigue
# artifact rather than a true length/velocity property.
segmented_all <- dplyr::bind_rows(
  if (length(isometric_steps_all) > 0L) dplyr::bind_rows(isometric_steps_all) |> dplyr::mutate(protocol_family = "isometric"),
  if (length(isovelocity_steps_all) > 0L) dplyr::bind_rows(isovelocity_steps_all) |> dplyr::mutate(protocol_family = "isovelocity")
)
if (nrow(segmented_all) > 0L) {
  for (fam in unique(segmented_all$protocol_family)) {
    p_fat <- build_fatigue_check_plot(
      dplyr::filter(segmented_all, protocol_family == fam),
      title = sprintf("Fatigue check (%s): muscle force vs. stimulation order", fam)
    )
    ggplot2::ggsave(file.path(SUMMARY_PLOT_DIR, sprintf("fatigue_check_%s.png", fam)), p_fat, width = 10, height = 6, dpi = 150)
  }
  # ADDITIONAL interpolated-baseline fatigue check (isometric only, same
  # substitution pattern as the FL summary above) -- original
  # fatigue_check_isometric.png is untouched.
  if ("isometric" %in% segmented_all$protocol_family) {
    p_fat_interp <- build_fatigue_check_plot(
      dplyr::filter(segmented_all, protocol_family == "isometric") |>
        dplyr::mutate(muscle_force_Nm = .data$muscle_force_Nm_interp),
      title = "Fatigue check (isometric) [interpolated baseline]: muscle force vs. stimulation order"
    )
    ggplot2::ggsave(file.path(SUMMARY_PLOT_DIR, "fatigue_check_isometric_interpbaseline.png"), p_fat_interp,
                    width = 10, height = 6, dpi = 150)
  }
}

# Measured (E6 encoder) vs. predicted (commanded-angle) strain -- two levels:
#  - "_commanded": CONTINUOUS per-sample strain_pct (commanded) vs.
#    strain_measured_pct (encoder), actively-stimulated samples only, one
#    pooled panel per protocol_family (dynamic/isovelocity/isometric).
#    frequency_sweep excluded -- passive-only, no actively-stimulated rows.
#  - "_measured": STEP-LEVEL shortening_strain_pct (the FL curve's OWN
#    x-axis, from the step's commanded operating_point) vs. the new
#    strain_measured_step_pct (encoder, same sign fold) -- ISOMETRIC ONLY.
#    Isovelocity's shortening_strain_pct is a strain-RATE (%/s, from a
#    commanded velocity), not a position, so it is not unit-comparable to
#    positional encoder strain here (see build_step_strain_validation_plot()
#    docstring); dynamic/frequency_sweep have no discrete-step
#    operating_point structure. Uses segmented_all, built above for the
#    fatigue-check plots.
if (length(strain_check_all) > 0L) {
  strain_check_df <- dplyr::bind_rows(strain_check_all)
  p_strain <- build_measured_vs_predicted_strain_plot(strain_check_df)
  ggplot2::ggsave(file.path(SUMMARY_PLOT_DIR, "strain_validation_commanded.png"), p_strain, width = 12, height = 5, dpi = 150)
  cli::cli_alert_info("strain_validation_commanded.png: {nrow(strain_check_df)} actively-stimulated samples pooled across {dplyr::n_distinct(strain_check_df$trial_id)} trial(s)")
} else {
  cli::cli_alert_warning("No actively-stimulated samples found for strain_validation_commanded.png -- skipped")
}

if (nrow(segmented_all) > 0L && "strain_measured_step_pct" %in% names(segmented_all)) {
  step_strain_df <- dplyr::filter(segmented_all, is.finite(.data$shortening_strain_pct), is.finite(.data$strain_measured_step_pct))
  if (nrow(step_strain_df) > 0L) {
    p_step_strain <- build_step_strain_validation_plot(step_strain_df)
    ggplot2::ggsave(file.path(SUMMARY_PLOT_DIR, "strain_validation_measured.png"), p_step_strain, width = 10, height = 5, dpi = 150)
    cli::cli_alert_info("strain_validation_measured.png: {nrow(step_strain_df)} step(s) pooled across {dplyr::n_distinct(step_strain_df$trial_id)} trial(s)")
  } else {
    cli::cli_alert_warning("No finite step-level strain pairs found for strain_validation_measured.png -- skipped")
  }
} else {
  cli::cli_alert_warning("No segmented (isometric/isovelocity) steps found for strain_validation_measured.png -- skipped")
}

# PI-requested (2026-07-16): three additional validation figures, ALL
# samples (no active-stim restriction -- these are purely mechanical/
# geometric checks), one panel per protocol category where data exists:
#   1) angle_validation.png -- measured (E6 encoder) vs. commanded angle.
#   2/3) strain_validation_sono_vs_{encoder,commanded}.png -- measured
#      (sonomicrometry, RIGHT muscle -- the only wired sono channel in this
#      rig) vs. predicted (curvature geometry, right-folded) strain, from
#      each of the two candidate angle sources. See
#      plot_angle_sono_validation.R header for full scope/sign-convention
#      rationale.
if (length(angle_check_all) > 0L) {
  angle_check_df <- dplyr::bind_rows(angle_check_all)
  p_angle <- build_angle_validation_plot(angle_check_df)
  ggplot2::ggsave(file.path(SUMMARY_PLOT_DIR, "angle_validation.png"), p_angle, width = 12, height = 8, dpi = 150)
  cli::cli_alert_info("angle_validation.png: {nrow(angle_check_df)} samples pooled across {dplyr::n_distinct(angle_check_df$trial_id)} trial(s)")
} else {
  cli::cli_alert_warning("No angle data found for angle_validation.png -- skipped")
}

if (length(sono_strain_check_all) > 0L) {
  sono_check_df <- dplyr::bind_rows(sono_strain_check_all)
  if (all(is.na(sono_check_df$strain_sono_pct))) {
    cli::cli_alert_warning("sono_right channel/calibration unavailable in every trial -- sono validation figures skipped")
  } else {
    p_sono_enc <- build_sono_strain_validation_plot(sono_check_df, "strain_pred_encoder_right_pct", "encoder")
    ggplot2::ggsave(file.path(SUMMARY_PLOT_DIR, "strain_validation_sono_vs_encoder.png"), p_sono_enc, width = 12, height = 8, dpi = 150)
    p_sono_cmd <- build_sono_strain_validation_plot(sono_check_df, "strain_pred_commanded_right_pct", "commanded")
    ggplot2::ggsave(file.path(SUMMARY_PLOT_DIR, "strain_validation_sono_vs_commanded.png"), p_sono_cmd, width = 12, height = 8, dpi = 150)
    cli::cli_alert_info(
      "strain_validation_sono_vs_{{encoder,commanded}}.png: {sum(is.finite(sono_check_df$strain_sono_pct))} samples with valid sono data, pooled across {dplyr::n_distinct(dplyr::filter(sono_check_df, is.finite(strain_sono_pct))$trial_id)} trial(s)"
    )
  }
} else {
  cli::cli_alert_warning("No sono-strain data found -- strain_validation_sono_vs_{{encoder,commanded}}.png skipped")
}

# Muscle force vs. time, pooled across trials, one figure per protocol
# category -- x-axis truncated to [pre-stim baseline, longest stim duration +
# relaxation window] (see plot_force_vs_time.R header for rationale).
cli::cli_h1("Building force-vs-time plots")
for (fam in names(force_ts_all)) {
  ts_list <- force_ts_all[[fam]]
  if (length(ts_list) == 0L) {
    cli::cli_alert_warning("No force-vs-time data for category '{fam}' -- skipped")
    next
  }
  ts_df <- dplyr::bind_rows(ts_list)

  # Isovelocity: color by muscle_side (same left/right convention as the FV
  # summary plot) and facet by contraction_mode -- the earlier single
  # trial_id-colored pooled line mixed concentric (force typically LOWER
  # than the velocity-matched no-stim baseline) and eccentric (typically
  # higher) steps from BOTH sides into one trace, which visually looked
  # "unadjusted" even though force_sign was already applied per-step. This
  # view lets left vs. right be compared directly WITHIN the same
  # contraction mode, where side-adjustment should make them track together.
  facet_var <- if (fam == "isovelocity") "contraction_mode" else NULL
  color_var <- if (fam == "isovelocity") "muscle_side" else "trial_id"

  p_ts <- build_force_vs_time_plot(ts_df, title = sprintf("Muscle force vs. time (%s, pooled across trials)", fam),
                                    facet_var = facet_var, color_var = color_var)
  ggplot2::ggsave(file.path(SUMMARY_PLOT_DIR, sprintf("force_vs_time_%s.png", fam)), p_ts,
                  width = if (!is.null(facet_var)) 13 else 10, height = 6, dpi = 150)
  cli::cli_alert_info("force_vs_time_{fam}.png: {dplyr::n_distinct(ts_df$unit_id)} stim event(s)/step(s) across {dplyr::n_distinct(ts_df$trial_id)} trial(s)")

  # Dynamic-only: additionally break the pooled plot into one figure PER
  # commanded variable (frequency, amplitude, duty, phase), each faceted by
  # that variable's distinct levels -- dynamic trials mix multiple commanded
  # conditions (see detect_trial_type()'s freq_amp_combo/constant_multi_param
  # layouts), and a single pooled line averages across all of them together.
  if (fam == "dynamic") {
    dyn_breakdowns <- list(
      frequency = list(col = "freq_hz", label = "frequency (Hz)"),
      amplitude = list(col = "amp_deg", label = "amplitude (deg)"),
      phase     = list(col = "phase",   label = "phase"),
      duty      = list(col = "duty",    label = "duty")
    )
    for (bname in names(dyn_breakdowns)) {
      b <- dyn_breakdowns[[bname]]
      if (!b$col %in% names(ts_df) || dplyr::n_distinct(ts_df[[b$col]], na.rm = TRUE) < 1L) {
        cli::cli_alert_warning("force_vs_time_dynamic_by_{bname}.png: no finite '{b$col}' values -- skipped")
        next
      }
      p_dyn <- build_force_vs_time_plot(
        ts_df, title = sprintf("Muscle force vs. time (dynamic, grouped by %s)", b$label),
        facet_var = b$col, color_var = "muscle_side"
      )
      ggplot2::ggsave(file.path(SUMMARY_PLOT_DIR, sprintf("force_vs_time_dynamic_by_%s.png", bname)), p_dyn,
                      width = 13, height = 8, dpi = 150)
      cli::cli_alert_info("force_vs_time_dynamic_by_{bname}.png: {dplyr::n_distinct(ts_df[[b$col]], na.rm = TRUE)} distinct {bname} level(s)")
    }
  }
}

cli::cli_h1("Pipeline complete")
print(trial_report_tbl |> dplyr::select(trial_id, category, status, fit_status, note))
