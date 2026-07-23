# plot_sono_strain_validation_pooled.R
# PI-requested pooled companion to the per-specimen <specimen>_strainValidSonoEnc.png
# figures (run_fv_fl_power_pipeline.R + plot_angle_sono_validation.R): ONE
# figure PER PROTOCOL FAMILY (isometric, isovelocity, dynamic, frequency_sweep),
# each pooling measured-(sonomicrometry)-vs-predicted-(encoder-angle) RIGHT-
# muscle strain across ALL THREE specimens (bass16/17/18), split into two
# panels: ACTIVE vs PASSIVE.
#
# Definition of active/passive (PI-directed, 2026-07-22): the sono-
# instrumented channel is the RIGHT muscle only (see plot_angle_sono_
# validation.R header) -- "active" therefore means the RIGHT muscle is the
# one being commanded to contract (right stim); "passive" is everything else
# for the right muscle's own mechanical state, i.e. no stim OR left stim
# (during a left-stim burst the right muscle itself is not being driven,
# even though the specimen is bending).
#   active_passive = "active (right stim)"       if stim_state == "right stim"
#                   = "passive (no/left stim)"    otherwise (no stim, left stim)
#
# Data recipe per trial is IDENTICAL to run_fv_fl_power_pipeline.R's
# sono_check construction for each category (verified against that file,
# 2026-07-22): load_bender_flat() -> attach_predicted_strain() ->
# attach_measured_strain() -> attach_sono_strain() (the vetted sign
# convention -- force_sign_right folds ONLY the predicted side, never the
# sono side; see attach_sono_strain() docstring). The pipeline additionally
# routes isometric/isovelocity through analyze_isometric()/analyze_isovelocity()
# first, but that step-segmentation is NOT needed here: attach_predicted_strain()'s
# active_mask argument only ever feeds the (unused-by-this-figure)
# strain_active_pct/strain_passive_pct columns, and attach_measured_strain()
# depends only on enc.deg/dclamp.m/muscle_strain_r_m, none of which require
# step segmentation. Confirmed by reading both functions (muscle_geometry.R,
# plot_strain_validation.R) before writing this script -- this is a
# deliberate simplification, not an oversight.
#
# Sono oversampling fix (2026-07-22): the DS3's true internal update rate is
# ~241-247 Hz (daq_sono_internal_sample_rate_hertz) but its analog output is
# digitized on the 1 kHz AI clock, so every real sono reading was previously
# counted (and plotted) ~4x. .process_trial() now decimates to one sample
# per true DS3 update interval (round(1000 Hz / ds3_rate_hz) -- NOT an
# exact-value-repeat dedup; checked empirically that the raw signal doesn't
# hold a flat step between updates, see .process_trial()'s own comment),
# applied AFTER stim-state labeling so no stim pulse is missed.
#
# Output:
#   figs_summary/pooled_strainValidSonoEnc_<protocol_family>.png -- one file
#     per family that has any valid sono data across the corpus (ALL trials).
#   figs_summary/pooled_strainValidSonoEnc_dynamic_later.png -- dynamic only,
#     restricted to each specimen's "later (stable)" trials (excludes early
#     tissue-preconditioning trials, dynamic_trial_precondition.R).
#   figs_summary/pooled_strainValidSonoEnc_allProtocols_later.png -- ONE
#     multi-panel figure, ALL FOUR protocol families x active/passive,
#     restricted to "later (stable)" trials for every family (the session-
#     chronology cutoff applied protocol-agnostically via
#     classify_session_precondition() -- see that function's docstring for
#     why this is a deliberate extension, not dynamic-only).
#
# Run with:  Rscript R/plot_sono_strain_validation_pooled.R

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2); library(cli); library(rhdf5); library(purrr)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")
src("00_load_bender_flat.R")
src("01_calibrate.R")
src("muscle_geometry.R")
src("plot_strain_validation.R")     # attach_measured_strain
src("plot_force_vs_time.R")         # .detect_stim_events -- REQUIRED by .stim_window_state_label();
                                     # without it, .stim_window_state_label()'s internal tryCatch
                                     # silently swallows the "function not found" error and labels
                                     # EVERY sample "no stim" with no visible error (caught 2026-07-22
                                     # via a manual stim-state audit -- see script run history).
src("plot_angle_sono_validation.R") # attach_sono_strain, .stim_window_state_label, STIM_STATE_LEVELS
src("parse_trial_filename.R")

OUT_DIR      <- FIGS_SUMMARY_DIR
DATA_OUT_DIR <- file.path(.crittergripper_root(), "02_processed", "data_processed")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SPECIMEN_SUBFOLDERS <- c(
  bass16 = BASS16_RAW_SUBFOLDER,
  bass17 = BASS17_RAW_SUBFOLDER,
  bass18 = BASS18_RAW_SUBFOLDER
)
PROTOCOL_FAMILY_LEVELS <- c("isometric", "isovelocity", "dynamic", "frequency_sweep")
SPECIMEN_COLORS <- c(bass16 = "#1b9e77", bass17 = "#d95f02", bass18 = "#7570b3")
ACTIVE_PASSIVE_LEVELS <- c("passive (no/left stim)", "active (right stim)")
STEP_ACTIVITY_LEVELS  <- c("purely passive (no stim in step/cycle)", "contains stimulation")

# =============================================================================
# Geometry reader (same pattern as diag_sono_*.R scripts)
# =============================================================================
.read_geom_and_lidx <- function(h5path) {
  h5 <- H5Fopen(h5path, "H5F_ACC_RDONLY")
  on.exit(try(H5Fclose(h5), silent = TRUE), add = TRUE)
  m_attrs <- tryCatch(h5readAttributes(h5, "/metadata"), error = function(e) list())
  dbl1 <- function(key, default = NA_real_) {
    v <- suppressWarnings(as.numeric(m_attrs[[key]][1L]))
    if (length(v) == 0L || is.na(v)) default else v
  }
  list(
    width_mm       = dbl1("measurement_specimen_local_body_width_millimeter"),
    depth_mm       = dbl1("measurement_target_muscle_depth_millimeter"),
    lidx_pos_motor = dbl1("daq_specimen_lateral_index_on_positive_motor_side"),
    lidx_right     = dbl1("daq_specimen_side_index_right"),
    ds3_rate_hz    = dbl1("daq_sono_internal_sample_rate_hertz", default = 247.0)
  )
}

AI_SAMPLE_RATE_HZ <- 1000.0  # daq_ai_sample_rate_hz, all rig configs (see .cursorrules)

# =============================================================================
# Per-trial sono-check extraction (mirrors run_fv_fl_power_pipeline.R's
# sono_check construction for every category -- see header note above)
# =============================================================================
.process_trial <- function(specimen, h5path, category) {
  td <- tryCatch(load_bender_flat(h5path, do_filter = TRUE, loadtorques = "x"), error = function(e) NULL)
  if (is.null(td) || nrow(td) == 0L) return(tibble())

  geom <- .read_geom_and_lidx(h5path)
  if (!is.finite(geom$width_mm) || !is.finite(geom$lidx_right) || !is.finite(geom$lidx_pos_motor)) return(tibble())

  # Matches run_fv_fl_power_pipeline.R exactly: frequency_sweep is passive-
  # only (active_mask always FALSE); isometric/isovelocity/dynamic use the
  # real stim column. Harmless either way for THIS figure (see header note),
  # kept faithful for anyone reusing strain_active_pct/strain_passive_pct.
  active_mask <- if (category == "frequency_sweep") rep(FALSE, nrow(td)) else (as.character(td$stim) != "0")

  td <- attach_predicted_strain(td, local_body_width_mm = geom$width_mm, measured_muscle_depth_mm = geom$depth_mm,
                                 active_mask = active_mask)
  td <- attach_measured_strain(td)
  td <- attach_sono_strain(td, h5path, geom$lidx_right, geom$lidx_pos_motor)

  # Step/cycle-level activity flag (PI-directed, 2026-07-23): "step_number"
  # exists only for segmented_finite (isometric/isovelocity, one designed
  # activation window per step); single_finite (dynamic/frequency_sweep) has
  # no step_number column at all, but does have "cycle_index". Either way,
  # this identifies the natural experimental unit (one designed step, or one
  # bending cycle) so we can ask "did stim fire ANYWHERE in this whole
  # unit", independent of the WINDOWED per-sample stim_state below -- a
  # step/cycle that had a real stim burst but has since relaxed can still
  # carry residual force/strain into samples the windowed label would call
  # "no stim", which is exactly the contamination this flag is meant to
  # separate out. Computed on the FULL-RATE td (not the decimated `out`)
  # so a short burst can't be missed by decimation.
  step_or_cycle_id <- if (category %in% c("isometric", "isovelocity")) td$step_number else td$cycle_index
  raw_stim_active  <- as.character(td$stim) != "0"
  step_has_stim    <- as.logical(stats::ave(as.integer(raw_stim_active), step_or_cycle_id, FUN = max))

  stim_state <- as.character(.stim_window_state_label(td))
  out <- tibble(
    specimen = specimen, trial_id = tools::file_path_sans_ext(basename(h5path)), protocol_family = category,
    strain_pred_encoder_right_pct = td$strain_pred_encoder_right_pct,
    strain_sono_pct               = td$strain_sono_pct,
    step_has_stim                 = step_has_stim,
    stim_state                    = stim_state
  )

  # De-duplicate sono oversampling (2026-07-22 fix): the sono channel is
  # digitized on the 1 kHz AI clock, but the DS3's own true internal update
  # rate is ~241-247 Hz (daq_sono_internal_sample_rate_hertz), so every real
  # sono reading was previously being counted (and plotted) ~4x. Checked
  # empirically first (2026-07-22): the raw calibrated mm signal does NOT
  # show an exact-value zero-order-hold staircase sample-to-sample (only
  # ~2% exact-repeat diffs in a spot check) -- consecutive AI samples differ
  # by small, continuously-varying amounts, consistent with the DS3's analog
  # output stage smoothing/reconstructing its own ~241-247 Hz updates before
  # the 1 kHz ADC digitizes it, rather than holding a flat step. Exact-value
  # dedup is therefore the WRONG tool here (it would keep ~98% of "duplicate"
  # oversampled rows). Instead, DECIMATE to one sample per true DS3 update
  # interval (round(AI_rate / ds3_rate)), which is robust regardless of
  # whether the true hardware behavior is a hard step or a smoothed ramp --
  # either way, only the DS3's own true update rate carries independent
  # information. Applied to `out` (not the full-rate td used for stim-state
  # labeling above), so .detect_stim_events() still sees every row and
  # cannot miss a single-sample stim pulse.
  decim_n <- max(1L, round(AI_SAMPLE_RATE_HZ / geom$ds3_rate_hz))
  n_before <- nrow(out)
  out <- out[seq(1L, n_before, by = decim_n), , drop = FALSE]
  cli::cli_alert_success("{basename(h5path)} ({category}): {sum(is.finite(out$strain_pred_encoder_right_pct) & is.finite(out$strain_sono_pct))} valid sono samples ({n_before} -> {nrow(out)}, 1-in-{decim_n} decimated to the true ~{round(geom$ds3_rate_hz)} Hz DS3 update rate)")
  out
}

# =============================================================================
# Batch: discover + classify every trial per specimen (same parser the main
# pipeline uses), process, pool
# =============================================================================
cli::cli_h1("Pooling sono-strain validation across bass16/17/18, all protocol families")

all_rows <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  d <- raw_source_dir(subfolder)
  manifest <- tryCatch(parse_trial_directory(d), error = function(e) { cli::cli_alert_danger("{specimen}: {conditionMessage(e)}"); tibble() })
  if (nrow(manifest) == 0L) return(tibble())
  manifest <- dplyr::filter(manifest, .data$protocol %in% PROTOCOL_FAMILY_LEVELS)
  tibble(specimen = specimen, h5path = manifest$fullpath, category = manifest$protocol)
})

pooled <- purrr::pmap_dfr(all_rows, function(specimen, h5path, category) {
  tryCatch(.process_trial(specimen, h5path, category),
           error = function(e) { cli::cli_alert_danger("{basename(h5path)}: {conditionMessage(e)}"); tibble() })
})

src("dynamic_trial_precondition.R")

pooled <- pooled |>
  dplyr::filter(is.finite(.data$strain_pred_encoder_right_pct), is.finite(.data$strain_sono_pct)) |>
  dplyr::mutate(
    protocol_family = factor(.data$protocol_family, levels = PROTOCOL_FAMILY_LEVELS),
    active_passive = factor(
      dplyr::if_else(.data$stim_state == "right stim", "active (right stim)", "passive (no/left stim)"),
      levels = ACTIVE_PASSIVE_LEVELS
    ),
    trial_num    = extract_bender_trial_num(.data$trial_id),
    precondition = classify_session_precondition(.data$specimen, .data$trial_num),
    step_activity = factor(
      dplyr::if_else(.data$step_has_stim, "contains stimulation", "purely passive (no stim in step/cycle)"),
      levels = STEP_ACTIVITY_LEVELS
    )
  )

write.csv(pooled, file.path(DATA_OUT_DIR, "sono_strain_validation_pooled_samples.csv"), row.names = FALSE)
cli::cli_alert_success("Saved {nrow(pooled)} pooled samples -> data_processed/sono_strain_validation_pooled_samples.csv")

# =============================================================================
# One figure per protocol family: 2 panels (passive / active), pooled
# specimens colored, r/RMSE per panel
# =============================================================================
.build_family_plot <- function(fam, df = NULL, title_suffix = "") {
  if (is.null(df)) df <- dplyr::filter(pooled, .data$protocol_family == fam)
  if (nrow(df) < 20L) return(NULL)

  lims <- range(c(df$strain_pred_encoder_right_pct, df$strain_sono_pct), na.rm = TRUE)
  ref_df <- tibble(x = lims, y = lims)

  r2_labels <- df |>
    dplyr::group_by(.data$active_passive) |>
    dplyr::summarise(
      r    = suppressWarnings(cor(.data$strain_pred_encoder_right_pct, .data$strain_sono_pct, use = "complete.obs")),
      rmse = sqrt(mean((.data$strain_pred_encoder_right_pct - .data$strain_sono_pct)^2, na.rm = TRUE)),
      n    = dplyr::n(), .groups = "drop"
    ) |>
    dplyr::mutate(label = sprintf("r=%.3f, RMSE=%.2f, n=%s", .data$r, .data$rmse, format(.data$n, big.mark = ",")))

  n_specimens <- dplyr::n_distinct(df$specimen)

  ggplot(df, aes(x = .data$strain_pred_encoder_right_pct, y = .data$strain_sono_pct, color = .data$specimen)) +
    geom_point(shape = 1, size = 0.8, alpha = 0.2, stroke = 0.3) +
    geom_line(data = ref_df, aes(x = .data$x, y = .data$y), inherit.aes = FALSE,
              linetype = "dashed", color = "black", linewidth = 0.6) +
    geom_text(data = r2_labels, aes(x = -Inf, y = Inf, label = .data$label),
              inherit.aes = FALSE, hjust = -0.05, vjust = 1.3, size = 3.4, color = "black") +
    facet_wrap(~active_passive) +
    scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen", drop = FALSE) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
    coord_equal(xlim = lims, ylim = lims) +
    labs(title = sprintf("Pooled measured (sonomicrometry, right muscle) vs. predicted (encoder angle) strain -- %s%s", fam, title_suffix),
         subtitle = sprintf("Pooled across %d specimen(s) (bass16/17/18); RIGHT muscle only; dashed = 1:1 reference.\nActive = right stim (the sono-instrumented muscle itself contracting); passive = no stim OR left stim.", n_specimens),
         x = "Predicted strain (%, from encoder angle, right-folded)",
         y = "Measured strain (%, sonomicrometry, right muscle)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
}

for (fam in PROTOCOL_FAMILY_LEVELS) {
  p <- .build_family_plot(fam)
  if (is.null(p)) {
    cli::cli_alert_warning("{fam}: insufficient pooled sono data -- skipped")
    next
  }
  fout <- file.path(OUT_DIR, sprintf("pooled_strainValidSonoEnc_%s.png", fam))
  ggplot2::ggsave(fout, p, width = 11, height = 6, dpi = 150)
  cli::cli_alert_success("Saved {fout}")
}

# =============================================================================
# Dynamic, "later (stable)" trials only -- excludes each specimen's early
# tissue-preconditioning trials using the hard cutoff decided 2026-07-22
# (dynamic_trial_precondition.R). Purpose: confirm the pooled dynamic-active
# panel above recovers isometric/isovelocity-like tightness once each
# specimen's own preconditioning trials are excluded, rather than showing a
# dynamic-protocol-wide artifact.
# =============================================================================
dynamic_later <- pooled |>
  dplyr::filter(.data$protocol_family == "dynamic", .data$precondition == "later (stable)")

write.csv(dynamic_later, file.path(DATA_OUT_DIR, "sono_strain_validation_pooled_dynamic_later_samples.csv"), row.names = FALSE)
cli::cli_alert_success("Saved {nrow(dynamic_later)} 'later (stable)' dynamic samples -> data_processed/sono_strain_validation_pooled_dynamic_later_samples.csv")

p_dyn_later <- .build_family_plot("dynamic", df = dynamic_later,
                                   title_suffix = " [LATER (stable) trials only -- preconditioning trials excluded]")
if (is.null(p_dyn_later)) {
  cli::cli_alert_warning("dynamic (later-only): insufficient pooled sono data -- skipped")
} else {
  fout_later <- file.path(OUT_DIR, "pooled_strainValidSonoEnc_dynamic_later.png")
  ggplot2::ggsave(fout_later, p_dyn_later, width = 11, height = 6, dpi = 150)
  cli::cli_alert_success("Saved {fout_later}")
}

# =============================================================================
# ALL protocol families in ONE multi-panel figure, "later (stable)" trials
# only -- style of the original per-specimen <specimen>_strainValidSonoEnc.png
# (facet by protocol_family), extended with an active/passive facet dimension
# and pooled across all 3 specimens (PI-directed, 2026-07-22).
#
# The "later (stable)" filter applies to EVERY protocol family here (not just
# dynamic) via classify_session_precondition() -- see that function's
# docstring for why the dynamic-derived cutoff is applied protocol-agnostic-
# ally (a session-chronology tissue-settling effect, not dynamic-specific).
# For isometric/isovelocity/frequency_sweep, this typically removes little or
# nothing (those protocols mostly run AFTER a specimen's first few dynamic
# trials in every session to date), but is applied uniformly on principle
# rather than special-cased per family.
# =============================================================================
cli::cli_h1("All-protocol-family combined figure, 'later (stable)' trials only")

all_later <- dplyr::filter(pooled, .data$precondition == "later (stable)")
write.csv(all_later, file.path(DATA_OUT_DIR, "sono_strain_validation_pooled_allprotocols_later_samples.csv"), row.names = FALSE)
cli::cli_alert_success("Saved {nrow(all_later)} 'later (stable)' samples across all protocol families -> data_processed/sono_strain_validation_pooled_allprotocols_later_samples.csv")

n_excluded_by_family <- pooled |>
  dplyr::group_by(.data$protocol_family) |>
  dplyr::summarise(n_total = dplyr::n(), n_early = sum(.data$precondition == "early (preconditioning)", na.rm = TRUE), .groups = "drop")
print(n_excluded_by_family)

#' @param col_var Column name (string) to facet the grid's COLUMNS on --
#'   "active_passive" (per-SAMPLE windowed stim state) by default, or
#'   "step_activity" (per-STEP/CYCLE any-stim-anywhere) for the alternate
#'   version below.
#' @param col_caption Short phrase describing col_var, substituted into the
#'   subtitle so the two versions of this figure don't share a misleading
#'   caption.
.build_allprotocols_plot <- function(df, col_var = "active_passive", col_caption = NULL) {
  r2_labels <- df |>
    dplyr::group_by(.data$protocol_family, .data[[col_var]]) |>
    dplyr::summarise(
      r    = suppressWarnings(cor(.data$strain_pred_encoder_right_pct, .data$strain_sono_pct, use = "complete.obs")),
      rmse = sqrt(mean((.data$strain_pred_encoder_right_pct - .data$strain_sono_pct)^2, na.rm = TRUE)),
      n    = dplyr::n(), .groups = "drop"
    ) |>
    dplyr::mutate(label = sprintf("r=%.3f, RMSE=%.2f, n=%s", .data$r, .data$rmse, format(.data$n, big.mark = ",")))

  n_specimens <- dplyr::n_distinct(df$specimen)
  col_caption <- col_caption %||% ""

  ggplot(df, aes(x = .data$strain_pred_encoder_right_pct, y = .data$strain_sono_pct, color = .data$specimen)) +
    geom_point(shape = 1, size = 0.6, alpha = 0.15, stroke = 0.25) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_text(data = r2_labels, aes(x = -Inf, y = Inf, label = .data$label),
              inherit.aes = FALSE, hjust = -0.05, vjust = 1.3, size = 2.9, color = "black") +
    facet_grid(stats::reformulate(col_var, "protocol_family"), scales = "free") +
    scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen", drop = FALSE) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
    labs(title = "Pooled measured (sonomicrometry, right muscle) vs. predicted (encoder angle) strain -- all protocol families",
         subtitle = sprintf("Pooled across %d specimen(s) (bass16/17/18); RIGHT muscle only; dashed = 1:1 reference; sono decimated to one point per true ~241-247 Hz DS3 update (no oversampling duplicates).\n\"LATER (stable)\" trials only -- each specimen's early tissue-preconditioning trials excluded (dynamic_trial_precondition.R), applied to every protocol family.%s", n_specimens, col_caption),
         x = "Predicted strain (%, from encoder angle, right-folded)",
         y = "Measured strain (%, sonomicrometry, right muscle)") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", strip.text = element_text(size = 9))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

p_all_later <- .build_allprotocols_plot(all_later, col_var = "active_passive")
fout_all_later <- file.path(OUT_DIR, "pooled_strainValidSonoEnc_allProtocols_later.png")
ggplot2::ggsave(fout_all_later, p_all_later, width = 11, height = 12, dpi = 150)
cli::cli_alert_success("Saved {fout_all_later}")

# =============================================================================
# ALTERNATE version of the figure above: columns split by STEP/CYCLE-level
# activity (did stim fire ANYWHERE in this designed step or bending cycle?)
# rather than the per-SAMPLE windowed stim_state (PI-directed, 2026-07-23).
#
# Motivation: the "passive (no/left stim)" bucket in the sample-level figure
# can still include samples that sit AFTER a right-stim burst within the same
# step/cycle (relaxation tail, incomplete mechanical recovery) -- i.e.
# residual force from stimulation earlier in that same unit, mislabeled as
# equivalent to a step/cycle that was never stimulated at all. This version
# instead asks, per whole step (isometric/isovelocity) or bending cycle
# (dynamic/frequency_sweep): "purely passive (no stim in step/cycle)" if
# stim NEVER fired anywhere in that unit, vs. "contains stimulation" if it
# fired at any point (any side) -- computed from the raw stim column on the
# full-rate signal (`.process_trial()`'s step_has_stim), not the windowed
# per-sample label, so it can't miss a brief burst to decimation or
# window-boundary effects.
# =============================================================================
cli::cli_h1("All-protocol-family combined figure, STEP/CYCLE-level activity split")

p_all_later_step <- .build_allprotocols_plot(
  all_later, col_var = "step_activity",
  col_caption = "\nColumns split by STEP/CYCLE-level activity (any stim anywhere in that step/cycle), not per-sample windowed stim state -- isolates residual post-stimulation force/strain from truly-never-stimulated units."
)
fout_all_later_step <- file.path(OUT_DIR, "pooled_strainValidSonoEnc_allProtocols_later_stepActivity.png")
ggplot2::ggsave(fout_all_later_step, p_all_later_step, width = 11, height = 12, dpi = 150)
cli::cli_alert_success("Saved {fout_all_later_step}")

cli::cli_alert_success("plot_sono_strain_validation_pooled.R complete")
