# diag_precondition_power_check.R
# Does the early-trial tissue-preconditioning effect (dynamic_trial_
# precondition.R) show up in muscle POWER OUTPUT too, for comparable bending
# (curvature.invm) and stimulation (freq.Hz, active cycles only)? PI-directed
# follow-up, 2026-07-22, after the sono-strain-offset cutoff was decided:
# "It would be interesting to see whether this translates to differences in
# muscle power output (for comparable bending and stimulation)."
#
# Two data sources:
#   (a) sono-strain offset per trial -- reused DIRECTLY from the already-
#       vetted pooled CSV (plot_sono_strain_validation_pooled.R's
#       sono_strain_validation_pooled_samples.csv), NOT recomputed here, so
#       this script cannot silently drift from the numbers that set the
#       cutoff in the first place.
#   (b) per-cycle muscle power for dynamic trials -- computed fresh via the
#       SAME calc_muscle_torque()/summarize_muscle_cycles() recipe
#       run_fv_fl_power_pipeline.R uses (that pipeline is plot-only, no CSV
#       export exists to reuse -- see compare_specimen_specific_properties.R
#       for the established duplication pattern this script follows).
#
# Run with:  Rscript R/diag_precondition_power_check.R
# Outputs -> figures: 02_processed/figs_diagnostic/ (FIGS_DIAGNOSTIC_DIR)
#            data:    02_processed/data_processed/

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(ggplot2); library(cli); library(rhdf5); library(readr)
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
src("plot_force_vs_time.R")  # .detect_stim_events(), RELAXATION_WINDOW_S
src("dynamic_trial_precondition.R")

OUT_DIR      <- FIGS_DIAGNOSTIC_DIR
DATA_OUT_DIR <- file.path(.crittergripper_root(), "02_processed", "data_processed")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SPECIMEN_SUBFOLDERS <- c(bass16 = BASS16_RAW_SUBFOLDER, bass17 = BASS17_RAW_SUBFOLDER, bass18 = BASS18_RAW_SUBFOLDER)
SPECIMEN_COLORS     <- c(bass16 = "#1b9e77", bass17 = "#d95f02", bass18 = "#7570b3")

# =============================================================================
# Part A: sono-strain offset per trial (reused, not recomputed)
# =============================================================================
cli::cli_h1("Part A: sono-strain offset by trial order (reusing pooled sono CSV)")

pooled_path <- file.path(DATA_OUT_DIR, "sono_strain_validation_pooled_samples.csv")
pooled <- readr::read_csv(pooled_path, show_col_types = FALSE)

offset_by_trial <- pooled |>
  dplyr::filter(.data$protocol_family == "dynamic", .data$active_passive == "active (right stim)",
                is.finite(.data$strain_pred_encoder_right_pct), is.finite(.data$strain_sono_pct)) |>
  dplyr::mutate(trial_num = extract_bender_trial_num(.data$trial_id)) |>
  dplyr::filter(is.finite(.data$trial_num)) |>
  dplyr::group_by(.data$specimen, .data$trial_id, .data$trial_num) |>
  dplyr::summarise(n = dplyr::n(),
                    offset_pct = mean(.data$strain_sono_pct) - mean(.data$strain_pred_encoder_right_pct),
                    r = suppressWarnings(cor(.data$strain_pred_encoder_right_pct, .data$strain_sono_pct)),
                    .groups = "drop") |>
  dplyr::filter(.data$n >= 30L) |>
  dplyr::mutate(precondition = classify_dynamic_precondition(.data$specimen, .data$trial_num)) |>
  dplyr::arrange(.data$specimen, .data$trial_num)

write.csv(offset_by_trial, file.path(DATA_OUT_DIR, "dynamic_precondition_offset_by_trial.csv"), row.names = FALSE)
print(offset_by_trial, n = 100)

cutoff_df <- tibble(specimen = names(DYNAMIC_PRECONDITION_CUTOFF_TRIALNUM),
                     cutoff = DYNAMIC_PRECONDITION_CUTOFF_TRIALNUM)

p_offset <- ggplot(offset_by_trial, aes(x = trial_num, y = offset_pct, color = precondition)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(data = cutoff_df, aes(xintercept = cutoff - 0.5), linetype = "dotted", color = "black") +
  geom_point(size = 2.4) +
  facet_wrap(~specimen, scales = "free_x") +
  scale_color_manual(values = c("early (preconditioning)" = "#d95f02", "later (stable)" = "#1b9e77"), drop = FALSE) +
  labs(title = "Dynamic, active (right-stim) trials: sono-strain offset by chronological trial order",
       subtitle = "Offset = mean(sono strain) - mean(encoder-predicted strain), active-window samples, %L0. Dotted line = the specimen-specific hard cutoff.",
       x = "Trial number (bender_NN, chronological)", y = "Mean offset (percentage points of L0)", color = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
fout_offset <- file.path(OUT_DIR, "dynamic_precondition_offset_vs_trialorder.png")
ggplot2::ggsave(fout_offset, p_offset, width = 10, height = 4.5, dpi = 150)
cli::cli_alert_success("Saved {fout_offset}")

# =============================================================================
# Part B: per-cycle muscle power for dynamic trials (freshly computed)
# =============================================================================
cli::cli_h1("Part B: per-cycle muscle power output, all dynamic trials bass16/17/18")

`%||%` <- function(x, y) if (is.null(x)) y else x

# Duplicated from run_fv_fl_power_pipeline.R / compare_specimen_specific_
# properties.R (both already note this is a deliberate duplication -- no
# importable module boundary exists yet short of running run_fv_fl_power_
# pipeline.R's whole top-level script). Keep all 3 copies in sync if any
# changes.
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

#' cycletype == "act" filtering below is now OBSERVED-stim-based (03_analyze.R
#' set_cycle_types() fix, 2026-07-24) -- previously relied on a design flag
#' (is_active_by_cycle) that used a different, unsynchronized cycle counter
#' and mislabelled the active window on every dynamic trial. summarize_
#' muscle_cycles() also no longer splits one physical cycle into a
#' left-window row and a right-window row (see that function's docstring) --
#' each row below is one whole active cycle, both muscles' contribution
#' included.
.attach_dynamic_muscle_force <- function(td, torque_col, lidx_pos_motor, lidx_left, lidx_right,
                                          relaxation_s = RELAXATION_WINDOW_S) {
  td$.row_id <- seq_len(nrow(td))
  td2 <- set_cycle_types(td)
  msc <- tryCatch(
    calc_muscle_torque(td2, torque_col = torque_col, include_all_active_samples = TRUE, phase_match_all_rows = TRUE),
    error = function(e) { cli::cli_warn("calc_muscle_torque failed: {conditionMessage(e)}"); NULL }
  )
  if (is.null(msc) || nrow(msc) == 0L || !(".row_id" %in% names(msc))) return(NULL)
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
  rec_lidx <- dplyr::case_when(row_side == "L" ~ lidx_left, row_side == "R" ~ lidx_right, .default = NA_real_)
  force_sign <- rec_lidx * lidx_pos_motor
  if (!all(!is.finite(force_sign))) {
    msc$muscle_torque.Nm <- force_sign[msc$.row_id] * msc$muscle_torque.Nm
    msc$stim <- row_side[msc$.row_id]
  }
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

.collect_dynamic_cycles <- function(specimen, source_dir) {
  manifest <- parse_trial_directory(source_dir)
  files <- manifest$fullpath[manifest$protocol == "dynamic"]
  purrr::map_dfr(files, function(f) {
    trial_id <- tools::file_path_sans_ext(basename(f))
    td <- tryCatch(.load_one(f), error = function(e) NULL)
    if (is.null(td) || all(as.character(td$stim) == "0")) return(tibble())
    geom_f <- .read_file_level_geometry(f)
    msc <- tryCatch(
      .attach_dynamic_muscle_force(td, "torque_inertia_corrected_Nm", geom_f$lidx_pos_motor, geom_f$lidx_left, geom_f$lidx_right),
      error = function(e) { cli::cli_alert_danger("{trial_id}: {conditionMessage(e)}"); NULL }
    )
    if (is.null(msc) || nrow(msc) == 0L) return(tibble())
    msc$muscle_mass.kg <- geom_f$muscle$muscle_mass_kg
    cyc <- summarize_muscle_cycles(msc)
    if (nrow(cyc) == 0L) return(tibble())
    cyc$specimen  <- specimen
    cyc$trial_id  <- trial_id
    cyc$trial_num <- extract_bender_trial_num(trial_id)
    cli::cli_alert_success("{trial_id}: {nrow(cyc)} active cycles")
    cyc
  })
}

dyn_cycles <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  .collect_dynamic_cycles(specimen, raw_source_dir(subfolder))
})
dyn_cycles <- dyn_cycles |>
  dplyr::filter(is.finite(.data$trial_num)) |>
  dplyr::mutate(precondition = classify_dynamic_precondition(.data$specimen, .data$trial_num))

write.csv(dyn_cycles, file.path(DATA_OUT_DIR, "dynamic_precondition_power_percycle.csv"), row.names = FALSE)
cli::cli_alert_success("Saved {nrow(dyn_cycles)} active cycles -> data_processed/dynamic_precondition_power_percycle.csv")

# =============================================================================
# Power vs. trial order
# =============================================================================
trial_power <- dyn_cycles |>
  dplyr::filter(is.finite(.data$avg_power.W)) |>
  dplyr::group_by(.data$specimen, .data$trial_id, .data$trial_num, .data$precondition) |>
  dplyr::summarise(n_cycles = dplyr::n(),
                    mean_avg_power_W = mean(.data$avg_power.W, na.rm = TRUE),
                    mean_curvature_invm = mean(abs(.data$curvature.invm), na.rm = TRUE),
                    mean_freq_hz = mean(.data$freq.Hz, na.rm = TRUE),
                    .groups = "drop") |>
  dplyr::arrange(.data$specimen, .data$trial_num)
write.csv(trial_power, file.path(DATA_OUT_DIR, "dynamic_precondition_power_by_trial.csv"), row.names = FALSE)
print(trial_power, n = 100)

p_power_trial <- ggplot(trial_power, aes(x = trial_num, y = mean_avg_power_W, color = precondition)) +
  geom_vline(data = cutoff_df, aes(xintercept = cutoff - 0.5), linetype = "dotted", color = "black") +
  geom_point(size = 2.4) +
  facet_wrap(~specimen, scales = "free") +
  scale_color_manual(values = c("early (preconditioning)" = "#d95f02", "later (stable)" = "#1b9e77"), drop = FALSE) +
  labs(title = "Dynamic trials: trial-mean cycle-averaged muscle power by chronological trial order",
       subtitle = "Dotted line = the specimen-specific preconditioning cutoff. Same trial-order axis as the offset plot above.",
       x = "Trial number (bender_NN, chronological)", y = "Trial-mean avg_power.W", color = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
fout_power_trial <- file.path(OUT_DIR, "dynamic_precondition_power_vs_trialorder.png")
ggplot2::ggsave(fout_power_trial, p_power_trial, width = 10, height = 4.5, dpi = 150)
cli::cli_alert_success("Saved {fout_power_trial}")

# =============================================================================
# Power vs. curvature (bending amplitude), controlling for comparable bending
# and stim by plotting against curvature.invm directly rather than pooling
# across amplitudes
# =============================================================================
p_power_curve <- ggplot(dplyr::filter(dyn_cycles, is.finite(.data$avg_power.W), is.finite(.data$curvature.invm)),
                         aes(x = abs(curvature.invm), y = avg_power.W, color = precondition)) +
  geom_point(alpha = 0.4, size = 1.2) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~specimen, scales = "free") +
  scale_color_manual(values = c("early (preconditioning)" = "#d95f02", "later (stable)" = "#1b9e77"), drop = FALSE) +
  labs(title = "Dynamic, active cycles: muscle power vs. bending amplitude, early vs. later trials",
       subtitle = "Comparable-bending check: do early and later trials differ in power at the SAME curvature?",
       x = "|Curvature| (1/m)", y = "Cycle-averaged muscle power (W)", color = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
fout_power_curve <- file.path(OUT_DIR, "dynamic_precondition_power_vs_curvature.png")
ggplot2::ggsave(fout_power_curve, p_power_curve, width = 10, height = 4.5, dpi = 150)
cli::cli_alert_success("Saved {fout_power_curve}")

# =============================================================================
# Regression: does precondition status predict power, controlling for
# curvature and frequency (the "comparable bending and stimulation" ask)?
# =============================================================================
cli::cli_h1("Regression: avg_power.W ~ precondition + curvature + freq, per specimen")
for (sp in unique(dyn_cycles$specimen)) {
  dsp <- dplyr::filter(dyn_cycles, .data$specimen == sp, is.finite(.data$avg_power.W),
                        is.finite(.data$curvature.invm), is.finite(.data$freq.Hz), !is.na(.data$precondition))
  if (dplyr::n_distinct(dsp$precondition) < 2L || nrow(dsp) < 10L) {
    cli::cli_alert_warning("{sp}: insufficient data/levels for regression (n={nrow(dsp)}, levels={dplyr::n_distinct(dsp$precondition)})")
    next
  }
  fit <- lm(avg_power.W ~ precondition + abs(curvature.invm) + freq.Hz, dsp)
  cli::cli_h2(sp)
  print(summary(fit)$coefficients)
}

cli::cli_alert_success("diag_precondition_power_check.R complete -- outputs in {OUT_DIR}/ and {DATA_OUT_DIR}/")
