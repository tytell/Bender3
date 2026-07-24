# diag_precondition_sono_length_activeVsPassive.R
# PI direction, 2026-07-24: "ground truth [the geometric-vs-sono power
# divergence] using strictly the difference between active and passive."
#
# calc_muscle_torque() already does exactly this kind of active-minus-passive
# differencing for FORCE (phase-matched act torque minus averaged passive
# torque at the same phase/cycle-position) -- that machinery is what made the
# torque signal in diag_sono_vs_geometric_dynamic_power.R trustworthy. But
# the LENGTH/velocity side of that script was never treated the same way:
# both the geometric (commanded-angle) and sono velocity traces are used as
# raw ABSOLUTE signals, never differenced against a phase-matched passive
# baseline. This script closes that asymmetry by reusing calc_muscle_torque()
# itself -- generically written around a "torque_col" -- with
# torque_col = the 40 Hz-filtered sono strain instead of a force channel.
# Output is directly analogous to muscle_torque.Nm: "sono_strain_excess_pct"
# = (active-cycle sono strain at a given phase) - (phase-matched AVERAGE
# passive-cycle sono strain at that same phase, same trial). A muscle that is
# actively shortening beyond whatever the passive/motor-driven tissue would
# do at that same phase should show a negative excess (shortening = strain
# increase by this pipeline's sign convention... see note below) during the
# stim window; if the early-trial "excess shortening" flagged throughout this
# investigation is REAL active muscle behavior (not a geometric artifact),
# this ground-truthed excess should ALSO show the same early > later pattern.
#
# Sign convention note: strain_sono_pct is defined shortening-positive
# (attach_sono_strain(): strain = -(mm - L0)/L0*100). So MORE shortening =
# HIGHER (more positive) strain, and "excess active shortening" shows up as
# sono_strain_excess_pct > 0 (active strain HIGHER than the passive baseline
# at the same phase) -- opposite intuition from a raw mm delta.
#
# Scope: dynamic trials, RIGHT-STIM windows only (sono's only wired
# channel), ALL trials (early + later -- testing whether the pattern
# tracks trial order the same way the raw offset/power findings did).
#
# Run with:  Rscript R/diag_precondition_sono_length_activeVsPassive.R
# Outputs -> figures: 02_processed/figs_diagnostic/ (FIGS_DIAGNOSTIC_DIR)
#            data:    02_processed/data_processed/

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr); library(ggplot2); library(cli); library(rhdf5); library(readr); library(signal)
})

# Same 40 Hz zero-phase Butterworth LP recipe as diag_sono_vs_geometric_
# dynamic_power.R -- duplicated (that script is a standalone analysis, not a
# sourceable library). See that script's header for the filter-choice
# rationale (from diag_sono_smoothing.R's comparison).
.butter_lp_sono <- function(x, cutoff_hz, sample_rate_hz, order = 4L) {
  nyq <- sample_rate_hz / 2.0
  if (cutoff_hz >= nyq) return(x)
  ok <- is.finite(x)
  if (sum(ok) < 20L) return(x)
  filt <- signal::butter(order, cutoff_hz / nyq, type = "low")
  out  <- x
  out[ok] <- signal::filtfilt(filt, x[ok])
  out
}
SONO_LOWPASS_CUTOFF_HZ <- 40.0

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")
src("00_load_bender_flat.R")
src("01_calibrate.R")
src("02_deconvolve.R")
src("muscle_geometry.R")
src("03_analyze.R")
src("parse_trial_filename.R")
src("plot_force_vs_time.R")           # .detect_stim_events(), RELAXATION_WINDOW_S
src("plot_strain_validation.R")       # attach_measured_strain()
src("plot_angle_sono_validation.R")   # attach_sono_strain(), .sono_reference_length_mm()
src("dynamic_trial_precondition.R")

OUT_DIR      <- FIGS_DIAGNOSTIC_DIR
SUMMARY_DIR  <- FIGS_SUMMARY_DIR   # cross-fish (3-specimen) version lives here too
DATA_OUT_DIR <- file.path(.crittergripper_root(), "02_processed", "data_processed")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUMMARY_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SPECIMEN_SUBFOLDERS <- c(bass16 = BASS16_RAW_SUBFOLDER, bass17 = BASS17_RAW_SUBFOLDER, bass18 = BASS18_RAW_SUBFOLDER)
SPECIMEN_COLORS     <- c(bass16 = "#1b9e77", bass17 = "#d95f02", bass18 = "#7570b3")

# =============================================================================
# Per-file geometry (same deliberate-duplication pattern as the other
# diag_precondition_*.R / diag_sono_vs_geometric_dynamic_power.R scripts).
# =============================================================================
.read_file_level_geometry <- function(filename) {
  h5 <- rhdf5::H5Fopen(filename, "H5F_ACC_RDONLY")
  on.exit(try(rhdf5::H5Fclose(h5), silent = TRUE), add = TRUE)
  m_attrs <- tryCatch(rhdf5::h5readAttributes(h5, "/metadata"), error = function(e) list())
  dbl1 <- function(v, default = NA_real_) {
    v <- suppressWarnings(as.numeric(v[1L]))
    if (length(v) == 0L || is.na(v)) default else v
  }
  list(local_body_width_mm = dbl1(m_attrs[["measurement_specimen_local_body_width_millimeter"]]),
       measured_depth_mm   = dbl1(m_attrs[["measurement_target_muscle_depth_millimeter"]]),
       ai_rate_hz          = dbl1(m_attrs[["daq_ai_sample_rate_hz"]], default = 1000.0),
       lidx_pos_motor      = dbl1(m_attrs[["daq_specimen_lateral_index_on_positive_motor_side"]]),
       lidx_left           = dbl1(m_attrs[["daq_specimen_side_index_left"]]),
       lidx_right          = dbl1(m_attrs[["daq_specimen_side_index_right"]]))
}

# =============================================================================
# Per-trial: reuse calc_muscle_torque()'s phase-matched active-minus-passive
# machinery, but feed it the 40 Hz-filtered sono STRAIN instead of a torque
# channel. Restrict to RIGHT-STIM windows afterward (sono only instruments
# the right muscle).
# =============================================================================
.collect_one_trial <- function(f, specimen) {
  trial_id <- tools::file_path_sans_ext(basename(f))
  td <- tryCatch(load_bender_flat(f, do_filter = TRUE, loadtorques = "x"), error = function(e) NULL)
  if (is.null(td) || all(as.character(td$stim) == "0")) return(tibble())

  geom <- .read_file_level_geometry(f)
  td <- tryCatch(attach_predicted_strain(td, geom$local_body_width_mm, geom$measured_depth_mm),
                 error = function(e) td)
  td <- tryCatch(attach_sono_strain(td, f, geom$lidx_right, geom$lidx_pos_motor), error = function(e) td)
  if (!("sono_right_mm" %in% names(td)) || all(!is.finite(td$sono_right_mm))) return(tibble())

  # Condition sono BEFORE any row restriction, on the full continuous trial
  # (see diag_sono_vs_geometric_dynamic_power.R header for the filtfilt-
  # continuity rationale).
  L0 <- .sono_reference_length_mm(td$angle.deg, td$sono_right_mm)
  if (!is.finite(L0) || L0 <= 0) return(tibble())
  sono_filt_mm <- .butter_lp_sono(td$sono_right_mm, SONO_LOWPASS_CUTOFF_HZ, geom$ai_rate_hz)
  td$strain_sono_filt_pct <- -(sono_filt_mm - L0) / L0 * 100.0

  td$.row_id <- seq_len(nrow(td))
  td2 <- set_cycle_types(td)

  # Reuse calc_muscle_torque() verbatim, just pointed at the sono-strain
  # column instead of a torque column -- same phase-matched act-minus-
  # (averaged)-passive logic already vetted for force.
  msc <- tryCatch(
    calc_muscle_torque(td2, torque_col = "strain_sono_filt_pct",
                        include_all_active_samples = TRUE, phase_match_all_rows = TRUE),
    error = function(e) { cli::cli_warn("{trial_id}: calc_muscle_torque failed: {conditionMessage(e)}"); NULL }
  )
  if (is.null(msc) || nrow(msc) == 0L || !(".row_id" %in% names(msc))) return(tibble())

  events <- tryCatch(.detect_stim_events(td), error = function(e) tibble::tibble())
  row_side <- rep(NA_character_, nrow(td))
  if (nrow(events) > 0L) {
    events <- dplyr::arrange(events, .data$stim_t0_s)
    for (i in seq_len(nrow(events))) {
      e <- events[i, ]
      in_win <- td$t.s >= e$stim_t0_s & td$t.s <= (e$stim_t1_s + RELAXATION_WINDOW_S)
      row_side[in_win] <- e$muscle_side
    }
  }
  # RIGHT-STIM ONLY -- sono only instruments the right muscle.
  msc$row_side <- row_side[msc$.row_id]
  msc <- dplyr::filter(msc, .data$row_side == "R")
  if (nrow(msc) == 0L) return(tibble())
  if ("cycletype" %in% names(msc)) msc <- dplyr::filter(msc, .data$cycletype == "act")
  if (nrow(msc) == 0L) return(tibble())

  # muscle_torque.Nm here is actually a STRAIN delta (pct) -- rename for
  # clarity; keep the raw filtered strain (.torque) and the passive
  # baseline (avg_pass_torque.Nm) alongside for inspection.
  tibble::tibble(
    specimen                = specimen,
    trial_id                = trial_id,
    trial_num               = extract_bender_trial_num(trial_id),
    cycle                   = msc$cycle,
    t.s                     = msc$t.s,
    sono_strain_active_pct  = msc$.torque,
    sono_strain_passive_baseline_pct = msc$avg_pass_torque.Nm,
    sono_strain_excess_pct  = msc$muscle_torque.Nm
  )
}

cli::cli_h1("Collecting RIGHT-STIM dynamic samples, active-minus-passive sono strain, bass16/17/18")
all_rows <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  manifest <- parse_trial_directory(raw_source_dir(subfolder))
  files <- manifest$fullpath[manifest$protocol == "dynamic"]
  purrr::map_dfr(files, function(f) {
    out <- tryCatch(.collect_one_trial(f, specimen), error = function(e) {
      cli::cli_alert_danger("{basename(f)}: {conditionMessage(e)}"); tibble()
    })
    if (nrow(out) > 0L) cli::cli_alert_success("{tools::file_path_sans_ext(basename(f))}: {nrow(out)} right-stim active samples")
    out
  })
})
all_rows <- all_rows |>
  dplyr::filter(is.finite(.data$trial_num), is.finite(.data$sono_strain_excess_pct)) |>
  dplyr::mutate(precondition = classify_dynamic_precondition(.data$specimen, .data$trial_num))

write.csv(all_rows, file.path(DATA_OUT_DIR, "sono_length_activeVsPassive_persample.csv"), row.names = FALSE)
cli::cli_alert_success("Saved {nrow(all_rows)} samples -> data_processed/sono_length_activeVsPassive_persample.csv")

# =============================================================================
# Per-cycle, then per-trial aggregation of the excess.
# =============================================================================
cyc <- all_rows |>
  dplyr::group_by(.data$specimen, .data$trial_id, .data$trial_num, .data$precondition, .data$cycle) |>
  dplyr::summarise(n_samples = dplyr::n(),
                    mean_excess_pct = mean(.data$sono_strain_excess_pct, na.rm = TRUE),
                    peak_excess_pct = max(abs(.data$sono_strain_excess_pct), na.rm = TRUE),
                    .groups = "drop")

trial_excess <- cyc |>
  dplyr::group_by(.data$specimen, .data$trial_id, .data$trial_num, .data$precondition) |>
  dplyr::summarise(n_cycles = dplyr::n(),
                    mean_excess_pct = mean(.data$mean_excess_pct, na.rm = TRUE),
                    .groups = "drop") |>
  dplyr::arrange(.data$specimen, .data$trial_num)

write.csv(trial_excess, file.path(DATA_OUT_DIR, "sono_length_activeVsPassive_bytrial.csv"), row.names = FALSE)
cli::cli_h1("Trial-level mean active-minus-passive sono strain excess, by trial order")
print(trial_excess, n = 60)

# =============================================================================
# Plot 1: trial-mean excess by trial order, per specimen -- does the
# phase-matched active-vs-passive GROUND-TRUTHED excess also show the same
# early > later pattern as the raw offset/power findings?
# =============================================================================
cutoff_df <- tibble::tibble(specimen = names(DYNAMIC_PRECONDITION_CUTOFF_TRIALNUM),
                             cutoff = DYNAMIC_PRECONDITION_CUTOFF_TRIALNUM)

p1 <- ggplot(trial_excess, aes(x = .data$trial_num, y = .data$mean_excess_pct, color = .data$specimen, shape = .data$precondition)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(data = cutoff_df, aes(xintercept = cutoff - 0.5), linetype = "dotted", color = "black", inherit.aes = FALSE) +
  geom_point(size = 2.8) +
  facet_wrap(~specimen, scales = "free_x") +
  scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen") +
  scale_shape_manual(values = c("early (preconditioning)" = 17, "later (stable)" = 16), name = NULL) +
  labs(title = "Ground-truthed active-vs-passive sono strain excess, by trial order (dynamic, right-stim)",
       subtitle = "Excess = active-cycle sono strain minus PHASE-MATCHED AVERAGE passive-cycle sono strain, same trial (reuses\ncalc_muscle_torque()'s act-minus-passive logic, same recipe already vetted for force). Positive = active cycle shortens\nMORE than the passive baseline at that phase (shortening-positive convention). Dotted line = existing hard cutoff.",
       x = "Trial number (bender_NN, chronological)", y = "Trial-mean active-minus-passive sono strain excess (pct L0)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
fout1 <- file.path(OUT_DIR, "dynamic_precondition_sonoLengthExcess_vs_trialorder.png")
ggplot2::ggsave(fout1, p1, width = 11, height = 5, dpi = 150)
cli::cli_alert_success("Saved {fout1}")

# Cross-fish (3-specimen) summary copy -- genuinely pooled content, belongs
# in figs_summary/ per the figure-placement rules.
fsum <- file.path(SUMMARY_DIR, "sonoLengthExcess_activeVsPassive_byTrialOrder.png")
ggplot2::ggsave(fsum, p1, width = 11, height = 5, dpi = 150)
cli::cli_alert_success("Saved {fsum}")

# =============================================================================
# Plot 2: cycle-level boxplot, early vs later -- distribution summary.
# =============================================================================
n_lab <- cyc |>
  dplyr::count(.data$precondition, name = "n_cycles") |>
  dplyr::mutate(x_label = sprintf("%s\n(n=%d cycles)", .data$precondition, .data$n_cycles))
cyc2 <- dplyr::left_join(cyc, dplyr::select(n_lab, "precondition", "x_label"), by = "precondition")
cyc2$x_label <- factor(cyc2$x_label, levels = dplyr::arrange(n_lab, .data$precondition)$x_label)

p2 <- ggplot(cyc2, aes(x = .data$x_label, y = .data$mean_excess_pct)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, fill = "grey95", color = "grey40") +
  geom_jitter(aes(color = .data$specimen), width = 0.12, size = 2.0, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen") +
  labs(title = "Dynamic, right-stim cycles: ground-truthed active-vs-passive sono\nstrain excess, early vs. later",
       subtitle = "Each point = one active cycle's mean excess (active sono strain minus\nphase-matched passive baseline, same trial).",
       x = NULL, y = "Cycle-mean active-minus-passive sono strain excess (pct L0)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
fout2 <- file.path(OUT_DIR, "dynamic_precondition_sonoLengthExcess_boxplot.png")
ggplot2::ggsave(fout2, p2, width = 8, height = 5.5, dpi = 150)
cli::cli_alert_success("Saved {fout2}")

# =============================================================================
# Summary: early vs later median excess (does it echo the offset/power
# early>later pattern, ground-truthed against passive baseline?)
# =============================================================================
gap_tbl <- cyc |>
  dplyr::group_by(.data$precondition) |>
  dplyr::summarise(n_cycles = dplyr::n(), median_excess_pct = median(.data$mean_excess_pct, na.rm = TRUE),
                    mean_excess_pct = mean(.data$mean_excess_pct, na.rm = TRUE), .groups = "drop")
cli::cli_h1("Early vs. later median/mean active-minus-passive sono strain excess (pct L0)")
print(gap_tbl)
write.csv(gap_tbl, file.path(DATA_OUT_DIR, "sono_length_activeVsPassive_earlyLaterGap.csv"), row.names = FALSE)

cli::cli_alert_success("diag_precondition_sono_length_activeVsPassive.R complete -- outputs in {OUT_DIR}/ and {DATA_OUT_DIR}/")
