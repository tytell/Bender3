# diag_sono_amp_ratio_multitrial.R
# Diagnostic (not a pipeline deliverable): extends diag_sono_timing.R's
# amplitude-ratio check (measured sono p2p / encoder-predicted p2p, per
# cycle) from bass16 trials 17/19 only to EVERY dynamic (single_finite)
# trial across bass16, bass17, and bass18 -- and, unlike diag_sono_timing.R,
# splits the ratio by per-cycle DOMINANT stim state (no stim / left stim /
# right stim) so we can test whether the excess-shortening ("tall, not
# wide") pattern in strainValidSonoEnc.png's dynamic panel requires active
# muscle stimulation (bio/slip hypothesis) or is present even in passive
# cycles (sensor/geometry artifact hypothesis).
#
# Motivating question (see strainValidSonoEnc.png, dynamic panel): the
# measured-vs-predicted strain cloud is TALLER than it is WIDE, i.e. sono
# reports more shortening than the encoder-derived curvature geometrically
# implies. A pure fixed-ms timing/phase lag (diag_sono_timing.R's original
# hypothesis) produces a loop with roughly EQUAL x/y spread at a given lag;
# it does not by itself explain an amplitude ratio that departs from 1. An
# amp_ratio_mean > 1 that is concentrated in stimulated cycles points at
# something physical happening under active load (clamp slip, tendon/SEE
# compliance, crystal movement relative to fixed axis); an amp_ratio > 1
# that is equally present in unstimulated (passive) cycles points at a
# sensor/geometry artifact independent of muscle activity.
#
# Geometry is read PER-FILE from that file's own /metadata attributes
# (measurement_specimen_local_body_width_millimeter,
# measurement_target_muscle_depth_millimeter,
# daq_specimen_lateral_index_on_positive_motor_side/side_index_left/right)
# -- NOT hardcoded per specimen -- because body width differs across
# bass16/17/18 (confirmed 42 / 38 / 44 mm) even though the lateral-index
# sign convention is the same rig-wide.
#
# Run with:  Rscript R/diag_sono_amp_ratio_multitrial.R
# Outputs -> figures: 02_processed/figs_diagnostic/ (FIGS_DIAGNOSTIC_DIR,
#            paths_config.R -- same canonical location as the original
#            per-specimen strainValidSonoEnc.png figures this diagnoses).
#            data tables: 02_processed/data_processed/ (code-generated
#            derived data, per .cursorrules file-placement rules).

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2); library(cli); library(rhdf5)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")
src("00_load_bender_flat.R")
src("01_calibrate.R")
src("muscle_geometry.R")
src("plot_strain_validation.R")     # attach_measured_strain
src("plot_angle_sono_validation.R") # attach_sono_strain, .stim_window_state_label, STIM_STATE_LEVELS
src("plot_force_vs_time.R")         # .detect_stim_events (used by .stim_window_state_label)

OUT_DIR      <- FIGS_DIAGNOSTIC_DIR
DATA_OUT_DIR <- file.path(.crittergripper_root(), "02_processed", "data_processed")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SAMPLE_RATE_HZ <- 1000.0

SPECIMEN_SUBFOLDERS <- c(
  bass16 = BASS16_RAW_SUBFOLDER,
  bass17 = BASS17_RAW_SUBFOLDER,
  bass18 = BASS18_RAW_SUBFOLDER
)

# Minimum predicted per-cycle p2p (%) before an amplitude ratio is trusted --
# near-zero predicted amplitude (near-quiescent cycles) makes the ratio blow
# up on noise alone and is not informative about slip/artifact.
MIN_P2P_PRED_PCT <- 0.5

# =============================================================================
# Helpers
# =============================================================================

.dynamic_trial_files <- function(subfolder) {
  d <- raw_source_dir(subfolder)
  sort(list.files(d, pattern = "_dynamic\\.h5$", full.names = TRUE))
}

#' Per-file geometry + lateral-index read (NOT hardcoded -- see header).
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
    lidx_left      = dbl1("daq_specimen_side_index_left"),
    lidx_right     = dbl1("daq_specimen_side_index_right")
  )
}

#' Read cycle_index and per-cycle frequency_hertz (single_finite only).
.read_cycle_freq <- function(h5path) {
  h5 <- H5Fopen(h5path, "H5F_ACC_RDONLY")
  on.exit(try(H5Fclose(h5), silent = TRUE), add = TRUE)
  ci_raw  <- tryCatch(h5read(h5, "/timeseries/cycle_index"), error = function(e) NULL)
  freq_ds <- tryCatch(h5read(h5, "/metadata/index_cycle_frequency_hertz"), error = function(e) NULL)
  list(ci = as.integer(ci_raw), freq_by_cycle = as.numeric(freq_ds))
}

#' Per-cycle peak-to-peak amplitude.
.p2p_per_cycle <- function(x, ci_vec, cycles) {
  vapply(cycles, function(c) {
    v <- x[ci_vec == c]
    if (sum(is.finite(v)) < 4L) return(NA_real_)
    max(v, na.rm = TRUE) - min(v, na.rm = TRUE)
  }, numeric(1L))
}

#' Modal (most frequent) stim-window state within one cycle's samples.
.dominant_stim_state <- function(state_chr) {
  state_chr <- state_chr[!is.na(state_chr)]
  if (length(state_chr) == 0L) return(NA_character_)
  tab <- table(state_chr)
  names(tab)[which.max(tab)]
}

#' Process one dynamic trial: load, attach strains, attach stim state,
#' return one row per cycle (amp_ratio + metadata). Returns empty tibble on
#' any load/geometry/sono failure (logged, not fatal to the batch).
.process_trial <- function(specimen, h5path) {
  trial_label <- tools::file_path_sans_ext(basename(h5path))

  td <- tryCatch(load_bender_flat(h5path, do_filter = FALSE, loadtorques = "x"),
                 error = function(e) { cli::cli_alert_warning("{trial_label}: load failed -- {conditionMessage(e)}"); NULL })
  if (is.null(td) || nrow(td) == 0L) {
    cli::cli_alert_warning("{trial_label}: empty/failed load, skip")
    return(tibble())
  }

  geom <- .read_geom_and_lidx(h5path)
  if (!is.finite(geom$width_mm) || !is.finite(geom$lidx_right) || !is.finite(geom$lidx_pos_motor)) {
    cli::cli_alert_warning("{trial_label}: missing geometry attrs, skip")
    return(tibble())
  }

  info <- .read_cycle_freq(h5path)
  if (is.null(info$ci) || length(info$ci) != nrow(td)) {
    cli::cli_alert_warning("{trial_label}: cycle_index missing/length mismatch ({length(info$ci)} vs {nrow(td)} rows), skip")
    return(tibble())
  }
  td$cycle_idx <- info$ci
  td$freq_cycle <- dplyr::if_else(
    td$cycle_idx >= 1L & td$cycle_idx <= length(info$freq_by_cycle),
    info$freq_by_cycle[pmax(1L, td$cycle_idx)], NA_real_
  )

  td <- attach_predicted_strain(
    td,
    local_body_width_mm      = geom$width_mm,
    measured_muscle_depth_mm = geom$depth_mm,
    active_mask               = rep(FALSE, nrow(td))
  )
  td <- attach_measured_strain(td)
  td <- attach_sono_strain(td, h5path, geom$lidx_right, geom$lidx_pos_motor)

  if (!any(is.finite(td$strain_sono_pct))) {
    cli::cli_alert_warning("{trial_label}: no valid sono samples (missing calibration or channel), skip")
    return(tibble())
  }

  stim_state_fac <- tryCatch(.stim_window_state_label(td), error = function(e) NULL)
  if (is.null(stim_state_fac)) {
    cli::cli_alert_warning("{trial_label}: stim state labeling failed, skip")
    return(tibble())
  }
  td$stim_state <- as.character(stim_state_fac)

  cyc_list <- sort(unique(td$cycle_idx[td$cycle_idx >= 1L]))
  if (length(cyc_list) == 0L) {
    cli::cli_alert_warning("{trial_label}: no valid cycle_idx>=1, skip")
    return(tibble())
  }

  p2p_pred  <- .p2p_per_cycle(td$strain_pred_encoder_right_pct, td$cycle_idx, cyc_list)
  p2p_meas  <- .p2p_per_cycle(td$strain_sono_pct,               td$cycle_idx, cyc_list)
  freq_by_c <- vapply(cyc_list, function(c) {
    v <- td$freq_cycle[td$cycle_idx == c]; v <- v[is.finite(v)]
    if (length(v) == 0L) NA_real_ else v[1L]
  }, numeric(1L))
  dom_state <- vapply(cyc_list, function(c) .dominant_stim_state(td$stim_state[td$cycle_idx == c]), character(1L))
  stim_frac <- vapply(cyc_list, function(c) {
    s <- td$stim_state[td$cycle_idx == c]
    if (length(s) == 0L) return(NA_real_)
    mean(s != "no stim", na.rm = TRUE)
  }, numeric(1L))

  cli::cli_alert_success("{trial_label}: {length(cyc_list)} cycles, {sum(is.finite(p2p_pred) & is.finite(p2p_meas))} with valid p2p")

  tibble(
    specimen = specimen, trial = trial_label, cycle = cyc_list,
    freq_hz = freq_by_c, stim_state = dom_state, stim_frac = stim_frac,
    p2p_pred_pct = p2p_pred, p2p_meas_pct = p2p_meas
  )
}

# =============================================================================
# Batch over every dynamic trial, every specimen
# =============================================================================
cli::cli_h1("Discovering dynamic trials across bass16/17/18")

all_files <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  files <- .dynamic_trial_files(subfolder)
  cli::cli_alert_info("{specimen}: {length(files)} dynamic trial file(s)")
  tibble(specimen = specimen, h5path = files)
})
cli::cli_alert_info("Total dynamic trials found: {nrow(all_files)}")

results <- purrr::pmap_dfr(all_files, function(specimen, h5path) {
  tryCatch(.process_trial(specimen, h5path),
           error = function(e) {
             cli::cli_alert_danger("{basename(h5path)}: unexpected error -- {conditionMessage(e)}")
             tibble()
           })
})

if (nrow(results) == 0L) {
  cli::cli_abort("No usable cycles extracted from any trial -- nothing to analyze.")
}

results <- results |>
  dplyr::mutate(
    amp_ratio = .data$p2p_meas_pct / .data$p2p_pred_pct,
    ok_amp    = is.finite(.data$amp_ratio) & is.finite(.data$p2p_pred_pct) & .data$p2p_pred_pct >= MIN_P2P_PRED_PCT,
    stim_state = factor(.data$stim_state, levels = c("no stim", "left stim", "right stim"))
  )

write.csv(results, file.path(DATA_OUT_DIR, "amp_ratio_all_dynamic_trials.csv"), row.names = FALSE)
cli::cli_alert_success("Saved per-cycle table -> {DATA_OUT_DIR}/amp_ratio_all_dynamic_trials.csv ({nrow(results)} cycles, {sum(results$ok_amp)} usable)")

# =============================================================================
# Summary: amp_ratio by specimen x stim_state
# =============================================================================
cli::cli_h1("Summary: amplitude ratio (measured p2p / predicted p2p) by specimen x stim state")

summary_tbl <- results |>
  dplyr::filter(.data$ok_amp) |>
  dplyr::group_by(.data$specimen, .data$stim_state) |>
  dplyr::summarise(
    n_cycles   = dplyr::n(),
    n_trials   = dplyr::n_distinct(.data$trial),
    amp_ratio_mean   = mean(.data$amp_ratio, na.rm = TRUE),
    amp_ratio_median = median(.data$amp_ratio, na.rm = TRUE),
    amp_ratio_sd     = sd(.data$amp_ratio, na.rm = TRUE),
    pct_gt_1 = mean(.data$amp_ratio > 1, na.rm = TRUE) * 100.0,
    .groups = "drop"
  ) |>
  dplyr::arrange(.data$specimen, .data$stim_state)

print(summary_tbl, n = 100)
write.csv(summary_tbl, file.path(DATA_OUT_DIR, "amp_ratio_summary_by_specimen_stimstate.csv"), row.names = FALSE)

overall_tbl <- results |>
  dplyr::filter(.data$ok_amp) |>
  dplyr::group_by(.data$stim_state) |>
  dplyr::summarise(
    n_cycles = dplyr::n(),
    amp_ratio_mean   = mean(.data$amp_ratio, na.rm = TRUE),
    amp_ratio_median = median(.data$amp_ratio, na.rm = TRUE),
    amp_ratio_sd     = sd(.data$amp_ratio, na.rm = TRUE),
    .groups = "drop"
  )
cli::cli_h2("Pooled across all specimens/trials:")
print(overall_tbl)

# =============================================================================
# Plots
# =============================================================================
cli::cli_h1("Plots")

df_plot <- results |> dplyr::filter(.data$ok_amp, !is.na(.data$stim_state))

p_box <- ggplot(df_plot, aes(x = .data$stim_state, y = .data$amp_ratio, fill = .data$stim_state)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_boxplot(outlier.alpha = 0.15, width = 0.6) +
  facet_wrap(~specimen) +
  scale_fill_manual(values = c("no stim" = "grey80", "left stim" = "#1d4ed8", "right stim" = "#b91c1c")) +
  labs(title = "Per-cycle amplitude ratio (sono p2p / encoder-predicted p2p), by dominant stim state",
       subtitle = "ratio > 1 means the muscle shortened MORE than joint kinematics implies for that cycle\nAll dynamic trials, bass16/17/18, cycles with predicted p2p >= 0.5%",
       x = NULL, y = "Amplitude ratio (measured / predicted)") +
  theme_bw(base_size = 11) + theme(legend.position = "none")
ggplot2::ggsave(file.path(OUT_DIR, "amp_ratio_by_stimstate_boxplot.png"), p_box, width = 9, height = 4.5, dpi = 150)

p_freq <- ggplot(df_plot, aes(x = .data$freq_hz, y = .data$amp_ratio, color = .data$stim_state)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_jitter(width = 0.05, height = 0, alpha = 0.35, size = 0.9) +
  facet_wrap(~specimen) +
  scale_color_manual(values = c("no stim" = "grey60", "left stim" = "#1d4ed8", "right stim" = "#b91c1c")) +
  labs(title = "Amplitude ratio vs. cycle frequency, by dominant stim state",
       subtitle = "If real slip/bio: ratio should be highest/most variable for stimulated cycles, and may scale with frequency",
       x = "Cycle frequency (Hz)", y = "Amplitude ratio (measured / predicted)", color = "Stim state") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggplot2::ggsave(file.path(OUT_DIR, "amp_ratio_vs_freq_by_stimstate.png"), p_freq, width = 10, height = 4.5, dpi = 150)

cli::cli_alert_success("Saved amp_ratio_by_stimstate_boxplot.png and amp_ratio_vs_freq_by_stimstate.png")
cli::cli_alert_success("diag_sono_amp_ratio_multitrial.R complete -- outputs in {OUT_DIR}/")
