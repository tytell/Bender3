# diag_sono_dynamic_by_amplitude.R
# Diagnostic (not a pipeline deliverable): re-renders the "dynamic" panel of
# strainValidSonoEnc.png (measured sono strain vs. encoder-predicted strain,
# ALL samples) split into SEPARATE panels by that sample's CYCLE'S predicted
# peak-to-peak amplitude, each with its own (equal-aspect, symmetric) axis
# limits -- so the low-amplitude cycles (which pool together with
# high-amplitude cycles in the single combined dynamic panel, making the
# "tall not wide" departure from 1:1 hard to see at a glance) get their own
# zoomed-in view.
#
# Motivation (see diag_sono_amp_ratio_multitrial.R): per-cycle amplitude
# ratio (measured p2p / predicted p2p) is ~5x at predicted p2p < 2%, falling
# to ~1.3x at predicted p2p 20-40%, and is essentially IDENTICAL between
# stim and no-stim cycles at matched amplitude -- i.e. the excess is an
# amplitude-dependent additive-ish effect, not a stim-driven one. This
# figure makes that visible directly in the measured-vs-predicted sample
# cloud (not just the per-cycle summary numbers), pooled across every
# dynamic trial in bass16/17/18.
#
# Bin edges match diag_sono_amp_ratio_multitrial.R's MIN_P2P_PRED_PCT /
# amp_bin cut points for direct comparison with that script's summary table.
#
# Run with:  Rscript R/diag_sono_dynamic_by_amplitude.R
# Outputs -> 02_processed/figs_diagnostic/ (FIGS_DIAGNOSTIC_DIR, paths_config.R
#            -- same canonical location as the original per-specimen
#            strainValidSonoEnc.png figures this diagnoses).

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2); library(cli); library(rhdf5); library(patchwork)
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

OUT_DIR <- FIGS_DIAGNOSTIC_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SPECIMEN_SUBFOLDERS <- c(
  bass16 = BASS16_RAW_SUBFOLDER,
  bass17 = BASS17_RAW_SUBFOLDER,
  bass18 = BASS18_RAW_SUBFOLDER
)

# Same bin edges as diag_sono_amp_ratio_multitrial.R's ad hoc summary cut,
# now promoted to a named constant since two scripts share it.
AMP_BIN_BREAKS <- c(0, 2, 5, 10, 20, 40)
AMP_BIN_LABELS <- c("0-2%", "2-5%", "5-10%", "10-20%", "20-40%")

STIM_COLORS <- c("no stim" = "grey75", "left stim" = "#1d4ed8", "right stim" = "#b91c1c")

# =============================================================================
# Helpers (same as diag_sono_amp_ratio_multitrial.R)
# =============================================================================

.dynamic_trial_files <- function(subfolder) {
  d <- raw_source_dir(subfolder)
  sort(list.files(d, pattern = "_dynamic\\.h5$", full.names = TRUE))
}

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

.read_cycle_freq <- function(h5path) {
  h5 <- H5Fopen(h5path, "H5F_ACC_RDONLY")
  on.exit(try(H5Fclose(h5), silent = TRUE), add = TRUE)
  ci_raw  <- tryCatch(h5read(h5, "/timeseries/cycle_index"), error = function(e) NULL)
  freq_ds <- tryCatch(h5read(h5, "/metadata/index_cycle_frequency_hertz"), error = function(e) NULL)
  list(ci = as.integer(ci_raw), freq_by_cycle = as.numeric(freq_ds))
}

.p2p_per_cycle <- function(x, ci_vec, cycles) {
  vapply(cycles, function(c) {
    v <- x[ci_vec == c]
    if (sum(is.finite(v)) < 4L) return(NA_real_)
    max(v, na.rm = TRUE) - min(v, na.rm = TRUE)
  }, numeric(1L))
}

#' Process one dynamic trial -> per-SAMPLE tibble (not per-cycle), each row
#' tagged with its own cycle's predicted-amplitude bin. Returns empty tibble
#' on any load/geometry/sono failure (logged, not fatal to the batch).
.process_trial_samples <- function(specimen, h5path) {
  trial_label <- tools::file_path_sans_ext(basename(h5path))

  td <- tryCatch(load_bender_flat(h5path, do_filter = FALSE, loadtorques = "x"),
                 error = function(e) { cli::cli_alert_warning("{trial_label}: load failed -- {conditionMessage(e)}"); NULL })
  if (is.null(td) || nrow(td) == 0L) return(tibble())

  geom <- .read_geom_and_lidx(h5path)
  if (!is.finite(geom$width_mm) || !is.finite(geom$lidx_right) || !is.finite(geom$lidx_pos_motor)) return(tibble())

  info <- .read_cycle_freq(h5path)
  if (is.null(info$ci) || length(info$ci) != nrow(td)) return(tibble())
  td$cycle_idx <- info$ci

  td <- attach_predicted_strain(
    td,
    local_body_width_mm      = geom$width_mm,
    measured_muscle_depth_mm = geom$depth_mm,
    active_mask               = rep(FALSE, nrow(td))
  )
  td <- attach_measured_strain(td)
  td <- attach_sono_strain(td, h5path, geom$lidx_right, geom$lidx_pos_motor)
  if (!any(is.finite(td$strain_sono_pct))) return(tibble())

  stim_state_fac <- tryCatch(.stim_window_state_label(td), error = function(e) NULL)
  if (is.null(stim_state_fac)) return(tibble())
  td$stim_state <- as.character(stim_state_fac)

  cyc_list <- sort(unique(td$cycle_idx[td$cycle_idx >= 1L]))
  if (length(cyc_list) == 0L) return(tibble())
  p2p_pred <- .p2p_per_cycle(td$strain_pred_encoder_right_pct, td$cycle_idx, cyc_list)
  bin_by_cycle <- cut(p2p_pred, breaks = AMP_BIN_BREAKS, labels = AMP_BIN_LABELS)
  cycle_bin_lut <- setNames(as.character(bin_by_cycle), as.character(cyc_list))

  cli::cli_alert_success("{trial_label}: {nrow(td)} samples, {length(cyc_list)} cycles")

  td |>
    dplyr::filter(.data$cycle_idx >= 1L,
                  is.finite(.data$strain_pred_encoder_right_pct),
                  is.finite(.data$strain_sono_pct)) |>
    dplyr::transmute(
      specimen = specimen, trial = trial_label,
      pred_pct = .data$strain_pred_encoder_right_pct,
      meas_pct = .data$strain_sono_pct,
      stim_state = .data$stim_state,
      amp_bin = cycle_bin_lut[as.character(.data$cycle_idx)]
    )
}

# =============================================================================
# Batch over every dynamic trial, every specimen
# =============================================================================
cli::cli_h1("Loading all dynamic trials (bass16/17/18) for per-sample amplitude-bin plot")

all_files <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  tibble(specimen = specimen, h5path = .dynamic_trial_files(subfolder))
})
cli::cli_alert_info("Total dynamic trials: {nrow(all_files)}")

samples <- purrr::pmap_dfr(all_files, function(specimen, h5path) {
  tryCatch(.process_trial_samples(specimen, h5path),
           error = function(e) {
             cli::cli_alert_danger("{basename(h5path)}: unexpected error -- {conditionMessage(e)}")
             tibble()
           })
})

samples <- samples |>
  dplyr::filter(!is.na(.data$amp_bin)) |>
  dplyr::mutate(
    amp_bin = factor(.data$amp_bin, levels = AMP_BIN_LABELS),
    stim_state = factor(.data$stim_state, levels = names(STIM_COLORS))
  )
cli::cli_alert_success("Pooled {nrow(samples)} samples across {dplyr::n_distinct(samples$trial)} trials")

# =============================================================================
# One equal-aspect, symmetric-limit panel per amplitude bin
# =============================================================================
cli::cli_h1("Building per-amplitude-bin panels")

.build_bin_panel <- function(bin_label) {
  df <- dplyr::filter(samples, .data$amp_bin == bin_label)
  if (nrow(df) < 20L) return(NULL)

  lim <- max(abs(range(c(df$pred_pct, df$meas_pct), na.rm = TRUE))) * 1.05
  r    <- suppressWarnings(cor(df$pred_pct, df$meas_pct, use = "complete.obs"))
  rmse <- sqrt(mean((df$pred_pct - df$meas_pct)^2, na.rm = TRUE))

  ggplot(df, aes(x = .data$pred_pct, y = .data$meas_pct, color = .data$stim_state)) +
    geom_point(shape = 1, size = 0.5, alpha = 0.15, stroke = 0.25) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
    scale_color_manual(values = STIM_COLORS, drop = FALSE) +
    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(title = sprintf("Predicted p2p %s (n=%s)", bin_label, format(nrow(df), big.mark = ",")),
         subtitle = sprintf("r=%.3f, RMSE=%.2f%%", r, rmse),
         x = "Predicted strain (%, encoder)", y = "Measured strain (%, sono)") +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", plot.title = element_text(size = 9.5), plot.subtitle = element_text(size = 9))
}

panels <- purrr::compact(purrr::map(AMP_BIN_LABELS, .build_bin_panel))

p_combined <- patchwork::wrap_plots(panels, ncol = 3) +
  patchwork::plot_annotation(
    title = "Dynamic trials only: measured (sono) vs. predicted (encoder) strain, split by cycle amplitude",
    subtitle = "Pooled across all bass16/17/18 dynamic trials, RIGHT muscle only. Each panel has its own equal-aspect, symmetric axis limits.\nLow-amplitude cycles depart much further from the dashed 1:1 line (proportionally) than high-amplitude cycles.",
    caption = "Color: grey = no stim, blue = left stim, red = right stim (per-sample stim window state)",
    theme = theme(plot.title = element_text(size = 13), plot.subtitle = element_text(size = 10),
                  plot.caption = element_text(size = 9, hjust = 0))
  )

fout <- file.path(OUT_DIR, "dynamic_by_amplitude_bin.png")
ggplot2::ggsave(fout, p_combined, width = 13, height = 9, dpi = 150)
cli::cli_alert_success("Saved {fout}")
cli::cli_alert_success("diag_sono_dynamic_by_amplitude.R complete -- outputs in {OUT_DIR}/")
