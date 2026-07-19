# diag_torque_smoothing.R
# Diagnostic (not a pipeline deliverable): compares RAW calibrated torque vs.
# inertia-corrected torque vs. several inertia-corrected + low-pass-filtered
# variants, for every isometric step, so the PI can visually pick a smoothing
# cutoff -- then that choice gets applied uniformly across the pipeline
# (build_segmented_force_timeseries()/build_isovelocity_force_timeseries()'s
# NOISE_FILTER_HZ, and .smooth_trace_display_only()'s DISPLAY_SMOOTH_HZ).
#
# Also produces a pooled power-spectral-density plot of the UNSMOOTHED
# inertia-corrected signal during the active stim window, to directly test
# the ~32 Hz wave observed by eye in forceTime_isometric_legacy.png.
#
# Run with:  Rscript R/diag_torque_smoothing.R

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(ggplot2); library(cli); library(rhdf5)
})

SOURCE_DIR <- "/Users/yjimenez/Library/CloudStorage/OneDrive-ProvidenceCollege/01_JimenezLab/01_PermanentArchive/bender_crittergripper/2026-07-14_bass16_bender"
OUTPUT_DIR <- "/Users/yjimenez/Library/CloudStorage/OneDrive-ProvidenceCollege/01_JimenezLab/02_ResearchHub/proj_crittergripper/figures/bass16_figures"
SUMMARY_PLOT_DIR <- file.path(OUTPUT_DIR, "summary_plots")

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("00_load_bender_flat.R")
src("02_deconvolve.R")
src("muscle_geometry.R")
src("fit_fv_fl.R")
src("03_analyze.R")
src("parse_trial_filename.R")
src("plot_force_vs_time.R")  # for .lowpass_filtfilt(), BASELINE_PAD_S, RELAXATION_WINDOW_S

CUTOFFS_HZ <- c(30, 15, 8)  # candidate smoothing options, widest to tightest

manifest <- parse_trial_directory(SOURCE_DIR)
iso_files <- dplyr::filter(manifest, .data$protocol == "isometric")
cli::cli_h1("Torque smoothing diagnostic: {nrow(iso_files)} isometric trial(s)")

timeseries_all <- list()
psd_all        <- list()

for (i in seq_len(nrow(iso_files))) {
  row <- iso_files[i, ]
  tid <- row$trial_id

  td <- load_bender_flat(row$fullpath, do_filter = TRUE, loadtorques = "x")
  comp <- deconvolve_bender(row$fullpath, hub_path = NULL, verbose = FALSE, components = TRUE)
  if (is.null(td) || is.null(comp)) { cli::cli_alert_warning("{tid}: load/deconvolve failed, skipped"); next }

  N <- min(nrow(td), length(comp$tau_corrected))
  td <- td[seq_len(N), , drop = FALSE]
  td$tau_raw_Nm       <- comp$tau_raw[seq_len(N)]
  td$tau_corrected_Nm <- comp$tau_corrected[seq_len(N)]

  geom <- .read_segmented_step_geometry(row$fullpath)
  steps <- dplyr::filter(geom$steps, is.finite(.data$stim_t0_s), is.finite(.data$stim_t1_s))

  for (j in seq_len(nrow(steps))) {
    s <- steps[j, ]
    step_rows <- td$step_number == s$step_number
    if (!any(step_rows)) next
    t_step   <- td$t.s[step_rows]
    raw_step <- td$tau_raw_Nm[step_rows]
    cor_step <- td$tau_corrected_Nm[step_rows]

    win <- t_step >= (s$stim_t0_s - BASELINE_PAD_S) & t_step <= (s$stim_t1_s + RELAXATION_WINDOW_S)
    if (!any(win, na.rm = TRUE)) next
    t_rel <- t_step[win] - s$stim_t0_s
    unit_id <- paste0(tid, "_step", s$step_number)

    base <- tibble::tibble(
      trial_id = tid, unit_id = unit_id, t_rel = t_rel,
      raw = raw_step[win], corrected_unsmoothed = cor_step[win]
    )
    for (fc in CUTOFFS_HZ) {
      base[[paste0("corrected_", fc, "Hz")]] <- .lowpass_filtfilt(cor_step, cutoff_hz = fc)[win]
    }
    timeseries_all[[unit_id]] <- base

    # PSD of the UNSMOOTHED corrected signal over the ACTIVE stim window only
    # (where the ~32 Hz wave was reported), one periodogram per step.
    act_win <- t_step >= s$stim_t0_s & t_step <= s$stim_t1_s
    y <- cor_step[act_win]
    y <- y[is.finite(y)]
    if (length(y) >= 30L) {
      y <- y - mean(y)
      sp <- tryCatch(stats::spec.pgram(y, taper = 0.1, detrend = TRUE, demean = TRUE, plot = FALSE),
                     error = function(e) NULL)
      if (!is.null(sp)) {
        psd_all[[paste0(unit_id, "_psd")]] <- tibble::tibble(
          trial_id = tid, unit_id = unit_id, freq_hz = sp$freq * 1000, power = sp$spec
        )
      }
    }
  }
}

ts_df  <- dplyr::bind_rows(timeseries_all)
psd_df <- dplyr::bind_rows(psd_all)


# =============================================================================
# Plot 1: time-domain comparison -- raw / unsmoothed-corrected / N Hz options
# =============================================================================

ts_long <- ts_df |>
  tidyr::pivot_longer(
    cols = c("raw", "corrected_unsmoothed", paste0("corrected_", CUTOFFS_HZ, "Hz")),
    names_to = "signal", values_to = "torque_Nm"
  ) |>
  dplyr::mutate(signal = factor(.data$signal,
    levels = c("raw", "corrected_unsmoothed", paste0("corrected_", CUTOFFS_HZ, "Hz")),
    labels = c("raw calibrated (pre-inertial-correction)", "inertia-corrected, UNsmoothed",
               paste0("inertia-corrected, ", CUTOFFS_HZ, " Hz low-pass"))
  ))

p_ts <- ggplot(ts_long, aes(x = .data$t_rel, y = .data$torque_Nm, color = .data$signal)) +
  geom_line(linewidth = 0.4, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  facet_wrap(~unit_id, scales = "free_y", ncol = 4) +
  scale_color_brewer(palette = "Dark2") +
  labs(title = "Torque smoothing diagnostic: raw vs. inertia-corrected vs. low-pass options",
       subtitle = "Every isometric step; PI picks ONE cutoff to apply uniformly pipeline-wide",
       x = "Time relative to stimulus onset (s)", y = "Torque (N*m)", color = "Signal") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(size = 7))

ggplot2::ggsave(file.path(SUMMARY_PLOT_DIR, "diag_torque_smoothing_timedomain.png"), p_ts,
                width = 16, height = 12, dpi = 150, limitsize = FALSE)


# =============================================================================
# Plot 2: pooled power spectrum of the UNSMOOTHED corrected signal (active
# window only) -- tests the ~32 Hz wave hypothesis directly.
# =============================================================================

if (nrow(psd_df) > 0L) {
  psd_mean <- psd_df |>
    dplyr::mutate(freq_bin = round(.data$freq_hz / 2) * 2) |>
    dplyr::group_by(.data$freq_bin) |>
    dplyr::summarise(power = mean(.data$power, na.rm = TRUE), n = dplyr::n(), .groups = "drop") |>
    dplyr::filter(.data$freq_bin <= 100, .data$freq_bin > 0)

  peak_freq <- psd_mean$freq_bin[which.max(psd_mean$power)]

  p_psd <- ggplot(psd_mean, aes(x = .data$freq_bin, y = .data$power)) +
    geom_line(color = "#1d4ed8", linewidth = 0.8) +
    geom_vline(xintercept = 32, linetype = "dashed", color = "red") +
    annotate("text", x = 32, y = max(psd_mean$power), label = "32 Hz", color = "red",
             hjust = -0.1, vjust = 1, size = 3.5) +
    scale_y_log10() +
    labs(title = "Pooled power spectrum: inertia-corrected torque, active stim window only",
         subtitle = sprintf("n = %d isometric steps pooled; largest peak observed at %.0f Hz (red dashed = 32 Hz reference)",
                            dplyr::n_distinct(psd_df$unit_id), peak_freq),
         x = "Frequency (Hz)", y = "Power spectral density (log scale)") +
    theme_bw(base_size = 12)

  ggplot2::ggsave(file.path(SUMMARY_PLOT_DIR, "diag_torque_smoothing_spectrum.png"), p_psd,
                  width = 9, height = 6, dpi = 150)
  cli::cli_alert_info("Pooled spectrum peak at {round(peak_freq)} Hz (n = {dplyr::n_distinct(psd_df$unit_id)} steps)")
} else {
  cli::cli_alert_warning("No PSD data computed -- spectrum plot skipped")
}

cli::cli_h1("Diagnostic complete: diag_torque_smoothing_timedomain.png, diag_torque_smoothing_spectrum.png")
