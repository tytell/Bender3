# diag_torque_smoothing.R
# Diagnostic (not a pipeline deliverable): compares RAW calibrated torque vs.
# inertia-corrected torque vs. several inertia-corrected + low-pass-filtered
# variants, for every isometric AND isovelocity step (bass16 only -- "a few
# examples," not a full 3-fish sweep), so the PI can visually pick a smoothing
# cutoff -- then that choice gets applied uniformly across the pipeline
# (build_segmented_force_timeseries()/build_isovelocity_force_timeseries()'s
# NOISE_FILTER_HZ, and .smooth_trace_display_only()'s DISPLAY_SMOOTH_HZ).
#
# EXTENDED 2026-07-22 (PI question: "might need to adjust per experiment type
# or angular acceleration") to ALSO cover isovelocity steps, not just
# isometric holds -- isometric has NO commanded motion during the stim window
# (bend angle is static), so any high-frequency content there is pure
# mechanical/electrical noise; isovelocity steps are ACTIVELY MOVING for the
# whole window, so the same high-frequency content could ALSO be a genuine
# inertial/motion-coupled transient (e.g. at ramp start/stop) rather than
# noise -- a single pipeline-wide cutoff picked from isometric-only data could
# be wrong for isovelocity. `protocol` is now carried through both plots
# (color-split PSD, protocol-prefixed facet labels) specifically so the two
# can be compared side by side, not pooled together. This diagnostic does NOT
# yet compute angular acceleration itself (bend-angle 2nd derivative) --
# that's a `angular_acceleration_check` follow-up, not done here.
#
# Also produces a pooled power-spectral-density plot of the UNSMOOTHED
# inertia-corrected signal during the active stim window, to directly test
# the ~32 Hz wave observed by eye in forceTime_isometric_legacy.png.
#
# EXTENDED AGAIN 2026-07-22 (PI: "run on bass17/bass18, include the
# broadband-power gap and examples of different filters for different
# samples") --
#   (1) parameterized BASS_ID/SOURCE_DIR/OUTPUT_DIR (env-var pattern, same as
#       export_snr_summary_figures.R) so this now runs per-specimen, not
#       hardcoded to bass16.
#   (2) the isometric-vs-isovelocity broadband power gap from the first
#       extension is now QUANTIFIED, not just visually apparent: the pooled
#       spectrum plot's subtitle reports the median isovelocity/isometric
#       power ratio across the shared 0-100 Hz band directly, per specimen.
#   (3) new Plot 3 (`diag_torque_smoothing_examples.png`): a SMALL, readable
#       set of example steps (2 isometric + 2 isovelocity, the median-power
#       step of each protocol so it's representative, not cherry-picked) at
#       large panel size, all filter variants overlaid -- the full facet grid
#       in Plot 1 has every step but is too dense to actually SEE how much
#       each cutoff changes a given trace; this is the "how different do the
#       filters actually look" companion.
#
# Run with:  Rscript R/diag_torque_smoothing.R
#   (optionally: BENDER3_BASS_ID=bass17 Rscript R/diag_torque_smoothing.R,
#   same env-var override pattern as export_snr_summary_figures.R)

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(ggplot2); library(cli); library(rhdf5)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")

# Defaults come from paths_config.R (single source of truth) -- see that
# file if the OneDrive folder layout ever moves again. BASS_ID/SOURCE_DIR/
# OUTPUT_DIR overridable via env vars (same pattern as
# export_snr_summary_figures.R) so this can be re-run per specimen without
# editing the file.
BASS_ID    <- Sys.getenv("BENDER3_BASS_ID", "bass16")
.raw_subfolder_for <- function(id) switch(id,
  bass16 = BASS16_RAW_SUBFOLDER, bass17 = BASS17_RAW_SUBFOLDER,
  bass18 = BASS18_RAW_SUBFOLDER, cli::cli_abort("Unknown BASS_ID {id}"))
SOURCE_DIR <- Sys.getenv("BENDER3_SOURCE_DIR", raw_source_dir(.raw_subfolder_for(BASS_ID)))
OUTPUT_DIR <- Sys.getenv("BENDER3_OUTPUT_DIR", figs_dir(BASS_ID))
# Flat, no summary_plots subfolder (2026-07-21, matches run_fv_fl_power_pipeline.R).
SUMMARY_PLOT_DIR <- OUTPUT_DIR

src("00_load_bender_flat.R")
src("02_deconvolve.R")
src("muscle_geometry.R")
src("fit_fv_fl.R")
src("03_analyze.R")
src("parse_trial_filename.R")
src("plot_force_vs_time.R")  # for .lowpass_filtfilt(), BASELINE_PAD_S, RELAXATION_WINDOW_S

CUTOFFS_HZ <- c(30, 15, 8)  # candidate smoothing options, widest to tightest

manifest <- parse_trial_directory(SOURCE_DIR)
stim_files <- dplyr::filter(manifest, .data$protocol %in% c("isometric", "isovelocity"))
cli::cli_h1("Torque smoothing diagnostic [{BASS_ID}]: {sum(stim_files$protocol=='isometric')} isometric + {sum(stim_files$protocol=='isovelocity')} isovelocity trial(s)")

timeseries_all <- list()
psd_all        <- list()

for (i in seq_len(nrow(stim_files))) {
  row <- stim_files[i, ]
  tid <- row$trial_id
  proto <- row$protocol

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
    unit_id <- paste0(proto, "_", tid, "_step", s$step_number)

    base <- tibble::tibble(
      trial_id = tid, unit_id = unit_id, protocol = proto, t_rel = t_rel,
      raw = raw_step[win], corrected_unsmoothed = cor_step[win]
    )
    for (fc in CUTOFFS_HZ) {
      base[[paste0("corrected_", fc, "Hz")]] <- .lowpass_filtfilt(cor_step, cutoff_hz = fc)[win]
    }
    timeseries_all[[unit_id]] <- base

    # PSD of the UNSMOOTHED corrected signal over the ACTIVE stim window only
    # (where the ~32 Hz wave was reported), one periodogram per step. Kept
    # per-protocol (not pooled) so isometric (static hold, no commanded
    # motion) and isovelocity (actively moving the whole window) can be
    # compared directly -- a peak that's isometric-only is noise; a peak
    # that's isovelocity-only (or much larger there) is more likely
    # motion/inertial-coupling, not measurement noise.
    act_win <- t_step >= s$stim_t0_s & t_step <= s$stim_t1_s
    y <- cor_step[act_win]
    y <- y[is.finite(y)]
    if (length(y) >= 30L) {
      y <- y - mean(y)
      sp <- tryCatch(stats::spec.pgram(y, taper = 0.1, detrend = TRUE, demean = TRUE, plot = FALSE),
                     error = function(e) NULL)
      if (!is.null(sp)) {
        psd_all[[paste0(unit_id, "_psd")]] <- tibble::tibble(
          trial_id = tid, unit_id = unit_id, protocol = proto, freq_hz = sp$freq * 1000, power = sp$spec
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
  labs(title = sprintf("Torque smoothing diagnostic [%s]: raw vs. inertia-corrected vs. low-pass options", BASS_ID),
       subtitle = "Every isometric + isovelocity step, facet label prefixed isometric_/isovelocity_; PI picks ONE cutoff to apply uniformly pipeline-wide (or confirms it needs to differ by protocol) -- see diag_torque_smoothing_examples.png for a legible close-up",
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
    dplyr::group_by(.data$protocol, .data$freq_bin) |>
    dplyr::summarise(power = mean(.data$power, na.rm = TRUE), n = dplyr::n(), .groups = "drop") |>
    dplyr::filter(.data$freq_bin <= 100, .data$freq_bin > 0)

  peak_by_protocol <- psd_mean |>
    dplyr::group_by(.data$protocol) |>
    dplyr::slice_max(.data$power, n = 1) |>
    dplyr::ungroup()
  for (r in seq_len(nrow(peak_by_protocol))) {
    cli::cli_alert_info("{peak_by_protocol$protocol[r]} spectrum peak at {round(peak_by_protocol$freq_bin[r])} Hz (n = {dplyr::n_distinct(dplyr::filter(psd_df, protocol==peak_by_protocol$protocol[r])$unit_id)} steps)")
  }

  # QUANTIFIED broadband power gap (PI-requested 2026-07-22, was visual-only
  # before): median isovelocity/isometric power ratio across every SHARED
  # frequency bin in the 0-100 Hz band, so "2-3 orders of magnitude" (the
  # earlier bass16-only visual read) becomes a number reported per specimen,
  # not eyeballed off a log-scale plot each time.
  gap_wide <- psd_mean |>
    dplyr::select(freq_bin, protocol, power) |>
    tidyr::pivot_wider(names_from = "protocol", values_from = "power") |>
    dplyr::filter(is.finite(.data$isometric), is.finite(.data$isovelocity), .data$isometric > 0)
  broadband_ratio <- if (nrow(gap_wide) > 0L) stats::median(gap_wide$isovelocity / gap_wide$isometric, na.rm = TRUE) else NA_real_
  cli::cli_alert_info("[{BASS_ID}] BROADBAND POWER GAP: median isovelocity/isometric power ratio across 0-100 Hz = {round(broadband_ratio, 1)}x")

  p_psd <- ggplot(psd_mean, aes(x = .data$freq_bin, y = .data$power, color = .data$protocol)) +
    geom_line(linewidth = 0.8) +
    geom_vline(xintercept = 32, linetype = "dashed", color = "red") +
    annotate("text", x = 32, y = max(psd_mean$power), label = "32 Hz", color = "red",
             hjust = -0.1, vjust = 1, size = 3.5) +
    scale_y_log10() +
    scale_color_manual(values = c(isometric = "#1d4ed8", isovelocity = "#ea580c"), name = "protocol") +
    labs(title = sprintf("Pooled power spectrum [%s]: inertia-corrected torque, active stim window only, BY PROTOCOL", BASS_ID),
         subtitle = sprintf("isometric n=%d steps (static hold, NO commanded motion) vs isovelocity n=%d steps (actively moving) -- BROADBAND POWER GAP: median isovelocity/isometric ratio across 0-100 Hz = %.1fx (a peak PRESENT ONLY in isovelocity, or this much larger there, is more likely motion/inertial-coupling than measurement noise); red dashed = 32 Hz reference",
                            dplyr::n_distinct(dplyr::filter(psd_df, protocol=="isometric")$unit_id),
                            dplyr::n_distinct(dplyr::filter(psd_df, protocol=="isovelocity")$unit_id),
                            broadband_ratio),
         x = "Frequency (Hz)", y = "Power spectral density (log scale)") +
    theme_bw(base_size = 12)

  ggplot2::ggsave(file.path(SUMMARY_PLOT_DIR, "diag_torque_smoothing_spectrum.png"), p_psd,
                  width = 9, height = 6.5, dpi = 150)
} else {
  cli::cli_alert_warning("No PSD data computed -- spectrum plot skipped")
}


# =============================================================================
# Plot 3: legible EXAMPLES -- 2 isometric + 2 isovelocity representative steps
# (the median-total-power step per protocol, not cherry-picked), large panels,
# all filter variants overlaid, so how much each cutoff actually changes a
# given trace is visible (Plot 1's full facet grid has every step but is too
# dense to read this off).
# =============================================================================

if (nrow(ts_df) > 0L) {
  step_power <- ts_df |>
    dplyr::group_by(.data$unit_id, .data$protocol) |>
    dplyr::summarise(total_var = stats::var(.data$corrected_unsmoothed, na.rm = TRUE), .groups = "drop")
  example_units <- step_power |>
    dplyr::group_by(.data$protocol) |>
    dplyr::slice_min(abs(.data$total_var - stats::median(.data$total_var, na.rm = TRUE)), n = 2, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::pull(.data$unit_id)

  ex_long <- ts_long |> dplyr::filter(.data$unit_id %in% example_units)
  p_ex <- ggplot(ex_long, aes(x = .data$t_rel, y = .data$torque_Nm, color = .data$signal)) +
    geom_line(linewidth = 0.7, alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    facet_wrap(~unit_id, scales = "free_y", ncol = 2) +
    scale_color_brewer(palette = "Dark2") +
    labs(title = sprintf("Torque smoothing EXAMPLES [%s]: median-power step per protocol, legible close-up", BASS_ID),
         subtitle = "2 isometric + 2 isovelocity steps (median total variance within their own protocol, so representative, not cherry-picked) -- compare how much each low-pass option changes the trace",
         x = "Time relative to stimulus onset (s)", y = "Torque (N*m)", color = "Signal") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", strip.text = element_text(size = 9))

  ggplot2::ggsave(file.path(SUMMARY_PLOT_DIR, "diag_torque_smoothing_examples.png"), p_ex,
                  width = 12, height = 9, dpi = 150)
} else {
  cli::cli_alert_warning("No timeseries data -- examples plot skipped")
}

cli::cli_h1("Diagnostic complete [{BASS_ID}]: diag_torque_smoothing_timedomain.png, diag_torque_smoothing_spectrum.png, diag_torque_smoothing_examples.png")
