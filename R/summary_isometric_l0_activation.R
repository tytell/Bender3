# summary_isometric_l0_activation.R
# PI request 2026-07-24: a normalised force-vs-time figure restricted to
# ISOMETRIC NEAR-ZERO (L0) contractions, pooled across every source in which
# they occur, plus activation/relaxation-time boxplots compared to the red
# (slow) muscle values of Coughlin & Carroll (2006, Comp. Biochem. Physiol.
# A 145:533-539, Fig. 2).
#
# Sources of L0 (commanded ~0 deg) stimulated contractions:
#   1. isometric protocol -- the near-zero-strain step(s).
#   2. isovelocity protocol -- the V=0 (operating_point == 0) step.
#   3. dynamic protocol -- the pre- and post-cycling L0 bookends
#      (detect_dynamic_l0_bookends(); left+right, before & after the sweep).
#
# TOP row: normalised muscle force (0..1, peak-normalised, baseline-subtracted,
#   sign-folded contraction-positive) vs. time relative to stim onset. Thin
#   lines = individual contractions (one per step/trial/bookend), coloured by
#   specimen; thick line = per-specimen mean trend on a common time grid.
#
# BOTTOM row: boxplots (+ raw points) of activation time and relaxation time
#   per specimen.
#   UPDATE 2026-07-24 (PI direction): "Activation/relaxation times should be
#   determined using the 10 to 90% (activation) and 90 to 10% (relaxation)
#   window. That's what Coughlin used." Definitions (changed from the prior
#   stim-onset/half-decay convention):
#     activation time = t(rise to 90% of peak) - t(rise to 10% of peak).
#     relaxation time = t(fall to 10% of peak, after the peak) -
#                        t(fall to 90% of peak, after the peak).
#   Both are pure RISE/FALL times relative to the contraction's OWN peak --
#   no longer anchored to stim onset/offset -- matching Coughlin's method.
#   Shaded band = Coughlin & Carroll (2006) red-muscle mean +/- SD read from
#   their Fig. 2 (TA ~78 ms, TR ~150 ms) -- a DIFFERENT paper (kinetics-
#   focused) than Coughlin (2000) (the power/work/tension paper whose 10-90/
#   90-10 CONVENTION we now match). C&C (2006)'s own TA/TR use THEIR paper's
#   stim-onset/half-decay definition, not this 10-90/90-10 one, so the
#   numeric comparison to the C&C band remains qualitative (order-of-
#   magnitude fast-vs-slow) for that reason ON TOP OF the existing tetanic-
#   vs-twitch caveat below.
#
# Run: Rscript R/summary_isometric_l0_activation.R  -> figs_summary/

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr); library(ggplot2); library(patchwork); library(cli); library(signal)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R"); src("00_load_bender_flat.R"); src("01_calibrate.R"); src("02_deconvolve.R")
src("muscle_geometry.R"); src("03_analyze.R"); src("parse_trial_filename.R")
src("plot_force_vs_time.R"); src("fit_fv_fl.R"); src("extract_dynamic_l0_bookends.R")

OUT_DIR      <- FIGS_SUMMARY_DIR
DATA_OUT_DIR <- file.path(.crittergripper_root(), "02_processed", "data_processed")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
SPECIMEN_SUBFOLDERS <- c(bass16 = BASS16_RAW_SUBFOLDER, bass17 = BASS17_RAW_SUBFOLDER, bass18 = BASS18_RAW_SUBFOLDER)
SPECIMEN_COLORS     <- c(bass16 = "#1b9e77", bass17 = "#d95f02", bass18 = "#7570b3")

NEARZERO_DEG   <- 0.5    # |operating_point| below this = "near-zero" isometric hold
PAD_S          <- 0.2
RELAX_S        <- 1.0
BOOKEND_RELAX_S <- 0.35   # bookends are ~54 ms twitches -> short decay tail

# Coughlin & Carroll (2006) Fig. 2 values (approx read-off, ms).
# RED = slow swimming muscle; WHITE/FAST = the fast feeding muscles
# (sternohyoideus + epaxial), which are the white/fast counterpart.
COUGHLIN <- list(activation = c(mean = 78, sd = 22), relaxation = c(mean = 150, sd = 25))
COUGHLIN_WHITE <- list(activation = c(lo = 10, hi = 20, mid = 15),
                       relaxation = c(lo = 28, hi = 45, mid = 37))

.deconv_td <- function(f) {
  td <- tryCatch(load_bender_flat(f, do_filter = TRUE, loadtorques = "x"), error = function(e) NULL)
  if (is.null(td)) return(NULL)
  tau <- tryCatch(deconvolve_bender(f, hub_path = NULL, verbose = FALSE), error = function(e) NULL)
  if (is.null(tau)) return(NULL)
  N <- min(nrow(td), length(tau)); td <- td[seq_len(N), , drop = FALSE]
  td$torque_inertia_corrected_Nm <- tau[seq_len(N)]
  attr(td, "Filename") <- f; td
}

# force timeseries for dynamic L0 bookends (mirrors build_segmented_force_timeseries
# columns; sign & baseline derived directly from the burst since bookends carry
# no force_sign/passive_force_Nm metadata).
.bookend_timeseries <- function(td, bookends, trial_id) {
  if (nrow(bookends) == 0L) return(tibble::tibble())
  bookends <- dplyr::arrange(bookends, .data$stim_t0_s)
  # cap each burst's window BEFORE the NEXT stim event of ANY kind (the other
  # bookend of the pair OR the first cycling burst, which can start <1s after
  # the pre-right bookend) -- otherwise the full RELAX_S tail swallows it and
  # creates a spurious late peak / wrong normalisation.
  all_ev <- tryCatch(.detect_stim_events(td), error = function(e) NULL)
  onsets <- if (!is.null(all_ev) && nrow(all_ev)) sort(all_ev$stim_t0_s) else numeric(0)
  purrr::map_dfr(seq_len(nrow(bookends)), function(i) {
    b <- bookends[i, ]
    base_win <- td$t.s >= b$t_pre_baseline_start_s & td$t.s <= b$t_pre_baseline_end_s
    act_win  <- td$t.s >= b$stim_t0_s & td$t.s <= b$stim_t1_s
    if (!any(base_win) || !any(act_win)) return(NULL)
    baseline <- mean(td$torque_inertia_corrected_Nm[base_win], na.rm = TRUE)
    sgn <- sign(mean(td$torque_inertia_corrected_Nm[act_win], na.rm = TRUE) - baseline)
    if (!is.finite(sgn) || sgn == 0) sgn <- 1
    nxt <- onsets[onsets > b$stim_t1_s + 1e-6]
    next_onset <- if (length(nxt)) min(nxt) else Inf
    # bookends are ~54 ms twitches: a short post-stim tail fully captures their
    # rise+decay; a long tail only invites pre-cycling motor-motion artifacts
    # (whose torque can exceed the small twitch peak). Cap at 0.35 s AND before
    # the next stim event of any kind.
    win_end <- min(b$stim_t1_s + BOOKEND_RELAX_S, next_onset - 0.05)
    win <- td$t.s >= (b$stim_t0_s - PAD_S) & td$t.s <= win_end
    if (!any(win)) return(NULL)
    tibble::tibble(
      trial_id = trial_id, unit_id = paste0(trial_id, "_bookend", i),
      muscle_side = b$muscle_side, source = "dynamic L0 bookend",
      t_rel = td$t.s[win] - b$stim_t0_s,
      muscle_force_Nm = sgn * (td$torque_inertia_corrected_Nm[win] - baseline),
      stim_duration_s = b$stim_t1_s - b$stim_t0_s)
  })
}

cli::cli_h1("Collecting isometric near-zero (L0) contractions")
units <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  manifest <- parse_trial_directory(raw_source_dir(subfolder))
  purrr::pmap_dfr(list(manifest$fullpath, manifest$protocol), function(f, proto) {
    trial_id <- tools::file_path_sans_ext(basename(f))
    if (proto %in% c("isometric", "isovelocity")) {
      td <- .deconv_td(f); if (is.null(td)) return(tibble())
      res <- tryCatch(if (proto == "isometric") analyze_isometric(td) else analyze_isovelocity(td),
                      error = function(e) NULL)
      if (is.null(res) || is.null(res$step_summary) || nrow(res$step_summary) == 0L) return(tibble())
      tdA <- res$td
      if (!("torque_inertia_corrected_Nm" %in% names(tdA)) && nrow(tdA) == nrow(td))
        tdA$torque_inertia_corrected_Nm <- td$torque_inertia_corrected_Nm
      ss0 <- dplyr::filter(res$step_summary, abs(.data$operating_point) < NEARZERO_DEG)
      if (nrow(ss0) == 0L) return(tibble())
      out <- build_segmented_force_timeseries(tdA, ss0, trial_id)
      if (nrow(out) == 0L) return(tibble())
      out$source <- if (proto == "isometric") "isometric L0 step" else "isovelocity V=0 step"
      out$specimen <- specimen; out
    } else if (proto == "dynamic") {
      td <- .deconv_td(f); if (is.null(td)) return(tibble())
      be <- tryCatch(detect_dynamic_l0_bookends(td), error = function(e) NULL)
      if (is.null(be) || nrow(be) == 0L) return(tibble())
      out <- .bookend_timeseries(td, be, trial_id)
      if (nrow(out) == 0L) return(tibble())
      out$specimen <- specimen; out
    } else tibble()
  })
})
cli::cli_alert_info("Collected {dplyr::n_distinct(units$unit_id)} L0 contractions across {dplyr::n_distinct(units$specimen)} specimens")

# ---- per-unit: filter, normalise, compute activation/relaxation times ----
.process_unit <- function(u) {
  u <- dplyr::arrange(u, .data$t_rel)
  f <- .lowpass_filtfilt(u$muscle_force_Nm, DISPLAY_SMOOTH_HZ, fs_hz = 1000)
  t   <- u$t_rel
  dur <- u$stim_duration_s[1L]
  # peak is the contraction peak (rise+plateau of THIS burst), not any late
  # residual -- search [0, stim_offset + 0.15s].
  cwin <- t >= 0 & t <= (dur + 0.15)
  if (!any(cwin)) return(NULL)
  cwin_idx <- which(cwin)
  peak_rel_i <- which.max(f[cwin_idx])
  peak <- f[cwin_idx[peak_rel_i]]
  if (!is.finite(peak) || peak <= 0) return(NULL)
  fn <- f / peak
  peak_t <- t[cwin_idx[peak_rel_i]]

  # ACTIVATION (PI direction, 2026-07-24 -- matches Coughlin (2000)'s own
  # convention): rise time from the FIRST crossing of 10% peak to the FIRST
  # SUBSEQUENT crossing of 90% peak, searched within the rise phase only
  # (contraction window up to the peak) -- no longer anchored to stim onset.
  rise_win <- cwin & t <= peak_t
  idx10r <- which(rise_win & fn >= 0.10)
  t10_rise <- if (length(idx10r)) t[idx10r[1L]] else NA_real_
  ta_ms <- NA_real_; t90_rise <- NA_real_
  if (is.finite(t10_rise)) {
    idx90r <- which(rise_win & t >= t10_rise & fn >= 0.90)
    if (length(idx90r)) { t90_rise <- t[idx90r[1L]]; ta_ms <- (t90_rise - t10_rise) * 1000 }
  }

  # RELAXATION (PI direction, 2026-07-24): fall time from the FIRST
  # post-peak crossing below 90% peak to the FIRST SUBSEQUENT crossing below
  # 10% peak -- no longer anchored to stim offset. Searched over the whole
  # post-peak trace available for this unit (already time-bounded upstream
  # by BOOKEND_RELAX_S / RELAX_S, so no extra cap needed here).
  fall_win <- t > peak_t
  idx90f <- which(fall_win & fn <= 0.90)
  tr_ms <- NA_real_; t90_fall <- NA_real_; t10_fall <- NA_real_
  if (length(idx90f)) {
    t90_fall <- t[idx90f[1L]]
    idx10f <- which(fall_win & t >= t90_fall & fn <= 0.10)
    if (length(idx10f)) { t10_fall <- t[idx10f[1L]]; tr_ms <- (t10_fall - t90_fall) * 1000 }
  }

  u$force_norm <- fn
  attr(u, "ta_ms") <- ta_ms; attr(u, "tr_ms") <- tr_ms
  attr(u, "t10_rise") <- t10_rise; attr(u, "t90_rise") <- t90_rise
  attr(u, "t90_fall") <- t90_fall; attr(u, "t10_fall") <- t10_fall
  u
}

by_unit <- split(units, units$unit_id)
proc <- purrr::compact(lapply(by_unit, .process_unit))
traces <- dplyr::bind_rows(proc)

# DIAGNOSTIC: which units carry high force late (t_rel > 0.6s)? -> late-peak source
late_dbg <- traces |>
  dplyr::group_by(.data$unit_id, .data$specimen, .data$source) |>
  dplyr::summarise(late_hi = any(.data$force_norm > 0.8 & .data$t_rel > 0.6, na.rm = TRUE), .groups = "drop") |>
  dplyr::filter(.data$late_hi)
cli::cli_h2("units with force_norm>0.8 at t_rel>0.6 (late-peak diagnostic)")
print(dplyr::count(late_dbg, .data$specimen, .data$source))
times <- purrr::imap_dfr(proc, function(u, id) tibble::tibble(
  unit_id = id, specimen = u$specimen[1L], source = u$source[1L],
  stim_duration_s = u$stim_duration_s[1L],
  activation_ms = attr(u, "ta_ms"), relaxation_ms = attr(u, "tr_ms"),
  t10_rise_s = attr(u, "t10_rise"), t90_rise_s = attr(u, "t90_rise"),
  t90_fall_s = attr(u, "t90_fall"), t10_fall_s = attr(u, "t10_fall")))
cli::cli_h2("stim duration (s) and count by source")
print(times |> dplyr::group_by(.data$source) |>
  dplyr::summarise(n = dplyr::n(), dur_med = round(median(.data$stim_duration_s), 3),
                    dur_min = round(min(.data$stim_duration_s), 3),
                    dur_max = round(max(.data$stim_duration_s), 3), .groups = "drop"))
print(times |> dplyr::count(.data$specimen, .data$source))
write.csv(times, file.path(DATA_OUT_DIR, "isometric_l0_activation_times.csv"), row.names = FALSE)
# Raw per-sample traces (unit_id, specimen, source, t_rel, force_norm) so
# downstream figures (e.g. summary_coughlin2000_bass_comparison.R) can reuse
# the same collected/processed bookend twitches without recomputing them.
write.csv(dplyr::select(traces, "unit_id", "specimen", "source", "t_rel", "force_norm"),
          file.path(DATA_OUT_DIR, "isometric_l0_activation_traces.csv"), row.names = FALSE)
cli::cli_alert_info("activation/relaxation available for n={sum(is.finite(times$activation_ms))}/{sum(is.finite(times$relaxation_ms))} contractions")

# ---- box panel with BOTH Coughlin bands (red = slow, gray = white/fast) ----
.red_band <- function(c) data.frame(lo = c[["mean"]] - c[["sd"]], hi = c[["mean"]] + c[["sd"]], mid = c[["mean"]])
.box_panel <- function(times_in, col, red, white, title, ylab) {
  d <- times_in[is.finite(times_in[[col]]), ]
  rb <- .red_band(red)
  ggplot(d, aes(x = specimen, y = .data[[col]])) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = white[["lo"]], ymax = white[["hi"]], fill = "grey45", alpha = 0.18) +
    geom_hline(yintercept = white[["mid"]], color = "grey35", linetype = "dashed", linewidth = 0.5) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = rb$lo, ymax = rb$hi, fill = "#b30000", alpha = 0.12) +
    geom_hline(yintercept = rb$mid, color = "#b30000", linetype = "dashed", linewidth = 0.5) +
    geom_boxplot(aes(color = specimen), fill = NA, outlier.shape = NA, width = 0.55) +
    geom_jitter(aes(color = specimen), width = 0.12, height = 0, size = 1.8, alpha = 0.7) +
    annotate("text", x = 0.6, y = red[["mean"]], label = "C&C red (slow)", hjust = 0, vjust = -0.5, size = 2.7, color = "#b30000") +
    annotate("text", x = 0.6, y = white[["mid"]], label = "C&C white (fast)", hjust = 0, vjust = -0.5, size = 2.7, color = "grey30") +
    scale_color_manual(values = SPECIMEN_COLORS, guide = "none") +
    labs(title = title, x = NULL, y = ylab) +
    theme_bw(base_size = 12)
}

# ---- build + save one full kinetics figure from a traces/times subset ----
build_kinetics_fig <- function(traces_in, times_in, fout, top_source_note) {
  grid <- seq(-PAD_S, 1.2, by = 0.005)
  mean_trend <- purrr::map_dfr(split(traces_in, traces_in$specimen), function(df) {
    units_here <- dplyr::n_distinct(df$unit_id)
    m <- purrr::map(split(df, df$unit_id), function(u)
      stats::approx(u$t_rel, u$force_norm, xout = grid, rule = 1)$y)
    M <- do.call(cbind, m)
    tibble::tibble(specimen = df$specimen[1L], t_rel = grid,
                   force_norm = rowMeans(M, na.rm = TRUE),
                   n = rowSums(is.finite(M)),
                   min_support = pmax(10, 0.50 * units_here))
  }) |> dplyr::filter(.data$n >= .data$min_support)

  pTop <- ggplot() +
    geom_hline(yintercept = c(0, 1), linetype = "dotted", color = "grey85") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_line(data = traces_in, aes(x = t_rel, y = force_norm, group = unit_id, color = specimen),
              alpha = 0.20, linewidth = 0.3) +
    geom_line(data = mean_trend, aes(x = t_rel, y = force_norm, color = specimen), linewidth = 1.4) +
    scale_color_manual(values = SPECIMEN_COLORS, name = "individual") +
    coord_cartesian(xlim = c(-PAD_S, 1.2), ylim = c(-0.15, 1.1)) +
    labs(title = "A. Isometric near-zero (L0) contractions: normalised force vs. time",
         subtitle = sprintf("Peak-normalised, baseline-subtracted, sign-folded. Thin = individual contractions (%s); thick = per-specimen mean (>=50%% support). n=%d contractions.", top_source_note, dplyr::n_distinct(traces_in$unit_id)),
         x = "Time relative to stim onset (s)", y = "Normalised force (F / peak)") +
    theme_bw(base_size = 12)

  pTA <- .box_panel(times_in, "activation_ms", COUGHLIN$activation, COUGHLIN_WHITE$activation,
                    "B. Activation time (10% to 90% of peak)", "Activation time (ms)")
  pTR <- .box_panel(times_in, "relaxation_ms", COUGHLIN$relaxation, COUGHLIN_WHITE$relaxation,
                    "C. Relaxation time (90% to 10% of peak)", "Relaxation time (ms)")

  fig <- pTop / (pTA | pTR) +
    patchwork::plot_layout(heights = c(1.25, 1)) +
    patchwork::plot_annotation(
      title = "Isometric L0 contractile kinetics across individuals vs. Coughlin & Carroll (2006) red (slow) + white (fast) muscle",
      subtitle = "Ours = 10-90% rise / 90-10% fall times (Coughlin (2000) convention, PI-directed 2026-07-24); C&C (2006) band uses THEIR OWN stim-onset/half-decay\ndefinition on single twitches, not this one -- comparison is qualitative (fast-vs-slow), not a matched-definition benchmark.",
      theme = theme(plot.title = element_text(face = "bold", size = 12)))
  ggplot2::ggsave(fout, fig, width = 11, height = 9, dpi = 150)
  cli::cli_alert_success("Saved {fout}")
}

# (1) all sources
build_kinetics_fig(traces, times,
  file.path(OUT_DIR, "isometric_L0_activation_kinetics.png"),
  "isometric L0 steps + isovelocity V=0 + dynamic pre/post L0 bookends")

# (2) dynamic bookend twitches ONLY -- so the force traces and the kinetics
#     boxplots are a single comparable contraction type (matches the
#     type-controlled early-vs-later figure).
be_ids <- unique(times$unit_id[times$source == "dynamic L0 bookend"])
build_kinetics_fig(dplyr::filter(traces, .data$unit_id %in% be_ids),
  dplyr::filter(times, .data$source == "dynamic L0 bookend"),
  file.path(OUT_DIR, "isometric_L0_activation_kinetics_bookendsOnly.png"),
  "dynamic pre/post L0 bookend twitches ONLY")
