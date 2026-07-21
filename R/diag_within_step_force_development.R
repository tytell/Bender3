# diag_within_step_force_development.R
# READ-ONLY DIAGNOSTIC (2026-07-16, PI-requested). Does NOT modify any
# pipeline analysis: it sources the existing functions and only reads data +
# writes PNGs. Purpose: within-step force development for bass17 -- does
# muscle force reach a plateau during the 300 ms stimulus or is it still
# rising at stim-off, and does time-to-plateau differ by length (isometric)
# / velocity (isovelocity)? Also marks WHERE the pipeline samples its
# "active" value (mean over [stim_t0, stim_t1 + DEACTIVATION_WINDOW_S]).
#
# Run: Rscript R/diag_within_step_force_development.R

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(stringr)
  library(ggplot2); library(cli); library(patchwork); library(signal)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")

# Defaults come from paths_config.R (single source of truth) -- see that
# file if the OneDrive folder layout ever moves again.
SRC_DIR <- raw_source_dir(BASS17_RAW_SUBFOLDER)
OUT_DIR <- FIGS_SUMMARY_DIR

src("00_load_bender_flat.R")
src("01_calibrate.R")
src("02_deconvolve.R")
src("muscle_geometry.R")
src("fit_fv_fl.R")
src("03_analyze.R")
src("parse_trial_filename.R")
src("plot_force_vs_time.R")  # build_segmented_force_timeseries, build_isovelocity_force_timeseries,
                             # .smooth_trace_display_only, DEACTIVATION not here (03_analyze)

`%||%` <- function(x, y) if (is.null(x)) y else x
wrap_sub <- function(s, w = 88) paste(strwrap(s, width = w), collapse = "\n")
DEACTIVATION_WINDOW_S <- 0.5  # matches build_segmented_step_summary default -- the active-mean tail
BASE_PAD_S  <- 0.2            # matches build_segmented_force_timeseries()
RELAX_S     <- 1.0            # matches build_segmented_force_timeseries()
LEVEL_COLORS <- c("#2c7fb8", "#addd8e", "#d95f0e")  # low / mid / high

load_one <- function(f) {
  td <- load_bender_flat(f, do_filter = TRUE, loadtorques = "x")
  tau <- deconvolve_bender(f, hub_path = NULL, verbose = FALSE)
  N <- min(nrow(td), length(tau)); td <- td[seq_len(N), , drop = FALSE]
  td$torque_inertia_corrected_Nm <- tau[seq_len(N)]
  attr(td, "Filename") <- f
  td
}

#' Raw inertia-corrected torque windowed to [stim_t0 - BASE_PAD, stim_t1 + RELAX]
#' for one step, with t_rel = t.s - stim_t0. Read-only slice of the analyzed td.
raw_step_window <- function(td, step_row) {
  sn <- step_row$step_number
  rows <- which(td$step_number == sn)
  if (length(rows) == 0L) return(tibble::tibble())
  sub <- td[rows, ]
  t_rel <- sub$t.s - step_row$stim_t0_s
  keep <- t_rel >= -BASE_PAD_S & t_rel <= (step_row$stim_t1_s - step_row$stim_t0_s + RELAX_S)
  tibble::tibble(step_number = sn, t_rel = t_rel[keep],
                 torque_Nm = sub$torque_inertia_corrected_Nm[keep])
}

#' Peak / stim-offset / time-to-plateau report on a baseline-subtracted,
#' sign-corrected muscle-force development trace (positive = active pull).
#' All timing relative to stim onset; stim-off at stim_dur_s.
step_report <- function(t_rel, force_Nm, stim_dur_s) {
  y <- .smooth_trace_display_only(force_Nm)  # light 8 Hz, display+metric consistency
  in_stim <- t_rel >= 0 & t_rel <= stim_dur_s
  full <- is.finite(y)
  peak_full     <- if (any(full)) max(y[full]) else NA_real_
  peak_in_stim  <- if (any(in_stim & is.finite(y))) max(y[in_stim & is.finite(y)]) else NA_real_
  off_i         <- which.min(abs(t_rel - stim_dur_s))
  force_at_off  <- if (length(off_i)) y[off_i] else NA_real_
  # time-to-95%-of-in-stim-peak (proxy for time-to-plateau)
  t_plateau <- NA_real_
  if (is.finite(peak_in_stim) && peak_in_stim > 0) {
    hit <- which(in_stim & is.finite(y) & y >= 0.95 * peak_in_stim)
    if (length(hit)) t_plateau <- min(t_rel[hit])
  }
  list(peak_full = peak_full, peak_in_stim = peak_in_stim,
       force_at_off = force_at_off, t_to_plateau_s = t_plateau,
       still_rising = is.finite(force_at_off) && is.finite(peak_in_stim) &&
         force_at_off >= 0.98 * peak_in_stim)
}

# =============================================================================
# ISOMETRIC (bass17 has a single isometric trial -> "strongest" is trivial)
# =============================================================================

manifest <- parse_trial_directory(SRC_DIR)
iso_files <- manifest$fullpath[manifest$protocol == "isometric"]
cli::cli_h1("Isometric: {length(iso_files)} trial(s) -> using {basename(iso_files[1])}")

iso_f  <- iso_files[1]
iso_td <- load_one(iso_f)
iso    <- analyze_isometric(iso_td, filename = iso_f)
iso_ss <- iso$step_summary

# Block 1 = left muscle / left bend (steps 1-8). Steps 1 (0%), 5 (~16%), 8 (~28%).
iso_pick_steps <- c(1L, 5L, 8L)
iso_sel <- dplyr::filter(iso_ss, step_number %in% iso_pick_steps) |>
  dplyr::arrange(abs(shortening_strain_pct))
iso_sel$level <- factor(round(iso_sel$shortening_strain_pct), levels = round(sort(iso_sel$shortening_strain_pct)))
iso_sel$level_lab <- sprintf("%.0f%% strain (step %d)", iso_sel$shortening_strain_pct, iso_sel$step_number)

# Subtracted (active) development trace via the pipeline's own function.
iso_ts <- build_segmented_force_timeseries(iso$td, iso_ss, trial_id = "bass17_iso")
iso_ts$step_number <- as.integer(stringr::str_extract(iso_ts$unit_id, "(?<=_step)\\d+"))
iso_ts_sel <- dplyr::filter(iso_ts, step_number %in% iso_pick_steps) |>
  dplyr::left_join(dplyr::select(iso_sel, step_number, level_lab, shortening_strain_pct), by = "step_number")

# Raw torque windows + per-step markers (passive baseline, active-window mean).
iso_raw <- purrr::map_dfr(iso_pick_steps, function(sn) raw_step_window(iso_td, iso_ss[iso_ss$step_number == sn, ][1, ]))
iso_raw <- dplyr::left_join(iso_raw, dplyr::select(iso_sel, step_number, level_lab), by = "step_number")

iso_markers <- purrr::map_dfr(iso_pick_steps, function(sn) {
  s <- iso_ss[iso_ss$step_number == sn, ][1, ]
  dur <- s$stim_t1_s - s$stim_t0_s
  # active-mean window [0, dur + deactivation]; active mean = mean raw torque there
  rw <- raw_step_window(iso_td, s)
  active_win <- rw$t_rel >= 0 & rw$t_rel <= (dur + DEACTIVATION_WINDOW_S)
  tibble::tibble(
    step_number = sn, stim_off = dur,
    active_win_end = dur + DEACTIVATION_WINDOW_S,
    passive_mean = s$passive_force_Nm,
    active_mean  = mean(rw$torque_Nm[active_win], na.rm = TRUE),
    level_lab = iso_sel$level_lab[iso_sel$step_number == sn]
  )
})

cli::cli_h2("Isometric within-step report (subtracted active muscle force, N*m)")
iso_rep <- purrr::map_dfr(iso_pick_steps, function(sn) {
  d <- dplyr::filter(iso_ts_sel, step_number == sn)
  dur <- iso_markers$stim_off[iso_markers$step_number == sn]
  r <- step_report(d$t_rel, d$muscle_force_Nm, dur)
  tibble::tibble(step = sn,
                 strain_pct = round(iso_sel$shortening_strain_pct[iso_sel$step_number == sn], 1),
                 peak_active_Nm = signif(r$peak_in_stim, 3),
                 force_at_stimoff_Nm = signif(r$force_at_off, 3),
                 time_to_plateau_s = signif(r$t_to_plateau_s, 3),
                 still_rising_at_off = r$still_rising)
})
print(iso_rep)


# =============================================================================
# ISOVELOCITY (pick strongest of 3 trials by peak |muscle_force_Nm|)
# =============================================================================

isov_files <- manifest$fullpath[manifest$protocol == "isovelocity"]
cli::cli_h1("Isovelocity: {length(isov_files)} trial(s) -- selecting strongest by peak |torque|")

isov_analyzed <- purrr::map(isov_files, function(f) {
  td <- load_one(f)
  list(f = f, td = td, res = analyze_isovelocity(td, filename = f))
})
isov_strength <- purrr::map_dbl(isov_analyzed, function(a) {
  max(abs(a$res$step_summary$muscle_force_Nm), na.rm = TRUE)
})
for (i in seq_along(isov_analyzed)) {
  cli::cli_alert("{basename(isov_analyzed[[i]]$f)}: peak |muscle force| = {signif(isov_strength[i],3)} N*m")
}
pick <- which.max(isov_strength)
cli::cli_alert_success("Selected {basename(isov_analyzed[[pick]]$f)} (strongest)")

isov_f  <- isov_analyzed[[pick]]$f
isov_td <- isov_analyzed[[pick]]$td
isov    <- isov_analyzed[[pick]]$res
isov_ss <- isov$step_summary

# Block 1 = left/left (steps 1-4): velocities 0, ~95, ~189, ~284 %/s.
# Span with steps 1 (0), 3 (~189), 4 (~284 %/s).
isov_pick_steps <- c(1L, 3L, 4L)
isov_sel <- dplyr::filter(isov_ss, step_number %in% isov_pick_steps) |>
  dplyr::arrange(abs(shortening_strain_pct))
isov_sel$level_lab <- sprintf("%.0f %%/s (step %d)", isov_sel$shortening_strain_pct, isov_sel$step_number)

# Velocity-matched subtracted development trace (pipeline function). May be
# empty for a step lacking a no-stim velocity match -- handled below.
isov_ts <- tryCatch(build_isovelocity_force_timeseries(isov$td, isov_ss, trial_id = "bass17_isov"),
                    error = function(e) tibble::tibble())
if (nrow(isov_ts) > 0) {
  isov_ts$step_number <- as.integer(stringr::str_extract(isov_ts$unit_id, "(?<=_step)\\d+"))
  isov_ts_sel <- dplyr::filter(isov_ts, step_number %in% isov_pick_steps) |>
    dplyr::left_join(dplyr::select(isov_sel, step_number, level_lab, shortening_strain_pct), by = "step_number")
} else {
  isov_ts_sel <- tibble::tibble()
}

isov_raw <- purrr::map_dfr(isov_pick_steps, function(sn) raw_step_window(isov_td, isov_ss[isov_ss$step_number == sn, ][1, ]))
isov_raw <- dplyr::left_join(isov_raw, dplyr::select(isov_sel, step_number, level_lab), by = "step_number")

isov_markers <- purrr::map_dfr(isov_pick_steps, function(sn) {
  s <- isov_ss[isov_ss$step_number == sn, ][1, ]
  dur <- s$stim_t1_s - s$stim_t0_s
  rw <- raw_step_window(isov_td, s)
  active_win <- rw$t_rel >= 0 & rw$t_rel <= (dur + DEACTIVATION_WINDOW_S)
  tibble::tibble(step_number = sn, stim_off = dur, active_win_end = dur + DEACTIVATION_WINDOW_S,
                 active_mean = mean(rw$torque_Nm[active_win], na.rm = TRUE),
                 level_lab = isov_sel$level_lab[isov_sel$step_number == sn])
})

if (nrow(isov_ts_sel) > 0) {
  cli::cli_h2("Isovelocity within-step report (velocity-matched active force, N*m)")
  isov_rep <- purrr::map_dfr(isov_pick_steps, function(sn) {
    d <- dplyr::filter(isov_ts_sel, step_number == sn)
    if (nrow(d) == 0) return(tibble::tibble(step = sn, note = "no velocity-matched baseline"))
    dur <- isov_markers$stim_off[isov_markers$step_number == sn]
    r <- step_report(d$t_rel, d$muscle_force_Nm, dur)
    tibble::tibble(step = sn,
                   strain_rate_pct_s = round(isov_sel$shortening_strain_pct[isov_sel$step_number == sn], 1),
                   peak_active_Nm = signif(r$peak_in_stim, 3),
                   force_at_stimoff_Nm = signif(r$force_at_off, 3),
                   time_to_plateau_s = signif(r$t_to_plateau_s, 3),
                   still_rising_at_off = r$still_rising)
  })
  print(isov_rep)
} else {
  cli::cli_alert_warning("No velocity-matched subtracted traces available -- isovelocity report from RAW torque only")
}


# =============================================================================
# Plot builders
# =============================================================================

stim_vlines <- function(stim_dur) list(
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30"),
  geom_vline(xintercept = stim_dur, linetype = "dashed", color = "grey30")
)

# --- Overlay: subtracted active force, common y, colored by level ---
build_overlay <- function(ts_sel, stim_dur, title, sub, ylab) {
  ts_sel <- ts_sel |> dplyr::arrange(unit_id, t_rel) |>
    dplyr::group_by(unit_id) |>
    dplyr::mutate(y = .smooth_trace_display_only(muscle_force_Nm)) |> dplyr::ungroup()
  ggplot(ts_sel, aes(t_rel, y, color = level_lab)) +
    annotate("rect", xmin = 0, xmax = stim_dur, ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.10) +
    annotate("rect", xmin = stim_dur, xmax = stim_dur + DEACTIVATION_WINDOW_S, ymin = -Inf, ymax = Inf, fill = "grey50", alpha = 0.10) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    stim_vlines(stim_dur) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = LEVEL_COLORS, name = NULL) +
    labs(title = title, subtitle = wrap_sub(sub), x = "Time relative to stim onset (s)", y = ylab) +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
}

# --- Small multiples: RAW inertia-corrected torque + extraction markers ---
build_facets <- function(raw_df, markers, stim_dur, title, sub, show_passive = TRUE) {
  raw_df <- raw_df |> dplyr::arrange(step_number, t_rel) |>
    dplyr::group_by(step_number) |>
    dplyr::mutate(y = .smooth_trace_display_only(torque_Nm)) |> dplyr::ungroup()
  p <- ggplot(raw_df, aes(t_rel, y)) +
    annotate("rect", xmin = 0, xmax = stim_dur, ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.10) +
    geom_rect(data = markers, aes(xmin = 0, xmax = active_win_end, ymin = -Inf, ymax = Inf),
              inherit.aes = FALSE, fill = "grey50", alpha = 0.08) +
    geom_line(aes(y = torque_Nm), color = "grey75", linewidth = 0.3) +
    geom_line(color = "black", linewidth = 0.7) +
    stim_vlines(stim_dur) +
    geom_hline(data = markers, aes(yintercept = active_mean), color = "#d95f0e", linetype = "solid", linewidth = 0.7) +
    facet_wrap(~ level_lab, scales = "free_y", ncol = 1) +
    labs(title = title, subtitle = wrap_sub(sub),
         x = "Time relative to stim onset (s)", y = "Inertia-corrected torque (N*m)") +
    theme_bw(base_size = 11)
  if (show_passive && "passive_mean" %in% names(markers)) {
    p <- p + geom_hline(data = markers, aes(yintercept = passive_mean), color = "#1d4ed8", linetype = "dotted", linewidth = 0.7)
  }
  p
}

iso_dur  <- iso_markers$stim_off[1]
isov_dur <- isov_markers$stim_off[1]

p_iso_overlay <- build_overlay(
  iso_ts_sel, iso_dur,
  "bass17 isometric: within-step active force development (baseline-subtracted)",
  "Common y | orange = 300 ms stim, grey = 0.5 s deactivation tail | colored by held strain",
  "Active muscle force (N*m, shortening +)")

p_iso_facets <- build_facets(
  iso_raw, iso_markers, iso_dur,
  "bass17 isometric: raw inertia-corrected torque + extraction markers",
  "Free y | orange line = active-window mean (pipeline value) | blue dotted = pre-stim passive baseline | grey band = [stim_t0, stim_t1+0.5s] averaging window",
  show_passive = TRUE)

ggplot2::ggsave(file.path(OUT_DIR, "diag_within_step_force_isometric_overlay.png"),
                p_iso_overlay, width = 9, height = 5.5, dpi = 150)
ggplot2::ggsave(file.path(OUT_DIR, "diag_within_step_force_isometric_facets.png"),
                p_iso_facets, width = 8, height = 9, dpi = 150)

# Isovelocity: overlay from subtracted trace if present, else raw overlay.
if (nrow(isov_ts_sel) > 0) {
  p_isov_overlay <- build_overlay(
    isov_ts_sel, isov_dur,
    "bass17 isovelocity: within-step active force (velocity-matched subtraction)",
    "Common y | orange = 300 ms stim, grey = 0.5 s deactivation tail | colored by strain rate",
    "Active muscle force (N*m)")
} else {
  isov_raw2 <- isov_raw |> dplyr::group_by(step_number) |>
    dplyr::mutate(y = .smooth_trace_display_only(torque_Nm)) |> dplyr::ungroup()
  p_isov_overlay <- ggplot(isov_raw2, aes(t_rel, y, color = level_lab)) +
    annotate("rect", xmin = 0, xmax = isov_dur, ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.10) +
    stim_vlines(isov_dur) + geom_line(linewidth = 0.9) +
    scale_color_manual(values = LEVEL_COLORS, name = NULL) +
    labs(title = "bass17 isovelocity: raw inertia-corrected torque across velocities",
         subtitle = "No velocity-matched baseline available -- RAW torque shown (dominated by velocity-dependent passive ramp)",
         x = "Time relative to stim onset (s)", y = "Inertia-corrected torque (N*m)") +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
}

p_isov_facets <- build_facets(
  isov_raw, isov_markers, isov_dur,
  "bass17 isovelocity: raw inertia-corrected torque + extraction markers",
  "Free y | orange line = active-window mean (pipeline value) | grey band = [stim_t0, stim_t1+0.5s] averaging window",
  show_passive = FALSE)

ggplot2::ggsave(file.path(OUT_DIR, "diag_within_step_force_isovelocity_overlay.png"),
                p_isov_overlay, width = 9, height = 5.5, dpi = 150)
ggplot2::ggsave(file.path(OUT_DIR, "diag_within_step_force_isovelocity_facets.png"),
                p_isov_facets, width = 8, height = 9, dpi = 150)

cli::cli_alert_success("Saved 4 diagnostic figures to {OUT_DIR}")
