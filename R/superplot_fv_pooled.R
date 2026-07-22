# superplot_fv_pooled.R
# Pooled Force-Velocity SUPERPLOT across individuals (bass16/17/18) --
# PI-requested 2026-07-22, "similar to FL superplots, raw and normalized."
#
# X-AXIS CONVENTION: shortening_strain_pct, the SAME muscle-centric column
# analyze_isovelocity()/resolve_step_contraction() already computes for
# every step -- for isovelocity this is a STRAIN RATE (%/s: the same
# curvature -> strain geometry as FL's shortening_strain_pct, run on
# operating_point when it's a velocity instead of a position -- see
# 03_analyze.R's own header comment, "FL x-axis / FV x-axis use predicted
# MUSCLE strain / strain-rate"). Crucially this is SIGNED BY CONTRACTION
# MODE (concentric = +, eccentric = -), NOT by raw commanded operating_point
# sign -- two steps at the same |op| on OPPOSITE sides can have opposite raw
# op signs yet the SAME contraction_mode (e.g. left concentric at op=+112
# and right concentric at op=-112, mirrored motor directions for a
# symmetric protocol), so pooling on raw op would mix concentric-here with
# eccentric-there at the same nominal velocity. shortening_strain_pct
# already resolves this (verified empirically 2026-07-22: same |op|, same
# contraction_mode, opposite raw op sign -> IDENTICAL signed
# shortening_strain_pct) -- using it here instead of operating_point is the
# ONLY change from a naive axis-swap of the FL superplot's code.
#
# DESIGN (mirrors superplot_fl_pooled.R's V=0-only rule, axis-swapped): FL
# pooled every V=0 contraction across STRAIN; FV pools every STRAIN=0-ish
# (L0) contraction across STRAIN RATE, PLUS the genuinely-moving isovelocity
# ramps (which sample a RANGE of position/strain during the ramp -- there is
# no single-strain subset of a moving ramp without the sono L0-crossing
# machinery, see NOTE below). Three strain-rate-axis sources are pooled:
#   - ISOMETRIC L0 REP : strain_pct ~= 0 steps only (V=0 anchor, "pure" --
#                        same L0 restriction FL used for isovelocity/dynamic,
#                        just swapped to the other axis here).
#   - ISOVELOCITY V=0   : that trial's own embedded V=0 holds (V=0 anchor,
#                        same points superplot_fl_pooled.R's
#                        extract_isovelocity_zero_points() already computes).
#   - ISOVELOCITY MOVING: the actual concentric/eccentric ramps -- this IS
#                        the velocity-dependent part of the curve. Uses
#                        compute_isovelocity_vector_batch()'s REAL
#                        angle-matched passive (2026-07-22 correction, see
#                        diag_ytorque_sign_examples.R's header) -- angle_
#                        matched(/_cross_trial) sources ONLY; static_
#                        baseline_fallback rows are excluded (same
#                        motion-vs-motionless problem as that correction).
#   - DYNAMIC L0 BOOKEND: V=0 anchor, same bookend points as the FL
#                        superplot's extract_dynamic_l0_points().
#
# NOTE on strain during moving ramps (limitation, not a bug): the per-step
# scalar force plotted here is MEAN-over-the-active-window (the pipeline's
# standard extraction, see bass16_forceextractionmethod.png), sampled over
# whatever RANGE of strain the ramp swept through during that window -- NOT
# force-at-a-fixed-L0 the way a textbook FV curve is defined. A stricter
# alternative exists in the codebase (.mfv_fv_l0_crossing() /
# compute_isovelocity_vector_batch()'s `fv_l0` output -- force sampled AT the
# sono-confirmed L0 crossing) but is right-muscle-only and sono-dependent,
# which would drop most of the data for a first pass. This superplot uses
# the SAME force convention as every other plot in this pipeline instead, so
# it is at least internally consistent; `fv_l0` is a candidate refinement,
# not silently substituted in here.
#
# TRIAL-TO-TRIAL FORCE-SCALE NORMALIZATION: same F/F0 mechanic as the FL
# superplot, axis-swapped: F0 = that SAME trial+side's own mean V=0 force
# (isometric: its own L0 rep; isovelocity: its own V=0 hold; dynamic: its
# own bookend, pre+post averaged) -- physically the SAME F0 quantity FL
# normalizes by (isometric/near-isometric force at the reference state), so
# a V=0 point in the normalized FV plot should land near 1.0 by construction
# (sanity check, not circular: the MOVING points are the real test of
# whether normalization helps).
#
# Sign convention: same as FL (2026-07-22 revision) -- muscle_force_vector_N
# / muscle_force_vector_geom_N reported RAW, no per-step correction. See
# muscle_force_vector.R's module header.
#
# Superplot layers (NO model fit -- connect-the-dots only, per cursorrule),
# identical convention to the FL superplot:
#   - individual points for each trial (per-trial, per-velocity-bin mean) = dots
#   - one THIN line per trial connecting its dots
#   - one line per INDIVIDUAL (fish) = mean of its trials per bin
#   - one THICKER line = GROUP mean across individuals per bin
#
# Run with:  Rscript R/superplot_fv_pooled.R

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(stringr); library(fs); library(ggplot2); library(rhdf5); library(cli)
})

.root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
for (f in c("paths_config.R", "00_load_bender_flat.R", "01_calibrate.R", "02_deconvolve.R",
            "muscle_geometry.R", "fit_fv_fl.R", "03_analyze.R",
            "parse_trial_filename.R", "plot_strain_validation.R",
            "plot_angle_sono_validation.R", "muscle_force_vector.R",
            "plot_force_vs_time.R",
            "extract_dynamic_l0_bookends.R")) {
  source(file.path(.root, f))
}
# .dynamic_l0_trial_rows() (dynamic L0 bookend -> vector-force row, via
# attach_vector_muscle_force()'s td6_override path) is defined in
# superplot_fl_pooled.R -- reused verbatim here rather than re-derived, so
# there is exactly one implementation of "how a dynamic bookend becomes a
# vector-force row" in the codebase. Cut the source at "assemble" so this
# script's OWN assembly/plotting code below (which differs entirely from the
# FL superplot's) is not double-defined by sourcing the whole file.
.fl_pooled_lines <- readLines(file.path(.root, "superplot_fl_pooled.R"))
.cut_at <- grep("^# -+ assemble -+", .fl_pooled_lines)[1]
if (is.na(.cut_at)) cli::cli_abort("superplot_fv_pooled.R: could not find the assemble-section marker in superplot_fl_pooled.R -- has its structure changed?")
.fl_pooled_defs_file <- tempfile(fileext = ".R")
writeLines(.fl_pooled_lines[seq_len(.cut_at - 1L)], .fl_pooled_defs_file)
source(.fl_pooled_defs_file)
unlink(.fl_pooled_defs_file)

SPECIMEN_DIRS <- c(
  bass16 = raw_source_dir(BASS16_RAW_SUBFOLDER),
  bass17 = raw_source_dir(BASS17_RAW_SUBFOLDER),
  bass18 = raw_source_dir(BASS18_RAW_SUBFOLDER)
)
OUTPUT_DIR <- FIGS_DIAGNOSTIC_DIR
fs::dir_create(OUTPUT_DIR, recurse = TRUE)

STRAIN_RATE_BIN_PCT_S <- 40.0  # common strain-rate-bin width for the mean lines / dots
SNR_MIN <- MFV_UHAT_SNR_MIN
L0_STRAIN_EPSILON_PCT <- 0.5   # same tolerance as superplot_fl_pooled.R

.load_trial_td <- function(fullpath) {
  td <- load_bender_flat(fullpath, do_filter = TRUE, loadtorques = "x")
  if (is.null(td) || nrow(td) == 0L) stop("load_bender_flat returned no data")
  tau <- deconvolve_bender(fullpath, hub_path = NULL, verbose = FALSE)
  if (is.null(tau)) stop("deconvolve_bender failed")
  N <- min(nrow(td), length(tau))
  td <- td[seq_len(N), , drop = FALSE]
  td$torque_inertia_corrected_Nm <- tau[seq_len(N)]
  attr(td, "Filename") <- fullpath
  td
}

# =============================================================================
# V=0 anchor sources (strain_rate_pct_s = 0 by construction)
# =============================================================================

#' Isometric trials' own L0 (strain ~= 0) rep(s) -- V=0 anchor, "pure" (single
#' representative length), matching the L0 restriction FL used for its
#' isovelocity/dynamic V=0 points.
extract_isometric_l0_points <- function(fish, manifest) {
  iso <- dplyr::filter(manifest, protocol == "isometric")
  out <- list()
  for (i in seq_len(nrow(iso))) {
    tid <- iso$trial_id[i]; fp <- iso$fullpath[i]
    res <- tryCatch({ td <- .load_trial_td(fp); r <- analyze_isometric(td, filename = fp); r$trial_id <- tid; r },
                    error = function(e) { cli::cli_warn("iso load {tid}: {conditionMessage(e)}"); NULL })
    if (is.null(res)) next
    vec <- tryCatch(attach_vector_muscle_force(res, fp, "isometric"),
                    error = function(e) { cli::cli_warn("iso vec {tid}: {conditionMessage(e)}"); NULL })
    if (is.null(vec)) next
    ss0 <- vec$step_summary |>
      dplyr::filter(.data$muscle_side %in% c("left", "right"),
                    is.finite(.data$shortening_strain_pct),
                    abs(.data$shortening_strain_pct) < L0_STRAIN_EPSILON_PCT)
    if (nrow(ss0) == 0L) next
    out[[tid]] <- ss0 |>
      dplyr::transmute(fish = fish, trial_id = tid, protocol = "isometric", muscle_side = .data$muscle_side,
                       strain_rate_pct_s = 0.0,
                       force_emp_N = .data$muscle_force_vector_N, force_geom_N = .data$muscle_force_vector_geom_N,
                       activation_snr = .data$activation_snr,
                       baseline_force_noise_N = .data$baseline_force_noise_N)
  }
  if (length(out) == 0L) tibble::tibble() else dplyr::bind_rows(out)
}

#' Isovelocity trials' own embedded V=0 holds -- V=0 anchor (same source
#' superplot_fl_pooled.R's extract_isovelocity_zero_points() computes, but
#' here we only need the STATIC (correctly-computed) V=0 rows -- see that
#' script's per-fish batch loop below for why the MOVING rows need the
#' separate angle-matched batch call instead.
extract_isovelocity_v0_points <- function(fish, manifest) {
  isv <- dplyr::filter(manifest, protocol == "isovelocity")
  out <- list()
  for (i in seq_len(nrow(isv))) {
    tid <- isv$trial_id[i]; fp <- isv$fullpath[i]
    res <- tryCatch({ td <- .load_trial_td(fp); r <- analyze_isovelocity(td, filename = fp); r$trial_id <- tid; r },
                    error = function(e) { cli::cli_warn("isv load {tid}: {conditionMessage(e)}"); NULL })
    if (is.null(res)) next
    vec <- tryCatch(attach_vector_muscle_force(res, fp, "isovelocity"),
                    error = function(e) { cli::cli_warn("isv vec {tid}: {conditionMessage(e)}"); NULL })
    if (is.null(vec)) next
    ss0 <- vec$step_summary |>
      dplyr::filter(.data$muscle_side %in% c("left", "right"), .data$contraction_mode == "isometric_zero")
    if (nrow(ss0) == 0L) next
    out[[tid]] <- ss0 |>
      dplyr::transmute(fish = fish, trial_id = tid, protocol = "isovelocity", muscle_side = .data$muscle_side,
                       strain_rate_pct_s = 0.0,
                       force_emp_N = .data$muscle_force_vector_N, force_geom_N = .data$muscle_force_vector_geom_N,
                       activation_snr = .data$activation_snr,
                       baseline_force_noise_N = .data$baseline_force_noise_N)
  }
  if (length(out) == 0L) tibble::tibble() else dplyr::bind_rows(out)
}

#' Dynamic trials' V=0 L0 bookends -- V=0 anchor, identical mechanism to
#' superplot_fl_pooled.R's .dynamic_l0_trial_rows()/extract_dynamic_l0_points()
#' (reused verbatim rather than re-derived -- sourced from that file below).
extract_dynamic_v0_points <- function(fish, manifest) {
  dyn <- dplyr::filter(manifest, protocol == "dynamic")
  out <- list()
  for (i in seq_len(nrow(dyn))) {
    tid <- dyn$trial_id[i]; fp <- dyn$fullpath[i]
    rows <- tryCatch(.dynamic_l0_trial_rows(fp, tid),
                     error = function(e) { cli::cli_warn("dyn load {tid}: {conditionMessage(e)}"); NULL })
    if (is.null(rows) || nrow(rows) == 0L) next
    out[[tid]] <- rows |>
      dplyr::transmute(fish = fish, trial_id = tid, protocol = "dynamic", muscle_side = .data$muscle_side,
                       strain_rate_pct_s = 0.0,
                       force_emp_N = .data$muscle_force_vector_N, force_geom_N = .data$muscle_force_vector_geom_N,
                       activation_snr = .data$activation_snr,
                       baseline_force_noise_N = .data$baseline_force_noise_N)
  }
  if (length(out) == 0L) tibble::tibble() else dplyr::bind_rows(out)
}

# =============================================================================
# Moving (velocity-dependent) source: isovelocity concentric/eccentric ramps,
# via the REAL angle-matched batch (2026-07-22 correction).
# =============================================================================

#' All of one fish's isovelocity trials' MOVING (concentric/eccentric) steps,
#' via compute_isovelocity_vector_batch() (angle-matched passive, matching
#' FV_isovelocity_uhatBoth.png's own path exactly -- NOT the static-baseline
#' shortcut extract_isovelocity_zero_points() uses for its V=0-only return).
#' Batches ALL of this fish's isovelocity trials together (required for
#' cross-trial passive borrowing) -- cannot be called per-trial.
extract_isovelocity_moving_points <- function(fish, manifest) {
  isv <- dplyr::filter(manifest, protocol == "isovelocity")
  if (nrow(isv) == 0L) return(tibble::tibble())
  iso_inputs <- list()
  for (i in seq_len(nrow(isv))) {
    tid <- isv$trial_id[i]; fp <- isv$fullpath[i]
    res <- tryCatch({ td <- .load_trial_td(fp); r <- analyze_isovelocity(td, filename = fp); r$trial_id <- tid; r },
                    error = function(e) { cli::cli_warn("isv load {tid}: {conditionMessage(e)}"); NULL })
    if (!is.null(res)) iso_inputs[[tid]] <- list(trial_id = tid, filename = fp, res = res)
  }
  if (length(iso_inputs) == 0L) return(tibble::tibble())
  batch <- tryCatch(compute_isovelocity_vector_batch(iso_inputs),
                    error = function(e) { cli::cli_warn("{fish} isovelocity batch: {conditionMessage(e)}"); NULL })
  if (is.null(batch)) return(tibble::tibble())
  out <- list()
  for (tid in names(batch$step_summaries)) {
    ss <- batch$step_summaries[[tid]]
    mv <- ss |>
      dplyr::filter(.data$muscle_side %in% c("left", "right"),
                    .data$contraction_mode %in% c("concentric", "eccentric"),
                    # angle-matched(/_cross_trial) only -- static_baseline_fallback has the
                    # same motion-vs-motionless problem as the ytorquesignexamples correction.
                    grepl("^angle_matched", .data$passive_source))
    if (nrow(mv) == 0L) next
    out[[tid]] <- mv |>
      dplyr::transmute(fish = fish, trial_id = tid, protocol = "isovelocity", muscle_side = .data$muscle_side,
                       strain_rate_pct_s = .data$shortening_strain_pct,
                       force_emp_N = .data$muscle_force_vector_N, force_geom_N = .data$muscle_force_vector_geom_N,
                       activation_snr = .data$activation_snr,
                       baseline_force_noise_N = .data$baseline_force_noise_N)
  }
  if (length(out) == 0L) tibble::tibble() else dplyr::bind_rows(out)
}


# ---------------------------------------------------------------- assemble ---

cli::cli_h1("Pooled FV superplot: bass16/17/18 (V=0 anchor + isovelocity moving ramps)")
all_points <- list()
for (fish in names(SPECIMEN_DIRS)) {
  cli::cli_h2(fish)
  manifest <- parse_trial_directory(SPECIMEN_DIRS[[fish]])
  iso_pts <- extract_isometric_l0_points(fish, manifest)
  isv0_pts <- extract_isovelocity_v0_points(fish, manifest)
  isvm_pts <- extract_isovelocity_moving_points(fish, manifest)
  dyn_pts <- extract_dynamic_v0_points(fish, manifest)
  cli::cli_alert_info(
    "{fish}: isometric L0 = {nrow(iso_pts)}; isovelocity V=0 hold = {nrow(isv0_pts)}; isovelocity MOVING = {nrow(isvm_pts)}; dynamic L0 bookend = {nrow(dyn_pts)}")
  all_points[[fish]] <- dplyr::bind_rows(iso_pts, isv0_pts, isvm_pts, dyn_pts)
}
pooled <- dplyr::bind_rows(all_points)
if (nrow(pooled) == 0L) cli::cli_abort("No pooled FV points extracted.")

# F0 (this trial+side's own mean V=0 force, SNR-gated) -- same F/F0 mechanic
# as the FL superplot, computed AFTER pooling so it can draw on ALL of a
# trial's V=0 contributions (its own hold/bookend/L0-rep rows) at once.
#
# TRIAL-LEVEL F0 FALLS BACK TO FISH+SIDE (added 2026-07-22, empirical finding
# during this build): isovelocity's OWN embedded V=0 holds almost never clear
# SNR_MIN (observed range ~0.1-1.9 across all 3 fish, vs. isometric's usual
# several-fold-higher SNR) -- these holds are brief, embedded reps, not
# isometric's dedicated sustained contractions. An isovelocity trial's
# MOVING-ramp rows would therefore get F0 = NA under a strict trial-only
# match, discarding ALL 72 moving points from the normalized plot. Fall back
# to that SAME fish+side's POOLED F0 across every OTHER SNR-passing V=0
# source (isometric L0 reps, dynamic bookends -- i.e. borrowing the same
# physiological reference point across trials of the SAME individual+side,
# not inventing a new one) when the trial's own V=0 doesn't qualify.
# SNR + MAGNITUDE gated (mfv_gate_f0, 2026-07-22): F0 must ALSO exceed its own
# baseline force-noise floor, not just clear SNR (same fix as the FL superplot).
f0_by_trial <- pooled |>
  dplyr::filter(.data$strain_rate_pct_s == 0) |>
  dplyr::group_by(.data$trial_id, .data$muscle_side) |>
  dplyr::summarise(f0_emp  = mfv_gate_f0(.data$force_emp_N,  .data$activation_snr, .data$baseline_force_noise_N),
                   f0_geom = mfv_gate_f0(.data$force_geom_N, .data$activation_snr, .data$baseline_force_noise_N),
                   .groups = "drop")
f0_by_fish <- pooled |>
  dplyr::filter(.data$strain_rate_pct_s == 0) |>
  dplyr::group_by(.data$fish, .data$muscle_side) |>
  dplyr::summarise(f0_emp_fish  = mfv_gate_f0(.data$force_emp_N,  .data$activation_snr, .data$baseline_force_noise_N),
                   f0_geom_fish = mfv_gate_f0(.data$force_geom_N, .data$activation_snr, .data$baseline_force_noise_N),
                   .groups = "drop")
pooled <- pooled |>
  dplyr::left_join(f0_by_trial, by = c("trial_id", "muscle_side")) |>
  dplyr::left_join(f0_by_fish, by = c("fish", "muscle_side")) |>
  dplyr::mutate(f0_emp = dplyr::coalesce(.data$f0_emp, .data$f0_emp_fish),
               f0_geom = dplyr::coalesce(.data$f0_geom, .data$f0_geom_fish),
               force_emp_norm = .data$force_emp_N / .data$f0_emp,
               force_geom_norm = .data$force_geom_N / .data$f0_geom) |>
  dplyr::select(-"f0_emp_fish", -"f0_geom_fish")

method_levels <- c(empirical = "empirical u_hat", geometric = "geometric / longitudinal u_hat")
n_iso <- dplyr::n_distinct(dplyr::filter(pooled, protocol == "isometric")$trial_id)
n_isv <- dplyr::n_distinct(dplyr::filter(pooled, protocol == "isovelocity")$trial_id)
n_dyn <- dplyr::n_distinct(dplyr::filter(pooled, protocol == "dynamic")$trial_id)
n_moving <- sum(pooled$strain_rate_pct_s != 0, na.rm = TRUE)
fish_cols <- c(bass16 = "#1b9e77", bass17 = "#d95f02", bass18 = "#7570b3")

#' Shared FV-superplot builder (mirrors superplot_fl_pooled.R's
#' .build_fl_superplot(), x-axis swapped to velocity).
.build_fv_superplot <- function(pooled, emp_col, geom_col, title, subtitle, y_lab) {
  long <- pooled |>
    tidyr::pivot_longer(c(dplyr::all_of(emp_col), dplyr::all_of(geom_col)),
                        names_to = "method_raw", values_to = "force_N") |>
    dplyr::mutate(method = factor(dplyr::if_else(.data$method_raw == emp_col,
                                                 method_levels[["empirical"]], method_levels[["geometric"]]),
                                  levels = method_levels)) |>
    dplyr::filter(is.finite(.data$force_N), is.finite(.data$strain_rate_pct_s))
  long$velocity_bin <- round(long$strain_rate_pct_s / STRAIN_RATE_BIN_PCT_S) * STRAIN_RATE_BIN_PCT_S

  per_trial <- long |>
    dplyr::group_by(method, fish, trial_id, velocity_bin) |>
    dplyr::summarise(force_N = mean(force_N, na.rm = TRUE), .groups = "drop")
  per_fish <- per_trial |>
    dplyr::group_by(method, fish, velocity_bin) |>
    dplyr::summarise(force_N = mean(force_N, na.rm = TRUE), .groups = "drop")
  per_group <- per_fish |>
    dplyr::group_by(method, velocity_bin) |>
    dplyr::summarise(force_N = mean(force_N, na.rm = TRUE), n_fish = dplyr::n(), .groups = "drop")

  ggplot(mapping = aes(x = velocity_bin, y = force_N)) +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3, linetype = "dashed") +
    geom_line(data = per_trial, aes(group = trial_id, colour = fish), linewidth = 0.3, alpha = 0.35) +
    geom_point(data = per_trial, aes(colour = fish), size = 0.8, alpha = 0.4) +
    geom_line(data = per_fish, aes(group = fish, colour = fish), linewidth = 1.0, alpha = 0.95) +
    geom_point(data = per_fish, aes(colour = fish), size = 2.1, alpha = 0.95) +
    geom_line(data = per_group, aes(group = 1), colour = "black", linewidth = 1.8) +
    geom_point(data = per_group, colour = "black", size = 2.4) +
    facet_wrap(~method, ncol = 2) +
    scale_colour_manual(values = fish_cols, name = "Individual (fish)") +
    labs(title = title, subtitle = subtitle,
        x = "Muscle shortening strain rate (%/s; + = concentric/shortening, - = eccentric/lengthening; 0 = isometric/V=0 anchor)",
        y = y_lab) +
    theme_bw(12) +
    theme(legend.position = "right", plot.subtitle = element_text(size = 8), strip.text = element_text(face = "bold"))
}

p <- .build_fv_superplot(
  pooled, "force_emp_N", "force_geom_N",
  title = "Pooled Force-Velocity superplot (6-axis line-of-action muscle force)",
  subtitle = sprintf(
    "V=0 anchor (isometric L0 reps + isovelocity's own V=0 holds + dynamic L0 bookends) + isovelocity's actual moving ramps (angle-matched passive only, muscle-centric shortening_strain_pct, NOT raw commanded op sign), pooled across %d individuals | %d isometric + %d isovelocity + %d dynamic trial(s), %d moving-ramp points | %g%%/s strain-rate bins | connect-the-dots, NO fit\nthin = per trial, medium = per individual mean, thick black = group mean | RAW absolute force -- see _normalized companion for F/F0 | force is MEAN-over-active-window (see bass16_forceextractionmethod.png), NOT sampled at a fixed L0 crossing -- see this script's header NOTE",
    dplyr::n_distinct(pooled$fish), n_iso, n_isv, n_dyn, n_moving, STRAIN_RATE_BIN_PCT_S),
  y_lab = "Muscle force along u_hat (N, raw sign, no per-step correction)")

n_norm_avail <- dplyr::n_distinct(dplyr::filter(pooled, is.finite(.data$force_emp_norm))$trial_id)
p_norm <- .build_fv_superplot(
  pooled, "force_emp_norm", "force_geom_norm",
  title = "Pooled Force-Velocity superplot, NORMALIZED (F / trial's own V=0 force)",
  subtitle = sprintf(
    "Same pooling as the RAW superplot, each point divided by that TRIAL+SIDE's own mean V=0 force (F/F0) to correct cross-trial force-scale differences before pooling | %d/%d trials contributed a usable F0 | %d individuals | %g%%/s strain-rate bins | connect-the-dots, NO fit\nPROTOTYPE (2026-07-22) -- compare against the RAW file before treating this as canonical",
    n_norm_avail, dplyr::n_distinct(pooled$trial_id), dplyr::n_distinct(pooled$fish), STRAIN_RATE_BIN_PCT_S),
  y_lab = "Muscle force / trial's own V=0 force (F/F0, dimensionless)")

outfile <- file.path(OUTPUT_DIR, "FVsuperplot_isovelocity_pooled.png")
ggplot2::ggsave(outfile, p, width = 13, height = 6.5, dpi = 150)
cli::cli_alert_success("Wrote {outfile}")

outfile_norm <- file.path(OUTPUT_DIR, "FVsuperplot_isovelocity_pooled_normalized.png")
ggplot2::ggsave(outfile_norm, p_norm, width = 13, height = 6.5, dpi = 150)
cli::cli_alert_success("Wrote {outfile_norm}")

cli::cli_alert_info("Pooled points: {nrow(pooled)} rows ({n_moving} moving) | fish: {paste(sort(unique(pooled$fish)), collapse=', ')} | F0 available for {n_norm_avail}/{dplyr::n_distinct(pooled$trial_id)} trials")
