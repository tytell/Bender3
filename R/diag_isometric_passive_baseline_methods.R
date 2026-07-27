# diag_isometric_passive_baseline_methods.R
# PI question (2026-07-25): "What are some alternatives for calculating the
# passive baseline used in isometric specific-tension?" -- follow-on to the
# specific-tension investigation in analysis_muscle_force_vector_log.md's
# 2026-07-25 "specific tension" addendum, which found that muscle_force_Nm
# (active-minus-passive) is a tiny (~3%) residual on top of a large, still-
# decaying passive viscoelastic torque, so the CHOICE of passive baseline
# has an outsized effect on the resulting specific tension.
#
# This script computes SIX candidate passive-baseline methods per isometric
# step from the SAME raw trace (none of them touch production code -- this
# is a read-only diagnostic; build_segmented_step_summary() (03_analyze.R)
# is called unmodified and only its already-exposed columns/outputs are
# reused or extended):
#
#   pre_static          Production default (passive_force_Nm_static): mean
#                        torque over the FULL pre-stim baseline window
#                        [t_pre_baseline_start_s, t_pre_baseline_end_s].
#   post_static         Mean torque over the FULL post-stim baseline window.
#                        Already computed by build_segmented_step_summary()
#                        but unused by the production static-baseline path --
#                        included here as a comparison anchor (same
#                        construction as pre_static, just sampled later
#                        along the same decay).
#   interp_linear       Production ALTERNATE (passive_force_Nm_interp,
#                        PI-directed 2026-07-16): linear interpolation of
#                        pre_static/post_static, evaluated at the active
#                        window's own mean time. Cancels a LOCALLY LINEAR
#                        relaxation trend instead of assuming the passive
#                        torque is flat across the whole hold.
#   pre_last0.3s        NEW: mean torque over just the LAST 0.3s of the
#                        pre-stim window (closest in time to stim onset,
#                        minimizing both (a) contamination from ramp-
#                        settling ringing early in the window and (b) the
#                        time-gap the static/interp methods must
#                        extrapolate/interpolate across).
#   interp_closest      NEW: same linear-interpolation CONSTRUCTION as
#                        interp_linear, but using narrow closest-in-time
#                        sub-windows (pre_last0.3s and a symmetric
#                        post_first0.3s: mean over the FIRST 0.3s of the
#                        post-stim window) as the two anchors instead of the
#                        full pre/post window means. Shortens the
#                        interpolation span, which should reduce residual
#                        error from the true (nonlinear/exponential) decay
#                        being approximated as locally linear.
#   raw_no_subtraction  NOT a passive-baseline method (zero baseline) --
#                        included ONLY as a diagnostic upper-bound anchor,
#                        per the prior investigation's finding that raw
#                        active-window torque alone falls in Coughlin's
#                        range. Do not treat as a real candidate fix.
#
# Restricted to muscle_side == "right" steps (the sono-instrumented muscle),
# matching the production isometric-tension script's side convention (see
# diag_precondition_tension_vs_offset_isometric.R header) -- all 5 isometric
# trials across bass16/17/18, ALL preconditioning phases (not just "later
# (stable)") to maximize n for this method-comparison (precondition is kept
# as a plotting facet/color, not a filter).
#
# Run with:  Rscript R/diag_isometric_passive_baseline_methods.R
# Outputs -> figures: 02_processed/figs_diagnostic/ (FIGS_DIAGNOSTIC_DIR)
#            data:    02_processed/data_processed/

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(ggplot2); library(cli); library(readr)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")
src("00_load_bender_flat.R")
src("01_calibrate.R")
src("02_deconvolve.R")
src("muscle_geometry.R")
src("plot_force_vs_time.R")  # .smooth_trace_display_only() -- required by build_segmented_step_summary()
src("03_analyze.R")
src("parse_trial_filename.R")
src("dynamic_trial_precondition.R")  # extract_bender_trial_num(), classify_session_precondition()

OUT_DIR      <- FIGS_DIAGNOSTIC_DIR
DATA_OUT_DIR <- file.path(.crittergripper_root(), "02_processed", "data_processed")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SPECIMEN_SUBFOLDERS <- c(bass16 = BASS16_RAW_SUBFOLDER, bass17 = BASS17_RAW_SUBFOLDER, bass18 = BASS18_RAW_SUBFOLDER)
SPECIMEN_COLORS     <- c(bass16 = "#1b9e77", bass17 = "#d95f02", bass18 = "#7570b3")
CLOSEST_WINDOW_S    <- 0.3  # width of the "closest-in-time" sub-windows (pre_last / post_first)

# Coughlin (2000) reference band -- same values as summary_coughlin2000_bass_
# comparison.R::COUGHLIN2000$tension_kNm2 (mean=180, sd=33.6 kN/m^2).
COUGHLIN_TENSION_MEAN_KNM2 <- 180
COUGHLIN_TENSION_SD_KNM2   <- 33.6

# =============================================================================
# Part A: load every isometric trial, compute all 6 baseline methods per step
# =============================================================================
cli::cli_h1("Part A: isometric passive-baseline method comparison, bass16/17/18")

.load_one <- function(f) {
  td <- load_bender_flat(f, do_filter = TRUE, loadtorques = "x")
  tau <- deconvolve_bender(f, hub_path = NULL, verbose = FALSE)
  N <- min(nrow(td), length(tau)); td <- td[seq_len(N), , drop = FALSE]
  td$torque_inertia_corrected_Nm <- tau[seq_len(N)]
  attr(td, "Filename") <- f
  td
}

#' Mean torque over [t0, t1] for one step's samples in td.
.window_mean <- function(td, step_number, t0, t1) {
  if (!is.finite(t0) || !is.finite(t1)) return(NA_real_)
  rows <- td$step_number == step_number & td$t.s >= t0 & td$t.s <= t1
  if (!any(rows, na.rm = TRUE)) return(NA_real_)
  mean(td$torque_inertia_corrected_Nm[rows], na.rm = TRUE)
}

.collect_isometric_steps <- function(specimen, source_dir) {
  manifest <- parse_trial_directory(source_dir)
  files <- manifest$fullpath[manifest$protocol == "isometric"]
  purrr::map_dfr(files, function(f) {
    trial_id <- tools::file_path_sans_ext(basename(f))
    td0 <- tryCatch(.load_one(f), error = function(e) NULL)
    if (is.null(td0)) return(tibble())
    built <- tryCatch(build_segmented_step_summary(td0, f), error = function(e) {
      cli::cli_alert_danger("{trial_id}: build_segmented_step_summary failed: {conditionMessage(e)}"); NULL
    })
    if (is.null(built) || nrow(built$step_summary) == 0L) return(tibble())
    if (!is.finite(built$r_m) || built$r_m <= 0) {
      cli::cli_alert_warning("{trial_id}: missing/invalid r_m -- skipped")
      return(tibble())
    }
    steps <- built$step_summary
    td    <- built$td

    # --- NEW candidate baseline #1: pre_last0.3s (closest-to-onset anchor) ---
    steps$passive_force_Nm_pre_last <- vapply(seq_len(nrow(steps)), function(i) {
      s <- steps[i, ]
      .window_mean(td, s$step_number, s$t_pre_baseline_end_s - CLOSEST_WINDOW_S, s$t_pre_baseline_end_s)
    }, numeric(1L))

    # --- NEW candidate baseline #2 input: post_first0.3s (symmetric anchor) ---
    post_first <- vapply(seq_len(nrow(steps)), function(i) {
      s <- steps[i, ]
      .window_mean(td, s$step_number, s$t_post_baseline_start_s, s$t_post_baseline_start_s + CLOSEST_WINDOW_S)
    }, numeric(1L))

    # --- NEW candidate baseline #2: interp_closest -- same linear-interp
    # construction as production's passive_force_Nm_interp, but anchored on
    # the narrow closest-in-time sub-windows above instead of the full pre/
    # post window means (shorter interpolation span).
    t_pre_last_mid   <- steps$t_pre_baseline_end_s - CLOSEST_WINDOW_S / 2
    t_post_first_mid <- steps$t_post_baseline_start_s + CLOSEST_WINDOW_S / 2
    t_active_mid     <- (steps$stim_t0_s + (steps$stim_t1_s + 0.5)) / 2
    steps$passive_force_Nm_interp_closest <- dplyr::if_else(
      is.finite(post_first) & is.finite(steps$passive_force_Nm_pre_last) & (t_post_first_mid != t_pre_last_mid),
      steps$passive_force_Nm_pre_last +
        (post_first - steps$passive_force_Nm_pre_last) *
          (t_active_mid - t_pre_last_mid) / (t_post_first_mid - t_pre_last_mid),
      steps$passive_force_Nm_pre_last
    )

    steps$specimen  <- specimen
    steps$trial_id  <- trial_id
    steps$trial_num <- extract_bender_trial_num(trial_id)
    steps$r_m       <- built$r_m
    cli::cli_alert_success("{trial_id}: {sum(steps$muscle_side == 'right', na.rm = TRUE)} right-side steps")
    steps
  })
}

iso_steps <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  .collect_isometric_steps(specimen, raw_source_dir(subfolder))
})
iso_steps$precondition <- classify_session_precondition(iso_steps$specimen, iso_steps$trial_num)

right_steps <- iso_steps |> dplyr::filter(.data$muscle_side == "right")
cli::cli_alert_info("{nrow(right_steps)} right-side isometric steps across {dplyr::n_distinct(right_steps$trial_id)} trials (all preconditioning phases)")

# =============================================================================
# Part B: fold into long format, one row per (step x method)
# =============================================================================
METHOD_LEVELS <- c("raw_no_subtraction", "pre_static", "post_static",
                    "pre_last0.3s", "interp_linear", "interp_closest")
METHOD_LABELS <- c(
  raw_no_subtraction = "raw (no subtraction)",
  pre_static          = "pre-static (production)",
  post_static         = "post-static",
  "pre_last0.3s"      = "pre, last 0.3s",
  interp_linear       = "interp (full pre/post) [prod. alt.]",
  interp_closest      = "interp (closest 0.3s anchors)"
)

.tension_from_passive <- function(steps, passive_Nm) {
  muscle_force_Nm <- steps$force_sign * (steps$active_force_Nm - passive_Nm)
  force_N <- muscle_force_Nm / steps$r_m
  force_N / MEASURED_RED_MUSCLE_CSA_CM2
}

long <- purrr::map_dfr(METHOD_LEVELS, function(m) {
  passive_Nm <- switch(m,
    raw_no_subtraction = 0,
    pre_static          = right_steps$passive_force_Nm_static,
    post_static         = right_steps$post_force_Nm_static,
    "pre_last0.3s"      = right_steps$passive_force_Nm_pre_last,
    interp_linear       = right_steps$passive_force_Nm_interp,
    interp_closest      = right_steps$passive_force_Nm_interp_closest
  )
  tibble::tibble(
    specimen = right_steps$specimen, trial_id = right_steps$trial_id,
    step_number = right_steps$step_number, precondition = right_steps$precondition,
    shortening_strain_pct = right_steps$shortening_strain_pct,
    method = factor(m, levels = METHOD_LEVELS),
    specific_tension_kNm2 = .tension_from_passive(right_steps, passive_Nm) * 10  # N/cm^2 -> kN/m^2
  )
})

write.csv(long, file.path(DATA_OUT_DIR, "isometric_passive_baseline_method_comparison.csv"), row.names = FALSE)
cli::cli_alert_success("Saved {nrow(long)} (step x method) rows -> data_processed/isometric_passive_baseline_method_comparison.csv")

summary_tbl <- long |>
  dplyr::group_by(.data$method) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_abs_kNm2   = mean(abs(.data$specific_tension_kNm2), na.rm = TRUE),
    median_abs_kNm2 = median(abs(.data$specific_tension_kNm2), na.rm = TRUE),
    max_abs_kNm2    = max(abs(.data$specific_tension_kNm2), na.rm = TRUE),
    .groups = "drop"
  )
cli::cli_h2("Per-method summary (|specific tension|, kN/m^2), all right-side steps pooled")
print(summary_tbl, n = 20)
cli::cli_alert_info("Coughlin (2000) reference: {COUGHLIN_TENSION_MEAN_KNM2} +/- {COUGHLIN_TENSION_SD_KNM2} kN/m^2")

# per-trial max (same aggregation as the production Coughlin panel B), by method
trial_max <- long |>
  dplyr::group_by(.data$specimen, .data$trial_id, .data$precondition, .data$method) |>
  dplyr::summarise(max_abs_kNm2 = max(abs(.data$specific_tension_kNm2), na.rm = TRUE), .groups = "drop")
write.csv(trial_max, file.path(DATA_OUT_DIR, "isometric_passive_baseline_method_trialmax.csv"), row.names = FALSE)
cli::cli_h2("Per-trial max |specific tension| (kN/m^2), by method")
print(trial_max |> tidyr::pivot_wider(names_from = "method", values_from = "max_abs_kNm2"), n = 20)

# =============================================================================
# Part C: plots -> figs_diagnostic/
# =============================================================================
long$method <- factor(long$method, levels = METHOD_LEVELS, labels = METHOD_LABELS[METHOD_LEVELS])

# Plot 1: per-step tension vs. strain, one connected "ladder" per step across
# methods, faceted by specimen -- shows how much a SINGLE step's tension
# moves just from switching baseline method.
p1 <- ggplot(long, aes(x = abs(.data$shortening_strain_pct), y = .data$specific_tension_kNm2,
                        group = interaction(.data$trial_id, .data$step_number))) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = COUGHLIN_TENSION_MEAN_KNM2 - COUGHLIN_TENSION_SD_KNM2,
           ymax = COUGHLIN_TENSION_MEAN_KNM2 + COUGHLIN_TENSION_SD_KNM2,
           fill = "grey70", alpha = 0.3) +
  geom_hline(yintercept = COUGHLIN_TENSION_MEAN_KNM2, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
  geom_line(color = "grey70", alpha = 0.5, linewidth = 0.3) +
  geom_point(aes(color = .data$method), size = 1.8, alpha = 0.85) +
  facet_wrap(~specimen, scales = "free_x") +
  labs(title = "Isometric specific tension: passive-baseline method comparison, per step",
       subtitle = paste0("Right-muscle steps, all 5 isometric trials (all preconditioning phases), all 6 candidate baseline methods.\n",
                          "Grey ribbon/dashed line = Coughlin (2000) tetanic reference, ",
                          COUGHLIN_TENSION_MEAN_KNM2, " +/- ", COUGHLIN_TENSION_SD_KNM2, " kN/m^2.\n",
                          "Thin grey lines connect the SAME step across methods -- shows within-step sensitivity to baseline choice."),
       x = "|Predicted muscle shortening strain| (%L0)", y = "Specific tension (kN/m^2)", color = "Baseline method") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
fout1 <- file.path(OUT_DIR, "isometric_passive_baseline_methods_perstep.png")
ggplot2::ggsave(fout1, p1, width = 11, height = 6.5, dpi = 150)
cli::cli_alert_success("Saved {fout1}")

# Plot 2: distribution of |specific tension| by method, pooled across all
# right-side steps/trials -- shows the overall shift in scale each method
# produces relative to Coughlin's band.
p2 <- ggplot(long, aes(x = .data$method, y = abs(.data$specific_tension_kNm2))) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = COUGHLIN_TENSION_MEAN_KNM2 - COUGHLIN_TENSION_SD_KNM2,
           ymax = COUGHLIN_TENSION_MEAN_KNM2 + COUGHLIN_TENSION_SD_KNM2,
           fill = "grey70", alpha = 0.3) +
  geom_hline(yintercept = COUGHLIN_TENSION_MEAN_KNM2, linetype = "dashed", color = "grey40") +
  geom_boxplot(aes(fill = .data$method), outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_jitter(width = 0.12, size = 1, alpha = 0.4) +
  scale_y_log10() +
  labs(title = "Isometric |specific tension| distribution by passive-baseline method (log scale)",
       subtitle = paste0("All right-muscle steps pooled, bass16/17/18, all 5 isometric trials. Grey ribbon/dashed line = Coughlin (2000), ",
                          COUGHLIN_TENSION_MEAN_KNM2, " +/- ", COUGHLIN_TENSION_SD_KNM2, " kN/m^2."),
       x = NULL, y = "|Specific tension| (kN/m^2, log scale)") +
  theme_bw(base_size = 11) + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
fout2 <- file.path(OUT_DIR, "isometric_passive_baseline_methods_distribution.png")
ggplot2::ggsave(fout2, p2, width = 9, height = 6, dpi = 150)
cli::cli_alert_success("Saved {fout2}")

# Plot 3: per-trial max tension by method (same aggregation as the
# production Coughlin panel B) -- which method's TRIAL-MAX lands closest to
# Coughlin's band at the aggregate level actually used downstream.
trial_max$method <- factor(trial_max$method, levels = METHOD_LEVELS, labels = METHOD_LABELS[METHOD_LEVELS])
p3 <- ggplot(trial_max, aes(x = .data$method, y = .data$max_abs_kNm2, group = .data$trial_id)) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = COUGHLIN_TENSION_MEAN_KNM2 - COUGHLIN_TENSION_SD_KNM2,
           ymax = COUGHLIN_TENSION_MEAN_KNM2 + COUGHLIN_TENSION_SD_KNM2,
           fill = "grey70", alpha = 0.3) +
  geom_hline(yintercept = COUGHLIN_TENSION_MEAN_KNM2, linetype = "dashed", color = "grey40") +
  geom_line(aes(color = .data$specimen), alpha = 0.6, linewidth = 0.6) +
  geom_point(aes(color = .data$specimen, shape = .data$precondition), size = 2.8) +
  scale_y_log10() +
  scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen") +
  labs(title = "Isometric trial-max |specific tension| by passive-baseline method (log scale)",
       subtitle = paste0("One line per trial (n=5 trials total, all preconditioning phases) tracking how the SAME trial's headline\n",
                          "tension number (max over right-side steps -- the production Coughlin-panel aggregation) moves across methods.\n",
                          "Grey ribbon/dashed line = Coughlin (2000), ", COUGHLIN_TENSION_MEAN_KNM2, " +/- ", COUGHLIN_TENSION_SD_KNM2, " kN/m^2."),
       x = NULL, y = "Trial-max |specific tension| (kN/m^2, log scale)", shape = "Preconditioning") +
  theme_bw(base_size = 11) + theme(legend.position = "right", axis.text.x = element_text(angle = 30, hjust = 1))
fout3 <- file.path(OUT_DIR, "isometric_passive_baseline_methods_trialmax.png")
ggplot2::ggsave(fout3, p3, width = 10, height = 6, dpi = 150)
cli::cli_alert_success("Saved {fout3}")

cli::cli_alert_success("diag_isometric_passive_baseline_methods.R complete -- outputs in {OUT_DIR}/ and {DATA_OUT_DIR}/")
