# diag_isometric_strongest_contraction_baseline_visual.R
# PI ask (2026-07-25): "Which of the alternative subtraction methods is most
# likely to resolve true muscle forces? Provide a visual of active traces and
# what the passive calculation might look like -- one 3-panel figure, the
# strongest muscle contraction."
#
# RECOMMENDATION: of the 3 alternatives discussed in chat (pooled/model-based
# passive fit, matched no-stim control step, signature/onset-kinetics
# detection), the MATCHED-CONTROL-STEP method is the gold standard WHERE
# AVAILABLE (already used for isovelocity's real non-stim "bilateral_
# simultaneous" calibration reps) -- but ISOMETRIC trials have NO such rep
# (verified 2026-07-25: every isometric step is stimulated, left_/right_
# unilateral only). That option is unavailable for this protocol as
# currently run. Of the two that remain buildable, this script implements
# the POOLED/MODEL-BASED passive-relaxation fit (alternative #1): instead of
# one step's own noisy 2-window (pre/post) mean, regress passive torque on
# BOTH elapsed session time (captures viscoelastic relaxation across the
# whole trial) AND commanded operating_point (captures the position-
# dependent passive-torque level, which a naive time-only pooled model
# CONFLATES with time -- caught empirically while building this: a first
# time-only-loess draft gave a wildly wrong answer because isometric steps'
# operating_point changes systematically with elapsed time in this ramped
# protocol, and time-only smoothing rides straight through that confound).
#
# "Strongest muscle contraction" = the step with the single largest RAW
# (no-subtraction) |specific tension| across the corpus (isometric_passive_
# baseline_method_comparison.csv, 2026-07-25): bass17, trial
# 2026-07-15_bass17_bender_15_isometric, step 16 (28.4% strain, ~321 kN/m^2
# raw). Same step used in the initial specific-tension investigation's
# /tmp/iso_step16_trace.png.
#
# Run with:  Rscript R/diag_isometric_strongest_contraction_baseline_visual.R
# Output -> 02_processed/figs_diagnostic/isometric_strongest_contraction_
#           baseline_methods_visual.png

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2); library(cli); library(patchwork)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")
src("00_load_bender_flat.R")
src("01_calibrate.R")
src("02_deconvolve.R")
src("muscle_geometry.R")
src("plot_force_vs_time.R")
src("03_analyze.R")
src("parse_trial_filename.R")

OUT_DIR <- FIGS_DIAGNOSTIC_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

TARGET_TRIAL_FILE <- "2026-07-15_bass17_bender_15_isometric.h5"
TARGET_STEP       <- 16L
CLOSEST_WINDOW_S  <- 0.3
COUGHLIN_TENSION_MEAN_KNM2 <- 180
COUGHLIN_TENSION_SD_KNM2   <- 33.6

f <- file.path(raw_source_dir(BASS17_RAW_SUBFOLDER), TARGET_TRIAL_FILE)
td <- load_bender_flat(f, do_filter = TRUE, loadtorques = "x")
tau <- deconvolve_bender(f, hub_path = NULL, verbose = FALSE)
N <- min(nrow(td), length(tau)); td <- td[seq_len(N), , drop = FALSE]
td$torque_inertia_corrected_Nm <- tau[seq_len(N)]
attr(td, "Filename") <- f

# IMPORTANT: t.s is PER-STEP LOCAL time (resets to 0 at every step_number,
# confirmed 2026-07-25: step 1 and step 9 both span t.s in [0, 11.0]) --
# NOT a continuous whole-trial clock. Any cross-step aggregation below joins
# on step_number explicitly and/or uses wall_clock_start for real elapsed
# time; never compares raw t.s values across different steps.
built <- build_segmented_step_summary(td, f)
steps <- built$step_summary
td2   <- built$td
r_m   <- built$r_m

s <- steps[steps$step_number == TARGET_STEP, ]
cli::cli_alert_info("Target step {TARGET_STEP}: side={s$muscle_side}, strain={round(s$shortening_strain_pct,1)}%, force_sign={s$force_sign}")

# ==========================================================================
# Candidate passive-baseline methods for the target step (single-step,
# step_number-scoped windows -- matches diag_isometric_passive_baseline_
# methods.R exactly).
# ==========================================================================
.step_window_mean <- function(td, step_number, t0, t1) {
  rows <- td$step_number == step_number & td$t.s >= t0 & td$t.s <= t1
  if (!any(rows, na.rm = TRUE)) return(NA_real_)
  mean(td$torque_inertia_corrected_Nm[rows], na.rm = TRUE)
}

pre_last_t0   <- s$t_pre_baseline_end_s - CLOSEST_WINDOW_S
post_first_t1 <- s$t_post_baseline_start_s + CLOSEST_WINDOW_S
pre_last      <- .step_window_mean(td2, TARGET_STEP, pre_last_t0, s$t_pre_baseline_end_s)
post_first    <- .step_window_mean(td2, TARGET_STEP, s$t_post_baseline_start_s, post_first_t1)
t_active_mid     <- (s$stim_t0_s + (s$stim_t1_s + 0.5)) / 2
t_pre_mid        <- (s$t_pre_baseline_start_s + s$t_pre_baseline_end_s) / 2
t_post_mid       <- (s$t_post_baseline_start_s + s$t_post_baseline_end_s) / 2
t_pre_last_mid   <- pre_last_t0 + CLOSEST_WINDOW_S / 2
t_post_first_mid <- s$t_post_baseline_start_s + CLOSEST_WINDOW_S / 2

interp_closest <- pre_last + (post_first - pre_last) *
  (t_active_mid - t_pre_last_mid) / (t_post_first_mid - t_pre_last_mid)

# ==========================================================================
# NEW: pooled/model-based passive fit -- regress EVERY step's pre- AND
# post-baseline window MEAN (already computed by build_segmented_step_
# summary() for all 32 steps) on elapsed REAL session time (wall_clock_
# start-based, since t.s resets per step and cannot be compared across
# steps) AND commanded operating_point (needed because this ramped-FL
# protocol's operating_point changes systematically with elapsed time --
# omitting it lets the model conflate position-dependent passive level with
# time-dependent relaxation; see header). Quadratic operating_point term
# since passive torque vs. bend angle is not expected to be linear through
# zero.
# ==========================================================================
wc_posix <- suppressWarnings(as.POSIXct(steps$wall_clock_start, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC"))
t0_session <- min(wc_posix, na.rm = TRUE)
steps$elapsed_s <- as.numeric(difftime(wc_posix, t0_session, units = "secs"))
steps$.t_pre_mid_s  <- (steps$t_pre_baseline_start_s + steps$t_pre_baseline_end_s) / 2
steps$.t_post_mid_s <- (steps$t_post_baseline_start_s + steps$t_post_baseline_end_s) / 2

train <- dplyr::bind_rows(
  tibble::tibble(elapsed_s = steps$elapsed_s + steps$.t_pre_mid_s, operating_point = steps$operating_point,
                torque = steps$passive_force_Nm_static, anchor = "pre"),
  tibble::tibble(elapsed_s = steps$elapsed_s + steps$.t_post_mid_s, operating_point = steps$operating_point,
                torque = steps$post_force_Nm_static, anchor = "post")
) |> dplyr::filter(is.finite(.data$elapsed_s), is.finite(.data$operating_point), is.finite(.data$torque))

pooled_fit <- stats::lm(torque ~ elapsed_s + operating_point + I(operating_point^2), data = train)
target_elapsed_active <- steps$elapsed_s[steps$step_number == TARGET_STEP] + t_active_mid
pooled_pred_at_active <- as.numeric(predict(pooled_fit, newdata = data.frame(
  elapsed_s = target_elapsed_active, operating_point = s$operating_point
)))
cli::cli_alert_info("Pooled regression: R^2={round(summary(pooled_fit)$r.squared,3)}, predicted passive at step {TARGET_STEP} = {round(pooled_pred_at_active,4)} N*m")

active_Nm <- s$active_force_Nm  # legacy peak-window mean, same as production

methods <- tibble::tibble(
  method = factor(c("raw_no_subtraction", "pre_static", "post_static",
                    "interp_linear", "interp_closest", "pooled_regression"),
                  levels = c("raw_no_subtraction", "pre_static", "post_static",
                            "interp_linear", "interp_closest", "pooled_regression")),
  passive_Nm = c(0, s$passive_force_Nm_static, s$post_force_Nm_static,
                s$passive_force_Nm_interp, interp_closest, pooled_pred_at_active)
)
methods$muscle_force_Nm <- s$force_sign * (active_Nm - methods$passive_Nm)
methods$force_N <- methods$muscle_force_Nm / r_m
methods$specific_tension_kNm2 <- (methods$force_N / MEASURED_RED_MUSCLE_CSA_CM2) * 10
print(methods)

# ==========================================================================
# Panel A: WHY the pooled model needs a position term -- all 32 steps' own
# pre/post baseline means vs. commanded operating_point, colored by elapsed
# session time. Target step highlighted.
# ==========================================================================
op_grid <- seq(min(train$operating_point), max(train$operating_point), length.out = 100)
fit_curve <- tibble::tibble(
  operating_point = op_grid,
  torque = as.numeric(predict(pooled_fit, newdata = data.frame(
    elapsed_s = target_elapsed_active, operating_point = op_grid)))
)

pA <- ggplot(train, aes(x = .data$operating_point, y = .data$torque)) +
  geom_point(aes(color = .data$elapsed_s, shape = .data$anchor), size = 2.4, alpha = 0.85) +
  geom_line(data = fit_curve, aes(x = .data$operating_point, y = .data$torque),
            inherit.aes = FALSE, color = "darkorchid", linewidth = 0.8) +
  geom_vline(xintercept = s$operating_point, linetype = "dotted", color = "firebrick") +
  scale_color_viridis_c(name = "Elapsed\nsession\ntime (s)") +
  labs(title = "Why the pooled model needs a POSITION term: passive torque vs. commanded angle",
       subtitle = "All 32 steps' own pre-/post-baseline window means (bass17, this trial). Purple = pooled fit\nat step 16's elapsed time. Red dotted = step 16's own operating_point.",
       x = "Commanded operating_point (deg)", y = "Passive-window torque (N*m)", shape = "Window") +
  theme_bw(base_size = 10)

# ==========================================================================
# Panel B: zoomed to the target step -- raw trace + all candidate passive
# estimates as anchor points, active window shaded.
# ==========================================================================
zoom_lo <- s$t_pre_baseline_start_s - 0.3
zoom_hi <- s$t_post_baseline_end_s + 0.3
td_zoom <- td2[td2$step_number == TARGET_STEP & td2$t.s >= zoom_lo & td2$t.s <= zoom_hi, ]

anchor_pts <- tibble::tibble(
  t = c(t_pre_mid, t_post_mid, t_pre_last_mid, t_post_first_mid, t_active_mid, t_active_mid),
  v = c(s$passive_force_Nm_static, s$post_force_Nm_static, pre_last, post_first,
        pooled_pred_at_active, active_Nm),
  label = c("pre_static", "post_static", "pre_last0.3s", "post_first0.3s", "pooled_regression", "active (raw)")
)

pB <- ggplot() +
  annotate("rect", xmin = s$stim_t0_s, xmax = s$stim_t1_s + 0.5, ymin = -Inf, ymax = Inf,
           fill = "firebrick", alpha = 0.12) +
  geom_line(data = td_zoom, aes(x = .data$t.s, y = .data$torque_inertia_corrected_Nm),
            color = "grey30", linewidth = 0.4) +
  geom_segment(aes(x = t_pre_mid, xend = t_post_mid, y = s$passive_force_Nm_static, yend = s$post_force_Nm_static),
              color = "orange", linewidth = 0.6) +
  geom_hline(yintercept = pooled_pred_at_active, color = "darkorchid", linewidth = 0.6, linetype = "dashed") +
  geom_point(data = anchor_pts, aes(x = .data$t, y = .data$v, color = .data$label), size = 3) +
  scale_color_brewer(palette = "Dark2", name = "Passive estimate / active") +
  coord_cartesian(ylim = c(min(anchor_pts$v) - 0.05, max(anchor_pts$v) + 0.05)) +
  labs(title = sprintf("Zoomed: step %d (strongest contraction, %.1f%% strain)", TARGET_STEP, s$shortening_strain_pct),
       subtitle = "Grey = raw torque (this step only, y-axis clipped to the passive/active comparison range --\nthe 75 Hz stim burst's own artifact spikes go far outside this window, see raw trace note below).\nPurple dashed = pooled position+time regression. Orange = interp_linear's pre-post secant.\nRed band = stim + deactivation tail.",
       x = "Step-local time (s)", y = "Torque (N*m)") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")

# ==========================================================================
# Panel C: resulting specific tension by method (log scale), Coughlin ref.
# ==========================================================================
pC <- ggplot(methods, aes(x = .data$method, y = abs(.data$specific_tension_kNm2))) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = COUGHLIN_TENSION_MEAN_KNM2 - COUGHLIN_TENSION_SD_KNM2,
           ymax = COUGHLIN_TENSION_MEAN_KNM2 + COUGHLIN_TENSION_SD_KNM2,
           fill = "grey70", alpha = 0.3) +
  geom_hline(yintercept = COUGHLIN_TENSION_MEAN_KNM2, linetype = "dashed", color = "grey40") +
  geom_col(aes(fill = .data$method), width = 0.6) +
  geom_text(aes(label = signif(abs(.data$specific_tension_kNm2), 3)), vjust = -0.4, size = 3) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  scale_y_log10() +
  labs(title = "Resulting |specific tension| by method",
       subtitle = "Grey ribbon/dashed = Coughlin (2000), 180 +/- 33.6 kN/m^2",
       x = NULL, y = "|Specific tension| (kN/m^2, log scale)") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 25, hjust = 1))

combined <- pA / pB / pC +
  patchwork::plot_annotation(
    title = "Isometric strongest contraction: raw active trace vs. candidate passive-baseline calculations",
    subtitle = sprintf("bass17, %s, step %d -- largest raw specific tension in the corpus (%.0f kN/m^2 unsubtracted)",
                       TARGET_TRIAL_FILE, TARGET_STEP, abs(methods$specific_tension_kNm2[methods$method == "raw_no_subtraction"]))
  )

fout <- file.path(OUT_DIR, "isometric_strongest_contraction_baseline_methods_visual.png")
ggplot2::ggsave(fout, combined, width = 9, height = 13, dpi = 150)
cli::cli_alert_success("Saved {fout}")
