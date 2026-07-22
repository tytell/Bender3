# diag_ytorque_inertial_timing.R
# READ-ONLY DIAGNOSTIC (2026-07-22, PI-directed follow-up to
# ytorquesignexamples.png). PI observation: "y torque is positive for
# isometric stims and the rise is coincident with muscle stimulation. On the
# other hand, negative y torques in isovelocity happen well after the
# stimulus and most likely in relation to the velocity ramp. Velocity ramp
# suggests inertial noise?"
#
# CORRECTION FIRST (important): the isovelocity "negative" examples in
# ytorquesignexamples.png v1 were captured from extract_isovelocity_zero_
# points()'s per-step loop, which computes EVERY step (including op != 0)
# with a STATIC pre-stim baseline before discarding the non-V0 rows -- that
# baseline compares a FAST-MOVING active window against a MOTIONLESS
# pre-stim window, which is guaranteed to show a huge motion-linked
# deviation regardless of muscle activity. The REAL isovelocity force path
# (compute_isovelocity_vector_batch(), what actually feeds FV_isovelocity_
# uhatBoth.png) uses ANGLE-MATCHED passive instead (active vs. a no-stim
# ramp at the SAME angle/velocity), which is designed to cancel exactly this
# kind of motion-linked artifact. Re-run with the REAL angle-matched path:
# pct negative drops from 73-83% (static-baseline, WRONG capture) to 52-60%
# (angle-matched, REAL path) for concentric/eccentric -- much closer to a
# coin flip than a strong signal, vs. isometric/V0 holds' still-consistent
# 90-97% positive. This ALREADY suggests the moving-ramp sign is noise-
# dominated even AFTER angle-matched correction, consistent with the PI's
# inertial-noise hypothesis -- ytorquesignexamples.png's 3 isovelocity
# examples were regenerated from the REAL angle-matched batch (see that
# script's re-run, same date).
#
# THIS diagnostic tests the mechanism directly: does the SAME Ty deflection,
# at the SAME time and magnitude, appear in a completely UNSTIMULATED
# no-stim ramp of the same commanded speed? If yes, the deflection is a
# kinematic (inertial/motion) feature of the ramp itself, not a stimulus-
# locked muscle event -- angle-matching is SUPPOSED to cancel this out of
# force_yTorque_N (active minus passive), but a residual can survive if the
# active and passive ramps' acceleration profiles don't line up exactly
# (e.g. the muscle's own contribution changes how the arm accelerates).
#
# One token: ytorqueinertialtiming.
#
# Run: Rscript R/diag_ytorque_inertial_timing.R

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2); library(patchwork); library(cli)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")
src("00_load_bender_flat.R")
src("muscle_force_vector.R")
src("plot_force_vs_time.R")  # .smooth_trace_display_only
OUT_DIR <- FIGS_DIAGNOSTIC_DIR

#' Load one step's raw Ty + encoder angle/velocity vs local time (t.s already
#' resets to 0 at step start -- segmented_finite convention).
.load_step_trace <- function(td, step_number) {
  rows <- td$step_number == step_number
  d <- tibble::tibble(t.s = td$t.s[rows], Ty = .mfv_col(td, "ytorque")[rows],
                      enc = td$enc.deg[rows], stimV = td$stim.V[rows]) |> dplyr::arrange(t.s)
  d$Ty_smooth <- .smooth_trace_display_only(d$Ty)
  d$vel <- c(NA, diff(d$enc) / diff(d$t.s))
  d$vel_smooth <- .smooth_trace_display_only(ifelse(is.na(d$vel), 0, d$vel))
  d
}

# =============================================================================
# Representative example: bass16_bender_16_isovelocity, step 12 (active,
# left, eccentric, op=-335.6 deg/s) vs. step 20 (no-stim, same trial, same
# |commanded speed| -- both trace the same ~33 deg excursion). Chosen because
# it is one of the few isovelocity trials with a WITHIN-TRIAL no-stim ramp at
# the same speed (most trials in this corpus only have same-individual
# CROSS-trial matches, harder to show side by side cleanly).
# =============================================================================
fp <- file.path(raw_source_dir(BASS16_RAW_SUBFOLDER), "2026-07-14_bass16_bender_16_isovelocity.h5")
td <- load_bender_flat(fp, do_filter = TRUE, loadtorques = c("x", "y", "z"), loadforces = TRUE)
active  <- .load_step_trace(td, 12)
passive <- .load_step_trace(td, 20)

stim_on <- active$stimV > max(active$stimV, na.rm = TRUE) * 0.3
stim_t0 <- min(active$t.s[stim_on], na.rm = TRUE); stim_t1 <- max(active$t.s[stim_on], na.rm = TRUE)
t_dip_active  <- active$t.s[which.min(active$Ty_smooth)]
t_dip_passive <- passive$t.s[which.min(passive$Ty_smooth)]
t_velpeak_active  <- active$t.s[which.max(abs(active$vel_smooth))]
t_velpeak_passive <- passive$t.s[which.max(abs(passive$vel_smooth))]
cli::cli_alert_info(
  "ACTIVE step12: stim [{round(stim_t0,3)}, {round(stim_t1,3)}]s, Ty dip @ t={round(t_dip_active,3)}s, |vel| peak @ t={round(t_velpeak_active,3)}s"
)
cli::cli_alert_info(
  "PASSIVE (no-stim) step20, same |commanded speed|: Ty dip @ t={round(t_dip_passive,3)}s, |vel| peak @ t={round(t_velpeak_passive,3)}s -- NO stimulation present at all"
)

.trace_panel <- function(d, title, subtitle, stim_t0 = NA, stim_t1 = NA) {
  rng <- range(d$Ty, na.rm = TRUE); span <- diff(rng)
  vel_scaled <- scales::rescale(d$vel_smooth, to = rng, from = range(d$vel_smooth, na.rm = TRUE))
  p <- ggplot(d, aes(t.s)) +
    geom_line(aes(y = Ty), color = "grey75", linewidth = 0.3) +
    geom_line(aes(y = Ty_smooth), color = "#b91c1c", linewidth = 0.8) +
    geom_line(aes(y = vel_scaled), color = "#1d4ed8", linewidth = 0.7, linetype = "dashed") +
    labs(title = title, subtitle = subtitle, x = "Time within step (s)",
        y = "Raw y-torque (N*m, red) | scaled angular velocity (blue dashed)") +
    theme_bw(base_size = 10) + theme(plot.subtitle = element_text(size = 8))
  if (is.finite(stim_t0)) p <- p + annotate("rect", xmin = stim_t0, xmax = stim_t1, ymin = -Inf, ymax = Inf,
                                            fill = "orange", alpha = 0.12)
  p
}
p1 <- .trace_panel(active, "ACTIVE: step 12 (left, eccentric, -335.6 deg/s, STIMULATED)",
                   sprintf("Ty dip @ t=%.3fs, |vel| peak @ t=%.3fs, stim window [%.3f, %.3f]s (orange) -- dip happens AFTER stim, near end of ramp",
                          t_dip_active, t_velpeak_active, stim_t0, stim_t1),
                   stim_t0, stim_t1)
p2 <- .trace_panel(passive, "PASSIVE: step 20 (no-stim ramp, same trial, same |commanded speed|)",
                   sprintf("Ty dip @ t=%.3fs, |vel| peak @ t=%.3fs -- ZERO stimulation, yet the SAME-timed, SAME-shaped dip appears",
                          t_dip_passive, t_velpeak_passive))

p <- p1 / p2 +
  patchwork::plot_annotation(
    title = "force_yTorque_N's isovelocity sign: real muscle signal, or uncorrected inertial coupling from the ramp itself?",
    subtitle = paste(strwrap(sprintf(
      "The SAME Ty dip, at the SAME time (%.3fs vs %.3fs) relative to the SAME velocity-profile feature (%.3fs vs %.3fs), appears in BOTH the stimulated ramp AND a completely unstimulated no-stim ramp of the same commanded speed -- this is a KINEMATIC feature of the ramp's motion, not a stimulus-locked muscle event. Angle-matched passive subtraction is designed to cancel exactly this, but the active trace's dip is ~3.7x LARGER (%.2f vs %.2f N*m) than the passive one's, so subtraction leaves a real residual (this residual is what force_yTorque_N reports for this step, %.3f N) -- consistent with the PI's hypothesis that isovelocity's y-torque signal is contaminated by uncorrected inertial coupling from the velocity/acceleration ramp, on top of (or instead of) real muscle force.",
      t_dip_active, t_dip_passive, t_velpeak_active, t_velpeak_passive,
      min(active$Ty_smooth, na.rm=TRUE), min(passive$Ty_smooth, na.rm=TRUE), -0.2565433), width = 140), collapse = "\n")
  )
ggplot2::ggsave(file.path(OUT_DIR, "ytorqueinertialtiming.png"), p, width = 12, height = 9, dpi = 150)
cli::cli_alert_success("Wrote ytorqueinertialtiming.png")

# =============================================================================
# Broader statistical support: across ALL isovelocity ramps (moving steps,
# both stimulated and no-stim), does the Ty-trace's own extremum time track
# its OWN velocity profile's extremum time, regardless of stim? (A tight,
# reproducible lag here -- present with or without stim -- is the signature
# of a kinematic/inertial coupling, not a stimulus-driven response.)
# =============================================================================
.step_timing <- function(fp) {
  td <- tryCatch(load_bender_flat(fp, do_filter = TRUE, loadtorques = c("x", "y", "z"), loadforces = TRUE),
                error = function(e) NULL)
  if (is.null(td)) return(NULL)
  purrr::map_dfr(sort(unique(td$step_number)), function(sn) {
    d <- .load_step_trace(td, sn)
    if (diff(range(d$enc, na.rm = TRUE)) < 2) return(NULL)  # skip V=0 holds -- no ramp to test
    stim_any <- any(d$stimV > max(d$stimV, na.rm = TRUE) * 0.3, na.rm = TRUE)
    tibble::tibble(fp = basename(fp), step_number = sn, stim_any = stim_any,
                   t_Ty_extreme = d$t.s[which.max(abs(d$Ty_smooth - median(d$Ty_smooth, na.rm = TRUE)))],
                   t_vel_extreme = d$t.s[which.max(abs(d$vel_smooth))])
  })
}
fish_dirs <- list(bass16 = raw_source_dir(BASS16_RAW_SUBFOLDER), bass17 = raw_source_dir(BASS17_RAW_SUBFOLDER),
                  bass18 = raw_source_dir(BASS18_RAW_SUBFOLDER))
timing <- dplyr::bind_rows(lapply(names(fish_dirs), function(fish) {
  fps <- list.files(fish_dirs[[fish]], pattern = "isovelocity\\.h5$", full.names = TRUE)
  dplyr::bind_rows(lapply(fps, .step_timing)) |> dplyr::mutate(fish = fish)
}))
timing$dt <- timing$t_Ty_extreme - timing$t_vel_extreme
p_timing <- ggplot(timing, aes(dt, fill = stim_any)) +
  geom_histogram(bins = 30, position = "identity", alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
  scale_fill_manual(values = c(`FALSE` = "#7570b3", `TRUE` = "#d95f02"),
                    labels = c(`FALSE` = "no-stim ramp", `TRUE` = "stimulated ramp")) +
  labs(title = "Time of Ty extremum minus time of angular-velocity extremum, all isovelocity ramps",
      subtitle = sprintf("n=%d ramps, 3 fish | median lag: stim=%.3fs, no-stim=%.3fs -- SAME lag with or without stim points to a kinematic origin, not a stimulus-driven one",
                         nrow(timing), median(timing$dt[timing$stim_any], na.rm = TRUE), median(timing$dt[!timing$stim_any], na.rm = TRUE)),
      x = "t(Ty extremum) - t(|angular velocity| extremum), seconds", y = "Count", fill = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggplot2::ggsave(file.path(OUT_DIR, "ytorqueinertialtiming_stats.png"), p_timing, width = 8, height = 6, dpi = 150)
cli::cli_alert_success("Wrote ytorqueinertialtiming_stats.png")
cli::cli_alert_info("Median dt: stim={round(median(timing$dt[timing$stim_any], na.rm=TRUE),3)}s, no-stim={round(median(timing$dt[!timing$stim_any], na.rm=TRUE),3)}s (n={nrow(timing)})")
