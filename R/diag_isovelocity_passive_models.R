# diag_isovelocity_passive_models.R
# READ-ONLY DIAGNOSTIC (2026-07-22, PI-requested). Sources the pipeline and only
# READS data + WRITES PNGs -- modifies no production analysis.
#
# Question: after fixing the ISOMETRIC passive (relaxation-fit, M2), how does the
# ISOVELOCITY passive baseline compare, and is it right? The PRE-FIX production
# (compute_isovelocity_vector_batch) subtracted an ANGLE+signed-VELOCITY-matched
# no-stim ramp -- the correct raw material (it carries angle-elastic +
# velocity-viscous + inertial together) -- BUT collapsed it to a scalar
# window-MEAN, then subtracted that from the Method-D (peak) active value. This
# diagnostic shows why that mean is wrong and what pointwise (angle-matched,
# sample-by-sample) subtraction changes. RESOLVED: pointwise angle-matched
# subtraction (.mfv_ramp_passive_pointwise) was adopted in production 2026-07-22
# (PI-approved). This script still computes BOTH the mean-collapse and the
# pointwise result itself (independent of the pipeline), so it remains the
# canonical BEFORE/AFTER comparison regardless of the production default.
#
# LOGIC vs ISOMETRIC (the PI's question):
#   ISOMETRIC   -- fixed angle; passive varies only in TIME (viscoelastic
#                  relaxation). 1 d.o.f., bracketed by quiescent pre/post
#                  samples -> fit vs time, subtract pointwise. SOLVED (M2).
#   ISOVELOCITY -- MOVING; passive varies in ANGLE (elastic, LARGE -- up to
#                  ~6 N across one active window here), plus with velocity
#                  (viscous) and direction (hysteresis). The active window
#                  itself sweeps angle, so a window-MEAN passive is a poor
#                  stand-in for the passive at the Method-D peak's OWN angle.
#                  The relaxation-fit-vs-time idea does NOT transfer; the fix is
#                  to subtract the angle-matched ramp POINTWISE (by angle) then
#                  take Method D on the delta -- the direct analog of the
#                  isometric pointwise fix.
#
# Three panels:
#   1 rampshape   -- representative active ramps (bass18, one trial): projected
#                    active g(t) with the POINTWISE angle-matched passive and the
#                    current window-MEAN passive (flat) overlaid; stim shaded.
#                    Shows the passive sweeping across the window vs the flat mean.
#   2 fvpayoff    -- FV: force vs strain rate under PRODUCTION (Method-D active -
#                    MEAN passive) vs POINTWISE, all fish. How the FV curve moves.
#   3 rampstruct  -- no-stim passive g vs ANGLE for +v vs -v ramps (same |v|):
#                    the passive is NOT a clean single-valued elastic curve --
#                    it carries large INERTIAL transients at the ramp turnarounds
#                    and differs by sweep direction, so angle-ALONE matching
#                    leaves a residual even done pointwise (flagged 2nd-order limit).
#
# FINDING (all 3 fish): passive range across one active window is median
# 2.1/3.0/3.1 N (bass16/17/18) and up to ~6 N, growing with |v|/sweep;
# pointwise-minus-production muscle force differs by median 0.2/0.4/0.6 N and up
# to +-4.6 N -- larger than the force itself. Under production (Method-D active -
# window-MEAN passive) the low-force fish bass16/17 show a CONCAVE-UP FV (force
# rising at both velocity extremes) -- the SAME artifact as the FL concave-up;
# pointwise FLATTENS it to ~0 (appropriate: those fish are at the noise floor).
# bass18 (real force) goes from a flat production FV to a plausible bell under
# pointwise, but OVERSHOOTS negative at high |v| -- residual inertial-transient /
# angle-alignment error (Panel 3), a flagged 2nd-order limit. So the window-MEAN
# passive manufactures a velocity-dependent artifact just as the isometric static
# mean did; pointwise is the direct analog fix -- ADOPTED in production 2026-07-22
# (.mfv_ramp_passive_pointwise). This script keeps computing both mean and
# pointwise itself, so it remains the canonical BEFORE/AFTER record.
#
# Run: Rscript R/diag_isovelocity_passive_models.R
# Canon token: isovpassivemodels -- see FIGURES_README.md.

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(stringr)
  library(ggplot2); library(tidyr); library(cli)
})

.root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
for (f in c("paths_config.R","00_load_bender_flat.R","01_calibrate.R","02_deconvolve.R",
            "muscle_geometry.R","fit_fv_fl.R","03_analyze.R","parse_trial_filename.R",
            "plot_force_vs_time.R","muscle_force_vector.R")) source(file.path(.root, f))

OUT_DIR <- FIGS_DIAGNOSTIC_DIR
DIRS <- list(bass16 = raw_source_dir(BASS16_RAW_SUBFOLDER),
             bass17 = raw_source_dir(BASS17_RAW_SUBFOLDER),
             bass18 = raw_source_dir(BASS18_RAW_SUBFOLDER))
DEACT <- 0.5; PW <- 0.15

load_one <- function(f) {
  td <- load_bender_flat(f, do_filter = TRUE, loadtorques = "x")
  tau <- deconvolve_bender(f, hub_path = NULL, verbose = FALSE)
  N <- min(nrow(td), length(tau)); td <- td[seq_len(N), , drop = FALSE]
  td$torque_inertia_corrected_Nm <- tau[seq_len(N)]; attr(td, "Filename") <- f; td
}
method_d <- function(t, g, t0, t1) {
  active <- t >= t0 & t <= (t1 + DEACT); search <- t >= t0 & t <= t1
  fm <- mean(g[active], na.rm = TRUE); ts <- t[search]; vs <- g[search]
  if (sum(is.finite(vs)) < 3 || diff(range(ts, na.rm = TRUE)) < PW) return(fm)
  vsm <- .smooth_trace_display_only(vs); use_max <- abs(max(vsm) - fm) >= abs(min(vsm) - fm)
  tp <- ts[if (use_max) which.max(vsm) else which.min(vsm)]
  nb <- active & t >= (tp - PW / 2) & t <= (tp + PW / 2); mean(g[nb], na.rm = TRUE)
}

trace_rows <- list(); summ <- list(); hyst_rows <- list()
for (fish in names(DIRS)) {
  man <- parse_trial_directory(DIRS[[fish]])
  isv <- man$fullpath[man$protocol == "isovelocity"]
  for (ti in seq_along(isv)) {
    fp <- isv[ti]
    res <- tryCatch(analyze_isovelocity(load_one(fp), filename = fp), error = function(e) NULL); if (is.null(res)) next
    td6 <- tryCatch(mfv_load_six_axis(fp, res$td), error = function(e) NULL); if (is.null(td6)) next
    geom <- .mfv_read_geom(fp); arms <- resolve_muscle_moment_arms(geom$width_mm, geom$depth_mm_raw, geom$dvert_mm)
    ss <- res$step_summary
    ns <- .mfv_no_stim_steps(td6, unique(td6$step_number))
    r_v <- c(0, arms$r_m_m, arms$d_m)
    projg <- function(rows, op) {
      uh <- uhat_geometric(op)$uhat; rxu <- .mfv_cross(r_v, uh); eff <- sqrt(sum(rxu^2))
      Tx <- .mfv_col(td6, "xtorque"); Ty <- .mfv_col(td6, "ytorque"); Tz <- .mfv_col(td6, "ztorque")
      list(g = (Tx[rows] * rxu[1] + Ty[rows] * rxu[2] + Tz[rows] * rxu[3]) / eff^2, rxu = rxu, eff = eff)
    }
    for (j in seq_len(nrow(ss))) {
      s <- ss[j, ]
      if (!s$muscle_side %in% c("left", "right")) next
      if (!s$contraction_mode %in% c("concentric", "eccentric")) next
      if (!is.finite(s$stim_t0_s) || !is.finite(s$stim_t1_s)) next
      step_rows <- td6$step_number == s$step_number; if (!any(step_rows)) next
      op <- s$operating_point
      cand <- ns[abs(ss$operating_point[match(ns, ss$step_number)] - op) < MFV_VELOCITY_MATCH_TOL]
      if (length(cand) == 0L) next
      ramp <- .mfv_ramp_from_step(td6, td6$step_number == cand[1L])
      srows <- which(step_rows); t_step <- td6$t.s[srows]; ang_step <- td6$enc.deg[srows]
      pj <- projg(step_rows, op); g_step <- pj$g; rxu <- pj$rxu; eff <- pj$eff
      interp_full <- function(vp) {
        ok <- is.finite(ramp$ang) & is.finite(vp); if (sum(ok) < 2L) return(rep(NA_real_, length(ang_step)))
        ag <- aggregate(vp[ok], by = list(a = ramp$ang[ok]), FUN = mean)
        approx(ag$a, ag$x, xout = ang_step, rule = 2)$y
      }
      pfx <- interp_full(ramp$T[[1L]]); pfy <- interp_full(ramp$T[[2L]]); pfz <- interp_full(ramp$T[[3L]])
      g_pass_full <- (pfx * rxu[1L] + pfy * rxu[2L] + pfz * rxu[3L]) / eff^2
      act_win <- t_step >= s$stim_t0_s & t_step <= s$stim_t1_s
      mean_pass <- mean(g_pass_full[act_win], na.rm = TRUE)
      md_act <- method_d(t_step, g_step, s$stim_t0_s, s$stim_t1_s)
      prod_muscle <- md_act - mean_pass
      pw_muscle <- method_d(t_step, g_step - g_pass_full, s$stim_t0_s, s$stim_t1_s)
      summ[[length(summ) + 1L]] <- tibble(
        fish = fish, trial = basename(fp), step = s$step_number,
        contraction_mode = s$contraction_mode, muscle_side = s$muscle_side,
        strain_rate = s$shortening_strain_pct, op = op,
        passive_mean = mean_pass, passive_range = diff(range(g_pass_full[act_win], na.rm = TRUE)),
        prod_muscle = prod_muscle, pw_muscle = pw_muscle)
      # panel-1 traces: bass18 first trial, right side, positive velocities
      if (fish == "bass18" && ti == 1L && s$muscle_side == "right" && op > 0) {
        disp <- t_step >= (s$stim_t0_s - 0.15) & t_step <= (s$stim_t1_s + 0.3)
        trace_rows[[length(trace_rows) + 1L]] <- tibble(
          step = s$step_number, vlab = sprintf("%s v=%+.0f deg/s", s$contraction_mode, op),
          trel = t_step[disp] - s$stim_t0_s, g_act = g_step[disp], g_pass = g_pass_full[disp],
          mean_pass = mean_pass, dur = s$stim_t1_s - s$stim_t0_s)
      }
    }
    # panel-3 hysteresis: no-stim passive g vs angle, one +v and one -v pair (bass18 t1)
    if (fish == "bass18" && ti == 1L) {
      for (sn in ns) {
        op <- ss$operating_point[match(sn, ss$step_number)]
        if (!is.finite(op) || abs(abs(op) - 213.6) > 5) next   # a single mid |v|
        rows <- td6$step_number == sn; pj <- projg(rows, op)
        hyst_rows[[length(hyst_rows) + 1L]] <- tibble(
          dir = if (op > 0) "+v (lengthening sweep)" else "-v (shortening sweep)",
          ang = td6$enc.deg[rows], g = pj$g)
      }
    }
  }
}
S <- bind_rows(summ); TR <- bind_rows(trace_rows); HY <- bind_rows(hyst_rows)

# ---- PANEL 1: ramp shape + pointwise vs mean passive ----
p1 <- ggplot(TR, aes(trel)) +
  annotate("rect", xmin = 0, xmax = median(TR$dur), ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.12) +
  geom_line(aes(y = g_act, color = "active g(t)"), linewidth = 0.7) +
  geom_line(aes(y = g_pass, color = "pointwise angle-matched passive"), linewidth = 0.7) +
  geom_line(aes(y = mean_pass, color = "current window-MEAN passive"), linewidth = 0.7, linetype = "dashed") +
  facet_wrap(~vlab, scales = "free_y") +
  scale_color_manual(name = NULL, values = c("active g(t)" = "#111111",
    "pointwise angle-matched passive" = "#059669", "current window-MEAN passive" = "#b91c1c")) +
  labs(title = "Panel 1 -- Isovelocity: active vs POINTWISE angle-matched passive vs current window-MEAN (bass18, one trial)",
       subtitle = "The passive SWEEPS across the active window (muscle is moving through angle); the flat red MEAN is a poor stand-in for the passive at the Method-D peak's own angle.",
       x = "Time relative to stim onset (s)", y = "Projected force g(t) (N)") +
  theme_bw(10) + theme(legend.position = "bottom", plot.subtitle = element_text(size = 8))
ggsave(file.path(OUT_DIR, "isovpassivemodels_1_rampshape.png"), p1, width = 12, height = 7, dpi = 130)

# ---- PANEL 2: FV payoff ----
fv <- S %>% transmute(fish, strain_rate,
  `production (MethodD - MEAN passive)` = prod_muscle,
  `pointwise (angle-matched)` = pw_muscle) %>%
  pivot_longer(-c(fish, strain_rate), names_to = "method", values_to = "F")
p2 <- ggplot(fv, aes(strain_rate, F, color = method)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.8, span = 1.0) +
  geom_point(size = 1.0, alpha = 0.5) +
  facet_wrap(~fish, scales = "free_y") +
  scale_color_manual(values = c("production (MethodD - MEAN passive)" = "#b91c1c", "pointwise (angle-matched)" = "#059669")) +
  labs(title = "Panel 2 -- FV payoff: muscle force vs strain rate under MEAN vs POINTWISE passive",
       subtitle = "Window-MEAN passive leaves a velocity-growing residual (red rises with |rate|); pointwise removes most of it but can overshoot NEGATIVE at high |rate| (residual hysteresis).",
       x = "Strain rate (%/s, signed)", y = "Muscle force (N)") +
  theme_bw(11) + theme(legend.position = "bottom", plot.subtitle = element_text(size = 8))
ggsave(file.path(OUT_DIR, "isovpassivemodels_2_fvpayoff.png"), p2, width = 10, height = 4.5, dpi = 140)

# ---- PANEL 3: hysteresis ----
if (nrow(HY) > 0) {
  p3 <- ggplot(HY, aes(ang, g, color = dir)) +
    geom_path(linewidth = 0.7) +
    scale_color_manual(name = NULL, values = c("+v (lengthening sweep)" = "#2563eb", "-v (shortening sweep)" = "#ea580c")) +
    labs(title = "Panel 3 -- No-stim passive g vs ANGLE, +v vs -v sweep (same |v|, bass18)",
         subtitle = "The passive is NOT a clean single-valued elastic curve: large INERTIAL transients at the ramp turnarounds + direction dependence -> angle-alone matching leaves a residual even done pointwise.",
         x = "Joint angle (deg)", y = "Projected passive g (N)") +
    theme_bw(11) + theme(legend.position = "bottom", plot.subtitle = element_text(size = 8))
  ggsave(file.path(OUT_DIR, "isovpassivemodels_3_hysteresis.png"), p3, width = 8, height = 5, dpi = 140)
}

cli::cli_h2("Isovelocity passive: production (MethodD - MEAN) vs pointwise")
for (fh in unique(S$fish)) {
  d <- S[S$fish == fh, ]
  cli::cli_alert("{fh}: n={nrow(d)}  |passive range across window| median={sprintf('%.2f', median(d$passive_range, na.rm=TRUE))} N  |pointwise-prod| median={sprintf('%.2f', median(abs(d$pw_muscle-d$prod_muscle), na.rm=TRUE))} N  max={sprintf('%.2f', max(abs(d$pw_muscle-d$prod_muscle), na.rm=TRUE))} N")
}
cli::cli_alert_success("Saved isovpassivemodels_1..3 to {OUT_DIR}")
