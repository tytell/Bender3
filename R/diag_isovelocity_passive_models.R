# diag_isovelocity_passive_models.R
# READ-ONLY DIAGNOSTIC (2026-07-22, PI-requested; corpus-gap fix + reimplementation
# bug fix added same day). Sources the pipeline and only READS data + WRITES
# PNGs -- modifies no production analysis.
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
# (PI-approved). This script still computes the pointwise result itself, calling
# the SAME production helper (not a reimplementation), so it remains the
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
# CORPUS GAP (found 2026-07-22, PI-directed follow-up to the "flag 1" negative-
# overshoot investigation): the exact-velocity stim-off library ONLY covers a
# SUBSET of commanded velocities per fish (e.g. bass18 has stim-off ramps at
# 107/214/320 deg/s but its OTHER trials move at 142/285/427 deg/s, with NO
# matching stim-off ramp anywhere in the corpus). Production originally fell
# back to a STATIC single-angle baseline for those steps -- meaningless for a
# swept ramp, and the actual source of flag 1's negative overshoot (confirmed:
# 100% of strongly-negative points were `static_baseline_fallback`, ALL
# `angle_matched` points were positive). Fix: production now falls back to the
# NEAREST same-sign-velocity stim-off ramp (still angle-matched, `.mfv_ramp_
# passive_pointwise`'s overlap guard applies) before giving up to the static
# baseline (`passive_source == "angle_matched_nearest_v"`). This script mirrors
# that exact fallback chain via a fish-wide no-stim-ramp library (production
# builds this per-specimen in `compute_isovelocity_vector_batch`; this script
# builds an equivalent one locally so it stays independent of the pipeline).
#
# REIMPLEMENTATION BUG (found + fixed same day): this script's OWN active-vs-
# passive projection (`projg()`) was calling `uhat_geometric(op)` with `op` = the
# step's VELOCITY (deg/s), not a bend angle -- `uhat_geometric()` expects a bend
# ANGLE in degrees, so at high |v| (e.g. 427) it silently computed cos/sin of a
# nonsense "angle" and flipped sign, manufacturing a spurious negative dip in
# Panel 2 that was NOT present in production. `.mfv_finalize_step()` only calls
# `uhat_geometric(s$operating_point)` for category == "isometric" (where
# operating_point genuinely IS the angle); for isovelocity its reference u_hat
# (-> muscle_force_vector_geom_N) is the FIXED longitudinal c(1,0,0). `projg()`
# now uses that same fixed u_hat, matching production exactly.
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
#                    it carries INERTIAL transients at the ramp turnarounds and
#                    differs by sweep direction, a 2nd-order limitation of
#                    angle-alone matching (does not produce sign flips; see
#                    FINDING below for the corrected magnitude).
#
# FINDING (all 3 fish, POST corpus-gap + reimplementation-bug fix): passive
# range across one active window is median 1.8/2.2/2.4 N (bass16/17/18), up to
# ~5.6 N, growing with |v|/sweep; pointwise-minus-production muscle force
# differs by median 0.7/1.2/1.2 N and up to +-5.6 N -- still larger than the
# force itself for the low-force fish. Under production (Method-D active -
# window-MEAN passive) bass16/17 show a CONCAVE-UP FV (force rising at both
# velocity extremes) -- the SAME artifact as the FL concave-up; pointwise
# FLATTENS it to ~0 (appropriate: those fish are at the noise floor). bass18
# (real force) goes from a flat production FV to a genuine, CLEANLY POSITIVE
# curve under pointwise across the FULL velocity range including the
# nearest-velocity-fallback trials (142-427 deg/s) -- NO negative overshoot
# anywhere once (a) the corpus-gap fallback and (b) the uhat_geometric(velocity)
# bug are both fixed. TARGET SHAPE CORRECTION (PI-directed 2026-07-22): FV is
# NOT supposed to be bell-shaped (that's the FL target) -- the physiological
# target is a Hill hyperbola, monotonic-DECREASING with increasing shortening
# velocity (eccentric > isometric > concentric), not a symmetric peak at V=0.
# Checked against that target (SNR-passing right-side points, grouped by
# |velocity|): bass18 shows eccentric > concentric at 127/255 %/s (1.68 vs 1.37
# N; 1.65 vs 1.03 N) -- the correct qualitative ordering -- breaking down only
# at the single highest velocity tested (382 %/s: 1.76 vs 1.83 N, concentric
# slightly exceeds eccentric, worth a closer look since it's the
# nearest-velocity-fallback point). bass16/17 show NO consistent
# eccentric-vs-concentric ordering at any velocity (both near/below the noise
# floor) -- i.e. bass18 pointwise is not just the best-LOOKING curve, it is the
# only one reproducing the correct Hill-type SIGN relationship. The earlier
# "overshoot negative at high |v| -- residual inertial-transient" conclusion
# was an artifact of the two bugs above, not a real limitation of pointwise
# angle-matched passive subtraction. Window-MEAN
# passive still manufactures a velocity-dependent artifact just as the
# isometric static mean did; pointwise (with the corpus-gap fallback) is the
# direct analog fix -- ADOPTED in production 2026-07-22 (.mfv_ramp_passive_
# pointwise + nearest-velocity fallback, PI-approved).
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

  # -- pass 1: load every trial once, build a FISH-WIDE no-stim ramp library
  # (mirrors compute_isovelocity_vector_batch()'s passive_library -- required
  # so a trial with NO no-stim steps of its own, e.g. bass18's 142/285/427
  # series, can still borrow a same-sign-velocity ramp from a DIFFERENT trial
  # of the same fish, exactly as production does).
  trial_cache <- list(); fish_lib <- list()
  for (ti in seq_along(isv)) {
    fp <- isv[ti]
    res <- tryCatch(analyze_isovelocity(load_one(fp), filename = fp), error = function(e) NULL); if (is.null(res)) next
    td6 <- tryCatch(mfv_load_six_axis(fp, res$td), error = function(e) NULL); if (is.null(td6)) next
    geom <- .mfv_read_geom(fp); arms <- resolve_muscle_moment_arms(geom$width_mm, geom$depth_mm_raw, geom$dvert_mm)
    ss <- res$step_summary
    ns <- .mfv_no_stim_steps(td6, unique(td6$step_number))
    trial_cache[[length(trial_cache) + 1L]] <- list(fp = fp, ti = ti, td6 = td6, ss = ss, arms = arms, ns = ns)
    for (sn in ns) {
      op_ns <- ss$operating_point[match(sn, ss$step_number)]
      if (is.finite(op_ns)) {
        fish_lib[[length(fish_lib) + 1L]] <- list(op = op_ns, ramp = .mfv_ramp_from_step(td6, td6$step_number == sn))
      }
    }
  }
  lib_ops <- vapply(fish_lib, function(l) l$op, numeric(1))

  for (tc in trial_cache) {
    fp <- tc$fp; ti <- tc$ti; td6 <- tc$td6; ss <- tc$ss; arms <- tc$arms; ns <- tc$ns
    ns_ops <- ss$operating_point[match(ns, ss$step_number)]
    r_v <- c(0, arms$r_m_m, arms$d_m)
    # u_hat = FIXED longitudinal (1,0,0), NOT uhat_geometric(op): op here is a
    # VELOCITY (deg/s), and uhat_geometric() expects a bend ANGLE (deg) --
    # calling it on a velocity silently computes cos/sin of a nonsense
    # "angle" (e.g. 427 "deg") and flips sign at high |v|, manufacturing a
    # spurious negative dip that ISN'T in production. .mfv_finalize_step()
    # itself only calls uhat_geometric(s$operating_point) for category ==
    # "isometric" (where operating_point genuinely IS the bend angle); for
    # isovelocity/dynamic its `uhat_ref` (-> muscle_force_vector_geom_N) is
    # the fixed c(1,0,0) longitudinal cross-check. Match that here so this
    # diagnostic's own force projection is bit-for-bit consistent with
    # production instead of a divergent reimplementation (found 2026-07-22
    # while reconciling this panel against compute_isovelocity_vector_batch()).
    projg <- function(rows, op = NULL) {
      uh <- c(1.0, 0.0, 0.0); rxu <- .mfv_cross(r_v, uh); eff <- sqrt(sum(rxu^2))
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
      srows <- which(step_rows); t_step <- td6$t.s[srows]; ang_step <- td6$enc.deg[srows]
      act_win_srows <- t_step >= s$stim_t0_s & t_step <= s$stim_t1_s
      ang_act <- ang_step[act_win_srows]
      # Uses the ACTUAL production helper (.mfv_ramp_passive_pointwise), not a
      # hand-rolled reimplementation, so this diagnostic is bit-for-bit
      # consistent with compute_isovelocity_vector_batch() -- including its
      # angle-overlap guard (returns NULL, not an extrapolated result, if the
      # candidate ramp doesn't cover >=50% of this step's active-window
      # angles). Exact-velocity match (within-trial) first; else loop
      # same-sign-velocity candidates from the FISH-WIDE library in order of
      # increasing |Δv| until one clears the overlap guard (mirrors the
      # production fallback chain added 2026-07-22 after this diagnostic's
      # own within-trial-only, single-candidate matching was found to hide
      # the trials whose stim-off ramps didn't cover their moving steps'
      # exact velocities -- see analysis_muscle_force_vector_log.md).
      exact <- ns[abs(ns_ops - op) < MFV_VELOCITY_MATCH_TOL]
      passive_pw <- NULL
      if (length(exact) > 0L) {
        passive_pw <- .mfv_ramp_passive_pointwise(td6, step_rows, ang_act,
                        .mfv_ramp_from_step(td6, td6$step_number == exact[1L]))
      }
      if (is.null(passive_pw)) {
        same_sign <- which(is.finite(lib_ops) & sign(lib_ops) == sign(op))
        if (length(same_sign) > 0L) {
          ord <- same_sign[order(abs(lib_ops[same_sign] - op))]
          for (li in ord) {
            passive_pw <- .mfv_ramp_passive_pointwise(td6, step_rows, ang_act, fish_lib[[li]]$ramp)
            if (!is.null(passive_pw)) break
          }
        }
      }
      if (is.null(passive_pw)) next
      pj <- projg(step_rows, op); g_step <- pj$g; rxu <- pj$rxu; eff <- pj$eff
      pfx <- passive_pw$xtorque[srows]; pfy <- passive_pw$ytorque[srows]; pfz <- passive_pw$ztorque[srows]
      g_pass_full <- (pfx * rxu[1L] + pfy * rxu[2L] + pfz * rxu[3L]) / eff^2
      act_win <- act_win_srows
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
       subtitle = "Window-MEAN passive leaves a velocity-growing residual (red rises with |rate|); pointwise (using every fish's stim-off ramps, nearest-velocity fallback where no exact match exists) removes it cleanly -- no negative overshoot at any |rate|.",
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
