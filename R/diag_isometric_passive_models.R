# diag_isometric_passive_models.R
# READ-ONLY DIAGNOSTIC / PROTOTYPE (2026-07-22, PI-requested). Sources the
# pipeline and only READS data + WRITES PNGs -- modifies no production analysis.
#
# Question: the isometric FL concave-up survives drift correction, fatigue, and
# stim (all previously ruled out). It is a PASSIVE-SUBTRACTION residual. This
# prototypes improved isometric passive subtraction and VISUALIZES the reasoning
# in 4 plots, comparing three passive models (projected onto the geometric u_hat,
# so the quantity IS the production muscle_force_vector_geom_N):
#   M0 -- static pre-stim window mean (current production).
#   M1 -- pre->post linear interpolation (existing passive_force_Nm_interp idea),
#         subtracted POINTWISE (at each active sample's own time), then Method D.
#   M2 -- viscoelastic relaxation fit (loess over the quiescent pre+post samples)
#         subtracted POINTWISE, then Method D. This is the pointwise passive.
#
# Four plots (each settles one decision):
#   1 relaxshape        -- do the quiescent samples show a relaxation curve M2
#                          fits and M0/M1 miss? (bass17, one trial, per step)
#   2 leverage          -- passive(model) - passive(M0) at the active time vs
#                          |bend|: can a better baseline touch the arms at all?
#   3 zeromusclecontrol -- a NO-activation quiescent window minus M0 (should be
#                          ~0); if it scales with |bend| like the real force, the
#                          concave-up is subtraction residual (model-free proof).
#   4 flshape           -- payoff: muscle force vs signed strain under M0/M1/M2,
#                          annotated with cor(F,|strain|).
#
# FINDINGS (bass16/17/18): M0 pre-stim baseline is STALE (~0.5-1 s behind the
# active window on a relaxing passive). The zero-muscle control reproduces
# 121-169% of the bass16/17 |bend| slope with NO muscle -> the concave-up is
# largely subtraction artifact for the low-force fish. Pointwise M2 removes it
# without a spurious bell: cor(F,|strain|) M0->M2 bass16 +0.19->+0.03, bass17
# +0.93->+0.25, bass18 +0.57->+0.42 (bass18's monotonic rise is genuine, not the
# artifact). CAVEAT: M2 interpolates the passive across the ~0.3 s stim gap
# (unobservable), assuming the contraction does not discontinuously perturb the
# passive -- true FL cannot be resolved below the passive-drift floor for the
# low-force fish, which must be magnitude/SNR-gated, not trusted as a shape.
#
# Run: Rscript R/diag_isometric_passive_models.R
# Canon token: isopassivemodels -- see FIGURES_README.md.

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
RELAX <- 1.0; DEACT <- 0.5; PW <- 0.15  # post-stim relax margin, deactivation tail, Method-D narrow window (s)

load_one <- function(f) {
  td <- load_bender_flat(f, do_filter = TRUE, loadtorques = "x")
  tau <- deconvolve_bender(f, hub_path = NULL, verbose = FALSE)
  N <- min(nrow(td), length(tau)); td <- td[seq_len(N), , drop = FALSE]
  td$torque_inertia_corrected_Nm <- tau[seq_len(N)]; attr(td, "Filename") <- f; td
}

# projected torque-force g(t) onto the geometric u_hat (basis of muscle_force_vector_geom_N)
proj_g <- function(td6, rows, s, arms) {
  uh <- uhat_geometric(s$operating_point)$uhat; r <- c(0, arms$r_m_m, arms$d_m)
  rxu <- .mfv_cross(r, uh); eff <- sqrt(sum(rxu^2))
  Tx <- .mfv_col(td6, "xtorque"); Ty <- .mfv_col(td6, "ytorque"); Tz <- .mfv_col(td6, "ztorque")
  g <- (Tx * rxu[1] + Ty * rxu[2] + Tz * rxu[3]) / eff^2
  tibble(t = td6$t.s[rows], g = g[rows])
}
# Method D on a (possibly baseline-subtracted) trace
method_d_val <- function(t, g, t0, t1) {
  active <- t >= t0 & t <= (t1 + DEACT); search <- t >= t0 & t <= t1
  fm <- mean(g[active], na.rm = TRUE); ts <- t[search]; vs <- g[search]
  if (sum(is.finite(vs)) < 3 || diff(range(ts, na.rm = TRUE)) < PW) return(fm)
  vsm <- .smooth_trace_display_only(vs); use_max <- abs(max(vsm) - fm) >= abs(min(vsm) - fm)
  tp <- ts[if (use_max) which.max(vsm) else which.min(vsm)]
  nb <- active & t >= (tp - PW / 2) & t <= (tp + PW / 2); mean(g[nb], na.rm = TRUE)
}
# loess viscoelastic-relaxation fit over quiescent samples, predicted at tnew
m2_predict <- function(tq, gq, tnew) {
  ok <- is.finite(tq) & is.finite(gq); tq <- tq[ok]; gq <- gq[ok]
  if (length(tq) < 8) return(rep(NA_real_, length(tnew)))
  idx <- seq(1, length(tq), by = max(1, floor(length(tq) / 400)))
  df <- data.frame(t = tq[idx], g = gq[idx])
  fit <- tryCatch(loess(g ~ t, data = df, span = 0.6, degree = 2,
                        control = loess.control(surface = "direct")), error = function(e) NULL)
  if (is.null(fit)) return(rep(NA_real_, length(tnew)))
  as.numeric(tryCatch(predict(fit, newdata = data.frame(t = tnew)), error = function(e) rep(NA_real_, length(tnew))))
}

trace_rows <- list(); summ <- list()
for (fish in names(DIRS)) {
  man <- parse_trial_directory(DIRS[[fish]])
  iso <- man$fullpath[man$protocol == "isometric"]
  for (ti in seq_along(iso)) {
    fp <- iso[ti]
    res <- tryCatch(analyze_isometric(load_one(fp), filename = fp), error = function(e) NULL); if (is.null(res)) next
    td6 <- tryCatch(mfv_load_six_axis(fp, res$td), error = function(e) NULL); if (is.null(td6)) next
    geom <- .mfv_read_geom(fp); arms <- resolve_muscle_moment_arms(geom$width_mm, geom$depth_mm_raw, geom$dvert_mm)
    ss <- res$step_summary
    for (j in seq_len(nrow(ss))) {
      s <- ss[j, ]
      if (!s$muscle_side %in% c("left", "right") || !is.finite(s$stim_t0_s)) next
      if (!is.finite(s$t_post_baseline_start_s)) next
      rows <- td6$step_number == s$step_number; if (!any(rows)) next
      G <- proj_g(td6, rows, s, arms); t <- G$t; g <- G$g; t0 <- s$stim_t0_s; t1 <- s$stim_t1_s
      pre  <- t >= s$t_pre_baseline_start_s & t <= s$t_pre_baseline_end_s
      post <- t >= s$t_post_baseline_start_s & t <= s$t_post_baseline_end_s
      quies <- (t >= s$t_pre_baseline_start_s & t <= (t0 - 0.02)) | (t >= (t1 + RELAX) & t <= s$t_post_baseline_end_s)
      M0 <- mean(g[pre], na.rm = TRUE); postv <- mean(g[post], na.rm = TRUE)
      pre_mid  <- (s$t_pre_baseline_start_s + s$t_pre_baseline_end_s) / 2
      post_mid <- (s$t_post_baseline_start_s + s$t_post_baseline_end_s) / 2
      act_mid  <- (t0 + t1 + DEACT) / 2
      M1   <- M0 + (postv - M0) * (act_mid - pre_mid) / (post_mid - pre_mid)   # scalar @ act_mid (leverage)
      M1_t <- M0 + (postv - M0) * (t - pre_mid) / (post_mid - pre_mid)         # pointwise line
      gL_t <- m2_predict(t[quies], g[quies], t)                               # pointwise relaxation fit
      M2   <- m2_predict(t[quies], g[quies], act_mid); if (!is.finite(M2)) M2 <- M1
      av   <- method_d_val(t, g, t0, t1)
      mM0 <- method_d_val(t, g - M0,   t0, t1)
      mM1 <- method_d_val(t, g - M1_t, t0, t1)
      mM2 <- if (all(is.finite(gL_t))) method_d_val(t, g - gL_t, t0, t1) else mM1
      pw0 <- t >= (t1 + RELAX) & t <= (t1 + RELAX + PW); pseudo <- mean(g[pw0], na.rm = TRUE)
      summ[[length(summ) + 1L]] <- tibble(
        fish = fish, trial = basename(fp), step = s$step_number,
        strain = s$shortening_strain_pct, abs_strain = abs(s$shortening_strain_pct),
        M0 = M0, M1 = M1, M2 = M2, active = av, pseudo = pseudo,
        muscle_M0 = mM0, muscle_M1 = mM1, muscle_M2 = mM2,
        pseudo_force = pseudo - M0, artifact_M2 = M2 - M0)
      if (fish == "bass17" && ti == 1) {
        keep <- t >= (s$t_pre_baseline_start_s - 0.3) & t <= s$t_post_baseline_end_s
        M1L <- M0 + (postv - M0) * (t - pre_mid) / (post_mid - pre_mid)
        trace_rows[[length(trace_rows) + 1L]] <- tibble(
          step = s$step_number, strain = round(s$shortening_strain_pct),
          t = t[keep], trel = t[keep] - t0, g = g[keep], M0 = M0, M1 = M1L[keep],
          quies = quies[keep], t0 = t0, t1 = t1,
          slab = sprintf("%+d%% (step %d)", round(s$shortening_strain_pct), s$step_number))
      }
    }
  }
}
S <- bind_rows(summ); TR <- bind_rows(trace_rows)

# ---- PLOT 1: relaxation shape + model fit (bass17 one trial, per step) ----
TR2 <- TR %>% group_by(step) %>% group_modify(function(d, k) {
  q <- d[d$quies, ]
  if (nrow(q) >= 8) {
    fit <- tryCatch(loess(g ~ t, data = q, span = 0.6, degree = 2, control = loess.control(surface = "direct")), error = function(e) NULL)
    d$M2 <- if (!is.null(fit)) as.numeric(predict(fit, newdata = d)) else d$M1
  } else d$M2 <- d$M1
  d
}) %>% ungroup()
dur <- median(TR$t1 - TR$t0)
p1 <- ggplot(TR2, aes(trel, g)) +
  annotate("rect", xmin = 0, xmax = dur, ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.12) +
  geom_point(data = ~subset(.x, quies), color = "grey55", size = 0.25, alpha = 0.5) +
  geom_line(aes(y = M0, color = "M0 static"), linewidth = 0.6) +
  geom_line(aes(y = M1, color = "M1 linear interp"), linewidth = 0.6) +
  geom_line(aes(y = M2, color = "M2 relaxation fit"), linewidth = 0.6) +
  scale_color_manual(name = "passive model", values = c("M0 static" = "#1d4ed8", "M1 linear interp" = "#059669", "M2 relaxation fit" = "#b91c1c")) +
  facet_wrap(~slab, scales = "free_y", ncol = 4) +
  labs(title = "Plot 1 -- Isometric passive relaxation shape & model fit (bass17, one trial)",
       subtitle = "grey = quiescent (no-activation) projected torque-force; orange = stim gap (passive unobserved). Does M2 track the relaxation across the gap where M0/M1 miss it?",
       x = "Time relative to stim onset (s)", y = "Projected force g(t) (N)") +
  theme_bw(10) + theme(legend.position = "bottom", plot.subtitle = element_text(size = 8))
ggsave(file.path(OUT_DIR, "isopassivemodels_1_relaxshape.png"), p1, width = 13, height = 7, dpi = 130)

# ---- PLOT 2: baseline leverage vs |bend| ----
lev <- S %>% transmute(fish, abs_strain, `M1-M0` = M1 - M0, `M2-M0` = M2 - M0) %>%
  pivot_longer(c(`M1-M0`, `M2-M0`), names_to = "model", values_to = "delta")
p2 <- ggplot(lev, aes(abs_strain, delta, color = model)) +
  geom_hline(yintercept = 0, color = "grey70") + geom_point(size = 1.1, alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.7) +
  facet_wrap(~fish, scales = "free_y") +
  scale_color_manual(values = c("M1-M0" = "#059669", "M2-M0" = "#b91c1c")) +
  labs(title = "Plot 2 -- Baseline leverage: how much a better model changes subtracted passive vs |bend|",
       subtitle = "passive estimate at the active-window time, model minus static M0. Growing with |bend| = it can remove that much of the arms.",
       x = "|muscle strain| (%)", y = "passive(model) - passive(M0)  (N)") +
  theme_bw(11) + theme(plot.subtitle = element_text(size = 8))
ggsave(file.path(OUT_DIR, "isopassivemodels_2_leverage.png"), p2, width = 10, height = 4.5, dpi = 140)

# ---- PLOT 3: zero-muscle control ----
z <- S %>% transmute(fish, abs_strain, `real muscle force (active-M0)` = abs(muscle_M0),
  `zero-muscle control (quiescent-M0)` = abs(pseudo_force)) %>%
  pivot_longer(-c(fish, abs_strain), names_to = "kind", values_to = "val")
p3 <- ggplot(z, aes(abs_strain, val, color = kind)) +
  geom_point(size = 1.1, alpha = 0.7) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.7) +
  facet_wrap(~fish, scales = "free_y") +
  scale_color_manual(values = c("real muscle force (active-M0)" = "#111111", "zero-muscle control (quiescent-M0)" = "#e11d48"), name = NULL) +
  labs(title = "Plot 3 -- Zero-muscle control: does subtracting one quiescent window from the baseline leak with |bend|?",
       subtitle = "red = a no-activation window minus M0 (should be ~0 if subtraction were perfect). red ~ black at large |bend| => the concave-up is subtraction residual.",
       x = "|muscle strain| (%)", y = "|force| (N)") +
  theme_bw(11) + theme(legend.position = "bottom", plot.subtitle = element_text(size = 8))
ggsave(file.path(OUT_DIR, "isopassivemodels_3_zeromusclecontrol.png"), p3, width = 10, height = 5, dpi = 140)

# ---- PLOT 4: payoff, FL shape under each method ----
fl <- S %>% transmute(fish, strain, M0 = muscle_M0, M1 = muscle_M1, M2 = muscle_M2) %>%
  pivot_longer(c(M0, M1, M2), names_to = "method", values_to = "F")
p4 <- ggplot(fl, aes(strain, F, color = method)) +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.8, span = 1.0) +
  geom_point(size = 0.9, alpha = 0.5) +
  facet_wrap(~fish, scales = "free_y") +
  scale_color_manual(values = c("M0" = "#1d4ed8", "M1" = "#059669", "M2" = "#b91c1c")) +
  labs(title = "Plot 4 -- Payoff: FL shape (muscle force vs signed strain) under each passive model",
       subtitle = "M0 static | M1 pointwise linear interp | M2 pointwise relaxation fit. cor(F,|strain|): ~0 = flat; + = concave-up; - = bell.",
       x = "Muscle strain (%)", y = "Muscle force (N)") +
  theme_bw(11) + theme(plot.subtitle = element_text(size = 8))
ggsave(file.path(OUT_DIR, "isopassivemodels_4_flshape.png"), p4, width = 10, height = 4.5, dpi = 140)

ann <- fl %>% group_by(fish, method) %>% summarise(r = cor(F, abs(strain), use = "complete.obs"), .groups = "drop")
cli::cli_h2("Payoff: cor(F,|strain|) by passive model")
for (fh in unique(ann$fish)) {
  d <- ann[ann$fish == fh, ]
  cli::cli_alert("{fh}: M0={sprintf('%+.2f', d$r[d$method=='M0'])}  M1={sprintf('%+.2f', d$r[d$method=='M1'])}  M2={sprintf('%+.2f', d$r[d$method=='M2'])}")
}
cli::cli_h2("Zero-muscle control: |force|~|strain| slope (real vs no-muscle)")
for (fh in unique(S$fish)) {
  d <- S[S$fish == fh, ]
  rs <- coef(lm(abs(muscle_M0) ~ abs_strain, d))[2]; cs <- coef(lm(abs(pseudo_force) ~ abs_strain, d))[2]
  cli::cli_alert("{fh}: real slope={sprintf('%.4f', rs)}  control slope={sprintf('%.4f', cs)}  (control = {sprintf('%.0f%%', 100*cs/rs)} of real)")
}
cli::cli_alert_success("Saved isopassivemodels_1..4 to {OUT_DIR}")
