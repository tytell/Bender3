# diag_isometric_baseline_drift.R
# READ-ONLY DIAGNOSTIC (2026-07-22, PI-requested). Sources the pipeline and only
# READS data + WRITES a PNG -- modifies no analysis.
#
# Purpose (isometric baseline-DRIFT test): the isometric hold is ~6 s long, and
# the raw inertia-corrected torque CREEPS monotonically across the whole hold
# (viscoelastic stress relaxation of the bent specimen, amplitude scaling with
# |bend|). The production baseline is the pre-stim window mean (STATIC), sampled
# early in that creep, so active - static drifts upward for seconds after
# stim-off -- not an activation twitch. This figure overlays, per isometric step
# (all fish, faceted), the muscle force under the CURRENT static baseline (left)
# vs a pre->post INTERPOLATED baseline (right; the same linear pre/post scheme as
# 03_analyze.R's passive_force_Nm_interp). If the drift is a linear relaxation
# trend, the interpolated column returns to ~0 after the stim transient.
#
# NOTE: worked in the legacy single-axis inertia-corrected torque domain (which
# already carries pre/post baseline windows + the interp reference), NOT the
# 6-axis vector projection the FL superplot samples. The drift MECHANISM and its
# correction are torque-level; porting interp into the vector path is a separate
# (production) change, deliberately not done here.
#
# Run: Rscript R/diag_isometric_baseline_drift.R
# Canon token: isometricbaselinedrift -- see FIGURES_README.md.

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(stringr)
  library(ggplot2); library(tidyr); library(cli)
})

.root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.root, f))
for (f in c("paths_config.R","00_load_bender_flat.R","01_calibrate.R","02_deconvolve.R",
            "muscle_geometry.R","fit_fv_fl.R","03_analyze.R","parse_trial_filename.R",
            "plot_force_vs_time.R")) src(f)

OUT_DIR <- FIGS_DIAGNOSTIC_DIR
DIRS <- list(bass16 = raw_source_dir(BASS16_RAW_SUBFOLDER),
             bass17 = raw_source_dir(BASS17_RAW_SUBFOLDER),
             bass18 = raw_source_dir(BASS18_RAW_SUBFOLDER))
PAD <- 0.2

load_one <- function(f) {
  td <- load_bender_flat(f, do_filter = TRUE, loadtorques = "x")
  tau <- deconvolve_bender(f, hub_path = NULL, verbose = FALSE)
  N <- min(nrow(td), length(tau)); td <- td[seq_len(N), , drop = FALSE]
  td$torque_inertia_corrected_Nm <- tau[seq_len(N)]
  attr(td, "Filename") <- f
  td
}

rows <- list(); pk <- list()
for (fish in names(DIRS)) {
  man <- parse_trial_directory(DIRS[[fish]])
  for (fp in man$fullpath[man$protocol == "isometric"]) {
    r <- tryCatch(analyze_isometric(load_one(fp), filename = fp), error = function(e) NULL)
    if (is.null(r)) next
    td <- r$td; ss <- r$step_summary
    for (j in seq_len(nrow(ss))) {
      s <- ss[j, ]
      if (!s$muscle_side %in% c("left", "right") || !is.finite(s$stim_t0_s)) next
      if (!is.finite(s$t_pre_baseline_start_s) || !is.finite(s$t_post_baseline_start_s)) next
      sr <- which(td$step_number == s$step_number); if (!length(sr)) next
      sub <- td[sr, ]
      keep <- sub$t.s >= (s$t_pre_baseline_start_s - PAD) & sub$t.s <= (s$t_post_baseline_end_s + PAD)
      t_rel <- sub$t.s[keep] - s$stim_t0_s
      raw   <- sub$torque_inertia_corrected_Nm[keep]
      pre_mean <- s$passive_force_Nm_static; post_mean <- s$post_force_Nm_static
      pre_mid  <- (s$t_pre_baseline_start_s  + s$t_pre_baseline_end_s)  / 2 - s$stim_t0_s
      post_mid <- (s$t_post_baseline_start_s + s$t_post_baseline_end_s) / 2 - s$stim_t0_s
      interp_line <- pre_mean + (post_mean - pre_mean) * (t_rel - pre_mid) / (post_mid - pre_mid)
      f_static <- raw - pre_mean; f_interp <- raw - interp_line
      dur <- s$stim_t1_s - s$stim_t0_s
      instim <- t_rel >= 0 & t_rel <= dur
      post_lo <- s$t_post_baseline_start_s - s$stim_t0_s; post_hi <- s$t_post_baseline_end_s - s$stim_t0_s
      orient <- sign(mean(f_static[instim], na.rm = TRUE)); if (!is.finite(orient) || orient == 0) orient <- 1
      grp <- paste(basename(fp), s$step_number)
      rows[[length(rows) + 1L]] <- tibble(
        fish = fish, grp = grp, strain = s$shortening_strain_pct, t_rel = t_rel,
        static = orient * .smooth_trace_display_only(f_static),
        interp = orient * .smooth_trace_display_only(f_interp), dur = dur)
      pk[[length(pk) + 1L]] <- tibble(
        fish = fish, strain = s$shortening_strain_pct,
        peak_static = max(orient * f_static[instim], na.rm = TRUE),
        peak_interp = max(orient * f_interp[instim], na.rm = TRUE),
        resid_post_static = orient * mean(f_static[t_rel >= post_lo & t_rel <= post_hi], na.rm = TRUE),
        resid_post_interp = orient * mean(f_interp[t_rel >= post_lo & t_rel <= post_hi], na.rm = TRUE))
    }
  }
}

D <- bind_rows(rows); PK <- bind_rows(pk)
DL <- D %>% pivot_longer(c(static, interp), names_to = "baseline", values_to = "y")
DL$baseline <- factor(DL$baseline, levels = c("static", "interp"),
  labels = c("STATIC pre-stim baseline (current)", "pre->post INTERPOLATED baseline"))
dur <- median(D$dur)

p <- ggplot(DL, aes(t_rel, y, group = grp, color = strain)) +
  annotate("rect", xmin = 0, xmax = dur, ymin = -Inf, ymax = Inf, fill = "orange", alpha = 0.12) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_line(linewidth = 0.45, alpha = 0.85) +
  scale_color_gradient2(low = "#2166ac", mid = "grey80", high = "#b2182b", midpoint = 0, name = "strain (%)") +
  facet_grid(fish ~ baseline, scales = "free_y") +
  labs(title = "Isometric baseline-DRIFT test: does pre->post interpolation remove the post-stim creep?",
       subtitle = "Oriented so in-stim deflection is +; orange = stim | LEFT = current static baseline, RIGHT = interpolated | hold is ~6 s; if creep is linear, RIGHT returns to ~0",
       x = "Time relative to stim onset (s)", y = "Muscle force (N*m, oriented +)") +
  theme_bw(11) + theme(legend.position = "right", plot.subtitle = element_text(size = 8))
ggsave(file.path(OUT_DIR, "isometricbaselinedrift.png"), p, width = 11, height = 9, dpi = 140)

shp <- function(v, s) {
  ok <- is.finite(v) & is.finite(s); v <- v[ok]; s <- s[ok]
  ra <- cor(v, abs(s)); rs <- cor(v, s)
  verdict <- if (is.finite(ra) && ra < -0.2) "CONCAVE-DOWN (bell, max at L0)"
             else if (is.finite(ra) && ra > abs(rs) + 0.1) "concave-up (U)"
             else if (is.finite(rs) && abs(rs) > ra + 0.1) "monotonic (signed)"
             else "flat/mixed"
  sprintf("cor(F,|strain|)=%+.2f cor(F,signed)=%+.2f -> %s", ra, rs, verdict)
}
cli::cli_h2("Per-fish shape: peak in-stim force vs strain, STATIC vs INTERP")
for (fh in unique(PK$fish)) {
  d <- PK[PK$fish == fh, ]
  cli::cli_alert("{fh} (n={nrow(d)})")
  cli::cli_alert("  STATIC: {shp(d$peak_static, d$strain)}")
  cli::cli_alert("  INTERP: {shp(d$peak_interp, d$strain)}")
  cli::cli_alert("  mean |resid at post-baseline|: static={sprintf('%.3f', mean(abs(d$resid_post_static), na.rm=TRUE))} interp={sprintf('%.3f', mean(abs(d$resid_post_interp), na.rm=TRUE))}")
}
cli::cli_alert_success("Saved isometricbaselinedrift.png to {OUT_DIR}")
