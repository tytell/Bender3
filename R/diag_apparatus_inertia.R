# diag_apparatus_inertia.R
# READ-ONLY DIAGNOSTIC (2026-07-22, PI-directed). Pre-wiring evidence for
# generalizing the Bender inertial correction from uniaxial (zTorque only) to
# multi-axial, and for fixing Gap 1 (apparatus MOI is a silent no-op in every
# real trial: calibration_inertia_apparatus_moi_gram_millimeter_squared = NaN,
# so 02_deconvolve.R drops the apparatus term entirely).
#
# This script CHANGES NOTHING in the production path. It reads the raw archive
# and the committed fit artifact, computes evidence, and writes two figures +
# a printed numeric summary to figs_diagnostic/. Per the repo rule "reproduce/
# diagnose before changing", it must run and be reviewed BEFORE 02_deconvolve.R
# is touched.
#
# Two tokens:
#   apparatusinertiafit       -- empty-apparatus corpus: per-channel (all 6
#                                F/T axes) empirical inertia I vs angular
#                                acceleration alpha, how well geometry (aor,
#                                width) explains each channel's I (F4 vs F5),
#                                and validation of the stored zTorque F4 fit.
#                                Answers PI decisions 1 (form) and 3 (per-
#                                channel vs coupled: this single-bending-axis
#                                corpus can only identify a per-channel I-vs-
#                                alpha VECTOR, not a full cross-axis tensor --
#                                the script demonstrates that limitation).
#   multiaxialinertiacompare  -- real specimen trials: per-axis raw torque vs
#                                (raw - sign*I*alpha), with the sign auto-chosen
#                                on the most inertia-dominated window (PI
#                                decision 2), variance-reduction quantified per
#                                axis, and the yTorque residual re-examined
#                                against ytorqueinertialtiming (PI decision 4).
#
# Run: Rscript R/diag_apparatus_inertia.R

suppressPackageStartupMessages({
  library(rhdf5); library(dplyr); library(tibble); library(purrr)
  library(ggplot2); library(patchwork); library(cli)
})

.root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.root, f))
src("paths_config.R")
src("apparatus_inertia_fit.R")
src("inertia_diag_common.R")   # CHANNELS, unit consts, CACHE_DIR, .read_trial, .alpha_rad_s2, .per_channel_inertia

OUT_DIR <- FIGS_DIAGNOSTIC_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

APPARATUS_SUBFOLDER <- Sys.getenv("BENDER3_APPARATUS_CALIB_SUBFOLDER",
                                  "2026-07-06_apparatus_inertial_calibration")
FIT_JSON <- Sys.getenv("BENDER3_APPARATUS_FIT_JSON",
                       "2026-07-06_apparatus_inertia_calibration.json")


# =============================================================================
# PART A: empty-apparatus corpus -- per-channel inertia + geometry model
# =============================================================================
cli::cli_h1("PART A: empty-apparatus corpus (per-channel inertia vs geometry)")

art <- read_apparatus_inertia_fit(resolve_apparatus_fit_path(FIT_JSON))
json_trials <- art$trials  # list of per-trial dicts (file, aor, width, I zTorque, r2)
corpus_dir <- raw_source_dir(APPARATUS_SUBFOLDER)
excluded <- unlist(art$excluded_files)

corpus <- purrr::map_dfr(json_trials, function(tr) {
  fp <- file.path(corpus_dir, tr$file)
  if (!file.exists(fp)) { cli::cli_alert_warning("missing corpus file: {tr$file}"); return(NULL) }
  cli::cli_alert_info("reading {tr$file}")
  d <- .read_trial(fp)
  if (is.null(d)) return(NULL)
  alpha <- .alpha_rad_s2(d$t_s, d$angle_deg)
  pc <- .per_channel_inertia(d$ft, alpha)
  pc$file <- tr$file
  pc$aor_mm <- as.numeric(tr$aor_millimeter)
  pc$width_mm <- as.numeric(tr$width_millimeter)
  pc$excluded <- tr$file %in% excluded
  pc
})

fit_corpus <- corpus |> dplyr::filter(!excluded, is.finite(I_gmm2))

# Per-channel geometry model: how well does (aor, width) explain each channel's I?
# F4 = a + b*aor^2 + c*width^2 ; F5 = a + b*(aor^2 + (width/2)^2). Report R^2.
geom_fit <- purrr::map_dfr(CHANNELS, function(ch) {
  d <- fit_corpus |> dplyr::filter(channel == ch)
  if (nrow(d) < 4L) return(tibble(channel = ch, f4_r2 = NA_real_, f5_r2 = NA_real_, n = nrow(d)))
  f4 <- stats::lm(I_gmm2 ~ I(aor_mm^2) + I(width_mm^2), data = d)
  f5 <- stats::lm(I_gmm2 ~ I(aor_mm^2 + (width_mm / 2)^2), data = d)
  tibble(channel = ch, f4_r2 = summary(f4)$r.squared, f5_r2 = summary(f5)$r.squared, n = nrow(d),
         mean_r2_vs_alpha = mean(d$r2, na.rm = TRUE),
         mean_I_gmm2 = mean(d$I_gmm2, na.rm = TRUE))
})

cli::cli_h2("Per-channel inertia-vs-alpha and geometry-model quality (corpus)")
print(as.data.frame(geom_fit), row.names = FALSE, digits = 4)

# Validate the stored zTorque F4 fit against this script's own extraction.
zt <- fit_corpus |> dplyr::filter(channel == "zTorque")
zt$I_f4_pred <- vapply(seq_len(nrow(zt)), function(i)
  apparatus_inertia_from_fit(art, zt$aor_mm[i], zt$width_mm[i])$i_gmm2, numeric(1L))
zt$I_f5_pred <- vapply(seq_len(nrow(zt)), function(i)
  apparatus_inertia_from_fit(art, zt$aor_mm[i], zt$width_mm[i], form_id = "F5")$i_gmm2, numeric(1L))
zt$I_json <- vapply(zt$file, function(f) {
  tr <- json_trials[[which(vapply(json_trials, function(x) x$file, "") == f)]]
  as.numeric(tr$i_gram_millimeter_squared)
}, numeric(1L))

cli::cli_h2("zTorque: this script's I vs JSON-stored I vs F4/F5 predictions (g*mm^2)")
print(as.data.frame(zt[, c("file", "aor_mm", "width_mm", "I_gmm2", "I_json", "I_f4_pred", "I_f5_pred", "r2")]),
      row.names = FALSE, digits = 5)

# ---- Part A figure ---------------------------------------------------------
pA1 <- ggplot(fit_corpus, aes(channel, r2)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.15, height = 0, alpha = 0.6, color = "#1d4ed8") +
  labs(title = "(a) Per-channel inertia signal: R^2 of channel vs angular acceleration",
       subtitle = "High R^2 = that axis carries a real I*alpha term (needs correcting); low = negligible",
       x = "F/T channel", y = "R^2 (channel ~ alpha), per corpus trial") +
  theme_bw(base_size = 10)

pA2 <- ggplot(dplyr::filter(geom_fit, is.finite(f4_r2)),
              aes(x = channel)) +
  geom_col(aes(y = f4_r2), fill = "#b91c1c", width = 0.4,
           position = position_nudge(x = -0.2)) +
  geom_col(aes(y = f5_r2), fill = "#0369a1", width = 0.4,
           position = position_nudge(x =  0.2)) +
  labs(title = "(b) How well geometry (aor, width) explains each channel's I -- F4 (red) vs F5 (blue)",
       subtitle = "F4 = a + b*aor^2 + c*width^2 ; F5 = a + b*(aor^2 + (width/2)^2)",
       x = "F/T channel", y = "geometry-model R^2") +
  theme_bw(base_size = 10)

pA3 <- ggplot(zt, aes(I_json, I_f4_pred)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = "F4"), size = 2) +
  geom_point(aes(y = I_f5_pred, color = "F5"), size = 2) +
  scale_color_manual(values = c(F4 = "#b91c1c", F5 = "#0369a1"), name = "form") +
  labs(title = "(c) zTorque: predicted vs stored empirical I (g*mm^2)",
       subtitle = "Points near the dashed y=x line predict well at real trial geometries",
       x = "empirical I (JSON per-trial)", y = "predicted I") +
  theme_bw(base_size = 10)

pA <- (pA1 / pA2 / pA3) +
  patchwork::plot_annotation(
    title = "apparatusinertiafit: empty-apparatus corpus, per-channel inertia and geometry model",
    subtitle = "PI decisions 1 (F4 vs F5) and 3 (which channels need correcting; single bending axis => per-channel vector, not a full tensor)")
ggplot2::ggsave(file.path(OUT_DIR, "apparatusinertiafit.png"), pA, width = 10, height = 12, dpi = 150)
cli::cli_alert_success("Wrote apparatusinertiafit.png")

# Per-channel F4 geometry coefficients, reused by Part B to predict apparatus I
# at a specimen's geometry for EACH axis (mini Phase-1 fit, for the compare only).
channel_f4_models <- setNames(lapply(CHANNELS, function(ch) {
  d <- fit_corpus |> dplyr::filter(channel == ch)
  if (nrow(d) < 4L) return(NULL)
  stats::lm(I_gmm2 ~ I(aor_mm^2) + I(width_mm^2), data = d)
}), CHANNELS)

predict_channel_I_gmm2 <- function(ch, aor_mm, width_mm) {
  m <- channel_f4_models[[ch]]
  if (is.null(m)) return(NA_real_)
  as.numeric(stats::predict(m, newdata = data.frame(aor_mm = aor_mm, width_mm = width_mm)))
}


# =============================================================================
# PART B: real specimen trials -- per-axis corrected vs uncorrected + sign
# =============================================================================
cli::cli_h1("PART B: real specimen trials (multi-axial corrected vs uncorrected)")

# bass16: one isovelocity (moving ramps -> inertia-active) + one isometric
# (motionless -> alpha ~ 0 control). Extendable via env.
specimen_files <- local({
  d16 <- raw_source_dir(BASS16_RAW_SUBFOLDER)
  iso  <- list.files(d16, pattern = "isovelocity\\.h5$", full.names = TRUE)
  isom <- list.files(d16, pattern = "isometric\\.h5$", full.names = TRUE)
  c(if (length(iso)  > 0) iso[[1]],
    if (length(isom) > 0) isom[[1]])
})

.window_stats <- function(raw, alpha, I_kg_m2) {
  # variance of raw vs both correction signs, on this window.
  cand_minus <- raw - I_kg_m2 * alpha   # sign = +1 : M - I*alpha
  cand_plus  <- raw + I_kg_m2 * alpha   # sign = -1 : M + I*alpha
  v0 <- stats::var(raw, na.rm = TRUE)
  vm <- stats::var(cand_minus, na.rm = TRUE)
  vp <- stats::var(cand_plus, na.rm = TRUE)
  best_sign <- if (vp < vm) -1L else 1L
  best_var  <- min(vp, vm)
  list(var_raw = v0, var_minus = vm, var_plus = vp,
       best_sign = best_sign, var_reduction_pct = 100 * (1 - best_var / v0))
}

partB_rows <- list()
partB_plots <- list()

for (fp in specimen_files) {
  cli::cli_alert_info("specimen: {basename(fp)}")
  d <- .read_trial(fp)
  if (is.null(d)) next
  alpha <- .alpha_rad_s2(d$t_s, d$angle_deg)
  aor <- d$aor_mm; width <- d$width_mm
  spec_moi <- if (is.finite(d$specimen_moi_gmm2)) d$specimen_moi_gmm2 else 0.0

  # Most inertia-dominated run: highest RMS |alpha| AND low |corr(zTorque, angle)|.
  run_start <- c(TRUE, diff(d$t_s) <= 0); run <- cumsum(run_start)
  run_tbl <- purrr::map_dfr(unique(run), function(r) {
    idx <- which(run == r)
    if (length(idx) < 50L) return(NULL)
    zt <- d$ft[idx, 6L]
    ok <- is.finite(zt) & is.finite(d$angle_deg[idx]) & is.finite(alpha[idx])
    if (sum(ok) < 50L) return(NULL)
    cor_ta <- suppressWarnings(abs(stats::cor(zt[ok], d$angle_deg[idx][ok])))
    tibble(run = r, n = length(idx), alpha_rms = sqrt(mean(alpha[idx]^2, na.rm = TRUE)),
           abs_corr_torque_angle = cor_ta)
  })
  if (nrow(run_tbl) == 0L) next
  # inertia-dominated = strong motion, weak torque~angle (elastic) coupling.
  run_tbl$score <- run_tbl$alpha_rms * (1 - pmin(run_tbl$abs_corr_torque_angle, 1))
  best_run <- run_tbl$run[which.max(run_tbl$score)]
  win_idx <- which(run == best_run)

  for (j in seq_len(6L)) {
    ch <- CHANNELS[j]
    I_app_gmm2 <- predict_channel_I_gmm2(ch, aor, width)
    I_gmm2 <- if (is.finite(I_app_gmm2)) I_app_gmm2 else 0.0
    if (ch == "zTorque") I_gmm2 <- I_gmm2 + spec_moi  # specimen bends about z
    I_kg_m2 <- I_gmm2 * kg_m2_per_g_mm2
    raw <- d$ft[win_idx, j]
    ws <- .window_stats(raw, alpha[win_idx], I_kg_m2)
    partB_rows[[length(partB_rows) + 1L]] <- tibble(
      file = basename(fp), channel = ch, aor_mm = aor, width_mm = width,
      I_apparatus_gmm2 = I_app_gmm2, I_total_gmm2 = I_gmm2,
      var_raw = ws$var_raw, var_minus = ws$var_minus, var_plus = ws$var_plus,
      best_sign = ws$best_sign, var_reduction_pct = ws$var_reduction_pct,
      artifact_sign = as.integer(art$correction_slope_sign))
  }

  # figure: per-axis raw vs corrected(best sign) over the chosen window
  sub <- purrr::map_dfr(seq_len(6L), function(j) {
    ch <- CHANNELS[j]
    I_app_gmm2 <- predict_channel_I_gmm2(ch, aor, width)
    I_gmm2 <- if (is.finite(I_app_gmm2)) I_app_gmm2 else 0.0
    if (ch == "zTorque") I_gmm2 <- I_gmm2 + spec_moi
    I_kg_m2 <- I_gmm2 * kg_m2_per_g_mm2
    raw <- d$ft[win_idx, j]
    ws <- .window_stats(raw, alpha[win_idx], I_kg_m2)
    corr <- raw - ws$best_sign * I_kg_m2 * alpha[win_idx]
    t_rel <- d$t_s[win_idx] - min(d$t_s[win_idx], na.rm = TRUE)
    dplyr::bind_rows(
      tibble(t_rel = t_rel, value = raw, series = "raw", channel = ch),
      tibble(t_rel = t_rel, value = corr, series = "corrected", channel = ch))
  })
  p <- ggplot(sub, aes(t_rel, value, color = series)) +
    geom_line(linewidth = 0.4) +
    facet_wrap(~channel, scales = "free_y", ncol = 2) +
    scale_color_manual(values = c(raw = "grey60", corrected = "#b91c1c")) +
    labs(title = paste0("multiaxialinertiacompare: ", basename(fp)),
         subtitle = "Per-axis raw vs (raw - sign*I*alpha), sign auto-chosen on the most inertia-dominated window",
         x = "time within window (s)", y = "force (N) / torque (N*m)") +
    theme_bw(base_size = 9) + theme(legend.position = "bottom")
  partB_plots[[basename(fp)]] <- p
}

partB <- dplyr::bind_rows(partB_rows)
cli::cli_h2("Per-axis correction sign + variance reduction on inertia-dominated window")
print(as.data.frame(partB[, c("file","channel","I_total_gmm2","var_reduction_pct","best_sign","artifact_sign")]),
      row.names = FALSE, digits = 4)

if (length(partB_plots) > 0L) {
  combined <- patchwork::wrap_plots(partB_plots, ncol = 1)
  ggplot2::ggsave(file.path(OUT_DIR, "multiaxialinertiacompare.png"), combined,
                  width = 11, height = 6 * length(partB_plots), dpi = 150, limitsize = FALSE)
  cli::cli_alert_success("Wrote multiaxialinertiacompare.png")
}

cli::cli_alert_success("diag_apparatus_inertia.R complete")
