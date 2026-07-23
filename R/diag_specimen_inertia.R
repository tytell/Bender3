# diag_specimen_inertia.R
# READ-ONLY DIAGNOSTIC (2026-07-22, PI-directed follow-up to diag_apparatus_
# inertia.R). Question from the PI: run "a similar analysis on specimen inertia
# using lxwxh frustum dimensions in metadata and provided tissue density."
#
# WHY: diag_apparatus_inertia.R showed the APPARATUS carries almost no yTorque
# inertia (R^2 of yTorque-vs-alpha ~0.006), so the yTorque residual that
# ytorqueinertialtiming found is NOT apparatus inertia. The two remaining
# candidates are (a) SPECIMEN inertia (the clamped body mass, offset dorso-
# ventrally by the clamp offset d), and (b) the active-vs-passive acceleration-
# profile mismatch. This script quantifies (a) directly.
#
# Two complementary estimates of specimen inertia:
#   PART 1 (ANALYTIC): from the body cross-section (local_body_height x
#     local_body_width), the bent span (clamp_separation) and the tissue density
#     in metadata, model the segment as a solid elliptical cylinder and compute
#     its full inertia tensor about the SENSOR ORIGIN (parallel-axis with the
#     dorsoventral clamp offset d). Validate the bending term (I_zz) against the
#     stored calibration_inertia_specimen_moi. Crucially, the yTorque term is
#     the product of inertia I_yz = -m*y_cm*d, which is ZERO for a mediolaterally
#     centered specimen -- so the analysis reports how large a mediolateral CoM
#     offset y_cm would be needed to explain the observed yTorque, rather than
#     assuming one.
#   PART 2 (EMPIRICAL): on the specimen's own PASSIVE (no-stim) ramps, regress
#     each channel on [angle + alpha] to SEPARATE the elastic term (angle) from
#     the inertial term (alpha). The alpha coefficient is the TOTAL (apparatus +
#     specimen) per-channel inertia measured on the real loaded body; subtract
#     the apparatus prediction (from the 11-trial corpus fit) to isolate the
#     SPECIMEN contribution per channel. For yTorque this is the direct test:
#     is there an alpha-linear term beyond the apparatus, and does removing it
#     flatten the ramp-correlated residual?
#
# One token: specimeninertiacompare. Changes nothing in the production path.
#
# Run: Rscript R/diag_specimen_inertia.R

suppressPackageStartupMessages({
  library(rhdf5); library(dplyr); library(tibble); library(purrr)
  library(ggplot2); library(patchwork); library(cli)
})

.root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.root, f))
src("paths_config.R")
src("apparatus_inertia_fit.R")
src("inertia_diag_common.R")

OUT_DIR <- FIGS_DIAGNOSTIC_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

APPARATUS_SUBFOLDER <- Sys.getenv("BENDER3_APPARATUS_CALIB_SUBFOLDER",
                                  "2026-07-06_apparatus_inertial_calibration")
FIT_JSON <- Sys.getenv("BENDER3_APPARATUS_FIT_JSON",
                       "2026-07-06_apparatus_inertia_calibration.json")


# =============================================================================
# PART 0: apparatus per-channel inertia model (for empirical subtraction)
# =============================================================================
cli::cli_h1("PART 0: apparatus per-channel inertia model (from corpus)")

art <- read_apparatus_inertia_fit(resolve_apparatus_fit_path(FIT_JSON))
corpus_dir <- raw_source_dir(APPARATUS_SUBFOLDER)
excluded <- unlist(art$excluded_files)

corpus <- purrr::map_dfr(art$trials, function(tr) {
  fp <- file.path(corpus_dir, tr$file)
  if (!file.exists(fp) || tr$file %in% excluded) return(NULL)
  d <- .read_trial(fp)  # no-stim cache key: reuses diag_apparatus_inertia.R's cache
  if (is.null(d)) return(NULL)
  alpha <- .alpha_rad_s2(d$t_s, d$angle_deg)
  pc <- .per_channel_inertia(d$ft, alpha)
  pc$aor_mm <- as.numeric(tr$aor_millimeter); pc$width_mm <- as.numeric(tr$width_millimeter)
  pc
})

# Per-channel F4 geometry model of I (g*mm^2) + the channel's typical slope sign.
channel_app_model <- setNames(lapply(CHANNELS, function(ch) {
  dch <- corpus |> dplyr::filter(channel == ch, is.finite(I_gmm2))
  if (nrow(dch) < 4L) return(NULL)
  list(model = stats::lm(I_gmm2 ~ I(aor_mm^2) + I(width_mm^2), data = dch),
       sign = stats::median(dch$slope_sign, na.rm = TRUE))
}), CHANNELS)

#' Signed apparatus inertia (kg*m^2) predicted for a channel at (aor, width).
apparatus_inertia_signed_kg_m2 <- function(ch, aor_mm, width_mm) {
  m <- channel_app_model[[ch]]
  if (is.null(m)) return(NA_real_)
  i_gmm2 <- as.numeric(stats::predict(m$model, newdata = data.frame(aor_mm = aor_mm, width_mm = width_mm)))
  i_gmm2 * kg_m2_per_g_mm2 * m$sign
}


# =============================================================================
# PART 1: ANALYTIC specimen inertia tensor from body geometry + density
# =============================================================================
cli::cli_h1("PART 1: analytic specimen inertia (elliptical cylinder from metadata)")

specimen_subfolders <- list(bass16 = BASS16_RAW_SUBFOLDER, bass17 = BASS17_RAW_SUBFOLDER,
                            bass18 = BASS18_RAW_SUBFOLDER)

# One representative isovelocity trial per specimen (has geometry + moving ramps).
specimen_trials <- purrr::imap(specimen_subfolders, function(sub, sid) {
  fps <- list.files(raw_source_dir(sub), pattern = "isovelocity\\.h5$", full.names = TRUE)
  if (length(fps) == 0L) return(NULL)
  list(sid = sid, fp = fps[[1L]])
}) |> purrr::compact()

#' Solid elliptical-cylinder inertia tensor (about centroid, then sensor origin).
#' Sensor frame: X longitudinal (length), Y mediolateral (body width), Z
#' dorsoventral (body height). Bending is about Z. r_cm offset = (0, y_cm, d).
specimen_tensor <- function(h_mm, w_mm, L_mm, rho, d_mm, y_cm = 0.0) {
  bY <- w_mm / 2.0   # mediolateral semi-axis
  bZ <- h_mm / 2.0   # dorsoventral semi-axis
  m  <- rho * pi * bY * bZ * L_mm            # g
  I_zz_c <- m * (bY^2 / 4 + L_mm^2 / 12)     # bend axis (dorsoventral)
  I_yy_c <- m * (bZ^2 / 4 + L_mm^2 / 12)
  I_xx_c <- m * (bY^2 + bZ^2) / 4
  list(mass_g = m, bY = bY, bZ = bZ,
       I_xx_c = I_xx_c, I_yy_c = I_yy_c, I_zz_c = I_zz_c,
       I_zz_o = I_zz_c + m * y_cm^2,          # zTorque (adds specimen MOI)
       I_yy_o = I_yy_c + m * d_mm^2,
       I_yz_o = -m * y_cm * d_mm)             # yTorque product of inertia
}

part1 <- purrr::map_dfr(specimen_trials, function(tt) {
  cli::cli_alert_info("specimen geometry: {tt$sid} ({basename(tt$fp)})")
  d <- .read_trial(tt$fp, with_stim = TRUE)
  if (is.null(d)) return(NULL)
  ten <- specimen_tensor(d$body_height_mm, d$body_width_mm, d$clamp_sep_mm,
                         d$density_g_mm3, d$clamp_offset_vert_mm, y_cm = 0.0)
  tibble(specimen = tt$sid,
         body_height_mm = d$body_height_mm, body_width_mm = d$body_width_mm,
         clamp_sep_mm = d$clamp_sep_mm, density_g_mm3 = d$density_g_mm3,
         d_vert_mm = d$clamp_offset_vert_mm,
         mass_g = ten$mass_g,
         I_zz_analytic_gmm2 = ten$I_zz_o, I_zz_stored_gmm2 = d$specimen_moi_gmm2,
         I_yy_analytic_gmm2 = ten$I_yy_o, I_xx_analytic_gmm2 = ten$I_xx_c,
         I_yz_per_ycm = -ten$mass_g * d$clamp_offset_vert_mm)  # I_yz per 1 mm of y_cm
})

cli::cli_h2("Analytic specimen inertia vs stored bending MOI (g*mm^2)")
print(as.data.frame(part1[, c("specimen","mass_g","I_zz_analytic_gmm2","I_zz_stored_gmm2",
                              "I_yy_analytic_gmm2","I_xx_analytic_gmm2")]),
      row.names = FALSE, digits = 5)
cli::cli_alert_info(paste0(
  "yTorque product-of-inertia I_yz = -m*y_cm*d is ZERO for a centered specimen; ",
  "per 1 mm of mediolateral CoM offset it is {round(mean(part1$I_yz_per_ycm),0)} g*mm^2 (bass mean)."))


# =============================================================================
# PART 2: EMPIRICAL specimen inertia from passive ramps ([angle + alpha])
# =============================================================================
cli::cli_h1("PART 2: empirical specimen inertia from passive no-stim ramps")

part2_rows <- list(); part2_ytrace <- list()

for (tt in specimen_trials) {
  d <- .read_trial(tt$fp, with_stim = TRUE)
  if (is.null(d)) next
  kin <- .kinematics(d$t_s, d$angle_deg)
  aor <- d$aor_mm; width <- d$width_mm

  # PASSIVE + MOVING samples only: no stim, and |velocity| above a small floor so
  # the ramps (which carry alpha) dominate over motionless holds.
  passive <- !d$active & is.finite(kin$alpha_rad_s2) & abs(kin$vel_deg_s) > 5
  if (sum(passive, na.rm = TRUE) < 200L) {
    cli::cli_alert_warning("{tt$sid}: too few passive moving samples ({sum(passive, na.rm=TRUE)})")
    next
  }

  emp <- .per_channel_inertia_vs_elastic(d$ft[passive, , drop = FALSE],
                                         d$angle_deg[passive], kin$alpha_rad_s2[passive])
  emp$specimen <- tt$sid
  emp$apparatus_kg_m2 <- vapply(emp$channel, function(ch) apparatus_inertia_signed_kg_m2(ch, aor, width), numeric(1L))
  emp$specimen_kg_m2  <- emp$inertia_kg_m2 - emp$apparatus_kg_m2
  emp$total_I_gmm2 <- emp$inertia_kg_m2 * g_mm2_per_kg_m2
  emp$apparatus_I_gmm2 <- emp$apparatus_kg_m2 * g_mm2_per_kg_m2
  emp$specimen_I_gmm2  <- emp$specimen_kg_m2 * g_mm2_per_kg_m2
  part2_rows[[tt$sid]] <- emp

  # yTorque trace test on the most inertia-dominated passive run: raw vs
  # apparatus-corrected vs total(empirical)-corrected -- does the ramp feature go?
  run_start <- c(TRUE, diff(d$t_s) <= 0); run <- cumsum(run_start)
  run_tbl <- purrr::map_dfr(unique(run), function(r) {
    idx <- which(run == r & passive)
    if (length(idx) < 50L) return(NULL)
    tibble(run = r, arms = sqrt(mean(kin$alpha_rad_s2[idx]^2, na.rm = TRUE)), n = length(idx))
  })
  if (nrow(run_tbl) > 0L) {
    best <- run_tbl$run[which.max(run_tbl$arms)]
    idx <- which(run == best)
    yt <- d$ft[idx, 5L]  # yTorque
    a  <- kin$alpha_rad_s2[idx]
    I_app <- apparatus_inertia_signed_kg_m2("yTorque", aor, width)
    I_tot <- emp$inertia_kg_m2[emp$channel == "yTorque"]
    t_rel <- d$t_s[idx] - min(d$t_s[idx], na.rm = TRUE)
    part2_ytrace[[tt$sid]] <- dplyr::bind_rows(
      tibble(specimen = tt$sid, t_rel = t_rel, value = yt, series = "raw"),
      tibble(specimen = tt$sid, t_rel = t_rel, value = yt - I_app * a, series = "apparatus-corrected"),
      tibble(specimen = tt$sid, t_rel = t_rel, value = yt - I_tot * a, series = "empirical-total-corrected"))
  }
}

part2 <- dplyr::bind_rows(part2_rows)
cli::cli_h2("Empirical per-channel inertia on passive ramps: total vs apparatus vs specimen (g*mm^2)")
print(as.data.frame(part2[, c("specimen","channel","total_I_gmm2","apparatus_I_gmm2",
                              "specimen_I_gmm2","alpha_partial_r2")]),
      row.names = FALSE, digits = 4)

# Implied mediolateral CoM offset that would produce the empirical yTorque specimen inertia.
yt_rows <- part2 |> dplyr::filter(channel == "yTorque")
if (nrow(yt_rows) > 0L) {
  yt_rows <- yt_rows |> dplyr::left_join(
    part1[, c("specimen", "I_yz_per_ycm")], by = "specimen") |>
    dplyr::mutate(implied_y_cm_mm = specimen_I_gmm2 / I_yz_per_ycm)
  cli::cli_h2("yTorque: implied mediolateral CoM offset to explain the empirical specimen term")
  print(as.data.frame(yt_rows[, c("specimen","specimen_I_gmm2","alpha_partial_r2","implied_y_cm_mm")]),
        row.names = FALSE, digits = 4)
}


# =============================================================================
# Figure
# =============================================================================
p1 <- ggplot(part2, aes(channel)) +
  geom_col(aes(y = abs(total_I_gmm2), fill = "total (app+spec)"), width = 0.28, position = position_nudge(x = -0.28)) +
  geom_col(aes(y = abs(apparatus_I_gmm2), fill = "apparatus"), width = 0.28) +
  geom_col(aes(y = abs(specimen_I_gmm2), fill = "specimen (=total-app)"), width = 0.28, position = position_nudge(x = 0.28)) +
  facet_wrap(~specimen, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("total (app+spec)" = "grey40", "apparatus" = "#0369a1", "specimen (=total-app)" = "#b91c1c"), name = NULL) +
  labs(title = "(a) Empirical per-channel inertia on passive ramps (|I|, g*mm^2)",
       subtitle = "specimen = total (measured on loaded body) minus apparatus (corpus fit)",
       x = "F/T channel", y = "|inertia| (g*mm^2)") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")

pY <- if (length(part2_ytrace) > 0L) {
  yd <- dplyr::bind_rows(part2_ytrace)
  ggplot(yd, aes(t_rel, value, color = series)) +
    geom_line(linewidth = 0.4) +
    facet_wrap(~specimen, ncol = 1, scales = "free") +
    scale_color_manual(values = c(raw = "grey60", `apparatus-corrected` = "#0369a1",
                                  `empirical-total-corrected` = "#b91c1c")) +
    labs(title = "(b) yTorque on the most inertia-dominated passive ramp",
         subtitle = "Does removing an I*alpha term flatten the ramp-correlated residual? apparatus alone vs empirical total",
         x = "time within ramp (s)", y = "yTorque (N*m)") +
    theme_bw(base_size = 9) + theme(legend.position = "bottom")
} else NULL

fig <- if (!is.null(pY)) (p1 | pY) else p1
fig <- fig + patchwork::plot_annotation(
  title = "specimeninertiacompare: is the yTorque residual specimen inertia?",
  subtitle = "Analytic tensor (Part 1, printed) + empirical passive-ramp separation (Part 2). PRE-WIRING diagnostic.")
ggplot2::ggsave(file.path(OUT_DIR, "specimeninertiacompare.png"), fig,
                width = 13, height = 10, dpi = 150, limitsize = FALSE)
cli::cli_alert_success("Wrote specimeninertiacompare.png")
cli::cli_alert_success("diag_specimen_inertia.R complete")
