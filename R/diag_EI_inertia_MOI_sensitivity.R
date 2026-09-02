# diag_EI_inertia_MOI_sensitivity.R
#
# Sensitivity of passive EI-vs-frequency (high-confidence inertial correction
# only: bass16/17/18, July 6 F4 fit at each trial's own AOR/plate + analytic
# I_spec) to a +/-50% mis-estimate of total MOI.
#
# Prefix diagnostic: correction-sensitivity check, not a paper candidate.
#   -> figs_diagnostic/diagnostic_EI_vs_frequency_inertia_MOI_sensitivity.png
#   -> data_processed/diagnostic_EI_inertia_MOI_sensitivity_per_cycle.csv
#
# Run: Rscript R/diag_EI_inertia_MOI_sensitivity.R

suppressPackageStartupMessages({
  library(rhdf5); library(dplyr); library(tibble); library(purrr); library(tidyr)
  library(ggplot2); library(cli); library(readr)
})

.root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.root, f))
src("paths_config.R")
src("00_load_bender_flat.R")
src("03_analyze.R")
src("apparatus_inertia_fit.R")
src("inertia_diag_common.R")

APPARATUS_FIT <- read_apparatus_inertia_fit(
  resolve_apparatus_fit_path("2026-07-06_apparatus_inertia_calibration.json")
)
DENSITY_DEFAULT_G_MM3 <- 0.00106
I_SCALES <- c(0.5, 1.0, 1.5)
GRID_FREQ <- c(1, 2, 2.75, 3, 5, 7)
GRID_CURV <- c(2, 3.36, 4, 4.47, 5.59, 6, 8)

OUT_FIG  <- file.path(FIGS_DIAGNOSTIC_DIR, "diagnostic_EI_vs_frequency_inertia_MOI_sensitivity.png")
OUT_CSV  <- file.path(.crittergripper_root(), "02_processed", "data_processed",
                      "diagnostic_EI_inertia_MOI_sensitivity_per_cycle.csv")
dir.create(FIGS_DIAGNOSTIC_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)

ARCHIVE <- .permanent_archive_root()
dyn_files <- c(
  file.path(ARCHIVE, "2026-07-14_bass16_bender",
            paste0("2026-07-14_bass16_bender_", c("17", "18", "19"), "_dynamic.h5")),
  file.path(ARCHIVE, "2026-07-15_bass17_bender",
            paste0("2026-07-15_bass17_bender_", c("17", "18", "19", "20"), "_dynamic.h5")),
  file.path(ARCHIVE, "2026-07-16_bass18_bender",
            paste0("2026-07-16_bass18_bender_", c("14", "15"), "_dynamic.h5"))
)
dyn_files <- dyn_files[file.exists(dyn_files)]

.snap_to_grid <- function(x, grid, abs_tol = 0.15) {
  vapply(x, function(v) {
    if (!is.finite(v) || length(grid) == 0L) return(NA_real_)
    d <- abs(grid - v)
    j <- which.min(d)
    if (d[j] <= abs_tol) grid[j] else NA_real_
  }, numeric(1L))
}

.dbl1_attr <- function(a, k) {
  v <- suppressWarnings(as.numeric(a[[k]][1L]))
  if (length(v) == 0L) NA_real_ else v
}

.read_inertia_geom <- function(path) {
  a <- tryCatch(h5readAttributes(path, "/metadata"), error = function(e) list())
  list(
    width_mm  = .dbl1_attr(a, "measurement_specimen_local_body_width_millimeter"),
    height_mm = .dbl1_attr(a, "measurement_specimen_local_body_height_millimeter"),
    clamp_mm  = .dbl1_attr(a, "measurement_clamp_separation_millimeter"),
    density_g_mm3 = .dbl1_attr(a, "measurement_specimen_density_gram_per_cubic_millimeter"),
    aor_mm    = .dbl1_attr(a, "calibration_inertia_apparatus_aor_to_clamp_millimeter"),
    plate_mm  = .dbl1_attr(a, "calibration_inertia_apparatus_plate_to_plate_millimeter")
  )
}

.specimen_inertia_lateral_kg_m2 <- function(width_mm, height_mm, clamp_mm, density_g_mm3) {
  bY <- width_mm / 2.0; bZ <- height_mm / 2.0; L <- clamp_mm; rho <- density_g_mm3
  if (!all(is.finite(c(bY, bZ, L, rho))) || rho <= 0) return(NA_real_)
  m_g <- rho * pi * bY * bZ * L
  (m_g * (bZ^2 / 4 + L^2 / 12)) * kg_m2_per_g_mm2
}

.I_total_signed <- function(path) {
  geom <- .read_inertia_geom(path)
  stopifnot(is.finite(geom$aor_mm), is.finite(geom$plate_mm))
  app <- apparatus_inertia_from_fit(APPARATUS_FIT, geom$aor_mm, geom$plate_mm)
  dens <- geom$density_g_mm3
  if (!is.finite(dens)) dens <- DENSITY_DEFAULT_G_MM3
  I_spec <- .specimen_inertia_lateral_kg_m2(geom$width_mm, geom$height_mm, geom$clamp_mm, dens)
  if (!is.finite(I_spec)) stop("I_spec missing: ", basename(path))
  I_app_signed <- app$i_gmm2 * kg_m2_per_g_mm2 * APPARATUS_FIT$correction_slope_sign
  I_app_signed + sign(I_app_signed) * I_spec
}

.tau_col <- function(td) {
  hit <- intersect(c("ztorque.Nm", "ztorque0.Nm", "tau_Nm"), names(td))
  if (length(hit) == 0L) NA_character_ else hit[[1L]]
}

process_file <- function(path) {
  cli::cli_alert_info("{basename(path)}")
  td <- load_bender_flat(path, do_filter = TRUE, loadtorques = "z", loadforces = FALSE)
  if (is.null(td) || nrow(td) == 0L) return(tibble::tibble())
  tq <- .tau_col(td)
  if (is.na(tq)) return(tibble::tibble())
  I_tot <- .I_total_signed(path)
  ang <- if ("enc.deg" %in% names(td) && any(is.finite(td$enc.deg))) td$enc.deg else td$angle.deg
  td$alpha <- .kinematics(td$t.s, ang)$alpha_rad_s2
  td$tau_raw <- td[[tq]]
  fish <- unique(stats::na.omit(td$fishcode))[1L]
  td <- set_cycle_types(td) |>
    dplyr::filter(.data$cycletype == "pass", is.finite(.data$cycle), .data$cycle >= 0,
                  is.finite(.data$tau_raw), is.finite(.data$curve.invm), is.finite(.data$alpha))
  if (nrow(td) < 10L) return(tibble::tibble())

  purrr::map_dfr(I_SCALES, function(k) {
    tau_k <- td$tau_raw - k * I_tot * td$alpha
    td$tau_k <- tau_k
    td |>
      dplyr::group_by(.data$cycle) |>
      dplyr::summarise(
        n = dplyr::n(),
        period_s = diff(range(.data$t.s, na.rm = TRUE)),
        curv_amp = (max(.data$curve.invm, na.rm = TRUE) - min(.data$curve.invm, na.rm = TRUE)) / 2,
        EI_Nm2 = tryCatch(calc_avg_stiffness(.data$curve.invm, .data$tau_k), error = function(e) NA_real_),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        fishcode = fish, filename = basename(path),
        I_scale = k, I_total_signed_kgm2 = I_tot,
        freq_meas = dplyr::if_else(.data$period_s > 0, 1.0 / .data$period_s, NA_real_),
        freq_nominal = .snap_to_grid(.data$freq_meas, GRID_FREQ, abs_tol = 0.12),
        curv_nominal = .snap_to_grid(.data$curv_amp, GRID_CURV, abs_tol = 0.25)
      ) |>
      dplyr::filter(is.finite(.data$EI_Nm2), .data$EI_Nm2 != 0,
                    is.finite(.data$freq_nominal), is.finite(.data$curv_nominal),
                    .data$freq_nominal >= 1, .data$n >= 10L)
  })
}

cli::cli_h1("MOI sensitivity +/-50% on {length(dyn_files)} high-confidence dynamic files")
rows <- purrr::map_dfr(dyn_files, process_file)

# Keep only cycles present at all three scales, then drop singleton cells.
keep <- rows |>
  dplyr::count(.data$fishcode, .data$filename, .data$cycle) |>
  dplyr::filter(.data$n == length(I_SCALES)) |>
  dplyr::select(-"n")
rows <- dplyr::semi_join(rows, keep, by = c("fishcode", "filename", "cycle"))
rows <- rows |>
  dplyr::group_by(.data$fishcode, .data$freq_nominal, .data$curv_nominal, .data$I_scale) |>
  dplyr::filter(dplyr::n() >= 3L) |>
  dplyr::ungroup()

readr::write_csv(rows, OUT_CSV)
cli::cli_alert_success("Wrote {nrow(rows)} rows -> {OUT_CSV}")

tab <- rows |>
  dplyr::group_by(.data$freq_nominal, .data$I_scale) |>
  dplyr::summarise(EI_med = median(abs(.data$EI_Nm2)), n = dplyr::n(), .groups = "drop") |>
  tidyr::pivot_wider(names_from = "I_scale", values_from = c("EI_med", "n"))
# pivot names: EI_med_0.5 may become EI_med_0.5
cli::cli_h2("Median |EI| (N m^2), kappa pooled")
print(as.data.frame(tab), row.names = FALSE, digits = 3)

base <- rows |>
  dplyr::filter(.data$I_scale == 1) |>
  dplyr::select(.data$fishcode, .data$filename, .data$cycle, EI_1 = .data$EI_Nm2)
cmp <- rows |>
  dplyr::filter(.data$I_scale != 1) |>
  dplyr::inner_join(base, by = c("fishcode", "filename", "cycle")) |>
  dplyr::mutate(pct = 100 * (abs(.data$EI_Nm2) - abs(.data$EI_1)) / abs(.data$EI_1))
cli::cli_h2("Percent change vs 1.0x MOI (median, kappa pooled)")
print(as.data.frame(
  cmp |>
    dplyr::group_by(.data$freq_nominal, .data$I_scale) |>
    dplyr::summarise(pct_med = median(.data$pct, na.rm = TRUE),
                     pct_q10 = stats::quantile(.data$pct, 0.1, na.rm = TRUE),
                     pct_q90 = stats::quantile(.data$pct, 0.9, na.rm = TRUE),
                     .groups = "drop")
), row.names = FALSE, digits = 2)

rows$I_lab <- factor(
  sprintf("%.1fx MOI", rows$I_scale),
  levels = sprintf("%.1fx MOI", I_SCALES)
)
rows$kap_lab <- paste0("kappa = ", rows$curv_nominal, " /m")

p <- ggplot(rows, aes(.data$freq_nominal, abs(.data$EI_Nm2), color = .data$I_lab)) +
  geom_point(size = 1.3, alpha = 0.22, position = position_jitter(width = 0.07, height = 0)) +
  stat_summary(fun = mean, geom = "line", aes(group = .data$I_lab), linewidth = 0.85) +
  stat_summary(fun = mean, geom = "point", size = 2.3, shape = 21, fill = "white") +
  facet_wrap(~ kap_lab, nrow = 2) +
  scale_color_manual(values = c("0.5x MOI" = "#4daf4a", "1.0x MOI" = "#377eb8", "1.5x MOI" = "#e41a1c"),
                     name = "I_total scale") +
  labs(title = "Passive |EI| vs frequency: +/-50% total MOI",
       subtitle = "bass16/17/18 only; July 6 I_app at trial geometry + analytic I_spec; 0.5x / 1.0x / 1.5x I_total",
       x = "Frequency (Hz)", y = "|EI| (N m^2)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank())

ggplot2::ggsave(OUT_FIG, p, width = 10, height = 7.5, dpi = 150)
cli::cli_alert_success("Wrote {OUT_FIG}")
cli::cli_h1("diag_EI_inertia_MOI_sensitivity.R complete")
