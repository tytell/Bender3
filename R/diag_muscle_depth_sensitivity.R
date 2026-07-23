# diag_muscle_depth_sensitivity.R
# Diagnostic (not a pipeline deliverable): tests whether mis-estimating
# muscle depth (measurement_target_muscle_depth_millimeter, currently 0/
# unset in every bass16/17/18 file -> DEFAULT_MUSCLE_DEPTH_MM = 1.0 mm
# fallback, muscle_geometry.R) could explain why measured (sono) strain
# systematically EXCEEDS encoder-predicted strain in dynamic trials (see
# diag_sono_amp_ratio_multitrial.R / diag_sono_amp_ratio_destaircased.R).
#
# Geometry: predicted strain = curvature * r_m, r_m = (half_body_width -
# muscle_depth) / 1000 [[meters]] (compute_predicted_strain(),
# muscle_geometry.R). Predicted strain is therefore DIRECTLY PROPORTIONAL
# to r_m -- so the "gain" needed to make predicted strain catch up to
# measured strain is exactly r_m_needed / r_m_current, and the muscle depth
# that would produce it is depth_needed = half_width - r_m_needed.
#
# Since r_m INCREASES as depth DECREASES (muscle closer to the surface,
# farther from the neutral mid-width axis), the depth has a hard physical
# ceiling at 0 mm (muscle exactly at the surface) -- there is no shallower
# option. That ceiling bounds the MAXIMUM possible predicted-strain gain
# achievable by adjusting depth alone, independent of what the true depth
# turns out to be.
#
# Run with:  Rscript R/diag_muscle_depth_sensitivity.R
# Outputs -> figures: 02_processed/figs_diagnostic/ (FIGS_DIAGNOSTIC_DIR)
#            data:    02_processed/data_processed/

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2); library(cli)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")
src("muscle_geometry.R")  # DEFAULT_MUSCLE_DEPTH_MM

OUT_DIR      <- FIGS_DIAGNOSTIC_DIR
DATA_OUT_DIR <- file.path(.crittergripper_root(), "02_processed", "data_processed")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Confirmed identical across every archived bass16/17/18 trial (all
# protocol families) -- see diag_sono_amp_ratio_multitrial.R's morphometrics
# audit. Body width differs by specimen; depth assumption is rig-wide.
BODY_WIDTH_MM <- c(bass16 = 42.0, bass17 = 38.0, bass18 = 44.0)
DEPTH_CURRENT_MM <- DEFAULT_MUSCLE_DEPTH_MM  # 1.0 mm, from muscle_geometry.R

# Observed measured/predicted amplitude gains this sensitivity analysis is
# checking against (from diag_sono_amp_ratio_multitrial.R /
# diag_sono_amp_ratio_destaircased.R's no-stim linear fits and per-bin
# ratios -- both raw-staircase and de-staircased, since neither changed the
# slope/gain component).
TARGET_GAINS <- c(
  "De-staircased linear-fit slope"      = 1.15,
  "Raw high-amplitude bin (20-40%)"     = 1.32,
  "Raw dominant bin (5-10%, most data)" = 1.64,
  "Raw pooled overall"                  = 1.71
)

# =============================================================================
# 1) What's the MAXIMUM gain physically achievable by varying depth, given
#    depth cannot go below 0 mm (muscle at the surface)?
# =============================================================================
half_width_mm  <- BODY_WIDTH_MM / 2.0
r_m_current_mm <- half_width_mm - DEPTH_CURRENT_MM
r_m_ceiling_mm <- half_width_mm - 0.0  # depth = 0, the physical limit
max_gain_achievable <- r_m_ceiling_mm / r_m_current_mm

ceiling_tbl <- tibble(
  specimen = names(BODY_WIDTH_MM),
  half_width_mm = half_width_mm,
  depth_current_assumed_mm = DEPTH_CURRENT_MM,
  r_m_current_mm = r_m_current_mm,
  r_m_ceiling_mm_at_depth0 = r_m_ceiling_mm,
  max_gain_achievable_by_depth_alone = round(max_gain_achievable, 4)
)
cli::cli_h1("Physical ceiling: max predicted-strain gain achievable via depth alone (depth >= 0)")
print(ceiling_tbl)
write.csv(ceiling_tbl, file.path(DATA_OUT_DIR, "muscle_depth_sensitivity_ceiling.csv"), row.names = FALSE)

# =============================================================================
# 2) What depth WOULD be required to hit each observed gain target?
# =============================================================================
req_rows <- list()
for (tn in names(TARGET_GAINS)) {
  g <- TARGET_GAINS[[tn]]
  r_m_needed_mm   <- r_m_current_mm * g
  depth_needed_mm <- half_width_mm - r_m_needed_mm
  req_rows[[tn]] <- tibble(
    specimen = names(BODY_WIDTH_MM), target = tn, target_gain = g,
    r_m_needed_mm = round(r_m_needed_mm, 2), depth_needed_mm = round(depth_needed_mm, 2),
    physically_possible = depth_needed_mm >= 0
  )
}
req_tbl <- dplyr::bind_rows(req_rows)
cli::cli_h1("Depth that WOULD be required to fully explain each observed gain")
print(req_tbl, n = 100)
write.csv(req_tbl, file.path(DATA_OUT_DIR, "muscle_depth_sensitivity_required_depth.csv"), row.names = FALSE)

# =============================================================================
# Plot: gain-vs-depth curve per specimen, target gains as reference lines,
#       physical ceiling (depth=0) marked.
# =============================================================================
depth_grid <- seq(-15, 10, by = 0.25)
curve_df <- purrr::map_dfr(names(BODY_WIDTH_MM), function(sp) {
  hw <- half_width_mm[[sp]]; r_cur <- r_m_current_mm[[sp]]
  tibble(specimen = sp, depth_mm = depth_grid, r_m_mm = hw - depth_grid, gain = (hw - depth_grid) / r_cur)
})

p <- ggplot(curve_df, aes(x = .data$depth_mm, y = .data$gain, color = .data$specimen)) +
  annotate("rect", xmin = -15, xmax = 0, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 1.0, linetype = "solid", color = "grey40") +
  geom_hline(data = tibble(target = names(TARGET_GAINS), gain = unname(TARGET_GAINS)),
             aes(yintercept = .data$gain), linetype = "dashed", color = "black", inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  geom_vline(xintercept = DEPTH_CURRENT_MM, linetype = "dotted", color = "#b91c1c") +
  annotate("text", x = -15, y = max(unlist(TARGET_GAINS)) + 0.05, hjust = 0, vjust = 0,
           label = "physically impossible\n(depth < 0, muscle outside body)", size = 3.0, color = "grey30") +
  annotate("text", x = 0.3, y = 1.02, hjust = 0, label = "depth = 0\n(surface, ceiling)", size = 2.8) +
  annotate("text", x = DEPTH_CURRENT_MM + 0.3, y = 1.0, hjust = 0, vjust = 1.3,
           label = "current\nassumption\n(1 mm)", size = 2.8, color = "#b91c1c") +
  scale_y_continuous(breaks = seq(0.8, 1.8, 0.1)) +
  coord_cartesian(xlim = c(-15, 10), ylim = c(0.85, 1.8)) +
  labs(title = "Muscle-depth sensitivity: can mis-estimated muscle depth explain measured > predicted strain?",
       subtitle = paste0(
         "Predicted strain is proportional to r_m = half-body-width - muscle depth, so gain = r_m(depth) / r_m(1mm, current assumption).\n",
         "Dashed black lines = OBSERVED measured/predicted amplitude gains (de-staircased slope 1.15 up to raw pooled 1.71).\n",
         "Grey band = physically impossible (depth < 0, i.e. muscle would sit outside the body surface)."),
       x = "Assumed muscle depth (mm)", y = "Predicted-strain gain relative to current 1mm assumption",
       color = "Specimen") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

fout <- file.path(OUT_DIR, "muscle_depth_sensitivity.png")
ggplot2::ggsave(fout, p, width = 9.5, height = 6.5, dpi = 150)
cli::cli_alert_success("Saved {fout}")
cli::cli_alert_success("diag_muscle_depth_sensitivity.R complete")
