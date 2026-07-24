# summary_precondition_power_tension_earlyVsLater.R
# PI question, 2026-07-24: "I agree on the exclusion, but doesn't that
# leave me with very few muscle power (W/kg) and tension data? Add a
# mass-specific power and tension plot to figs_summary, making sure to
# compare the outputs of the early vs later groups."
#
# Pools the trial-level power/tension tables that the three protocol-
# specific diagnostics already computed (no raw H5 reprocessing) and shows,
# side by side, exactly how many trials survive the early (preconditioning)
# vs. later (stable) split in each protocol family, plus whether the
# power/tension values themselves differ between the two groups.
#
# Reuses:
#   dynamic_precondition_power_vs_offset_by_trial.csv     (diag_precondition_
#     power_vs_offset.R)          -- has precondition already
#   isovelocity_power_vs_offset_by_trial.csv               (diag_precondition_
#     power_vs_offset_isovelocity.R) -- trial_num only, classify here
#   isometric_tension_vs_offset_by_trial.csv               (diag_precondition_
#     tension_vs_offset_isometric.R) -- trial_num only, classify here
#
# Run with:  Rscript R/summary_precondition_power_tension_earlyVsLater.R
# Outputs -> figures: 02_processed/figs_summary/ (FIGS_SUMMARY_DIR)
#            data:    02_processed/data_processed/

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(ggplot2); library(cli); library(readr)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")
src("dynamic_trial_precondition.R")

OUT_DIR      <- FIGS_SUMMARY_DIR
DATA_OUT_DIR <- file.path(.crittergripper_root(), "02_processed", "data_processed")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SPECIMEN_COLORS <- c(bass16 = "#1b9e77", bass17 = "#d95f02", bass18 = "#7570b3")

# =============================================================================
# Load the three trial-level tables and tidy into one long data frame:
# specimen, trial_id, protocol_family, precondition, metric, value
# =============================================================================
dyn <- readr::read_csv(file.path(DATA_OUT_DIR, "dynamic_precondition_power_vs_offset_by_trial.csv"), show_col_types = FALSE) |>
  dplyr::transmute(.data$specimen, .data$trial_id, .data$trial_num, .data$precondition,
                    protocol_family = "dynamic",
                    `mean power (W/kg)` = .data$mean_avg_power_Wkg,
                    `max power (W/kg)`  = .data$max_peak_power_Wkg)

isovel <- readr::read_csv(file.path(DATA_OUT_DIR, "isovelocity_power_vs_offset_by_trial.csv"), show_col_types = FALSE) |>
  dplyr::mutate(precondition = classify_session_precondition(.data$specimen, .data$trial_num)) |>
  dplyr::transmute(.data$specimen, .data$trial_id, .data$trial_num, .data$precondition,
                    protocol_family = "isovelocity",
                    `mean power (W/kg)` = .data$mean_power_Wkg,
                    `max power (W/kg)`  = .data$max_power_Wkg)

isomet <- readr::read_csv(file.path(DATA_OUT_DIR, "isometric_tension_vs_offset_by_trial.csv"), show_col_types = FALSE) |>
  dplyr::mutate(precondition = classify_session_precondition(.data$specimen, .data$trial_num)) |>
  dplyr::transmute(.data$specimen, .data$trial_id, .data$trial_num, .data$precondition,
                    protocol_family = "isometric",
                    `mean tension (N/cm^2)` = .data$mean_tension_Ncm2,
                    `max tension (N/cm^2)`  = .data$max_tension_Ncm2)

pooled <- dplyr::bind_rows(dyn, isovel, isomet) |>
  tidyr::pivot_longer(cols = c(dplyr::starts_with("mean "), dplyr::starts_with("max ")),
                       names_to = "metric", values_to = "value") |>
  dplyr::filter(is.finite(.data$value), !is.na(.data$precondition)) |>
  dplyr::mutate(
    stat        = ifelse(startsWith(.data$metric, "mean"), "mean (trial-average)", "max (trial-peak)"),
    panel_label = sprintf("%s: %s", tools::toTitleCase(.data$protocol_family), .data$metric)
  )

write.csv(pooled, file.path(DATA_OUT_DIR, "precondition_power_tension_earlyVsLater_pooled.csv"), row.names = FALSE)

# =============================================================================
# Sample-size table -- directly answers "doesn't exclusion leave very few
# data points?"
# =============================================================================
n_tbl <- pooled |>
  dplyr::distinct(.data$specimen, .data$trial_id, .data$protocol_family, .data$precondition) |>
  dplyr::count(.data$protocol_family, .data$precondition, name = "n_trials") |>
  tidyr::pivot_wider(names_from = "precondition", values_from = "n_trials", values_fill = 0)
cli::cli_h1("Trial counts surviving early/later split, by protocol")
print(n_tbl)
write.csv(n_tbl, file.path(DATA_OUT_DIR, "precondition_power_tension_earlyVsLater_trialcounts.csv"), row.names = FALSE)

n_label_tbl <- pooled |>
  dplyr::distinct(.data$specimen, .data$trial_id, .data$protocol_family, .data$precondition) |>
  dplyr::count(.data$protocol_family, .data$precondition, name = "n_trials") |>
  dplyr::mutate(x_label = sprintf("%s\n(n=%d trials)", .data$precondition, .data$n_trials))

pooled <- pooled |>
  dplyr::left_join(dplyr::select(n_label_tbl, .data$protocol_family, .data$precondition, .data$x_label),
                    by = c("protocol_family", "precondition"))

# =============================================================================
# Plot: one panel per protocol x metric (mean/max), x = early vs later
# (n-labeled), points colored by specimen, boxplot for the group median/IQR.
# Panels use independent y-axes (scales="free") since power (W/kg) and
# tension (N/cm^2) are different physical quantities -- axis labels make the
# unit explicit per panel, so this is not implying cross-panel comparability.
# =============================================================================
panel_order <- c("Dynamic: mean power (W/kg)", "Isovelocity: mean power (W/kg)", "Isometric: mean tension (N/cm^2)",
                  "Dynamic: max power (W/kg)", "Isovelocity: max power (W/kg)", "Isometric: max tension (N/cm^2)")
pooled$panel_label <- factor(pooled$panel_label, levels = panel_order)
pooled$precondition <- factor(pooled$precondition, levels = DYNAMIC_PRECONDITION_LEVELS)
x_lab_order <- pooled |> dplyr::distinct(.data$precondition, .data$x_label) |> dplyr::arrange(.data$precondition)
pooled$x_label <- factor(pooled$x_label, levels = x_lab_order$x_label)

p <- ggplot(pooled, aes(x = .data$x_label, y = .data$value)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, fill = "grey95", color = "grey40") +
  geom_jitter(aes(color = .data$specimen), width = 0.12, size = 2.4, alpha = 0.85) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  facet_wrap(~panel_label, scales = "free", ncol = 3) +
  scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen") +
  labs(title = "Mass-specific muscle power and specific tension: early (preconditioning) vs. later (stable) trials",
       subtitle = "Each point = one trial (dynamic_trial_precondition.R specimen-specific cutoff). Top row = trial-mean; bottom row = trial-max/peak.\nDynamic and isovelocity report power (W/kg); isometric reports specific tension (N/cm^2, no external power by design). Independent y-axes per panel.",
       x = NULL, y = "Value (unit per panel)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", strip.text = element_text(size = 9),
        axis.text.x = element_text(size = 8))

fout <- file.path(OUT_DIR, "precondition_power_tension_earlyVsLater.png")
ggplot2::ggsave(fout, p, width = 12, height = 7.5, dpi = 150)
cli::cli_alert_success("Saved {fout}")

cli::cli_alert_success("summary_precondition_power_tension_earlyVsLater.R complete -- outputs in {OUT_DIR}/ and {DATA_OUT_DIR}/")
