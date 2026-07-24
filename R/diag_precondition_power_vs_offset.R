# diag_precondition_power_vs_offset.R
# Direct follow-up to diag_precondition_power_check.R's trial-order plots
# (PI, 2026-07-24): "dynamic_precondition_offset_vs_trialorder.png and
# dynamic_precondition_power_vs_trialorder.png are somewhat disconcerting.
# The implication is that high power measurements (a good thing, presumably)
# are collected only when the muscle is slipping and measured strain offset
# is highest." Rather than inferring this indirectly from two separate
# trial-order plots, this script correlates normalized muscle power DIRECTLY
# against the sono-strain offset, per trial.
#
# Reuses the two CSVs diag_precondition_power_check.R already wrote (no raw
# H5 reprocessing needed):
#   dynamic_precondition_power_percycle.csv  -- per-cycle avg_power.Wkg /
#     peak_power.Wkg (mass-normalized, from summarize_muscle_cycles())
#   dynamic_precondition_offset_by_trial.csv -- per-trial sono-strain offset
#     (mean(sono) - mean(encoder-predicted), % of L0) and pointwise r
#
# Aggregated to TRIAL level (not per-cycle) because the offset metric itself
# is only meaningful as a trial-level summary (see dynamic_trial_
# precondition.R docstring) -- a per-cycle join would require reconciling
# two different cycle-index conventions (calc_muscle_torque's `cycle` =
# floor(t.norm) vs the HDF5 cycle_index dataset), which is unnecessary
# complexity for this question.
#
# Run with:  Rscript R/diag_precondition_power_vs_offset.R
# Outputs -> figures: 02_processed/figs_diagnostic/ (FIGS_DIAGNOSTIC_DIR)
#            data:    02_processed/data_processed/

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2); library(cli); library(readr)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")

OUT_DIR      <- FIGS_DIAGNOSTIC_DIR
DATA_OUT_DIR <- file.path(.crittergripper_root(), "02_processed", "data_processed")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SPECIMEN_COLORS <- c(bass16 = "#1b9e77", bass17 = "#d95f02", bass18 = "#7570b3")

# =============================================================================
# Load + aggregate to trial level
# =============================================================================
percycle <- readr::read_csv(file.path(DATA_OUT_DIR, "dynamic_precondition_power_percycle.csv"), show_col_types = FALSE)
offset   <- readr::read_csv(file.path(DATA_OUT_DIR, "dynamic_precondition_offset_by_trial.csv"), show_col_types = FALSE)

trial_power <- percycle |>
  dplyr::filter(is.finite(.data$avg_power.Wkg), is.finite(.data$peak_power.Wkg)) |>
  dplyr::group_by(.data$specimen, .data$trial_id) |>
  dplyr::summarise(
    n_cycles           = dplyr::n(),
    mean_avg_power_Wkg = mean(.data$avg_power.Wkg, na.rm = TRUE),
    max_peak_power_Wkg = max(.data$peak_power.Wkg, na.rm = TRUE),
    .groups = "drop"
  )

joined <- offset |>
  dplyr::select(.data$specimen, .data$trial_id, .data$trial_num, .data$precondition, .data$offset_pct, .data$r) |>
  dplyr::inner_join(trial_power, by = c("specimen", "trial_id"))

write.csv(joined, file.path(DATA_OUT_DIR, "dynamic_precondition_power_vs_offset_by_trial.csv"), row.names = FALSE)
cli::cli_alert_success("Joined {nrow(joined)} trials (power x offset) -> data_processed/dynamic_precondition_power_vs_offset_by_trial.csv")
print(joined, n = 100)

.cor_label <- function(x, y) {
  ct <- suppressWarnings(cor.test(x, y))
  sprintf("pooled: r=%.3f, p=%.3g, n=%d", unname(ct$estimate), ct$p.value, sum(is.finite(x) & is.finite(y)))
}

# =============================================================================
# Plot 1: trial-mean normalized power vs. trial sono-strain offset (pooled)
# =============================================================================
lab1 <- .cor_label(joined$offset_pct, joined$mean_avg_power_Wkg)
cli::cli_alert_info("Mean power vs offset -- {lab1}")

p1 <- ggplot(joined, aes(x = offset_pct, y = mean_avg_power_Wkg)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7, fill = "grey80") +
  geom_point(aes(color = specimen, shape = precondition), size = 2.8, stroke = 0.9) +
  scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen") +
  scale_shape_manual(values = c("early (preconditioning)" = 17, "later (stable)" = 16), name = NULL) +
  annotate("text", x = -Inf, y = Inf, label = lab1, hjust = -0.05, vjust = 1.4, size = 3.4) +
  labs(title = "Trial-mean normalized muscle power vs. trial-mean sono-strain offset",
       subtitle = "Each point = one dynamic trial. Offset = mean(sono strain) - mean(encoder-predicted strain), active-window, %L0.\nTriangle = early (preconditioning) trial, circle = later (stable) trial.",
       x = "Trial-mean sono-strain offset (percentage points of L0)",
       y = "Trial-mean cycle-averaged power (W/kg)") +
  theme_bw(base_size = 11) + theme(legend.position = "right")
fout1 <- file.path(OUT_DIR, "dynamic_precondition_meanpower_vs_offset.png")
ggplot2::ggsave(fout1, p1, width = 8, height = 5.5, dpi = 150)
cli::cli_alert_success("Saved {fout1}")

# =============================================================================
# Plot 2: trial-max normalized peak power vs. trial sono-strain offset (pooled)
# =============================================================================
lab2 <- .cor_label(joined$offset_pct, joined$max_peak_power_Wkg)
cli::cli_alert_info("Max peak power vs offset -- {lab2}")

p2 <- ggplot(joined, aes(x = offset_pct, y = max_peak_power_Wkg)) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7, fill = "grey80") +
  geom_point(aes(color = specimen, shape = precondition), size = 2.8, stroke = 0.9) +
  scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen") +
  scale_shape_manual(values = c("early (preconditioning)" = 17, "later (stable)" = 16), name = NULL) +
  annotate("text", x = -Inf, y = Inf, label = lab2, hjust = -0.05, vjust = 1.4, size = 3.4) +
  labs(title = "Trial-max normalized PEAK muscle power vs. trial-mean sono-strain offset",
       subtitle = "Each point = one dynamic trial (max cycle peak_power.Wkg reached in that trial).\nTriangle = early (preconditioning) trial, circle = later (stable) trial.",
       x = "Trial-mean sono-strain offset (percentage points of L0)",
       y = "Trial-max peak instantaneous power (W/kg)") +
  theme_bw(base_size = 11) + theme(legend.position = "right")
fout2 <- file.path(OUT_DIR, "dynamic_precondition_maxpower_vs_offset.png")
ggplot2::ggsave(fout2, p2, width = 8, height = 5.5, dpi = 150)
cli::cli_alert_success("Saved {fout2}")

# =============================================================================
# Plot 3: same as Plot 1, faceted by specimen -- checks whether the pooled
# relationship is a real within-specimen effect or purely a between-specimen
# confound (e.g. one specimen simply produces more power AND happens to have
# the largest early-trial offsets).
# =============================================================================
per_specimen_cor <- joined |>
  dplyr::group_by(.data$specimen) |>
  dplyr::summarise(
    n = dplyr::n(),
    r_mean = suppressWarnings(cor(.data$offset_pct, .data$mean_avg_power_Wkg, use = "complete.obs")),
    r_max  = suppressWarnings(cor(.data$offset_pct, .data$max_peak_power_Wkg, use = "complete.obs")),
    .groups = "drop"
  )
cli::cli_h2("Within-specimen correlations (offset vs. mean/max normalized power)")
print(per_specimen_cor)
write.csv(per_specimen_cor, file.path(DATA_OUT_DIR, "dynamic_precondition_power_vs_offset_within_specimen_cor.csv"), row.names = FALSE)

lab3 <- per_specimen_cor |>
  dplyr::mutate(label = sprintf("r=%.3f, n=%d", .data$r_mean, .data$n))

p3 <- ggplot(joined, aes(x = offset_pct, y = mean_avg_power_Wkg, color = specimen)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  geom_point(aes(shape = precondition), size = 2.6, stroke = 0.9) +
  geom_text(data = lab3, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = -0.05, vjust = 1.4, size = 3.2, color = "black") +
  facet_wrap(~specimen, scales = "free") +
  scale_color_manual(values = SPECIMEN_COLORS, guide = "none") +
  scale_shape_manual(values = c("early (preconditioning)" = 17, "later (stable)" = 16), name = NULL) +
  labs(title = "Trial-mean normalized power vs. sono-strain offset, WITHIN each specimen",
       subtitle = "Does the pooled relationship hold within a specimen, or is it a between-specimen confound? Triangle = early trial, circle = later trial.",
       x = "Trial-mean sono-strain offset (percentage points of L0)",
       y = "Trial-mean cycle-averaged power (W/kg)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
fout3 <- file.path(OUT_DIR, "dynamic_precondition_power_vs_offset_byspecimen.png")
ggplot2::ggsave(fout3, p3, width = 11, height = 4.5, dpi = 150)
cli::cli_alert_success("Saved {fout3}")

# =============================================================================
# Plot 4: same as Plot 3, but using trial-max normalized PEAK power instead
# of trial-mean cycle-averaged power (PI-directed, 2026-07-24: "perform a
# similar analysis using peak power instead of cycle average").
# =============================================================================
lab4 <- per_specimen_cor |>
  dplyr::mutate(label = sprintf("r=%.3f, n=%d", .data$r_max, .data$n))

p4 <- ggplot(joined, aes(x = offset_pct, y = max_peak_power_Wkg, color = specimen)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  geom_point(aes(shape = precondition), size = 2.6, stroke = 0.9) +
  geom_text(data = lab4, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = -0.05, vjust = 1.4, size = 3.2, color = "black") +
  facet_wrap(~specimen, scales = "free") +
  scale_color_manual(values = SPECIMEN_COLORS, guide = "none") +
  scale_shape_manual(values = c("early (preconditioning)" = 17, "later (stable)" = 16), name = NULL) +
  labs(title = "Trial-max normalized PEAK power vs. sono-strain offset, WITHIN each specimen",
       subtitle = "Peak-power counterpart of the mean-power by-specimen figure above. Triangle = early trial, circle = later trial.",
       x = "Trial-mean sono-strain offset (percentage points of L0)",
       y = "Trial-max peak instantaneous power (W/kg)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
fout4 <- file.path(OUT_DIR, "dynamic_precondition_maxpower_vs_offset_byspecimen.png")
ggplot2::ggsave(fout4, p4, width = 11, height = 4.5, dpi = 150)
cli::cli_alert_success("Saved {fout4}")

cli::cli_alert_success("diag_precondition_power_vs_offset.R complete -- outputs in {OUT_DIR}/ and {DATA_OUT_DIR}/")
