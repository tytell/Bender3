# diag_precondition_passive_sono_fidelity.R
# PI question, 2026-07-24: "Is the sono fidelity in passive trials ALSO
# better than earlier ones?" -- i.e. does the early (preconditioning) ->
# later (stable) improvement seen in the ACTIVE (right-stim) dynamic samples
# (dynamic_trial_precondition.R cutoff; r 0.293 -> 0.905, RMSE 6.34 -> 0.96,
# offset 4.92 -> 0.43 pct-points) ALSO show up in PASSIVE (no/left stim)
# dynamic samples, where there is no muscle tension actively pulling on the
# grips?
#
# Purely a re-slice of the already-vetted pooled sono CSV
# (sono_strain_validation_pooled_samples.csv, from plot_sono_strain_
# validation_pooled.R) by precondition status -- no new data collection or
# processing, just the PASSIVE half of the same table the ACTIVE analysis
# already used (diag_precondition_power_check.R Part A).
#
# Run with:  Rscript R/diag_precondition_passive_sono_fidelity.R
# Outputs -> figures: 02_processed/figs_diagnostic/ (FIGS_DIAGNOSTIC_DIR)
#            data:    02_processed/data_processed/

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(ggplot2); library(cli); library(readr)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")
src("dynamic_trial_precondition.R")

OUT_DIR      <- FIGS_DIAGNOSTIC_DIR
DATA_OUT_DIR <- file.path(.crittergripper_root(), "02_processed", "data_processed")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SPECIMEN_COLORS <- c(bass16 = "#1b9e77", bass17 = "#d95f02", bass18 = "#7570b3")

pooled_path <- file.path(DATA_OUT_DIR, "sono_strain_validation_pooled_samples.csv")
pooled <- readr::read_csv(pooled_path, show_col_types = FALSE)

d <- pooled |>
  dplyr::filter(.data$protocol_family == "dynamic") |>
  dplyr::mutate(
    trial_num    = extract_bender_trial_num(.data$trial_id),
    precondition = classify_dynamic_precondition(.data$specimen, .data$trial_num)
  ) |>
  dplyr::filter(!is.na(.data$precondition))

.fidelity_table <- function(df, group_cols) {
  df |>
    dplyr::filter(is.finite(.data$strain_pred_encoder_right_pct), is.finite(.data$strain_sono_pct)) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      n          = dplyr::n(),
      r          = suppressWarnings(cor(.data$strain_pred_encoder_right_pct, .data$strain_sono_pct)),
      rmse       = sqrt(mean((.data$strain_pred_encoder_right_pct - .data$strain_sono_pct)^2)),
      offset_pct = mean(.data$strain_sono_pct) - mean(.data$strain_pred_encoder_right_pct),
      .groups = "drop"
    )
}

cli::cli_h1("Dynamic PASSIVE (no/left stim) samples, by precondition")
passive_tbl <- .fidelity_table(dplyr::filter(d, .data$active_passive == "passive (no/left stim)"), "precondition")
print(passive_tbl)
passive_by_specimen_tbl <- .fidelity_table(dplyr::filter(d, .data$active_passive == "passive (no/left stim)"), c("specimen", "precondition"))
print(passive_by_specimen_tbl, n = 20)

cli::cli_h1("Dynamic ACTIVE (right stim) samples, by precondition (for comparison)")
active_tbl <- .fidelity_table(dplyr::filter(d, .data$active_passive == "active (right stim)"), "precondition")
print(active_tbl)

combined_tbl <- dplyr::bind_rows(
  dplyr::mutate(passive_tbl, active_passive = "passive (no/left stim)"),
  dplyr::mutate(active_tbl, active_passive = "active (right stim)")
)
write.csv(combined_tbl, file.path(DATA_OUT_DIR, "dynamic_precondition_sono_fidelity_active_vs_passive.csv"), row.names = FALSE)
write.csv(passive_by_specimen_tbl, file.path(DATA_OUT_DIR, "dynamic_precondition_sono_fidelity_passive_by_specimen.csv"), row.names = FALSE)
cli::cli_alert_success("Saved fidelity tables -> data_processed/dynamic_precondition_sono_fidelity_*.csv")

# =============================================================================
# Plot: passive dynamic sono validation, early vs. later, side by side
# (mirrors the active-only version already used throughout this
# investigation, e.g. plot_sono_strain_validation_pooled.R's family plot).
# =============================================================================
df_plot <- dplyr::filter(d, .data$active_passive == "passive (no/left stim)",
                          is.finite(.data$strain_pred_encoder_right_pct), is.finite(.data$strain_sono_pct))
lims <- range(c(df_plot$strain_pred_encoder_right_pct, df_plot$strain_sono_pct), na.rm = TRUE)
ref_df <- tibble(x = lims, y = lims)

label_tbl <- passive_tbl |>
  dplyr::mutate(label = sprintf("r=%.3f, RMSE=%.2f, offset=%.3f pct-pts, n=%s",
                                 .data$r, .data$rmse, .data$offset_pct, format(.data$n, big.mark = ",")))

p <- ggplot(df_plot, aes(x = .data$strain_pred_encoder_right_pct, y = .data$strain_sono_pct, color = .data$specimen)) +
  geom_point(shape = 1, size = 0.6, alpha = 0.12, stroke = 0.25) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_text(data = label_tbl, aes(x = -Inf, y = Inf, label = .data$label), inherit.aes = FALSE,
            hjust = -0.05, vjust = 1.4, size = 3.2, color = "black") +
  facet_wrap(~precondition) +
  scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen", drop = FALSE) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2))) +
  coord_equal(xlim = lims, ylim = lims) +
  labs(title = "Dynamic, PASSIVE (no/left stim) samples: sono strain vs. encoder-predicted strain, early vs. later trials",
       subtitle = "Right muscle sono channel; dashed = 1:1 reference. Same specimen-specific early/later cutoff used throughout this investigation\n(dynamic_trial_precondition.R), applied here to PASSIVE samples for comparison against the ACTIVE (right-stim) version already built.",
       x = "Predicted strain (%, from encoder angle, right-folded)",
       y = "Measured strain (%, sonomicrometry, right muscle)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
fout <- file.path(OUT_DIR, "dynamic_precondition_passive_sonoValidation_earlyVsLater.png")
ggplot2::ggsave(fout, p, width = 10, height = 5.5, dpi = 150)
cli::cli_alert_success("Saved {fout}")

cli::cli_alert_success("diag_precondition_passive_sono_fidelity.R complete -- outputs in {OUT_DIR}/ and {DATA_OUT_DIR}/")
