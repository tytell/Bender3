# diag_precondition_calibration_gain_check.R
# PI question, 2026-07-24: "Based on the offset, can we just add a
# calibration factor to truly get the desired strain? If so, does the
# active data suggest such a factor would be consistent across conditions
# or variable (and therefore kinda useless)?"
#
# Tests this directly: for every dynamic ACTIVE (right-stim) trial, fits
# strain_sono_pct ~ strain_pred_encoder_right_pct (per trial) and reports
# the slope (a multiplicative GAIN factor) and intercept (an additive
# OFFSET factor) separately, then asks whether either is consistent enough
# across trials/specimens to serve as a fixed correction.
#
# Purely a re-slice/re-fit of the already-vetted pooled sono CSV
# (sono_strain_validation_pooled_samples.csv) -- no new data collection.
#
# Run with:  Rscript R/diag_precondition_calibration_gain_check.R
# Outputs -> figures: 02_processed/figs_diagnostic/ (FIGS_DIAGNOSTIC_DIR)
#            data:    02_processed/data_processed/

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(ggplot2); library(cli); library(readr)
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
  dplyr::filter(.data$protocol_family == "dynamic", .data$active_passive == "active (right stim)") |>
  dplyr::mutate(
    trial_num    = extract_bender_trial_num(.data$trial_id),
    precondition = classify_dynamic_precondition(.data$specimen, .data$trial_num)
  ) |>
  dplyr::filter(!is.na(.data$precondition), is.finite(.data$strain_pred_encoder_right_pct), is.finite(.data$strain_sono_pct))

# =============================================================================
# Per-trial linear fit: strain_sono_pct = intercept + slope * strain_pred_pct
# slope = a hypothetical multiplicative GAIN correction (e.g. from mis-
#   estimated muscle depth / body width); intercept = a hypothetical
#   additive OFFSET correction. If either were a real fixed calibration
#   error, it should be roughly CONSTANT across trials (at least within a
#   specimen, since geometry doesn't change within a session).
# =============================================================================
.fit_one_trial <- function(df) {
  if (nrow(df) < 20L || stats::sd(df$strain_pred_encoder_right_pct) < 1e-6) {
    return(tibble::tibble(n = nrow(df), slope = NA_real_, intercept = NA_real_, r = NA_real_))
  }
  fit <- stats::lm(strain_sono_pct ~ strain_pred_encoder_right_pct, data = df)
  tibble::tibble(
    n         = nrow(df),
    slope     = unname(stats::coef(fit)[2]),
    intercept = unname(stats::coef(fit)[1]),
    r         = suppressWarnings(stats::cor(df$strain_pred_encoder_right_pct, df$strain_sono_pct))
  )
}

per_trial <- d |>
  dplyr::group_by(.data$specimen, .data$trial_id, .data$trial_num, .data$precondition) |>
  dplyr::group_modify(~.fit_one_trial(.x)) |>
  dplyr::ungroup() |>
  dplyr::arrange(.data$specimen, .data$trial_num)

write.csv(per_trial, file.path(DATA_OUT_DIR, "dynamic_precondition_calibration_gain_by_trial.csv"), row.names = FALSE)

cli::cli_h1("Per-trial gain (slope) + offset (intercept), dynamic ACTIVE (right-stim)")
print(per_trial, n = 40)

cli::cli_h1("Consistency summary, by precondition (pooled across specimens)")
consistency_by_precondition <- per_trial |>
  dplyr::filter(is.finite(.data$slope)) |>
  dplyr::group_by(.data$precondition) |>
  dplyr::summarise(
    n_trials         = dplyr::n(),
    median_slope     = median(.data$slope), sd_slope = sd(.data$slope),
    slope_range      = sprintf("%.2f to %.2f", min(.data$slope), max(.data$slope)),
    median_intercept = median(.data$intercept), sd_intercept = sd(.data$intercept),
    intercept_range  = sprintf("%.2f to %.2f", min(.data$intercept), max(.data$intercept)),
    .groups = "drop"
  )
print(consistency_by_precondition)

cli::cli_h1("Consistency summary, LATER (stable) trials only, BY SPECIMEN")
consistency_later_by_specimen <- per_trial |>
  dplyr::filter(is.finite(.data$slope), .data$precondition == "later (stable)") |>
  dplyr::group_by(.data$specimen) |>
  dplyr::summarise(
    n_trials         = dplyr::n(),
    median_slope     = median(.data$slope), sd_slope = sd(.data$slope),
    median_intercept = median(.data$intercept), sd_intercept = sd(.data$intercept),
    .groups = "drop"
  )
print(consistency_later_by_specimen)

write.csv(consistency_by_precondition, file.path(DATA_OUT_DIR, "dynamic_precondition_calibration_gain_consistency_by_precondition.csv"), row.names = FALSE)
write.csv(consistency_later_by_specimen, file.path(DATA_OUT_DIR, "dynamic_precondition_calibration_gain_consistency_later_by_specimen.csv"), row.names = FALSE)

# =============================================================================
# Plot: per-trial gain (slope) vs. trial order, faceted by specimen --
# directly visualizes whether a single calibration factor could work.
# =============================================================================
cutoff_df <- tibble::tibble(specimen = names(DYNAMIC_PRECONDITION_CUTOFF_TRIALNUM),
                             cutoff = DYNAMIC_PRECONDITION_CUTOFF_TRIALNUM)

p_slope <- ggplot(per_trial, aes(x = trial_num, y = slope, color = precondition)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  geom_vline(data = cutoff_df, aes(xintercept = cutoff - 0.5), linetype = "dotted", color = "black") +
  geom_point(size = 2.4) +
  facet_wrap(~specimen, scales = "free_x") +
  scale_color_manual(values = c("early (preconditioning)" = "#d95f02", "later (stable)" = "#1b9e77"), drop = FALSE) +
  labs(title = "Per-trial GAIN (regression slope, sono strain vs. predicted strain) by trial order",
       subtitle = "Dashed line = slope of 1 (perfect gain, no correction needed). Dotted line = the specimen-specific preconditioning cutoff.\nA usable single calibration factor would need this to be roughly FLAT and consistent -- it is not, in early trials.",
       x = "Trial number (bender_NN, chronological)", y = "Gain (regression slope, dimensionless)", color = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
fout_slope <- file.path(OUT_DIR, "dynamic_precondition_calibration_gain_vs_trialorder.png")
ggplot2::ggsave(fout_slope, p_slope, width = 10, height = 4.5, dpi = 150)
cli::cli_alert_success("Saved {fout_slope}")

p_intercept <- ggplot(per_trial, aes(x = trial_num, y = intercept, color = precondition)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(data = cutoff_df, aes(xintercept = cutoff - 0.5), linetype = "dotted", color = "black") +
  geom_point(size = 2.4) +
  facet_wrap(~specimen, scales = "free_x") +
  scale_color_manual(values = c("early (preconditioning)" = "#d95f02", "later (stable)" = "#1b9e77"), drop = FALSE) +
  labs(title = "Per-trial OFFSET (regression intercept, sono strain vs. predicted strain) by trial order",
       subtitle = "Dashed line = intercept of 0 (no offset needed). Dotted line = the specimen-specific preconditioning cutoff.",
       x = "Trial number (bender_NN, chronological)", y = "Offset (regression intercept, %L0)", color = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
fout_intercept <- file.path(OUT_DIR, "dynamic_precondition_calibration_offset_vs_trialorder.png")
ggplot2::ggsave(fout_intercept, p_intercept, width = 10, height = 4.5, dpi = 150)
cli::cli_alert_success("Saved {fout_intercept}")

cli::cli_alert_success("diag_precondition_calibration_gain_check.R complete -- outputs in {OUT_DIR}/ and {DATA_OUT_DIR}/")
