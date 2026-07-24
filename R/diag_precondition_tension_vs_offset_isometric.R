# diag_precondition_tension_vs_offset_isometric.R
# Isometric counterpart of diag_precondition_power_vs_offset.R (dynamic) /
# diag_precondition_power_vs_offset_isovelocity.R (isovelocity), using
# SPECIFIC TENSION (N/cm^2) instead of power (PI-directed, 2026-07-24:
# "Please run analysis for isometric using tension" -- isometric contractions
# do no external mechanical work by design, since the motor does not move
# during the activation window, so POWER is ~zero there regardless of
# muscle effort; TENSION/FORCE is the physically correct isometric analog
# of "how hard is the muscle working", see analysis_muscle_force_vector_log
# .md's 2026-07-24 entry for the full reasoning given to the PI).
#
# Tension conversion is the SAME, already-vetted formula
# muscle_geometry.R::add_specific_properties_to_fit() applies to a FV/FL
# fit's single F0 anchor point (Force0_N = F0_Nm / r_m; specific_tension_
# Ncm2 = Force0_N / csa_muscle_cm2) -- applied here PER STEP (each step's own
# measured muscle_force_Nm, a TORQUE despite the "_force" name) instead of
# only at one fit-level F0. No new physics, a finer-grained application of
# an already-vetted conversion.
#
# Restricted to muscle_side == "right" steps (the sono-instrumented muscle),
# matching the dynamic/isovelocity scripts' side convention.
#
# Data sources:
#   (a) sono-strain offset per trial -- reused from the already-vetted pooled
#       CSV (sono_strain_validation_pooled_samples.csv), filtered to
#       protocol_family == "isometric" & active (right stim).
#   (b) per-step muscle tension -- computed fresh via analyze_isometric()
#       (03_analyze.R), which this script does NOT duplicate/reimplement.
#
# Run with:  Rscript R/diag_precondition_tension_vs_offset_isometric.R
# Outputs -> figures: 02_processed/figs_diagnostic/ (FIGS_DIAGNOSTIC_DIR)
#            data:    02_processed/data_processed/

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(ggplot2); library(cli); library(readr)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")
src("00_load_bender_flat.R")
src("01_calibrate.R")
src("02_deconvolve.R")
src("muscle_geometry.R")
src("plot_force_vs_time.R")  # .smooth_trace_display_only() -- REQUIRED by build_segmented_step_summary()'s
                              # .legacy_peak_window_mean() (same dependency gap found + fixed in the
                              # isovelocity script, 2026-07-24).
src("03_analyze.R")
src("parse_trial_filename.R")
src("dynamic_trial_precondition.R")  # extract_bender_trial_num()

OUT_DIR      <- FIGS_DIAGNOSTIC_DIR
DATA_OUT_DIR <- file.path(.crittergripper_root(), "02_processed", "data_processed")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SPECIMEN_SUBFOLDERS <- c(bass16 = BASS16_RAW_SUBFOLDER, bass17 = BASS17_RAW_SUBFOLDER, bass18 = BASS18_RAW_SUBFOLDER)
SPECIMEN_COLORS     <- c(bass16 = "#1b9e77", bass17 = "#d95f02", bass18 = "#7570b3")

# =============================================================================
# Part A: per-step isometric specific tension, right muscle only
# =============================================================================
cli::cli_h1("Part A: per-step isometric specific tension, right muscle, bass16/17/18")

.load_one <- function(f) {
  td <- load_bender_flat(f, do_filter = TRUE, loadtorques = "x")
  tau <- deconvolve_bender(f, hub_path = NULL, verbose = FALSE)
  N <- min(nrow(td), length(tau)); td <- td[seq_len(N), , drop = FALSE]
  td$torque_inertia_corrected_Nm <- tau[seq_len(N)]
  attr(td, "Filename") <- f
  td
}

.collect_isometric_steps <- function(specimen, source_dir) {
  manifest <- parse_trial_directory(source_dir)
  files <- manifest$fullpath[manifest$protocol == "isometric"]
  purrr::map_dfr(files, function(f) {
    trial_id <- tools::file_path_sans_ext(basename(f))
    td <- tryCatch(.load_one(f), error = function(e) NULL)
    if (is.null(td)) return(tibble())
    res <- tryCatch(analyze_isometric(td, f), error = function(e) {
      cli::cli_alert_danger("{trial_id}: analyze_isometric failed: {conditionMessage(e)}"); NULL
    })
    if (is.null(res) || nrow(res$step_summary) == 0L) return(tibble())
    if (!is.finite(res$r_m) || res$r_m <= 0 ||
        !is.finite(res$muscle$csa_muscle_cm2) || res$muscle$csa_muscle_cm2 <= 0) {
      cli::cli_alert_warning("{trial_id}: missing/invalid r_m or muscle CSA -- skipped")
      return(tibble())
    }

    # Same conversion muscle_geometry.R::add_specific_properties_to_fit()
    # already applies to a fit's single F0: Force0_N = F0_Nm / r_m (lever-arm
    # torque -> force), specific_tension_Ncm2 = Force0_N / csa_muscle_cm2.
    # Applied per-STEP here (each step's own muscle_force_Nm), not just at
    # one fit-level F0.
    steps <- res$step_summary
    steps$force_N               <- steps$muscle_force_Nm / res$r_m
    steps$specific_tension_Ncm2 <- steps$force_N / res$muscle$csa_muscle_cm2

    steps$specimen  <- specimen
    steps$trial_id  <- trial_id
    steps$trial_num <- extract_bender_trial_num(trial_id)
    cli::cli_alert_success("{trial_id}: {sum(steps$muscle_side == 'right', na.rm = TRUE)} right-side steps")
    steps
  })
}

iso_steps <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  .collect_isometric_steps(specimen, raw_source_dir(subfolder))
})

write.csv(iso_steps, file.path(DATA_OUT_DIR, "isometric_tension_percycle.csv"), row.names = FALSE)
cli::cli_alert_success("Saved {nrow(iso_steps)} isometric steps (all sides) -> data_processed/isometric_tension_percycle.csv")

trial_tension <- iso_steps |>
  dplyr::filter(.data$muscle_side == "right", is.finite(.data$specific_tension_Ncm2)) |>
  dplyr::group_by(.data$specimen, .data$trial_id, .data$trial_num) |>
  dplyr::summarise(
    n_steps           = dplyr::n(),
    mean_tension_Ncm2 = mean(.data$specific_tension_Ncm2, na.rm = TRUE),
    max_tension_Ncm2  = max(abs(.data$specific_tension_Ncm2), na.rm = TRUE),
    .groups = "drop"
  )
print(trial_tension, n = 50)

# =============================================================================
# Part B: sono-strain offset per trial (reused, not recomputed)
# =============================================================================
cli::cli_h1("Part B: sono-strain offset by trial, isometric, active (right stim)")

pooled_path <- file.path(DATA_OUT_DIR, "sono_strain_validation_pooled_samples.csv")
pooled <- readr::read_csv(pooled_path, show_col_types = FALSE)

offset_by_trial <- pooled |>
  dplyr::filter(.data$protocol_family == "isometric", .data$active_passive == "active (right stim)",
                is.finite(.data$strain_pred_encoder_right_pct), is.finite(.data$strain_sono_pct)) |>
  dplyr::group_by(.data$specimen, .data$trial_id) |>
  dplyr::summarise(n = dplyr::n(),
                    offset_pct = mean(.data$strain_sono_pct) - mean(.data$strain_pred_encoder_right_pct),
                    r = suppressWarnings(cor(.data$strain_pred_encoder_right_pct, .data$strain_sono_pct)),
                    .groups = "drop")
print(offset_by_trial, n = 50)

# =============================================================================
# Join + plot
# =============================================================================
joined <- offset_by_trial |>
  dplyr::inner_join(trial_tension, by = c("specimen", "trial_id"))

write.csv(joined, file.path(DATA_OUT_DIR, "isometric_tension_vs_offset_by_trial.csv"), row.names = FALSE)
cli::cli_alert_success("Joined {nrow(joined)} isometric trials (tension x offset) -> data_processed/isometric_tension_vs_offset_by_trial.csv")

if (nrow(joined) < 4L) {
  cli::cli_alert_warning("Only {nrow(joined)} joined isometric trial(s) (only 5 isometric trials exist across all 3 specimens combined, and not all have right-stim sono coverage) -- too few for a correlation/regression line. Plotting the raw points only, NO fitted trend line.")
} else {
  cli::cli_alert_info("n={nrow(joined)} -- proceeding with correlation + trend line")
}

.cor_label <- function(x, y) {
  n_ok <- sum(is.finite(x) & is.finite(y))
  if (n_ok < 4L) return(sprintf("n=%d (too few for a correlation)", n_ok))
  ct <- suppressWarnings(cor.test(x, y))
  sprintf("pooled: r=%.3f, p=%.3g, n=%d", unname(ct$estimate), ct$p.value, n_ok)
}

lab1 <- .cor_label(joined$offset_pct, joined$mean_tension_Ncm2)
cli::cli_alert_info("Isometric mean tension vs offset -- {lab1}")
p1 <- ggplot(joined, aes(x = offset_pct, y = mean_tension_Ncm2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60")
if (nrow(joined) >= 4L) p1 <- p1 + geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7, fill = "grey80")
p1 <- p1 +
  geom_point(aes(color = specimen), size = 3.2) +
  scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen") +
  annotate("text", x = -Inf, y = Inf, label = lab1, hjust = -0.05, vjust = 1.4, size = 3.4) +
  labs(title = "Isometric: trial-mean specific tension vs. trial-mean sono-strain offset",
       subtitle = "Right-muscle steps only. Tension = (muscle_force_Nm / r_m) / muscle_csa_cm2 (lever-arm torque -> force -> area-specific tension), per step.\nOffset = mean(sono strain) - mean(encoder-predicted strain), active-window, %L0. Only 5 isometric trials exist in the whole corpus -- interpret cautiously.",
       x = "Trial-mean sono-strain offset (percentage points of L0)",
       y = "Trial-mean specific tension (N/cm^2)") +
  theme_bw(base_size = 11) + theme(legend.position = "right")
fout1 <- file.path(OUT_DIR, "isometric_meantension_vs_offset.png")
ggplot2::ggsave(fout1, p1, width = 8, height = 5.5, dpi = 150)
cli::cli_alert_success("Saved {fout1}")

lab2 <- .cor_label(joined$offset_pct, joined$max_tension_Ncm2)
cli::cli_alert_info("Isometric max tension vs offset -- {lab2}")
p2 <- ggplot(joined, aes(x = offset_pct, y = max_tension_Ncm2))
if (nrow(joined) >= 4L) p2 <- p2 + geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7, fill = "grey80")
p2 <- p2 +
  geom_point(aes(color = specimen), size = 3.2) +
  scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen") +
  annotate("text", x = -Inf, y = Inf, label = lab2, hjust = -0.05, vjust = 1.4, size = 3.4) +
  labs(title = "Isometric: trial-max specific tension vs. trial-mean sono-strain offset",
       subtitle = "Right-muscle steps only (max |specific_tension_Ncm2| across steps in the trial). Same formula as above.",
       x = "Trial-mean sono-strain offset (percentage points of L0)",
       y = "Trial-max |specific tension| (N/cm^2)") +
  theme_bw(base_size = 11) + theme(legend.position = "right")
fout2 <- file.path(OUT_DIR, "isometric_maxtension_vs_offset.png")
ggplot2::ggsave(fout2, p2, width = 8, height = 5.5, dpi = 150)
cli::cli_alert_success("Saved {fout2}")

cli::cli_alert_success("diag_precondition_tension_vs_offset_isometric.R complete -- outputs in {OUT_DIR}/ and {DATA_OUT_DIR}/")
