# diag_precondition_power_vs_offset_isovelocity.R
# Extends diag_precondition_power_vs_offset.R's "does high measured power
# coincide with high sono-strain offset (slip)?" check to ISOVELOCITY trials
# (PI-directed, 2026-07-24: "perform a peak/avg power analysis on the other
# protocol types as well").
#
# Why ISOVELOCITY only (not isometric, not frequency_sweep) -- see chat
# summary for the full explanation given to the PI:
#   - isometric: the motor does not move during the activation window (that
#     is the definition of "isometric"), so mechanical power there is ~zero
#     by design, not by measurement -- there is no meaningful power signal
#     to correlate against offset. FORCE is the natural isometric analog;
#     not computed here (separate ask).
#   - frequency_sweep: passive-only (active_mask always FALSE in
#     plot_sono_strain_validation_pooled.R / run_fv_fl_power_pipeline.R --
#     no stimulation is ever delivered), so calc_muscle_torque()'s
#     cycletype == "act" filter (the existing per-cycle power framework)
#     returns zero rows for every frequency_sweep trial. There is no
#     muscle-driven power to analyze.
#   - isovelocity: genuinely produces mechanical power (force x commanded
#     velocity) and the codebase already has a VETTED conversion from the
#     per-step torque/strain-rate columns to real Watts -- see
#     muscle_geometry.R's add_specific_properties_to_fit() header:
#     "real muscle power at that operating point is peak_power_composite *
#     dclamp_m / (100 * r_m)" where peak_power_composite = F0(muscle_force_Nm,
#     itself a TORQUE despite the "_Nm" force-like name) * Vmax(shortening_
#     strain_pct, %/s for isovelocity per 03_analyze.R's build_segmented_
#     step_summary() comment). That formula is torque x angular_velocity =
#     power; this script applies the IDENTICAL, already-vetted conversion
#     PER STEP (using each step's own measured muscle_force_Nm and
#     shortening_strain_pct) instead of only at the Hill fit's single
#     extrapolated Vmax -- no new physics, just a finer-grained application
#     of the same formula.
#
# Restricted to muscle_side == "right" steps (the sono-instrumented muscle),
# matching diag_precondition_power_vs_offset.R's dynamic analysis convention
# of comparing power on the SAME muscle whose strain is being validated.
#
# Data sources:
#   (a) sono-strain offset per trial -- reused from the already-vetted pooled
#       CSV (sono_strain_validation_pooled_samples.csv), filtered to
#       protocol_family == "isovelocity" & active (right stim), exactly the
#       same recipe diag_precondition_power_check.R used for dynamic.
#   (b) per-step muscle power -- computed fresh via analyze_isovelocity()
#       (03_analyze.R), which this script does NOT duplicate/reimplement.
#
# Run with:  Rscript R/diag_precondition_power_vs_offset_isovelocity.R
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
                              # .legacy_peak_window_mean(); without it every step's active_force_Nm
                              # (and hence muscle_force_Nm/power) errors out (caught 2026-07-24).
src("fit_fv_fl.R")           # fit_force_velocity_curve() -- REQUIRED by analyze_isovelocity() (caught 2026-07-24).
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
# Part A: per-step isovelocity power, right muscle only (freshly computed)
# =============================================================================
cli::cli_h1("Part A: per-step isovelocity muscle power, right muscle, bass16/17/18")

.load_one <- function(f) {
  td <- load_bender_flat(f, do_filter = TRUE, loadtorques = "x")
  tau <- deconvolve_bender(f, hub_path = NULL, verbose = FALSE)
  N <- min(nrow(td), length(tau)); td <- td[seq_len(N), , drop = FALSE]
  td$torque_inertia_corrected_Nm <- tau[seq_len(N)]
  attr(td, "Filename") <- f
  td
}

.collect_isovelocity_steps <- function(specimen, source_dir) {
  manifest <- parse_trial_directory(source_dir)
  files <- manifest$fullpath[manifest$protocol == "isovelocity"]
  purrr::map_dfr(files, function(f) {
    trial_id <- tools::file_path_sans_ext(basename(f))
    td <- tryCatch(.load_one(f), error = function(e) NULL)
    if (is.null(td)) return(tibble())
    res <- tryCatch(analyze_isovelocity(td, f), error = function(e) {
      cli::cli_alert_danger("{trial_id}: analyze_isovelocity failed: {conditionMessage(e)}"); NULL
    })
    if (is.null(res) || nrow(res$step_summary) == 0L) return(tibble())

    # Same conversion muscle_geometry.R::add_specific_properties_to_fit()
    # already applies to the Hill fit's composite peak_power (Nm*(%/s) ->
    # W): torque(Nm) x strain_rate(%/s) x dclamp_m/(100*r_m) = torque x
    # angular_velocity(rad/s) = physical power (W). Applied per-STEP here.
    power_conv <- res$dclamp_m / (100.0 * res$r_m)
    steps <- res$step_summary
    steps$power_W   <- steps$muscle_force_Nm * steps$shortening_strain_pct * power_conv
    steps$power_Wkg <- steps$power_W / res$muscle$muscle_mass_kg

    steps$specimen  <- specimen
    steps$trial_id  <- trial_id
    steps$trial_num <- extract_bender_trial_num(trial_id)
    cli::cli_alert_success("{trial_id}: {sum(steps$muscle_side == 'right', na.rm = TRUE)} right-side steps")
    steps
  })
}

isovel_steps <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  .collect_isovelocity_steps(specimen, raw_source_dir(subfolder))
})

write.csv(isovel_steps, file.path(DATA_OUT_DIR, "isovelocity_power_percycle.csv"), row.names = FALSE)
cli::cli_alert_success("Saved {nrow(isovel_steps)} isovelocity steps (all sides) -> data_processed/isovelocity_power_percycle.csv")

trial_power <- isovel_steps |>
  dplyr::filter(.data$muscle_side == "right", is.finite(.data$power_Wkg)) |>
  dplyr::group_by(.data$specimen, .data$trial_id, .data$trial_num) |>
  dplyr::summarise(
    n_steps        = dplyr::n(),
    mean_power_Wkg = mean(.data$power_Wkg, na.rm = TRUE),
    max_power_Wkg  = max(abs(.data$power_Wkg), na.rm = TRUE),
    .groups = "drop"
  )
print(trial_power, n = 50)

# =============================================================================
# Part B: sono-strain offset per trial (reused, not recomputed)
# =============================================================================
cli::cli_h1("Part B: sono-strain offset by trial, isovelocity, active (right stim)")

pooled_path <- file.path(DATA_OUT_DIR, "sono_strain_validation_pooled_samples.csv")
pooled <- readr::read_csv(pooled_path, show_col_types = FALSE)

offset_by_trial <- pooled |>
  dplyr::filter(.data$protocol_family == "isovelocity", .data$active_passive == "active (right stim)",
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
  dplyr::inner_join(trial_power, by = c("specimen", "trial_id"))

write.csv(joined, file.path(DATA_OUT_DIR, "isovelocity_power_vs_offset_by_trial.csv"), row.names = FALSE)
cli::cli_alert_success("Joined {nrow(joined)} isovelocity trials (power x offset) -> data_processed/isovelocity_power_vs_offset_by_trial.csv")

if (nrow(joined) < 4L) {
  cli::cli_alert_warning("Fewer than 4 joined isovelocity trials -- skipping plots (insufficient sono coverage or step data)")
} else {
  .cor_label <- function(x, y) {
    ct <- suppressWarnings(cor.test(x, y))
    sprintf("pooled: r=%.3f, p=%.3g, n=%d", unname(ct$estimate), ct$p.value, sum(is.finite(x) & is.finite(y)))
  }

  lab1 <- .cor_label(joined$offset_pct, joined$mean_power_Wkg)
  cli::cli_alert_info("Isovelocity mean power vs offset -- {lab1}")
  p1 <- ggplot(joined, aes(x = offset_pct, y = mean_power_Wkg)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7, fill = "grey80") +
    geom_point(aes(color = specimen), size = 2.8) +
    scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen") +
    annotate("text", x = -Inf, y = Inf, label = lab1, hjust = -0.05, vjust = 1.4, size = 3.4) +
    labs(title = "Isovelocity: trial-mean normalized muscle power vs. trial-mean sono-strain offset",
         subtitle = "Right-muscle steps only. Power = muscle_force_Nm x shortening_strain_pct x dclamp_m/(100*r_m) (torque x angular velocity), per step.\nOffset = mean(sono strain) - mean(encoder-predicted strain), active-window, %L0. No hard early/later cutoff exists for isovelocity (dynamic-only, see dynamic_trial_precondition.R).",
         x = "Trial-mean sono-strain offset (percentage points of L0)",
         y = "Trial-mean per-step power (W/kg)") +
    theme_bw(base_size = 11) + theme(legend.position = "right")
  fout1 <- file.path(OUT_DIR, "isovelocity_meanpower_vs_offset.png")
  ggplot2::ggsave(fout1, p1, width = 8, height = 5.5, dpi = 150)
  cli::cli_alert_success("Saved {fout1}")

  lab2 <- .cor_label(joined$offset_pct, joined$max_power_Wkg)
  cli::cli_alert_info("Isovelocity max power vs offset -- {lab2}")
  p2 <- ggplot(joined, aes(x = offset_pct, y = max_power_Wkg)) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7, fill = "grey80") +
    geom_point(aes(color = specimen), size = 2.8) +
    scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen") +
    annotate("text", x = -Inf, y = Inf, label = lab2, hjust = -0.05, vjust = 1.4, size = 3.4) +
    labs(title = "Isovelocity: trial-max normalized peak muscle power vs. trial-mean sono-strain offset",
         subtitle = "Right-muscle steps only (max |power_Wkg| across steps in the trial). Same power formula as above.",
         x = "Trial-mean sono-strain offset (percentage points of L0)",
         y = "Trial-max |per-step power| (W/kg)") +
    theme_bw(base_size = 11) + theme(legend.position = "right")
  fout2 <- file.path(OUT_DIR, "isovelocity_maxpower_vs_offset.png")
  ggplot2::ggsave(fout2, p2, width = 8, height = 5.5, dpi = 150)
  cli::cli_alert_success("Saved {fout2}")
}

cli::cli_alert_success("diag_precondition_power_vs_offset_isovelocity.R complete -- outputs in {OUT_DIR}/ and {DATA_OUT_DIR}/")
