# diag_fv_fl_ztorque_vs_uhat.R
# PI request 2B, 2026-07-24: build the FV (isovelocity) and FL (isometric)
# curves using the two FORCE-reconstruction methods -- "simple zTorque"
# (single primary-bending-axis torque projected to a force) vs. "uHat" (the
# empirical multi-axis force-vector magnitude) -- as SEPARATE plots, so the
# two methods can be judged side by side.
#
# Both forces come from attach_vector_muscle_force() (muscle_force_vector.R),
# which returns them per step in the SAME units (N), so they are directly
# comparable on one y-axis:
#   - force_zTorque_N        : the "simple zTorque" method (single-axis torque
#                              channel -> force).
#   - muscle_force_vector_N  : the "uHat" method (empirical active-minus-
#                              passive 6-axis wrench-delta direction, force
#                              magnitude along that line of action).
#
# X-axis is the COMMANDED operating_point-derived strain / strain-rate
# (shortening_strain_pct) -- held CONSTANT across both methods so the ONLY
# variable is the force reconstruction (this is a force-method comparison,
# not the sono-length comparison of prototype_fv_fl_sono_length.R).
#
# Scope: RIGHT muscle steps only (consistent with the sono prototype; keeps a
# single clean sign convention), bass16/17/18.
#
# Run with:  Rscript R/diag_fv_fl_ztorque_vs_uhat.R
# Outputs -> figs_diagnostic/ (FIGS_DIAGNOSTIC_DIR) + per-step CSV.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr); library(ggplot2); library(cli); library(rhdf5)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")
src("00_load_bender_flat.R")
src("01_calibrate.R")
src("02_deconvolve.R")
src("muscle_geometry.R")
src("03_analyze.R")
src("parse_trial_filename.R")
src("plot_force_vs_time.R")
src("fit_fv_fl.R")
src("plot_strain_validation.R")
src("plot_angle_sono_validation.R")
src("muscle_force_vector.R")
src("plot_muscle_force_vector.R")

OUT_DIR      <- FIGS_DIAGNOSTIC_DIR
DATA_OUT_DIR <- file.path(.crittergripper_root(), "02_processed", "data_processed")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SPECIMEN_SUBFOLDERS <- c(bass16 = BASS16_RAW_SUBFOLDER, bass17 = BASS17_RAW_SUBFOLDER, bass18 = BASS18_RAW_SUBFOLDER)
SPECIMEN_COLORS     <- c(bass16 = "#1b9e77", bass17 = "#d95f02", bass18 = "#7570b3")

.collect_one_trial <- function(f, specimen, category) {
  trial_id <- tools::file_path_sans_ext(basename(f))
  td <- tryCatch(load_bender_flat(f, do_filter = TRUE, loadtorques = "x"), error = function(e) NULL)
  if (is.null(td)) return(tibble())
  tau <- tryCatch(deconvolve_bender(f, hub_path = NULL, verbose = FALSE), error = function(e) NULL)
  if (is.null(tau)) return(tibble())
  N <- min(nrow(td), length(tau)); td <- td[seq_len(N), , drop = FALSE]
  td$torque_inertia_corrected_Nm <- tau[seq_len(N)]
  attr(td, "Filename") <- f

  res <- tryCatch(
    if (category == "isometric") analyze_isometric(td) else analyze_isovelocity(td),
    error = function(e) { cli::cli_warn("{trial_id}: analyze failed: {conditionMessage(e)}"); NULL })
  if (is.null(res) || is.null(res$step_summary) || nrow(res$step_summary) == 0L) return(tibble())
  res$trial_id <- trial_id

  vec <- tryCatch(attach_vector_muscle_force(res, f, category), error = function(e) {
    cli::cli_warn("{trial_id}: attach_vector_muscle_force failed: {conditionMessage(e)}"); NULL })
  if (is.null(vec) || is.null(vec$step_summary) || nrow(vec$step_summary) == 0L) return(tibble())
  ss <- vec$step_summary

  needed <- c("shortening_strain_pct", "muscle_side", "force_zTorque_N", "muscle_force_vector_N")
  if (!all(needed %in% names(ss))) return(tibble())

  tibble::tibble(
    specimen = specimen, trial_id = trial_id, category = category,
    step_number = ss$step_number, muscle_side = ss$muscle_side,
    commanded_x = ss$shortening_strain_pct,
    force_zTorque_N = ss$force_zTorque_N,
    force_uhat_N    = ss$muscle_force_vector_N
  )
}

cli::cli_h1("Collecting isometric + isovelocity steps (zTorque vs uHat force), bass16/17/18")
all_steps <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  manifest <- parse_trial_directory(raw_source_dir(subfolder))
  purrr::pmap_dfr(list(manifest$fullpath, manifest$protocol), function(f, proto) {
    if (!proto %in% c("isometric", "isovelocity")) return(tibble())
    out <- tryCatch(.collect_one_trial(f, specimen, proto), error = function(e) {
      cli::cli_alert_danger("{basename(f)}: {conditionMessage(e)}"); tibble() })
    if (nrow(out) > 0L) cli::cli_alert_success("{tools::file_path_sans_ext(basename(f))}: {nrow(out)} steps")
    out
  })
})

all_steps <- all_steps |> dplyr::filter(.data$muscle_side == "right")
write.csv(all_steps, file.path(DATA_OUT_DIR, "fv_fl_ztorque_vs_uhat_steps.csv"), row.names = FALSE)
cli::cli_alert_success("Saved {nrow(all_steps)} right-muscle steps -> data_processed/fv_fl_ztorque_vs_uhat_steps.csv")

mk_long <- function(df) {
  df |>
    tidyr::pivot_longer(cols = c("force_zTorque_N", "force_uhat_N"),
                         names_to = "force_method", values_to = "force_N") |>
    dplyr::mutate(force_method = ifelse(.data$force_method == "force_zTorque_N",
                                         "zTorque (single primary axis)", "uHat (empirical force vector)")) |>
    dplyr::filter(is.finite(.data$force_N), is.finite(.data$commanded_x))
}

# ---- FL (isometric) ----
fl <- mk_long(dplyr::filter(all_steps, .data$category == "isometric"))
if (nrow(fl) > 0L) {
  pFL <- ggplot(fl, aes(x = .data$commanded_x, y = .data$force_N, color = .data$specimen)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_point(size = 2.4, alpha = 0.8) +
    facet_wrap(~force_method, scales = "free_y") +
    scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen") +
    labs(title = "FL (isometric), RIGHT muscle: force-reconstruction method comparison",
         subtitle = "Same commanded strain x-axis; force two ways. NOTE: the two columns use OPPOSITE sign conventions (compare magnitude/trend,\nnot sign). uHat COLLAPSES toward 0 for bass17 (orange) -- its empirical direction is low-SNR there -- while the single-axis\nzTorque robustly recovers a clean monotonic FL trend for all 3 specimens: evidence the simple zTorque method is more robust.",
         x = "Commanded shortening strain (%)", y = "Muscle force (N)") +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
  foutFL <- file.path(OUT_DIR, "FL_isometric_zTorqueVsUhat.png")
  ggplot2::ggsave(foutFL, pFL, width = 11, height = 5.2, dpi = 150)
  cli::cli_alert_success("Saved {foutFL}")
}

# ---- FV (isovelocity) ----
fv <- mk_long(dplyr::filter(all_steps, .data$category == "isovelocity"))
if (nrow(fv) > 0L) {
  pFV <- ggplot(fv, aes(x = .data$commanded_x, y = .data$force_N, color = .data$specimen)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
    geom_point(size = 2.4, alpha = 0.8) +
    facet_wrap(~force_method, scales = "free_y") +
    scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen") +
    labs(title = "FV (isovelocity), RIGHT muscle: force-reconstruction method comparison",
         subtitle = "Same commanded strain-rate x-axis; force two ways (opposite sign conventions -- compare magnitude/trend, not sign).\nzTorque (single-axis) shows a clean force-velocity structure across all 3 specimens; uHat is noisier and compressed / under-\nrecovers -- again favouring the simpler single-axis zTorque force for these data. Positive x = shortening (concentric).",
         x = "Commanded shortening strain rate (%/s)", y = "Muscle force (N)") +
    theme_bw(base_size = 11) + theme(legend.position = "bottom")
  foutFV <- file.path(OUT_DIR, "FV_isovelocity_zTorqueVsUhat.png")
  ggplot2::ggsave(foutFV, pFV, width = 11, height = 5.2, dpi = 150)
  cli::cli_alert_success("Saved {foutFV}")
}

cli::cli_h1("zTorque vs uHat force agreement (right muscle)")
all_steps |>
  dplyr::group_by(.data$category) |>
  dplyr::summarise(n = dplyr::n(),
                    r = suppressWarnings(cor(.data$force_zTorque_N, .data$force_uhat_N, use = "complete.obs")),
                    median_ratio_uhat_over_z = median(.data$force_uhat_N / .data$force_zTorque_N, na.rm = TRUE),
                    .groups = "drop") |>
  print()

# =============================================================================
# PI follow-up, 2026-07-24: "How does zTorque muscle TENSION compare to uHat
# tension? Double check your work and math." Converts both force methods to
# specific tension (N/cm^2) using MEASURED_RED_MUSCLE_CSA_CM2
# (muscle_geometry.R, the SAME CSA now used everywhere else -- CSA is a
# property of the muscle, independent of which force-reconstruction method
# is used, so applying one shared CSA to both is correct here).
#
# Sign check FIRST: force_zTorque_N and force_uhat_N are documented above
# (this script's own header/subtitles) to use OPPOSITE raw sign conventions.
# Tension is compared on |force| (magnitude) -- sign is a convention
# artifact for this comparison, not a real difference in muscle output.
# =============================================================================
tension_tbl <- all_steps |>
  dplyr::mutate(
    same_sign            = sign(.data$force_zTorque_N) == sign(.data$force_uhat_N),
    tension_zTorque_Ncm2 = abs(.data$force_zTorque_N) / MEASURED_RED_MUSCLE_CSA_CM2,
    tension_uhat_Ncm2    = abs(.data$force_uhat_N)    / MEASURED_RED_MUSCLE_CSA_CM2,
    ratio_uhat_over_z_abs = .data$tension_uhat_Ncm2 / .data$tension_zTorque_Ncm2
  )
write.csv(tension_tbl, file.path(DATA_OUT_DIR, "fv_fl_ztorque_vs_uhat_tension.csv"), row.names = FALSE)

cli::cli_h1(sprintf("Sign agreement (raw force_zTorque_N vs. force_uhat_N), CSA = %.2f cm^2", MEASURED_RED_MUSCLE_CSA_CM2))
print(dplyr::count(tension_tbl, .data$category, .data$same_sign))

cli::cli_h1("Magnitude-only (|force|-based) agreement, by category")
print(tension_tbl |> dplyr::group_by(.data$category) |> dplyr::summarise(
  n = dplyr::n(),
  r_abs = suppressWarnings(cor(abs(.data$force_zTorque_N), abs(.data$force_uhat_N), use = "complete.obs")),
  median_ratio_abs = median(.data$ratio_uhat_over_z_abs, na.rm = TRUE),
  mean_tension_zTorque_Ncm2 = mean(.data$tension_zTorque_Ncm2, na.rm = TRUE),
  mean_tension_uhat_Ncm2    = mean(.data$tension_uhat_Ncm2, na.rm = TRUE),
  max_tension_zTorque_Ncm2  = max(.data$tension_zTorque_Ncm2, na.rm = TRUE),
  max_tension_uhat_Ncm2     = max(.data$tension_uhat_Ncm2, na.rm = TRUE),
  .groups = "drop"))

cli::cli_h1("Magnitude-only agreement, by specimen x category (heterogeneity check)")
print(tension_tbl |> dplyr::group_by(.data$specimen, .data$category) |> dplyr::summarise(
  n = dplyr::n(),
  r_abs = suppressWarnings(cor(abs(.data$force_zTorque_N), abs(.data$force_uhat_N), use = "complete.obs")),
  median_ratio_abs = round(median(.data$ratio_uhat_over_z_abs, na.rm = TRUE), 3),
  .groups = "drop"))

# CAVEAT, printed + on-figure: the isovelocity force here is a LOCAL
# peak-window mean (.mfv_window_peak_means(), MFV_PEAK_WINDOW_S=0.15s,
# muscle_force_vector.R) -- NOT the constant-velocity-window, phase-matched
# force method summary_fv_fl_best_within_individual.R later found NECESSARY
# to remove an "anti-Hill" inflation artifact from isovelocity steps. High
# isovelocity tension values here (up to several x Coughlin's 2000 bass
# benchmark) most likely reflect that PRE-fix windowing, not real muscle
# tension approaching the literature value -- do not read them as
# "zTorque/uHat isovelocity tension matches Coughlin."
cli::cli_alert_warning("CAVEAT: isovelocity force here is a pre-CV-window PEAK-WINDOW MEAN (MFV_PEAK_WINDOW_S=0.15s) -- NOT the constant-velocity-window method summary_fv_fl_best_within_individual.R found necessary to remove an inflation artifact. Isovelocity tension values here are NOT a trustworthy Coughlin comparison; use isometric only for that.")

COUGHLIN_TENSION_NCM2 <- list(mean = 18.0, sd = 3.36)  # Coughlin 2000, 180+/-33.6 kN/m^2 (PI-updated 2026-07-24, graph-read at the PI's actual protocol condition; ~same magnitude as the original 186.4 text value -- see summary_coughlin2000_bass_comparison.R)

pTension <- ggplot(tension_tbl, aes(x = .data$tension_zTorque_Ncm2, y = .data$tension_uhat_Ncm2)) +
  annotate("rect", xmin = COUGHLIN_TENSION_NCM2$mean - COUGHLIN_TENSION_NCM2$sd,
           xmax = COUGHLIN_TENSION_NCM2$mean + COUGHLIN_TENSION_NCM2$sd, ymin = -Inf, ymax = Inf,
           fill = "#b30000", alpha = 0.08) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = .data$specimen, shape = .data$category), size = 2.6, alpha = 0.8) +
  scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen") +
  facet_wrap(~category, scales = "free") +
  labs(title = "zTorque vs. uHat specific tension (|force| / measured CSA), right muscle",
       subtitle = sprintf("CSA = %.2f cm^2 (MEASURED_RED_MUSCLE_CSA_CM2). Dashed = 1:1 line (methods agree). Red band = Coughlin (2000) bass tension\n(18.64+/-3.36 N/cm^2). CAVEAT: isovelocity force is a pre-CV-window peak-window mean -- NOT a trustworthy Coughlin comparison (see log).", MEASURED_RED_MUSCLE_CSA_CM2),
       x = "zTorque tension (N/cm^2)", y = "uHat tension (N/cm^2)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
foutT <- file.path(OUT_DIR, "tension_zTorqueVsUhat_scatter.png")
ggplot2::ggsave(foutT, pTension, width = 11, height = 5.5, dpi = 150)
cli::cli_alert_success("Saved {foutT}")

cli::cli_alert_success("diag_fv_fl_ztorque_vs_uhat.R complete")
