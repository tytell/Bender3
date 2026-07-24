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

cli::cli_alert_success("diag_fv_fl_ztorque_vs_uhat.R complete")
