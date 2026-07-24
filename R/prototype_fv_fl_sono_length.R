# prototype_fv_fl_sono_length.R
# PI request 2A, 2026-07-24: prototype FV (isovelocity) and FL (isometric)
# curves computed from the REAL muscle length (sonomicrometry) instead of the
# commanded operating_point, to test whether the actual muscle
# length/velocity gives cleaner curves than the commanded target.
#
# What changes vs. the production pipeline (run_fv_fl_power_pipeline.R):
#   - FORCE (y-axis) is UNCHANGED: muscle_force_Nm from build_segmented_step_
#     summary() (active-minus-passive windowed torque on the single primary-
#     bending-axis channel -- see the z-torque/uHat clarification in the
#     decision log; this is the "z-torque" method). 2B compares force methods.
#   - LENGTH / VELOCITY (x-axis) is recomputed from sono:
#       * FL (isometric): x = mean 40 Hz-filtered sono strain during the
#         step's active window, sign-folded by contraction_mode EXACTLY as
#         attach_step_measured_strain() folds the encoder strain -- directly
#         comparable to the commanded shortening_strain_pct (%).
#       * FV (isovelocity): x = sono strain RATE (%/s) = OLS slope of the
#         filtered sono strain vs. time across the steady middle 60% of the
#         active window (trimming accel/decel ramps), sign-folded the same way.
#
# KNOWN PROTOTYPE LIMITATION (FV only): the sono strain rate rank-agrees with
# the commanded strain rate (r~0.97) but comes out ~10x SMALLER in magnitude.
# Part of this is expected/real (muscle-tendon-grip series compliance means
# the muscle does not achieve the full rigid-geometry commanded excursion
# rate), but the factor is too large to be compliance alone: the middle-60%-
# of-stim-window heuristic almost certainly still includes non-ramp samples
# (pre/post holds), diluting the slope. Isolating the true constant-velocity
# ramp segment (e.g. from the commanded angle's own velocity plateau) is the
# needed refinement before any sono-based Vmax/FV fit should be trusted. FL
# (isometric) does NOT have this issue -- it is a position, not a rate, and
# sono vs commanded there agree at r~0.99 in BOTH value and scale.
#
# Scope: RIGHT-muscle steps only -- sono instruments only the right muscle,
# so a sono x-axis exists only for right-side steps (left steps keep their
# commanded value in the pipeline, untouched here). Pools bass16/17/18.
#
# Run with:  Rscript R/prototype_fv_fl_sono_length.R
# Outputs -> figs_summary/ (FIGS_SUMMARY_DIR) + per-step CSV.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr); library(ggplot2); library(cli); library(rhdf5); library(signal)
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

OUT_DIR      <- FIGS_SUMMARY_DIR
DATA_OUT_DIR <- file.path(.crittergripper_root(), "02_processed", "data_processed")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SPECIMEN_SUBFOLDERS <- c(bass16 = BASS16_RAW_SUBFOLDER, bass17 = BASS17_RAW_SUBFOLDER, bass18 = BASS18_RAW_SUBFOLDER)
SPECIMEN_COLORS     <- c(bass16 = "#1b9e77", bass17 = "#d95f02", bass18 = "#7570b3")
DEACT_WIN_S         <- 0.5
SONO_LOWPASS_CUTOFF_HZ <- 40.0

.butter_lp_sono <- function(x, cutoff_hz, sample_rate_hz, order = 4L) {
  nyq <- sample_rate_hz / 2.0
  if (cutoff_hz >= nyq) return(x)
  ok <- is.finite(x)
  if (sum(ok) < 20L) return(x)
  filt <- signal::butter(order, cutoff_hz / nyq, type = "low")
  out  <- x; out[ok] <- signal::filtfilt(filt, x[ok]); out
}

.read_sono_geom <- function(filename) {
  h5 <- rhdf5::H5Fopen(filename, "H5F_ACC_RDONLY"); on.exit(try(rhdf5::H5Fclose(h5), silent = TRUE), add = TRUE)
  m <- tryCatch(rhdf5::h5readAttributes(h5, "/metadata"), error = function(e) list())
  d1 <- function(v, def = NA_real_) { v <- suppressWarnings(as.numeric(v[1L])); if (length(v)==0L||is.na(v)) def else v }
  list(ai_rate_hz = d1(m[["daq_ai_sample_rate_hz"]], 1000.0),
       lidx_pos_motor = d1(m[["daq_specimen_lateral_index_on_positive_motor_side"]]),
       lidx_right = d1(m[["daq_specimen_side_index_right"]]))
}

# Sign-fold to match attach_step_measured_strain()'s convention exactly.
.fold_by_mode <- function(raw, mode) {
  dplyr::case_when(
    mode == "concentric"     ~ abs(raw),
    mode == "eccentric"      ~ -abs(raw),
    mode == "isometric_zero" ~ 0.0,
    .default = NA_real_
  )
}

# =============================================================================
# One trial -> step rows with BOTH commanded and sono x-axis.
# =============================================================================
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
  ss  <- res$step_summary
  tdA <- res$td

  # sono strain, filtered, on the analysis td
  geom <- .read_sono_geom(f)
  tdA <- tryCatch(attach_sono_strain(tdA, f, geom$lidx_right, geom$lidx_pos_motor), error = function(e) tdA)
  if (!("sono_right_mm" %in% names(tdA)) || all(!is.finite(tdA$sono_right_mm))) return(tibble())
  L0 <- .sono_reference_length_mm(tdA$angle.deg, tdA$sono_right_mm)
  if (!is.finite(L0) || L0 <= 0) return(tibble())
  sono_filt_mm <- .butter_lp_sono(tdA$sono_right_mm, SONO_LOWPASS_CUTOFF_HZ, geom$ai_rate_hz)
  tdA$strain_sono_filt_pct <- -(sono_filt_mm - L0) / L0 * 100.0

  # per-step sono x-axis
  ss$sono_x <- vapply(seq_len(nrow(ss)), function(i) {
    s <- ss[i, ]
    if (is.na(s$contraction_mode) || !is.finite(s$stim_t0_s) || !is.finite(s$stim_t1_s)) return(NA_real_)
    win <- tdA$step_number == s$step_number & tdA$t.s >= s$stim_t0_s & tdA$t.s <= (s$stim_t1_s + DEACT_WIN_S)
    if (!any(win, na.rm = TRUE)) return(NA_real_)
    strain <- tdA$strain_sono_filt_pct[win]; tt <- tdA$t.s[win]
    ok <- is.finite(strain) & is.finite(tt)
    if (sum(ok) < 5L) return(NA_real_)
    if (category == "isometric") {
      raw <- mean(strain[ok])                 # % (position-like strain the muscle was held at)
    } else {
      # strain RATE (%/s): OLS slope over the steady middle 60% of the window
      tt2 <- tt[ok]; st2 <- strain[ok]
      q <- stats::quantile(tt2, c(0.2, 0.8), names = FALSE)
      mid <- tt2 >= q[1] & tt2 <= q[2]
      if (sum(mid) < 5L) mid <- rep(TRUE, length(tt2))
      raw <- unname(stats::coef(stats::lm(st2[mid] ~ tt2[mid]))[2L])
    }
    .fold_by_mode(raw, s$contraction_mode)
  }, numeric(1L))

  tibble::tibble(
    specimen = specimen, trial_id = trial_id, category = category,
    step_number = ss$step_number, muscle_side = ss$muscle_side,
    contraction_mode = ss$contraction_mode,
    commanded_x = ss$shortening_strain_pct,
    sono_x = ss$sono_x,
    muscle_force_Nm = ss$muscle_force_Nm
  )
}

cli::cli_h1("Collecting isometric + isovelocity steps (commanded vs sono x-axis), bass16/17/18")
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

# RIGHT muscle only (sono limitation), finite force
all_steps <- all_steps |>
  dplyr::filter(.data$muscle_side == "right", is.finite(.data$muscle_force_Nm))
write.csv(all_steps, file.path(DATA_OUT_DIR, "fv_fl_sono_vs_commanded_steps.csv"), row.names = FALSE)
cli::cli_alert_success("Saved {nrow(all_steps)} right-muscle steps -> data_processed/fv_fl_sono_vs_commanded_steps.csv")

# =============================================================================
# Long form: one row per step per x-axis source, for faceted overlay.
# =============================================================================
mk_long <- function(df) {
  df |>
    tidyr::pivot_longer(cols = c("commanded_x", "sono_x"), names_to = "x_source", values_to = "x_val") |>
    dplyr::mutate(x_source = ifelse(.data$x_source == "commanded_x",
                                     "commanded (operating_point)", "sono (measured muscle length)")) |>
    dplyr::filter(is.finite(.data$x_val))
}

# ---- FL (isometric): force vs strain ----
fl <- mk_long(dplyr::filter(all_steps, .data$category == "isometric"))
pFL <- ggplot(fl, aes(x = .data$x_val, y = .data$muscle_force_Nm, color = .data$specimen)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(size = 2.4, alpha = 0.8) +
  facet_wrap(~x_source, scales = "free_x") +
  scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen") +
  labs(title = "FL (isometric), RIGHT muscle: force vs. commanded target vs. real (sono) muscle strain",
       subtitle = "Same force (active-minus-passive primary-axis torque); x-axis swapped from commanded operating_point to the 40 Hz-filtered\nsono strain actually held during each step. Cleaner clustering under sono = the commanded target was not the true length.",
       x = "Shortening strain (%, sign-folded by contraction mode)", y = "Muscle force (N*m)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
foutFL <- file.path(OUT_DIR, "FL_isometric_sonoVsCommanded.png")
ggplot2::ggsave(foutFL, pFL, width = 11, height = 5.2, dpi = 150)
cli::cli_alert_success("Saved {foutFL}")

# ---- FV (isovelocity): force vs strain rate ----
fv <- mk_long(dplyr::filter(all_steps, .data$category == "isovelocity"))
pFV <- ggplot(fv, aes(x = .data$x_val, y = .data$muscle_force_Nm, color = .data$specimen)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
  geom_point(size = 2.4, alpha = 0.8) +
  facet_wrap(~x_source, scales = "free_x") +
  scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen") +
  labs(title = "FV (isovelocity), RIGHT muscle: force vs. commanded velocity vs. real (sono) shortening rate",
       subtitle = "Same force; x-axis swapped from commanded angular velocity to the sono strain RATE (OLS slope over the steady middle 60%% of\neach step). PROTOTYPE CAVEAT: sono and commanded rank-agree (r=0.97) but sono magnitudes are ~10x SMALLER -- partly real\nseries-compliance, partly because this window heuristic does not yet isolate the constant-velocity ramp precisely (see header).",
       x = "Shortening strain rate (%/s, sign-folded by contraction mode)", y = "Muscle force (N*m)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
foutFV <- file.path(OUT_DIR, "FV_isovelocity_sonoVsCommanded.png")
ggplot2::ggsave(foutFV, pFV, width = 11, height = 5.2, dpi = 150)
cli::cli_alert_success("Saved {foutFV}")

# Quick fidelity: commanded vs sono x correlation per protocol
cli::cli_h1("commanded-x vs sono-x agreement (right muscle)")
all_steps |>
  dplyr::group_by(.data$category) |>
  dplyr::summarise(n = dplyr::n(),
                    r = suppressWarnings(cor(.data$commanded_x, .data$sono_x, use = "complete.obs")),
                    .groups = "drop") |>
  print()

cli::cli_alert_success("prototype_fv_fl_sono_length.R complete")
