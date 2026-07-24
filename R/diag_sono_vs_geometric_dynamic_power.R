# diag_sono_vs_geometric_dynamic_power.R
# PI direction, 2026-07-24 (reframing the preconditioning investigation):
# "The early dynamic trials are showing REAL muscle behavior, just at
# lengths that were not prescribed by the motor and geometry. But we can
# still extract the muscle mechanics (power) using the lengths measured
# using sonomicrometry for those. So why don't we run a sono-based vs
# encoder-based calculation of muscle power for dynamic trials."
#
# Directly tests that hypothesis: computes per-cycle dynamic muscle power
# TWO ways from the SAME muscle_torque.Nm (force is a torque-sensor
# measurement, independent of any length assumption):
#   (1) GEOMETRIC (the pipeline's current method, unchanged) -- velocity
#       from the COMMANDED-angle kinematics (pos.rad = angle.deg*pi/180,
#       the same quantity 03_analyze.R's add_muscle_instantaneous() uses).
#   (2) SONO -- velocity from the DIRECTLY MEASURED muscle length
#       (sono_right_mm, attach_sono_strain()), CONDITIONED before
#       differentiating (PI-directed, 2026-07-24, following up on the
#       first pass of this script which differentiated the raw decimated
#       signal and got visibly noisy peak-power estimates): zero-phase
#       Butterworth low-pass at 40 Hz (order 4, filtfilt -- same recipe
#       already compared/vetted in diag_sono_smoothing.R, "below DS3
#       Nyquist, preserves muscle motion"; every dynamic cycling frequency
#       in this corpus is well under 10 Hz) applied to the FULL-TRIAL
#       sono_right_mm at its native 1 kHz AI rate (zero-phase = no lag,
#       and filtering the continuous full-trial series avoids the edge
#       artifacts a filtfilt would produce if applied separately to each
#       disjoint stim-window slice), THEN decimated to the true DS3 update
#       rate (daq_sono_internal_sample_rate_hertz, ~241-247 Hz -- "decimate,
#       don't dedupe", the fix already adopted in plot_sono_strain_
#       validation_pooled.R) before differentiating. This removes both the
#       1 kHz AI staircase/oversampling AND genuine sensor/ADC noise above
#       the real muscle-motion bandwidth, which the raw-signal derivative
#       in the first pass was amplifying into spurious peak-power spikes.
#       Force_N = muscle_torque.Nm / r_m (SAME lever arm used everywhere
#       else in this pipeline -- unavoidable, since there is no muscle-
#       specific force transducer, only a whole-body torque sensor).
#
# Scope: RIGHT-STIM cycles ONLY (ALL trials, early AND later -- the whole
# point is to test whether sono recovers valid mechanics from the
# "early/preconditioning" trials that dynamic_trial_precondition.R
# currently excludes). Sono only instruments the right muscle, so
# left-stim cycles have no independent length measurement and are dropped
# here (they still exist, unaffected, in the geometric-only analyses
# elsewhere in this pipeline).
#
# Run with:  Rscript R/diag_sono_vs_geometric_dynamic_power.R
# Outputs -> figures: 02_processed/figs_diagnostic/ (FIGS_DIAGNOSTIC_DIR)
#            data:    02_processed/data_processed/

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr); library(ggplot2); library(cli); library(rhdf5); library(readr); library(signal)
})

# Zero-phase Butterworth LP filter via filtfilt -- same recipe as
# diag_sono_smoothing.R's .butter_lp(), duplicated here (that script is a
# standalone analysis, not a sourceable library).
.butter_lp_sono <- function(x, cutoff_hz, sample_rate_hz, order = 4L) {
  nyq <- sample_rate_hz / 2.0
  if (cutoff_hz >= nyq) return(x)
  ok <- is.finite(x)
  if (sum(ok) < 20L) return(x)
  filt <- signal::butter(order, cutoff_hz / nyq, type = "low")
  out  <- x
  out[ok] <- signal::filtfilt(filt, x[ok])
  out
}
SONO_LOWPASS_CUTOFF_HZ <- 40.0

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")
src("00_load_bender_flat.R")
src("01_calibrate.R")
src("02_deconvolve.R")
src("muscle_geometry.R")
src("03_analyze.R")
src("parse_trial_filename.R")
src("plot_force_vs_time.R")           # .detect_stim_events(), RELAXATION_WINDOW_S
src("plot_strain_validation.R")       # attach_measured_strain()
src("plot_angle_sono_validation.R")   # attach_sono_strain(), .read_sono_right_mm_aligned()
src("dynamic_trial_precondition.R")

OUT_DIR      <- FIGS_DIAGNOSTIC_DIR
DATA_OUT_DIR <- file.path(.crittergripper_root(), "02_processed", "data_processed")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SPECIMEN_SUBFOLDERS <- c(bass16 = BASS16_RAW_SUBFOLDER, bass17 = BASS17_RAW_SUBFOLDER, bass18 = BASS18_RAW_SUBFOLDER)
SPECIMEN_COLORS     <- c(bass16 = "#1b9e77", bass17 = "#d95f02", bass18 = "#7570b3")

# =============================================================================
# Per-file geometry (duplicated from run_fv_fl_power_pipeline.R /
# diag_precondition_power_check.R -- same deliberate-duplication pattern
# those two note; extended here with local_body_width_mm/muscle_depth_mm
# since attach_predicted_strain()/attach_sono_strain() need them too, not
# just compute_muscle_mass_and_csa()).
# =============================================================================
.read_file_level_geometry_full <- function(filename) {
  h5 <- rhdf5::H5Fopen(filename, "H5F_ACC_RDONLY")
  on.exit(try(rhdf5::H5Fclose(h5), silent = TRUE), add = TRUE)
  m_attrs <- tryCatch(rhdf5::h5readAttributes(h5, "/metadata"), error = function(e) list())
  dbl1 <- function(v, default = NA_real_) {
    v <- suppressWarnings(as.numeric(v[1L]))
    if (length(v) == 0L || is.na(v)) default else v
  }
  local_body_width_mm  <- dbl1(m_attrs[["measurement_specimen_local_body_width_millimeter"]])
  measured_depth_mm    <- dbl1(m_attrs[["measurement_target_muscle_depth_millimeter"]])
  dclamp_mm            <- dbl1(m_attrs[["measurement_clamp_separation_millimeter"]])
  local_body_height_mm <- dbl1(m_attrs[["measurement_specimen_local_body_height_millimeter"]])
  density_g_per_mm3    <- dbl1(m_attrs[["measurement_specimen_density_gram_per_cubic_millimeter"]])
  sono_rate_hz         <- dbl1(m_attrs[["daq_sono_internal_sample_rate_hertz"]], default = 241.0)
  ai_rate_hz           <- dbl1(m_attrs[["daq_ai_sample_rate_hz"]], default = 1000.0)
  muscle <- compute_muscle_mass_and_csa(local_body_width_mm, local_body_height_mm, dclamp_mm, density_g_per_mm3)
  list(muscle = muscle,
       local_body_width_mm = local_body_width_mm,
       measured_depth_mm   = measured_depth_mm,
       sono_rate_hz        = sono_rate_hz,
       ai_rate_hz          = ai_rate_hz,
       lidx_pos_motor = dbl1(m_attrs[["daq_specimen_lateral_index_on_positive_motor_side"]]),
       lidx_left      = dbl1(m_attrs[["daq_specimen_side_index_left"]]),
       lidx_right     = dbl1(m_attrs[["daq_specimen_side_index_right"]]))
}

.load_one <- function(f) {
  td <- load_bender_flat(f, do_filter = TRUE, loadtorques = "x")
  tau <- deconvolve_bender(f, hub_path = NULL, verbose = FALSE)
  N <- min(nrow(td), length(tau)); td <- td[seq_len(N), , drop = FALSE]
  td$torque_inertia_corrected_Nm <- tau[seq_len(N)]
  attr(td, "Filename") <- f
  td
}

# =============================================================================
# Per-trial: RIGHT-STIM active-cycle samples, force-signed torque, +
# sono_right_mm carried along row-for-row, decimated to the true DS3 rate.
# =============================================================================
.collect_one_trial <- function(f, specimen) {
  trial_id <- tools::file_path_sans_ext(basename(f))
  td <- tryCatch(.load_one(f), error = function(e) NULL)
  if (is.null(td) || all(as.character(td$stim) == "0")) return(tibble())

  geom <- .read_file_level_geometry_full(f)
  strain_res <- attach_predicted_strain(td, geom$local_body_width_mm, geom$measured_depth_mm)
  r_m <- strain_res$muscle_strain_r_m[1L]
  td  <- strain_res
  td  <- tryCatch(attach_sono_strain(td, f, geom$lidx_right, geom$lidx_pos_motor), error = function(e) td)
  if (!("sono_right_mm" %in% names(td)) || all(!is.finite(td$sono_right_mm))) return(tibble())
  if (!is.finite(r_m) || r_m <= 0) return(tibble())

  # Condition sono BEFORE any row restriction, on the full continuous
  # trial (filtfilt is zero-phase but assumes a continuous time series --
  # filtering disjoint stim-window slices separately would introduce edge
  # artifacts at every slice boundary). See module header for rationale.
  td$sono_right_mm_filt <- .butter_lp_sono(td$sono_right_mm, SONO_LOWPASS_CUTOFF_HZ, geom$ai_rate_hz)

  td$.row_id <- seq_len(nrow(td))
  td2 <- set_cycle_types(td)
  msc <- tryCatch(
    calc_muscle_torque(td2, torque_col = "torque_inertia_corrected_Nm",
                        include_all_active_samples = TRUE, phase_match_all_rows = TRUE),
    error = function(e) { cli::cli_warn("{trial_id}: calc_muscle_torque failed: {conditionMessage(e)}"); NULL }
  )
  if (is.null(msc) || nrow(msc) == 0L || !(".row_id" %in% names(msc))) return(tibble())

  events <- tryCatch(.detect_stim_events(td), error = function(e) tibble::tibble())
  row_side <- rep(NA_character_, nrow(td))
  if (nrow(events) > 0L) {
    events <- dplyr::arrange(events, .data$stim_t0_s)
    for (i in seq_len(nrow(events))) {
      e <- events[i, ]
      in_win <- td$t.s >= e$stim_t0_s & td$t.s <= (e$stim_t1_s + RELAXATION_WINDOW_S)
      row_side[in_win] <- e$muscle_side
    }
  }
  # RIGHT-STIM ONLY -- sono only instruments the right muscle, no
  # left-side length measurement exists to compare against.
  force_sign_right <- geom$lidx_right * geom$lidx_pos_motor
  msc$row_side <- row_side[msc$.row_id]
  msc <- dplyr::filter(msc, .data$row_side == "R")
  if (nrow(msc) == 0L) return(tibble())
  msc$muscle_torque.Nm <- force_sign_right * msc$muscle_torque.Nm
  if ("cycletype" %in% names(msc)) msc <- dplyr::filter(msc, .data$cycletype == "act")
  if (nrow(msc) == 0L) return(tibble())

  msc$sono_right_mm      <- td$sono_right_mm[msc$.row_id]
  msc$sono_right_mm_filt <- td$sono_right_mm_filt[msc$.row_id]
  msc$specimen  <- specimen
  msc$trial_id  <- trial_id
  msc$trial_num <- extract_bender_trial_num(trial_id)
  msc$r_m       <- r_m
  msc$muscle_mass.kg <- geom$muscle$muscle_mass_kg
  msc$sono_rate_hz   <- geom$sono_rate_hz
  msc$ai_rate_hz     <- geom$ai_rate_hz
  dplyr::select(msc, dplyr::any_of(c(
    "specimen", "trial_id", "trial_num", "cycle", "t.s", "pos.rad", "muscle_torque.Nm",
    "sono_right_mm", "sono_right_mm_filt", "r_m", "muscle_mass.kg", "sono_rate_hz", "ai_rate_hz",
    "curvature.invm", "freq.Hz"
  )))
}

cli::cli_h1("Collecting RIGHT-STIM dynamic cycles, all trials (early + later), bass16/17/18")
all_rows <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  manifest <- parse_trial_directory(raw_source_dir(subfolder))
  files <- manifest$fullpath[manifest$protocol == "dynamic"]
  purrr::map_dfr(files, function(f) {
    out <- tryCatch(.collect_one_trial(f, specimen), error = function(e) {
      cli::cli_alert_danger("{basename(f)}: {conditionMessage(e)}"); tibble()
    })
    if (nrow(out) > 0L) cli::cli_alert_success("{tools::file_path_sans_ext(basename(f))}: {nrow(out)} right-stim active samples")
    out
  })
})
all_rows <- all_rows |>
  dplyr::filter(is.finite(.data$trial_num)) |>
  dplyr::mutate(precondition = classify_dynamic_precondition(.data$specimen, .data$trial_num))

# =============================================================================
# Decimate to the true DS3 update rate (per trial -- rate is file-specific,
# ~241-247 Hz across this corpus), THEN differentiate within each cycle to
# get both velocity estimates on the same time base.
# =============================================================================
.decimate_and_diff <- function(df) {
  decim_n <- max(1L, round(df$ai_rate_hz[1L] / df$sono_rate_hz[1L]))
  idx <- seq(1L, nrow(df), by = decim_n)
  d <- df[idx, , drop = FALSE]
  d <- dplyr::arrange(d, .data$cycle, .data$t.s)
  d |>
    dplyr::group_by(.data$cycle) |>
    dplyr::mutate(
      dt_s      = .data$t.s - dplyr::lag(.data$t.s),
      dpos_rad  = .data$pos.rad - dplyr::lag(.data$pos.rad),
      dsono_mm  = .data$sono_right_mm_filt - dplyr::lag(.data$sono_right_mm_filt)
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(is.finite(.data$dt_s), .data$dt_s > 0) |>
    dplyr::mutate(
      omega_geom_rad_s = .data$dpos_rad / .data$dt_s,
      v_sono_m_s        = -( .data$dsono_mm / 1000.0) / .data$dt_s,   # shortening-positive
      force_N           = .data$muscle_torque.Nm / .data$r_m,
      insta_power_geom_W = .data$muscle_torque.Nm * .data$omega_geom_rad_s,
      insta_power_sono_W = .data$force_N * .data$v_sono_m_s
    )
}

decim_df <- all_rows |>
  dplyr::group_by(.data$specimen, .data$trial_id) |>
  dplyr::group_modify(~.decimate_and_diff(.x)) |>
  dplyr::ungroup()

write.csv(decim_df, file.path(DATA_OUT_DIR, "sono_vs_geometric_dynamic_power_persample.csv"), row.names = FALSE)
cli::cli_alert_success("Saved {nrow(decim_df)} decimated right-stim samples -> data_processed/sono_vs_geometric_dynamic_power_persample.csv")

# =============================================================================
# Per-cycle summary (both methods), then per-trial aggregation
# =============================================================================
cyc <- decim_df |>
  dplyr::filter(is.finite(.data$insta_power_geom_W), is.finite(.data$insta_power_sono_W)) |>
  dplyr::group_by(.data$specimen, .data$trial_id, .data$trial_num, .data$precondition, .data$cycle) |>
  dplyr::summarise(
    n_samples          = dplyr::n(),
    muscle_mass.kg     = dplyr::first(.data$muscle_mass.kg),
    avg_power_geom_W   = mean(.data$insta_power_geom_W, na.rm = TRUE),
    peak_power_geom_W  = max(abs(.data$insta_power_geom_W), na.rm = TRUE),
    avg_power_sono_W   = mean(.data$insta_power_sono_W, na.rm = TRUE),
    peak_power_sono_W  = max(abs(.data$insta_power_sono_W), na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    avg_power_geom_Wkg  = .data$avg_power_geom_W  / .data$muscle_mass.kg,
    peak_power_geom_Wkg = .data$peak_power_geom_W / .data$muscle_mass.kg,
    avg_power_sono_Wkg  = .data$avg_power_sono_W  / .data$muscle_mass.kg,
    peak_power_sono_Wkg = .data$peak_power_sono_W / .data$muscle_mass.kg
  )

write.csv(cyc, file.path(DATA_OUT_DIR, "sono_vs_geometric_dynamic_power_percycle.csv"), row.names = FALSE)
cli::cli_alert_success("Saved {nrow(cyc)} cycles -> data_processed/sono_vs_geometric_dynamic_power_percycle.csv")

trial_power <- cyc |>
  dplyr::group_by(.data$specimen, .data$trial_id, .data$trial_num, .data$precondition) |>
  dplyr::summarise(
    n_cycles                = dplyr::n(),
    mean_avg_power_geom_Wkg = mean(.data$avg_power_geom_Wkg, na.rm = TRUE),
    mean_avg_power_sono_Wkg = mean(.data$avg_power_sono_Wkg, na.rm = TRUE),
    max_peak_power_geom_Wkg = max(.data$peak_power_geom_Wkg, na.rm = TRUE),
    max_peak_power_sono_Wkg = max(.data$peak_power_sono_Wkg, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::arrange(.data$specimen, .data$trial_num)

write.csv(trial_power, file.path(DATA_OUT_DIR, "sono_vs_geometric_dynamic_power_bytrial.csv"), row.names = FALSE)
cli::cli_h1("Trial-level mean power: geometric vs. sono, by trial order")
print(trial_power, n = 60)

# =============================================================================
# Plot 1: trial-mean power, geometric vs. sono, by trial order -- does the
# early-trial elevation the geometric method shows also appear (or shrink)
# once real (sono) muscle length is used?
# =============================================================================
cutoff_df <- tibble::tibble(specimen = names(DYNAMIC_PRECONDITION_CUTOFF_TRIALNUM),
                             cutoff = DYNAMIC_PRECONDITION_CUTOFF_TRIALNUM)

trial_long <- trial_power |>
  tidyr::pivot_longer(cols = c("mean_avg_power_geom_Wkg", "mean_avg_power_sono_Wkg"),
                       names_to = "method", values_to = "mean_avg_power_Wkg") |>
  dplyr::mutate(method = ifelse(.data$method == "mean_avg_power_geom_Wkg",
                                 "geometric (commanded-angle kinematics)", "sono (measured muscle length)"))

p1 <- ggplot(trial_long, aes(x = .data$trial_num, y = .data$mean_avg_power_Wkg, color = .data$method, shape = .data$precondition)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(data = cutoff_df, aes(xintercept = cutoff - 0.5), linetype = "dotted", color = "black", inherit.aes = FALSE) +
  geom_line(aes(group = interaction(.data$trial_id, .data$method)), alpha = 0) +
  geom_point(size = 2.6) +
  facet_wrap(~specimen, scales = "free") +
  scale_color_manual(values = c("geometric (commanded-angle kinematics)" = "#7570b3", "sono (measured muscle length)" = "#1b9e77"), name = NULL) +
  scale_shape_manual(values = c("early (preconditioning)" = 17, "later (stable)" = 16), name = NULL) +
  labs(title = "Dynamic, RIGHT-STIM trials: trial-mean mass-specific power, geometric vs. sono-based",
       subtitle = "Same measured torque, two different velocity/length sources. Dotted line = the existing hard early/later cutoff -- shown for\nreference only; this comparison spans ALL trials (early + later), testing whether sono-based power removes the early-trial gap.",
       x = "Trial number (bender_NN, chronological)", y = "Trial-mean cycle-averaged power (W/kg)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
fout1 <- file.path(OUT_DIR, "dynamic_sonoVsGeometric_power_vs_trialorder.png")
ggplot2::ggsave(fout1, p1, width = 11, height = 5, dpi = 150)
cli::cli_alert_success("Saved {fout1}")

# =============================================================================
# Plot 2: per-cycle scatter, geometric vs. sono avg power, colored by
# precondition -- direct agreement check, cycle-by-cycle.
# =============================================================================
lims <- range(c(cyc$avg_power_geom_Wkg, cyc$avg_power_sono_Wkg), na.rm = TRUE)
ref_df <- tibble::tibble(x = lims, y = lims)
cor_lab <- cyc |>
  dplyr::group_by(.data$precondition) |>
  dplyr::summarise(r = suppressWarnings(cor(.data$avg_power_geom_Wkg, .data$avg_power_sono_Wkg, use = "complete.obs")),
                    n = dplyr::n(), .groups = "drop") |>
  dplyr::mutate(label = sprintf("r=%.3f, n=%d", .data$r, .data$n))

p2 <- ggplot(cyc, aes(x = .data$avg_power_geom_Wkg, y = .data$avg_power_sono_Wkg, color = .data$specimen)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_point(size = 1.6, alpha = 0.6) +
  geom_text(data = cor_lab, aes(x = -Inf, y = Inf, label = .data$label), inherit.aes = FALSE,
            hjust = -0.05, vjust = 1.4, size = 3.2) +
  facet_wrap(~precondition) +
  scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen") +
  coord_equal(xlim = lims, ylim = lims) +
  labs(title = "Per-cycle avg power: geometric (commanded-angle) vs. sono (measured length)",
       subtitle = "Each point = one active cycle, right-stim only. Dashed = 1:1 reference.",
       x = "Cycle-averaged power, GEOMETRIC method (W/kg)",
       y = "Cycle-averaged power, SONO method (W/kg)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
fout2 <- file.path(OUT_DIR, "dynamic_sonoVsGeometric_power_scatter.png")
ggplot2::ggsave(fout2, p2, width = 10, height = 5.5, dpi = 150)
cli::cli_alert_success("Saved {fout2}")

# =============================================================================
# Summary table: does the early-vs-later GAP shrink under the sono method?
# =============================================================================
gap_tbl <- trial_power |>
  tidyr::pivot_longer(cols = c("mean_avg_power_geom_Wkg", "mean_avg_power_sono_Wkg"),
                       names_to = "method", values_to = "power_Wkg") |>
  dplyr::group_by(.data$method, .data$precondition) |>
  dplyr::summarise(n_trials = dplyr::n(), median_power_Wkg = median(.data$power_Wkg, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = "precondition", values_from = c("n_trials", "median_power_Wkg"))
cli::cli_h1("Early vs. later median power, geometric vs. sono method")
print(gap_tbl)
write.csv(gap_tbl, file.path(DATA_OUT_DIR, "sono_vs_geometric_dynamic_power_earlyLaterGap.csv"), row.names = FALSE)

# =============================================================================
# Plot 3: boxplot comparison, geometric vs. sono, at CYCLE level (PI-
# requested, 2026-07-24) -- avg AND peak power side by side (2x2: rows =
# {avg, peak}, columns = {geometric, sono}), x = precondition (n-labeled),
# points colored by specimen. cyc has far more rows than trial_power (36
# early / 89 later cycles vs. 14/16 trials), giving boxes with enough
# points to actually show a distribution rather than 2-4 dots.
# =============================================================================
box_long <- cyc |>
  tidyr::pivot_longer(
    cols = c("avg_power_geom_Wkg", "peak_power_geom_Wkg", "avg_power_sono_Wkg", "peak_power_sono_Wkg"),
    names_to = "metric", values_to = "power_Wkg"
  ) |>
  dplyr::mutate(
    method = ifelse(grepl("_geom_", .data$metric), "geometric (commanded-angle)", "sono (measured length, 40 Hz LP-filtered)"),
    stat   = ifelse(grepl("^avg_", .data$metric), "mean (cycle-averaged)", "peak (cycle-max |instantaneous|)")
  ) |>
  dplyr::filter(is.finite(.data$power_Wkg))

n_lab_box <- box_long |>
  dplyr::distinct(.data$specimen, .data$trial_id, .data$cycle, .data$precondition) |>
  dplyr::count(.data$precondition, name = "n_cycles") |>
  dplyr::mutate(x_label = sprintf("%s\n(n=%d cycles)", .data$precondition, .data$n_cycles))
box_long <- dplyr::left_join(box_long, dplyr::select(n_lab_box, "precondition", "x_label"), by = "precondition")
box_long$x_label <- factor(box_long$x_label, levels = dplyr::arrange(n_lab_box, .data$precondition)$x_label)
box_long$method  <- factor(box_long$method, levels = c("geometric (commanded-angle)", "sono (measured length, 40 Hz LP-filtered)"))
box_long$stat    <- factor(box_long$stat, levels = c("mean (cycle-averaged)", "peak (cycle-max |instantaneous|)"))

# PI follow-up, 2026-07-24: "some early trials have much higher cycle power,
# much closer to Coughlin (2000) -- look into that." Coughlin's bass power
# (14.4 +/- 1.9 W/kg, DERIVED as work x frequency at 0.572L, 2.4 L/s -- see
# summary_coughlin2000_bass_comparison.R) is overlaid on every facet so it's
# visible exactly WHERE the apparent "closer to Coughlin" early-trial power
# survives (mean, geometric only) vs. where it does NOT survive sono
# correction (mean, sono method -- median actually goes slightly NEGATIVE).
COUGHLIN_POWER_WKG <- list(mean = 14.4, sd = 1.9)

p3 <- ggplot(box_long, aes(x = .data$x_label, y = .data$power_Wkg)) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = COUGHLIN_POWER_WKG$mean - COUGHLIN_POWER_WKG$sd,
           ymax = COUGHLIN_POWER_WKG$mean + COUGHLIN_POWER_WKG$sd,
           fill = "#b30000", alpha = 0.12) +
  geom_hline(yintercept = COUGHLIN_POWER_WKG$mean, color = "#b30000", linetype = "dashed", linewidth = 0.4) +
  geom_boxplot(outlier.shape = NA, width = 0.5, fill = "grey95", color = "grey40") +
  geom_jitter(aes(color = .data$specimen), width = 0.12, size = 1.8, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.4) +
  facet_grid(stat ~ method, scales = "free_y") +
  scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen") +
  labs(title = "Dynamic, right-stim cycles: mass-specific power, geometric vs. sono method (early vs. later)",
       subtitle = sprintf("Sono velocity from a %.0f Hz zero-phase Butterworth LP filter on sono_right_mm (full trial), decimated to the true DS3\nrate before differentiating -- see module header. Each point = one active right-stim cycle. Red band = Coughlin (2000) bass\npower, derived (1 pt, 0.572L @ 2.4 L/s) -- see summary_coughlin2000_bass_comparison.R for provenance/caveats.", SONO_LOWPASS_CUTOFF_HZ),
       x = NULL, y = "Power (W/kg)") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", strip.text = element_text(size = 9))

fout3 <- file.path(OUT_DIR, "dynamic_sonoVsGeometric_power_boxplot.png")
ggplot2::ggsave(fout3, p3, width = 10, height = 7, dpi = 150)
cli::cli_alert_success("Saved {fout3}")

cli::cli_alert_success("diag_sono_vs_geometric_dynamic_power.R complete -- outputs in {OUT_DIR}/ and {DATA_OUT_DIR}/")
