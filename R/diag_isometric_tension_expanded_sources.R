# diag_isometric_tension_expanded_sources.R
# PI directive (2026-07-25): "You MUST include the isometric pre-post stim
# in dynamic as well as the V=0 form isovelocity [in the isometric specific-
# tension analysis], provided they pass all the other inclusion criteria."
#
# The isometric-protocol tension panel (diag_precondition_tension_vs_offset_
# isometric.R) only draws from n=5 isometric-protocol FILES (n=4 after the
# "later (stable)" precondition filter). This script adds TWO more genuine
# isometric-CONTRACTION sources that live inside other protocols' files:
#
#   isovelocity_v0     Every isovelocity trial commands operating_point==0
#                       deg/s ("hold still") for BOTH unilateral sides (plus
#                       a non-stim bilateral "calibration" step at v=0) --
#                       confirmed directly in the raw metadata, 11 files
#                       across bass16/17/18. A v=0 commanded step IS an
#                       isometric contraction in every physical sense (motor
#                       not moving); analyze_isovelocity() (03_analyze.R)
#                       ALREADY computes muscle_force_Nm for these rows with
#                       its normal (PI-validated, unswapped) muscle_side ->
#                       force_sign mapping -- reused here completely
#                       UNCHANGED, just filtered down to the v0 rows.
#   dynamic_l0_bookend  Every dynamic (single_finite) trial brackets its
#                       active cycling with up to 4 static L0 (commanded 0
#                       deg) stimulated bursts -- left+right BEFORE cycling
#                       starts, left+right AFTER it ends (the "isometric
#                       pre-/post-stim" the PI is referring to; "before
#                       t=0" / "after the final bending cycle"). Detected
#                       via the EXISTING detect_dynamic_l0_bookends()
#                       (extract_dynamic_l0_bookends.R, unmodified) --
#                       already used in production to feed plot_fatigue_
#                       timeline.R, just not previously folded into a
#                       tension comparison.
#
# SIGN-CONVENTION CAVEAT -- TESTED EMPIRICALLY, NOT ASSUMED (2026-07-25):
# the FIRST version of this script assumed detect_dynamic_l0_bookends()'s
# muscle_side ("left"/"right") needed the SAME row_side->lidx SWAP already
# confirmed/fixed for dynamic CYCLING bursts (analysis_muscle_force_vector_
# log.md, 2026-07-25 "root cause found" addendum) -- reasoning that both
# come from the same per-sample stim_side detection on the same dynamic
# file. Directly tested against 36 bookends across 9 files, all 3 specimens
# (bass16/17/18): the UNSWAPPED mapping (side=="left" -> lidx_left, SAME
# convention as genuine isometric/isovelocity's index_step_recruitment-based
# muscle_side) gives 36/36 POSITIVE muscle_force_Nm (physiologically
# correct -- a real stimulated contraction producing more torque than its
# own resting baseline); the SWAPPED mapping gives 36/36 NEGATIVE. This
# REFINES (narrows) the original bug finding: the confirmed stim_side
# mislabeling is specific to make_stimuli()'s CYCLING-phase sine-burst
# logic (which extremum of the oscillation a burst is centered on), NOT to
# the separate pre-/post-stim BOOKEND triggering (a static, non-cyclic
# firing with no phase-to-extremum computation involved) -- the two are
# apparently different code paths inside make_stimuli() with different
# (in)correctness. Bookend rows below therefore use the SAME UNSWAPPED
# mapping as isometric_v0/isovelocity_v0, not a special-cased swap.
#
# CORRECTION to a hypothesis raised while building this (see the empirical
# test above): muscle_force_vector.R::.mfv_finalize_step()'s force_zTorque_N
# side-correction, and run_fv_fl_power_pipeline.R's L0-bookend block
# (~line 548) feeding plot_fatigue_timeline.R, ALREADY use this SAME
# unswapped mapping for bookend rows -- now CONFIRMED correct, not a bug.
#
# METHOD CONSISTENCY: uses the SAME scalar zTorque static-pre-baseline
# method as the production isometric tension calc (passive_force_Nm_static
# equivalent) for ALL THREE sources, rather than mixing in the separate
# 6-axis vector-force method already computed for L0 bookends elsewhere --
# the 2026-07-25 baseline-method-sensitivity diagnostic (diag_isometric_
# passive_baseline_methods.R) already showed the specific WINDOW choice
# barely matters, so using one consistent, already-vetted method across all
# 3 sources avoids an apples-to-oranges mix of two different underlying
# force-extraction algorithms.
#
# Run with:  Rscript R/diag_isometric_tension_expanded_sources.R
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
src("plot_force_vs_time.R")           # .smooth_trace_display_only(), .detect_stim_events()
src("fit_fv_fl.R")                    # fit_force_velocity_curve() -- REQUIRED by analyze_isovelocity()
src("03_analyze.R")
src("extract_dynamic_l0_bookends.R")  # detect_dynamic_l0_bookends() -- REUSED UNCHANGED
src("parse_trial_filename.R")
src("dynamic_trial_precondition.R")   # extract_bender_trial_num(), classify_session_precondition()

OUT_DIR      <- FIGS_DIAGNOSTIC_DIR
DATA_OUT_DIR <- file.path(.crittergripper_root(), "02_processed", "data_processed")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SPECIMEN_SUBFOLDERS <- c(bass16 = BASS16_RAW_SUBFOLDER, bass17 = BASS17_RAW_SUBFOLDER, bass18 = BASS18_RAW_SUBFOLDER)
SPECIMEN_COLORS     <- c(bass16 = "#1b9e77", bass17 = "#d95f02", bass18 = "#7570b3")
COUGHLIN_TENSION_MEAN_KNM2 <- 180
COUGHLIN_TENSION_SD_KNM2   <- 33.6

.load_one <- function(f) {
  td <- load_bender_flat(f, do_filter = TRUE, loadtorques = "x")
  tau <- deconvolve_bender(f, hub_path = NULL, verbose = FALSE)
  N <- min(nrow(td), length(tau)); td <- td[seq_len(N), , drop = FALSE]
  td$torque_inertia_corrected_Nm <- tau[seq_len(N)]
  attr(td, "Filename") <- f
  td
}

# Same full geometry reader as run_fv_fl_power_pipeline.R / .read_segmented_
# step_geometry() (03_analyze.R) -- duplicated here per this repo's
# established per-script pattern (see those files' own copies) rather than
# exporting a shared one.
.read_file_level_geometry <- function(filename) {
  h5 <- rhdf5::H5Fopen(filename, "H5F_ACC_RDONLY")
  on.exit(try(rhdf5::H5Fclose(h5), silent = TRUE), add = TRUE)
  m_attrs <- tryCatch(rhdf5::h5readAttributes(h5, "/metadata"), error = function(e) list())
  dbl1 <- function(v, default = NA_real_) {
    v <- suppressWarnings(as.numeric(v[1L]))
    if (length(v) == 0L || is.na(v)) default else v
  }
  local_body_width_mm  <- dbl1(m_attrs[["measurement_specimen_local_body_width_millimeter"]])
  local_body_height_mm <- dbl1(m_attrs[["measurement_specimen_local_body_height_millimeter"]])
  dclamp_mm            <- dbl1(m_attrs[["measurement_clamp_separation_millimeter"]])
  density_g_per_mm3    <- dbl1(m_attrs[["measurement_specimen_density_gram_per_cubic_millimeter"]])
  muscle <- compute_muscle_mass_and_csa(local_body_width_mm, local_body_height_mm, dclamp_mm, density_g_per_mm3)
  list(
    local_body_width_mm  = local_body_width_mm,
    muscle_depth_mm_raw  = dbl1(m_attrs[["measurement_target_muscle_depth_millimeter"]]),
    dclamp_mm            = dclamp_mm,
    muscle               = muscle,
    lidx_pos_motor       = dbl1(m_attrs[["daq_specimen_lateral_index_on_positive_motor_side"]]),
    lidx_left            = dbl1(m_attrs[["daq_specimen_side_index_left"]]),
    lidx_right           = dbl1(m_attrs[["daq_specimen_side_index_right"]])
  )
}

.file_r_m <- function(geom) {
  depth <- resolve_muscle_depth_mm(geom$muscle_depth_mm_raw)
  compute_predicted_strain(1.0, geom$local_body_width_mm, depth$depth_mm)$r_m
}

.time_window_mean <- function(td, t0, t1) {
  if (!is.finite(t0) || !is.finite(t1)) return(NA_real_)
  rows <- td$t.s >= t0 & td$t.s <= t1
  if (!any(rows, na.rm = TRUE)) return(NA_real_)
  mean(td$torque_inertia_corrected_Nm[rows], na.rm = TRUE)
}

# =============================================================================
# Source 1: genuine isometric-protocol steps (right-side only, matching
# production's own convention) -- IDENTICAL collection to
# diag_precondition_tension_vs_offset_isometric.R, tagged source_protocol.
# =============================================================================
cli::cli_h1("Source 1: isometric-protocol steps")

.collect_isometric <- function(specimen, source_dir) {
  manifest <- parse_trial_directory(source_dir)
  files <- manifest$fullpath[manifest$protocol == "isometric"]
  purrr::map_dfr(files, function(f) {
    trial_id <- tools::file_path_sans_ext(basename(f))
    td <- tryCatch(.load_one(f), error = function(e) NULL)
    if (is.null(td)) return(tibble())
    res <- tryCatch(analyze_isometric(td, f), error = function(e) NULL)
    if (is.null(res) || nrow(res$step_summary) == 0L || !is.finite(res$r_m) || res$r_m <= 0) return(tibble())
    steps <- res$step_summary
    steps$force_N               <- steps$muscle_force_Nm / res$r_m
    steps$specific_tension_Ncm2 <- steps$force_N / MEASURED_RED_MUSCLE_CSA_CM2
    steps$specimen  <- specimen
    steps$trial_id  <- trial_id
    steps$trial_num <- extract_bender_trial_num(trial_id)
    steps$source_protocol <- "isometric"
    dplyr::select(steps, specimen, trial_id, trial_num, source_protocol, muscle_side,
                  shortening_strain_pct, specific_tension_Ncm2)
  })
}

iso_rows <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  .collect_isometric(specimen, raw_source_dir(subfolder))
})
cli::cli_alert_success("isometric: {nrow(iso_rows)} steps (all sides)")

# =============================================================================
# Source 2: isovelocity v=0 steps ("V=0 form isovelocity") -- muscle_side/
# force_sign UNCHANGED from analyze_isovelocity()'s normal (validated)
# computation; only filtered down to the v=0 rows.
# =============================================================================
cli::cli_h1("Source 2: isovelocity v=0 (\"isometric-form\") steps")

.collect_isovelocity_v0 <- function(specimen, source_dir) {
  manifest <- parse_trial_directory(source_dir)
  files <- manifest$fullpath[manifest$protocol == "isovelocity"]
  purrr::map_dfr(files, function(f) {
    trial_id <- tools::file_path_sans_ext(basename(f))
    td <- tryCatch(.load_one(f), error = function(e) NULL)
    if (is.null(td)) return(tibble())
    res <- tryCatch(analyze_isovelocity(td, f), error = function(e) NULL)
    if (is.null(res) || nrow(res$step_summary) == 0L || !is.finite(res$r_m) || res$r_m <= 0) return(tibble())
    steps <- res$step_summary
    steps <- dplyr::filter(steps, abs(.data$shortening_value) < 1e-6, .data$muscle_side %in% c("left", "right"))
    if (nrow(steps) == 0L) return(tibble())
    steps$force_N               <- steps$muscle_force_Nm / res$r_m
    steps$specific_tension_Ncm2 <- steps$force_N / MEASURED_RED_MUSCLE_CSA_CM2
    steps$specimen  <- specimen
    steps$trial_id  <- trial_id
    steps$trial_num <- extract_bender_trial_num(trial_id)
    steps$source_protocol <- "isovelocity_v0"
    steps$shortening_strain_pct <- 0.0  # v=0 by construction -- x-axis anchor point
    dplyr::select(steps, specimen, trial_id, trial_num, source_protocol, muscle_side,
                  shortening_strain_pct, specific_tension_Ncm2)
  })
}

isovel_rows <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  .collect_isovelocity_v0(specimen, raw_source_dir(subfolder))
})
cli::cli_alert_success("isovelocity_v0: {nrow(isovel_rows)} steps (all sides)")

# =============================================================================
# Source 3: dynamic L0 bookends ("isometric pre-/post-stim in dynamic") --
# NEW aggregation; force_sign uses the SWAPPED mapping (see header caveat).
# =============================================================================
cli::cli_h1("Source 3: dynamic L0 bookends (pre-cycling / post-cycling isometric bursts)")

.collect_dynamic_l0_bookends <- function(specimen, source_dir) {
  manifest <- parse_trial_directory(source_dir)
  files <- manifest$fullpath[manifest$protocol == "dynamic"]
  purrr::map_dfr(files, function(f) {
    trial_id <- tools::file_path_sans_ext(basename(f))
    td <- tryCatch(.load_one(f), error = function(e) NULL)
    if (is.null(td)) return(tibble())
    geom <- tryCatch(.read_file_level_geometry(f), error = function(e) NULL)
    if (is.null(geom) || !is.finite(geom$lidx_pos_motor) || !is.finite(geom$lidx_left) || !is.finite(geom$lidx_right)) {
      return(tibble())
    }
    r_m <- .file_r_m(geom)
    if (!is.finite(r_m) || r_m <= 0) return(tibble())

    bookends <- tryCatch(detect_dynamic_l0_bookends(td), error = function(e) tibble())
    if (nrow(bookends) == 0L) return(tibble())

    rows <- lapply(seq_len(nrow(bookends)), function(i) {
      b <- bookends[i, ]
      active_Nm  <- .time_window_mean(td, b$stim_t0_s, b$stim_t1_s + 0.5)
      passive_Nm <- .time_window_mean(td, b$t_pre_baseline_start_s, b$t_pre_baseline_end_s)
      if (!is.finite(active_Nm) || !is.finite(passive_Nm)) return(NULL)
      # UNSWAPPED mapping -- SAME convention as isometric/isovelocity
      # (resolve_step_contraction(), muscle_geometry.R). Empirically verified
      # (header comment) that bookends do NOT need the cycling-burst swap:
      # 36/36 test bookends gave physiologically-correct positive force with
      # this direct mapping, 0/36 with the swap.
      rec_lidx <- dplyr::case_when(
        b$muscle_side == "left"  ~ geom$lidx_left,
        b$muscle_side == "right" ~ geom$lidx_right,
        .default = NA_real_
      )
      force_sign <- rec_lidx * geom$lidx_pos_motor
      if (!is.finite(force_sign)) return(NULL)
      muscle_force_Nm <- force_sign * (active_Nm - passive_Nm)
      tibble::tibble(
        specimen = specimen, trial_id = trial_id, trial_num = extract_bender_trial_num(trial_id),
        source_protocol = "dynamic_l0_bookend", muscle_side = b$muscle_side,
        shortening_strain_pct = 0.0,
        specific_tension_Ncm2 = (muscle_force_Nm / r_m) / MEASURED_RED_MUSCLE_CSA_CM2
      )
    })
    dplyr::bind_rows(rows)
  })
}

dynbookend_rows <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  .collect_dynamic_l0_bookends(specimen, raw_source_dir(subfolder))
})
cli::cli_alert_success("dynamic_l0_bookend: {nrow(dynbookend_rows)} bursts (all sides)")

# =============================================================================
# Combine, filter (right side + "later (stable)" -- SAME criteria as the
# production isometric tension panel), aggregate, plot.
# =============================================================================
all_rows <- dplyr::bind_rows(iso_rows, isovel_rows, dynbookend_rows)
all_rows$precondition <- classify_session_precondition(all_rows$specimen, all_rows$trial_num)

write.csv(all_rows, file.path(DATA_OUT_DIR, "isometric_tension_expanded_sources_allsteps.csv"), row.names = FALSE)

qualifying <- all_rows |>
  dplyr::filter(.data$muscle_side == "right", .data$precondition == "later (stable)",
                is.finite(.data$specific_tension_Ncm2))

cli::cli_h2("Right-side, \"later (stable)\" steps/bursts by source")
print(qualifying |> dplyr::count(.data$source_protocol), n = 10)

trial_tension_expanded <- qualifying |>
  dplyr::group_by(.data$specimen, .data$trial_id, .data$trial_num, .data$source_protocol) |>
  dplyr::summarise(n = dplyr::n(),
                    mean_tension_Ncm2 = mean(.data$specific_tension_Ncm2, na.rm = TRUE),
                    max_tension_Ncm2  = max(abs(.data$specific_tension_Ncm2), na.rm = TRUE),
                    .groups = "drop")
write.csv(trial_tension_expanded, file.path(DATA_OUT_DIR, "isometric_tension_expanded_sources_trialtension.csv"), row.names = FALSE)
cli::cli_h2("Trial-level tension, expanded sources (right side, later/stable only)")
print(trial_tension_expanded, n = 60)

old_n  <- sum(trial_tension_expanded$source_protocol == "isometric")
new_n  <- nrow(trial_tension_expanded)
cli::cli_alert_info("n trial/burst-level tension points: isometric-only = {old_n} -> expanded (all 3 sources) = {new_n}")

# =============================================================================
# Plots -> figs_diagnostic/
# =============================================================================
SOURCE_LEVELS <- c("isometric", "isovelocity_v0", "dynamic_l0_bookend")
SOURCE_LABELS <- c(isometric = "isometric (production)",
                   isovelocity_v0 = "isovelocity, v=0 steps",
                   dynamic_l0_bookend = "dynamic, L0 pre/post bookends")
trial_tension_expanded$source_protocol <- factor(trial_tension_expanded$source_protocol,
                                                   levels = SOURCE_LEVELS, labels = SOURCE_LABELS[SOURCE_LEVELS])

p1 <- ggplot(trial_tension_expanded, aes(x = .data$max_tension_Ncm2 * 10, y = .data$source_protocol)) +
  annotate("rect", ymin = -Inf, ymax = Inf,
           xmin = COUGHLIN_TENSION_MEAN_KNM2 - COUGHLIN_TENSION_SD_KNM2,
           xmax = COUGHLIN_TENSION_MEAN_KNM2 + COUGHLIN_TENSION_SD_KNM2,
           fill = "grey70", alpha = 0.3) +
  geom_vline(xintercept = COUGHLIN_TENSION_MEAN_KNM2, linetype = "dashed", color = "grey40") +
  geom_jitter(aes(color = .data$specimen), height = 0.15, size = 2.6, alpha = 0.85) +
  scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen") +
  labs(title = "Isometric-form specific tension: expanded sources vs. isometric-protocol-only",
       subtitle = paste0("Trial/burst-max |specific tension|, right side, \"later (stable)\" only -- ",
                          "n=", old_n, " (isometric only) -> n=", new_n, " (all 3 sources).\n",
                          "Grey ribbon/dashed line = Coughlin (2000), ", COUGHLIN_TENSION_MEAN_KNM2,
                          " +/- ", COUGHLIN_TENSION_SD_KNM2, " kN/m^2."),
       x = "Trial/burst-max |specific tension| (kN/m^2)", y = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "right")
fout1 <- file.path(OUT_DIR, "isometric_tension_expanded_sources_trialmax.png")
ggplot2::ggsave(fout1, p1, width = 10, height = 5.5, dpi = 150)
cli::cli_alert_success("Saved {fout1}")

qualifying$source_protocol <- factor(qualifying$source_protocol, levels = SOURCE_LEVELS, labels = SOURCE_LABELS[SOURCE_LEVELS])
p2 <- ggplot(qualifying, aes(x = .data$source_protocol, y = abs(.data$specific_tension_Ncm2) * 10)) +
  annotate("rect", xmin = -Inf, xmax = Inf,
           ymin = COUGHLIN_TENSION_MEAN_KNM2 - COUGHLIN_TENSION_SD_KNM2,
           ymax = COUGHLIN_TENSION_MEAN_KNM2 + COUGHLIN_TENSION_SD_KNM2,
           fill = "grey70", alpha = 0.3) +
  geom_hline(yintercept = COUGHLIN_TENSION_MEAN_KNM2, linetype = "dashed", color = "grey40") +
  geom_boxplot(aes(fill = .data$source_protocol), outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_jitter(aes(color = .data$specimen), width = 0.15, size = 1.4, alpha = 0.6) +
  scale_fill_discrete(guide = "none") +
  scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen") +
  scale_y_log10() +
  labs(title = "Isometric-form specific tension distribution by source (log scale)",
       subtitle = paste0("All qualifying (right side, \"later (stable)\") steps/bursts pooled, bass16/17/18. Grey ribbon/dashed = Coughlin (2000), ",
                          COUGHLIN_TENSION_MEAN_KNM2, " +/- ", COUGHLIN_TENSION_SD_KNM2, " kN/m^2."),
       x = NULL, y = "|Specific tension| (kN/m^2, log scale)") +
  theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 15, hjust = 1))
fout2 <- file.path(OUT_DIR, "isometric_tension_expanded_sources_distribution.png")
ggplot2::ggsave(fout2, p2, width = 10, height = 6, dpi = 150)
cli::cli_alert_success("Saved {fout2}")

cli::cli_alert_success("diag_isometric_tension_expanded_sources.R complete -- outputs in {OUT_DIR}/ and {DATA_OUT_DIR}/")
