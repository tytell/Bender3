# compare_specimen_specific_properties.R
# Cross-specimen comparison of mass-/area-specific muscle properties
# (specific tension N/cm^2, mass-specific power W/kg) for the specimen
# corpora already processed individually by run_fv_fl_power_pipeline.R
# (PI-directed follow-up, 2026-07-16: "compare the normalized values of both
# bass for each protocol type"). READ-ONLY comparison -- does not write
# anything into either specimen's own per-specimen figures directory.
#
# Reuses the SAME calculation functions run_fv_fl_power_pipeline.R uses:
# compute_muscle_mass_and_csa()/add_specific_properties_to_fit()
# (muscle_geometry.R), analyze_isometric()/analyze_isovelocity()/
# summarize_muscle_cycles() (03_analyze.R), fit_force_velocity_curve()
# (fit_fv_fl.R). Isometric FL no longer fits a model by default (PI
# direction, 2026-07-16 -- see fit_fv_fl.R module header); this script
# reports isometric FL descriptively (n points, force range) instead.
# Numbers here are cross-checked against that pipeline's own per-specimen
# console output (2026-07-16 run) and match exactly.
#
# NOTE: .read_file_level_geometry() and .attach_dynamic_muscle_force() below
# are intentionally DUPLICATED from run_fv_fl_power_pipeline.R (that file
# defines them inline for its own per-trial loop and has no importable
# module boundary short of running its whole top-level pipeline). Keep both
# copies in sync if either changes; a follow-up could hoist them into a
# shared module (e.g. 03_analyze.R) if a third caller ever needs them.
#
# Run with: Rscript R/compare_specimen_specific_properties.R

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(stringr)
  library(fs); library(ggplot2); library(cli); library(rhdf5); library(patchwork)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")

# Default comes from paths_config.R (single source of truth) -- see that
# file if the OneDrive folder layout ever moves again.
OUTPUT_DIR <- Sys.getenv("BENDER3_COMPARISON_OUTPUT_DIR", FIGS_SUMMARY_DIR)
fs::dir_create(OUTPUT_DIR, recurse = TRUE)

src("00_load_bender_flat.R")
src("01_calibrate.R")
src("02_deconvolve.R")
src("muscle_geometry.R")
src("fit_fv_fl.R")
src("03_analyze.R")
src("parse_trial_filename.R")
src("plot_force_vs_time.R")  # .detect_stim_events(), RELAXATION_WINDOW_S

`%||%` <- function(x, y) if (is.null(x)) y else x

# Duplicated from run_fv_fl_power_pipeline.R -- see module header note above.
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
  muscle <- compute_muscle_mass_and_csa(local_body_width_mm, local_body_height_mm,
                                        dclamp_mm, density_g_per_mm3)
  list(
    local_body_width_mm  = local_body_width_mm,
    local_body_height_mm = local_body_height_mm,
    muscle_depth_mm_raw  = dbl1(m_attrs[["measurement_target_muscle_depth_millimeter"]]),
    dclamp_mm            = dclamp_mm,
    density_g_per_mm3    = density_g_per_mm3,
    muscle               = muscle,
    lidx_pos_motor       = dbl1(m_attrs[["daq_specimen_lateral_index_on_positive_motor_side"]]),
    lidx_left            = dbl1(m_attrs[["daq_specimen_side_index_left"]]),
    lidx_right           = dbl1(m_attrs[["daq_specimen_side_index_right"]])
  )
}

# Duplicated from run_fv_fl_power_pipeline.R -- see module header note above.
.attach_dynamic_muscle_force <- function(td, torque_col, lidx_pos_motor, lidx_left, lidx_right,
                                          relaxation_s = RELAXATION_WINDOW_S) {
  td$.row_id <- seq_len(nrow(td))
  td2 <- set_cycle_types(td)
  msc <- tryCatch(
    calc_muscle_torque(td2, torque_col = torque_col, include_all_active_samples = TRUE,
                       phase_match_all_rows = TRUE),
    error = function(e) { cli::cli_warn("calc_muscle_torque failed: {conditionMessage(e)}"); NULL }
  )
  td$muscle_force_Nm <- NA_real_
  if (!is.null(msc) && nrow(msc) > 0L && ".row_id" %in% names(msc)) {
    events <- tryCatch(.detect_stim_events(td), error = function(e) tibble::tibble())
    row_side <- rep(NA_character_, nrow(td))
    if (nrow(events) > 0L) {
      events <- dplyr::arrange(events, .data$stim_t0_s)
      for (i in seq_len(nrow(events))) {
        e <- events[i, ]
        in_win <- td$t.s >= e$stim_t0_s & td$t.s <= (e$stim_t1_s + relaxation_s)
        row_side[in_win] <- e$muscle_side
      }
    }
    rec_lidx <- dplyr::case_when(
      row_side == "L" ~ lidx_left,
      row_side == "R" ~ lidx_right,
      .default = NA_real_
    )
    force_sign <- rec_lidx * lidx_pos_motor
    if (!all(!is.finite(force_sign))) {
      msc$muscle_torque.Nm <- force_sign[msc$.row_id] * msc$muscle_torque.Nm
      msc$stim <- row_side[msc$.row_id]
    }
    td$muscle_force_Nm[msc$.row_id] <- msc$muscle_torque.Nm
  }
  msc_active <- if (!is.null(msc) && "cycletype" %in% names(msc)) dplyr::filter(msc, .data$cycletype == "act") else msc
  list(td = td, msc = msc_active)
}


# =============================================================================
# Per-specimen data collection
# =============================================================================

SPECIMENS <- list(
  bass16 = raw_source_dir(BASS16_RAW_SUBFOLDER),
  bass17 = raw_source_dir(BASS17_RAW_SUBFOLDER)
)

collect_specimen <- function(source_dir, label) {
  cli::cli_h1("Collecting: {label}")
  manifest <- parse_trial_directory(source_dir)
  geom0 <- .read_file_level_geometry(manifest$fullpath[1L])
  depth0 <- resolve_muscle_depth_mm(geom0$muscle_depth_mm_raw)
  r_m <- (geom0$local_body_width_mm / 2 - depth0$depth_mm) / 1000
  dclamp_m <- geom0$dclamp_mm / 1000

  load_one <- function(f) {
    td <- load_bender_flat(f, do_filter = TRUE, loadtorques = "x")
    tau <- deconvolve_bender(f, hub_path = NULL, verbose = FALSE)
    N <- min(nrow(td), length(tau)); td <- td[seq_len(N), , drop = FALSE]
    td$torque_inertia_corrected_Nm <- tau[seq_len(N)]
    attr(td, "Filename") <- f
    td
  }

  iso_steps <- purrr::map(manifest$fullpath[manifest$protocol == "isometric"], function(f) {
    analyze_isometric(load_one(f), filename = f)$step_summary
  }) |> dplyr::bind_rows()

  isov_steps <- purrr::map(manifest$fullpath[manifest$protocol == "isovelocity"], function(f) {
    analyze_isovelocity(load_one(f), filename = f)$step_summary
  }) |> dplyr::bind_rows()
  isov_fits <- purrr::map(c(left = "left", right = "right"), function(sd) {
    ss <- dplyr::filter(isov_steps, muscle_side == sd, is.finite(muscle_force_Nm))
    f0_row <- dplyr::filter(ss, abs(shortening_value) < 1e-6)
    f0_iso <- if (nrow(f0_row) > 0L) mean(f0_row$muscle_force_Nm, na.rm = TRUE) else NA_real_
    ss_conc <- dplyr::filter(ss, contraction_mode %in% c("concentric", "isometric_zero"))
    f <- fit_force_velocity_curve(ss_conc$shortening_strain_pct, ss_conc$muscle_force_Nm, side_label = sd, f0_isometric = f0_iso)
    add_specific_properties_to_fit(f, r_m, dclamp_m, geom0$muscle, kind = "FV")
  })

  dyn_cyc <- purrr::map(manifest$fullpath[manifest$protocol == "dynamic"], function(f) {
    td <- load_one(f)
    if (all(as.character(td$stim) == "0")) return(NULL)
    geom_f <- .read_file_level_geometry(f)
    bc <- .attach_dynamic_muscle_force(td, "torque_inertia_corrected_Nm", geom_f$lidx_pos_motor, geom_f$lidx_left, geom_f$lidx_right)
    if (is.null(bc$msc) || nrow(bc$msc) == 0L) return(NULL)
    bc$msc$muscle_mass.kg <- geom_f$muscle$muscle_mass_kg
    summarize_muscle_cycles(bc$msc)
  }) |> dplyr::bind_rows()

  list(label = label, geom = geom0, iso_steps = iso_steps, isov_fits = isov_fits, dyn_cyc = dyn_cyc)
}

results <- purrr::imap(SPECIMENS, ~ collect_specimen(.x, .y))


# =============================================================================
# Console comparison table
# =============================================================================

cli::cli_h1("Cross-specimen comparison")
lims <- coughlin_steady_state_power_limits()
for (r in results) {
  cli::cli_h2(r$label)
  cli::cli_alert_info("Muscle mass {signif(r$geom$muscle$muscle_mass_g,3)} g, CSA {signif(r$geom$muscle$csa_muscle_cm2,3)} cm^2")
  for (sd in c("left","right")) {
    ss <- dplyr::filter(r$iso_steps, muscle_side == sd, is.finite(muscle_force_Nm))
    cli::cli_alert(if (nrow(ss) > 0)
      "Isometric [{sd}]: n={nrow(ss)} points -- no model fit (connect-the-mean default); force range [{signif(min(ss$muscle_force_Nm),3)}, {signif(max(ss$muscle_force_Nm),3)}] N*m"
      else "Isometric [{sd}]: no finite points")
  }
  for (sd in c("left","right")) {
    f <- r$isov_fits[[sd]]
    cli::cli_alert("Isovelocity [{sd}]: {f$status}{if (identical(f$status,'ok')) sprintf(' -- %.3g N/cm^2, %.3g W/kg', f$specific_tension_Ncm2 %||% NA, f$mass_specific_peak_power_Wkg %||% NA) else ''}")
  }
  if (nrow(r$dyn_cyc) > 0) {
    cli::cli_alert("Dynamic: n={nrow(r$dyn_cyc)} cycles -- mean/cycle W/kg [{round(min(r$dyn_cyc$avg_power.Wkg,na.rm=TRUE),1)}, {round(max(r$dyn_cyc$avg_power.Wkg,na.rm=TRUE),1)}]; peak instantaneous W/kg [{round(min(r$dyn_cyc$peak_power.Wkg,na.rm=TRUE),1)}, {round(max(r$dyn_cyc$peak_power.Wkg,na.rm=TRUE),1)}], {sum(r$dyn_cyc$exceeds_coughlin_hi,na.rm=TRUE)}/{nrow(r$dyn_cyc)} ABOVE Coughlin band")
  }
}


# =============================================================================
# Comparison figure
# =============================================================================

wrap_text <- function(s, width) paste(strwrap(s, width = width), collapse = "\n")

side_colors <- c(bass16 = "#1d4ed8", bass17 = "#b91c1c")

dyn_all <- purrr::map_dfr(results, function(r) dplyr::mutate(r$dyn_cyc, specimen = r$label))

p_mean <- ggplot(dyn_all, aes(x = specimen, y = avg_power.Wkg, color = specimen)) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.8) +
  geom_boxplot(width = 0.35, outlier.shape = NA, fill = NA, linewidth = 0.8) +
  scale_color_manual(values = side_colors, guide = "none") +
  labs(title = "Dynamic: mean power per cycle", x = NULL, y = "Mean muscle power (W/kg)") +
  theme_bw(base_size = 11)

p_peak <- ggplot(dyn_all, aes(x = specimen, y = peak_power.Wkg, color = specimen)) +
  annotate("rect", xmin = 0.4, xmax = 2.6, ymin = lims$lo, ymax = lims$hi, fill = "grey40", alpha = 0.15) +
  geom_hline(yintercept = lims$mean, linetype = "dashed", color = "grey40") +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.8) +
  geom_boxplot(width = 0.35, outlier.shape = NA, fill = NA, linewidth = 0.8) +
  scale_color_manual(values = side_colors, guide = "none") +
  labs(title = "Dynamic: peak instantaneous power per cycle",
       subtitle = sprintf("Grey band = Coughlin et al. 1996 (%.0f-%.0f W/kg)", lims$lo, lims$hi),
       x = NULL, y = "Peak muscle power (W/kg)") +
  theme_bw(base_size = 11)

# A fit can report status "ok" (converged) yet be physically degenerate --
# e.g. bass17 left has Vmax ~ 0 (all points effectively at zero velocity),
# which makes F0/specific-tension/power numerically well-defined but not a
# real measurement. Exclude any fit with a near-zero Vmax from the
# comparison (quality gate specific to THIS comparison, not a change to the
# underlying fit_force_velocity_curve()/analyze_isovelocity() status logic).
MIN_MEANINGFUL_VMAX_PCT_PER_S <- 1e-3
isov_df <- purrr::map_dfr(results, function(r) {
  purrr::imap_dfr(r$isov_fits, function(f, sd) {
    if (!identical(f$status, "ok")) return(tibble::tibble())
    if (!is.finite(f$Vmax) || f$Vmax < MIN_MEANINGFUL_VMAX_PCT_PER_S) return(tibble::tibble())
    tibble::tibble(specimen = r$label, side = sd,
                   specific_tension_Ncm2 = f$specific_tension_Ncm2 %||% NA_real_,
                   mass_specific_peak_power_Wkg = f$mass_specific_peak_power_Wkg %||% NA_real_)
  })
})

p_isov_tension <- if (nrow(isov_df) > 0) {
  ggplot(isov_df, aes(x = interaction(specimen, side, sep = " "), y = specific_tension_Ncm2, fill = specimen)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = side_colors, guide = "none") +
    labs(title = "Isovelocity: specific tension", x = NULL, y = "Specific tension (N/cm^2)",
         caption = "Only fits with status \"ok\" shown (bass16 both sides + bass17 left FAILED/degenerate -- omitted)") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
} else NULL

p_isov_power <- if (nrow(isov_df) > 0) {
  ggplot(isov_df, aes(x = interaction(specimen, side, sep = " "), y = mass_specific_peak_power_Wkg, fill = specimen)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = side_colors, guide = "none") +
    labs(title = "Isovelocity: mass-specific peak power", x = NULL, y = "Peak power (W/kg)",
         caption = "Only fits with status \"ok\" shown") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
} else NULL

panels <- list(p_mean, p_peak, p_isov_tension, p_isov_power)
panels <- panels[!sapply(panels, is.null)]

p_combined <- patchwork::wrap_plots(panels, ncol = 2) +
  patchwork::plot_annotation(
    title = "bass16 vs. bass17: mass-/area-specific muscle properties by protocol",
    subtitle = wrap_text(
      "Isometric FL no longer fits a model by default (connect-the-mean only, PI direction 2026-07-16) -- no specific-tension panel for FL. Isovelocity: bass16 both sides FAILED to converge; bass17 left side degenerate (omitted).",
      95
    ),
    caption = wrap_text(sprintf(
      "Muscle mass: bass16 %.3g g, bass17 %.3g g | Muscle CSA: bass16 %.3g cm^2, bass17 %.3g cm^2 | est. from 3%% of an oval body cross-section x test-section length",
      results$bass16$geom$muscle$muscle_mass_g, results$bass17$geom$muscle$muscle_mass_g,
      results$bass16$geom$muscle$csa_muscle_cm2, results$bass17$geom$muscle$csa_muscle_cm2
    ), 100)
  )

out_path <- file.path(OUTPUT_DIR, "specimen_comparison_specific_properties.png")
ggplot2::ggsave(out_path, p_combined, width = 12, height = 9.5, dpi = 150)
cli::cli_alert_success("Saved comparison figure: {out_path}")
