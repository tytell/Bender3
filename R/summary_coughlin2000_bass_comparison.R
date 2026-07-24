# summary_coughlin2000_bass_comparison.R
# PI request 2026-07-24: compare OUR measured mass-specific power (W/kg),
# work per cycle (J/kg) and specific tension (kN/m^2) against published
# largemouth bass RED MUSCLE values from Coughlin (2000, J. Exp. Biol.
# 203:617-629) -- a DIFFERENT reference than the "Coughlin & Carroll (2006)"
# activation/relaxation-time band used in isometric_L0_activation_kinetics*.
# Style follows isometric_L0_activation_kinetics_bookendsOnly.png (twitch
# trace on top, boxplots below), per PI direction.
#
# Because the user's clarifying-question form was skipped, the RECOMMENDED
# defaults from that 10-question set were used (see chat for the full list):
#   1. ONE combined figure (not three).
#   2. Panel A = the SAME isometric/bookend twitch trace as before (twitch
#      data pairs with the tension metric; power/work come from a different,
#      dynamic protocol -- see caveats below), refined per the "except"
#      requests (tighter x-axis; explicit activation/relaxation windows).
#   3. Power & work computed from DYNAMIC trials only (true work-loop match
#      to Coughlin's method), right-muscle-stim cycles only (stim=="R"),
#      matching this repo's "active = right stim" convention throughout.
#   4. LATER (stable) trials only (dynamic_trial_precondition.R cutoffs) --
#      avoids the known early-trial sono/geometric power-inflation artifact.
#   5. The work/power benchmark is a SINGLE literature data point (no
#      position x speed range was extractable as text -- see below) and is
#      shown as such, explicitly labeled.
#   6. The huge tension gap (ours ~15-50x below Coughlin's) is annotated
#      with the most likely explanation (CSA method mismatch).
#   7. Twitch:tetanus ratio (0.279+/-0.030) NOT added (not requested).
#   8. Temperature caveat added: rig has NO logged bath/room temperature
#      anywhere in HDF5 metadata; assumed ambient ~20-22 C, uncontrolled.
#   9. Activation/relaxation windows shown as shaded bands directly on the
#      twitch trace (not a separate inset).
#  10. Existing isometric_L0_activation_kinetics*.png figures NOT retrofitted
#      -- this is a new, separate figure.
#
# =============================================================================
# Coughlin (2000) largemouth bass, 20 C bath -- values extractable as TEXT
# from the PDF (i.e., NOT read off a bar/scatter figure by eye):
#   - Tetanic isometric tension: 186.4 +/- 33.6 kN/m^2 (mean+/-S.E.M., N=18
#     bundles), POOLED across all 6 longitudinal positions -- position had
#     NO significant effect on tension (F=0.151, P=0.861), so there is no
#     single "50%L" tension value to target; the pooled value is the best
#     available number regardless of exact position.
#   - Work per oscillation cycle: 3.56 +/- 0.47 J/kg (mean+/-S.E.M., N=6
#     bundles) -- ONE data point, embedded as on-figure text in Fig. 6, for
#     a muscle bundle at 0.572L ("MID", closest available position to the
#     user's ~50%L; MID bundles were dissected from 0.55-0.70L) stimulated
#     under in vivo conditions for swimming at 2.4 L/s (4.05 Hz tailbeat --
#     the FASTEST bass swimming speed tested, not a full speed range).
#   - Power is NOT given directly for this point; it is DERIVED here using
#     Coughlin's own definition ("Power is the product of the net work per
#     cycle and the oscillation frequency"): 3.56 J/kg * 4.05 Hz = 14.4 W/kg.
#     The SD on power is likewise propagated (0.47 * 4.05 = 1.9 W/kg) --
#     this is an approximation, not an independently reported error.
#   - All other position x speed combinations for bass work/power (the full
#     grids in Figs 6, 7, 9) exist ONLY as bar/scatter graphics in the PDF;
#     PDF text extraction did not recover their values, so they are NOT
#     used here. If a fuller position x speed grid is needed, it would have
#     to be digitized from the original figures.
#
# CAVEATS shown on the figure:
#   - Apples-to-oranges risk: Coughlin's numbers are muscle stimulated under
#     SUBMAXIMAL, in-vivo SWIMMING conditions (not the animal's maximum,
#     W_max); our dynamic/isometric protocols are closer to maximal
#     characterization stimulation. This can push OUR power/tension values
#     either above or below the swimming-condition benchmark for reasons
#     unrelated to data quality.
#   - Tension gap: our isometric specific tension is ~15-50x BELOW
#     Coughlin's value. Likely explanation: our specific tension uses a
#     GEOMETRIC whole-body-oval CSA estimate (compute_muscle_mass_and_csa(),
#     muscle_geometry.R), while Coughlin's is a HISTOLOGICAL live-fiber-area
#     measurement on an isolated single-myotome bundle (Trypan
#     Blue/SDH-stained cross-sections) -- the geometric estimate very likely
#     overestimates true contractile CSA (includes non-muscle tissue,
#     connective tissue, and any submaximal recruitment in situ).
#   - No bath/room temperature is logged anywhere in this rig's HDF5
#     metadata (see .cursorrules: no such field exists in bender_h5_export.py
#     or bender_functions.py). Assumed ambient ~20-22 C, uncontrolled,
#     vs. Coughlin's temperature-controlled 20 C bath.
#
# Run with:  Rscript R/summary_coughlin2000_bass_comparison.R
# Outputs -> figures: 02_processed/figs_summary/ (FIGS_SUMMARY_DIR)
#            data:    02_processed/data_processed/

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(ggplot2); library(patchwork); library(cli); library(readr)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")
src("dynamic_trial_precondition.R")

OUT_DIR      <- FIGS_SUMMARY_DIR
DATA_OUT_DIR <- file.path(.crittergripper_root(), "02_processed", "data_processed")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SPECIMEN_COLORS <- c(bass16 = "#1b9e77", bass17 = "#d95f02", bass18 = "#7570b3")
PAD_S <- 0.2

# =============================================================================
# Coughlin (2000) bass reference values (see header for provenance/caveats)
# =============================================================================
COUGHLIN2000 <- list(
  tension_kNm2 = list(mean = 186.4, sd = 33.6, n = 18,
                       note = "tetanic force, pooled across all 6 positions (no position effect)"),
  work_Jkg     = list(mean = 3.56, sd = 0.47, n = 6,
                       note = "MID bundle, 0.572L, swimming at 2.4 L/s (4.05 Hz) -- single point, not a range"),
  power_Wkg    = list(mean = 3.56 * 4.05, sd = 0.47 * 4.05, n = 6,
                       note = "DERIVED (work x frequency); SD propagated, not independently reported")
)

# =============================================================================
# PANEL A: isometric bookend twitch trace (reuses collected traces/times from
# summary_isometric_l0_activation.R -- no raw H5 reprocessing).
# =============================================================================
traces_all <- readr::read_csv(file.path(DATA_OUT_DIR, "isometric_l0_activation_traces.csv"), show_col_types = FALSE)
times_all  <- readr::read_csv(file.path(DATA_OUT_DIR, "isometric_l0_activation_times.csv"), show_col_types = FALSE)

traces <- dplyr::filter(traces_all, .data$source == "dynamic L0 bookend")
times  <- dplyr::filter(times_all,  .data$source == "dynamic L0 bookend")

# representative activation/relaxation windows (pooled median, ms -> s) --
# same TA/TR definitions as summary_isometric_l0_activation.R:
#   activation = t(stim onset=0) -> t(rise to 90% of peak)
#   relaxation = t(stim offset)  -> t(fall to 50% of peak)
ta_med_s  <- stats::median(times$activation_ms, na.rm = TRUE) / 1000
tr_med_s  <- stats::median(times$relaxation_ms, na.rm = TRUE) / 1000
dur_med_s <- stats::median(times$stim_duration_s, na.rm = TRUE)
cli::cli_alert_info("Bookend medians: stim duration {round(dur_med_s*1000)} ms, activation {round(ta_med_s*1000)} ms, relaxation {round(tr_med_s*1000)} ms (n={nrow(times)})")

XLIM <- c(-PAD_S, 0.45)  # tightened -- bookends are ~54 ms twitches, fully decayed by ~0.4 s

grid <- seq(XLIM[1], XLIM[2], by = 0.002)
mean_trend <- purrr::map_dfr(split(traces, traces$specimen), function(df) {
  units_here <- dplyr::n_distinct(df$unit_id)
  m <- purrr::map(split(df, df$unit_id), function(u)
    stats::approx(u$t_rel, u$force_norm, xout = grid, rule = 1)$y)
  M <- do.call(cbind, m)
  tibble::tibble(specimen = df$specimen[1L], t_rel = grid,
                 force_norm = rowMeans(M, na.rm = TRUE),
                 n = rowSums(is.finite(M)),
                 min_support = pmax(10, 0.50 * units_here))
}) |> dplyr::filter(.data$n >= .data$min_support)

pTop <- ggplot() +
  annotate("rect", xmin = 0, xmax = dur_med_s, ymin = -Inf, ymax = Inf, fill = "grey70", alpha = 0.25) +
  annotate("rect", xmin = 0, xmax = ta_med_s, ymin = -Inf, ymax = Inf, fill = "#1b9e77", alpha = 0.15) +
  annotate("rect", xmin = dur_med_s, xmax = dur_med_s + tr_med_s, ymin = -Inf, ymax = Inf, fill = "#d95f02", alpha = 0.15) +
  geom_hline(yintercept = c(0, 0.9, 0.5, 1), linetype = "dotted", color = "grey80") +
  geom_vline(xintercept = c(0, dur_med_s), linetype = "dashed", color = "grey50") +
  geom_line(data = traces, aes(x = t_rel, y = force_norm, group = unit_id, color = specimen),
            alpha = 0.20, linewidth = 0.3) +
  geom_line(data = mean_trend, aes(x = t_rel, y = force_norm, color = specimen), linewidth = 1.4) +
  annotate("text", x = ta_med_s / 2, y = 1.08, label = sprintf("activation\n(median %.0f ms)", ta_med_s * 1000),
           size = 2.6, color = "#1b9e77", fontface = "italic") +
  annotate("text", x = dur_med_s + tr_med_s / 2, y = 1.08, label = sprintf("relaxation\n(median %.0f ms)", tr_med_s * 1000),
           size = 2.6, color = "#d95f02", fontface = "italic") +
  annotate("text", x = dur_med_s / 2, y = -0.13, label = "stim ON", size = 2.4, color = "grey30") +
  scale_color_manual(values = SPECIMEN_COLORS, name = "individual") +
  coord_cartesian(xlim = XLIM, ylim = c(-0.18, 1.15)) +
  labs(title = "A. Isometric L0 bookend twitches: normalised force vs. time",
       subtitle = sprintf("Same twitches as isometric_L0_activation_kinetics_bookendsOnly.png, x-axis tightened to the twitch itself.\nShaded: stim ON (grey), activation window (green, onset->90%% peak), relaxation window (orange, offset->50%% decay). n=%d twitches.", dplyr::n_distinct(traces$unit_id)),
       x = "Time relative to stim onset (s)", y = "Normalised force (F / peak)") +
  theme_bw(base_size = 12)

# =============================================================================
# PANEL B/C/D data: LATER (stable) trials only, right-muscle-stim cycles
# (stim=="R") for dynamic power/work, right-muscle isometric steps for
# tension (already restricted upstream).
# =============================================================================
tension <- readr::read_csv(file.path(DATA_OUT_DIR, "isometric_tension_vs_offset_by_trial.csv"), show_col_types = FALSE) |>
  dplyr::mutate(precondition = classify_session_precondition(.data$specimen, .data$trial_num)) |>
  dplyr::filter(.data$precondition == "later (stable)") |>
  dplyr::transmute(.data$specimen, .data$trial_id, metric = "Specific tension (kN/m^2)",
                    value = .data$max_tension_Ncm2 * 10)  # N/cm^2 -> kN/m^2 (x10)

percycle <- readr::read_csv(file.path(DATA_OUT_DIR, "dynamic_precondition_power_percycle.csv"), show_col_types = FALSE) |>
  dplyr::filter(.data$precondition == "later (stable)", .data$stim == "R")

power <- percycle |>
  dplyr::group_by(.data$specimen, .data$trial_id) |>
  dplyr::summarise(value = mean(.data$avg_power.Wkg, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(metric = "Mean cycle power (W/kg)")

work <- percycle |>
  dplyr::group_by(.data$specimen, .data$trial_id) |>
  dplyr::summarise(value = mean(.data$work.Jkg, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(metric = "Mean work per cycle (J/kg)")

pooled <- dplyr::bind_rows(tension, power, work) |>
  dplyr::filter(is.finite(.data$value)) |>
  dplyr::mutate(metric = factor(.data$metric, levels = c("Specific tension (kN/m^2)", "Mean cycle power (W/kg)", "Mean work per cycle (J/kg)")))

n_tbl <- pooled |> dplyr::count(.data$metric, name = "n_trials")
cli::cli_h2("n trials (later/stable) per metric")
print(n_tbl)
write.csv(pooled, file.path(DATA_OUT_DIR, "coughlin2000_bass_comparison_data.csv"), row.names = FALSE)

.coug_panel <- function(m, coug_lo, coug_hi, coug_mid, coug_label, ylab, title) {
  d <- dplyr::filter(pooled, .data$metric == m)
  yr <- range(c(d$value, coug_lo, coug_hi), na.rm = TRUE)
  ytop <- coug_hi + 0.12 * diff(yr)
  ggplot(d, aes(x = 1, y = .data$value)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = coug_lo, ymax = coug_hi, fill = "#b30000", alpha = 0.12) +
    geom_hline(yintercept = coug_mid, color = "#b30000", linetype = "dashed", linewidth = 0.5) +
    annotate("text", x = 1, y = ytop, label = coug_label, hjust = 0.5, vjust = 1, size = 2.5, color = "#b30000", lineheight = 0.9) +
    geom_boxplot(width = 0.5, outlier.shape = NA, fill = NA, color = "grey40") +
    geom_jitter(aes(color = .data$specimen), width = 0.15, size = 2.4, alpha = 0.85) +
    scale_color_manual(values = SPECIMEN_COLORS, name = "individual", guide = "none") +
    scale_x_continuous(breaks = NULL, limits = c(0.4, 1.6)) +
    coord_cartesian(ylim = c(min(yr[1], 0), ytop + 0.08 * diff(yr)), clip = "off") +
    labs(title = title, x = NULL, y = ylab) +
    theme_bw(base_size = 12) + theme(plot.margin = margin(t = 14, r = 6, b = 5.5, l = 5.5))
}

t_c <- COUGHLIN2000$tension_kNm2
pB <- .coug_panel("Specific tension (kN/m^2)",
                   t_c$mean - t_c$sd, t_c$mean + t_c$sd, t_c$mean,
                   sprintf("Coughlin 2000 bass, 20C\n(n=%d bundles, all positions)", t_c$n),
                   "Specific tension (kN/m^2)", "B. Isometric tension (trial-max)")

p_c <- COUGHLIN2000$power_Wkg
pC <- .coug_panel("Mean cycle power (W/kg)",
                   p_c$mean - p_c$sd, p_c$mean + p_c$sd, p_c$mean,
                   "Coughlin 2000 bass, 20C\n(derived, 1 pt: 0.572L @ 2.4 L/s)",
                   "Mean cycle power (W/kg)", "C. Dynamic power (right-stim cycles)")

w_c <- COUGHLIN2000$work_Jkg
pD <- .coug_panel("Mean work per cycle (J/kg)",
                   w_c$mean - w_c$sd, w_c$mean + w_c$sd, w_c$mean,
                   sprintf("Coughlin 2000 bass, 20C\n(n=%d bundles, 1 pt: 0.572L @ 2.4 L/s)", w_c$n),
                   "Work per cycle (J/kg)", "D. Dynamic work (right-stim cycles)")

fig <- pTop / (pB | pC | pD) +
  patchwork::plot_layout(heights = c(1.15, 1)) +
  patchwork::plot_annotation(
    title = "Bass red muscle vs. Coughlin (2000): tension, power, work -- later/stable trials only",
    subtitle = paste0(
      "Later (stable) trials only (dynamic_trial_precondition.R cutoffs); dynamic power/work use right-stim cycles only.\n",
      "CAVEATS: (1) Coughlin's tension is position-pooled (no 50%L-specific value; position had no effect); work/power is a SINGLE data point at\n",
      "0.572L, 2.4 L/s (fastest speed tested), not a full range -- power is DERIVED (work x freq), not directly reported.\n",
      "(2) Coughlin's values are muscle under SUBMAXIMAL in-vivo swimming conditions; ours are closer to maximal characterization stimulation -- not a like-for-like comparison.\n",
      "(3) Our tension uses MEASURED_RED_MUSCLE_CSA_CM2 (image-analysis CSA, reference specimen 'bass07', assumed ~50%L; muscle_geometry.R), not bass16/17/18's own CSA (never imaged) -- still ~11-20x below Coughlin's, so the gap is NOT fully explained by the earlier geometric-CSA guess alone.\n",
      "(4) No bath/room temperature is logged on this rig; assumed ambient ~20-22C (uncontrolled) vs. Coughlin's controlled 20C bath."),
    theme = theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 8.3)))

fout <- file.path(OUT_DIR, "coughlin2000_bass_power_work_tension_comparison.png")
ggplot2::ggsave(fout, fig, width = 13, height = 9.5, dpi = 150)
cli::cli_alert_success("Saved {fout}")
