# superplot_fl_pooled_snr_passing.R
# SNR-filtered companion to superplot_fl_pooled.R: keeps ONLY the points/steps
# whose activation_snr >= MFV_UHAT_SNR_MIN ("GOOD"), omitting the rest
# ("BAD" -- indistinguishable from baseline force noise per the SNR survey).
# activation_snr is computed at STEP granularity (one value per isometric
# step, or per isovelocity active ramp) and broadcast to every point that
# step/ramp contributes -- so filtering drops entire steps/ramps, not
# individual samples within a kept step/ramp.
#
# Reuses ALL extraction logic from superplot_fl_pooled.R by sourcing it (this
# also regenerates the existing UNFILTERED diagnostic figure, unchanged) --
# then re-bins/re-plots only the SNR-passing subset of the same `pooled`
# tibble into a second, separate output for bass_summary_figures/.
#
# Run with:  Rscript R/superplot_fl_pooled_snr_passing.R

.root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
source(file.path(.root, "superplot_fl_pooled.R"))  # builds `pooled` (now incl. activation_snr)

SUMMARY_OUTPUT_DIR <- "/Users/yjimenez/Library/CloudStorage/OneDrive-ProvidenceCollege/01_JimenezLab/02_ResearchHub/proj_crittergripper/figures/bass_summary_figures"
fs::dir_create(SUMMARY_OUTPUT_DIR, recurse = TRUE)

n_total <- nrow(pooled)
passing <- dplyr::filter(pooled, is.finite(.data$activation_snr), .data$activation_snr >= SNR_MIN)
n_pass  <- nrow(passing)
cli::cli_h1("SNR-passing pooled FL superplot (activation_snr >= {SNR_MIN})")
cli::cli_alert_info("Kept {n_pass} / {n_total} points ({round(100*n_pass/n_total,1)}%) across {dplyr::n_distinct(passing$fish)} fish, {dplyr::n_distinct(passing$trial_id)} trials")

if (n_pass == 0L) cli::cli_abort("No SNR-passing points -- nothing to plot.")

long_pass <- passing |>
  tidyr::pivot_longer(c(force_emp_N, force_geom_N), names_to = "method_raw", values_to = "force_N") |>
  dplyr::mutate(method = factor(dplyr::if_else(method_raw == "force_emp_N",
                                               method_levels[["empirical"]], method_levels[["geometric"]]),
                                levels = method_levels)) |>
  dplyr::filter(is.finite(force_N), is.finite(shortening_strain_pct))
long_pass$strain_bin <- round(long_pass$shortening_strain_pct / STRAIN_BIN_PCT) * STRAIN_BIN_PCT

per_trial_p <- long_pass |>
  dplyr::group_by(method, fish, trial_id, strain_bin) |>
  dplyr::summarise(force_N = mean(force_N, na.rm = TRUE), .groups = "drop")
per_fish_p <- per_trial_p |>
  dplyr::group_by(method, fish, strain_bin) |>
  dplyr::summarise(force_N = mean(force_N, na.rm = TRUE), .groups = "drop")
per_group_p <- per_fish_p |>
  dplyr::group_by(method, strain_bin) |>
  dplyr::summarise(force_N = mean(force_N, na.rm = TRUE), n_fish = dplyr::n(), .groups = "drop")

n_iso_p <- dplyr::n_distinct(dplyr::filter(passing, protocol == "isometric")$trial_id)
n_isv_p <- dplyr::n_distinct(dplyr::filter(passing, protocol == "isovelocity")$trial_id)

p2 <- ggplot(mapping = aes(x = strain_bin, y = force_N)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_line(data = per_trial_p, aes(group = trial_id, colour = fish), linewidth = 0.3, alpha = 0.35) +
  geom_point(data = per_trial_p, aes(colour = fish), size = 0.8, alpha = 0.4) +
  geom_line(data = per_fish_p, aes(group = fish, colour = fish), linewidth = 1.0, alpha = 0.95) +
  geom_point(data = per_fish_p, aes(colour = fish), size = 2.1, alpha = 0.95) +
  geom_line(data = per_group_p, aes(group = 1), colour = "black", linewidth = 1.8) +
  geom_point(data = per_group_p, colour = "black", size = 2.4) +
  facet_wrap(~method, ncol = 2) +
  scale_colour_manual(values = fish_cols, name = "Individual (fish)") +
  labs(
    title = "Pooled Force-Length superplot -- SNR-PASSING ONLY (activation_snr >= 3.0)",
    subtitle = sprintf(
      "Same pooling as the unfiltered diagnostic superplot, but %s (%d/%d, %.1f%%) omitted for activation_snr < %g\n%d isometric + %d isovelocity trial(s) contribute at least one passing step/ramp | %g%% length bins | connect-the-dots, NO fit",
      "low-confidence steps/ramps", n_total - n_pass, n_total, 100 * (n_total - n_pass) / n_total, SNR_MIN,
      n_iso_p, n_isv_p, STRAIN_BIN_PCT),
    x = "Muscle shortening strain (%, length; + = shortened)",
    y = "Muscle force along u_hat (N, + = reinforces passive direction)") +
  theme_bw(12) +
  theme(legend.position = "right",
        plot.subtitle = element_text(size = 8),
        strip.text = element_text(face = "bold"))

outfile2 <- file.path(SUMMARY_OUTPUT_DIR, "FLsuperplot_isometric_isovelocity_pooled_snrPass.png")
ggplot2::ggsave(outfile2, p2, width = 13, height = 6.5, dpi = 150)
cli::cli_alert_success("Wrote {outfile2}")
