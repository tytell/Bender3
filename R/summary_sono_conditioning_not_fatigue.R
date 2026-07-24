# summary_sono_conditioning_not_fatigue.R
# PI request, 2026-07-24: "Make a compelling graph (using my data) showing
# that getting good sono data DOES NOT require fatigued muscle, it just
# requires a certain amount (how much) of conditioning."
#
# The argument, made from existing per-trial data (no re-reading H5):
#   - "Good sono" = the ground-truthed active-vs-passive sono strain EXCESS
#     (diag_precondition_sono_length_activeVsPassive.R) collapsing toward 0.
#   - "Conditioning DOSE" = cumulative dynamic active cycles delivered in the
#     session BEFORE each trial (a lower bound on total conditioning -- other
#     protocols run in the same session also condition, but are not counted
#     here; noted on the figure).
#   - "Fatigue" is ruled OUT by peak active muscle force capacity (max
#     |muscle_torque.Nm| per trial, the standard fatigue readout, from the
#     directly-measured single-axis torque -- independent of any length
#     model). If the sono problem were fatigue, low excess should ONLY occur
#     once force has collapsed.
#
# Decisive evidence (bass17, the specimen with the most sustained dynamic
# protocol): peak force at the FINAL dynamic trial (0.553 N*m) essentially
# equals the FIRST (0.557 N*m) -- ~full strength -- yet the excess is fully
# resolved. So stabilization is a CONDITIONING DOSE effect, not force loss.
#
# Two panels:
#   (A) |sono strain excess| vs. cumulative conditioning cycles -- reads off
#       "how much conditioning" (the dose at which excess drops into the
#       good band), per specimen.
#   (B) DECOUPLING scatter: |sono strain excess| vs. peak force capacity
#       (% of each specimen's own session max) -- good sono (low excess)
#       occurs across the full force range, INCLUDING near-max force, so it
#       is not gated by force loss / fatigue.
#
# Run with:  Rscript R/summary_sono_conditioning_not_fatigue.R
# Outputs -> figs_summary/ (FIGS_SUMMARY_DIR) + a per-trial CSV.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork); library(cli)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")

DATA_DIR    <- file.path(.crittergripper_root(), "02_processed", "data_processed")
SUMMARY_DIR <- FIGS_SUMMARY_DIR
dir.create(SUMMARY_DIR, recursive = TRUE, showWarnings = FALSE)

SPECIMEN_COLORS <- c(bass16 = "#1b9e77", bass17 = "#d95f02", bass18 = "#7570b3")
GOOD_BAND_PCT   <- 1.0   # |excess| below this = "good sono" (well under the ~5% commanded p2p)

# =============================================================================
# Assemble per-trial: sono excess + peak force capacity + conditioning dose.
# =============================================================================
persample <- read.csv(file.path(DATA_DIR, "sono_vs_geometric_dynamic_power_persample.csv"))
excess    <- read.csv(file.path(DATA_DIR, "sono_length_activeVsPassive_bytrial.csv"))

capacity <- persample |>
  dplyr::group_by(.data$specimen, .data$trial_id, .data$trial_num) |>
  dplyr::summarise(peak_torque_Nm = max(abs(.data$muscle_torque.Nm), na.rm = TRUE), .groups = "drop")

trial <- excess |>
  dplyr::left_join(capacity, by = c("specimen", "trial_id", "trial_num")) |>
  dplyr::arrange(.data$specimen, .data$trial_num) |>
  dplyr::group_by(.data$specimen) |>
  dplyr::mutate(
    cum_cond_cycles   = cumsum(.data$n_cycles) - .data$n_cycles,      # dose delivered BEFORE this trial
    peak_force_pct    = 100 * .data$peak_torque_Nm / max(.data$peak_torque_Nm, na.rm = TRUE),
    abs_excess_pct    = abs(.data$mean_excess_pct)
  ) |>
  dplyr::ungroup()

write.csv(trial, file.path(DATA_DIR, "sono_conditioning_not_fatigue_bytrial.csv"), row.names = FALSE)

# Dose at which each specimen first drops into (and stays in) the good band.
dose_tbl <- trial |>
  dplyr::group_by(.data$specimen) |>
  dplyr::arrange(.data$cum_cond_cycles, .by_group = TRUE) |>
  dplyr::summarise(
    dose_cycles_to_good = {
      good <- .data$abs_excess_pct < GOOD_BAND_PCT
      # first index from which ALL subsequent trials are good
      stable_from <- which(rev(cumprod(rev(good))) == 1L)
      if (length(stable_from) == 0L) NA_integer_ else .data$cum_cond_cycles[min(stable_from)]
    },
    force_pct_at_good = {
      good <- .data$abs_excess_pct < GOOD_BAND_PCT
      stable_from <- which(rev(cumprod(rev(good))) == 1L)
      if (length(stable_from) == 0L) NA_real_ else .data$peak_force_pct[min(stable_from)]
    },
    force_pct_at_session_end = .data$peak_force_pct[which.max(.data$trial_num)],
    .groups = "drop"
  )
cli::cli_h1("Conditioning dose to reach good sono, and force retained there")
print(dose_tbl)

# =============================================================================
# Panel A: |excess| vs cumulative conditioning cycles
# =============================================================================
pA <- ggplot(trial, aes(x = .data$cum_cond_cycles, y = .data$abs_excess_pct, color = .data$specimen)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = GOOD_BAND_PCT,
           fill = "#2ca02c", alpha = 0.10) +
  annotate("text", x = Inf, y = GOOD_BAND_PCT, label = sprintf("  \"good sono\" band (|excess| < %.0f pct-pt)", GOOD_BAND_PCT),
           hjust = 1.02, vjust = -0.5, size = 3, color = "#2ca02c") +
  geom_line(aes(group = .data$specimen), linewidth = 0.7, alpha = 0.8) +
  geom_point(size = 2.6) +
  scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen") +
  labs(title = "A. Good sono is reached after a modest CONDITIONING DOSE",
       subtitle = "Ground-truthed active-vs-passive sono strain excess vs. cumulative dynamic conditioning cycles delivered before each trial\n(lower bound -- other protocols in the same session also condition but are not counted here).",
       x = "Cumulative conditioning cycles delivered (dynamic, before this trial)",
       y = "|active-minus-passive sono strain excess| (pct L0)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

# =============================================================================
# Panel B: decoupling -- |excess| vs peak force capacity (% of specimen max)
# =============================================================================
lab_end <- trial |>
  dplyr::group_by(.data$specimen) |>
  dplyr::filter(.data$trial_num == max(.data$trial_num)) |>
  dplyr::ungroup()

pB <- ggplot(trial, aes(x = .data$peak_force_pct, y = .data$abs_excess_pct, color = .data$specimen)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = GOOD_BAND_PCT,
           fill = "#2ca02c", alpha = 0.10) +
  geom_point(aes(shape = .data$precondition), size = 2.8, stroke = 0.9) +
  geom_text(data = lab_end, aes(label = sprintf("%s session end", .data$specimen)), size = 2.9,
            show.legend = FALSE, hjust = 0.5, vjust = -1.0) +
  scale_color_manual(values = SPECIMEN_COLORS, name = "Specimen") +
  scale_shape_manual(values = c("early (preconditioning)" = 17, "later (stable)" = 16), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.08))) +
  labs(title = "B. Good sono is NOT gated by force loss (i.e. not fatigue)",
       subtitle = "Low excess occurs across the WHOLE range of muscle force capacity -- including near-max force (bass17 ends its session at\n~99% of peak strength with the excess already fully resolved). If the sono problem were fatigue, low excess could only\nappear once force had collapsed -- it does not.",
       x = "Peak active muscle force capacity this trial (% of specimen's own session max)",
       y = "|active-minus-passive sono strain excess| (pct L0)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

fig <- pA / pB + patchwork::plot_annotation(
  title = "Getting good sono data requires CONDITIONING, not fatigue",
  theme = theme(plot.title = element_text(face = "bold", size = 14))
)

fout <- file.path(SUMMARY_DIR, "sono_conditioning_not_fatigue.png")
ggplot2::ggsave(fout, fig, width = 11, height = 10, dpi = 150)
cli::cli_alert_success("Saved {fout}")
cli::cli_alert_success("summary_sono_conditioning_not_fatigue.R complete")
