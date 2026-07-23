# diag_isovelocity_hillcheck.R
# READ-ONLY DIAGNOSTIC (2026-07-22, PI-requested: "I'm not seeing how bass18
# resembles the Hill-type relationship... please produce the plot"). Sources
# superplot_fv_pooled.R and only re-aggregates its already-built `pooled` frame
# -- computes no new force, modifies no production analysis.
#
# QUESTION: does the pointwise-corrected isovelocity FV (flag-1 fix,
# analysis_muscle_force_vector_log.md "flag 1 ROOT-CAUSED + FIXED") actually
# reproduce the Hill-type relationship (monotonic-DECREASING force with
# increasing shortening/concentric velocity; eccentric >= isometric >=
# concentric), NOT a bell/symmetric-peak-at-V=0 (that's the FL target, not
# FV's -- see the same log entry's "TARGET SHAPE CORRECTION")?
#
# ONE PANEL, 3 fish facets: every point plotted, alpha-flagged by its 4-tier
# confidence (confident/confidently_small/unstable_magnitude/unconfirmable --
# see mfv_confidence_tier(), muscle_force_vector.R), never dropped per
# cursorrule; a per-(fish, side) SUMMARY LINE connects the MEAN force at each
# unique signed strain-rate value, computed from confident + confidently_small
# points only (the line is a derived summary, not raw data, so it is
# legitimate to restrict it to the trustworthy subset while still showing
# every raw point). A dashed reference arrow/annotation marks the Hill-expected
# direction (force falls as concentric velocity rises) for visual comparison.
#
# REVISED 2026-07-22 (PI-directed "SNR-based confidence gating audit" --
# implements analysis_muscle_force_vector_log.md's 2026-07-22 proposal): the
# `confident` filter feeding the summary line below used to be RATIO-only
# (`activation_snr >= snr_min`), which is exactly what produced the V=0 notch
# described below -- bass18's stronger trial (bender_03)'s V=0 reps (1.99/
# 1.31 N, consistent with its own eccentric plateau) got excluded from the
# mean purely because that trial's OWN noise floor happens to be higher, not
# because the force itself is unreal (confirmed directly:
# R/diag_snr_magnitude_audit.R shows those exact reps are "confidently_small"
# -- SNR fails, but |force| clears its own baseline_force_noise_N easily).
# The summary line now includes "confidently_small" points (real, just
# below-ratio) alongside "confident" ones -- ONLY "unstable_magnitude"/
# "unconfirmable" points are excluded from the mean now. Point-level alpha
# uses the full 4-tier factor (was a 2-level SNR-only TRUE/FALSE).
#
# FINDING (revised after actually building and looking at this plot -- the
# PI's skepticism was justified; an earlier verbal claim that bass18 "reproduces
# the correct Hill-type ordering" was based on comparing concentric-vs-eccentric
# PAIRS at each |v| and ignoring the V=0 isometric anchor entirely. Once V=0 is
# included, the full curve is NOT a clean monotonic decline):
# bass18 right, ORIGINAL ratio-only gate: eccentric plateau ~1.65-1.76 N (-127
# to -382 %/s) -> DROPS SHARPLY to 0.56 N at V=0 (n=3, ALL from the SAME
# weaker trial-set, bender_06/07/08/09/10) -> back UP to 1.37 N (127) -> DOWN
# to 1.03 N (255) -> back UP to 1.83 N (382). A "W", not a Hill hyperbola.
# RE-RUN 2026-07-22 under the magnitude-aware confidence tier (see REVISED
# note above): the V=0 mean is now 0.94 N (n=7, up from n=3) -- the stronger
# trial's (bender_03) V=0 reps are "confidently_small" (SNR 2.08/2.55 < 3.0,
# but |force| clears baseline_force_noise_N easily) and are now correctly
# included instead of silently dropped. This CONFIRMS those reps were real,
# not noise -- but the V=0 notch is only PARTIALLY resolved (0.56 -> 0.94 N),
# not eliminated: 0.94 N is still below the eccentric plateau (1.65-2.33 N)
# and most concentric points (1.37-3.63 N, though its own low end, 1.03 N at
# 255 %/s, is close). So the magnitude-aware fix rules out "pure ratio-gating
# artifact" as the FULL explanation, but a genuine partial dip near V=0 (or a
# residual cross-trial-fatigue confound -- see original root-cause note below)
# remains and is NOT settled by this fix alone. bass16/17 show a PEAK AT V=0
# with decline on BOTH sides (a tent, not Hill's plateau-then-monotonic-
# decline, and not flat noise either) -- unexplained, flagged.
# CONCLUSION: the moving-ramp-only pairwise comparison (eccentric > concentric
# at 127/255 %/s) remains directionally Hill-consistent for bass18, but the
# FULL curve including V=0 is NOT a clean monotonic decline even under the
# magnitude-aware gate -- this needs a within-trial-only test (each trial's
# OWN V=0 + moving steps, not pooled across trials) before it can be called
# settled.
#
# Run: Rscript R/diag_isovelocity_hillcheck.R
# Canon token: isovhillcheck -- see FIGURES_README.md.

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(cli)
})

.root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
source(file.path(.root, "superplot_fl_pooled.R")); pooled_fl <- pooled
source(file.path(.root, "superplot_fv_pooled.R")); pooled_fv <- pooled
OUT <- FIGS_DIAGNOSTIC_DIR

sides_pal <- c(left = "#2563eb", right = "#ea580c")
snr_min <- MFV_UHAT_SNR_MIN

# isovelocity's own V=0 hold + moving ramps -- protocol-consistent single
# source (NOT mixing in isometric L0 / dynamic bookends, so this is purely
# "does isovelocity's own passive-corrected force respect Hill ordering").
d <- pooled_fv |>
  dplyr::filter(.data$protocol == "isovelocity", is.finite(.data$force_geom_N),
                is.finite(.data$strain_rate_pct_s), is.finite(.data$activation_snr),
                is.finite(.data$baseline_force_noise_N)) |>
  dplyr::mutate(confidence_tier = mfv_confidence_tier(.data$force_geom_N, .data$activation_snr,
                                                       .data$baseline_force_noise_N, snr_min))

means <- d |>
  dplyr::filter(.data$confidence_tier %in% c("confident", "confidently_small")) |>
  dplyr::group_by(.data$fish, .data$muscle_side, .data$strain_rate_pct_s) |>
  dplyr::summarise(force_geom_N = mean(.data$force_geom_N, na.rm = TRUE), n = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(.data$fish, .data$muscle_side, .data$strain_rate_pct_s)

p <- ggplot() +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_vline(xintercept = 0, color = "grey80", linetype = "dashed") +
  geom_point(data = d, aes(strain_rate_pct_s, force_geom_N, color = muscle_side, alpha = confidence_tier), size = 1.3) +
  geom_line(data = means, aes(strain_rate_pct_s, force_geom_N, color = muscle_side), linewidth = 1.0) +
  geom_point(data = means, aes(strain_rate_pct_s, force_geom_N, color = muscle_side), size = 2.2, shape = 21, fill = "white") +
  scale_color_manual(values = sides_pal, name = "muscle side") +
  scale_alpha_manual(values = MFV_CONFIDENCE_ALPHA, name = "confidence tier (SNR x magnitude)", drop = FALSE) +
  facet_wrap(~fish, scales = "free_y") +
  labs(title = "Isovelocity Hill check: is force-velocity monotonic-DECREASING (eccentric >= isometric >= concentric), not a bell?",
       subtitle = "Solid line = mean of confident+confidently_small points at each velocity (isovelocity's own V=0 hold + moving ramps only, pointwise angle-matched passive) -- REVISED 2026-07-22 to be magnitude-aware, see R/diag_snr_magnitude_audit.R + analysis log; only unstable_magnitude/unconfirmable points are excluded from the mean now. All 4 tiers shown, never dropped (alpha only). Hill target: HIGH on the left (eccentric/lengthening), falling toward LOW on the right (concentric/high shortening velocity) -- a bell (high at V=0, low at both extremes) would be the WRONG shape here (that's FL's target).",
       x = "Muscle shortening strain rate (%/s; + = concentric/shortening, - = eccentric/lengthening, 0 = isometric V=0 anchor)",
       y = "Muscle force (N, geometric u_hat)") +
  theme_bw(10) + theme(plot.subtitle = element_text(size = 7.5), legend.position = "bottom")
ggsave(file.path(OUT, "isovhillcheck.png"), p, width = 12, height = 5, dpi = 140)

cli::cli_h2("Isovelocity Hill check -- confidence-tier-passing (confident+confidently_small) summary means by fish/side/velocity")
print(as.data.frame(means), digits = 3)
cli::cli_alert_success("Wrote isovhillcheck.png to {OUT}")
