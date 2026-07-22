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
# ONE PANEL, 3 fish facets: every SNR-passing (activation_snr >= MFV_UHAT_SNR_MIN)
# point plotted at full opacity; SNR-failing points shown faint (alpha-flagged,
# never dropped, per cursorrule); a per-(fish, side) SUMMARY LINE connects the
# MEAN force at each unique signed strain-rate value, computed from
# SNR-passing points only (the line is a derived summary, not raw data, so it
# is legitimate to restrict it to the trustworthy subset while still showing
# every raw point). A dashed reference arrow/annotation marks the Hill-expected
# direction (force falls as concentric velocity rises) for visual comparison.
#
# FINDING (revised after actually building and looking at this plot -- the
# PI's skepticism was justified; an earlier verbal claim that bass18 "reproduces
# the correct Hill-type ordering" was based on comparing concentric-vs-eccentric
# PAIRS at each |v| and ignoring the V=0 isometric anchor entirely. Once V=0 is
# included, the full curve is NOT a clean monotonic decline):
# bass18 right: eccentric plateau ~1.65-1.76 N (-127 to -382 %/s) -> DROPS
# SHARPLY to 0.56 N at V=0 -> back UP to 1.37 N (127) -> DOWN to 1.03 N (255)
# -> back UP to 1.83 N (382). A "W", not a Hill hyperbola. Root cause of the
# V=0 dip: its 3 SNR-passing contributors are ALL from the SAME weaker
# trial-set (bender_06/07/08/09/10, the 142-427 deg/s series) while the
# STRONGER trial (bender_03, the 107-320 deg/s series) contributes V=0 values
# of 1.69/2.41 N -- fully consistent with the eccentric plateau -- but FAILS
# the SNR gate (2.08/2.55 < 3.0, its V=0 hold is brief/embedded, not a
# sustained tetanus). So the V=0 notch is a CROSS-TRIAL FATIGUE CONFOUND
# (comparing V=0 from weak trials against moving points from a MIX of strong
# and weak trials), not evidence against Hill -- but it means this pooled-
# across-trials comparison cannot cleanly test the Hill relationship; a
# within-single-trial comparison (bender_03 alone: V=0 ~2.05 N vs its own
# eccentric ~2.03-2.33 N vs its own concentric 0.99-3.63 N, high-variance/
# low-SNR ~2.4-2.8) is suggestive but not itself SNR-clean either. bass16/17
# show a PEAK AT V=0 with decline on BOTH sides (a tent, not Hill's plateau-
# then-monotonic-decline, and not flat noise either) -- unexplained, flagged.
# CONCLUSION: the moving-ramp-only pairwise comparison (eccentric > concentric
# at 127/255 %/s) remains directionally Hill-consistent for bass18, but the
# FULL curve including V=0 is confounded by cross-trial fatigue and does NOT
# cleanly demonstrate the Hill relationship as a single monotonic curve --
# this needs a within-trial-only test (each trial's OWN V=0 + moving steps,
# not pooled across trials) before it can be called settled.
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
                is.finite(.data$strain_rate_pct_s), is.finite(.data$activation_snr)) |>
  dplyr::mutate(confident = .data$activation_snr >= snr_min)

means <- d |>
  dplyr::filter(.data$confident) |>
  dplyr::group_by(.data$fish, .data$muscle_side, .data$strain_rate_pct_s) |>
  dplyr::summarise(force_geom_N = mean(.data$force_geom_N, na.rm = TRUE), n = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(.data$fish, .data$muscle_side, .data$strain_rate_pct_s)

p <- ggplot() +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_vline(xintercept = 0, color = "grey80", linetype = "dashed") +
  geom_point(data = d, aes(strain_rate_pct_s, force_geom_N, color = muscle_side, alpha = confident), size = 1.3) +
  geom_line(data = means, aes(strain_rate_pct_s, force_geom_N, color = muscle_side), linewidth = 1.0) +
  geom_point(data = means, aes(strain_rate_pct_s, force_geom_N, color = muscle_side), size = 2.2, shape = 21, fill = "white") +
  scale_color_manual(values = sides_pal, name = "muscle side") +
  scale_alpha_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.2), name = "confident (SNR >= 3)") +
  facet_wrap(~fish, scales = "free_y") +
  labs(title = "Isovelocity Hill check: is force-velocity monotonic-DECREASING (eccentric >= isometric >= concentric), not a bell?",
       subtitle = "Solid line = mean of SNR-passing points at each velocity (isovelocity's own V=0 hold + moving ramps only, pointwise angle-matched passive). Faint dots = SNR-failing (shown, not dropped). Hill target: HIGH on the left (eccentric/lengthening), falling toward LOW on the right (concentric/high shortening velocity) -- a bell (high at V=0, low at both extremes) would be the WRONG shape here (that's FL's target).",
       x = "Muscle shortening strain rate (%/s; + = concentric/shortening, - = eccentric/lengthening, 0 = isometric V=0 anchor)",
       y = "Muscle force (N, geometric u_hat)") +
  theme_bw(10) + theme(plot.subtitle = element_text(size = 7.5), legend.position = "bottom")
ggsave(file.path(OUT, "isovhillcheck.png"), p, width = 12, height = 5, dpi = 140)

cli::cli_h2("Isovelocity Hill check -- SNR-passing summary means by fish/side/velocity")
print(as.data.frame(means), digits = 3)
cli::cli_alert_success("Wrote isovhillcheck.png to {OUT}")
