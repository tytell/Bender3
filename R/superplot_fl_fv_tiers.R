# superplot_fl_fv_tiers.R
# FL / FV at the THREE pooling tiers (PI-requested 2026-07-22), to sit ALONGSIDE
# the existing across-protocol pooled superplots and make the pooling level of
# every claim explicit. Sources the two pooled builders to reuse their exact
# extraction (angle+velocity-matched POINTWISE isovelocity passive, relaxation-
# fit isometric passive, SNR + magnitude-gated F0) -- this file only RE-AGGREGATES
# their already-built `pooled` frames at different granularities; it computes no
# new force.
#
# THE TIERS (see FIGURES_README.md "Individual/Summary"):
#   within-trial          -- each line = ONE trial's own step series. Cleanest
#                            shape: same session, same F0, same passive ramps, NO
#                            cross-trial fatigue/scale normalization. (NEW.)
#   within-protocol pool  -- pooled ACROSS fish but WITHOUT mixing protocols
#                            (isometric-only FL here). Lets each protocol's shape
#                            be read under its OWN settled passive before being
#                            combined. FV's within-protocol pool already IS the
#                            FVsuperplot (isovelocity moving points only). (NEW FL.)
#   across-protocol pool  -- the existing FLsuperplot/FVsuperplot (isometric +
#                            isovelocity V=0 + dynamic L0 folded together). Most
#                            assumption-laden; trust only once per-protocol
#                            passives are settled. (EXISTS -- not rebuilt here.)
#
# All force panels use the GEOMETRIC u_hat (stable; the empirical u_hat is noisy
# for the low-force moving steps -- see FVsuperplot).
#
# Run: Rscript R/superplot_fl_fv_tiers.R  (also refreshes the two pooled superplots)
# Canon token: fltiers / fvtiers -- see FIGURES_README.md.

suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(tibble) })

.root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
source(file.path(.root, "superplot_fl_pooled.R")); pooled_fl <- pooled
source(file.path(.root, "superplot_fv_pooled.R")); pooled_fv <- pooled
OUT <- FIGS_DIAGNOSTIC_DIR

sides_pal <- c(left = "#2563eb", right = "#ea580c")

# Point-level confidence tier (ADDED 2026-07-22, SNR-magnitude conflation
# audit -- this file previously had NO confidence flagging at all, plotting
# every finite point at the same fixed alpha regardless of activation_snr or
# force magnitude; the inconsistency noted in
# analysis_muscle_force_vector_log.md's 2026-07-22 entry). Lines stay at
# their existing fixed alpha (a per-segment tier isn't meaningful for a
# connect-the-dots line); only the point layer is tier-flagged.

# ---- TIER 1a: within-trial FL (isometric trials' own length series) ----
wt_fl <- pooled_fl |>
  dplyr::filter(.data$protocol == "isometric", is.finite(.data$force_geom_N),
                is.finite(.data$shortening_strain_pct), is.finite(.data$activation_snr),
                is.finite(.data$baseline_force_noise_N)) |>
  dplyr::mutate(confidence_tier = mfv_confidence_tier(.data$force_geom_N, .data$activation_snr,
                                                       .data$baseline_force_noise_N))
p1 <- ggplot(wt_fl, aes(shortening_strain_pct, force_geom_N,
                        group = interaction(trial_id, muscle_side), color = muscle_side)) +
  geom_line(linewidth = 0.4, alpha = 0.85) +
  geom_point(aes(alpha = confidence_tier), size = 0.8) +
  scale_alpha_manual(values = MFV_CONFIDENCE_ALPHA, name = "confidence tier (SNR x magnitude)", drop = FALSE) +
  facet_wrap(~fish, scales = "free_y") +
  scale_color_manual(values = sides_pal, name = "muscle side") +
  labs(title = "Within-trial FL (isometric) -- each line = ONE trial's own length series, geometric u_hat, RAW force (no cross-trial normalization)",
       subtitle = "Cleanest FL view: same session/F0/passive. Free y per fish (bass18 is much stronger). Spread across lines = between-trial variability, not shape. Point alpha = confidence tier (see mfv_confidence_tier()); lines shown for ALL tiers.",
       x = "Muscle shortening strain (%, + = shortened)", y = "Muscle force (N, geometric u_hat)") +
  theme_bw(10) + theme(plot.subtitle = element_text(size = 8))
ggsave(file.path(OUT, "fltiers_1_within_trial.png"), p1, width = 12, height = 4.5, dpi = 140)

# ---- TIER 1b: within-trial FV (isovelocity trials' own moving series) ----
wt_fv <- pooled_fv |>
  dplyr::filter(.data$protocol == "isovelocity", .data$strain_rate_pct_s != 0,
                is.finite(.data$force_geom_N), is.finite(.data$strain_rate_pct_s),
                is.finite(.data$activation_snr), is.finite(.data$baseline_force_noise_N)) |>
  dplyr::mutate(confidence_tier = mfv_confidence_tier(.data$force_geom_N, .data$activation_snr,
                                                       .data$baseline_force_noise_N))
p2 <- ggplot(wt_fv, aes(strain_rate_pct_s, force_geom_N,
                        group = interaction(trial_id, muscle_side), color = muscle_side)) +
  geom_hline(yintercept = 0, color = "grey75") +
  geom_line(linewidth = 0.4, alpha = 0.85) +
  geom_point(aes(alpha = confidence_tier), size = 0.8) +
  scale_alpha_manual(values = MFV_CONFIDENCE_ALPHA, name = "confidence tier (SNR x magnitude)", drop = FALSE) +
  facet_wrap(~fish, scales = "free_y") +
  scale_color_manual(values = sides_pal, name = "muscle side") +
  labs(title = "Within-trial FV (isovelocity) -- each line = ONE trial's own velocity series, geometric u_hat, POINTWISE angle-matched passive",
       subtitle = "+ = concentric/shortening, - = eccentric/lengthening. All lines strictly angle-matched (exact or nearest same-sign-velocity stim-off ramp, see analysis log) -- no static-baseline-fallback rows. Point alpha = confidence tier.",
       x = "Muscle shortening strain rate (%/s, signed)", y = "Muscle force (N, geometric u_hat)") +
  theme_bw(10) + theme(plot.subtitle = element_text(size = 8))
ggsave(file.path(OUT, "fvtiers_1_within_trial.png"), p2, width = 12, height = 4.5, dpi = 140)

# ---- TIER 2: within-protocol pooled FL (ISOMETRIC ONLY, across fish, F/F0) ----
# The across-protocol FLsuperplot folds isovelocity V=0 + dynamic L0 into these
# length points; here we keep ONLY the isometric length series so the isometric
# FL shape is visible under its own relaxation-fit passive, un-mixed.
wp_fl <- pooled_fl |>
  dplyr::filter(.data$protocol == "isometric", is.finite(.data$force_geom_norm),
                is.finite(.data$shortening_strain_pct), is.finite(.data$activation_snr),
                is.finite(.data$baseline_force_noise_N)) |>
  dplyr::mutate(confidence_tier = mfv_confidence_tier(.data$force_geom_N, .data$activation_snr,
                                                       .data$baseline_force_noise_N))
bins <- wp_fl |>
  dplyr::mutate(bin = round(.data$shortening_strain_pct / 5) * 5) |>
  dplyr::group_by(.data$fish, .data$bin) |>
  dplyr::summarise(F = mean(.data$force_geom_norm, na.rm = TRUE), n = dplyr::n(), .groups = "drop")
overall <- bins |> dplyr::group_by(.data$bin) |>
  dplyr::summarise(F = mean(.data$F, na.rm = TRUE), .groups = "drop")
p3 <- ggplot() +
  geom_hline(yintercept = c(0, 1), color = "grey80") +
  geom_point(data = wp_fl, aes(shortening_strain_pct, force_geom_norm, color = fish, alpha = confidence_tier), size = 0.9) +
  geom_line(data = bins, aes(bin, F, color = fish), linewidth = 0.7) +
  geom_line(data = overall, aes(bin, F), color = "black", linewidth = 1.3) +
  scale_alpha_manual(values = MFV_CONFIDENCE_ALPHA, name = "confidence tier (SNR x magnitude)", drop = FALSE) +
  labs(title = "Within-protocol pooled FL -- ISOMETRIC ONLY, across fish, F/F0 (geometric u_hat)",
       subtitle = "Companion to the across-protocol FLsuperplot but with NO isovelocity-V=0 / dynamic-L0 mixing -- the isometric FL shape under its own relaxation-fit passive. Black = grand mean over ALL tiers. Point alpha = confidence tier.",
       x = "Muscle shortening strain (%, + = shortened)", y = "F / F0 (dimensionless)") +
  theme_bw(11) + theme(plot.subtitle = element_text(size = 8))
ggsave(file.path(OUT, "fltiers_2_within_protocol_isometriconly.png"), p3, width = 10, height = 5, dpi = 140)

cli::cli_alert_success("Wrote fltiers_1/fvtiers_1/fltiers_2 to {OUT}")
cli::cli_alert_info("within-trial FL rows={nrow(wt_fl)}; within-trial FV rows={nrow(wt_fv)}; within-protocol isometric-FL rows={nrow(wp_fl)}")
