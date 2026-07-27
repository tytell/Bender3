# diag_isometric_anchoredPooledRate_qc.R
# QUALITY-CONTROL CHECK (PI-directed, 2026-07-27) -- does NOT override
# `legacy` in any production figure. Tests a new candidate passive-baseline
# method, `anchoredPooledRate`, across ALL 5 isometric trials (bass16 x2,
# bass17 x1, bass18 x2), scoring it against `legacy` and `baselineInterp` on
# two criteria the PI named as important-but-not-non-negotiable: (1) sign-
# flip rate across methods, (2) concentric/eccentric symmetry.
#
# BACKGROUND (see analysis_muscle_force_vector_log.md, 2026-07-27 addenda):
# `legacy` (pre-stim-only baseline) was shown to have a DIRECTION-DEPENDENT
# bias -- it OVERESTIMATES |muscle_force| on concentric steps (14/14 in the
# bass17 example) and UNDERESTIMATES on eccentric steps (only 4/14) --
# because it ignores a real, measured, side-AGNOSTIC positional drift
# (passive torque relaxes toward zero throughout a held step, symmetric in
# both bend directions; the concentric/eccentric LABEL is just which side is
# scoring that same raw bend, not a property of the drift itself).
# `baselineInterp` (linear pre+post interpolation) already corrects for
# this drift and was NOT direction-biased in the same example.
# `pooledRegressionBaseline` (session-wide absolute-torque regression, see
# summary_fl_fv_tension_pooled_regression_baseline.R) disagreed in SIGN with
# both simpler methods on the one step stress-tested -- diagnosed as
# cross-step session drift contaminating one step's ABSOLUTE prediction.
#
# METHOD UNDER TEST: anchoredPooledRate. Keeps interp's anchor (each step's
# OWN pre-stim window value as the intercept -- this is what keeps interp
# reliable and legacy's bias small) but replaces the RATE of relaxation
# (normally each step's own noisy 2-point (post-pre)/(t_post-t_pre) slope,
# as used by baselineInterp) with a POOLED, same-trial, cross-step SMOOTHED
# rate-vs-operating_point model. Rationale: a rate estimated from just 2
# points per step is noisy; pooling many steps' rates (which should vary
# smoothly with bend angle, all measuring the SAME side-agnostic positional
# drift) should reduce that noise WITHOUT reintroducing pooled_regression's
# failure mode, because the ABSOLUTE passive level is still anchored
# locally (this step's own pre-window), never predicted from cross-step
# data directly.
#   local_rate_i    = (post_i - pre_i) / (t_post_mid_i - t_pre_mid_i)      [step's own, noisy]
#   rate_pooled(op) = lm(local_rate ~ operating_point + I(operating_point^2)), fit per TRIAL
#   passive_at_active_i = pre_i + rate_pooled(operating_point_i) * (t_active_mid_i - t_pre_mid_i)
#   muscle_force_Nm_anchoredPooledRate = force_sign_i * (active_force_Nm_i - passive_at_active_i)
#
# SCOPE (PI-directed 2026-07-27): isometric ONLY. Isovelocity is deferred --
# it already has a genuine no-stim `velocity_matched` comparison for most
# steps, a fundamentally different (self-contained, not extrapolated)
# problem; revisit separately.
#
# DECISION RULE (PI-directed): this script is diagnostic/QC ONLY. Overriding
# `legacy` in any production figure requires (a) this comparison plot,
# (b) a log addendum documenting the justification, BOTH before the switch
# -- not yet done here regardless of this script's result.
#
# Run with:  Rscript R/diag_isometric_anchoredPooledRate_qc.R
# Output -> figs_diagnostic/isometric_anchoredPooledRate_qc_{comparison,
#           signflip,symmetry}.png
#        -> data_processed/isometric_anchoredPooledRate_qc_allsteps.csv

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(ggplot2); library(cli); library(patchwork)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")
src("00_load_bender_flat.R")
src("01_calibrate.R")
src("02_deconvolve.R")
src("muscle_geometry.R")
src("plot_force_vs_time.R")
src("03_analyze.R")
src("parse_trial_filename.R")
src("dynamic_trial_precondition.R")

OUT_DIR      <- FIGS_DIAGNOSTIC_DIR
DATA_OUT_DIR <- file.path(.crittergripper_root(), "02_processed", "data_processed")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SPECIMEN_SUBFOLDERS <- c(bass16 = BASS16_RAW_SUBFOLDER, bass17 = BASS17_RAW_SUBFOLDER, bass18 = BASS18_RAW_SUBFOLDER)

.load_one <- function(f) {
  td <- load_bender_flat(f, do_filter = TRUE, loadtorques = "x")
  tau <- deconvolve_bender(f, hub_path = NULL, verbose = FALSE)
  N <- min(nrow(td), length(tau)); td <- td[seq_len(N), , drop = FALSE]
  td$torque_inertia_corrected_Nm <- tau[seq_len(N)]
  attr(td, "Filename") <- f
  td
}

# ==========================================================================
# anchoredPooledRate -- see module header for the full derivation. Keeps the
# intermediate .t_*_mid_s / local_rate_qc columns PERSISTENT (not deleted)
# so the cross-trial variant below (fit_anchored_crosstrial_rate()) can
# reuse the same per-step local-rate estimates when pooling across a
# specimen's multiple trials, instead of recomputing them.
# ==========================================================================
fit_anchored_pooled_rate <- function(steps, deactivation_window_s = 0.5) {
  steps$t_pre_mid_s    <- (steps$t_pre_baseline_start_s + steps$t_pre_baseline_end_s) / 2
  steps$t_post_mid_s   <- (steps$t_post_baseline_start_s + steps$t_post_baseline_end_s) / 2
  steps$t_active_mid_s <- (steps$stim_t0_s + (steps$stim_t1_s + deactivation_window_s)) / 2
  steps$local_rate_qc  <- (steps$post_force_Nm_static - steps$passive_force_Nm_static) /
    (steps$t_post_mid_s - steps$t_pre_mid_s)

  ok <- is.finite(steps$local_rate_qc) & is.finite(steps$operating_point)
  if (sum(ok) < 4L || dplyr::n_distinct(steps$operating_point[ok]) < 3L) {
    steps$muscle_force_Nm_anchoredPooledRate <- NA_real_
    steps$rate_pooled_qc <- NA_real_
  } else {
    fit <- tryCatch(
      stats::lm(local_rate_qc ~ operating_point + I(operating_point^2), data = steps[ok, ]),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      steps$muscle_force_Nm_anchoredPooledRate <- NA_real_
      steps$rate_pooled_qc <- NA_real_
    } else {
      steps$rate_pooled_qc <- as.numeric(predict(fit, newdata = data.frame(operating_point = steps$operating_point)))
      passive_at_active <- steps$passive_force_Nm_static +
        steps$rate_pooled_qc * (steps$t_active_mid_s - steps$t_pre_mid_s)
      steps$muscle_force_Nm_anchoredPooledRate <- steps$force_sign * (steps$active_force_Nm - passive_at_active)
    }
  }
  steps
}

# ==========================================================================
# anchoredPooledRate_crossTrial (PI directive, 2026-07-27 follow-up): SAME
# anchor (each step's own pre-stim window) as anchoredPooledRate, but the
# rate model now pools local_rate_qc ~ operating_point ACROSS ALL of one
# SPECIMEN's isometric trials (never across specimens/different fish) --
# more data again (bass16: 2 trials/24 steps, bass18: 2 trials/26 steps;
# bass17 has only 1 isometric trial, so this is a no-op there, which
# doubles as an internal check: bass17's cross-trial and per-trial numbers
# should come out IDENTICAL). Purpose: diagnose whether anchoredPooledRate's
# small residual eccentric-side bias (see log, 2026-07-27) is a per-trial,
# small-sample artifact that more data smooths out, or persists even with
# more pooling.
# ==========================================================================
fit_anchored_crosstrial_rate <- function(all_steps, deactivation_window_s = 0.5) {
  all_steps |>
    dplyr::group_by(.data$specimen) |>
    dplyr::group_modify(function(df, key) {
      ok <- is.finite(df$local_rate_qc) & is.finite(df$operating_point)
      if (sum(ok) < 4L || dplyr::n_distinct(df$operating_point[ok]) < 3L) {
        df$rate_pooled_crossTrial <- NA_real_
        df$muscle_force_Nm_anchoredPooledRate_crossTrial <- NA_real_
        return(df)
      }
      fit <- tryCatch(
        stats::lm(local_rate_qc ~ operating_point + I(operating_point^2), data = df[ok, ]),
        error = function(e) NULL
      )
      if (is.null(fit)) {
        df$rate_pooled_crossTrial <- NA_real_
        df$muscle_force_Nm_anchoredPooledRate_crossTrial <- NA_real_
        return(df)
      }
      df$fit_r2_crossTrial <- summary(fit)$r.squared
      df$rate_pooled_crossTrial <- as.numeric(predict(fit, newdata = data.frame(operating_point = df$operating_point)))
      passive_at_active <- df$passive_force_Nm_static +
        df$rate_pooled_crossTrial * (df$t_active_mid_s - df$t_pre_mid_s)
      df$muscle_force_Nm_anchoredPooledRate_crossTrial <- df$force_sign * (df$active_force_Nm - passive_at_active)
      df
    }) |>
    dplyr::ungroup()
}

# ==========================================================================
# Collect all 5 isometric trials, all 3 specimens
# ==========================================================================
cli::cli_h1("Collecting isometric steps, all specimens, 3-way method comparison")

.collect_isometric_qc <- function(specimen, source_dir) {
  manifest <- parse_trial_directory(source_dir)
  files <- manifest$fullpath[manifest$protocol == "isometric"]
  purrr::map_dfr(files, function(f) {
    trial_id <- tools::file_path_sans_ext(basename(f))
    td <- tryCatch(.load_one(f), error = function(e) NULL)
    if (is.null(td)) return(tibble())
    built <- tryCatch(build_segmented_step_summary(td, f), error = function(e) NULL)
    if (is.null(built) || nrow(built$step_summary) == 0L) return(tibble())
    steps <- fit_anchored_pooled_rate(built$step_summary)
    steps$specimen  <- specimen
    steps$trial_id  <- trial_id
    steps$trial_num <- extract_bender_trial_num(trial_id)
    steps
  })
}

all_steps <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  .collect_isometric_qc(specimen, raw_source_dir(subfolder))
})
all_steps$precondition <- classify_session_precondition(all_steps$specimen, all_steps$trial_num)

cli::cli_h1("Fitting anchoredPooledRate_crossTrial (per-specimen, pooled across that specimen's trials)")
all_steps <- fit_anchored_crosstrial_rate(all_steps)
r2_by_specimen <- all_steps |> dplyr::distinct(.data$specimen, .data$fit_r2_crossTrial)
print(r2_by_specimen)

# Restrict comparison to real per-side steps with all 4 methods available.
qc <- all_steps |>
  dplyr::filter(.data$muscle_side %in% c("left", "right"),
               .data$contraction_mode %in% c("concentric", "eccentric"),
               is.finite(.data$muscle_force_Nm), is.finite(.data$muscle_force_Nm_interp),
               is.finite(.data$muscle_force_Nm_anchoredPooledRate),
               is.finite(.data$muscle_force_Nm_anchoredPooledRate_crossTrial))
cli::cli_alert_info("{nrow(qc)} concentric/eccentric steps with all 4 methods available, across {dplyr::n_distinct(qc$trial_id)} trials")

readr::write_csv(qc, file.path(DATA_OUT_DIR, "isometric_anchoredPooledRate_qc_allsteps.csv"))

# ==========================================================================
# QC 1: sign-flip rate across the 4 methods, per trial
# ==========================================================================
cli::cli_h1("QC 1: sign-flip rate (legacy vs. baselineInterp vs. anchoredPooledRate vs. anchoredPooledRate_crossTrial)")

qc <- qc |>
  dplyr::mutate(
    sign_legacy = sign(.data$muscle_force_Nm),
    sign_interp = sign(.data$muscle_force_Nm_interp),
    sign_anchoredPooledRate = sign(.data$muscle_force_Nm_anchoredPooledRate),
    sign_anchoredPooledRate_crossTrial = sign(.data$muscle_force_Nm_anchoredPooledRate_crossTrial),
    n_pos = (.data$sign_legacy > 0) + (.data$sign_interp > 0) + (.data$sign_anchoredPooledRate > 0) +
      (.data$sign_anchoredPooledRate_crossTrial > 0),
    n_neg = (.data$sign_legacy < 0) + (.data$sign_interp < 0) + (.data$sign_anchoredPooledRate < 0) +
      (.data$sign_anchoredPooledRate_crossTrial < 0),
    sign_flip_any = .data$n_pos > 0 & .data$n_neg > 0
  )

signflip_by_trial <- qc |>
  dplyr::group_by(.data$specimen, .data$trial_id) |>
  dplyr::summarise(n_steps = dplyr::n(), n_flip = sum(.data$sign_flip_any), pct_flip = 100 * mean(.data$sign_flip_any),
                   .groups = "drop")
print(signflip_by_trial)

# ==========================================================================
# QC 2: concentric/eccentric asymmetry -- does anchoredPooledRate correct
# legacy's direction-dependent bias (PI-confirmed on bass17, 2026-07-27)?
# For each method, |muscle_force| relative to interp's own |muscle_force|
# (interp is the within-trial reference for "drift-corrected" -- see log),
# split by contraction_mode.
# ==========================================================================
cli::cli_h1("QC 2: concentric/eccentric asymmetry vs. baselineInterp reference")

asym <- qc |>
  dplyr::mutate(
    legacy_minus_interp = abs(.data$muscle_force_Nm) - abs(.data$muscle_force_Nm_interp),
    anchoredPooledRate_minus_interp = abs(.data$muscle_force_Nm_anchoredPooledRate) - abs(.data$muscle_force_Nm_interp),
    anchoredPooledRate_crossTrial_minus_interp = abs(.data$muscle_force_Nm_anchoredPooledRate_crossTrial) - abs(.data$muscle_force_Nm_interp)
  ) |>
  dplyr::group_by(.data$contraction_mode) |>
  dplyr::summarise(
    n = dplyr::n(),
    legacy_bigger_pct = 100 * mean(abs(.data$muscle_force_Nm) > abs(.data$muscle_force_Nm_interp)),
    legacy_mean_diff  = mean(.data$legacy_minus_interp),
    anchoredPooledRate_bigger_pct = 100 * mean(abs(.data$muscle_force_Nm_anchoredPooledRate) > abs(.data$muscle_force_Nm_interp)),
    anchoredPooledRate_mean_diff  = mean(.data$anchoredPooledRate_minus_interp),
    anchoredPooledRate_crossTrial_bigger_pct = 100 * mean(abs(.data$muscle_force_Nm_anchoredPooledRate_crossTrial) > abs(.data$muscle_force_Nm_interp)),
    anchoredPooledRate_crossTrial_mean_diff  = mean(.data$anchoredPooledRate_crossTrial_minus_interp),
    .groups = "drop"
  )
print(asym, width = Inf)

# ==========================================================================
# Comparison plot (required before any promotion decision, per PI directive)
# ==========================================================================
cli::cli_h1("Building comparison figure")

side_colors <- c(left = "#1d4ed8", right = "#b91c1c")

p1 <- ggplot(qc, aes(x = .data$shortening_strain_pct, y = .data$muscle_force_Nm, color = .data$muscle_side)) +
  geom_point(aes(shape = .data$trial_id), size = 1.8, alpha = 0.6) +
  scale_color_manual(values = side_colors, guide = "none") +
  labs(title = "legacy (pre-only)", x = NULL, y = "Muscle force (N*m)") +
  theme_bw(base_size = 10) + theme(legend.position = "none")

p2 <- ggplot(qc, aes(x = .data$shortening_strain_pct, y = .data$muscle_force_Nm_interp, color = .data$muscle_side)) +
  geom_point(aes(shape = .data$trial_id), size = 1.8, alpha = 0.6) +
  scale_color_manual(values = side_colors, guide = "none") +
  labs(title = "baselineInterp (pre+post linear)", x = NULL, y = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "none")

p3 <- ggplot(qc, aes(x = .data$shortening_strain_pct, y = .data$muscle_force_Nm_anchoredPooledRate, color = .data$muscle_side)) +
  geom_point(aes(shape = .data$trial_id), size = 1.8, alpha = 0.6) +
  scale_color_manual(values = side_colors, guide = "none") +
  labs(title = "anchoredPooledRate (own-pre anchor,\npooled rate PER TRIAL)", x = NULL, y = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "none")

p4 <- ggplot(qc, aes(x = .data$shortening_strain_pct, y = .data$muscle_force_Nm_anchoredPooledRate_crossTrial, color = .data$muscle_side)) +
  geom_point(aes(shape = .data$trial_id), size = 1.8, alpha = 0.6) +
  scale_color_manual(values = side_colors, name = "Muscle side") +
  labs(title = "anchoredPooledRate_crossTrial (own-pre\nanchor, pooled rate PER SPECIMEN)",
       x = "Muscle shortening strain (%)", y = NULL, shape = "Trial") +
  theme_bw(base_size = 10)

p_signflip <- ggplot(signflip_by_trial, aes(x = .data$trial_id, y = .data$pct_flip, fill = .data$specimen)) +
  geom_col() +
  geom_text(aes(label = sprintf("%d/%d", .data$n_flip, .data$n_steps)), vjust = -0.3, size = 3) +
  labs(title = "Sign-flip rate across 4 methods, by trial", x = NULL, y = "% steps with a sign disagreement") +
  coord_cartesian(ylim = c(0, max(signflip_by_trial$pct_flip, 10) * 1.25)) +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 30, hjust = 1))

p_asym <- ggplot(asym, aes(x = .data$contraction_mode)) +
  geom_col(aes(y = .data$legacy_mean_diff, fill = "legacy"), position = position_nudge(x = -0.3), width = 0.22) +
  geom_col(aes(y = .data$anchoredPooledRate_mean_diff, fill = "anchoredPooledRate"), position = position_nudge(x = 0), width = 0.22) +
  geom_col(aes(y = .data$anchoredPooledRate_crossTrial_mean_diff, fill = "anchoredPooledRate_crossTrial"), position = position_nudge(x = 0.3), width = 0.22) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c(legacy = "grey50", anchoredPooledRate = "#0f766e", anchoredPooledRate_crossTrial = "#7c3aed"), name = "Method") +
  labs(title = "Mean(|method| - |baselineInterp|) by contraction mode\n(0 = perfectly agrees with interp; symmetric bars = no CE/ECC bias)",
       x = "Contraction mode", y = "Mean |force| difference from interp (N*m)") +
  theme_bw(base_size = 10)

top_row <- p1 + p2 + p3 + p4 + patchwork::plot_layout(nrow = 1)
bottom_row <- p_signflip + p_asym + patchwork::plot_layout(nrow = 1, widths = c(1.3, 1))
p_full <- top_row / bottom_row +
  patchwork::plot_annotation(
    title = "Isometric passive-baseline method QC: legacy vs. baselineInterp vs. anchoredPooledRate vs. anchoredPooledRate_crossTrial (all 5 trials, bass16/17/18)",
    subtitle = "QUALITY-CONTROL CHECK ONLY -- does not override legacy in any production figure (see analysis_muscle_force_vector_log.md, 2026-07-27)",
    theme = theme(plot.title = element_text(size = 12, face = "bold"), plot.subtitle = element_text(size = 9))
  )

fout <- file.path(OUT_DIR, "isometric_anchoredPooledRate_qc_comparison.png")
ggplot2::ggsave(fout, p_full, width = 13, height = 9, dpi = 150)
cli::cli_alert_success("Saved {fout}")

cli::cli_alert_success("diag_isometric_anchoredPooledRate_qc.R complete.")
cli::cli_alert_info("Sign-flip summary:"); print(signflip_by_trial)
cli::cli_alert_info("Concentric/eccentric asymmetry summary:"); print(asym)
