# diag_isometric_anchoredRate_v2_dirTime.R
# PI directive (2026-07-27): refit the anchored rate model with the two
# terms the H1/H2 hypothesis tests said were missing, then re-run the QC
# comparison. See analysis_muscle_force_vector_log.md, addendum "H1/H2
# hypothesis tests: eccentric residual explained".
#
# DIAGNOSIS BEING TESTED. The eccentric-side residual of
# `anchoredPooledRate` was traced to TWO defects in its RATE model, neither
# of them a property of eccentric contraction:
#   (1) the real relaxation-rate-vs-bend-angle relationship is
#       DIRECTION-ASYMMETRIC (direction-split model beat the symmetric one
#       in all 3 specimens, dAIC -8.9/-64.8/-31.5), but the model was
#       symmetric-by-construction; and
#   (2) the relaxation RATE ITSELF decays over a session (~17% shrink,
#       27/34 matched pairs, p=1.6e-5), but the model had no time term --
#       and because the protocol runs all concentric blocks BEFORE all
#       eccentric blocks, that session decay masqueraded as an
#       "eccentric" bias.
#
# MODEL UNDER TEST: `anchoredRate_dirTime` (v2). Same anchor as every
# other anchored variant -- each step's OWN pre-stim window supplies the
# absolute passive level, so cross-step data can never contaminate it (the
# failure mode of `pooledRegressionBaseline`). Only the RATE model changes:
#
#   v1 (anchoredPooledRate_crossTrial):
#       local_rate ~ operating_point + I(operating_point^2)
#   v2 (anchoredRate_dirTime):
#       local_rate ~ bend_dir * (operating_point + I(operating_point^2))
#                    + elapsed_in_trial_s + trial_id
#
#   bend_dir  = negative/positive bend (op < 0 vs op >= 0), giving the
#               direction-split shape from H1.
#   elapsed_in_trial_s = seconds since that TRIAL's first step
#               (wall_clock_start based). The session-order decay found in
#               H2 is a WITHIN-trial effect (block order), so within-trial
#               elapsed time is the correct covariate for it.
#   trial_id  = per-trial intercept, absorbing between-trial state
#               differences when pooling a specimen's two trials. Drops out
#               automatically for bass17 (single isometric trial).
#
#   Fit PER SPECIMEN (pooled across that specimen's own isometric trials,
#   never across fish) -- v2 has more parameters, so it needs the larger
#   per-specimen n; the per-trial variant is not attempted for v2.
#
# SUCCESS CRITERION, and why it is not circular. Two scores are computed:
#
#   (a) REFERENCE-FREE: RMSE of (model-predicted rate - that step's OWN
#       measured local rate). This asks only "does the rate model
#       reproduce the rate actually measured at each step", needs no
#       ground truth and no reference method. Lower = better.
#
#   (b) DIVERGENCE FROM baselineInterp, split by contraction mode -- the
#       same metric tracked in the previous two addenda, so results stay
#       comparable. NOTE the expected direction: interp's rate for a step
#       is that step's OWN 2-point measured rate, so it already carries
#       both the direction asymmetry and the session decay. A rate model
#       that correctly captures both should therefore CONVERGE on interp
#       (residual -> 0 in BOTH modes, asymmetry collapsing), not diverge.
#       If v2 converges, the honest reading is that pooling buys noise
#       smoothing rather than bias correction -- so residual SCATTER (SD)
#       is reported alongside the mean, since that is where any remaining
#       benefit of pooling would show up.
#
# SCOPE: isometric only. DIAGNOSTIC ONLY -- promotes nothing, overrides
# nothing; `legacy` remains the production default.
#
# Run with:  Rscript R/diag_isometric_anchoredRate_v2_dirTime.R
# Output -> figs_diagnostic/isometric_anchoredRate_v2_dirTime.png
#        -> data_processed/isometric_anchoredRate_v2_dirTime_allsteps.csv

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
DEACTIVATION_WINDOW_S <- 0.5

.load_one <- function(f) {
  td <- load_bender_flat(f, do_filter = TRUE, loadtorques = "x")
  tau <- deconvolve_bender(f, hub_path = NULL, verbose = FALSE)
  N <- min(nrow(td), length(tau)); td <- td[seq_len(N), , drop = FALSE]
  td$torque_inertia_corrected_Nm <- tau[seq_len(N)]
  attr(td, "Filename") <- f
  td
}

.collect <- function(specimen, source_dir) {
  manifest <- parse_trial_directory(source_dir)
  files <- manifest$fullpath[manifest$protocol == "isometric"]
  purrr::map_dfr(files, function(f) {
    trial_id <- tools::file_path_sans_ext(basename(f))
    td <- tryCatch(.load_one(f), error = function(e) NULL)
    if (is.null(td)) return(tibble())
    built <- tryCatch(build_segmented_step_summary(td, f), error = function(e) NULL)
    if (is.null(built) || nrow(built$step_summary) == 0L) return(tibble())
    s <- built$step_summary
    s$t_pre_mid_s    <- (s$t_pre_baseline_start_s + s$t_pre_baseline_end_s) / 2
    s$t_post_mid_s   <- (s$t_post_baseline_start_s + s$t_post_baseline_end_s) / 2
    s$t_active_mid_s <- (s$stim_t0_s + (s$stim_t1_s + DEACTIVATION_WINDOW_S)) / 2
    s$local_rate_qc  <- (s$post_force_Nm_static - s$passive_force_Nm_static) /
      (s$t_post_mid_s - s$t_pre_mid_s)
    wc <- suppressWarnings(as.POSIXct(s$wall_clock_start, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC"))
    s$elapsed_in_trial_s <- as.numeric(difftime(wc, min(wc, na.rm = TRUE), units = "secs"))
    s$specimen <- specimen
    s$trial_id <- trial_id
    s
  })
}

cli::cli_h1("Collecting isometric steps (all specimens)")
all_steps <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(sub, sp) .collect(sp, raw_source_dir(sub)))
all_steps$bend_dir <- ifelse(all_steps$operating_point < 0, "negative bend", "positive bend")
cli::cli_alert_info("{nrow(all_steps)} steps across {dplyr::n_distinct(all_steps$trial_id)} isometric trials")

# ==========================================================================
# Fit v1 (symmetric, no time) and v2 (direction-split + session time),
# per specimen. Both anchored identically -- only the rate model differs.
# ==========================================================================
cli::cli_h1("Fitting rate models v1 (symmetric) and v2 (direction-split + session time)")

.apply_anchor <- function(df, rate) {
  df$force_sign * (df$active_force_Nm - (df$passive_force_Nm_static + rate * (df$t_active_mid_s - df$t_pre_mid_s)))
}

all_steps <- all_steps |>
  dplyr::group_by(.data$specimen) |>
  dplyr::group_modify(function(df, key) {
    ok <- is.finite(df$local_rate_qc) & is.finite(df$operating_point) & is.finite(df$elapsed_in_trial_s)
    df$rate_v1 <- NA_real_; df$rate_v2 <- NA_real_
    df$muscle_force_Nm_anchoredPooledRate_crossTrial <- NA_real_
    df$muscle_force_Nm_anchoredRate_dirTime <- NA_real_
    if (sum(ok) < 10L) return(df)

    d <- df[ok, ]
    m1 <- tryCatch(stats::lm(local_rate_qc ~ operating_point + I(operating_point^2), data = d),
                   error = function(e) NULL)
    # trial_id only enters when the specimen actually has >1 isometric trial
    f2 <- if (dplyr::n_distinct(d$trial_id) > 1L) {
      local_rate_qc ~ bend_dir * (operating_point + I(operating_point^2)) + elapsed_in_trial_s + trial_id
    } else {
      local_rate_qc ~ bend_dir * (operating_point + I(operating_point^2)) + elapsed_in_trial_s
    }
    m2 <- tryCatch(stats::lm(f2, data = d), error = function(e) NULL)

    if (!is.null(m1)) {
      df$rate_v1 <- as.numeric(predict(m1, newdata = df))
      df$muscle_force_Nm_anchoredPooledRate_crossTrial <- .apply_anchor(df, df$rate_v1)
    }
    if (!is.null(m2)) {
      df$rate_v2 <- as.numeric(predict(m2, newdata = df))
      df$muscle_force_Nm_anchoredRate_dirTime <- .apply_anchor(df, df$rate_v2)
    }
    df$r2_v1 <- if (!is.null(m1)) summary(m1)$r.squared else NA_real_
    df$r2_v2 <- if (!is.null(m2)) summary(m2)$r.squared else NA_real_
    df
  }) |>
  dplyr::ungroup()

# ==========================================================================
# Score (a): REFERENCE-FREE -- does the rate model reproduce each step's
# own measured local rate?
# ==========================================================================
cli::cli_h1("Score (a): rate-model accuracy vs. each step's OWN measured rate (reference-free)")

rate_acc <- all_steps |>
  dplyr::filter(is.finite(.data$local_rate_qc), is.finite(.data$rate_v1), is.finite(.data$rate_v2)) |>
  dplyr::group_by(.data$specimen) |>
  dplyr::summarise(
    n        = dplyr::n(),
    r2_v1    = dplyr::first(.data$r2_v1),
    r2_v2    = dplyr::first(.data$r2_v2),
    rmse_v1  = sqrt(mean((.data$rate_v1 - .data$local_rate_qc)^2)),
    rmse_v2  = sqrt(mean((.data$rate_v2 - .data$local_rate_qc)^2)),
    .groups  = "drop"
  ) |>
  dplyr::mutate(rmse_drop_pct = 100 * (1 - .data$rmse_v2 / .data$rmse_v1))
print(rate_acc, width = Inf)

# ==========================================================================
# Score (b): divergence from baselineInterp, by contraction mode
# ==========================================================================
cli::cli_h1("Score (b): divergence from baselineInterp, by contraction mode")

qc <- all_steps |>
  dplyr::filter(.data$muscle_side %in% c("left", "right"),
                .data$contraction_mode %in% c("concentric", "eccentric"),
                is.finite(.data$muscle_force_Nm), is.finite(.data$muscle_force_Nm_interp),
                is.finite(.data$muscle_force_Nm_anchoredPooledRate_crossTrial),
                is.finite(.data$muscle_force_Nm_anchoredRate_dirTime))
cli::cli_alert_info("{nrow(qc)} concentric/eccentric steps with all methods available")
readr::write_csv(qc, file.path(DATA_OUT_DIR, "isometric_anchoredRate_v2_dirTime_allsteps.csv"))

method_cols <- c(
  legacy                        = "muscle_force_Nm",
  anchoredPooledRate_crossTrial = "muscle_force_Nm_anchoredPooledRate_crossTrial",
  anchoredRate_dirTime          = "muscle_force_Nm_anchoredRate_dirTime"
)

diverg <- purrr::imap_dfr(method_cols, function(col, nm) {
  qc |>
    dplyr::mutate(d = abs(.data[[col]]) - abs(.data$muscle_force_Nm_interp)) |>
    dplyr::group_by(.data$contraction_mode) |>
    dplyr::summarise(n = dplyr::n(),
                     bigger_pct = 100 * mean(.data$d > 0),
                     mean_diff  = mean(.data$d),
                     sd_diff    = stats::sd(.data$d),
                     .groups = "drop") |>
    dplyr::mutate(method = nm)
}) |>
  dplyr::mutate(method = factor(.data$method, levels = names(method_cols)))
print(diverg, width = Inf)

asym_gap <- diverg |>
  dplyr::group_by(.data$method) |>
  dplyr::summarise(asym_gap_pct_points = abs(diff(.data$bigger_pct)),
                   asym_gap_mean_diff  = abs(diff(.data$mean_diff)), .groups = "drop")
cli::cli_alert_info("Concentric-vs-eccentric asymmetry GAP by method (smaller = more symmetric):")
print(asym_gap, width = Inf)

# Sign-flip rate including v2
qc <- qc |>
  dplyr::mutate(
    n_pos = (.data$muscle_force_Nm > 0) + (.data$muscle_force_Nm_interp > 0) +
      (.data$muscle_force_Nm_anchoredPooledRate_crossTrial > 0) + (.data$muscle_force_Nm_anchoredRate_dirTime > 0),
    n_neg = (.data$muscle_force_Nm < 0) + (.data$muscle_force_Nm_interp < 0) +
      (.data$muscle_force_Nm_anchoredPooledRate_crossTrial < 0) + (.data$muscle_force_Nm_anchoredRate_dirTime < 0),
    sign_flip_any = .data$n_pos > 0 & .data$n_neg > 0
  )
signflip <- qc |>
  dplyr::group_by(.data$specimen, .data$trial_id) |>
  dplyr::summarise(n_steps = dplyr::n(), n_flip = sum(.data$sign_flip_any),
                   pct_flip = 100 * mean(.data$sign_flip_any), .groups = "drop")
cli::cli_alert_info("Sign-flip rate (4 methods incl. v2):")
print(signflip)

# ==========================================================================
# Figure
# ==========================================================================
cli::cli_h1("Building figure")

rmse_long <- rate_acc |>
  dplyr::select("specimen", "rmse_v1", "rmse_v2") |>
  tidyr::pivot_longer(c("rmse_v1", "rmse_v2"), names_to = "model", values_to = "rmse") |>
  dplyr::mutate(model = ifelse(.data$model == "rmse_v1", "v1: symmetric, no time", "v2: direction-split + session time"))

pA <- ggplot(rmse_long, aes(x = .data$specimen, y = .data$rmse, fill = .data$model)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_text(data = rate_acc, aes(x = .data$specimen, y = pmax(.data$rmse_v1, .data$rmse_v2),
                                 label = sprintf("-%.0f%%", .data$rmse_drop_pct)),
            vjust = -0.5, size = 3.2, inherit.aes = FALSE) +
  scale_fill_manual(values = c("v1: symmetric, no time" = "#7c3aed",
                               "v2: direction-split + session time" = "#059669"), name = "Rate model") +
  labs(title = "(a) REFERENCE-FREE: does the rate model reproduce each step's OWN measured rate?",
       subtitle = "RMSE of (predicted rate - measured local rate). Lower = better. Label = v2's reduction vs v1",
       x = NULL, y = "RMSE (N*m/s)") +
  theme_bw(base_size = 10)

pB <- ggplot(diverg, aes(x = .data$contraction_mode, y = .data$mean_diff, fill = .data$method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  scale_fill_manual(values = c(legacy = "grey50", anchoredPooledRate_crossTrial = "#7c3aed",
                               anchoredRate_dirTime = "#059669"), name = "Method") +
  labs(title = "(b) THE VERDICT: divergence from baselineInterp, by contraction mode",
       subtitle = "Mean(|method| - |interp|). Bars at zero AND equal across modes = bias removed and asymmetry collapsed",
       x = NULL, y = "Mean |force| difference from interp (N*m)") +
  theme_bw(base_size = 10)

pC <- ggplot(diverg, aes(x = .data$method, y = .data$sd_diff, fill = .data$contraction_mode)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  scale_fill_manual(values = c(concentric = "#0ea5e9", eccentric = "#f97316"), name = "Mode") +
  labs(title = "(c) Scatter of that divergence (SD)",
       subtitle = "Where any noise-smoothing benefit of pooling would show up",
       x = NULL, y = "SD of |force| difference (N*m)") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 12, hjust = 1))

pD <- ggplot(qc, aes(x = .data$shortening_strain_pct, y = .data$muscle_force_Nm_anchoredRate_dirTime,
                     color = .data$muscle_side)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
  geom_point(aes(shape = .data$specimen), size = 2, alpha = 0.75) +
  scale_color_manual(values = c(left = "#1d4ed8", right = "#b91c1c"), name = "Side") +
  labs(title = "(d) Resulting force-length data, anchoredRate_dirTime (v2)",
       x = "Muscle shortening strain (%)", y = "Muscle force (N*m)", shape = "Specimen") +
  theme_bw(base_size = 10)

p_full <- (pA / pB / (pC | pD)) +
  patchwork::plot_annotation(
    title = "Isometric baseline v2: does adding a direction-split term + a session-time term to the rate model fix the eccentric residual?",
    subtitle = "DIAGNOSTIC ONLY -- promotes nothing; legacy remains the production default (see analysis_muscle_force_vector_log.md, 2026-07-27)",
    theme = theme(plot.title = element_text(size = 12, face = "bold"), plot.subtitle = element_text(size = 9))
  )

fout <- file.path(OUT_DIR, "isometric_anchoredRate_v2_dirTime.png")
ggplot2::ggsave(fout, p_full, width = 13, height = 13, dpi = 150)
cli::cli_alert_success("Saved {fout}")

cli::cli_h2("SUMMARY")
print(rate_acc, width = Inf)
print(diverg, width = Inf)
print(asym_gap, width = Inf)
