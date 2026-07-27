# diag_isometric_baseline_hypothesis_tests.R
# PI directive (2026-07-27): test the TWO hypotheses left open by the
# cross-trial pooling result (see analysis_muscle_force_vector_log.md,
# "cross-trial rate pooling test" addendum). Recap: anchoredPooledRate
# removes legacy's concentric bias almost completely, but leaves a small
# eccentric-side residual that got MORE consistent (88.2% -> 94.1%) when
# more data was pooled across trials -- i.e. it behaves like structure, not
# small-sample noise. Two candidate explanations were logged:
#
#   H1: the relaxation-RATE-vs-bend-angle relationship is genuinely
#       ASYMMETRIC between the two bend directions, and the current
#       `local_rate ~ operating_point + I(operating_point^2)` model cannot
#       represent that shape, so it under/over-corrects on one side.
#
#   H2: `baselineInterp` -- the REFERENCE every comparison above was scored
#       against -- is itself the biased one on eccentric steps, in which
#       case anchoredPooledRate isn't drifting away from truth, it's
#       drifting away from a biased yardstick.
#
# TEST DESIGN (all reference-free where possible -- the whole difficulty in
# this investigation is that no isometric no-stim ground truth exists, so
# these tests are built to avoid needing one):
#
#   Panel A (H1): local_rate vs operating_point per specimen, with the
#     CURRENT symmetric model overlaid against a DIRECTION-SPLIT model
#     (separate quadratic fit for negative-bend and positive-bend steps).
#     If the split model wins on AIC, the rate really is direction-
#     asymmetric and H1 is supported.
#
#   Panel B (H1): residuals of the CURRENT model, split by bend direction.
#     A correct model leaves residuals centred on zero in BOTH directions;
#     equal-and-opposite residual means = direction asymmetry the model is
#     absorbing incorrectly.
#
#   Panel C (H2): the MATCHED-BEND test, and the cleanest one available.
#     Passive drift is PI-confirmed to be purely positional and agnostic to
#     which muscle is recruited. The isometric protocol conveniently visits
#     the SAME signed operating_point twice per trial -- once with the side
#     whose label makes it "concentric", once with the side that makes it
#     "eccentric" (34 such matched pairs across the corpus). At matched
#     bend, local_rate MUST agree if drift is truly label-agnostic. It is
#     `baselineInterp` that consumes each step's own local_rate, so if
#     these matched pairs disagree systematically by label, interp inherits
#     a label-dependent bias -> H2 supported. CONFOUND, stated plainly: the
#     two members of a matched pair occur at different times in the
#     session, so a session-order effect is not separable from a true
#     label effect here; the pair ORDER is reported alongside.
#
#   Panel D (arbiter): reference-free REPEATABILITY. bass16 and bass18 each
#     ran TWO isometric trials at identical nominal operating points, so
#     the same (side, contraction_mode, strain) condition is measured
#     twice. Whichever method gives the most reproducible muscle force
#     across those repeats is extracting the most signal and the least
#     baseline artifact -- no ground truth needed. CAVEAT: real
#     between-trial fatigue also inflates this spread, but it inflates it
#     EQUALLY for all four methods (identical groups, paired comparison),
#     so the between-method ranking stays fair even though the absolute
#     level does not mean "method noise" alone. bass17 contributes nothing
#     here (one isometric trial only).
#
# SCOPE: isometric only (PI-directed). Diagnostic/QC ONLY -- promotes
# nothing, overrides nothing; `legacy` remains the production default.
#
# Run with:  Rscript R/diag_isometric_baseline_hypothesis_tests.R
# Output -> figs_diagnostic/isometric_baseline_hypothesis_tests.png
#        -> data_processed/isometric_baseline_hypothesis_tests_matchedbend.csv

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

# Same per-step quantities the QC script builds (deliberate duplication --
# the established pattern for diagnostics in this repo), plus all four
# candidate muscle-force columns so Panel D can score them side by side.
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
    s$specimen <- specimen
    s$trial_id <- trial_id
    s
  })
}

cli::cli_h1("Collecting isometric steps (all specimens)")
all_steps <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(sub, sp) .collect(sp, raw_source_dir(sub)))
cli::cli_alert_info("{nrow(all_steps)} steps across {dplyr::n_distinct(all_steps$trial_id)} isometric trials")

# Attach the two anchored variants (per-trial rate model, per-specimen rate
# model) so Panel D can score all four methods.
all_steps <- all_steps |>
  dplyr::group_by(.data$trial_id) |>
  dplyr::group_modify(function(df, key) {
    ok <- is.finite(df$local_rate_qc) & is.finite(df$operating_point)
    df$muscle_force_Nm_anchoredPooledRate <- NA_real_
    if (sum(ok) >= 4L && dplyr::n_distinct(df$operating_point[ok]) >= 3L) {
      fit <- tryCatch(stats::lm(local_rate_qc ~ operating_point + I(operating_point^2), data = df[ok, ]),
                      error = function(e) NULL)
      if (!is.null(fit)) {
        rate <- as.numeric(predict(fit, newdata = data.frame(operating_point = df$operating_point)))
        df$muscle_force_Nm_anchoredPooledRate <- df$force_sign *
          (df$active_force_Nm - (df$passive_force_Nm_static + rate * (df$t_active_mid_s - df$t_pre_mid_s)))
      }
    }
    df
  }) |>
  dplyr::ungroup() |>
  dplyr::group_by(.data$specimen) |>
  dplyr::group_modify(function(df, key) {
    ok <- is.finite(df$local_rate_qc) & is.finite(df$operating_point)
    df$muscle_force_Nm_anchoredPooledRate_crossTrial <- NA_real_
    df$resid_symmetric <- NA_real_
    if (sum(ok) >= 4L && dplyr::n_distinct(df$operating_point[ok]) >= 3L) {
      fit <- tryCatch(stats::lm(local_rate_qc ~ operating_point + I(operating_point^2), data = df[ok, ]),
                      error = function(e) NULL)
      if (!is.null(fit)) {
        rate <- as.numeric(predict(fit, newdata = data.frame(operating_point = df$operating_point)))
        df$muscle_force_Nm_anchoredPooledRate_crossTrial <- df$force_sign *
          (df$active_force_Nm - (df$passive_force_Nm_static + rate * (df$t_active_mid_s - df$t_pre_mid_s)))
        df$resid_symmetric <- df$local_rate_qc - rate
      }
    }
    df
  }) |>
  dplyr::ungroup()

# ==========================================================================
# H1 test: symmetric quadratic vs. direction-split quadratic, per specimen
# ==========================================================================
cli::cli_h1("H1: is the relaxation rate vs. bend-angle relationship direction-asymmetric?")

h1 <- purrr::map_dfr(unique(all_steps$specimen), function(sp) {
  d <- dplyr::filter(all_steps, .data$specimen == sp,
                     is.finite(.data$local_rate_qc), is.finite(.data$operating_point))
  if (nrow(d) < 8L) return(tibble())
  d$bend_dir <- ifelse(d$operating_point < 0, "negative bend", "positive bend")
  m_sym   <- stats::lm(local_rate_qc ~ operating_point + I(operating_point^2), data = d)
  m_split <- stats::lm(local_rate_qc ~ bend_dir * (operating_point + I(operating_point^2)), data = d)
  tibble::tibble(
    specimen   = sp,
    n          = nrow(d),
    r2_sym     = summary(m_sym)$r.squared,
    r2_split   = summary(m_split)$r.squared,
    aic_sym    = stats::AIC(m_sym),
    aic_split  = stats::AIC(m_split),
    delta_aic  = stats::AIC(m_split) - stats::AIC(m_sym),   # negative => split model wins
    p_split    = tryCatch(stats::anova(m_sym, m_split)$`Pr(>F)`[2], error = function(e) NA_real_)
  )
})
print(h1, width = Inf)

# Fitted curves for Panel A
curve_df <- purrr::map_dfr(unique(all_steps$specimen), function(sp) {
  d <- dplyr::filter(all_steps, .data$specimen == sp,
                     is.finite(.data$local_rate_qc), is.finite(.data$operating_point))
  if (nrow(d) < 8L) return(tibble())
  d$bend_dir <- ifelse(d$operating_point < 0, "negative bend", "positive bend")
  m_sym   <- stats::lm(local_rate_qc ~ operating_point + I(operating_point^2), data = d)
  m_split <- stats::lm(local_rate_qc ~ bend_dir * (operating_point + I(operating_point^2)), data = d)
  grid <- tibble::tibble(operating_point = seq(min(d$operating_point), max(d$operating_point), length.out = 120)) |>
    dplyr::mutate(bend_dir = ifelse(.data$operating_point < 0, "negative bend", "positive bend"))
  dplyr::bind_rows(
    dplyr::mutate(grid, rate = as.numeric(predict(m_sym, newdata = grid)),   model = "symmetric (current)"),
    dplyr::mutate(grid, rate = as.numeric(predict(m_split, newdata = grid)), model = "direction-split")
  ) |> dplyr::mutate(specimen = sp)
})

# ==========================================================================
# H2 test: matched-bend pairs -- same signed operating_point, same trial,
# one step labelled concentric and one eccentric.
# ==========================================================================
cli::cli_h1("H2: at MATCHED bend, does local_rate depend on the concentric/eccentric label?")

matched <- all_steps |>
  dplyr::filter(.data$contraction_mode %in% c("concentric", "eccentric"),
                is.finite(.data$local_rate_qc)) |>
  dplyr::group_by(.data$specimen, .data$trial_id, .data$operating_point) |>
  dplyr::filter(dplyr::n_distinct(.data$contraction_mode) == 2L) |>
  dplyr::summarise(
    rate_conc  = mean(.data$local_rate_qc[.data$contraction_mode == "concentric"], na.rm = TRUE),
    rate_ecc   = mean(.data$local_rate_qc[.data$contraction_mode == "eccentric"], na.rm = TRUE),
    step_conc  = dplyr::first(.data$step_number[.data$contraction_mode == "concentric"]),
    step_ecc   = dplyr::first(.data$step_number[.data$contraction_mode == "eccentric"]),
    .groups    = "drop"
  ) |>
  dplyr::mutate(
    rate_diff_ecc_minus_conc = .data$rate_ecc - .data$rate_conc,
    ecc_came_later = .data$step_ecc > .data$step_conc
  )

readr::write_csv(matched, file.path(DATA_OUT_DIR, "isometric_baseline_hypothesis_tests_matchedbend.csv"))
cli::cli_alert_info("{nrow(matched)} matched-bend pairs")
cli::cli_alert_info("mean(rate_ecc - rate_conc) = {signif(mean(matched$rate_diff_ecc_minus_conc, na.rm=TRUE),3)}")
h2_t <- tryCatch(stats::t.test(matched$rate_diff_ecc_minus_conc), error = function(e) NULL)
if (!is.null(h2_t)) cli::cli_alert_info("paired difference t-test p = {signif(h2_t$p.value,3)}")

# ==========================================================================
# Panel D: reference-free repeatability across the four methods
# ==========================================================================
cli::cli_h1("Arbiter: repeatability across repeated identical conditions (bass16/bass18)")

method_cols <- c(
  legacy                        = "muscle_force_Nm",
  baselineInterp                = "muscle_force_Nm_interp",
  anchoredPooledRate            = "muscle_force_Nm_anchoredPooledRate",
  anchoredPooledRate_crossTrial = "muscle_force_Nm_anchoredPooledRate_crossTrial"
)

repeat_base <- all_steps |>
  dplyr::filter(.data$muscle_side %in% c("left", "right"),
                .data$contraction_mode %in% c("concentric", "eccentric")) |>
  dplyr::mutate(cond = paste(.data$specimen, .data$muscle_side, .data$contraction_mode,
                             round(.data$shortening_strain_pct, 1), sep = "|")) |>
  dplyr::group_by(.data$cond) |>
  dplyr::filter(dplyr::n() >= 2L, dplyr::n_distinct(.data$trial_id) >= 2L) |>
  dplyr::ungroup()

repeatability <- purrr::imap_dfr(method_cols, function(col, nm) {
  d <- repeat_base |>
    dplyr::filter(is.finite(.data[[col]])) |>
    dplyr::group_by(.data$cond) |>
    dplyr::filter(dplyr::n() >= 2L) |>
    dplyr::summarise(sd_within = stats::sd(.data[[col]], na.rm = TRUE), .groups = "drop")
  tibble::tibble(method = nm, n_conditions = nrow(d), mean_sd_within = mean(d$sd_within, na.rm = TRUE))
}) |>
  dplyr::mutate(method = factor(.data$method, levels = names(method_cols)))
print(repeatability)

# ==========================================================================
# Figure
# ==========================================================================
cli::cli_h1("Building figure")

h1_lab <- h1 |>
  dplyr::mutate(lab = sprintf("%s: dAIC=%+.1f%s", .data$specimen, .data$delta_aic,
                              ifelse(is.finite(.data$p_split) & .data$p_split < 0.05, " *", "")))

pA <- ggplot(dplyr::filter(all_steps, is.finite(.data$local_rate_qc)),
             aes(x = .data$operating_point, y = .data$local_rate_qc)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
  geom_point(size = 1.6, alpha = 0.55, color = "grey30") +
  geom_line(data = curve_df, aes(y = .data$rate, color = .data$model, group = interaction(.data$model, .data$bend_dir)),
            linewidth = 0.9) +
  geom_text(data = h1_lab, aes(x = -Inf, y = Inf, label = .data$lab), hjust = -0.05, vjust = 1.4, size = 3, inherit.aes = FALSE) +
  facet_wrap(~ specimen, nrow = 1, scales = "free_y") +
  scale_color_manual(values = c("symmetric (current)" = "#0f766e", "direction-split" = "#dc2626"), name = "Rate model") +
  labs(title = "H1: is the relaxation rate direction-asymmetric?",
       subtitle = "Negative dAIC = the direction-split model fits better; * = significantly better (F-test p<0.05)",
       x = "Operating point (commanded bend, deg)", y = "Local relaxation rate (N*m/s)") +
  theme_bw(base_size = 10)

resid_df <- all_steps |>
  dplyr::filter(is.finite(.data$resid_symmetric)) |>
  dplyr::mutate(bend_dir = ifelse(.data$operating_point < 0, "negative bend", "positive bend"))
resid_means <- resid_df |>
  dplyr::group_by(.data$specimen, .data$bend_dir) |>
  dplyr::summarise(m = mean(.data$resid_symmetric), .groups = "drop")

pB <- ggplot(resid_df, aes(x = .data$bend_dir, y = .data$resid_symmetric, fill = .data$bend_dir)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_boxplot(alpha = 0.5, outlier.size = 0.8, width = 0.6) +
  geom_point(data = resid_means, aes(y = .data$m), shape = 23, size = 2.6, fill = "black", color = "black") +
  facet_wrap(~ specimen, nrow = 1) +
  scale_fill_manual(values = c("negative bend" = "#2563eb", "positive bend" = "#f59e0b"), guide = "none") +
  labs(title = "H1: residuals of the CURRENT (symmetric) rate model, by bend direction",
       subtitle = "Diamond = mean. Equal-and-opposite means would indicate the symmetric model is absorbing a real asymmetry",
       x = NULL, y = "Residual (N*m/s)") +
  theme_bw(base_size = 10)

pC <- ggplot(matched, aes(x = .data$operating_point, y = .data$rate_diff_ecc_minus_conc, color = .data$specimen)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_point(size = 2.2, alpha = 0.8) +
  labs(title = sprintf("H2: matched-bend test (n=%d pairs)", nrow(matched)),
       subtitle = paste0("Same signed bend, same trial, opposite labels. If drift is label-agnostic these sit on zero.\n",
                         sprintf("mean = %.2e N*m/s%s", mean(matched$rate_diff_ecc_minus_conc, na.rm = TRUE),
                                 if (!is.null(h2_t)) sprintf(", t-test p = %.3g", h2_t$p.value) else "")),
       x = "Operating point (deg)", y = "local_rate(eccentric) - local_rate(concentric)", color = "Specimen") +
  theme_bw(base_size = 10)

pD <- ggplot(repeatability, aes(x = .data$method, y = .data$mean_sd_within, fill = .data$method)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.2e", .data$mean_sd_within)), vjust = -0.4, size = 3) +
  scale_fill_manual(values = c(legacy = "grey50", baselineInterp = "#2563eb",
                               anchoredPooledRate = "#0f766e", anchoredPooledRate_crossTrial = "#7c3aed"),
                    guide = "none") +
  coord_cartesian(ylim = c(0, max(repeatability$mean_sd_within, na.rm = TRUE) * 1.25)) +
  labs(title = sprintf("Arbiter: repeatability across %d repeated conditions (bass16/bass18)", repeatability$n_conditions[1]),
       subtitle = "Mean within-condition SD of muscle force. LOWER = more reproducible. Reference-free (no ground truth needed)",
       x = NULL, y = "Mean within-condition SD (N*m)") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 12, hjust = 1))

p_full <- (pA / pB / (pC | pD)) +
  patchwork::plot_annotation(
    title = "Isometric passive-baseline: testing H1 (direction-asymmetric relaxation rate) and H2 (baselineInterp is the biased reference)",
    subtitle = "DIAGNOSTIC ONLY -- promotes nothing; legacy remains the production default (see analysis_muscle_force_vector_log.md, 2026-07-27)",
    theme = theme(plot.title = element_text(size = 12, face = "bold"), plot.subtitle = element_text(size = 9))
  )

fout <- file.path(OUT_DIR, "isometric_baseline_hypothesis_tests.png")
ggplot2::ggsave(fout, p_full, width = 13, height = 13, dpi = 150)
cli::cli_alert_success("Saved {fout}")

cli::cli_h2("SUMMARY")
cli::cli_alert_info("H1 (direction-asymmetric rate):"); print(h1, width = Inf)
cli::cli_alert_info("H2 (matched-bend label effect): mean = {signif(mean(matched$rate_diff_ecc_minus_conc, na.rm=TRUE),3)}, p = {if (!is.null(h2_t)) signif(h2_t$p.value,3) else NA}")
cli::cli_alert_info("Repeatability arbiter:"); print(repeatability)
