# summary_fv_fl_best_within_individual.R
# PI request 2026-07-24 (follow-up to prototype_fv_fl_sono_length.R and
# diag_fv_fl_ztorque_vs_uhat.R): build canonical-textbook-style FL and FV
# curves (see the attached Hill FL/FV reference) from the BEST-quality
# points, WITHIN A SINGLE INDIVIDUAL, in the standard sign conventions.
#
# Three fixes the PI called for, all implemented here:
#
#  1. FORCE sign convention. The "messy / negative force" in the earlier
#     zTorque-vs-uHat figure came from plotting the UNfolded force_zTorque_N.
#     Here force = muscle_force_Nm from build_segmented_step_summary() -- the
#     active-minus-passive torque already sign-FOLDED into the recruited
#     muscle's contraction-positive frame (force_sign) -- i.e. positive =
#     real active tension. Normalised to F/Fmax within the individual
#     (Fmax = that individual's isometric peak).
#
#  2. AXIS conventions (match the canonical plot):
#       FL: x = LENGTH strain, shorter muscle to the LEFT (negative), longer
#           to the RIGHT (positive). strain_len = (L - L0)/L0*100, so this is
#           the NEGATIVE of the shortening-positive sono strain.
#       FV: x = velocity, SHORTENING negative, LENGTHENING positive
#           (= d(strain_len)/dt). Normalised to V/Vmax.
#
#  3. FV constant-velocity window. Rather than a blind middle-60%-of-stim
#     heuristic (which diluted the sono slope ~10x), the constant-velocity
#     segment is detected from the COMMANDED angle (|angular velocity| >=
#     CV_FRAC of the step's peak |angular velocity|, longest contiguous run)
#     and that SAME time window is applied to the sono length to measure the
#     real muscle strain rate -- exactly as the PI directed.
#
# BEST-POINT SELECTION (principled, then post-hoc rationalised in the log):
#   - RIGHT muscle only (sono's only channel).
#   - genuine activation: muscle_force_Nm > 0 (active exceeded passive).
#   - FV only: the commanded ramp was actually realised in the real muscle
#     length -- sono-length-vs-time linear fit over the CV window has
#     R^2 >= FV_R2_MIN AND the window is >= FV_MIN_SAMPLES long.
#   - WITHIN-INDIVIDUAL: the individual with the best combined FL+FV coverage
#     is chosen for the figure (all individuals' points saved to CSV).
#
# Run with:  Rscript R/summary_fv_fl_best_within_individual.R
# Outputs -> figs_summary/ + per-step CSV in data_processed/.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr); library(ggplot2); library(patchwork); library(cli); library(rhdf5); library(signal)
})

.pipeline_root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.pipeline_root, f))
src("paths_config.R")
src("00_load_bender_flat.R")
src("01_calibrate.R")
src("02_deconvolve.R")
src("muscle_geometry.R")
src("03_analyze.R")
src("parse_trial_filename.R")
src("plot_force_vs_time.R")
src("fit_fv_fl.R")
src("plot_strain_validation.R")
src("plot_angle_sono_validation.R")

OUT_DIR      <- FIGS_SUMMARY_DIR
DATA_OUT_DIR <- file.path(.crittergripper_root(), "02_processed", "data_processed")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SPECIMEN_SUBFOLDERS <- c(bass16 = BASS16_RAW_SUBFOLDER, bass17 = BASS17_RAW_SUBFOLDER, bass18 = BASS18_RAW_SUBFOLDER)
SPECIMEN_COLORS     <- c(bass16 = "#1b9e77", bass17 = "#d95f02", bass18 = "#7570b3")

DEACT_WIN_S     <- 0.5
SONO_LP_HZ      <- 40.0
CV_FRAC         <- 0.6    # constant-velocity = |omega| >= CV_FRAC * peak |omega|
FV_R2_MIN       <- 0.90   # sono length must track a straight line over the CV window
FV_MIN_SAMPLES  <- 8L

.butter_lp <- function(x, cutoff_hz, sr_hz, order = 4L) {
  nyq <- sr_hz / 2.0
  if (cutoff_hz >= nyq) return(x)
  ok <- is.finite(x); if (sum(ok) < 20L) return(x)
  filt <- signal::butter(order, cutoff_hz / nyq, type = "low")
  out <- x; out[ok] <- signal::filtfilt(filt, x[ok]); out
}

.read_sono_geom <- function(filename) {
  h5 <- rhdf5::H5Fopen(filename, "H5F_ACC_RDONLY"); on.exit(try(rhdf5::H5Fclose(h5), silent = TRUE), add = TRUE)
  m <- tryCatch(rhdf5::h5readAttributes(h5, "/metadata"), error = function(e) list())
  d1 <- function(v, def = NA_real_) { v <- suppressWarnings(as.numeric(v[1L])); if (length(v)==0L||is.na(v)) def else v }
  list(ai_rate_hz = d1(m[["daq_ai_sample_rate_hz"]], 1000.0),
       lidx_pos_motor = d1(m[["daq_specimen_lateral_index_on_positive_motor_side"]]),
       lidx_right = d1(m[["daq_specimen_side_index_right"]]))
}

# Longest contiguous run of TRUE in a logical vector -> start/end indices.
.longest_run <- function(mask) {
  mask[is.na(mask)] <- FALSE
  if (!any(mask)) return(NULL)
  r <- rle(mask)
  ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1L
  true_idx <- which(r$values)
  best <- true_idx[which.max(r$lengths[true_idx])]
  list(start = starts[best], end = ends[best])
}

.collect_one_trial <- function(f, specimen, category) {
  trial_id <- tools::file_path_sans_ext(basename(f))
  td <- tryCatch(load_bender_flat(f, do_filter = TRUE, loadtorques = "x"), error = function(e) NULL)
  if (is.null(td)) return(tibble())
  tau <- tryCatch(deconvolve_bender(f, hub_path = NULL, verbose = FALSE), error = function(e) NULL)
  if (is.null(tau)) return(tibble())
  N <- min(nrow(td), length(tau)); td <- td[seq_len(N), , drop = FALSE]
  td$torque_inertia_corrected_Nm <- tau[seq_len(N)]
  attr(td, "Filename") <- f

  res <- tryCatch(
    if (category == "isometric") analyze_isometric(td) else analyze_isovelocity(td),
    error = function(e) { cli::cli_warn("{trial_id}: analyze failed: {conditionMessage(e)}"); NULL })
  if (is.null(res) || is.null(res$step_summary) || nrow(res$step_summary) == 0L) return(tibble())
  ss  <- res$step_summary
  tdA <- res$td
  if (!("torque_inertia_corrected_Nm" %in% names(tdA)) && nrow(tdA) == nrow(td))
    tdA$torque_inertia_corrected_Nm <- td$torque_inertia_corrected_Nm

  geom <- .read_sono_geom(f)
  tdA <- tryCatch(attach_sono_strain(tdA, f, geom$lidx_right, geom$lidx_pos_motor), error = function(e) tdA)
  if (!("sono_right_mm" %in% names(tdA)) || all(!is.finite(tdA$sono_right_mm))) return(tibble())
  L0 <- .sono_reference_length_mm(tdA$angle.deg, tdA$sono_right_mm)
  if (!is.finite(L0) || L0 <= 0) return(tibble())
  sono_filt <- .butter_lp(tdA$sono_right_mm, SONO_LP_HZ, geom$ai_rate_hz)
  # LENGTH strain (lengthening positive): shorter < 0, longer > 0.
  tdA$strain_len_pct <- (sono_filt - L0) / L0 * 100.0

  # per-step measures
  rows <- lapply(seq_len(nrow(ss)), function(i) {
    s <- ss[i, ]
    if (is.na(s$contraction_mode) || !is.finite(s$stim_t0_s) || !is.finite(s$stim_t1_s)) return(NULL)
    act <- tdA$step_number == s$step_number & tdA$t.s >= s$stim_t0_s & tdA$t.s <= (s$stim_t1_s + DEACT_WIN_S)
    if (!any(act, na.rm = TRUE)) return(NULL)
    tt <- tdA$t.s[act]; sl <- tdA$strain_len_pct[act]; ang <- tdA$angle.deg[act]
    ok <- is.finite(tt) & is.finite(sl)
    if (sum(ok) < 5L) return(NULL)

    force_sign <- if ("force_sign" %in% names(s)) s$force_sign else 1
    pass_Nm    <- if ("passive_force_Nm" %in% names(s)) s$passive_force_Nm else NA_real_
    tq <- if ("torque_inertia_corrected_Nm" %in% names(tdA)) tdA$torque_inertia_corrected_Nm[act] else rep(NA_real_, length(tt))

    if (category == "isometric") {
      x_val <- mean(sl[ok])          # length strain held (%, lengthening +)
      fit_r2 <- NA_real_; n_win <- sum(ok)
      force_Nm <- s$muscle_force_Nm  # mean active over hold is fine (no motion)
    } else {
      # constant-velocity window from COMMANDED angle, applied to BOTH sono and force
      om <- c(NA, diff(ang) / diff(tt))          # deg/s
      peak <- max(abs(om), na.rm = TRUE)
      cvmask <- is.finite(om) & abs(om) >= CV_FRAC * peak
      run <- .longest_run(cvmask)
      if (is.null(run) || (run$end - run$start + 1L) < FV_MIN_SAMPLES) return(NULL)
      w <- seq.int(run$start, run$end)
      tt_w <- tt[w]; sl_w <- sl[w]; tq_w <- tq[w]
      okw <- is.finite(tt_w) & is.finite(sl_w)
      if (sum(okw) < FV_MIN_SAMPLES) return(NULL)
      lmf <- stats::lm(sl_w[okw] ~ tt_w[okw])
      x_val <- unname(stats::coef(lmf)[2L])       # d(strain_len)/dt (%/s), lengthening +
      fit_r2 <- suppressWarnings(summary(lmf)$r.squared)
      n_win <- sum(okw)
      # active force measured WITHIN the constant-velocity window, phase-matched
      # to the velocity, minus the velocity-matched passive baseline, folded to
      # contraction-positive. Falls back to the whole-window mean if unavailable.
      force_Nm <- if (all(is.finite(tq_w[okw])) && is.finite(pass_Nm))
        force_sign * (mean(tq_w[okw], na.rm = TRUE) - pass_Nm) else s$muscle_force_Nm
    }
    tibble::tibble(
      specimen = specimen, trial_id = trial_id, category = category,
      step_number = s$step_number, muscle_side = s$muscle_side,
      contraction_mode = s$contraction_mode,
      x_val = x_val, fit_r2 = fit_r2, n_win = n_win,
      muscle_force_Nm = force_Nm,
      muscle_force_wholewin_Nm = s$muscle_force_Nm
    )
  })
  dplyr::bind_rows(rows)
}

cli::cli_h1("Collecting isometric (FL) + isovelocity (FV) steps, canonical conventions")
all_steps <- purrr::imap_dfr(SPECIMEN_SUBFOLDERS, function(subfolder, specimen) {
  manifest <- parse_trial_directory(raw_source_dir(subfolder))
  purrr::pmap_dfr(list(manifest$fullpath, manifest$protocol), function(f, proto) {
    if (!proto %in% c("isometric", "isovelocity")) return(tibble())
    out <- tryCatch(.collect_one_trial(f, specimen, proto), error = function(e) {
      cli::cli_alert_danger("{basename(f)}: {conditionMessage(e)}"); tibble() })
    if (nrow(out) > 0L) cli::cli_alert_success("{tools::file_path_sans_ext(basename(f))}: {nrow(out)} steps")
    out
  })
})

# ---- best-point gates ----
all_steps <- all_steps |> dplyr::filter(.data$muscle_side == "right", is.finite(.data$muscle_force_Nm))
sel <- all_steps |>
  dplyr::filter(.data$muscle_force_Nm > 0) |>
  dplyr::filter(.data$category == "isometric" | (.data$fit_r2 >= FV_R2_MIN & .data$n_win >= FV_MIN_SAMPLES))
write.csv(sel, file.path(DATA_OUT_DIR, "fv_fl_best_within_individual_steps.csv"), row.names = FALSE)

# ---- pick the best individual: most combined FL+FV points, must have both ----
coverage <- sel |>
  dplyr::group_by(.data$specimen) |>
  dplyr::summarise(n_fl = sum(.data$category == "isometric"),
                    n_fv = sum(.data$category == "isovelocity"),
                    has_both = n_fl > 0 & n_fv > 0,
                    score = ifelse(has_both, n_fl + n_fv, 0L), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(.data$score))
cli::cli_h1("Per-individual best-point coverage")
print(coverage)
best_specimen <- coverage$specimen[1L]
cli::cli_alert_info("Chosen individual for the canonical figure: {best_specimen}")

# =============================================================================
# OVERVIEW (diagnostic): every individual, points coloured by trial, so the
# best-individual choice can be judged directly. Force normalised PER
# specimen+panel to its own peak; FV velocity normalised per specimen.
# =============================================================================
ov <- sel |>
  dplyr::group_by(.data$specimen, .data$category) |>
  dplyr::mutate(F_norm = .data$muscle_force_Nm / max(.data$muscle_force_Nm, na.rm = TRUE)) |>
  dplyr::group_by(.data$specimen) |>
  dplyr::mutate(V_norm = .data$x_val / max(abs(.data$x_val[.data$category == "isovelocity"]), na.rm = TRUE)) |>
  dplyr::ungroup() |>
  dplyr::mutate(trial_short = sub("^\\d{4}-\\d{2}-\\d{2}_[^_]+_", "", .data$trial_id))

ov_fl <- dplyr::filter(ov, .data$category == "isometric")
ov_fv <- dplyr::filter(ov, .data$category == "isovelocity")

pOV_fl <- ggplot(ov_fl, aes(x = .data$x_val, y = .data$F_norm, color = .data$trial_short)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE, span = 1.0, color = "grey30", linewidth = 0.7) +
  geom_point(size = 2.2) +
  facet_wrap(~ specimen, nrow = 1) +
  labs(title = "Force-Length by individual (isometric)",
       x = "Muscle length strain (% L0)   <- shorter | longer ->",
       y = "Force / specimen+panel peak", color = "trial") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

pOV_fv <- ggplot(ov_fv, aes(x = .data$V_norm, y = .data$F_norm, color = .data$trial_short)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
  geom_smooth(aes(group = 1), method = "loess", se = FALSE, span = 1.1, color = "grey30", linewidth = 0.7) +
  geom_point(size = 2.2) +
  facet_wrap(~ specimen, nrow = 1) +
  labs(title = "Force-Velocity by individual (isovelocity)",
       x = expression(Normalised~velocity~V/V[max]~~"(shortening < 0 < lengthening)"),
       y = "Force / specimen+panel peak", color = "trial") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")

ov_fig <- pOV_fl / pOV_fv + patchwork::plot_annotation(
  title = "FL/FV candidate points -- all individuals, coloured by trial (best-point gates applied)",
  subtitle = "Gates: right muscle, F>0, FV sono-ramp R^2>=0.90. Grey line = per-specimen loess. Use to judge which individual gives the cleanest within-individual curves.",
  theme = theme(plot.title = element_text(face = "bold", size = 12)))
ov_out <- file.path(FIGS_DIAGNOSTIC_DIR, "FL_FV_candidates_all_individuals.png")
ggplot2::ggsave(ov_out, ov_fig, width = 13, height = 8, dpi = 150)
cli::cli_alert_success("Saved {ov_out}")

df <- dplyr::filter(sel, .data$specimen == best_specimen)
col_one <- unname(SPECIMEN_COLORS[best_specimen])

# ---- normalise WITHIN individual, PER PANEL to that panel's own peak ----
# NOTE: isometric F0 is normally the shared reference, but in this individual
# the isometric holds produced anomalously little force relative to the
# dynamic contractions (F0 ~= 1/17 of the eccentric peak -- a weak-isometric
# artifact, not a physiological FL/FV feature). To show the SHAPES faithfully
# we normalise each panel to its own within-individual maximum active tension.
# FL: pooling multiple isometric trials mixes a strongly- and a
# weakly-activated series (e.g. bass18 bender_04 vs bender_11), smearing the
# curve vertically. Use the SINGLE best-activated isometric trial (highest
# peak active tension) within the individual -- the truest length-tension set.
fl_all <- dplyr::filter(df, .data$category == "isometric")
best_fl_trial <- fl_all |>
  dplyr::group_by(.data$trial_id) |>
  dplyr::summarise(peak = max(.data$muscle_force_Nm, na.rm = TRUE), n = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(.data$peak)) |>
  dplyr::slice(1L) |> dplyr::pull(.data$trial_id)
cli::cli_alert_info("FL uses best-activated isometric trial: {best_fl_trial}")
fl <- dplyr::filter(fl_all, .data$trial_id == best_fl_trial)
fv <- dplyr::filter(df, .data$category == "isovelocity")
FL_max <- max(fl$muscle_force_Nm, na.rm = TRUE); if (!is.finite(FL_max) || FL_max <= 0) FL_max <- 1
FV_max <- max(fv$muscle_force_Nm, na.rm = TRUE); if (!is.finite(FV_max) || FV_max <= 0) FV_max <- 1
fl$F_norm <- fl$muscle_force_Nm / FL_max
fv$F_norm <- fv$muscle_force_Nm / FV_max
Vmax <- max(abs(fv$x_val), na.rm = TRUE); if (!is.finite(Vmax) || Vmax <= 0) Vmax <- 1
fv$V_norm <- fv$x_val / Vmax

# =============================================================================
# FL panel (canonical: shorter left, longer right)
# =============================================================================
pFL <- ggplot(fl, aes(x = .data$x_val, y = .data$F_norm)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
  { if (nrow(fl) >= 4) geom_smooth(method = "loess", se = FALSE, span = 1.0, color = col_one, linewidth = 0.9) } +
  geom_point(size = 3, color = col_one) +
  labs(title = "A. Force-Length (isometric holds)", x = "Muscle length strain (% L0)   <- shorter | longer ->",
       y = expression(Force~"/"~panel~peak)) +
  coord_cartesian(ylim = c(0, 1.05)) +
  theme_bw(base_size = 12)

# =============================================================================
# FV panel (canonical: shortening negative, lengthening positive)
# =============================================================================
pFV <- ggplot(fv, aes(x = .data$V_norm, y = .data$F_norm)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey70") +
  { if (nrow(fv) >= 4) geom_smooth(method = "loess", se = FALSE, span = 1.0, color = col_one, linewidth = 0.9) } +
  geom_point(size = 3, color = col_one) +
  annotate("text", x = -0.55, y = 1.02, label = "Shortening", hjust = 0.5, size = 3.4, color = "grey30") +
  annotate("text", x =  0.55, y = 1.02, label = "Lengthening", hjust = 0.5, size = 3.4, color = "grey30") +
  labs(title = "B. Force-Velocity (isovelocity ramps)", x = expression(Normalised~velocity~V/V[max]~~"(shortening < 0 < lengthening)"),
       y = expression(Force~"/"~panel~peak)) +
  coord_cartesian(ylim = c(0, 1.10)) +
  theme_bw(base_size = 12)

fig <- pFL + pFV + patchwork::plot_annotation(
  title = sprintf("Force-Length and Force-Velocity from best within-individual points (%s, RIGHT muscle)", best_specimen),
  subtitle = sprintf("Force = sign-folded active-minus-(velocity-matched)passive single-axis tension (positive); each panel normalised to its own within-individual peak.\nFL = the single best-activated isometric trial (%s) to avoid pooling weakly-activated holds. FV velocity AND force both read from the commanded\nconstant-velocity window (|omega| >= %.0f%% of peak) applied to sono; only clean sono ramps (R^2 >= %.2f) with genuine activation (F>0). n(FL)=%d, n(FV)=%d.",
                     sub("^\\d{4}-\\d{2}-\\d{2}_[^_]+_", "", best_fl_trial), 100*CV_FRAC, FV_R2_MIN, nrow(fl), nrow(fv)),
  theme = theme(plot.title = element_text(face = "bold", size = 13))
)

fout <- file.path(OUT_DIR, "FL_FV_best_within_individual.png")
ggplot2::ggsave(fout, fig, width = 12, height = 5.4, dpi = 150)
cli::cli_alert_success("Saved {fout}")
cli::cli_alert_success("summary_fv_fl_best_within_individual.R complete")
