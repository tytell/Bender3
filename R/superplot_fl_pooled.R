# superplot_fl_pooled.R
# Pooled Force-Length SUPERPLOT across individuals (bass16/17/18), combining
# BOTH protocol families (PI-directed, 2026-07-18):
#   - ISOMETRIC : each held step contributes one (length, force) point directly
#                 (shortening_strain_pct is a true positional strain there).
#   - ISOVELOCITY: each ACTIVE ramp is swept -- instantaneous vector force vs
#                 instantaneous ENCODER length, with an ANGLE-matched passive
#                 ramp subtracted POINTWISE (within-trial -> cross-trial SAME
#                 individual at the same signed velocity). Ramps with no
#                 angle-matched passive (static-baseline fallback) are DROPPED,
#                 because a static baseline leaves the large passive-bending
#                 slope in the trace and would swamp the small active signal.
#   NOTE (flagged, per user): isovelocity force is velocity-dependent, so its
#   FL contribution folds some force-VELOCITY behaviour into this force-LENGTH
#   plot. This is intentional pooling, not a bug.
#
# Muscle force is the 6-axis line-of-action (u_hat) tension from
# muscle_force_vector.R, standardized left == right by the PASSIVE-relative
# sign convention (REVISED 2026-07-18): positive/more extreme where the
# active deviation REINFORCES the passive torque's own direction at that
# angle, negative/less extreme where it OPPOSES it -- NOT uniformly positive.
# Two facets:
#   - empirical u_hat            (muscle_force_vector_N)
#   - geometric/longitudinal u_hat (muscle_force_vector_geom_N; for isovelocity
#                                   the geometric fallback IS the longitudinal
#                                   X estimate, matching .mfv_finalize_step()).
#
# Superplot layers (NO model fit -- connect-the-dots only, per cursorrule):
#   - individual points for each trial (per-trial, per strain-bin mean) = dots
#   - one THIN line per trial connecting its dots
#   - one line per INDIVIDUAL (fish) = mean of its trials per bin
#   - one THICKER line = GROUP mean across individuals per bin
#
# Run with:  Rscript R/superplot_fl_pooled.R

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(stringr); library(fs); library(ggplot2); library(rhdf5); library(cli)
})

.root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
for (f in c("paths_config.R", "00_load_bender_flat.R", "01_calibrate.R", "02_deconvolve.R",
            "muscle_geometry.R", "fit_fv_fl.R", "03_analyze.R",
            "parse_trial_filename.R", "plot_strain_validation.R",
            "plot_angle_sono_validation.R", "muscle_force_vector.R")) {
  source(file.path(.root, f))
}

# Defaults come from paths_config.R (single source of truth) -- see that
# file if the OneDrive folder layout ever moves again.
SPECIMEN_DIRS <- c(
  bass16 = raw_source_dir(BASS16_RAW_SUBFOLDER),
  bass17 = raw_source_dir(BASS17_RAW_SUBFOLDER),
  bass18 = raw_source_dir(BASS18_RAW_SUBFOLDER)
)
OUTPUT_DIR <- FIGS_DIAGNOSTIC_DIR
fs::dir_create(OUTPUT_DIR, recurse = TRUE)

STRAIN_BIN_PCT   <- 5.0   # common length-bin width for the mean lines / dots
DEACT_WIN_S      <- MFV_DEACTIVATION_WINDOW_S
VELOCITY_TOL     <- MFV_VELOCITY_MATCH_TOL
SNR_MIN          <- MFV_UHAT_SNR_MIN


# ------------------------------------------------------------------ helpers ---

#' Clamp separation (mm) -> test-section length for the strain conversion.
.read_dclamp_mm <- function(filename) {
  h5 <- rhdf5::H5Fopen(filename, "H5F_ACC_RDONLY")
  on.exit(try(rhdf5::H5Fclose(h5), silent = TRUE), add = TRUE)
  m <- tryCatch(rhdf5::h5readAttributes(h5, "/metadata"), error = function(e) list())
  v <- suppressWarnings(as.numeric(m[["measurement_clamp_separation_millimeter"]][1L]))
  if (length(v) == 0L || is.na(v)) NA_real_ else v
}

#' Load + inertia-deconvolve one trial into a td ready for analyze_*().
.load_trial_td <- function(fullpath) {
  td <- load_bender_flat(fullpath, do_filter = TRUE, loadtorques = "x")
  if (is.null(td) || nrow(td) == 0L) stop("load_bender_flat returned no data")
  tau <- deconvolve_bender(fullpath, hub_path = NULL, verbose = FALSE)
  if (is.null(tau)) stop("deconvolve_bender failed")
  N <- min(nrow(td), length(tau))
  td <- td[seq_len(N), , drop = FALSE]
  td$torque_inertia_corrected_Nm <- tau[seq_len(N)]
  attr(td, "Filename") <- fullpath
  td
}

#' Pointwise (not mean) angle-matched interpolation of a passive ramp's 6
#' channels onto EACH active-sample angle. Same >= 50% angle-overlap gate as
#' .mfv_interp_ramp_onto(); returns list(F=list(3 vecs), T=list(3 vecs)) each of
#' length(ang_act), or NULL if the overlap is insufficient.
.interp_ramp_pointwise <- function(ang_act, ramp) {
  if (length(ang_act) == 0L || all(!is.finite(ang_act))) return(NULL)
  ang_p <- ramp$ang
  if (sum(is.finite(ang_p)) < 2L) return(NULL)
  rng <- range(ang_p, na.rm = TRUE)
  if (mean(ang_act >= rng[1L] & ang_act <= rng[2L], na.rm = TRUE) < 0.5) return(NULL)
  interp1 <- function(vp) {
    ok <- is.finite(ang_p) & is.finite(vp)
    if (sum(ok) < 2L) return(rep(NA_real_, length(ang_act)))
    agg <- stats::aggregate(vp[ok], by = list(a = ang_p[ok]), FUN = mean)
    if (nrow(agg) < 2L) return(rep(NA_real_, length(ang_act)))
    tryCatch(stats::approx(agg$a, agg$x, xout = ang_act, rule = 2)$y,
             error = function(e) rep(NA_real_, length(ang_act)))
  }
  list(F = lapply(ramp$F, interp1), T = lapply(ramp$T, interp1))
}


# ------------------------------------------------- per-specimen extraction ---

#' Isometric (length, force) points for one specimen. Uses the per-step vector
#' force already tension-standardized in attach_vector_muscle_force().
extract_isometric_points <- function(fish, manifest) {
  iso <- dplyr::filter(manifest, protocol == "isometric")
  out <- list()
  for (i in seq_len(nrow(iso))) {
    tid <- iso$trial_id[i]; fp <- iso$fullpath[i]
    res <- tryCatch({
      td <- .load_trial_td(fp)
      r <- analyze_isometric(td, filename = fp); r$trial_id <- tid; r
    }, error = function(e) { cli::cli_warn("iso load {tid}: {conditionMessage(e)}"); NULL })
    if (is.null(res)) next
    vec <- tryCatch(attach_vector_muscle_force(res, fp, "isometric"),
                    error = function(e) { cli::cli_warn("iso vec {tid}: {conditionMessage(e)}"); NULL })
    if (is.null(vec)) next
    ss <- vec$step_summary
    out[[tid]] <- ss |>
      dplyr::filter(.data$muscle_side %in% c("left", "right"),
                    is.finite(.data$shortening_strain_pct)) |>
      dplyr::transmute(
        fish = fish, trial_id = tid, protocol = "isometric",
        muscle_side = .data$muscle_side,
        shortening_strain_pct = .data$shortening_strain_pct,
        force_emp_N  = .data$muscle_force_vector_N,
        force_geom_N = .data$muscle_force_vector_geom_N,
        activation_snr = .data$activation_snr)
  }
  if (length(out) == 0L) tibble::tibble() else dplyr::bind_rows(out)
}

#' Isovelocity instantaneous (length, force) SWEEP points for one specimen.
#' Angle-matched passive subtraction, pointwise, within-trial then cross-trial
#' (same individual, same signed velocity). Static-baseline ramps are dropped.
extract_isovelocity_sweep <- function(fish, manifest) {
  isv <- dplyr::filter(manifest, protocol == "isovelocity")
  if (nrow(isv) == 0L) return(list(points = tibble::tibble(), n_drop = 0L, n_keep = 0L))

  # load each isovelocity trial once (td6 + analyze res + geometry)
  trials <- list()
  for (i in seq_len(nrow(isv))) {
    tid <- isv$trial_id[i]; fp <- isv$fullpath[i]
    tr <- tryCatch({
      td  <- .load_trial_td(fp)
      res <- analyze_isovelocity(td, filename = fp); res$trial_id <- tid
      td6 <- mfv_load_six_axis(fp, res$td)
      if (is.null(td6)) stop("6-axis wrench unavailable")
      geom <- .mfv_read_geom(fp)
      arms <- resolve_muscle_moment_arms(geom$width_mm, geom$depth_mm_raw, geom$dvert_mm)
      dclamp_mm <- .read_dclamp_mm(fp)
      list(trial_id = tid, td6 = td6, ss = res$step_summary, arms = arms,
           strain_factor = (pi / 180) / (dclamp_mm / 1000) * arms$r_m_m * 100)
    }, error = function(e) { cli::cli_warn("isv load {tid}: {conditionMessage(e)}"); NULL })
    if (!is.null(tr)) trials[[tid]] <- tr
  }
  if (length(trials) == 0L) return(list(points = tibble::tibble(), n_drop = 0L, n_keep = 0L))

  # same-individual passive-ramp library (no-stim ramps), keyed by signed velocity
  passive_library <- list()
  for (tr in trials) {
    ns <- .mfv_no_stim_steps(tr$td6, unique(tr$td6$step_number))
    for (sn in ns) {
      op <- tr$ss$operating_point[match(sn, tr$ss$step_number)]
      if (!is.finite(op)) next
      passive_library[[length(passive_library) + 1L]] <- list(
        trial_id = tr$trial_id, operating_point = op,
        ramp = .mfv_ramp_from_step(tr$td6, tr$td6$step_number == sn))
    }
  }

  pts <- list(); n_drop <- 0L; n_keep <- 0L
  for (tr in trials) {
    td6 <- tr$td6; ss <- tr$ss; arms <- tr$arms
    no_stim_steps <- .mfv_no_stim_steps(td6, unique(td6$step_number))
    r_vec <- c(0.0, arms$r_m_m, arms$d_m)
    # longitudinal (u_hat = X) lever, for the geometric/longitudinal facet
    rxu_long <- .mfv_cross(r_vec, c(1, 0, 0)); a_long <- sum(rxu_long^2)

    for (i in seq_len(nrow(ss))) {
      s <- ss[i, ]
      if (!s$muscle_side %in% c("left", "right")) next
      if (!is.finite(s$stim_t0_s) || !is.finite(s$stim_t1_s)) next
      if (!is.finite(s$operating_point) || abs(s$operating_point) < 1e-6) next  # skip holds
      if (!is.finite(s$shortening_value) || s$shortening_value == 0) next
      step_rows <- td6$step_number == s$step_number
      sel <- step_rows & td6$t.s >= s$stim_t0_s & td6$t.s <= (s$stim_t1_s + DEACT_WIN_S)
      if (sum(sel) < 5L) next
      base_win <- td6$t.s >= s$t_pre_baseline_start_s & td6$t.s <= s$t_pre_baseline_end_s

      ang_act <- td6$enc.deg[sel]
      Fa <- list(.mfv_col(td6, "xforce")[sel],  .mfv_col(td6, "yforce")[sel],  .mfv_col(td6, "zforce")[sel])
      Ta <- list(.mfv_col(td6, "xtorque")[sel], .mfv_col(td6, "ytorque")[sel], .mfv_col(td6, "ztorque")[sel])

      # passive: within-trial angle match -> cross-trial (same fish) angle match
      pass <- NULL; psource <- NA_character_
      within_ns <- no_stim_steps[abs(
        ss$operating_point[match(no_stim_steps, ss$step_number)] - s$operating_point) < VELOCITY_TOL]
      if (length(within_ns) > 0L) {
        pass <- .interp_ramp_pointwise(ang_act, .mfv_ramp_from_step(td6, td6$step_number == within_ns[1L]))
        if (!is.null(pass)) psource <- "within_trial"
      }
      if (is.null(pass)) {
        for (lib in passive_library) {
          if (identical(lib$trial_id, tr$trial_id)) next
          if (abs(lib$operating_point - s$operating_point) >= VELOCITY_TOL) next
          cand <- .interp_ramp_pointwise(ang_act, lib$ramp)
          if (!is.null(cand)) { pass <- cand; psource <- "cross_trial"; break }
        }
      }
      if (is.null(pass)) { n_drop <- n_drop + 1L; next }  # static-fallback -> drop
      n_keep <- n_keep + 1L

      dF <- Map(function(a, p) a - p, Fa, pass$F)
      dT <- Map(function(a, p) a - p, Ta, pass$T)
      meandF <- vapply(dF, function(v) mean(v, na.rm = TRUE), numeric(1L))
      meandT <- vapply(dT, function(v) mean(v, na.rm = TRUE), numeric(1L))

      emp   <- uhat_empirical(meandF)
      noise <- .mfv_baseline_force_noise(td6, step_rows, base_win)
      snr   <- if (is.finite(noise) && noise > 0) emp$mag / noise else NA_real_
      uhat_emp <- if (all(is.finite(emp$uhat)) && is.finite(snr) && snr >= SNR_MIN) emp$uhat else c(1, 0, 0)

      # Sign convention (matches .mfv_finalize_step(), REVISED 2026-07-18 per PI
      # biomechanical correction): tension_sign is the sign of the RAW PASSIVE
      # torque projected onto the same r x u_hat axis, NOT an arbitrary
      # per-step "make this positive" rule -- so a sample whose active
      # deviation reinforces the passive direction reads positive/more
      # extreme, and one that opposes it reads negative/less extreme.
      # Ambiguous (passive projection indistinguishable from force noise, e.g.
      # near a zero-velocity/zero-angle window) falls back to the previous
      # always-positive convention.
      fr_mean <- solve_muscle_force_from_wrench(meandF, meandT, uhat_emp, arms$r_m_m, arms$d_m)
      rxu_e <- .mfv_cross(r_vec, uhat_emp); a_e <- sum(rxu_e^2)
      if (!is.finite(a_e) || a_e <= 0) next
      meanpassT <- vapply(pass$T, function(v) mean(v, na.rm = TRUE), numeric(1L))
      pass_moment <- sum(meanpassT * rxu_e) / a_e
      tsign <- if (!is.finite(pass_moment) || !is.finite(noise) || abs(pass_moment) <= noise) {
        if (is.finite(fr_mean$F_moment_N) && fr_mean$F_moment_N < 0) -1 else 1
      } else {
        sign(pass_moment)
      }
      f_emp  <- tsign * (dT[[1]] * rxu_e[1] + dT[[2]] * rxu_e[2] + dT[[3]] * rxu_e[3]) / a_e
      f_geom <- tsign * (dT[[1]] * rxu_long[1] + dT[[2]] * rxu_long[2] + dT[[3]] * rxu_long[3]) / a_long

      side_fold <- sign(s$shortening_value / s$operating_point)
      strain_inst <- ang_act * side_fold * tr$strain_factor

      pts[[length(pts) + 1L]] <- tibble::tibble(
        fish = fish, trial_id = tr$trial_id, protocol = "isovelocity",
        muscle_side = s$muscle_side, passive_source = psource,
        shortening_strain_pct = strain_inst,
        force_emp_N = f_emp, force_geom_N = f_geom,
        activation_snr = snr)
    }
  }
  list(points = if (length(pts) > 0L) dplyr::bind_rows(pts) else tibble::tibble(),
       n_drop = n_drop, n_keep = n_keep)
}


# ---------------------------------------------------------------- assemble ---

cli::cli_h1("Pooled FL superplot: bass16/17/18 (isometric + isovelocity sweep)")
all_points <- list(); drop_report <- list()
for (fish in names(SPECIMEN_DIRS)) {
  cli::cli_h2(fish)
  manifest <- parse_trial_directory(SPECIMEN_DIRS[[fish]])
  iso_pts <- extract_isometric_points(fish, manifest)
  isv     <- extract_isovelocity_sweep(fish, manifest)
  cli::cli_alert_info(
    "{fish}: isometric points = {nrow(iso_pts)}; isovelocity ramps kept = {isv$n_keep}, dropped (static-baseline) = {isv$n_drop}")
  all_points[[fish]] <- dplyr::bind_rows(iso_pts, isv$points)
}
pooled <- dplyr::bind_rows(all_points)
if (nrow(pooled) == 0L) cli::cli_abort("No pooled FL points extracted.")

# long form over the two u_hat methods, then bin length for dots + mean lines
method_levels <- c(empirical = "empirical u_hat",
                   geometric = "geometric / longitudinal u_hat")
long <- pooled |>
  tidyr::pivot_longer(c(force_emp_N, force_geom_N), names_to = "method_raw", values_to = "force_N") |>
  dplyr::mutate(method = factor(dplyr::if_else(method_raw == "force_emp_N",
                                               method_levels[["empirical"]], method_levels[["geometric"]]),
                                levels = method_levels)) |>
  dplyr::filter(is.finite(force_N), is.finite(shortening_strain_pct))

long$strain_bin <- round(long$shortening_strain_pct / STRAIN_BIN_PCT) * STRAIN_BIN_PCT

# per-trial per-bin mean (= the dots, one thin line per trial)
per_trial <- long |>
  dplyr::group_by(method, fish, trial_id, strain_bin) |>
  dplyr::summarise(force_N = mean(force_N, na.rm = TRUE), .groups = "drop")
# per-fish per-bin mean (individual mean) = mean of that fish's trial means
per_fish <- per_trial |>
  dplyr::group_by(method, fish, strain_bin) |>
  dplyr::summarise(force_N = mean(force_N, na.rm = TRUE), .groups = "drop")
# group mean per bin = mean across fish
per_group <- per_fish |>
  dplyr::group_by(method, strain_bin) |>
  dplyr::summarise(force_N = mean(force_N, na.rm = TRUE), n_fish = dplyr::n(), .groups = "drop")

n_iso  <- dplyr::n_distinct(dplyr::filter(pooled, protocol == "isometric")$trial_id)
n_isv  <- dplyr::n_distinct(dplyr::filter(pooled, protocol == "isovelocity")$trial_id)

fish_cols <- c(bass16 = "#1b9e77", bass17 = "#d95f02", bass18 = "#7570b3")

p <- ggplot(mapping = aes(x = strain_bin, y = force_N)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  # per-trial thin connect-the-dots lines + points (individual trials)
  geom_line(data = per_trial,
            aes(group = trial_id, colour = fish), linewidth = 0.3, alpha = 0.35) +
  geom_point(data = per_trial,
             aes(colour = fish), size = 0.8, alpha = 0.4) +
  # per-individual (fish) mean line + points
  geom_line(data = per_fish,
            aes(group = fish, colour = fish), linewidth = 1.0, alpha = 0.95) +
  geom_point(data = per_fish,
             aes(colour = fish), size = 2.1, alpha = 0.95) +
  # group mean across individuals (thicker, black)
  geom_line(data = per_group, aes(group = 1), colour = "black", linewidth = 1.8) +
  geom_point(data = per_group, colour = "black", size = 2.4) +
  facet_wrap(~method, ncol = 2) +
  scale_colour_manual(values = fish_cols, name = "Individual (fish)") +
  labs(
    title = "Pooled Force-Length superplot (6-axis line-of-action muscle force)",
    subtitle = sprintf(
      "isometric steps + isovelocity length-sweep (angle-matched passive), pooled across %d individuals | %d isometric + %d isovelocity trial(s) | %g%% length bins | connect-the-dots, NO fit\nthin = per trial, medium = per individual mean, thick black = group mean | isovelocity folds some force-velocity into force-length (flagged)",
      dplyr::n_distinct(pooled$fish), n_iso, n_isv, STRAIN_BIN_PCT),
    x = "Muscle shortening strain (%, length; + = shortened)",
    y = "Muscle force along u_hat (N, + = reinforces passive direction)") +
  theme_bw(12) +
  theme(legend.position = "right",
        plot.subtitle = element_text(size = 8),
        strip.text = element_text(face = "bold"))

outfile <- file.path(OUTPUT_DIR, "FLsuperplot_isometric_isovelocity_pooled.png")
ggplot2::ggsave(outfile, p, width = 13, height = 6.5, dpi = 150)
cli::cli_alert_success("Wrote {outfile}")
cli::cli_alert_info("Pooled points: {nrow(pooled)} rows | per-trial-bin dots: {nrow(per_trial)/2} | fish: {paste(sort(unique(pooled$fish)), collapse=', ')}")
