# superplot_fl_pooled.R
# Pooled Force-Length SUPERPLOT across individuals (bass16/17/18).
#
# V=0 ONLY RULE (PI-directed, 2026-07-21, supersedes the isovelocity-SWEEP
# design below): "the rule is that FL superplot only contains moments or
# steps where V = 0." An isovelocity active ramp has genuine nonzero
# velocity throughout its window -- pooling it into a force-LENGTH curve
# folds real force-VELOCITY behavior in regardless of how carefully its
# passive baseline is subtracted, which is exactly what this rule removes.
# Three V=0 sources are pooled, ALL genuinely static/isometric at the moment
# force is sampled:
#   - ISOMETRIC   : each held step (shortening_strain_pct is a true
#                   positional strain -- the whole point of an FL curve).
#   - ISOVELOCITY : NOT the sweep -- only that trial's own embedded V=0
#                   (operating_point == 0, contraction_mode ==
#                   "isometric_zero") holds, which repeat several times per
#                   block. These are genuine isometric-at-L0 contractions,
#                   just recorded inside an isovelocity-protocol file.
#   - DYNAMIC     : the pre-/post-cycling L0 stim bookends every dynamic
#                   trial brackets its sinusoidal sweep with (motor
#                   stationary at 0 deg, angle_sd ~ 0 -- see
#                   extract_dynamic_l0_bookends.R). Added 2026-07-21
#                   (PI-directed): "velocity = zero for those guys too and
#                   they are at L0."
# REMOVED 2026-07-21: the isovelocity continuous-ramp extraction
# (extract_isovelocity_sweep(), angle-matched pointwise passive subtraction)
# -- superseded by extract_isovelocity_zero_points() below. That sweep
# machinery is not used anywhere else in the codebase, so it was deleted
# outright rather than left dormant; see git history if it's ever needed
# again for a DIFFERENT (velocity-aware) plot.
#
# TRIAL-TO-TRIAL FORCE-SCALE NORMALIZATION (added 2026-07-21, PI-directed):
# pooling RAW absolute force across trials assumes every trial sits on the
# same force scale, which isn't guaranteed -- e.g. bass16 has TWO isometric
# trials (09, 15) and TWO isovelocity trials (14, 16) at the IDENTICAL
# nominal operating points, recorded at different points in the session, so
# whatever separates them (fatigue, clamp re-seating, etc.) shows up as pure
# scatter/inconsistency in a raw-force pooled FL curve -- likely a major
# contributor to this plot being "very unhelpful" (PI feedback, 2026-07-21).
# Each point is ALSO expressed as F/F0, F0 = that SAME trial+side's own mean
# force at L0 (isometric: shortening_strain_pct == 0; isovelocity: its own
# isometric_zero holds; dynamic: its own pre-/post-cycling L0 bookends,
# pre+post averaged -- so a real pre-vs-post fatigue drop within one dynamic
# trial shows up as a genuine deviation from F/F0 == 1, not hidden by
# always normalizing against itself) -- classic muscle-physiology
# normalization, corrects cross-trial scale differences while preserving
# each trial's own FL shape. Both RAW and NORMALIZED outputs are written
# side by side (separate files) so the PI can judge whether normalization
# actually helps before it replaces anything -- see the two ggsave calls at
# the bottom of this file.
#
# Muscle force is the 6-axis line-of-action (u_hat) tension from
# muscle_force_vector.R. SIGN CONVENTION (REVISED 2026-07-22, replaces the
# passive-relative reinforcing/opposing rule): muscle_force_vector_N and
# muscle_force_vector_geom_N are reported in their RAW, as-computed sign, no
# per-step correction -- confirmed empirically (328 real steps, 3 fish) NOT
# to be side- or L0-mirrored; mostly positive with a real negative minority
# concentrated in isovelocity, which is genuine signal, not an artifact to
# flip away. See muscle_force_vector.R's module header for the full
# rationale. Two facets:
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
            "plot_angle_sono_validation.R", "muscle_force_vector.R",
            "plot_force_vs_time.R",          # .detect_stim_events() -- extract_dynamic_l0_bookends.R depends on this
            "extract_dynamic_l0_bookends.R")) {
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
SNR_MIN          <- MFV_UHAT_SNR_MIN
# Tolerance (in %) around shortening_strain_pct == 0 (== operating_point == 0
# deg, L0) for identifying a trial's own L0 rep(s) to use as its F0
# normalizer -- added 2026-07-21 for the F/F0 pooling below. A tight
# floating-point tolerance, not a search band (same rationale as
# plot_fatigue_timeline.R's l0_angle_epsilon_deg).
L0_STRAIN_EPSILON_PCT <- 0.5


# ------------------------------------------------------------------ helpers ---

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

# ------------------------------------------------- per-specimen extraction ---

#' Isometric (length, force) points for one specimen. Uses the per-step vector
#' force already tension-standardized in attach_vector_muscle_force().
#'
#' Also attaches TRIAL-OWN-L0 NORMALIZED force (added 2026-07-21, PI-directed
#' -- see module header "TRIAL-TO-TRIAL FORCE-SCALE NORMALIZATION"):
#' force_emp_norm/force_geom_norm = force / that SAME trial+side's own mean
#' force at L0 (shortening_strain_pct == 0, i.e. operating_point == 0 deg).
#' This is the classic F/F0 muscle-physiology normalization (F0 = force at a
#' reference length, NOT necessarily this trial's overall peak) -- distinct
#' from, and NOT the same normalizer as, the "trial's own MAX force" tried
#' and abandoned for the fatigue-timeline diagnostic (plot_fatigue_timeline.R)
#' on 2026-07-21: that used MAX (any angle) as the reference and was tested
#' for TIME-based pooling; this uses L0 SPECIFICALLY (the biomechanical
#' reference length) as the reference and is for LENGTH-based (FL curve)
#' pooling -- same "divide by a per-trial reference force" mechanic, disjoint
#' purpose and disjoint reference point, do not conflate the two.
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
    ss <- vec$step_summary |>
      dplyr::filter(.data$muscle_side %in% c("left", "right"),
                    is.finite(.data$shortening_strain_pct))
    if (nrow(ss) == 0L) next

    # F0 (this trial+side's own mean L0 force) -- L0_EPSILON_PCT matches the
    # commanded-0-deg tolerance already used by plot_fatigue_timeline.R
    # (0.5 deg), converted to this fish's own strain-% units so it stays a
    # tight tolerance around the true commanded 0, not a search band.
    #
    # SNR-GATED (added 2026-07-21, after the first normalized superplot came
    # back with +-1000 spikes): F/F0 blows up whenever F0 itself is small/
    # noisy -- dividing by a near-noise-floor L0 force amplifies noise into a
    # huge ratio. Require the L0 rep(s) feeding F0 to ALSO pass the SAME
    # activation-SNR gate (SNR_MIN == MFV_UHAT_SNR_MIN) used everywhere else
    # in this codebase for confidence, rather than inventing a new threshold.
    # A trial+side with no SNR-passing L0 rep gets F0 = NA (normalized points
    # NA, raw points unaffected) instead of a numerically-unstable ratio.
    f0 <- ss |>
      dplyr::filter(abs(.data$shortening_strain_pct) < L0_STRAIN_EPSILON_PCT,
                    is.finite(.data$activation_snr), .data$activation_snr >= SNR_MIN) |>
      dplyr::group_by(.data$muscle_side) |>
      dplyr::summarise(f0_emp  = mean(.data$muscle_force_vector_N, na.rm = TRUE),
                       f0_geom = mean(.data$muscle_force_vector_geom_N, na.rm = TRUE),
                       .groups = "drop")

    out[[tid]] <- ss |>
      dplyr::left_join(f0, by = "muscle_side") |>
      dplyr::transmute(
        fish = fish, trial_id = tid, protocol = "isometric",
        muscle_side = .data$muscle_side,
        shortening_strain_pct = .data$shortening_strain_pct,
        force_emp_N  = .data$muscle_force_vector_N,
        force_geom_N = .data$muscle_force_vector_geom_N,
        force_emp_norm  = .data$muscle_force_vector_N / .data$f0_emp,
        force_geom_norm = .data$muscle_force_vector_geom_N / .data$f0_geom,
        activation_snr = .data$activation_snr)
  }
  if (length(out) == 0L) tibble::tibble() else dplyr::bind_rows(out)
}

#' Isovelocity's OWN V=0 (isometric_zero) hold points for one specimen --
#' NOT the continuous active-velocity sweep (REMOVED 2026-07-21, PI-directed
#' "V = 0 only" rule -- see module header). Each isovelocity trial embeds
#' repeated V=0 holds (operating_point == 0, contraction_mode ==
#' "isometric_zero") at multiple points in its block -- genuine isometric
#' contractions at L0, just recorded inside an isovelocity-protocol file, so
#' this now goes through the exact SAME attach_vector_muscle_force() path as
#' every other step-level force in this codebase (no bespoke sweep math
#' needed -- that machinery is gone).
extract_isovelocity_zero_points <- function(fish, manifest) {
  isv <- dplyr::filter(manifest, protocol == "isovelocity")
  out <- list()
  for (i in seq_len(nrow(isv))) {
    tid <- isv$trial_id[i]; fp <- isv$fullpath[i]
    res <- tryCatch({
      td <- .load_trial_td(fp)
      r <- analyze_isovelocity(td, filename = fp); r$trial_id <- tid; r
    }, error = function(e) { cli::cli_warn("isv load {tid}: {conditionMessage(e)}"); NULL })
    if (is.null(res)) next
    vec <- tryCatch(attach_vector_muscle_force(res, fp, "isovelocity"),
                    error = function(e) { cli::cli_warn("isv vec {tid}: {conditionMessage(e)}"); NULL })
    if (is.null(vec)) next
    ss0 <- vec$step_summary |>
      dplyr::filter(.data$muscle_side %in% c("left", "right"),
                    .data$contraction_mode == "isometric_zero")
    if (nrow(ss0) == 0L) next

    # F0 (this trial+side's own mean V=0 force, SNR-gated) -- same rationale
    # as extract_isometric_points()'s matching F0 block.
    f0 <- ss0 |>
      dplyr::filter(is.finite(.data$activation_snr), .data$activation_snr >= SNR_MIN) |>
      dplyr::group_by(.data$muscle_side) |>
      dplyr::summarise(f0_emp  = mean(.data$muscle_force_vector_N, na.rm = TRUE),
                       f0_geom = mean(.data$muscle_force_vector_geom_N, na.rm = TRUE),
                       .groups = "drop")

    out[[tid]] <- ss0 |>
      dplyr::left_join(f0, by = "muscle_side") |>
      dplyr::transmute(
        fish = fish, trial_id = tid, protocol = "isovelocity",
        muscle_side = .data$muscle_side,
        shortening_strain_pct = 0.0,
        force_emp_N  = .data$muscle_force_vector_N,
        force_geom_N = .data$muscle_force_vector_geom_N,
        force_emp_norm  = .data$muscle_force_vector_N / .data$f0_emp,
        force_geom_norm = .data$muscle_force_vector_geom_N / .data$f0_geom,
        activation_snr = .data$activation_snr)
  }
  if (length(out) == 0L) tibble::tibble() else dplyr::bind_rows(out)
}

#' Dynamic trials' V=0 L0 bookend contractions for one specimen (added
#' 2026-07-21, PI-directed -- see module header): the stim bursts bracketing
#' every dynamic trial's sinusoidal cycling (motor stationary at 0 deg
#' before it starts and after it ends) are genuine isometric holds at L0.
#' Detected via detect_dynamic_l0_bookends() (extract_dynamic_l0_bookends.R)
#' and routed through attach_vector_muscle_force() with the SAME synthetic-
#' step/td6_override pattern already used by run_fv_fl_power_pipeline.R's
#' fatigue-timeline dynamic integration -- one physics implementation,
#' reused, not re-derived here.
#'
#' F0 = that SAME trial+side's own mean bookend force (pre+post averaged,
#' SNR-gated) -- a dynamic trial only ever contributes L0 points, so this
#' lets a real pre-vs-post WITHIN-trial fatigue drop show up as a genuine
#' deviation from F/F0 == 1 (PI-directed choice, 2026-07-21) rather than
#' being trivially ~1 by normalizing each bookend against itself.
#' One dynamic trial's bookend rows, or NULL. A REAL function (not just a
#' brace-block passed to tryCatch) so its internal `return(NULL)` calls scope
#' to THIS function -- `return()` inside a bare `tryCatch({...})` block
#' returns from the ENCLOSING function instead, which silently aborted the
#' whole caller on the first trial with zero bookends (found + fixed
#' 2026-07-21: every fish's dynamic L0 bookend points were coming back
#' empty).
.dynamic_l0_trial_rows <- function(fp, tid) {
  td <- .load_trial_td(fp)
  bookends <- detect_dynamic_l0_bookends(td)
  if (nrow(bookends) == 0L) return(NULL)
  # attach_vector_muscle_force() filters td6$step_number == s$step_number --
  # dynamic/single_finite loads have no step_number column at all, so one
  # must be injected (constant -- every bookend on this trial shares the one
  # synthetic "step"). Load td6 ONCE per trial, reuse across all bookends.
  td_be <- td; td_be$step_number <- 1L
  td6_be <- mfv_load_six_axis(fp, td_be)
  if (is.null(td6_be)) return(NULL)
  td6_be$step_number <- 1L
  be_rows <- lapply(seq_len(nrow(bookends)), function(bi) {
    b <- bookends[bi, ]
    # category = "isometric" is physically correct here (a real fixed-angle
    # stimulated contraction), matching the fatigue-timeline dynamic
    # integration's exact same choice.
    ss_one <- tibble::tibble(
      step_number = 1L, muscle_side = b$muscle_side, operating_point = b$operating_point,
      contraction_mode = b$contraction_mode,
      stim_t0_s = b$stim_t0_s, stim_t1_s = b$stim_t1_s,
      t_pre_baseline_start_s = b$t_pre_baseline_start_s, t_pre_baseline_end_s = b$t_pre_baseline_end_s
    )
    vec_be <- tryCatch(
      attach_vector_muscle_force(list(step_summary = ss_one, trial_id = tid, td = td_be),
                                 fp, "isometric", td6_override = td6_be),
      error = function(e) NULL)
    if (is.null(vec_be) || nrow(vec_be$step_summary) == 0L) return(NULL)
    vec_be$step_summary[1L, ]
  })
  be_rows <- Filter(Negate(is.null), be_rows)
  if (length(be_rows) == 0L) NULL else dplyr::bind_rows(be_rows)
}

extract_dynamic_l0_points <- function(fish, manifest) {
  dyn <- dplyr::filter(manifest, protocol == "dynamic")
  out <- list()
  for (i in seq_len(nrow(dyn))) {
    tid <- dyn$trial_id[i]; fp <- dyn$fullpath[i]
    rows <- tryCatch(.dynamic_l0_trial_rows(fp, tid),
                     error = function(e) { cli::cli_warn("dyn load {tid}: {conditionMessage(e)}"); NULL })
    if (is.null(rows) || nrow(rows) == 0L) next

    f0 <- rows |>
      dplyr::filter(is.finite(.data$activation_snr), .data$activation_snr >= SNR_MIN) |>
      dplyr::group_by(.data$muscle_side) |>
      dplyr::summarise(f0_emp  = mean(.data$muscle_force_vector_N, na.rm = TRUE),
                       f0_geom = mean(.data$muscle_force_vector_geom_N, na.rm = TRUE),
                       .groups = "drop")

    out[[tid]] <- rows |>
      dplyr::left_join(f0, by = "muscle_side") |>
      dplyr::transmute(
        fish = fish, trial_id = tid, protocol = "dynamic",
        muscle_side = .data$muscle_side,
        shortening_strain_pct = 0.0,
        force_emp_N  = .data$muscle_force_vector_N,
        force_geom_N = .data$muscle_force_vector_geom_N,
        force_emp_norm  = .data$muscle_force_vector_N / .data$f0_emp,
        force_geom_norm = .data$muscle_force_vector_geom_N / .data$f0_geom,
        activation_snr = .data$activation_snr)
  }
  if (length(out) == 0L) tibble::tibble() else dplyr::bind_rows(out)
}


# ---------------------------------------------------------------- assemble ---

cli::cli_h1("Pooled FL superplot: bass16/17/18 (V=0 only -- isometric + isovelocity holds + dynamic bookends)")
all_points <- list()
for (fish in names(SPECIMEN_DIRS)) {
  cli::cli_h2(fish)
  manifest <- parse_trial_directory(SPECIMEN_DIRS[[fish]])
  iso_pts <- extract_isometric_points(fish, manifest)
  isv_pts <- extract_isovelocity_zero_points(fish, manifest)
  dyn_pts <- extract_dynamic_l0_points(fish, manifest)
  cli::cli_alert_info(
    "{fish}: isometric points = {nrow(iso_pts)}; isovelocity V=0 hold points = {nrow(isv_pts)}; dynamic L0 bookend points = {nrow(dyn_pts)}")
  all_points[[fish]] <- dplyr::bind_rows(iso_pts, isv_pts, dyn_pts)
}
pooled <- dplyr::bind_rows(all_points)
if (nrow(pooled) == 0L) cli::cli_abort("No pooled FL points extracted.")

method_levels <- c(empirical = "empirical u_hat",
                   geometric = "geometric / longitudinal u_hat")
n_iso  <- dplyr::n_distinct(dplyr::filter(pooled, protocol == "isometric")$trial_id)
n_isv  <- dplyr::n_distinct(dplyr::filter(pooled, protocol == "isovelocity")$trial_id)
n_dyn  <- dplyr::n_distinct(dplyr::filter(pooled, protocol == "dynamic")$trial_id)
fish_cols <- c(bass16 = "#1b9e77", bass17 = "#d95f02", bass18 = "#7570b3")

#' Shared FL-superplot builder (2026-07-21 refactor -- was inlined once; now
#' called twice, RAW and NORMALIZED, so both share identical binning/plotting
#' logic and can't silently drift apart).
#' @param emp_col,geom_col Which columns of `pooled` to pivot into `force_N`
#'   ("force_emp_N"/"force_geom_N" for raw, "force_emp_norm"/"force_geom_norm"
#'   for F/F0-normalized).
.build_fl_superplot <- function(pooled, emp_col, geom_col, title, subtitle, y_lab) {
  long <- pooled |>
    tidyr::pivot_longer(c(dplyr::all_of(emp_col), dplyr::all_of(geom_col)),
                        names_to = "method_raw", values_to = "force_N") |>
    dplyr::mutate(method = factor(dplyr::if_else(.data$method_raw == emp_col,
                                                 method_levels[["empirical"]], method_levels[["geometric"]]),
                                  levels = method_levels)) |>
    dplyr::filter(is.finite(.data$force_N), is.finite(.data$shortening_strain_pct))
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

  ggplot(mapping = aes(x = strain_bin, y = force_N)) +
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
    labs(title = title, subtitle = subtitle,
        x = "Muscle shortening strain (%, length; + = shortened)",
        y = y_lab) +
    theme_bw(12) +
    theme(legend.position = "right",
          plot.subtitle = element_text(size = 8),
          strip.text = element_text(face = "bold"))
}

p <- .build_fl_superplot(
  pooled, "force_emp_N", "force_geom_N",
  title = "Pooled Force-Length superplot (6-axis line-of-action muscle force)",
  subtitle = sprintf(
    "V=0 ONLY (isometric holds + isovelocity's own V=0 holds + dynamic L0 bookends), pooled across %d individuals | %d isometric + %d isovelocity + %d dynamic trial(s) | %g%% length bins | connect-the-dots, NO fit\nthin = per trial, medium = per individual mean, thick black = group mean | RAW absolute force -- see _normalized companion for F/F0",
    dplyr::n_distinct(pooled$fish), n_iso, n_isv, n_dyn, STRAIN_BIN_PCT),
  y_lab = "Muscle force along u_hat (N, raw sign, no per-step correction -- see 2026-07-22 sign convention note)")

# NORMALIZED companion (added 2026-07-21, PI-directed -- see module header):
# same pooling, each point expressed as F / that trial+side's own L0 force.
n_norm_avail <- dplyr::n_distinct(dplyr::filter(pooled, is.finite(.data$force_emp_norm))$trial_id)
p_norm <- .build_fl_superplot(
  pooled, "force_emp_norm", "force_geom_norm",
  title = "Pooled Force-Length superplot, NORMALIZED (F / trial's own L0 force)",
  subtitle = sprintf(
    "Same pooling as the RAW superplot, each point divided by that TRIAL+SIDE's own mean L0 force (F/F0) to correct cross-trial force-scale differences (e.g. fatigue between repeat trials) before pooling | %d/%d trials contributed a usable F0 | %d individuals | %g%% length bins | connect-the-dots, NO fit\nPROTOTYPE (2026-07-21) -- compare against the RAW file before treating this as canonical",
    n_norm_avail, dplyr::n_distinct(pooled$trial_id), dplyr::n_distinct(pooled$fish), STRAIN_BIN_PCT),
  y_lab = "Muscle force / trial's own L0 force (F/F0, dimensionless)")

outfile <- file.path(OUTPUT_DIR, "FLsuperplot_isometric_isovelocity_pooled.png")
ggplot2::ggsave(outfile, p, width = 13, height = 6.5, dpi = 150)
cli::cli_alert_success("Wrote {outfile}")

outfile_norm <- file.path(OUTPUT_DIR, "FLsuperplot_isometric_isovelocity_pooled_normalized.png")
ggplot2::ggsave(outfile_norm, p_norm, width = 13, height = 6.5, dpi = 150)
cli::cli_alert_success("Wrote {outfile_norm}")

cli::cli_alert_info("Pooled points: {nrow(pooled)} rows | fish: {paste(sort(unique(pooled$fish)), collapse=', ')} | F0 available for {n_norm_avail}/{dplyr::n_distinct(pooled$trial_id)} trials")
