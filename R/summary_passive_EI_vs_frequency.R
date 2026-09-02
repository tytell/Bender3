# summary_passive_EI_vs_frequency.R
#
# Passive flexural stiffness (EI) vs. frequency for largemouth bass, 2025+2026.
# Two panels: (1) dynamic freq x curvature grids that are 3x3 or 4x4 within a
# file, dropping the isolated 0.5 Hz 3x3 step (those fish already have 1-7 Hz);
# (2) frequency_sweep chirps faceted by amplitude_frequency_exponent (alpha in
# theta ~ f^alpha). Passive cycles only (observed stim). No remounts / FAILED.
#
# Inertial correction (PI 2026-09-01): discard the 2026-07-23 empty-apparatus
# OLS (R^2 ~ 0). Use the 2026-07-06 dedicated geometry-sweep fit
# (2026-07-06_apparatus_inertia_calibration.json, built from
# 2026-07-06_apparatus_inertial_calibration/*.h5). I_app is evaluated at each
# trial's AOR and plate-to-plate span; I_spec is the lateral elliptical-
# cylinder term. Artifact correction_slope_sign is applied so the correction
# reduces (does not inflate) high-frequency apparent EI. 2025 Bender2 files
# are a different apparatus and stay uncorrected.
#
# Frequency / curvature binning (FIXED 2026-09-01, repeating the 2026-07
# pattern): do NOT plot on index_cycle_frequency_hertz / operating_point
# design labels. On bass16/17/18 combo files those labels lag the real
# motion by ~one step after the first amplitude block (same class of bug as
# R/analyze_bass18_multiaxial_0p5hz_20260728.R: 2.75 Hz motion labelled
# 5 Hz). Physical cycle windows (floor(t.norm)) already track the waveform;
# bin on the measured period (1/period_s) and measured curvature amplitude,
# then snap to the known design grid -- same as
# R/diag_rod40a_dynamic_ei_by_condition.R (zero-crossing frequency rounded
# to nearest known step). Design-table columns are kept in the CSV as
# freq_cmd / curv_cmd for QA only.
#
# Prefix summary_ (PI, 2026-09-01): pooled across individuals.
#   -> figs_summary/summary_EI_vs_frequency_passive_inertia_highconf.png
#      (July 6 I_app at the trial's own AOR/plate + analytic I_spec)
#   -> figs_summary/summary_EI_vs_frequency_passive_other.png
#      (2025 Bender2 uncorrected; 2026 files with defaulted clamp geometry)
#   -> data_processed/summary_EI_vs_frequency_passive_per_cycle.csv
#
# Run: Rscript R/summary_passive_EI_vs_frequency.R

suppressPackageStartupMessages({
  library(rhdf5); library(dplyr); library(tibble); library(purrr)
  library(ggplot2); library(patchwork); library(readr); library(cli)
})

.root <- if (nzchar(Sys.getenv("BENDER3_R_ROOT"))) Sys.getenv("BENDER3_R_ROOT") else "R"
src <- function(f) source(file.path(.root, f))
src("paths_config.R")
src("00_load_bender_flat.R")
src("03_analyze.R")
src("analyze_frequency_sweep.R")
src("apparatus_inertia_fit.R")
src("inertia_diag_common.R")

APPARATUS_FIT <- read_apparatus_inertia_fit(
  resolve_apparatus_fit_path("2026-07-06_apparatus_inertia_calibration.json")
)
# Tissue density used when the file has body width/height/clamp but density NA
# (same override as analyze_bass17_multiaxial_20260723.R / bass18 0.5 Hz).
DENSITY_DEFAULT_G_MM3 <- 0.00106
# Clamp geometry locked on later 2026 laterals; used only when a 2026
# CritterGripper file is missing aor / plate-to-plate attrs (bass14).
AOR_DEFAULT_MM <- 20.5
PLATE_DEFAULT_MM <- 55.0

DATA_DIR <- file.path(.crittergripper_root(), "02_processed", "data_processed")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGS_SUMMARY_DIR, recursive = TRUE, showWarnings = FALSE)

ARCHIVE <- .permanent_archive_root()
BENDER2_2025 <- file.path(dirname(ARCHIVE), "BenderData 2025-11-26")
# PermanentArchive root is .../bender_crittergripper; 2025 lives one level up.
if (!dir.exists(BENDER2_2025)) {
  BENDER2_2025 <- file.path(dirname(ARCHIVE), "01_PermanentArchive", "BenderData 2025-11-26")
}
if (!dir.exists(BENDER2_2025)) {
  BENDER2_2025 <- "/Users/yjimenez/Library/CloudStorage/OneDrive-ProvidenceCollege/01_JimenezLab/01_PermanentArchive/BenderData 2025-11-26"
}

`%||%` <- function(x, y) if (is.null(x) || (length(x) == 1L && is.na(x))) y else x

.top_keys <- function(path) {
  h5 <- H5Fopen(path, "H5F_ACC_RDONLY")
  on.exit(try(H5Fclose(h5), silent = TRUE), add = TRUE)
  h5ls(h5, recursive = FALSE)$name
}

.chr1 <- function(v, default = NA_character_) {
  if (is.null(v) || length(v) == 0L) return(default)
  x <- v[[1L]]
  if (is.raw(x) || is.character(x)) {
    s <- paste(as.character(x), collapse = "")
    if (!nzchar(s) || s %in% c("NA", "NaN")) default else s
  } else {
    as.character(x)[1L]
  }
}

.read_amp_freq_exponent <- function(path) {
  # amplitude_frequency_exponent (α): commanded amplitude scales as f^α.
  # α = 0 constant amplitude; α = -1 amplitude ∝ 1/f (peak angular velocity nearer constant).
  keys <- tryCatch(.top_keys(path), error = function(e) character(0))
  grab <- function(attrs, nms) {
    for (nm in nms) {
      if (!is.null(attrs[[nm]])) {
        v <- suppressWarnings(as.numeric(attrs[[nm]][1L]))
        if (is.finite(v)) return(v)
      }
    }
    NA_real_
  }
  if ("NominalStimulus" %in% keys) {
    a <- tryCatch(h5readAttributes(path, "/NominalStimulus"), error = function(e) list())
    return(grab(a, "AmplitudeFrequencyExponent"))
  }
  if ("metadata" %in% keys) {
    a <- tryCatch(h5readAttributes(path, "/metadata"), error = function(e) list())
    return(grab(a, c("protocol_amplitude_frequency_exponent", "amplitude_frequency_exponent")))
  }
  if ("01_Metadata" %in% keys) {
    for (g in c("/01_Metadata/protocol_metadata", "/01_Metadata/bender_settings", "/01_Metadata")) {
      a <- tryCatch(h5readAttributes(path, g), error = function(e) NULL)
      if (is.null(a)) next
      v <- grab(a, c("amplitude_frequency_exponent", "protocol_amplitude_frequency_exponent"))
      if (is.finite(v)) return(v)
    }
  }
  NA_real_
}

.snap_to_grid <- function(x, grid, abs_tol = 0.15) {
  vapply(x, function(v) {
    if (!is.finite(v) || length(grid) == 0L) return(NA_real_)
    d <- abs(grid - v)
    j <- which.min(d)
    if (d[j] <= abs_tol) grid[j] else NA_real_
  }, numeric(1L))
}

.dbl1_attr <- function(a, k) {
  v <- suppressWarnings(as.numeric(a[[k]][1L]))
  if (length(v) == 0L) NA_real_ else v
}

.read_inertia_geom <- function(path) {
  keys <- tryCatch(.top_keys(path), error = function(e) character(0))
  grab <- function(grp) tryCatch(h5readAttributes(path, grp), error = function(e) list())
  a <- list()
  if ("metadata" %in% keys) {
    a <- grab("/metadata")
  } else if ("01_Metadata" %in% keys) {
    a <- c(grab("/01_Metadata"), grab("/01_Metadata/bender_settings"),
           grab("/01_Metadata/inertial_calibration_profile"))
  }
  list(
    width_mm  = .dbl1_attr(a, "measurement_specimen_local_body_width_millimeter"),
    height_mm = .dbl1_attr(a, "measurement_specimen_local_body_height_millimeter"),
    clamp_mm  = .dbl1_attr(a, "measurement_clamp_separation_millimeter"),
    density_g_mm3 = .dbl1_attr(a, "measurement_specimen_density_gram_per_cubic_millimeter"),
    aor_mm    = .dbl1_attr(a, "calibration_inertia_apparatus_aor_to_clamp_millimeter"),
    plate_mm  = .dbl1_attr(a, "calibration_inertia_apparatus_plate_to_plate_millimeter")
  )
}

.specimen_inertia_lateral_kg_m2 <- function(width_mm, height_mm, clamp_mm, density_g_mm3) {
  bY <- width_mm / 2.0; bZ <- height_mm / 2.0; L <- clamp_mm; rho <- density_g_mm3
  if (!all(is.finite(c(bY, bZ, L, rho))) || rho <= 0) return(NA_real_)
  m_g <- rho * pi * bY * bZ * L
  (m_g * (bZ^2 / 4 + L^2 / 12)) * kg_m2_per_g_mm2
}

# July 6 F4 fit at this trial's clamp geometry + analytic fish MOI.
# I_app_signed uses artifact correction_slope_sign (empty-apparatus tau vs
# alpha was negative). tau_corr = tau - I_total_signed * alpha_rad.
.apply_july06_inertia <- function(td, path) {
  td$I_app_kgm2 <- NA_real_
  td$I_spec_kgm2 <- NA_real_
  td$I_total_signed_kgm2 <- NA_real_
  td$inertia_applied <- FALSE
  td$inertia_geom_defaulted <- FALSE
  if (is.null(td) || nrow(td) == 0L) return(td)
  is_cg <- grepl("bender_crittergripper", path, ignore.case = TRUE)
  if (!is_cg) {
    cli::cli_alert_info("{basename(path)}: not CritterGripper -- inertia left uncorrected (2026-07-06 fit is rig-specific)")
    return(td)
  }
  geom <- .read_inertia_geom(path)
  aor <- geom$aor_mm
  plate <- geom$plate_mm
  used_default_geom <- FALSE
  if (!is.finite(aor) || !is.finite(plate)) {
    aor <- AOR_DEFAULT_MM
    plate <- PLATE_DEFAULT_MM
    used_default_geom <- TRUE
  }
  app <- apparatus_inertia_from_fit(APPARATUS_FIT, aor, plate)
  if (!isTRUE(app$in_domain)) {
    cli::cli_alert_warning("{basename(path)}: apparatus fit EXTRAPOLATED (aor={aor} mm, plate={plate} mm)")
  }
  dens <- geom$density_g_mm3
  if (!is.finite(dens)) dens <- DENSITY_DEFAULT_G_MM3
  I_spec <- .specimen_inertia_lateral_kg_m2(geom$width_mm, geom$height_mm, geom$clamp_mm, dens)
  if (!is.finite(I_spec)) I_spec <- 0
  I_app_signed <- app$i_gmm2 * kg_m2_per_g_mm2 * APPARATUS_FIT$correction_slope_sign
  I_total_signed <- I_app_signed + sign(I_app_signed) * I_spec

  ang <- if ("enc.deg" %in% names(td) && any(is.finite(td$enc.deg))) td$enc.deg else td$angle.deg
  alpha <- .kinematics(td$t.s, ang)$alpha_rad_s2
  n <- min(nrow(td), length(alpha))
  td$tau_Nm[seq_len(n)] <- td$tau_Nm[seq_len(n)] - I_total_signed * alpha[seq_len(n)]
  td$I_app_kgm2 <- I_app_signed
  td$I_spec_kgm2 <- I_spec
  td$I_total_signed_kgm2 <- I_total_signed
  td$inertia_applied <- TRUE
  td$inertia_geom_defaulted <- used_default_geom
  extra <- if (used_default_geom) " [aor/plate defaulted 20.5/55 mm]" else ""
  cli::cli_alert_success(
    "{basename(path)}: I_app={signif(I_app_signed, 4)} kg*m^2, I_spec={signif(I_spec, 4)} kg*m^2{extra}"
  )
  td
}

.tau_col <- function(td) {
  hit <- intersect(c("tau_Nm", "ztorque.Nm", "ztorque0.Nm", "xtorque.Nm",
                     "xtorque0.Nm", "torque.Nm"), names(td))
  if (length(hit) == 0L) NA_character_ else hit[[1L]]
}

.finish_td <- function(td, dclamp_m, fish, path, protocol) {
  if (is.null(td) || nrow(td) == 0L) return(NULL)
  if (!"curve.invm" %in% names(td) && is.finite(dclamp_m) && "angle.deg" %in% names(td)) {
    td$curve.invm <- td$angle.deg * pi / 180.0 / dclamp_m
  }
  if ("curve.invm" %in% names(td) && !"curverate.invms" %in% names(td)) {
    td$curverate.invms <- (dplyr::lead(td$curve.invm) - dplyr::lag(td$curve.invm)) /
      (dplyr::lead(td$t.s) - dplyr::lag(td$t.s))
  }
  if (!"cycle" %in% names(td) && "t.norm" %in% names(td)) {
    td$cycle <- floor(td$t.norm)
  }
  if (!"stim" %in% names(td)) {
    if ("stim_side" %in% names(td)) {
      ss <- tolower(trimws(as.character(td$stim_side)))
      td$stim <- factor(dplyr::case_match(ss, "left" ~ "L", "right" ~ "R",
                                          "both" ~ "B", .default = "0"),
                        levels = c("0", "L", "R", "B"))
    } else if ("stim.V" %in% names(td)) {
      td$stim <- factor(ifelse(is.finite(td$stim.V) & td$stim.V > 0.1, "L", "0"),
                        levels = c("0", "L", "R", "B"))
    } else {
      td$stim <- factor("0", levels = c("0", "L", "R", "B"))
    }
  }
  td$fishcode <- fish
  td$fullpathname <- path
  td$filename <- basename(path)
  td$protocol <- protocol
  td$dclamp.m <- dclamp_m
  if (!"freq.Hz" %in% names(td)) td$freq.Hz <- NA_real_
  if (!"amp.deg" %in% names(td)) td$amp.deg <- NA_real_
  if (!"curvature.invm" %in% names(td)) td$curvature.invm <- NA_real_
  td
}

.load_flat <- function(path, protocol) {
  td <- load_bender_flat(path, do_filter = TRUE, loadtorques = "z", loadforces = FALSE)
  if (is.null(td) || nrow(td) == 0L) return(NULL)
  fish <- unique(stats::na.omit(td$fishcode))[1L]
  if (length(fish) == 0L || is.na(fish)) {
    m <- regmatches(basename(path), regexpr("bass[0-9]+", basename(path), ignore.case = TRUE))
    fish <- if (length(m)) tolower(m) else "unknown"
  }
  dclamp <- unique(td$dclamp.m)[1L]
  tq <- .tau_col(td)
  if (is.na(tq)) return(NULL)
  td$tau_Nm <- td[[tq]]
  .finish_td(td, dclamp, fish, path, protocol)
}

.load_v1_01_metadata <- function(path, protocol) {
  h5 <- H5Fopen(path, "H5F_ACC_RDONLY")
  on.exit(try(H5Fclose(h5), silent = TRUE), add = TRUE)
  a <- tryCatch(h5readAttributes(h5, "/01_Metadata"), error = function(e) list())
  bs <- tryCatch(h5readAttributes(h5, "/01_Metadata/bender_settings"), error = function(e) list())
  fish <- .chr1(a[["specimen_id"]], NA_character_)
  if (is.na(fish)) {
    m <- regmatches(basename(path), regexpr("bass[0-9]+", basename(path), ignore.case = TRUE))
    fish <- if (length(m)) tolower(m) else "unknown"
  }
  dclamp_mm <- suppressWarnings(as.numeric((a[["measurement_clamp_separation_millimeter"]] %||%
                                              bs[["dclamp"]] %||% bs[["test_segment_length_mm"]])[1L]))
  dclamp_m <- if (is.finite(dclamp_mm)) dclamp_mm / 1000.0 else NA_real_

  ts_info <- h5ls(h5, recursive = TRUE)
  trial <- sort(ts_info$name[ts_info$group == "/02_TimeSeries" & grepl("^trial_", ts_info$name)])[1L]
  if (is.na(trial) || !nzchar(trial)) return(NULL)
  base <- paste0("/02_TimeSeries/", trial, "/")
  rd <- function(nm) {
    p <- paste0(base, nm)
    if (H5Lexists(h5, p)) tryCatch(h5read(h5, p), error = function(e) NULL) else NULL
  }
  t_s <- rd("time_second"); if (is.null(t_s)) t_s <- rd("t")
  if (is.null(t_s)) return(NULL)
  t_norm <- rd("time_normalized"); if (is.null(t_norm)) t_norm <- rd("tnorm")
  ang <- rd("angle_commanded_degree"); if (is.null(ang)) ang <- rd("angle_cmd")
  stim_side <- rd("stim_side")
  tau <- rd("primary_torque_raw")
  if (is.null(tau)) {
    ft <- rd("forcetorque_raw")
    if (is.null(ft)) ft <- rd("forcetorque")
    if (!is.null(ft)) {
      if (!is.matrix(ft)) ft <- as.matrix(ft)
      if (nrow(ft) == 6L && ncol(ft) != 6L) ft <- t(ft)
      tau <- ft[, 6L]
    }
  }
  if (is.null(tau) || is.null(ang)) return(NULL)
  N <- min(length(t_s), length(ang), length(tau))
  td <- tibble::tibble(
    t.s = as.numeric(t_s)[seq_len(N)],
    t.norm = if (is.null(t_norm)) NA_real_ else as.numeric(t_norm)[seq_len(N)],
    angle.deg = as.numeric(ang)[seq_len(N)],
    tau_Nm = as.numeric(tau)[seq_len(N)],
    stim_side = if (is.null(stim_side)) "none" else as.character(stim_side)[seq_len(N)]
  )
  .finish_td(td, dclamp_m, fish, path, protocol)
}

.load_bender2 <- function(path, protocol) {
  h5 <- H5Fopen(path, "H5F_ACC_RDONLY")
  on.exit(try(H5Fclose(h5), silent = TRUE), add = TRUE)
  root_a <- h5readAttributes(h5, "/")
  ns_a <- tryCatch(h5readAttributes(h5, "/NominalStimulus"), error = function(e) list())
  fish <- .chr1(root_a[["FishCode"]], "unknown")
  dclamp_mm <- suppressWarnings(as.numeric(root_a[["ClampDistance_mm"]][1L]))
  dclamp_m <- if (is.finite(dclamp_mm)) dclamp_mm / 1000.0 else NA_real_
  t_s <- as.numeric(h5read(h5, "/NominalStimulus/t"))
  t_norm <- as.numeric(h5read(h5, "/NominalStimulus/tnorm"))
  ang <- tryCatch(as.numeric(h5read(h5, "/NominalStimulus/Position")), error = function(e) NULL)
  if (is.null(ang)) {
    ang <- tryCatch(as.numeric(h5read(h5, "/Calibrated/Encoder")), error = function(e) NULL)
  }
  if (is.null(ang)) ang <- tryCatch(as.numeric(h5read(h5, "/RawInput/Encoder")), error = function(e) NULL)
  tau <- as.numeric(h5read(h5, "/Calibrated/xTorque"))
  stim.V <- tryCatch(as.numeric(h5read(h5, "/RawInput/Left stim")), error = function(e) NULL)
  if (is.null(stim.V)) stim.V <- tryCatch(as.numeric(h5read(h5, "/RawInput/activation_monitor")), error = function(e) NULL)
  N <- min(length(t_s), length(ang), length(tau))
  # Frequency-sweep t.norm is in radians in this schema (load_bender_data.R).
  ns_type <- .chr1(ns_a[["Type"]], "")
  if (identical(ns_type, "Frequency Sweep")) {
    t_norm <- t_norm / (2 * pi)
  }
  td <- tibble::tibble(
    t.s = t_s[seq_len(N)],
    t.norm = t_norm[seq_len(N)],
    angle.deg = ang[seq_len(N)],
    tau_Nm = tau[seq_len(N)],
    stim.V = if (is.null(stim.V)) 0 else stim.V[seq_len(N)]
  )
  # Bias-subtract like the 2026 loader (pre-t=0 window).
  pre <- is.finite(td$t.s) & td$t.s < -0.1
  if (any(pre)) {
    b <- mean(td$tau_Nm[pre], na.rm = TRUE)
    if (is.finite(b)) td$tau_Nm <- td$tau_Nm - b
  }
  td <- .finish_td(td, dclamp_m, fish, path, protocol)
  # Attach per-cycle design when present (Frequency Amplitude Combo).
  fbc <- ns_a[["FrequencyByCycle"]]
  abc <- ns_a[["AmplitudeByCycle"]]
  if (!is.null(fbc) && length(fbc) > 1L) {
    fbc <- as.numeric(fbc); abc <- as.numeric(abc)
    ci <- td$cycle
    ok <- is.finite(ci) & ci >= 0 & (ci + 1L) <= length(fbc)
    td$freq.Hz <- NA_real_
    td$amp.deg <- NA_real_
    td$freq.Hz[ok] <- fbc[ci[ok] + 1L]
    if (length(abc) == length(fbc)) td$amp.deg[ok] <- abc[ci[ok] + 1L]
    if (is.finite(dclamp_m) && "amp.deg" %in% names(td)) {
      td$curvature.invm <- td$amp.deg * pi / 180.0 / dclamp_m
    }
  }
  td
}

.load_any <- function(path, protocol) {
  keys <- tryCatch(.top_keys(path), error = function(e) character(0))
  td <- NULL
  if ("metadata" %in% keys && "timeseries" %in% keys) {
    td <- .load_flat(path, protocol)
  } else if ("01_Metadata" %in% keys) {
    td <- .load_v1_01_metadata(path, protocol)
  } else if ("Calibrated" %in% keys && "NominalStimulus" %in% keys) {
    td <- .load_bender2(path, protocol)
  } else {
    cli::cli_alert_warning("unrecognized schema: {basename(path)} [{paste(keys, collapse=', ')}]")
    return(NULL)
  }
  if (is.null(td) || !"tau_Nm" %in% names(td)) return(td)
  .apply_july06_inertia(td, path)
}

.ei_dynamic_cycles <- function(td) {
  tq <- "tau_Nm"
  td2 <- set_cycle_types(td) |>
    dplyr::filter(.data$cycletype == "pass", is.finite(.data$cycle), .data$cycle >= 0,
                  is.finite(.data[[tq]]), is.finite(.data$curve.invm))
  if (nrow(td2) < 10L) return(tibble::tibble())
  td2 |>
    dplyr::group_by(.data$fishcode, .data$filename, .data$cycle) |>
    dplyr::summarise(
      n = dplyr::n(),
      freq_cmd = dplyr::first(.data$freq.Hz),
      curv_cmd = dplyr::first(.data$curvature.invm),
      amp_deg = dplyr::first(.data$amp.deg),
      period_s = diff(range(.data$t.s, na.rm = TRUE)),
      curv_amp = (max(.data$curve.invm, na.rm = TRUE) - min(.data$curve.invm, na.rm = TRUE)) / 2,
      EI_Nm2 = tryCatch(calc_avg_stiffness(.data$curve.invm, .data[[tq]]), error = function(e) NA_real_),
      I_app_kgm2 = dplyr::first(.data$I_app_kgm2),
      I_spec_kgm2 = dplyr::first(.data$I_spec_kgm2),
      inertia_applied = dplyr::first(.data$inertia_applied),
      inertia_geom_defaulted = dplyr::first(.data$inertia_geom_defaulted),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      freq_meas = dplyr::if_else(.data$period_s > 0, 1.0 / .data$period_s, NA_real_),
      # Measured period is authoritative; design-table freq_cmd lags real
      # motion by ~one step on bass16/17/18 combo files (see header).
      freq.Hz = dplyr::coalesce(.data$freq_meas, .data$freq_cmd)
    )
}

.ei_sweep_cycles <- function(td) {
  td2 <- set_cycle_types(td)
  # Chirps are passive-only; still drop any cycle window that contains stim.
  cyc <- analyze_frequency_sweep(td2, torque_col = "tau_Nm")
  if (nrow(cyc) == 0L) return(tibble::tibble())
  pass_ok <- vapply(seq_len(nrow(cyc)), function(i) {
    win <- td2$t.s >= cyc$t_start[i] & td2$t.s < cyc$t_end[i]
    all(as.character(td2$stim[win]) == "0", na.rm = TRUE)
  }, logical(1L))
  cyc <- cyc[pass_ok, , drop = FALSE]
  if (nrow(cyc) == 0L) return(tibble::tibble())
  cyc |>
    dplyr::transmute(
      fishcode = td2$fishcode[1L], filename = td2$filename[1L],
      cycle = .data$cycle_index, n = NA_integer_,
      freq_cmd = NA_real_, curv_cmd = NA_real_, amp_deg = NA_real_,
      period_s = .data$t_end - .data$t_start,
      curv_amp = .data$mean_curvature_amplitude_invm,
      EI_Nm2 = .data$EI1.Nm2, freq_meas = .data$frequency_hz,
      freq.Hz = .data$frequency_hz, status = .data$status,
      amp_freq_exp = if ("amp_freq_exp" %in% names(td2)) td2$amp_freq_exp[1L] else NA_real_,
      I_app_kgm2 = if ("I_app_kgm2" %in% names(td2)) td2$I_app_kgm2[1L] else NA_real_,
      I_spec_kgm2 = if ("I_spec_kgm2" %in% names(td2)) td2$I_spec_kgm2[1L] else NA_real_,
      inertia_applied = if ("inertia_applied" %in% names(td2)) isTRUE(td2$inertia_applied[1L]) else FALSE,
      inertia_geom_defaulted = if ("inertia_geom_defaulted" %in% names(td2)) isTRUE(td2$inertia_geom_defaulted[1L]) else FALSE
    )
}

# =============================================================================
# Manifest (main sessions only; no remount / FAILED)
# =============================================================================
dyn_files <- c(
  file.path(BENDER2_2025, "06122025_bass01", c("bass01_002.h5", "bass01_003.h5", "bass01_004.h5")),
  file.path(ARCHIVE, "2026-06-15_bass14_bender", "2026-06-15_bass14_bender_13_dynamic.h5"),
  file.path(ARCHIVE, "2026-07-14_bass16_bender",
            paste0("2026-07-14_bass16_bender_", c("17", "18", "19"), "_dynamic.h5")),
  file.path(ARCHIVE, "2026-07-15_bass17_bender",
            paste0("2026-07-15_bass17_bender_", c("17", "18", "19", "20"), "_dynamic.h5")),
  file.path(ARCHIVE, "2026-07-16_bass18_bender",
            paste0("2026-07-16_bass18_bender_", c("14", "15"), "_dynamic.h5"))
)

sweep_dirs <- c(
  file.path(BENDER2_2025, c("2025-07-08_bass06", "2025-07-09_bass07", "2025-07-09_bass08",
                            "2025-07-17_bass09", "2025-07-17_bass10", "2025-07-21_bass11",
                            "2025-07-21_bass12")),
  file.path(ARCHIVE, c("2026-06-04_bass13_bender", "2026-07-14_bass16_bender",
                       "2026-07-15_bass17_bender", "2026-07-16_bass18_bender"))
)
sweep_files <- unique(unlist(lapply(sweep_dirs, function(d) {
  if (!dir.exists(d)) return(character(0))
  c(list.files(d, pattern = "frequency_sweep\\.h5$", full.names = TRUE),
    list.files(d, pattern = "_sweep\\.h5$", full.names = TRUE))
})))
# 2025 chirps are numbered, not suffix-tagged -- keep those whose Type is Frequency Sweep.
sweep_2025_extra <- unlist(lapply(
  file.path(BENDER2_2025, c("2025-07-08_bass06", "2025-07-09_bass07", "2025-07-09_bass08",
                            "2025-07-17_bass09", "2025-07-17_bass10", "2025-07-21_bass11",
                            "2025-07-21_bass12")),
  function(d) list.files(d, pattern = "\\.h5$", full.names = TRUE)
))
is_freq_sweep <- function(p) {
  a <- tryCatch(h5readAttributes(p, "/NominalStimulus"), error = function(e) NULL)
  identical(.chr1(a[["Type"]], ""), "Frequency Sweep")
}
sweep_files <- unique(c(sweep_files, sweep_2025_extra[vapply(sweep_2025_extra, is_freq_sweep, logical(1L))]))
sweep_files <- sweep_files[file.exists(sweep_files)]
dyn_files <- dyn_files[file.exists(dyn_files)]

cli::cli_h1("Passive EI vs frequency: {length(dyn_files)} grid files, {length(sweep_files)} sweep files")

# 0.5 Hz is a real 3x3 step (bass16/17/18) but sits as an isolated vertical
# cluster below the 1-7 Hz 4x4; those fish already have the richer 4x4 range.
GRID_FREQ <- c(1, 2, 2.75, 3, 5, 7)
GRID_CURV <- c(2, 3.36, 4, 4.47, 5.59, 6, 8)

process_dynamic <- function(path) {
  cli::cli_alert_info("dynamic {basename(path)}")
  td <- tryCatch(.load_any(path, "dynamic"), error = function(e) {
    cli::cli_alert_danger("{basename(path)}: {conditionMessage(e)}")
    NULL
  })
  if (is.null(td)) return(tibble::tibble())
  rows <- .ei_dynamic_cycles(td)
  if (nrow(rows) == 0L) return(tibble::tibble())
  rows |>
    dplyr::mutate(
      panel = "dynamic", amp_freq_exp = NA_real_,
      freq_nominal = .snap_to_grid(.data$freq.Hz, GRID_FREQ, abs_tol = 0.12),
      # Same lag hits operating_point / amp labels; snap measured p2p/2.
      curv_nominal = .snap_to_grid(dplyr::coalesce(.data$curv_amp, .data$curv_cmd), GRID_CURV, abs_tol = 0.25)
    ) |>
    dplyr::filter(is.finite(.data$EI_Nm2), .data$EI_Nm2 != 0,
                  is.finite(.data$freq_nominal), is.finite(.data$curv_nominal),
                  .data$freq_nominal >= 1,
                  .data$n >= 10L)
}

process_sweep <- function(path) {
  cli::cli_alert_info("sweep {basename(path)}")
  td <- tryCatch(.load_any(path, "frequency_sweep"), error = function(e) {
    cli::cli_alert_danger("{basename(path)}: {conditionMessage(e)}")
    NULL
  })
  if (is.null(td) || !"tau_Nm" %in% names(td)) return(tibble::tibble())
  rows <- tryCatch(.ei_sweep_cycles(td), error = function(e) {
    cli::cli_alert_danger("{basename(path)} sweep EI: {conditionMessage(e)}")
    tibble::tibble()
  })
  if (nrow(rows) == 0L) return(tibble::tibble())
  expn <- .read_amp_freq_exponent(path)
  rows |>
    dplyr::mutate(panel = "sweep", freq_nominal = NA_real_, curv_nominal = NA_real_,
                  amp_freq_exp = expn) |>
    dplyr::filter(is.finite(.data$EI_Nm2), .data$EI_Nm2 != 0, is.finite(.data$freq.Hz))
}

dyn_rows <- purrr::map_dfr(dyn_files, process_dynamic)
swp_rows <- purrr::map_dfr(sweep_files, process_sweep)

cli::cli_h2("Dynamic cycles kept: {nrow(dyn_rows)} across {dplyr::n_distinct(dyn_rows$fishcode)} fish")
if (nrow(dyn_rows) > 0L) {
  dyn_rows <- dyn_rows |>
    dplyr::group_by(.data$fishcode, .data$freq_nominal, .data$curv_nominal) |>
    dplyr::filter(dplyr::n() >= 3L) |>
    dplyr::ungroup()
  cli::cli_alert_info("After dropping singleton transition cells: {nrow(dyn_rows)} cycles")
  n_mis <- sum(is.finite(dyn_rows$freq_cmd) & is.finite(dyn_rows$freq_meas) &
                 abs(dyn_rows$freq_cmd - dyn_rows$freq_meas) > 0.3, na.rm = TRUE)
  cli::cli_alert_info("Design-label mismatches (|freq_cmd - freq_meas| > 0.3 Hz): {n_mis}/{nrow(dyn_rows)} cycles (rebinned onto measured freq/curvature)")
  print(as.data.frame(dyn_rows |> dplyr::count(.data$fishcode, .data$freq_nominal, .data$curv_nominal) |> dplyr::arrange(.data$fishcode)), row.names = FALSE)
}
cli::cli_h2("Sweep cycles kept: {nrow(swp_rows)} across {dplyr::n_distinct(swp_rows$fishcode)} fish")
if (nrow(swp_rows) > 0L) print(as.data.frame(swp_rows |> dplyr::count(.data$fishcode)), row.names = FALSE)

.mark_inertia_conf <- function(d) {
  if (nrow(d) == 0L) return(d)
  d |>
    dplyr::mutate(
      inertia_high_conf = .data$inertia_applied %in% TRUE &
        !(.data$inertia_geom_defaulted %in% TRUE) &
        is.finite(.data$I_spec_kgm2) & .data$I_spec_kgm2 > 0
    )
}
dyn_rows <- .mark_inertia_conf(dyn_rows)
swp_rows <- .mark_inertia_conf(swp_rows)

all_rows <- dplyr::bind_rows(dyn_rows, swp_rows)
csv_out <- file.path(DATA_DIR, "summary_EI_vs_frequency_passive_per_cycle.csv")
readr::write_csv(all_rows, csv_out)
cli::cli_alert_success("Wrote {csv_out}")

# =============================================================================
# Figures (high-confidence inertia correction vs everything else)
# =============================================================================
EXP_LABS <- c(
  `0` = "alpha = 0  (constant amplitude)",
  `-0.5` = "alpha = -0.5",
  `-1` = "alpha = -1  (amplitude ~ 1/f)"
)

.fish_palette <- function(fish) {
  fish <- sort(unique(as.character(fish)))
  pal <- grDevices::colorRampPalette(
    c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#1f78b4")
  )(max(length(fish), 1L))
  names(pal) <- fish
  pal
}

.ei_freq_figure <- function(dyn, swp, dyn_title, swp_title, swp_sub) {
  fish_lvls <- sort(unique(c(as.character(dyn$fishcode), as.character(swp$fishcode))))
  pal <- .fish_palette(fish_lvls)
  if (nrow(dyn) > 0L) {
    dyn$fishcode <- factor(dyn$fishcode, levels = fish_lvls)
    p_dyn <- ggplot(dyn, aes(.data$freq_nominal, abs(.data$EI_Nm2), color = .data$fishcode)) +
      geom_point(size = 1.6, alpha = 0.35, position = position_jitter(width = 0.08, height = 0)) +
      stat_summary(fun = mean, geom = "line", aes(group = .data$fishcode), linewidth = 0.7) +
      stat_summary(fun = mean, geom = "point", size = 2.4, shape = 21, fill = "white") +
      facet_wrap(~ paste0("kappa = ", curv_nominal, " /m"), nrow = 2) +
      scale_color_manual(values = pal, name = "Specimen", drop = FALSE) +
      labs(title = dyn_title, x = "Frequency (Hz)", y = "|EI| (N m^2)") +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom", panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank())
  } else {
    p_dyn <- ggplot() + theme_void() + labs(title = dyn_title, subtitle = "no dynamic cycles")
  }
  if (nrow(swp) > 0L) {
    swp <- swp |>
      dplyr::mutate(
        fishcode = factor(.data$fishcode, levels = fish_lvls),
        exp_lbl = factor(
          as.character(.data$amp_freq_exp),
          levels = c("0", "-0.5", "-1"),
          labels = unname(EXP_LABS)
        )
      )
    p_swp <- ggplot(dplyr::filter(swp, !is.na(.data$exp_lbl)),
                    aes(.data$freq.Hz, abs(.data$EI_Nm2), color = .data$fishcode)) +
      geom_point(size = 0.8, alpha = 0.22) +
      geom_smooth(se = FALSE, linewidth = 0.8, method = "loess", span = 0.35, formula = y ~ x) +
      facet_wrap(~ exp_lbl, ncol = 3, scales = "free_y") +
      scale_color_manual(values = pal, name = "Specimen", drop = FALSE) +
      labs(title = swp_title, subtitle = swp_sub,
           x = "Instantaneous frequency (Hz)", y = "|EI| (N m^2)") +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom", panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            strip.text = element_text(size = 10))
  } else {
    p_swp <- ggplot() + theme_void() + labs(title = swp_title, subtitle = "no sweep cycles")
  }
  (p_dyn / p_swp) +
    patchwork::plot_layout(guides = "collect", heights = c(1.1, 1)) +
    patchwork::plot_annotation(
      tag_levels = "A",
      theme = theme(plot.tag = element_text(size = 12, face = "bold"))
    ) &
    theme(legend.position = "bottom")
}

dyn_hi <- dplyr::filter(dyn_rows, .data$inertia_high_conf)
swp_hi <- dplyr::filter(swp_rows, .data$inertia_high_conf)
dyn_lo <- dplyr::filter(dyn_rows, !.data$inertia_high_conf)
swp_lo <- dplyr::filter(swp_rows, !.data$inertia_high_conf)

cli::cli_h2("High-confidence inertia: {nrow(dyn_hi)} dynamic + {nrow(swp_hi)} sweep cycles ({paste(sort(unique(c(dyn_hi$fishcode, swp_hi$fishcode))), collapse=', ')})")
cli::cli_h2("Other (uncorrected or defaulted geometry): {nrow(dyn_lo)} dynamic + {nrow(swp_lo)} sweep cycles ({paste(sort(unique(c(dyn_lo$fishcode, swp_lo$fishcode))), collapse=', ')})")

fig_hi <- .ei_freq_figure(
  dyn_hi, swp_hi,
  dyn_title = "Dynamic (July 6 I_app at trial geometry + analytic I_spec)",
  swp_title = "Frequency sweep, split by amplitude-frequency exponent",
  swp_sub = "High-confidence inertial correction only (bass16/17/18)"
)
fig_lo <- .ei_freq_figure(
  dyn_lo, swp_lo,
  dyn_title = "Dynamic (no high-confidence inertial correction)",
  swp_title = "Frequency sweep, split by amplitude-frequency exponent",
  swp_sub = "2025 Bender2 uncorrected; 2026 files with defaulted AOR/plate and no I_spec"
)

fig_hi_out <- file.path(FIGS_SUMMARY_DIR, "summary_EI_vs_frequency_passive_inertia_highconf.png")
fig_lo_out <- file.path(FIGS_SUMMARY_DIR, "summary_EI_vs_frequency_passive_other.png")
ggplot2::ggsave(fig_hi_out, fig_hi, width = 12, height = 11, dpi = 150)
ggplot2::ggsave(fig_lo_out, fig_lo, width = 12, height = 11, dpi = 150)
cli::cli_alert_success("Wrote {fig_hi_out}")
cli::cli_alert_success("Wrote {fig_lo_out}")
cli::cli_h1("summary_passive_EI_vs_frequency.R complete")
