# 02_deconvolve.R
# Inertial deconvolution: subtract the modeled apparatus + specimen inertia
# contribution from the calibrated torque time series.
#
# tau_corrected = tau_calibrated - (I_total * alpha + bias)
#
# where:
#   tau_calibrated  calibrated torque on the bending axis (N*m), from
#                   derived/forcetorque_calibrated in the raw file or hub
#   I_total         sum of apparatus and specimen MOI, converted from
#                   g*mm^2 to N*m/(deg/s^2)  via * 1e-9 * (pi/180)
#   alpha           angular acceleration (deg/s^2) from the angle time series
#   bias            calibration_inertia_bias_newton_meter (constant offset)
#
# HDF5 layout note: all calibration_inertia_* values are stored as GROUP
# ATTRIBUTES of /metadata, not as datasets.
#
# Architecture (D11, D12):
#   - Reads calibration_inertia_* from /metadata attrs in the raw file.
#   - Reads calibrated torque from derived/forcetorque_calibrated.
#   - Writes tau_corrected to hub derived/torque_inertia_corrected.
#   - NEVER modifies the raw archive.
#
# Usage:
#   source("R/02_deconvolve.R")
#   deconvolve_bender(raw_path = "path/to/raw.h5", hub_path = "path/to/hub.h5")

library(rhdf5)
library(dplyr)

deconvolve_bender <- function(
    raw_path,
    hub_path       = NULL,
    return_series  = TRUE,
    verbose        = TRUE
) {
  result <- NULL

  tryCatch({
    h5r <- H5Fopen(raw_path, "H5F_ACC_RDONLY")
    on.exit(try(H5Fclose(h5r), silent = TRUE), add = TRUE)

    # Attribute accessor (all calibration_inertia_* and daq_* scalars are attrs)
    m_attrs <- tryCatch(h5readAttributes(h5r, "/metadata"), error = function(e) list())
    m_a  <- function(key) m_attrs[[key]]
    m_ds <- function(key) {
      path <- paste0("/metadata/", key)
      if (H5Lexists(h5r, path)) tryCatch(h5read(h5r, path), error = function(e) NULL) else NULL
    }

    inertia <- .dcv_read_inertia(m_a, verbose)
    tau_raw <- .dcv_read_calibrated_torque(h5r, m_a, inertia$axis_sensor, hub_path, verbose)
    if (is.null(tau_raw)) stop("calibrated torque unavailable; run 01_calibrate.R first")

    angle_ds <- .dcv_read_angle(h5r, verbose)
    if (is.null(angle_ds)) stop("angle time series unavailable")

    # Angular acceleration via central differences (deg/s^2)
    t_s   <- angle_ds$t_s
    ang   <- angle_ds$angle_deg
    dt    <- c(diff(t_s)[1L], diff(t_s))
    vel   <- c(diff(ang)[1L] / dt[1L], diff(ang) / diff(t_s))
    alpha_deg_s2 <- c(NA_real_, diff(vel) / diff(t_s))

    N          <- min(length(tau_raw), length(alpha_deg_s2))
    tau_raw    <- tau_raw[seq_len(N)]
    alpha_use  <- alpha_deg_s2[seq_len(N)]

    # g*mm^2 -> N*m/(deg/s^2): multiply by 1e-9 (kg*m^2) * (pi/180) (rad/deg)
    g_mm2_to_Nm_degss <- 1e-9 * (pi / 180.0)
    I_total <- 0.0
    if (is.finite(inertia$apparatus_moi_g_mm2)) {
      I_total <- I_total + inertia$apparatus_moi_g_mm2 * g_mm2_to_Nm_degss
    }
    if (is.finite(inertia$specimen_moi_g_mm2)) {
      I_total <- I_total + inertia$specimen_moi_g_mm2 * g_mm2_to_Nm_degss
    }
    bias_Nm <- if (is.finite(inertia$bias_newton_meter)) inertia$bias_newton_meter else 0.0

    if (I_total == 0.0) {
      cli::cli_warn(
        "deconvolve_bender: no MOI available (apparatus + specimen both absent/NA); correction = bias only"
      )
    }

    tau_corrected <- tau_raw - (I_total * alpha_use + bias_Nm)

    if (verbose) {
      cli::cli_alert_success(
        "Inertial correction applied: I_total = {signif(I_total, 4)} N*m/(deg/s^2), bias = {signif(bias_Nm, 4)} N*m"
      )
    }

    result <- tau_corrected
    if (!is.null(hub_path)) .dcv_write_hub(hub_path, tau_corrected, verbose)

  }, error = function(e) {
    cli::cli_alert_danger("deconvolve_bender: {conditionMessage(e)}")
  })

  if (isTRUE(return_series)) invisible(result) else invisible(NULL)
}


# =============================================================================
# Helpers
# =============================================================================

.dcv_read_inertia <- function(m_a, verbose) {
  app_moi  <- .dcv_dbl(m_a("calibration_inertia_apparatus_moi_gram_millimeter_squared"))
  spec_moi <- .dcv_dbl(m_a("calibration_inertia_specimen_moi_gram_millimeter_squared"))
  bias     <- .dcv_dbl(m_a("calibration_inertia_bias_newton_meter"))
  axis     <- .dcv_chr(m_a("calibration_inertia_axis_sensor"))
  if (verbose) {
    cli::cli_alert_info(
      "Inertia: apparatus {signif(app_moi, 4)} g*mm^2, specimen {signif(spec_moi, 4)} g*mm^2, bias {signif(bias, 4)} N*m, axis {axis}"
    )
  }
  list(apparatus_moi_g_mm2 = app_moi, specimen_moi_g_mm2 = spec_moi,
       bias_newton_meter = bias, axis_sensor = axis)
}


.dcv_read_calibrated_torque <- function(h5r, m_a, axis_sensor, hub_path, verbose) {
  # ft_calibrated columns: 1=xForce, 2=yForce, 3=zForce, 4=xTorque, 5=yTorque, 6=zTorque
  axis_to_col <- c(x = 4L, y = 5L, z = 6L)
  torque_col  <- NA_integer_

  if (!is.na(axis_sensor) && nchar(axis_sensor) > 0L) {
    ax <- tolower(substr(axis_sensor, 1L, 1L))
    if (ax %in% names(axis_to_col)) torque_col <- axis_to_col[[ax]]
  }
  if (is.na(torque_col)) {
    pbax <- .dcv_chr(m_a("daq_primary_bending_axis"))
    if (!is.na(pbax)) {
      ax <- tolower(substr(pbax, 1L, 1L))
      if (ax %in% names(axis_to_col)) torque_col <- axis_to_col[[ax]]
    }
  }
  if (is.na(torque_col)) {
    if (verbose) cli::cli_warn("deconvolve: bending axis not determined; defaulting to xTorque")
    torque_col <- 4L
  }

  .read_ft <- function(h5_or_path) {
    tryCatch({
      m <- h5read(h5_or_path, "/derived/forcetorque_calibrated")
      if (!is.matrix(m)) m <- as.matrix(m)
      if (nrow(m) == 6L && ncol(m) != 6L) m <- t(m)
      if (ncol(m) == 6L) m else NULL
    }, error = function(e) NULL)
  }

  ft <- NULL
  der_path <- "/derived/forcetorque_calibrated"
  if (H5Lexists(h5r, der_path)) ft <- .read_ft(h5r)
  if (is.null(ft) && !is.null(hub_path) && file.exists(hub_path)) ft <- .read_ft(hub_path)
  if (is.null(ft)) return(NULL)

  ft[, torque_col]
}


.dcv_read_angle <- function(h5r, verbose) {
  # try single_finite flat path, then first step subgroup
  for (base in c("/timeseries", "/timeseries/step_001")) {
    tp <- paste0(base, "/time_second")
    if (!H5Lexists(h5r, tp)) next
    ap <- paste0(base, "/angle_measured_degree")
    if (!H5Lexists(h5r, ap)) ap <- paste0(base, "/angle_commanded_degree")
    t_s <- tryCatch(as.numeric(h5read(h5r, tp)), error = function(e) NULL)
    ang <- tryCatch(as.numeric(h5read(h5r, ap)), error = function(e) NULL)
    if (!is.null(t_s)) return(list(t_s = t_s, angle_deg = ang))
  }
  NULL
}


.dcv_write_hub <- function(hub_path, tau_corrected, verbose) {
  tryCatch({
    if (!file.exists(hub_path)) h5createFile(hub_path)
    if (!H5Lexists(hub_path, "/derived")) h5createGroup(hub_path, "derived")
    p <- "/derived/torque_inertia_corrected"
    if (H5Lexists(hub_path, p)) h5delete(hub_path, p)
    h5write(tau_corrected, hub_path, p)
    if (verbose) cli::cli_alert_success("torque_inertia_corrected written to hub: {hub_path}")
  }, error = function(e) {
    cli::cli_alert_danger("dcv_write_hub: {conditionMessage(e)}")
  })
}


.dcv_dbl <- function(v) {
  v <- suppressWarnings(as.numeric(v[1L]))
  if (length(v) == 0L || is.na(v)) NA_real_ else v
}

.dcv_chr <- function(v) {
  v <- as.character(v[1L])
  if (length(v) == 0L || is.na(v) || v %in% c("NA", "")) NA_character_ else v
}
