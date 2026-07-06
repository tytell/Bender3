# 01_calibrate.R
# Apply embedded sensor calibration values to the raw aidata archive.
#
# Architecture (D11, M2b):
#   - Raw file holds the immutable 8-channel aidata + calibration VALUES embedded
#     in /metadata (most as GROUP ATTRIBUTES; matrices as DATASETS within the group).
#   - calibrated F/T is computed here and optionally written to a hub HDF5 file
#     under derived/forcetorque_calibrated.
#   - Sono segment lengths are likewise interpolated and written to hub derived/sono_*.
#   - This script NEVER modifies the raw archive file.
#
# Usage:
#   source("R/01_calibrate.R")
#   calibrate_bender(raw_path = "path/to/raw.h5", hub_path = "path/to/hub.h5")
#
# If hub_path is NULL, results are returned as a list rather than written.

library(rhdf5)
library(dplyr)
library(jsonlite)

calibrate_bender <- function(raw_path, hub_path = NULL, verbose = TRUE) {
  result <- list()

  tryCatch({
    h5r <- H5Fopen(raw_path, "H5F_ACC_RDONLY")
    on.exit(try(H5Fclose(h5r), silent = TRUE), add = TRUE)

    # Attribute accessor (most scalars are stored as /metadata GROUP ATTRIBUTES)
    m_attrs <- tryCatch(h5readAttributes(h5r, "/metadata"), error = function(e) list())
    m_a  <- function(key) m_attrs[[key]]

    # Dataset accessor (calibration matrices, sono breakpoints, channel lists)
    m_ds <- function(key) {
      path <- paste0("/metadata/", key)
      if (H5Lexists(h5r, path)) tryCatch(h5read(h5r, path), error = function(e) NULL) else NULL
    }

    # -- F/T calibration ---------------------------------------------------
    ft_result <- tryCatch(
      .cal_forcetorque(h5r, m_ds, m_a, verbose),
      error = function(e) {
        if (verbose) cli::cli_alert_warning("F/T calibration skipped: {conditionMessage(e)}")
        NULL
      }
    )
    if (!is.null(ft_result)) result$forcetorque_calibrated <- ft_result

    # -- sono calibration --------------------------------------------------
    sono_result <- tryCatch(
      .cal_sono(h5r, m_ds, m_a, verbose),
      error = function(e) {
        if (verbose) cli::cli_alert_warning("Sono calibration skipped: {conditionMessage(e)}")
        NULL
      }
    )
    if (!is.null(sono_result)) result <- c(result, sono_result)

    # -- inertia parameters for downstream (all stored as attrs) -----------
    result$inertia <- list(
      apparatus_moi_g_mm2 = .cal_dbl(m_a("calibration_inertia_apparatus_moi_gram_millimeter_squared")),
      specimen_moi_g_mm2  = .cal_dbl(m_a("calibration_inertia_specimen_moi_gram_millimeter_squared")),
      bias_newton_meter   = .cal_dbl(m_a("calibration_inertia_bias_newton_meter")),
      axis_sensor         = .cal_chr(m_a("calibration_inertia_axis_sensor"))
    )

  }, error = function(e) {
    cli::cli_alert_danger("calibrate_bender: {conditionMessage(e)}")
  })

  # -- write to hub or return --------------------------------------------
  if (!is.null(hub_path) && length(result) > 0L) {
    .cal_write_hub(hub_path, raw_path, result, verbose)
  }

  invisible(result)
}


# -- F/T calibration block -----------------------------------------------

.cal_forcetorque <- function(h5r, m_ds, m_a, verbose) {
  # calibration_forcetorque_matrix is a DATASET: [6 x 6] in Python/HDF5 C-order.
  # rhdf5 reverses dim order: the R matrix is transposed relative to Python.
  cal_raw <- m_ds("calibration_forcetorque_matrix")
  if (is.null(cal_raw)) stop("calibration_forcetorque_matrix absent (file not calibrated)")
  cal_mat_R <- as.matrix(cal_raw)
  if (!all(dim(cal_mat_R) == c(6L, 6L))) {
    stop("calibration_forcetorque_matrix is not 6x6 (got ", paste(dim(cal_mat_R), collapse = "x"), ")")
  }

  # daq_ai_channel_map is a GROUP ATTRIBUTE (JSON string)
  # format: {"0": "ai0:xForce", "1": "ai1:yForce", ...}
  chan_map_raw <- m_a("daq_ai_channel_map")
  if (is.null(chan_map_raw)) stop("daq_ai_channel_map attribute absent")
  chan_map <- tryCatch(fromJSON(as.character(chan_map_raw)), error = function(e) NULL)
  if (is.null(chan_map)) stop("daq_ai_channel_map could not be parsed as JSON")

  ft_names_ordered <- c("xForce", "yForce", "zForce", "xTorque", "yTorque", "zTorque")
  chan_names <- as.character(unname(chan_map))

  # find R column index (1-based) for each F/T channel
  ft_cols_R <- vapply(ft_names_ordered, function(nm) {
    idx <- which(grepl(nm, chan_names, fixed = TRUE))
    if (length(idx) == 0L) NA_integer_ else as.integer(names(chan_map)[idx[1L]]) + 1L
  }, integer(1L))

  if (any(is.na(ft_cols_R))) {
    stop("F/T channels not found in daq_ai_channel_map: ",
         paste(ft_names_ordered[is.na(ft_cols_R)], collapse = ", "))
  }

  # aidata is [8 x N] in Python/HDF5; rhdf5 gives [N x 8] in R (dims reversed)
  aidata_R <- h5read(h5r, "/timeseries/aidata")
  if (!is.matrix(aidata_R)) aidata_R <- as.matrix(aidata_R)
  if (ncol(aidata_R) != 8L && nrow(aidata_R) == 8L) aidata_R <- t(aidata_R)
  if (ncol(aidata_R) != 8L) stop("aidata does not have 8 channels (got ", ncol(aidata_R), ")")

  aidata_ft <- aidata_R[, ft_cols_R, drop = FALSE]  # [N x 6]

  # Python: ft_cal = (aidata_ft.T @ cal_arr_python).T  =>  equivalent in R:
  # cal_mat_R = t(cal_arr_python), so t(cal_mat_R) = cal_arr_python
  ft_calibrated <- aidata_ft %*% t(cal_mat_R)  # [N x 6]
  colnames(ft_calibrated) <- ft_names_ordered

  if (verbose) cli::cli_alert_success(
    "F/T calibration applied: {nrow(ft_calibrated)} samples, 6 channels"
  )
  ft_calibrated
}


# -- sono calibration block -----------------------------------------------

.cal_sono <- function(h5r, m_ds, m_a, verbose) {
  result <- list()
  # daq_ai_channel_map is an attr
  chan_map_raw <- m_a("daq_ai_channel_map")
  chan_map <- tryCatch(fromJSON(as.character(chan_map_raw)), error = function(e) NULL)

  for (side in c("left", "right")) {
    key <- paste0("calibration_sono_", side, "_volt_millimeter_breakpoints")
    bp  <- tryCatch(as.numeric(m_ds(key)), error = function(e) NULL)
    if (is.null(bp) || length(bp) < 4L || any(!is.finite(bp))) next

    if (is.null(chan_map)) next
    chan_names <- as.character(unname(chan_map))
    sono_label <- paste0("sono_", side)
    idx <- which(grepl(sono_label, chan_names, fixed = TRUE))
    if (length(idx) == 0L) next

    sono_col_R <- as.integer(names(chan_map)[idx[1L]]) + 1L

    aidata_R <- tryCatch({
      m <- h5read(h5r, "/timeseries/aidata")
      if (!is.matrix(m)) m <- as.matrix(m)
      if (ncol(m) != 8L && nrow(m) == 8L) m <- t(m)
      m
    }, error = function(e) NULL)
    if (is.null(aidata_R)) next

    sono_v  <- aidata_R[, sono_col_R]
    sono_mm <- .cal_sono_interp(sono_v, bp)
    if (!is.null(sono_mm)) {
      result[[paste0("sono_", side, "_mm")]] <- sono_mm
      if (verbose) cli::cli_alert_success(
        "Sono {side}: calibrated {length(sono_mm)} samples"
      )
    }
  }
  result
}


.cal_sono_interp <- function(v_raw, bp) {
  # bp layout: [Low_V, High_V, Low_mm, High_mm] for 2-point;
  # or [V1, V2, ..., Vn, mm1, mm2, ..., mmn] for multi-point
  if (length(bp) < 4L || length(bp) %% 2L != 0L) return(NULL)
  n_pts  <- length(bp) / 2L
  v_pts  <- bp[seq_len(n_pts)]
  mm_pts <- bp[seq(n_pts + 1L, 2L * n_pts)]
  if (n_pts == 2L) {
    slope <- (mm_pts[2L] - mm_pts[1L]) / (v_pts[2L] - v_pts[1L])
    v_raw * slope + (mm_pts[1L] - v_pts[1L] * slope)
  } else {
    approx(v_pts, mm_pts, xout = v_raw, rule = 2L)$y
  }
}


# -- hub write -----------------------------------------------------------

.cal_write_hub <- function(hub_path, raw_path, result, verbose) {
  tryCatch({
    if (!file.exists(hub_path)) {
      h5createFile(hub_path)
      h5r2 <- H5Fopen(raw_path, "H5F_ACC_RDONLY")
      sv   <- tryCatch(h5readAttributes(h5r2, "/")$schema_version, error = function(e) NULL)
      H5Fclose(h5r2)
      if (!is.null(sv)) h5writeAttribute(sv, hub_path, "/", "schema_version")
    }
    if (!H5Lexists(hub_path, "/derived")) h5createGroup(hub_path, "derived")
    if (!is.null(result$forcetorque_calibrated)) {
      der <- "/derived/forcetorque_calibrated"
      if (H5Lexists(hub_path, der)) h5delete(hub_path, der)
      h5write(result$forcetorque_calibrated, hub_path, der)
    }
    for (nm in grep("^sono_", names(result), value = TRUE)) {
      p <- paste0("/derived/", nm)
      if (H5Lexists(hub_path, p)) h5delete(hub_path, p)
      h5write(result[[nm]], hub_path, p)
    }
    if (verbose) cli::cli_alert_success("Calibrated output written to hub: {hub_path}")
  }, error = function(e) {
    cli::cli_alert_danger("cal_write_hub: {conditionMessage(e)}")
  })
}


# -- scalar helpers -------------------------------------------------------

.cal_dbl <- function(v) {
  v <- suppressWarnings(as.numeric(v[1L]))
  if (length(v) == 0L || is.na(v)) NA_real_ else v
}

.cal_chr <- function(v) {
  v <- as.character(v[1L])
  if (length(v) == 0L || is.na(v) || v %in% c("NA", "")) NA_character_ else v
}
