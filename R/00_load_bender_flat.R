# load_bender_flat.R
# Flat-schema loader for Bender Phase-1 HDF5 files (metadata/ + timeseries/).
#
# Replaces load_bender_data.R (Bender2-era /Calibrated, /NominalStimulus, /RawInput).
# Hard cut: reads ONLY the new flat schema. Old files require the legacy loader.
#
# HDF5 layout (as written by bender_h5_export.py post-M1):
#   /metadata/        -- most scalars as GROUP ATTRIBUTES; matrices + index arrays
#                         as DATASETS within the group
#   /timeseries/      -- flat channel datasets (single_finite) or step_NNN/
#                         subgroups (segmented_finite)
#   /derived/         -- written by exporter for RA inspection (not ground truth)
#
# Returns a tibble with the analysis-column names expected by 03_analyze.R:
#   t.s, t.norm, angle.deg, anglevel.degps, enc.deg, stim, cycle, halfcycle,
#   curve.invm, curvature.invm, curverate.invms, freq.Hz, amp.deg, duty, phase,
#   is_active_by_cycle, dclamp.m, fishlength.m, fishmass.g, fishcode,
#   stimulus_type, filename, fullpathname, trial
# Torque columns when derived/forcetorque_calibrated is present:
#   xtorque0.Nm (raw pre-bias), xtorque.Nm (bias-sub + filtered), and y/z variants.
#
# Tibble attrs: SampleFrequency.Hz, Filename, protocol_sampling_mode.

library(rhdf5)
library(dplyr)
library(tibble)

load_bender_flat <- function(
    filename,
    do_filter     = TRUE,
    cutoffmult    = 8.0,
    cutoff        = NA_real_,
    loadtorques   = c("x"),
    loadforces    = FALSE
) {
  loadtorques <- .bfl_validate_axes(loadtorques)
  loadforces  <- .bfl_validate_axes(loadforces)

  td <- NULL

  tryCatch({
    h5 <- H5Fopen(filename, "H5F_ACC_RDONLY")
    on.exit(try(H5Fclose(h5), silent = TRUE), add = TRUE)

    # -- metadata access -------------------------------------------------
    # Most scalars live as GROUP ATTRIBUTES of /metadata.
    # Index arrays and calibration matrices live as DATASETS within /metadata.
    m_attrs <- tryCatch(h5readAttributes(h5, "/metadata"), error = function(e) list())

    # Attribute accessor (returns NULL when absent)
    m_a <- function(key) m_attrs[[key]]

    # Dataset accessor (for index_cycle_*, index_step_*, calibration matrices, etc.)
    m_ds <- function(key) {
      path <- paste0("/metadata/", key)
      if (H5Lexists(h5, path)) tryCatch(h5read(h5, path), error = function(e) NULL) else NULL
    }

    # -- key scalars -------------------------------------------------------
    sampling_mode <- .bfl_chr(m_a("protocol_sampling_mode"), "single_finite")
    sampfreq      <- .bfl_dbl(m_a("daq_ai_sample_rate_hertz"), 1000.0)
    dclamp_mm     <- .bfl_dbl(m_a("measurement_clamp_separation_millimeter"), NA_real_)
    dclamp_m      <- if (is.finite(dclamp_mm)) dclamp_mm / 1000.0 else NA_real_
    fishcode      <- .bfl_chr(m_a("specimen_id"),  NA_character_)
    protocol_type <- .bfl_chr(m_a("protocol_type"), NA_character_)
    fl_mm <- .bfl_dbl(m_a("measurement_specimen_bodylength_millimeter"),
                      .bfl_dbl(m_a("measurement_specimen_standardlength_millimeter"), NA_real_))
    fishlength_m  <- if (is.finite(fl_mm)) fl_mm / 1000.0 else NA_real_
    fishmass_g    <- .bfl_dbl(m_a("measurement_specimen_body_mass_gram"), NA_real_)

    # -- read timeseries --------------------------------------------------
    td <- if (sampling_mode == "single_finite") {
      .bfl_read_single_finite(h5)
    } else {
      .bfl_read_segmented_finite(h5)
    }

    if (is.null(td) || nrow(td) == 0L) {
      message("[load_bender_flat] empty timeseries in ", basename(filename))
      return(NULL)
    }

    # -- stim: remap none/left/right/both -> factor 0/L/R/B --------------
    stim_levels <- c("0", "L", "R", "B")
    if ("stim_side" %in% names(td)) {
      td <- td |>
        mutate(stim = factor(
          dplyr::case_match(
            tolower(trimws(as.character(stim_side))),
            "left"  ~ "L",
            "right" ~ "R",
            "both"  ~ "B",
            .default = "0"
          ),
          levels = stim_levels
        )) |>
        select(-stim_side)
    } else {
      td$stim <- factor("0", levels = stim_levels)
    }

    # -- cycle / halfcycle ------------------------------------------------
    td <- td |>
      mutate(
        cycle     = floor(t.norm),
        halfcycle = floor(t.norm * 2.0) / 2.0
      )

    # -- per-cycle design (single_finite only) ----------------------------
    if (sampling_mode == "single_finite") {
      td <- .bfl_attach_cycle_design(h5, td, m_ds, m_a, dclamp_m)
    } else {
      td$freq.Hz        <- NA_real_
      td$amp.deg        <- NA_real_
      td$duty           <- NA_real_
      td$phase          <- NA_real_
      td$curvature.invm <- NA_real_
      td$is_active_by_cycle <- NA
    }

    # -- calibrated torque from derived/ ----------------------------------
    td <- .bfl_attach_torque(h5, td, loadtorques, loadforces)

    # -- metadata scalars to columns --------------------------------------
    fn_base <- basename(filename)
    td <- bind_cols(td, tibble(
      filename      = fn_base,
      fullpathname  = file.path(normalizePath(dirname(filename)), filename),
      fishcode      = fishcode,
      stimulus_type = protocol_type,
      dclamp.m      = dclamp_m,
      fishlength.m  = fishlength_m,
      fishmass.g    = fishmass_g,
      trial         = .bfl_trial_from_filename(fn_base)
    ))

    # -- curvature --------------------------------------------------------
    if (is.finite(dclamp_m) && "angle.deg" %in% names(td)) {
      td <- td |>
        mutate(
          curve.invm      = angle.deg * pi / 180.0 / dclamp_m,
          curverate.invms = (lead(curve.invm) - lag(curve.invm)) /
                              (lead(t.s) - lag(t.s))
        )
    }

    # -- bias subtraction -------------------------------------------------
    td <- .bfl_subtract_bias(td, loadtorques, loadforces)

    # -- Butterworth LP filter --------------------------------------------
    if (do_filter && requireNamespace("signal", quietly = TRUE)) {
      fc <- cutoff
      if (!is.finite(fc)) {
        max_freq <- suppressWarnings(max(td$freq.Hz, na.rm = TRUE))
        fc <- if (is.finite(max_freq) && max_freq > 0) cutoffmult * max_freq else 50.0
      }
      nyq <- 0.5 * sampfreq
      if (fc < nyq) {
        filt     <- signal::butter(9L, fc / nyq, type = "low")
        raw_cols <- grep("torque0\\.Nm$|force0\\.N$", names(td), value = TRUE)
        if (length(raw_cols) > 0L) {
          td <- td |>
            mutate(across(
              all_of(raw_cols),
              list(s = \(x) if (sum(is.finite(x)) > 20L) signal::filtfilt(filt, x) else x)
            ))
          td <- td |>
            rename_with(
              \(nm) sub("(\\w+)0\\.(Nm?)_s$", "\\1.\\2", nm),
              .cols = ends_with("_s")
            )
        }
      }
    }

    attr(td, "SampleFrequency.Hz")     <- sampfreq
    attr(td, "Filename")               <- filename
    attr(td, "protocol_sampling_mode") <- sampling_mode

  }, error = function(err) {
    cli::cli_alert_danger("load_bender_flat: {conditionMessage(err)}")
    td <<- NULL
  })

  td
}


# =============================================================================
# Internal: timeseries readers
# =============================================================================

.bfl_read_single_finite <- function(h5) {
  ts_ds <- function(nm) {
    path <- paste0("/timeseries/", nm)
    if (H5Lexists(h5, path)) tryCatch(h5read(h5, path), error = function(e) NULL) else NULL
  }

  t_s <- ts_ds("time_second"); if (is.null(t_s)) return(NULL)
  N   <- length(t_s)

  # HDF5 byte-string stim_side: rhdf5 reads |S8 as character, values "none"/"left"/"right"/"both"
  stim_raw <- ts_ds("stim_side")

  tibble(
    t.s           = as.numeric(t_s),
    t.norm        = .bfl_vec_dbl(ts_ds("time_normalized"),                              N),
    angle.deg     = .bfl_vec_dbl(ts_ds("angle_commanded_degree"),                       N),
    enc.deg       = .bfl_vec_dbl(ts_ds("angle_measured_degree"),                        N),
    anglevel.degps= .bfl_vec_dbl(ts_ds("angular_velocity_commanded_degree_per_second"), N),
    stim_side     = if (!is.null(stim_raw)) as.character(stim_raw) else rep("none", N),
    stim.V        = .bfl_vec_dbl(ts_ds("stim_channel1_command_volt"),                   N),
    cycle_index   = {
      ci <- ts_ds("cycle_index")
      if (!is.null(ci)) as.integer(ci) else rep(-1L, N)
    }
  )
}


.bfl_read_segmented_finite <- function(h5) {
  # enumerate step_NNN subgroups by listing the whole file (recursive=TRUE)
  # and filtering on group == "/timeseries". recursive=FALSE only yields root-level groups.
  ts_info <- tryCatch(
    h5ls(h5, recursive = TRUE),
    error = function(e) NULL
  )
  if (is.null(ts_info)) return(NULL)

  step_names <- sort(ts_info$name[
    ts_info$group == "/timeseries" & grepl("^step_\\d+$", ts_info$name)
  ])
  if (length(step_names) == 0L) return(NULL)

  steps <- vector("list", length(step_names))
  for (i in seq_along(step_names)) {
    sname    <- step_names[[i]]
    base     <- paste0("/timeseries/", sname, "/")
    step_num <- as.integer(sub("step_", "", sname))

    ds <- function(nm) {
      path <- paste0(base, nm)
      if (H5Lexists(h5, path)) tryCatch(h5read(h5, path), error = function(e) NULL) else NULL
    }

    t_s <- ds("time_second"); if (is.null(t_s)) next
    N   <- length(t_s)
    stim_raw <- ds("stim_side")

    steps[[i]] <- tibble(
      t.s            = as.numeric(t_s),
      t.norm         = .bfl_vec_dbl(ds("time_normalized"),                              N),
      angle.deg      = .bfl_vec_dbl(ds("angle_commanded_degree"),                       N),
      enc.deg        = .bfl_vec_dbl(ds("angle_measured_degree"),                        N),
      anglevel.degps = .bfl_vec_dbl(ds("angular_velocity_commanded_degree_per_second"), N),
      stim_side      = if (!is.null(stim_raw)) as.character(stim_raw) else rep("none", N),
      stim.V         = .bfl_vec_dbl(ds("stim_channel1_command_volt"),                   N),
      cycle_index    = rep(-1L, N),
      step_number    = step_num
    )
  }

  steps <- Filter(Negate(is.null), steps)
  if (length(steps) == 0L) return(NULL)
  bind_rows(steps)
}


# =============================================================================
# Internal: per-cycle design (single_finite)
# =============================================================================

.bfl_attach_cycle_design <- function(h5, td, m_ds, m_a, dclamp_m) {
  freq_cyc   <- tryCatch(as.numeric(m_ds("index_cycle_frequency_hertz")),      error = function(e) NULL)
  amp_cyc    <- tryCatch(as.numeric(m_ds("index_cycle_motor_amplitude_degree")), error = function(e) NULL)
  active_cyc <- tryCatch(m_ds("index_cycle_active"),                            error = function(e) NULL)
  duty_cyc   <- tryCatch(as.numeric(m_ds("index_cycle_activation_duty")),       error = function(e) NULL)
  phase_cyc  <- tryCatch(as.numeric(m_ds("index_cycle_activation_phase")),      error = function(e) NULL)
  op_cyc     <- tryCatch(as.numeric(m_ds("index_cycle_operating_point")),        error = function(e) NULL)
  op_units   <- tryCatch(as.character(m_ds("index_cycle_operating_point_units")), error = function(e) NULL)

  na_fill <- function() {
    td$freq.Hz <- NA_real_; td$amp.deg <- NA_real_; td$duty <- NA_real_
    td$phase <- NA_real_; td$curvature.invm <- NA_real_; td$is_active_by_cycle <- NA
    td
  }

  if (is.null(freq_cyc) || length(freq_cyc) == 0L) return(na_fill())
  C <- length(freq_cyc)

  # Curvature per cycle from operating_point (when units mention per_meter) or amplitude
  curv_cyc <- if (!is.null(op_cyc) && !is.null(op_units) &&
                    any(grepl("per_meter|curvature", op_units, ignore.case = TRUE))) {
    as.numeric(op_cyc)
  } else if (!is.null(amp_cyc) && is.finite(dclamp_m)) {
    as.numeric(amp_cyc) * pi / 180.0 / dclamp_m
  } else {
    rep(NA_real_, C)
  }

  # cycle_index in timeseries is 1-based per schema (matches index_cycle_index)
  ci       <- td$cycle_index
  in_range <- !is.na(ci) & ci >= 1L & ci <= C

  td$freq.Hz        <- NA_real_
  td$amp.deg        <- NA_real_
  td$duty           <- NA_real_
  td$phase          <- NA_real_
  td$curvature.invm <- NA_real_
  td$is_active_by_cycle <- NA

  if (any(in_range)) {
    idx <- ci[in_range]
    td$freq.Hz[in_range]         <- freq_cyc[idx]
    if (!is.null(amp_cyc))   td$amp.deg[in_range]        <- amp_cyc[idx]
    if (!is.null(duty_cyc))  td$duty[in_range]            <- duty_cyc[idx]
    if (!is.null(phase_cyc)) td$phase[in_range]           <- phase_cyc[idx]
    td$curvature.invm[in_range]  <- curv_cyc[idx]
    if (!is.null(active_cyc)) {
      al <- as.logical(active_cyc)
      td$is_active_by_cycle[in_range] <- al[idx]
    }
  }

  # blank design cols outside the bending window
  outside <- !is.na(td$t.s) & td$t.s <= 0
  td$freq.Hz[outside] <- NA_real_; td$amp.deg[outside] <- NA_real_
  td$duty[outside]    <- NA_real_; td$phase[outside]   <- NA_real_
  td
}


# =============================================================================
# Internal: calibrated torque
# =============================================================================

.bfl_attach_torque <- function(h5, td, loadtorques, loadforces) {
  # derived/forcetorque_calibrated: Python [6, N] -> rhdf5 gives matrix dim [N, 6]
  # Column order (matches exporter _ft_names):
  #   1=xForce, 2=yForce, 3=zForce, 4=xTorque, 5=yTorque, 6=zTorque
  der_path <- "/derived/forcetorque_calibrated"
  if (!H5Lexists(h5, der_path)) return(td)

  ft <- tryCatch({
    m <- h5read(h5, der_path)
    if (!is.matrix(m)) m <- as.matrix(m)
    # ensure [N x 6]
    if (nrow(m) == 6L && ncol(m) != 6L) m <- t(m)
    if (ncol(m) != 6L) stop("unexpected shape")
    m
  }, error = function(e) {
    cli::cli_warn("load_bender_flat: could not read derived/forcetorque_calibrated: {conditionMessage(e)}")
    NULL
  })
  if (is.null(ft)) return(td)

  # for segmented, derived/ is concatenated across all steps -- lengths must match
  if (nrow(ft) != nrow(td)) {
    cli::cli_warn(
      "load_bender_flat: derived/ length {nrow(ft)} != timeseries length {nrow(td)}; torque skipped"
    )
    return(td)
  }

  if ("x" %in% loadtorques) td$xtorque0.Nm <- ft[, 4L]
  if ("y" %in% loadtorques) td$ytorque0.Nm <- ft[, 5L]
  if ("z" %in% loadtorques) td$ztorque0.Nm <- ft[, 6L]
  if ("x" %in% loadforces)  td$xforce0.N   <- ft[, 1L]
  if ("y" %in% loadforces)  td$yforce0.N   <- ft[, 2L]
  if ("z" %in% loadforces)  td$zforce0.N   <- ft[, 3L]
  td
}


# =============================================================================
# Internal: bias subtraction
# =============================================================================

.bfl_subtract_bias <- function(td, loadtorques, loadforces) {
  raw_cols <- intersect(
    c(if ("x" %in% loadtorques) "xtorque0.Nm",
      if ("y" %in% loadtorques) "ytorque0.Nm",
      if ("z" %in% loadtorques) "ztorque0.Nm",
      if ("x" %in% loadforces)  "xforce0.N",
      if ("y" %in% loadforces)  "yforce0.N",
      if ("z" %in% loadforces)  "zforce0.N"),
    names(td)
  )
  if (length(raw_cols) == 0L || !"t.s" %in% names(td)) return(td)
  pre_stim <- !is.na(td$t.s) & td$t.s < -0.1
  if (!any(pre_stim)) return(td)
  for (col in raw_cols) {
    b <- mean(td[[col]][pre_stim], na.rm = TRUE)
    if (is.finite(b)) td[[col]] <- td[[col]] - b
  }
  td
}


# =============================================================================
# Internal: small helpers
# =============================================================================

.bfl_validate_axes <- function(ax) {
  if (is.logical(ax) && length(ax) == 1L) {
    if (isTRUE(ax)) c("x", "y", "z") else character(0)
  } else if (is.character(ax) && length(ax) == 1L) {
    switch(ax, all = c("x", "y", "z"), none = character(0), c(ax))
  } else {
    ax
  }
}

.bfl_dbl <- function(v, default) {
  v <- suppressWarnings(as.numeric(v[1L]))
  if (length(v) == 0L || is.na(v)) default else v
}

.bfl_chr <- function(v, default) {
  v <- as.character(v[1L])
  if (length(v) == 0L || is.na(v) || v %in% c("NA", "")) default else v
}

.bfl_vec_dbl <- function(v, N) {
  if (is.null(v)) rep(NA_real_, N) else as.numeric(v)
}

.bfl_trial_from_filename <- function(fn) {
  m <- sub(".*_bender_(\\d+)\\.h5$", "\\1", fn, ignore.case = TRUE)
  if (!identical(m, fn) && grepl("^[0-9]+$", m)) return(as.numeric(m))
  m2 <- sub(".+_(\\d+)\\.h5$", "\\1", fn, ignore.case = TRUE)
  if (!identical(m2, fn) && grepl("^[0-9]+$", m2)) return(as.numeric(m2))
  NA_real_
}
