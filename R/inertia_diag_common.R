# inertia_diag_common.R
# Shared low-level readers + constants for the inertial-correction DIAGNOSTICS
# (diag_apparatus_inertia.R, diag_specimen_inertia.R). READ-ONLY: these helpers
# open raw .h5 files, compute kinematics, and never write to the archive or the
# production path.
#
# Kept in one place so both diagnostics use the SAME channel order, alpha
# definition, and cache -- a divergence between them would make the apparatus
# vs specimen comparison meaningless.

suppressPackageStartupMessages({ library(rhdf5); library(dplyr); library(tibble); library(purrr) })

# ATI fixed channel order in derived/forcetorque_calibrated (see 02_deconvolve.R
# .dcv_read_calibrated_torque: 1=xForce ... 6=zTorque). daq_forcetorque_channel_
# names is absent on these files, so the fixed order is the only option.
CHANNELS <- c("xForce", "yForce", "zForce", "xTorque", "yTorque", "zTorque")

g_mm2_per_kg_m2 <- 1e9          # kg*m^2 -> g*mm^2
kg_m2_per_g_mm2 <- 1e-9         # g*mm^2 -> kg*m^2

# OneDrive reads are slow; cache extracted signals per file. Default (no-stim)
# reads keep the original "<md5>_<base>.rds" key so previously cached corpus
# entries stay valid; stim-enabled reads use an "S_" prefix (they carry the
# active mask + extra geometry the default records may predate).
CACHE_DIR <- file.path(Sys.getenv("TMPDIR", "/tmp"), "bender_inertia_diag_cache")
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

#' np.gradient equivalent for uniform spacing dt (central interior, one-sided edges).
.np_gradient <- function(y, dt) {
  n <- length(y)
  if (n < 2L) return(rep(NA_real_, n))
  g <- numeric(n)
  g[1L] <- (y[2L] - y[1L]) / dt
  g[n]  <- (y[n] - y[n - 1L]) / dt
  if (n > 2L) g[2L:(n - 1L)] <- (y[3L:n] - y[1L:(n - 2L)]) / (2 * dt)
  g
}

#' Angular velocity (deg/s) and acceleration (rad/s^2) from an angle-degree
#' series, computed WITHIN each contiguous run so segmented_finite step seams
#' (diff(t) <= 0) are never differenced across.
#' @return list(vel_deg_s, alpha_rad_s2)
.kinematics <- function(t_s, angle_deg) {
  dt <- suppressWarnings(stats::median(diff(t_s), na.rm = TRUE))
  if (!is.finite(dt) || dt <= 0) dt <- 1e-3
  run_start <- c(TRUE, diff(t_s) <= 0)
  run <- cumsum(run_start)
  vel <- rep(NA_real_, length(angle_deg))
  alpha <- rep(NA_real_, length(angle_deg))
  for (r in unique(run)) {
    idx <- which(run == r)
    if (length(idx) < 5L) next
    v <- .np_gradient(angle_deg[idx], dt)
    vel[idx] <- v
    alpha[idx] <- .np_gradient(v, dt) * (pi / 180.0)
  }
  list(vel_deg_s = vel, alpha_rad_s2 = alpha)
}

#' Backwards-compatible alpha-only helper (rad/s^2).
.alpha_rad_s2 <- function(t_s, angle_deg) .kinematics(t_s, angle_deg)$alpha_rad_s2

#' Read (t_s, angle_deg, ft[N x 6], step_number, active mask, specimen_moi, aor,
#' width, and a set of geometry attrs) from a raw .h5, handling flat
#' single_finite vs segmented_finite. Cached to RDS.
#' @param with_stim if TRUE, also read the stim command voltages and return an
#'   `active` logical (|stim| > 0) so callers can select PASSIVE samples.
.read_trial <- function(fp, with_stim = FALSE) {
  key <- file.path(CACHE_DIR, paste0(if (with_stim) "S_" else "",
                                     tools::md5sum(fp), "_", basename(fp), ".rds"))
  if (file.exists(key)) return(readRDS(key))

  out <- tryCatch({
    h5r <- H5Fopen(fp, "H5F_ACC_RDONLY")
    on.exit(try(H5Fclose(h5r), silent = TRUE), add = TRUE)
    m_attrs <- tryCatch(h5readAttributes(h5r, "/metadata"), error = function(e) list())
    dbl1 <- function(k) { v <- suppressWarnings(as.numeric(m_attrs[[k]][1L])); if (length(v) == 0L) NA_real_ else v }

    ft <- h5read(h5r, "/derived/forcetorque_calibrated")
    if (!is.matrix(ft)) ft <- as.matrix(ft)
    if (nrow(ft) == 6L && ncol(ft) != 6L) ft <- t(ft)
    stopifnot(ncol(ft) == 6L)

    ls <- tryCatch(h5ls(h5r, recursive = TRUE), error = function(e) NULL)
    read_stim_run <- function(base) {
      if (!with_stim) return(NULL)
      g <- function(nm) if (H5Lexists(h5r, paste0(base, nm))) as.numeric(h5read(h5r, paste0(base, nm))) else NULL
      v1 <- g("stim_channel1_command_volt"); v2 <- g("stim_channel2_command_volt")
      if (is.null(v1) && is.null(v2)) {
        ss <- g("stim_state"); if (is.null(ss)) return(NULL) else return(abs(ss) > 0)
      }
      n <- max(length(v1), length(v2))
      a1 <- if (is.null(v1)) rep(0, n) else v1; a2 <- if (is.null(v2)) rep(0, n) else v2
      (abs(a1) > 1e-9) | (abs(a2) > 1e-9)
    }

    flat_t <- "/timeseries/time_second"
    if (H5Lexists(h5r, flat_t)) {
      t_s <- as.numeric(h5read(h5r, flat_t))
      ap <- if (H5Lexists(h5r, "/timeseries/angle_measured_degree")) "/timeseries/angle_measured_degree" else "/timeseries/angle_commanded_degree"
      ang <- as.numeric(h5read(h5r, ap))
      step_number <- rep(1L, length(t_s))
      active <- read_stim_run("/timeseries/")
      if (is.null(active)) active <- rep(FALSE, length(t_s))
    } else {
      steps <- sort(ls$name[ls$group == "/timeseries" & grepl("^step_\\d+$", ls$name)])
      t_parts <- a_parts <- s_parts <- act_parts <- vector("list", length(steps))
      for (i in seq_along(steps)) {
        b <- paste0("/timeseries/", steps[i], "/")
        t_parts[[i]] <- as.numeric(h5read(h5r, paste0(b, "time_second")))
        ap <- if (H5Lexists(h5r, paste0(b, "angle_measured_degree"))) paste0(b, "angle_measured_degree") else paste0(b, "angle_commanded_degree")
        a_parts[[i]] <- as.numeric(h5read(h5r, ap))
        s_parts[[i]] <- rep(as.integer(sub("step_", "", steps[i])), length(t_parts[[i]]))
        ar <- read_stim_run(b)
        act_parts[[i]] <- if (is.null(ar)) rep(FALSE, length(t_parts[[i]])) else ar[seq_len(length(t_parts[[i]]))]
      }
      t_s <- unlist(t_parts); ang <- unlist(a_parts); step_number <- unlist(s_parts)
      active <- unlist(act_parts)
    }
    n <- min(nrow(ft), length(t_s), length(ang))
    if (length(active) < n) active <- c(active, rep(FALSE, n - length(active)))
    list(t_s = t_s[seq_len(n)], angle_deg = ang[seq_len(n)],
         ft = ft[seq_len(n), , drop = FALSE], step_number = step_number[seq_len(n)],
         active = as.logical(active[seq_len(n)]),
         specimen_moi_gmm2 = dbl1("calibration_inertia_specimen_moi_gram_millimeter_squared"),
         aor_mm = dbl1("calibration_inertia_apparatus_aor_to_clamp_millimeter"),
         width_mm = dbl1("calibration_inertia_apparatus_plate_to_plate_millimeter"),
         body_height_mm = dbl1("measurement_specimen_local_body_height_millimeter"),
         body_width_mm  = dbl1("measurement_specimen_local_body_width_millimeter"),
         clamp_sep_mm   = dbl1("measurement_clamp_separation_millimeter"),
         clamp_offset_vert_mm = dbl1("measurement_clamp_offset_vertical_millimeter"),
         density_g_mm3  = dbl1("measurement_specimen_density_gram_per_cubic_millimeter"))
  }, error = function(e) { cli::cli_alert_danger("read {basename(fp)}: {conditionMessage(e)}"); NULL })

  if (!is.null(out)) saveRDS(out, key)
  out
}

#' Per-channel empirical inertia from channel-vs-alpha OLS (all 6 axes).
#' I = |slope| in kg*m^2 (channel N*m per alpha rad/s^2), reported also in g*mm^2.
.per_channel_inertia <- function(ft, alpha) {
  purrr::map_dfr(seq_len(6L), function(j) {
    y <- ft[, j]; ok <- is.finite(y) & is.finite(alpha)
    if (sum(ok) < 20L) return(tibble(channel = CHANNELS[j], slope_kg_m2 = NA_real_,
                                     I_gmm2 = NA_real_, r2 = NA_real_, slope_sign = NA_real_))
    fit <- stats::lm(y[ok] ~ alpha[ok])
    sl <- unname(stats::coef(fit)[2L])
    tibble(channel = CHANNELS[j], slope_kg_m2 = sl, I_gmm2 = abs(sl) * g_mm2_per_kg_m2,
           r2 = summary(fit)$r.squared, slope_sign = sign(sl))
  })
}

#' Per-channel inertia SEPARATED from elasticity: channel ~ angle + alpha.
#' The alpha coefficient is the inertial term (kg*m^2); the angle coefficient is
#' the elastic term (N*m/deg). Use on PASSIVE (no-stim) samples so muscle force
#' does not contaminate the fit. Returns both coefficients + partial-R^2 of alpha.
.per_channel_inertia_vs_elastic <- function(ft, angle_deg, alpha) {
  purrr::map_dfr(seq_len(6L), function(j) {
    y <- ft[, j]; ok <- is.finite(y) & is.finite(alpha) & is.finite(angle_deg)
    if (sum(ok) < 50L) return(tibble(channel = CHANNELS[j], inertia_kg_m2 = NA_real_,
                                     I_gmm2 = NA_real_, elastic_Nm_per_deg = NA_real_,
                                     r2_full = NA_real_, alpha_partial_r2 = NA_real_))
    yy <- y[ok]; aa <- alpha[ok]; gg <- angle_deg[ok]
    full <- stats::lm(yy ~ gg + aa)
    reduced <- stats::lm(yy ~ gg)  # elastic only
    sl_alpha <- unname(stats::coef(full)[3L])
    sl_ang   <- unname(stats::coef(full)[2L])
    ss_res_full <- sum(stats::resid(full)^2); ss_res_red <- sum(stats::resid(reduced)^2)
    partial_r2 <- if (ss_res_red > 0) (ss_res_red - ss_res_full) / ss_res_red else NA_real_
    tibble(channel = CHANNELS[j], inertia_kg_m2 = sl_alpha, I_gmm2 = abs(sl_alpha) * g_mm2_per_kg_m2,
           elastic_Nm_per_deg = sl_ang, r2_full = summary(full)$r.squared,
           alpha_partial_r2 = partial_r2)
  })
}
