# analyze_frequency_sweep.R
# Passive stiffness/damping vs. instantaneous frequency for frequency_sweep
# trials (its own protocol category, per PI decision 2026-07-15: structurally
# distinct from "Dynamic" -- 100% stim_type=="passive", no index_cycle_*
# per-cycle design grid, continuous chirp).
#
# The schema's instantaneous_frequency_hertz timeseries field (context_jlab_cg
# _h5schema.md \u00a73b) is NOT populated in this export (verified against the
# bass16 corpus, 2026-07-15), so instantaneous frequency is derived here from
# the commanded-angle waveform itself via zero-up-crossing timing -- one
# frequency + stiffness/damping estimate per bracketed oscillation cycle.
#
# No Force_muscle / FL / FV / power output here: frequency_sweep is
# passive-only (no active stimulation), so active-minus-passive muscle force
# is undefined (PI-confirmed scope, 2026-07-15) -- only passive mechanics
# (stiffness EI, damping etaI) vs. frequency is reported.

library(dplyr)

#' Zero-up-crossing cycle boundaries for a chirp waveform.
#' @return tibble(cycle_index, t_start, t_end, period_s, frequency_hz)
.detect_chirp_cycles <- function(t_s, angle_deg, min_amplitude_deg = 0.05) {
  ok <- is.finite(t_s) & is.finite(angle_deg)
  t_s <- t_s[ok]; angle_deg <- angle_deg[ok]
  ord <- order(t_s)
  t_s <- t_s[ord]; angle_deg <- angle_deg[ord]

  sgn <- sign(angle_deg)
  sgn[sgn == 0] <- NA_real_
  up_cross <- which(diff(sgn) == 2L) # -1 -> +1 transitions (zero-up-crossing)
  if (length(up_cross) < 2L) {
    return(tibble::tibble(cycle_index = integer(0), t_start = numeric(0),
                          t_end = numeric(0), period_s = numeric(0),
                          frequency_hz = numeric(0)))
  }

  t_cross <- t_s[up_cross]
  n_cyc <- length(t_cross) - 1L
  out <- tibble::tibble(
    cycle_index  = seq_len(n_cyc),
    t_start      = t_cross[seq_len(n_cyc)],
    t_end        = t_cross[seq_len(n_cyc) + 1L],
    period_s     = diff(t_cross),
    frequency_hz = 1.0 / diff(t_cross)
  )

  # drop cycles from the near-zero wait_before/wait_after bookends (amplitude
  # too small to trust a zero-crossing-derived frequency estimate).
  amp_ok <- vapply(seq_len(n_cyc), function(i) {
    win <- t_s >= out$t_start[i] & t_s <= out$t_end[i]
    max(abs(angle_deg[win]), na.rm = TRUE) >= min_amplitude_deg
  }, logical(1L))
  out[amp_ok, ]
}


#' Per-cycle passive stiffness (EI) / damping (etaI) vs. instantaneous
#' frequency for a loaded frequency_sweep td (with torque_inertia_corrected_Nm
#' attached, see deconvolve_bender()).
#'
#' @param td Tibble from load_bender_flat(), single_finite, with curve.invm /
#'   curverate.invms (computed by load_bender_flat when dclamp is available)
#'   and torque_inertia_corrected_Nm attached.
#' @return tibble: cycle_index, t_start, t_end, frequency_hz, EI1.Nm2,
#'   etaI1.Nm2s, peak_torque_Nm, mean_curvature_amplitude_invm, status
#'   ("ok"/"failed" per cycle -- e.g. degenerate/too-few-sample cycles).
analyze_frequency_sweep <- function(td, torque_col = "torque_inertia_corrected_Nm",
                                    min_amplitude_deg = 0.05, min_samples_per_cycle = 10L) {
  stopifnot(is.data.frame(td))
  req <- c("t.s", "angle.deg", "curve.invm", "curverate.invms", torque_col)
  missing <- setdiff(req, names(td))
  if (length(missing) > 0L) {
    cli::cli_abort("analyze_frequency_sweep: missing required column(s): {paste(missing, collapse=', ')}")
  }

  cycles <- .detect_chirp_cycles(td$t.s, td$angle.deg, min_amplitude_deg)
  if (nrow(cycles) == 0L) {
    cli::cli_warn("analyze_frequency_sweep: no chirp cycles detected (check amplitude / angle_commanded_degree)")
    return(tibble::tibble())
  }

  results <- vector("list", nrow(cycles))
  for (i in seq_len(nrow(cycles))) {
    win <- td$t.s >= cycles$t_start[i] & td$t.s < cycles$t_end[i]
    n_win <- sum(win, na.rm = TRUE)
    base_row <- list(
      cycle_index = cycles$cycle_index[i], t_start = cycles$t_start[i],
      t_end = cycles$t_end[i], frequency_hz = cycles$frequency_hz[i]
    )
    if (n_win < min_samples_per_cycle) {
      results[[i]] <- c(base_row, list(
        EI1.Nm2 = NA_real_, etaI1.Nm2s = NA_real_, peak_torque_Nm = NA_real_,
        mean_curvature_amplitude_invm = NA_real_,
        status = "failed", reason = sprintf("only %d samples in cycle (need >= %d)", n_win, min_samples_per_cycle)
      ))
      next
    }

    cv <- td$curve.invm[win]; cr <- td$curverate.invms[win]; tq <- td[[torque_col]][win]
    ok <- is.finite(cv) & is.finite(cr) & is.finite(tq)
    if (sum(ok) < min_samples_per_cycle) {
      results[[i]] <- c(base_row, list(
        EI1.Nm2 = NA_real_, etaI1.Nm2s = NA_real_, peak_torque_Nm = NA_real_,
        mean_curvature_amplitude_invm = NA_real_,
        status = "failed", reason = "too many non-finite samples within cycle"
      ))
      next
    }
    cv <- cv[ok]; cr <- cr[ok]; tq <- tq[ok]

    ei  <- tryCatch(.avg_mech(cv, tq), error = function(e) NA_real_)
    eta <- tryCatch(.avg_mech(cr, tq), error = function(e) NA_real_)
    status <- if (is.finite(ei) && is.finite(eta)) "ok" else "failed"
    reason <- if (status == "failed") "stiffness/damping regression returned non-finite" else NA_character_

    results[[i]] <- c(base_row, list(
      EI1.Nm2 = ei, etaI1.Nm2s = eta,
      peak_torque_Nm = max(abs(tq), na.rm = TRUE),
      mean_curvature_amplitude_invm = (max(cv, na.rm = TRUE) - min(cv, na.rm = TRUE)) / 2.0,
      status = status, reason = reason
    ))
  }

  dplyr::bind_rows(lapply(results, tibble::as_tibble))
}
