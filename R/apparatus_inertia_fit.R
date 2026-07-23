# apparatus_inertia_fit.R
# R-side evaluator for the versioned apparatus-inertia fit artifact
# (bender_apparatus_inertia_fit_v1 JSON, e.g.
# 2026-07-06_apparatus_inertia_calibration.json).
#
# WHY THIS EXISTS (Gap 1): the acquisition path stores raw torque + correction
# PARAMETERS only and leaves the apparatus MOI to be resolved POST-HOC in R
# (bender_functions.py ~line 469). That R-side step was never built:
# calibration_inertia_apparatus_moi_gram_millimeter_squared is NaN in every real
# trial file, so 02_deconvolve.R silently drops the apparatus term. This file
# ports bender_functions.py's apparatus_inertia_from_fit() (~line 6661) to R so
# the fit can be evaluated at each trial's OWN geometry.
#
# This module is a pure evaluator: it reads a REFERENCED JSON artifact (never
# embedded, matching the forcetorque_calibration_file convention) and returns a
# MOI. It performs NO correction and writes nothing. Nothing in the production
# pipeline calls it yet -- wiring into 02_deconvolve.R is a separate, gated step.
#
# UNITS (kept explicit, never folded into coefficients):
#   - artifact coefficients + returned I are in g*mm^2.
#   - the caller converts g*mm^2 -> N*m/(deg/s^2) via 1e-9 * (pi/180), exactly
#     as 02_deconvolve.R already does for the specimen MOI.

# Candidate fit forms, byte-for-byte mirror of bender_functions.py
# _APPARATUS_FIT_FORMS / _APPARATUS_FIT_TERMS (~line 6260). aor = axis-of-
# rotation-to-clamp distance (mm); width = clamp plate-to-plate span (mm).
.APPARATUS_FIT_TERMS <- list(
  F1 = c("a", "b"),
  F2 = c("a", "b", "c"),
  F3 = c("a", "b", "c"),
  F4 = c("a", "b", "c"),
  F5 = c("a", "b")
)

#' Design-matrix row for one (aor, width) under a fit form id. Mirrors
#' bender_functions.py _apparatus_fit_design_row() exactly.
#' @return numeric vector, same length/order as .APPARATUS_FIT_TERMS[[form_id]].
.apparatus_fit_design_row <- function(form_id, aor_mm, width_mm) {
  a <- as.numeric(aor_mm)
  w <- as.numeric(width_mm)
  switch(form_id,
    F1 = c(1.0, a * a),
    F2 = c(1.0, a, a * a),
    F3 = c(1.0, a * a, w),
    F4 = c(1.0, a * a, w * w),
    F5 = c(1.0, a * a + (w / 2.0)^2),
    stop("Unknown apparatus fit form: ", form_id)
  )
}

#' Resolve a referenced apparatus-inertia fit JSON to an absolute path.
#'
#' The trial file stores only a BASENAME in
#' calibration_inertia_apparatus_fit_file (e.g.
#' "2026-07-06_apparatus_inertia_calibration.json"); the artifact itself lives
#' in the repo (committed alongside the code). Resolution order:
#'   1. if fit_file is already an existing absolute/relative path, use it;
#'   2. $BENDER3_APPARATUS_FIT_DIR/<basename> (explicit override);
#'   3. <BENDER3_R_ROOT>/../<basename> (repo root, where the JSON is committed);
#'   4. <getwd()>/<basename>;
#'   5. alongside the raw .h5 file, if raw_path is supplied.
#' Returns the first hit, or NA_character_ if none exists (caller must warn).
resolve_apparatus_fit_path <- function(fit_file, raw_path = NULL) {
  if (is.null(fit_file) || is.na(fit_file) || !nzchar(fit_file)) return(NA_character_)
  if (file.exists(fit_file)) return(normalizePath(fit_file))

  base <- basename(fit_file)
  cands <- character(0)

  env_dir <- Sys.getenv("BENDER3_APPARATUS_FIT_DIR", unset = "")
  if (nzchar(env_dir)) cands <- c(cands, file.path(env_dir, base))

  r_root <- Sys.getenv("BENDER3_R_ROOT", unset = "R")
  cands <- c(cands, file.path(dirname(normalizePath(r_root, mustWork = FALSE)), base))
  cands <- c(cands, file.path(getwd(), base))
  if (!is.null(raw_path) && nzchar(raw_path)) cands <- c(cands, file.path(dirname(raw_path), base))

  hit <- cands[file.exists(cands)]
  if (length(hit) > 0L) normalizePath(hit[[1L]]) else NA_character_
}

#' Load + validate an apparatus-inertia fit artifact from JSON. Mirrors
#' bender_functions.py load_apparatus_inertia_fit(): rejects a wrong schema, a
#' blocked artifact, or one with no selected form/coefficients.
#' @return the parsed artifact (list), or stops with an actionable message.
read_apparatus_inertia_fit <- function(json_path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("read_apparatus_inertia_fit: the 'jsonlite' package is required")
  }
  if (is.na(json_path) || !file.exists(json_path)) {
    stop("read_apparatus_inertia_fit: artifact not found: ", json_path)
  }
  art <- jsonlite::fromJSON(json_path, simplifyVector = TRUE, simplifyDataFrame = FALSE)
  if (!is.list(art) || !identical(art$schema, "bender_apparatus_inertia_fit_v1")) {
    stop("read_apparatus_inertia_fit: not a bender_apparatus_inertia_fit_v1 artifact: ", json_path)
  }
  if (!is.null(art$blocked)) {
    stop("read_apparatus_inertia_fit: artifact is blocked (no usable fit): ", art$blocked)
  }
  if (is.null(art$fit_form_id) || !nzchar(art$fit_form_id) || is.null(art$fit_coefficients)) {
    stop("read_apparatus_inertia_fit: artifact has no selected fit form/coefficients.")
  }
  art
}

#' Evaluate apparatus MOI (g*mm^2) for one (aor, width) from a fit artifact.
#' Mirrors bender_functions.py apparatus_inertia_from_fit(~line 6661).
#'
#' @param artifact parsed artifact from read_apparatus_inertia_fit().
#' @param aor_mm,width_mm this trial's clamp geometry (mm). width_mm MUST be the
#'   clamp plate-to-plate span (attr calibration_inertia_apparatus_plate_to_plate
#'   _millimeter), i.e. the SAME quantity the calibration was fit against -- NOT
#'   the specimen body width (they differ: e.g. bass16 plate_to_plate 53 vs body
#'   width 42).
#' @param form_id optional override to evaluate a DIFFERENT candidate form than
#'   the artifact's selected one (diagnostic use: F4 vs F5). Coefficients are
#'   then taken from artifact$candidate_forms[[form_id]] when the override
#'   differs from the selected form.
#' @return list(i_gmm2, in_domain). in_domain is FALSE when (aor,width) is
#'   outside valid_domain; callers MUST warn and must not silently trust an
#'   extrapolated value.
apparatus_inertia_from_fit <- function(artifact, aor_mm, width_mm, form_id = NULL) {
  fid <- if (is.null(form_id)) artifact$fit_form_id else form_id
  if (is.null(fid) || !fid %in% names(.APPARATUS_FIT_TERMS)) {
    stop("apparatus_inertia_from_fit: unknown/absent fit form id: ", fid)
  }
  coef_src <- if (identical(fid, artifact$fit_form_id)) {
    artifact$fit_coefficients
  } else {
    cf <- artifact$candidate_forms[[fid]]
    if (is.null(cf) || is.null(cf$coefficients)) {
      stop("apparatus_inertia_from_fit: candidate form ", fid, " absent from artifact")
    }
    cf$coefficients
  }
  terms <- .APPARATUS_FIT_TERMS[[fid]]
  coef <- vapply(terms, function(nm) as.numeric(coef_src[[nm]]), numeric(1L))
  row  <- .apparatus_fit_design_row(fid, aor_mm, width_mm)
  i_gmm2 <- sum(row * coef)

  dom <- artifact$valid_domain
  a_rng <- if (!is.null(dom$aor_millimeter)) as.numeric(dom$aor_millimeter) else c(-Inf, Inf)
  w_rng <- if (!is.null(dom$width_millimeter)) as.numeric(dom$width_millimeter) else c(-Inf, Inf)
  in_domain <- is.finite(aor_mm) && is.finite(width_mm) &&
    aor_mm >= a_rng[1L] && aor_mm <= a_rng[2L] &&
    width_mm >= w_rng[1L] && width_mm <= w_rng[2L]

  list(i_gmm2 = i_gmm2, in_domain = as.logical(in_domain))
}
