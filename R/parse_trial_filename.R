# parse_trial_filename.R
# Filename parser for the LOCKED convention documented in
# context_jlab_cg_h5schema.md SS2 "Filename convention":
#   YYYY-MM-DD_<specimen_id>_bender_<NN>_<protocol>.h5
#   sim_YYYY-MM-DD_<specimen_id>_bender_<NN>_<protocol>.h5  (simulation-mode)
#
# NOT a reuse of R/00_load_bender_flat.R's .bfl_trial_from_filename(): that
# helper's regex assumes the trial number is the LAST underscore-delimited
# token before ".h5" (`..._bender_NN.h5`), which does not hold for this
# corpus -- every real filename here has a <protocol> suffix after <NN>
# (e.g. "..._bender_20_frequency_sweep.h5"), so that helper resolves trial
# to NA for all 21 files. Verified against the actual source directory
# listing (2026-07-15) before writing this parser, per task instructions.

library(stringr)

#' Parse one Bender raw-data filename into its documented components.
#'
#' @param path File path or bare filename.
#' @return A tibble (1 row) with session_date, specimen_id, trial_number,
#'   protocol, is_sim, trial_id (unique, human-readable: "<specimen>_bender_<NN>_<protocol>"),
#'   and ok (FALSE if the filename did not match the documented convention).
parse_bender_filename <- function(path) {
  fn <- basename(path)
  fn_noext <- sub("\\.h5$", "", fn, ignore.case = TRUE)

  is_sim <- grepl("^sim_", fn_noext)
  core <- sub("^sim_", "", fn_noext)

  # YYYY-MM-DD_<specimen_id>_bender_<NN>_<protocol>
  m <- str_match(core, "^(\\d{4}-\\d{2}-\\d{2})_([^_]+)_bender_(\\d+)_(.+)$")

  if (is.na(m[1, 1])) {
    return(tibble::tibble(
      filename = fn, session_date = NA_character_, specimen_id = NA_character_,
      trial_number = NA_integer_, protocol = NA_character_, is_sim = is_sim,
      trial_id = fn_noext, ok = FALSE
    ))
  }

  tibble::tibble(
    filename     = fn,
    session_date = m[1, 2],
    specimen_id  = m[1, 3],
    trial_number = as.integer(m[1, 4]),
    protocol     = m[1, 5],
    is_sim       = is_sim,
    trial_id     = sprintf("%s_bender_%02d_%s", m[1, 3], as.integer(m[1, 4]), m[1, 5]),
    ok           = TRUE
  )
}

#' Parse and validate a full directory of trial files; aborts loudly (rather
#' than silently dropping files) if any .h5 file fails to match the
#' documented convention or if trial_id collides across files.
parse_trial_directory <- function(dir_path, pattern = "\\.h5$") {
  files <- fs::dir_ls(dir_path, glob = paste0("*", sub("^\\\\", "", pattern)))
  files <- files[!grepl("\\.png$|\\.heic$", files, ignore.case = TRUE)]
  if (length(files) == 0L) cli::cli_abort("parse_trial_directory: no .h5 files found in {dir_path}")

  parsed <- purrr::map_dfr(files, parse_bender_filename)
  parsed$fullpath <- as.character(files)

  bad <- dplyr::filter(parsed, !.data$ok)
  if (nrow(bad) > 0L) {
    cli::cli_abort(c(
      "parse_trial_directory: {nrow(bad)} filename(s) did not match the documented convention:",
      setNames(bad$filename, rep("x", nrow(bad)))
    ))
  }

  dupes <- parsed$trial_id[duplicated(parsed$trial_id)]
  if (length(dupes) > 0L) {
    cli::cli_abort("parse_trial_directory: non-unique trial_id(s): {paste(unique(dupes), collapse=', ')}")
  }

  parsed
}
