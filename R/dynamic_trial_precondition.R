# dynamic_trial_precondition.R
# Specimen-specific HARD cutoff for dynamic-trial tissue preconditioning
# (decided 2026-07-22 -- see analysis_muscle_force_vector_log.md, "Dynamic-
# trial sono-strain gain traced to early-trial tissue preconditioning").
#
# Each specimen's earliest dynamic, ACTIVE (right-stim) trials show large
# excess-shortening bias -- sono-measured strain running far more negative
# than encoder-predicted strain -- that decays to near-zero within that same
# specimen's own session. This is consistent with a one-time tissue-
# preconditioning/settling effect, not a stable per-protocol artifact (the
# same degradation appears in isometric/isovelocity trials that happen to run
# EARLY in a session, e.g. bass18 trial 3, and is absent when those same
# protocols run LATER, after prior dynamic loading).
#
# The cutoff for each specimen is the first dynamic trial_num after which the
# active-window mean offset -- mean(strain_sono_pct) - mean(strain_pred_
# encoder_right_pct), i.e. a signed bias in percent of resting muscle length,
# NOT a correlation coefficient -- stays under 1.5 percentage points for
# EVERY later trial in that specimen. Offset (not pointwise r) was chosen
# deliberately: r is mathematically invariant to a constant additive bias, so
# it recovers to a "good" value BEFORE the mean offset itself has fully
# decayed (bass17: r > 0.85 by trial 6, offset < 1.5 pct-pt not until trial
# 9) -- offset is the stricter, safer gate for anything (e.g. power output)
# that depends on absolute strain accuracy rather than just tracking shape.
#
#   trial_num <  cutoff  ->  "early (preconditioning)"
#   trial_num >= cutoff  ->  "later (stable)"
#
# This ONLY labels dynamic, active(right-stim) trials. It does not touch
# isometric, isovelocity, frequency_sweep, or passive samples, and it does
# not delete or overwrite any existing data -- callers opt in by joining on
# specimen + trial_num (extract_bender_trial_num() below).

DYNAMIC_PRECONDITION_CUTOFF_TRIALNUM <- c(bass16 = 5L, bass17 = 9L, bass18 = 5L)

DYNAMIC_PRECONDITION_LEVELS <- c("early (preconditioning)", "later (stable)")

#' Parse the chronological bender_NN trial number out of a trial_id/filename
#' string, e.g. "2026-07-14_bass16_bender_05_dynamic" -> 5L. Every acquisition
#' run (dynamic, isovelocity, isometric, frequency_sweep) increments NN, so
#' this is a valid within-session chronological proxy even when applied only
#' to the "_dynamic" subset.
extract_bender_trial_num <- function(trial_id) {
  m <- regmatches(trial_id, regexpr("bender_(\\d+)", trial_id))
  suppressWarnings(as.integer(gsub("[^0-9]", "", m)))
}

#' Classify a dynamic trial as early(preconditioning)/later(stable) using the
#' hard specimen-specific cutoff above. Returns NA for specimens not in the
#' lookup or non-finite trial_num -- fails loud (NA), never silently assumes
#' a default cutoff for an unlisted specimen.
classify_dynamic_precondition <- function(specimen, trial_num) {
  cutoff <- DYNAMIC_PRECONDITION_CUTOFF_TRIALNUM[specimen]
  out <- ifelse(
    is.na(cutoff) | !is.finite(trial_num), NA_character_,
    ifelse(trial_num < cutoff, DYNAMIC_PRECONDITION_LEVELS[1], DYNAMIC_PRECONDITION_LEVELS[2])
  )
  factor(out, levels = DYNAMIC_PRECONDITION_LEVELS)
}

#' Same lookup/cutoff as classify_dynamic_precondition(), applied PROTOCOL-
#' AGNOSTICALLY (isometric/isovelocity/frequency_sweep trials included, not
#' just dynamic). This is deliberate, not a bug: trial_num is a per-specimen
#' SESSION-chronology counter that increments across every protocol type, and
#' the cutoff was derived (via dynamic-trial sono data) to mark a one-time
#' tissue-preconditioning/settling effect over the specimen's session -- a
#' property of the specimen's cumulative handling/loading history, not of the
#' dynamic protocol specifically. This is supported by the log's cross-
#' protocol observation (2026-07-22): isometric/isovelocity trials that
#' happen to run EARLY in a session (e.g. bass18 trial 3) show the SAME
#' degradation, and are clean when run later. Use this alias (rather than
#' calling classify_dynamic_precondition() directly) anywhere the exclusion
#' is being applied outside the dynamic protocol, so the intent is visible at
#' the call site.
classify_session_precondition <- function(specimen, trial_num) {
  classify_dynamic_precondition(specimen, trial_num)
}
