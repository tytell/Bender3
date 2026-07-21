# paths_config.R
# Single source of truth for every OneDrive-hosted path (raw .h5 sources +
# figure outputs) used across the R/ pipeline. Every script that needs one of
# these locations sources this file and calls its helpers -- NONE of them
# hardcode the OneDrive root independently. When the project folder structure
# moves (as it did 2026-07-19/2026-07-20: proj_crittergripper/figures/* ->
# proj_crittergripper/02_processed/figs_*), only THIS file needs to change.
#
# Every constant/helper can be overridden via BENDER3_* env vars (matching
# the existing BENDER3_SOURCE_DIR/BENDER3_OUTPUT_DIR override pattern already
# used by run_fv_fl_power_pipeline.R) so the pipeline still runs on another
# machine/mount point or against a relocated corpus without editing this
# file -- the hardcoded strings below are only the DEFAULTS.

.crittergripper_root <- function() {
  Sys.getenv(
    "BENDER3_CRITTERGRIPPER_ROOT",
    "/Users/yjimenez/Library/CloudStorage/OneDrive-ProvidenceCollege/01_JimenezLab/02_ResearchHub/proj_crittergripper"
  )
}

.permanent_archive_root <- function() {
  Sys.getenv(
    "BENDER3_ARCHIVE_ROOT",
    "/Users/yjimenez/Library/CloudStorage/OneDrive-ProvidenceCollege/01_JimenezLab/01_PermanentArchive/bender_crittergripper"
  )
}

#' Raw acquisition .h5 source folder for one dated trial block.
#'
#' bass16/17/18 (and other per-fish dated blocks) live under
#' 01_PermanentArchive/bender_crittergripper/{subfolder}/ -- this is the
#' archival source of truth and is NOT the same as
#' proj_crittergripper/01_inputs/data_raw/ (that flat folder holds a
#' different, not-yet-migrated subset of specimens: scup16, bass12,
#' snapper00, bass14, rod40a). Do not assume bass16/17/18 will appear there.
#'
#' @param subfolder Dated subfolder name, e.g. "2026-07-14_bass16_bender".
raw_source_dir <- function(subfolder) file.path(.permanent_archive_root(), subfolder)

# Known raw-source subfolders for the three vector-force specimens.
BASS16_RAW_SUBFOLDER <- "2026-07-14_bass16_bender"
BASS17_RAW_SUBFOLDER <- "2026-07-15_bass17_bender"
BASS18_RAW_SUBFOLDER <- "2026-07-16_bass18_bender"

#' Per-specimen figure output folder (holds trial_plots/ + summary_plots/
#' subdirs, created by the caller), e.g. figs_dir("bass16") ->
#' .../02_processed/figs_bass16.
figs_dir <- function(specimen_id) file.path(.crittergripper_root(), "02_processed", paste0("figs_", specimen_id))

# Cross-individual output folders (flat, not per-specimen).
# PI-directed rename 2026-07-21 (was figs_bass_summary) to match the
# .cursorrules "File placement" convention: cross-specimen/pooled summary
# figures -> 02_processed/figs_summary/.
FIGS_SUMMARY_DIR <- file.path(.crittergripper_root(), "02_processed", "figs_summary")
FIGS_DIAGNOSTIC_DIR   <- file.path(.crittergripper_root(), "02_processed", "figs_diagnostic")
# Test-run output only (not a deliverable -- see .gitignore); was
# figures/tests/ before the 2026-07-19/20 reorg, now under 02_processed/.
FIGS_TESTS_DIR        <- file.path(.crittergripper_root(), "02_processed", "figs_tests")

# Unrelated to the crittergripper reorg above -- incoming-file staging area,
# kept here only because plot_torque_vs_time_batch.R reads from it.
DROPZONE_DIR <- Sys.getenv(
  "BENDER3_DROPZONE",
  "/Users/yjimenez/Library/CloudStorage/OneDrive-ProvidenceCollege/01_JimenezLab/00_DropZone"
)
