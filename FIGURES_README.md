# Figure naming convention (adopted 2026-07-19)

Figures live outside the repo, on the shared OneDrive hub:
`proj_crittergripper/02_processed/`. This file documents the naming scheme so
a filename tells you what it is without opening it. It replaces ad hoc
suffixes (`_interpbaseline`, `_vector`, `_L0_sono`, `_snr_passing`) that had
grown inconsistently across scripts.

All folder paths below are single-sourced in `R/paths_config.R` -- if the
OneDrive layout moves again, update that file only (this doc + the
`R/*.R` scripts read from it, not from hardcoded strings).

## Folder layout
Renamed 2026-07-19/2026-07-20/2026-07-21 to match the `.cursorrules` "File
placement" convention (`figures/` -> `02_processed/figs_*`;
`bass_summary_figures/` -> `figs_summary/`, matching the "cross-specimen or
pooled summary figures" rule). Flattened 2026-07-21 (PI-directed): every
folder below is FLAT, no subfolders anywhere -- `trial_plots/`/
`summary_plots/` per-fish subfolders were removed; trial-plot and
summary-plot filenames already have distinct naming shapes
(`{bassID}_bender_{NN}_{protocol}...` vs. `{signal}_{protocol}...`), so
merging them into one folder creates no collisions.
- `figs_{bassID}/` — flat: one compound trial plot per trial file
  (`{bassID}_bender_{NN}_{protocol}[_{method}].png`) PLUS that fish's
  pooled summary plots, side by side in the same folder.
- `figs_diagnostic/` — flat: cross-individual pooled superplots, AND
  (2026-07-21) diagnostic/decision-making plots that compare filters or
  calculations rather than reporting a final result (see "Diagnostic vs.
  individual vs. summary" below) -- topic carried in the filename, not a
  subfolder, since these plots are per-decision, not per-fish.
- `figs_summary/` — flattened, curated copies (only SNR-passing
  material) pulled from all three fish for side-by-side browsing; every
  filename here is prefixed with the bass ID since the folder mixes fish.

## Diagnostic vs. individual vs. summary (adopted 2026-07-21)
Three tiers, not just two -- see `analysis_muscle_force_vector_log.md` for
the point-selection design this feeds:
- **Diagnostic** (`figs_diagnostic/`, flat, filename-tagged by topic, not
  per-fish): plots that compare filters/calculations to make a decision
  (e.g. empirical vs. geometric u_hat, legacy vs. vector force, own-step
  vs. interpolated baseline, smoothing cutoffs). Once a decision is locked,
  it's recorded in a manifest (`00_records/`, not yet built) rather than
  re-litigated per plot. Current tokens: `musclepullmethod`,
  `muscleforceestimate`, `passivebaselinemethod`, `torquesmoothingmethod`,
  `sonosmoothingmethod`, `sonotiminglag`, `forcedevtiming` (dual-tagged,
  see below), `fatiguetimeline` (dual-tagged), `snrfiltereffect`.
- **Individual** (`figs_{bassID}/`): trial-level plots, PLUS per-fish
  aggregated plots (all trials for that fish, trial identity kept visible
  via color/shape) using whichever method won its diagnostic decision.
  `forcedevtiming` and `fatiguetimeline` are dual-tagged -- diagnostic
  output that is ALSO a legitimate individual-tier visualization (within-
  step force development timing; near-L0 force vs. real elapsed session
  time), so they may appear in both places.
- **Summary** (`figs_summary/`): side-by-side per-individual panels (not
  pooled into one panel), raw / individual-mean / group-mean tiers.
  Individual-mean definition is still open (point-selection method design,
  see the decision log) -- do not assume a fitted curve vs. binned mean
  without checking that file first.

## Filename tokens
`{signal}_{protocol}[_{method}][_{filter}].png` inside per-fish folders;
`{bassID}_{signal}_{protocol}[_{method}][_{filter}].png` inside
`figs_summary/`. Tokens always appear in this order; omit a
bracketed token when it doesn't apply (e.g. `legacy` method is implicit and
dropped in a few original-workflow names — see table below for the exact
names in use).

- **protocol**: `isometric` | `isovelocity` | `dynamic` | `freqsweep`
- **signal**: what's plotted — `FL`, `FV`, `FVl0` (FV sampled at sono-confirmed
  L0 crossing), `forceTime`, `uhatCompare`, `fatigueCheck`, `angleValid`,
  `strainValidCmd`/`strainValidMeas`/`strainValidSonoEnc`/`strainValidSonoCmd`,
  `powerDynamic`/`powerDynamicMassSpec`, `stiffnessDamping`, `FLsuperplot`
- **method**: which force/baseline calculation —
  - `legacy` = original single-axis (zTorque only) calculation, no û
  - `baselineInterp` = isometric per-step baseline computed via
    cross-step interpolation instead of that step's own pre-stim window
  - `uhatEmp` / `uhatGeom` = 6-axis vector force via empirical /
    geometric line-of-action
  - `uhatBoth` = both û methods faceted side by side in one file
- **filter**: `snrPass` = regenerated with sub-threshold
  (`activation_snr < 3`) points/steps dropped entirely. Omitted = unfiltered
  (low-confidence points are still alpha-flagged in-plot where applicable,
  not dropped).

## Old name -> new name
Historical migration record (2026-07-19 renaming pass) -- the `trial_plots/`/
`summary_plots/` subfolder references below describe the layout AS IT WAS
THEN, not the current one (flattened 2026-07-21, see "Folder layout" above).
The filename tokens themselves are still current.

### Trial plots (was in `trial_plots/`, filename root = `{bassID}_bender_{NN}_{protocol}`, unchanged)
| Old | New |
|---|---|
| `{tid}.png` | unchanged |
| `{tid}_interpbaseline.png` | `{tid}_baselineInterp.png` |

### Summary plots (`summary_plots/`, per-fish, no bassID prefix needed — folder scopes it)
| Old | New |
|---|---|
| `summary_isometric_FL.png` | `FL_isometric_legacy.png` |
| `summary_isometric_FL_interpbaseline.png` | `FL_isometric_baselineInterp.png` |
| `summary_isometric_FL_vector.png` | `FL_isometric_uhatBoth.png` |
| `summary_isovelocity_FV.png` | `FV_isovelocity_legacy.png` |
| `summary_isovelocity_FV_vector.png` | `FV_isovelocity_uhatBoth.png` |
| `summary_isovelocity_FV_vector_L0_sono.png` | `FVl0_isovelocity_uhatBoth.png` |
| `uhat_empirical_vs_geometric.png` | `uhatCompare_empVsGeom.png` |
| `force_vs_time_isometric.png` | `forceTime_isometric_legacy.png` |
| `force_vs_time_isovelocity.png` | `forceTime_isovelocity_legacy.png` |
| `force_vs_time_isometric_vector.png` | `forceTime_isometric_uhatBoth.png` |
| `force_vs_time_isovelocity_vector.png` | `forceTime_isovelocity_uhatBoth.png` |
| `force_vs_time_dynamic_by_{bname}.png` | `forceTime_dynamic_by_{bname}.png` (bname = frequency/amplitude/phase/duty) |
| `fatigue_check_isometric.png` | `fatigueCheck_isometric_legacy.png` |
| `fatigue_check_isovelocity.png` | `fatigueCheck_isovelocity_legacy.png` |
| `fatigue_check_isometric_interpbaseline.png` | `fatigueCheck_isometric_baselineInterp.png` |
| `strain_validation_commanded.png` | `strainValidCmd.png` |
| `strain_validation_measured.png` | `strainValidMeas.png` |
| `strain_validation_sono_vs_encoder.png` | `strainValidSonoEnc.png` |
| `strain_validation_sono_vs_commanded.png` | `strainValidSonoCmd.png` |
| `angle_validation.png` | `angleValid.png` |
| `summary_dynamic_power.png` | `powerDynamic.png` |
| `summary_dynamic_power_massspecific.png` | `powerDynamicMassSpec.png` |
| `summary_frequency_sweep_stiffness_damping.png` | `stiffnessDamping_freqsweep.png` |

### Cross-individual superplots (`figs_diagnostic/`, flat)
| Old | New |
|---|---|
| `superplot_FL_pooled_bass16_17_18.png` | `FLsuperplot_isometric_isovelocity_pooled.png` |
| `superplot_FL_pooled_bass16_17_18_snr_passing.png` | `FLsuperplot_isometric_isovelocity_pooled_snrPass.png` |

### `figs_summary/` (flat, curated, bassID-prefixed; only SNR-passing material)
Trial plots copied as-is (already `bassID_bender_NN_protocol[_baselineInterp].png`).
Legacy summary copies (no filtering possible) and regenerated SNR-passing
summaries both follow `{bassID}_{new-summary-name-above}[_snrPass].png`, e.g.
`bass18_FL_isometric_uhatBoth_snrPass.png`,
`bass17_uhatCompare_empVsGeom_snrPass.png`.

## Adding a new figure
Pick tokens from the lists above in order; add a new token to this file in
the same change if you introduce a new `signal`/`method`/`filter` value.
Don't invent a new suffix style — extend this scheme.
