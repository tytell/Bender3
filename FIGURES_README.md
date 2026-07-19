# Figure naming convention (adopted 2026-07-19)

Figures live outside the repo, on the shared OneDrive hub:
`proj_crittergripper/figures/`. This file documents the naming scheme so a
filename tells you what it is without opening it. It replaces ad hoc
suffixes (`_interpbaseline`, `_vector`, `_L0_sono`, `_snr_passing`) that had
grown inconsistently across scripts.

## Folder layout (unchanged)
- `bass{ID}_figures/trial_plots/` — one compound diagnostic plot per trial file.
- `bass{ID}_figures/summary_plots/` — per-fish pooled summaries across trials.
- `diagnostic_figures/` — cross-individual pooled superplots (flat, bassID
  list embedded in the name since it spans specimens).
- `bass_summary_figures/` — flattened, curated copies (only SNR-passing
  material) pulled from all three fish for side-by-side browsing; every
  filename here is prefixed with the bass ID since the folder mixes fish.

## Filename tokens
`{signal}_{protocol}[_{method}][_{filter}].png` inside per-fish folders;
`{bassID}_{signal}_{protocol}[_{method}][_{filter}].png` inside
`bass_summary_figures/`. Tokens always appear in this order; omit a
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

### Trial plots (`trial_plots/`, filename root = `{bassID}_bender_{NN}_{protocol}`, unchanged)
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

### Cross-individual superplots (`diagnostic_figures/`, flat)
| Old | New |
|---|---|
| `superplot_FL_pooled_bass16_17_18.png` | `FLsuperplot_isometric_isovelocity_pooled.png` |
| `superplot_FL_pooled_bass16_17_18_snr_passing.png` | `FLsuperplot_isometric_isovelocity_pooled_snrPass.png` |

### `bass_summary_figures/` (flat, curated, bassID-prefixed; only SNR-passing material)
Trial plots copied as-is (already `bassID_bender_NN_protocol[_baselineInterp].png`).
Legacy summary copies (no filtering possible) and regenerated SNR-passing
summaries both follow `{bassID}_{new-summary-name-above}[_snrPass].png`, e.g.
`bass18_FL_isometric_uhatBoth_snrPass.png`,
`bass17_uhatCompare_empVsGeom_snrPass.png`.

## Adding a new figure
Pick tokens from the lists above in order; add a new token to this file in
the same change if you introduce a new `signal`/`method`/`filter` value.
Don't invent a new suffix style — extend this scheme.
