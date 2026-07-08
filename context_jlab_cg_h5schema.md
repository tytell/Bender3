# Bender / CritterGripper — HDF5 Schema

## Authority metadata

- **Purpose:** Defines the on-disk structure, group layout, key-naming taxonomy, and unit conventions for Bender HDF5 files.
- **Authority level:** Tier 2 (binding for all H5 write/read tasks — export routine, R pipeline, GUI→backend assignment).
- **Scope:** `bender_h5_export.py` (write), `01_calibration.R` / `02_correct.R` / `03_analyze.R` (read), GUI variable assignment in `bender_streamlit_gui.py`.
- **Owner:** PI (schema decisions are PI-only).
- **Version:** 2026-07-02 (amended 2026-07-02: complete M1 metadata contract — `schema_version` now ISO-8601 date, `note_bench` always written, `calibration_inertia_file` always written, `calibration_inertia_matrix` present-but-empty)
- **Last reviewed:** 2026-07-02

### Changelog


| Version        | Date       | Summary                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| -------------- | ---------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| v2.8           | 2026-07-08 | Apparatus-inertia fit provenance: exporter direct-writes three `calibration_inertia_apparatus_fit_*` keys from the Bender attrs `apparatus_inertia_calibration` (dict) + `apparatus_inertia_calibration_file` (path). `_fit_file` = provenance path, `_fit_form` = selected form equation, `_fit_json` = the FULL fit artifact serialized as a JSON string DATASET (coefficients, source_files, metrics, valid_domain, geometry_check, sign note, per-trial rows). Embeds VALUES (self-contained re-processing), mirroring the F/T calibration matrix, rather than a filename pointer. Fit built by `fit_apparatus_inertia.py`, loaded via the GUI Hardware-config section. String scalars follow the null-sentinel rule ("NA" when no fit loaded). Routing: the two Bender attrs are EXCLUDED, the three keys are MISSING_REQUIRED (computed, direct-write). |
| null-sentinel  | 2026-07-02 | Fixed-schema null-sentinel rule (see section 11): every canonical metadata scalar is ALWAYS emitted; absent values use a typed NA sentinel (float->NaN, string->"NA", int->float+NaN, bool->3-state string "true"/"false"/"NA", array->empty). Replaces omit-on-zero / skip-on-None. `required` fields must be present AND non-NA. Timeseries channels are exempt (absent = absent/empty dataset, never NaN-filled). Routing renames: `xsec_width` -> `measurement_specimen_local_body_width_millimeter`, `xsec_height` -> `measurement_specimen_local_body_height_millimeter` (retired `measurement_xsec_width_millimeter` / `measurement_xsec_height_millimeter`). |
| v2.7 (M2c)     | 2026-07-02 | D10 implementation (M2c): exporter writes the `single_finite` per-cycle design grid `index_cycle_*` (index/frequency/motor_amplitude/active/duty/phase/operating_point + units) from the `organize_cycles` per-cycle arrays, plus scalars `protocol_cycle_count` and `protocol_activation_start_cycle` (0-based). `segmented_finite` writes none of these. Routing key `n_end_cycles` renamed `protocol_end_cycles_count` → `protocol_end_cycle_count` to match this section. `index_step_*` untouched (D4 intact). |
| v2.7 (M2b)     | 2026-07-02 | D11 implementation (M2b): `forcetorque_raw` is NEVER written to `timeseries` (unconditional; simplified from the fallback-guard design). `calibration_forcetorque_matrix` is embedded ONLY when a REAL ATI matrix is loaded — the `np.eye(6)` identity fallback is REFUSED, so an absent matrix means "not calibrated" and `derived/forcetorque_calibrated` is likewise skipped. `derived/forcetorque_calibrated` re-anchored to `aidata` rows 0–5 via `daq_ai_channel_map`. Sono cal keys renamed to `calibration_sono_{left,right}_volt_millimeter_breakpoints`. Readers (`validate_plot_h5_batch.py`, `bender_h5_explore.py`) point at `derived/forcetorque_calibrated`. |
| v2.7           | 2026-06-30 | Frustum model raw inputs captured as first-class scalars (PI decision 2026-06-30). `tip_scale` DROPPED.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    |
| v2.7 (amended) | 2026-07-01 | (1) Rename four frustum scalars: `measurement_specimen_frustum_{height,width,length,density_*}` → `measurement_specimen_inertia_frustum_*` (adds "inertia" segment for consistency with `calibration_inertia_*` namespace). (2) Reclassify `clamp_offset_mm`/`specimen_profile_clamp_offset_mm` from draft `measurement_specimen_frustum_clamp_offset_millimeter` / `measurement_profile_clamp_offset_millimeter` (both retired) to `calibration_inertia_apparatus_aor_to_clamp_millimeter` (apparatus geometry, not specimen; many-to-one, first positive value wins, omit if absent). (3) Reclassify `clamp_plate_extension_mm` from `measurement_clamp_plate_extension_millimeter` (retired) to `calibration_inertia_apparatus_plate_to_plate_millimeter` (apparatus geometry; feeds Deferred flag 4). (4) `_hardware_inertia_baseline` / `CLAMP_DIM` superseded and deleted: `calibration_inertia_total_moi_gram_millimeter_squared` is specimen-MOI-only until Deferred flag 4 calibration-matrix workflow provides per-configuration apparatus term. `derived/torque_inertia_corrected` uses empirical I_est + specimen MOI (unchanged). (5) Governing Rule 10 added: every user-facing GUI input field must have a corresponding exported HDF5 key. |
| v2.6           | 2026-06-22 | Finalized Points 1-4 (folds in the standalone inertia-decisions plan). Point 1: raw geometry/mass/density as `measurement_specimen_*`, specimen MOI computed inline (analytic N-station integral, no `num_samples`), no provenance blob, drop `use_theoretical_inertial_correction`/`use_frustum_inertial_model` from schema. Point 2: `inertial_calibration_profile` flattened — `I_est` -> `calibration_inertia_apparatus_moi_gram_millimeter_squared` (`*(180/pi)*1e9`), `bias_est` -> `calibration_inertia_bias_newton_meter`, `axis_sensor` -> `calibration_inertia_axis_sensor`. Point 3: `calibration_inertia_system_moi` DROPPED (phantom); "apparatus" = I_est; `total_moi` confirmed = theoretical apparatus + specimen (Deferred flag 5). Point 4: exporter writes inspection-only `derived/forcetorque_calibrated` + `derived/torque_inertia_corrected` (R ignores). MOI keys carry `_gram_millimeter_squared`; no `attrs['units']`.                                                                                                                                                                                                                                                                                                           |
| v2.5           | 2026-06-22 | D10: `index_cycle_*`/`protocol_*` per-cycle design-grid arrays for `single_finite` (key names locked at M0); D11: exporter embeds calibration VALUES in `metadata`, removes calibrated `forcetorque_raw` from `timeseries` when matrix present (fallback guard), adds `calibration_forcetorque_file` provenance — supersedes §3d + Governing Rule 2 + F/T channel-order lock; D12: all MOI under `calibration_inertia_*`, drop `inertial_` prefix section, drop `calibration_inertia_used`/`_available`; D8 amendment: `03_analyze.R` is a PORT of validated corpus (not fresh); "inertial" → "inertia" normalization in all field names                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   |
| v2.4           | 2026-06-19 | Phase-1 spec freeze: flat `metadata`/`timeseries` groups, `calibration_inertial_*` flat keys, complete `index_step_*` column map, `metadata/99_Unrouted`, root-attr reduction to `schema_version` only (decisions D1-D9 per migration plan)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| v2.3           | 2026-06-15 | Phase 0: ledger-driven canonical exporter; `forcetorque_raw` `[6×N]` unsplit; per-sample index/stim arrays; `trial_index`/`step_order` dropped; `protocol_metadata` removed; `protocol_block_sequence` JSON; `99_Unrouted` fallback                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |


---

## PRECEDENCE DECLARATION (read first)

**This schema supersedes the prior `tier2_data_curation_sop` / `jlab_folder_architecture` self-documenting field-name requirements.** Previously-binding fields (`apparatus_id`, `instrumentation`, `specimen_id`, `session_date`, `protocol`, `analyst`) are folded into the prefix taxonomy below. Where this schema and any earlier SOP conflict, **this schema wins.**

Consequence (logged as backlog, not a blocker): the data-curation SOP and any provenance assertions referencing old key names must be updated to the new key names. See §8.

---

## 1. Variable naming convention (LOCKED)

**Pattern: `<description>_<unit>`**

- All lowercase. Underscores only — no hyphens, spaces, or special characters.
- Units are **fully spelled out** (no abbreviations). This is a locked decision.
- Products of units: join with underscore (`newton_meter` = N·m).
- Rates / ratios: use `_per_` (`degree_per_second`, `newton_per_millimeter`).
- Identifiers with no unit get **no suffix** (`specimen_id`, `segment`, `prep_condition`).
- Arrays keep their natural name; the unit suffix describes **one element** (`frustum_heights_millimeter` = array of heights, each in mm).

### Canonical unit spellings


| Physical quantity   | Spelled-out unit                                                          | Example field                 |
| ------------------- | ------------------------------------------------------------------------- | ----------------------------- |
| Length / distance   | `millimeter`                                                              | `clamp_separation_millimeter` |
| Mass                | `gram`                                                                    | `body_mass_gram`              |
| Temperature         | `celsius`                                                                 | `temperature_celsius`         |
| Angle               | `degree`                                                                  | `angle_commanded_degree`      |
| Angular velocity    | `degree_per_second`                                                       | `velocity_degree_per_second`  |
| Strain              | `percent` (positive = lengthening, negative = shortening)                 | `strain_percent`              |
| Time                | `second`                                                                  | `time_second`                 |
| Frequency           | `hertz`                                                                   | `sample_rate_hertz`           |
| Force               | `newton`                                                                  | `force_x_newton`              |
| Torque              | `newton_meter`                                                            | `torque_z_newton_meter`       |
| Stiffness           | `newton_per_millimeter` (or context-appropriate spelled-out derived unit) | —                             |
| Voltage             | `volt`                                                                    | `stim_command_volt`           |
| Dimensionless ratio | no suffix                                                                 | `time_normalized`             |


> **ATI force/torque units — LOCKED: newton (force) / newton_meter (torque).** `FT56491.cal` outputs SI units. All decoded channel names and calibration metadata keys use `newton` / `newton_meter`.

---

## 2. Top-level structure

Two groups in the raw file. **No numeric prefixes.**

```
<file>.h5
├── metadata      # flat, prefixed scalars/small arrays — see §4
└── timeseries    # continuous sample-indexed streams — see §3
```

The **research-hub copy** adds a third group (never written to the raw file):

```
└── derived       # corrected/computed outputs — HUB COPY ONLY
```

Group names are **identical** in both the raw file and the hub copy.

### Filename convention

Raw files follow the pattern:

```
YYYY-MM-DD_<specimen_id>_bender_<NN>_<protocol>.h5
```

Simulation-mode files are prefixed `sim_` and carry `metadata/session_simulated = True`:

```
sim_YYYY-MM-DD_<specimen_id>_bender_<NN>_<protocol>.h5
```

Analysis-pipeline rule: `filename starts with sim_` **OR** `session_simulated = True` → quarantine; absent flag → unknown/legacy.

### Root-level attributes (D6 — locked)

After Phase 1 migration (M1) only **one** attribute survives at the HDF5 file root:


| Root attr        | Value        | Notes                                                                                        |
| ---------------- | ------------ | -------------------------------------------------------------------------------------------- |
| `schema_version` | e.g. `"2026-07-02"` | File-level schema provenance; kept at root because it predates and outlives any group rename |


All other previous root attrs are moved to flat `metadata/` keys:


| Previous root attr | Disposition    | Target key                                                        |
| ------------------ | -------------- | ----------------------------------------------------------------- |
| `test_type`        | MOVE           | `metadata/protocol_type`                                          |
| `post-trial notes` | MOVE           | `metadata/note_bench` (+ `note_posthoc`)                          |
| `start_time_iso`   | MOVE           | `metadata/session_date`                                           |
| `filename`         | DROP from root | `metadata/filename` (already mirrored there; root copy redundant) |


After M1 the root contains exactly: `schema_version` attr + `metadata` group + `timeseries` group.
M1 must print the actual surviving root attr set during a sim export for PI eyeball before finalizing.

---

## 3. `timeseries`

Continuous streams, sample-indexed at the acquisition rate, sharing one `time_second` axis.

> **Phase 1 migration target (build target from v2.5).** Dataset names below are the canonical target. The on-disk group skeleton migrates from `02_TimeSeries/trial_XXXX` (zero-based, numeric prefix) to flat `timeseries/` channel arrays for `single_finite`, and `timeseries/step_NNN/` subgroups (one-based, 3-digit) for `segmented_finite`. The `01_`/`02_` numeric prefixes are dropped — hard cut, no dual-name reader (D2). The `index_step_*` flat parallel arrays (§4 `index_` section) replace the old `trial_index` subgroup in the same M2 commit. Per-sample index streams (§3b-i) concatenate cleanly into the flat timeline once the skeleton is collapsed.

### 3a. Motor and time axes


| Dataset                                        | Shape | Notes                                                                               |
| ---------------------------------------------- | ----- | ----------------------------------------------------------------------------------- |
| `time_second`                                  | `[N]` | Time axis in seconds. *(Was `t`; renamed in Phase 0.)*                              |
| `time_normalized`                              | `[N]` | Normalized time axis, dimensionless. *(Was `tnorm`.)*                               |
| `angle_commanded_degree`                       | `[N]` | Commanded motor angle (degrees). *(Was `angle_cmd`.)*                               |
| `angle_measured_degree`                        | `[N]` | Encoder angle, gearbox output / specimen frame (degrees). *(Was `angle_measured`.)* |
| `angular_velocity_commanded_degree_per_second` | `[N]` | Commanded angular velocity. *(Was `anglevel_cmd`.)*                                 |


### 3b. Stimulator command streams


| Dataset                         | Shape | Notes                                                                                     |
| ------------------------------- | ----- | ----------------------------------------------------------------------------------------- |
| `stim_channel1_command_volt`    | `[N]` | Stim command waveform, channel 1 (volts). *(Was `S1stimcmd`.)*                            |
| `stim_channel2_command_volt`    | `[N]` | Stim command waveform, channel 2 (volts). *(Was `S2stimcmd`.)*                            |
| `stim_side`                     | `[N]` | Per-sample categorical: `none` / `left` / `right` / `both` (handles co-contraction).      |
| `stim_state`                    | `[N]` | Per-sample stim phase: `passive` / `off` / `on` (3-value string).                         |
| `stim_type`                     | `[N]` | Per-sample activity: `active` / `passive`.                                                |
| `instantaneous_frequency_hertz` | `[N]` | Chirp instantaneous frequency (frequency_sweep only). *(Was `sweep_instantaneous_freq`.)* |


> **Stim fields — RESOLVED (PI).** `stim_state` stays a 3-value string (`passive`/`off`/`on`). `stim_type` stays per-sample (`active`/`passive`) — `active_passive` mixes both within one experiment. `stim_side` already exists as the per-sample categorical above (`none`/`left`/`right`/`both`). No code change.

### 3b-i. Per-sample index streams (LOCKED — converted from per-trial scalar attrs)

Each is a length-`[N]` array; within the current `trial_XXXX` skeleton the step-protocol fields are constant across a trial's samples (broadcast) and concatenate cleanly when the timeline is later flattened.


| Dataset            | Shape | Notes                                                                                                              |
| ------------------ | ----- | ------------------------------------------------------------------------------------------------------------------ |
| `cycle_index`      | `[N]` | Cycle number per sample; `-1` = not a numbered cycle (all step-protocol samples). *(Was `cycle_index_by_sample`.)* |
| `step_index`       | `[N]` | Ramp/step number at each sample (step protocols only).                                                             |
| `sequence_index`   | `[N]` | Pre-shuffle parameter index at each sample (step protocols only).                                                  |
| `block_index`      | `[N]` | Block number at each sample (step protocols only).                                                                 |
| `block_direction`  | `[N]` | Bending direction for the block (string; step protocols only).                                                     |
| `block_stim_sides` | `[N]` | Stim sides for the block (string; step protocols only).                                                            |


**Dropped:** `trial_index` (redundant with the file/trial identity + `sequence_index`) and `step_order` (reconstructable from per-sample `sequence_index` + `block_index`). Dynamic / frequency_sweep carry only `cycle_index` (no discrete shuffled steps/blocks).

### 3c. Raw analog input buffer


| Dataset  | Shape     | Notes                                                                                                                                                            |
| -------- | --------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `aidata` | `[8 × N]` | Raw analog input buffer. **IMMUTABLE** — written as-acquired, never decoded or modified at write time. Channel identities live in `metadata/daq_ai_channel_map`. |


### 3d. ATI 6-axis force/torque (D11 — supersedes prior §3d, Governing Rule 2, and F/T channel-order lock)

> **D11 (PI sign-off 2026-06-22) SUPERSEDES this section, Governing Rule 2, and the F/T channel-order lock.**
>
> **Architecture:** The raw file is a raw-voltage archive. The exporter does **not** compute or write a calibrated F/T array. Raw `aidata` rows 0–5 (F/T channels, identified via `metadata/daq_ai_channel_map`) are the immutable archive. The ATI 6×6 calibration matrix VALUES and sono calibration table VALUES are embedded in `metadata/` at collection time (see `calibration_forcetorque_matrix`, `calibration_forcetorque_file` in §4). Calibrated F/T output is computed in R and written only to hub `derived/`.
>
> **Removal of `forcetorque_raw` (M2b — simplified, no branching):** The `forcetorque_raw [6 × N]` timeseries dataset is **never written** by the exporter. `timeseries` holds only the immutable raw `aidata` archive (F/T = rows 0–5, identified via `daq_ai_channel_map`). The calibrated F/T copy lives ONLY in `derived/forcetorque_calibrated`, and that is written ONLY when a REAL `calibration_forcetorque_matrix` is embedded (identity fallback refused, see below). No real matrix → no embedded matrix AND no calibrated output. (Legacy pre-M2b files may still contain a `forcetorque_raw`; readers keep a fallback for those, but new writes never emit it.)
>
> **F/T channel-order re-anchor (D11):** Channel identity (row index → physical quantity) is a `daq_ai_channel_map` metadata fact, not a hard-coded R slice rule. Per-channel order (`0=Fx, 1=Fy, 2=Fz, 3=Tx, 4=Ty, 5=Tz`, newton / newton_meter) is declared in `daq_ai_channel_map` and read from there. The prior "LOCKED" rule asserting that the R pipeline slices `forcetorque_raw` directly by row index is **superseded**; R reads `daq_ai_channel_map` to identify channels.

**Not in `timeseries`:** calibrated F/T arrays, inertial-corrected torque series, `primary_torque_*`, and any computed derived quantities. All correction and calibration happens in R and lands in hub `derived/` only.

### 3e. Timeseries group structure (by `protocol_sampling_mode`)

`**single_finite`** — dynamic, sweep:

- `timeseries/` contains flat channel arrays directly.
- One continuous monotonic `time_second` axis.
- No subgroups.

`**segmented_finite`** — isometric, isovelocity, FL, FV:

- `timeseries/` contains subgroups `step_001/`, `step_002/`, `step_003/`, … (ONE-BASED, zero-padded 3-digit index; first step is `step_001`).
- Each subgroup holds the same channel arrays as the flat structure.
- `time_second` within each subgroup is **local** to that step, starting at 0.
- Subgroup prefix is always `step_` — never `trial_` or any other prefix.

The locked protocol→mode mapping is in §10.

---

## 4. `metadata`

Flat. **No subgroups.** Every key carries a descriptive prefix; alphabetical sort groups related fields.

### `calibration_` — sensor → physical-unit conversions and inertia parameters

```
calibration_forcetorque_matrix          [6×6]   ATI calibration matrix VALUES (raw ADC → newton / newton_meter);
                                                  embedded at collection (D11). Written ONLY when a REAL matrix is
                                                  loaded — the np.eye(6) identity fallback is REFUSED, so an ABSENT
                                                  key unambiguously means "not calibrated" (M2b). (Current Python
                                                  attr: calibration)
calibration_forcetorque_file            str     ATI calibration file name/path (provenance string); D11.
calibration_sono_left_volt_millimeter_breakpoints    [4 or 8]  Sono V→mm breakpoint table, left channel
                                                  [Low_V, High_V, Low_mm, High_mm] (or longer multi-point); VALUES
                                                  embedded at collection (D11, M2b). (Current: sono_cal_left)
calibration_sono_right_volt_millimeter_breakpoints   [4 or 8]  Sono V→mm breakpoint table, right channel; VALUES
                                                  embedded at collection (D11, M2b). (Current: sono_cal_right)

calibration_inertia_file                                  str    Path/name of the inertia calibration file; empty string if none
                                                                   (D3, D12 — flattened from calibration_link/calibration_file in M1).
calibration_inertia_apparatus_moi_gram_millimeter_squared float  Empirical apparatus (rig-alone) MOI in g*mm^2 from the empty-run
                                                                   calibration fit. Source is I_est (fit as torque[N*m] = I_est*alpha
                                                                   + bias with alpha in deg/s^2, so I_est is N*m/(deg/s^2)); converted
                                                                   at write time by *(180/pi) -> kg*m^2 then *1e9 -> g*mm^2 (D12 final,
                                                                   Point 2/3; replaces inertial_system_moi AND the raw I_est that lived
                                                                   in the inertial_calibration_profile subgroup).
calibration_inertia_bias_newton_meter                     float  Torque-offset term (bias_est, N*m) from the same fit; NOT an MOI, no
                                                                   unit conversion (D12 final, Point 2; flattened from
                                                                   inertial_calibration_profile/bias_est).
calibration_inertia_axis_sensor                           str    Sensor-frame torque axis the fit used (axis_sensor; categorical, no
                                                                   units) (D12 final, Point 2; flattened from
                                                                   inertial_calibration_profile/axis_sensor).
calibration_inertia_specimen_moi_gram_millimeter_squared  float  Specimen MOI in g*mm^2, computed inline at export from the N-station
                                                                   elliptical-frustum integration (set_specimen_geometry_inertial_model;
                                                                   exact per-segment integral, NO sample-count knob — the station
                                                                   height/width/position lists ARE the parameters) (D12 final, Point 1;
                                                                   replaces inertial_specimen_moi).
calibration_inertia_specimen_from_geometry                bool   Analytic specimen MOI available from the station geometry (D12;
                                                                   replaces inertial_specimen_from_geometry).
calibration_inertia_total_moi_gram_millimeter_squared     float  Total MOI in g*mm^2 = specimen MOI (i_total_system as of
                                                                   v2.7 amendment 2026-07-01). The hardcoded theoretical apparatus
                                                                   baseline (_hardware_inertia_baseline / CLAMP_DIM) is superseded and
                                                                   deleted; the per-configuration apparatus term will be supplied by
                                                                   the calibration-matrix workflow (Deferred flag 4) and added here
                                                                   when available. Until then this key holds specimen MOI only.
                                                                   (D12 final, Point 3; replaces inertial_total_moi.)
calibration_inertia_matrix                                [...]  Apparatus-MOI lookup table keyed by (clamp span,
                                                                   distance-from-AOR); empty array until apparatus-calibration
                                                                   workflow is implemented (D12, Deferred flag 4).

calibration_inertia_apparatus_aor_to_clamp_millimeter    float  Axial distance from the sensor face (the axis of rotation in the
                                                                   standard inline configuration) to the near edge of the rotating
                                                                   clamp body (mm). Source: specimen_profile_clamp_offset_mm (live
                                                                   GUI path, primary) or frustum_inputs['clamp_offset_mm'] (frustum
                                                                   path, fallback); many-to-one, first positive value wins. Key is
                                                                   OMITTED entirely if no positive source is present — never
                                                                   default-written as 0. Feeds apparatus-calibration-matrix lookup
                                                                   (Deferred flag 4). Replaces retired
                                                                   measurement_profile_clamp_offset_millimeter and draft
                                                                   measurement_specimen_frustum_clamp_offset_millimeter. (v2.7)
calibration_inertia_apparatus_plate_to_plate_millimeter  float  Plate-to-plate span between the two body plates within the moving
                                                                   clamp (mm). GUI field "Inter-clamp span". Source:
                                                                   clamp_plate_extension_mm. Feeds apparatus-calibration-matrix
                                                                   lookup (Deferred flag 4). Replaces retired
                                                                   measurement_clamp_plate_extension_millimeter. (v2.7)
calibration_inertia_apparatus_fit_file  str  Provenance path/name of the loaded apparatus-inertia fit artifact
                                                                   (apparatus_inertia_fit.json, built by fit_apparatus_inertia.py from
                                                                   the empty-apparatus calibration runs). 'NA' when no fit is loaded.
                                                                   Direct flat-key write in the exporter. (v2.8)
calibration_inertia_apparatus_fit_form  str  Selected fit-form equation for I_apparatus(aor, width) -- e.g. F4,
                                                                   the separable quadratic a + b*aor^2 + c*width^2. Human-readable
                                                                   mirror of the form stored inside _fit_json. 'NA' when unloaded. (v2.8)
calibration_inertia_apparatus_fit_json  str  The FULL apparatus-inertia fit artifact serialized as a JSON string
                                                                   DATASET (not an attr): fit_form_id + coefficients, source_files,
                                                                   excluded_files, LOO-RMSE, R2, valid_domain, geometry_check, sign
                                                                   note, aor_provenance, and per-trial (aor, width, I, aor_source)
                                                                   rows. Embeds VALUES (not a filename pointer) so the file is
                                                                   self-contained for re-processing, mirroring how the F/T calibration
                                                                   matrix embeds its 6x6 values. 'NA' when no fit is loaded. (v2.8)
```

> **Units convention (Point 1, locked):** MOI keys carry the spelled-out unit suffix
> `_gram_millimeter_squared` in the key name (matching Governing Rule 9). The exporter
> does NOT write a separate `attrs['units']` for these; the suffix is the single source
> of unit truth. `calibration_inertia_bias_newton_meter` and the `measurement_specimen_`*
> geometry/mass/density keys likewise carry their unit in the suffix, not in an attr.

> **D3 (locked).** `calibration_link` subgroup (`use_inertial_calibration`, `calibration_file`, `calibration_available`) is flattened in M1. `calibration_link/calibration_file` → `calibration_inertia_file` (D3, D12). `calibration_link/use_inertial_calibration` and `calibration_link/calibration_available` are **dropped** — no correction runs at export time, so availability and used-flag are meaningless (D3, D12).

> **D12 (locked, finalized 2026-06-22 with Points 1-3).** All computed MOI and inertia-correction parameters live under `calibration_inertia_`*; the `inertial_` prefix section is retired. Finalized decisions:
>
> - **Apparatus MOI = empirical I_est.** `inertial_calibration_profile/I_est` flattens to `calibration_inertia_apparatus_moi_gram_millimeter_squared`, unit-converted `*(180/pi)*1e9` at write time. `bias_est` -> `calibration_inertia_bias_newton_meter`; `axis_sensor` -> `calibration_inertia_axis_sensor`. The `inertial_calibration_profile` subgroup is fully flattened (Point 2).
> - `**calibration_inertia_system_moi` is DROPPED** (Point 3): it had no source attribute (phantom) and "system" was ambiguous; total already carries apparatus+specimen. "Apparatus" now means the rig alone (= I_est).
> - **No provenance blob** (Point 1): `calibration_inertia_moi_provenance` is dropped. Raw geometry/mass/density are first-class named `measurement_specimen`_* fields, so there is nothing left to bury in a JSON blob.
> - **No MOI discretization parameter** (Point 1): the live specimen-MOI path integrates the N stations with an exact per-segment integral; the station lists are the reproducibility inputs. The `num_samples` knobs belong to unused legacy builders and stay excluded.
> - `calibration_inertia_used` and `calibration_inertia_available` are **not written** (no correction at export). Raw geometry/mass/density stay under `measurement_specimen_`*.

### `daq_` — acquisition / DAQ configuration + rig hardware/wiring provenance

All config-sourced hardware/wiring/mechanics provenance is routed canonically under `daq_`
(PI decision: no separate `hardware_` group). The verbatim config-provenance writer is retired;
every config field now flows through the routing ledger except the legacy `units` / `unit_rules`
dicts (dropped as redundant with the spelled-out unit suffixes).

```
daq_ai_channel_map                      channel index → identity map (JSON; e.g. "0":"ai0:xForce" … "6":"ai6:sono_right")
daq_ai_sample_rate_hertz                AI + encoder sample clock (hertz) (Current Python attr: daq_ai_sample_rate_hz, a @property)
daq_ao_do_sample_rate_hertz             AO stim + DO motor sample clock (hertz)
daq_instrumentation                     list of active sensors this session (was: instrumentation)
daq_device_name                         NI device name (Current: device_name)
daq_motor_port                          motor DO port (Current: motor_port)
daq_encoder_channel                     encoder counter channel (Current: encoder_chan)
daq_stim_channels                       stim AO channel list (Current: stim_channels)
daq_stim_channel1 / daq_stim_channel2   per-side stim AO channel (Current: S1stim_chan / S2stim_chan)
daq_stim_channel1_side / _channel2_side stim side labels (Current: S1side / S2side)
daq_stim_monitor_channel / _name        stim-monitor AI channel + name (Current: stim_monitor_chan / stim_monitor_name)
daq_forcetorque_channels / _channel_names   6-axis F/T AI channels + names (Current: SG_chan / SG_name)
daq_sono_channels / _channel_names      sono AI channels + names (Current: sono_channel / sono_name)
daq_sono_enabled                        sono in use this session (Current: use_sono)
daq_sono_internal_sample_rate_hertz     DS3 internal update rate ~241–242 Hz (Current: sono_internal_samplefreq)
daq_motor_gear_ratio                    gearbox ratio, dimensionless (Current: motor_gear_ratio)
daq_motor_full_steps_per_revolution     motor steps/rev (Current: motor_full_steps_per_rev)
daq_encoder_pulses_per_revolution       encoder CPR (Current: encoder_pulses_per_rev)
daq_motor_axis_sensor                   motor rotation axis label (Current: motor_axis)
daq_bending_axis_sensor                 sensor-frame bending axis (Current: bending_axis_sensor)
daq_bending_axis_specimen               specimen-frame bending axis (Current: bending_axis_specimen)
daq_primary_bending_axis                preferred torque axis for QC/correction (Current: primary_bending_axis)
daq_positive_motor_direction            "left" / "right" (Current: positive_motor_direction)
daq_specimen_lateral_index_on_positive_motor_side   signed lateral index (Current: same name)
daq_specimen_side_index_left / _right   derived per-side signed indices
daq_motor_positive_bend_lateral_index   signed lateral index the positive motor direction bends toward (run-computed; isometric/isovelocity)
daq_collection_type                     acquisition type string written by the exporter (D5 — kept alongside protocol_sampling_mode;
                                          both are present; neither is deleted)
```

### `index_` — structural indices into the continuous `timeseries`

Parallel arrays, **one element per step**. Replaces old per-trial subgroups (`trial_index`). Written
by M2 in the same commit that removes `trial_index`. Source of truth: per-step `entry` dicts in the
exporter, enumerated from `bender_h5_export.py:501-531` (D4 — no invented fields).

> **D4 (locked).** Column set = code, not invention. Only fields present in the per-step `entry`
> dicts collected by the acquisition engine appear here. Motor-position reference fields added by
> `_record_motor_position_reference()` are included as a block below; M2 must enumerate their exact
> names from the code before writing.

#### All protocols (shared)

```
index_step_number                       [n_steps] int   1-based step number (matches step_NNN subgroup; source: step_index)
index_step_wall_clock_start             [n_steps] str   ISO-8601 real-world start of this step (no unit suffix)
index_step_duration_second              [n_steps] float Recorded acquisition duration of this step
index_step_rest_before_second           [n_steps] float Unrecorded rest gap before this step; 0 for step 1
index_step_operating_point              [n_steps] float Independent variable value for this step (protocol-dependent)
index_step_operating_point_units        [n_steps] str   Units of operating_point (e.g. "millimeter", "bodylength_per_second")
index_step_recruitment                  [n_steps]       Recruitment condition for this step
index_step_unilateral_posture_lateral_index  [n_steps]  Signed lateral index for unilateral posture
index_step_motor_positive_bend_lateral_index [n_steps]  Signed lateral index the positive motor direction bends toward
index_step_bilateral_mirror_motor       [n_steps] bool  Whether bilateral mirror-motor mode was active
index_step_stim_t0_second              [n_steps] float Stim onset time within the step (source: stim_t0)
index_step_stim_t1_second              [n_steps] float Stim offset time within the step (source: stim_t1)
index_step_t_pre_baseline_start_second  [n_steps] float Pre-stim baseline window start (source: t_pre_baseline_start)
index_step_t_pre_baseline_end_second    [n_steps] float Pre-stim baseline window end   (source: t_pre_baseline_end)
index_step_t_active_start_second        [n_steps] float Active window start             (source: t_active_start)
index_step_t_active_end_second          [n_steps] float Active window end               (source: t_active_end)
index_step_t_post_baseline_start_second [n_steps] float Post-stim baseline window start (source: t_post_baseline_start)
index_step_t_post_baseline_end_second   [n_steps] float Post-stim baseline window end   (source: t_post_baseline_end)
```

#### All protocols (structural — sample offsets into flat/concatenated timeseries)

```
index_step_sample_start                 [n_steps] int   First sample index of this step in the concatenated timeline
index_step_sample_end                   [n_steps] int   Last sample index (inclusive) of this step
index_step_type                         [n_steps] str   Step type label
index_step_target_value_native          [n_steps] float Target value in protocol-native units
index_step_curvature_per_meter          [n_steps] float Step curvature (source: curvature_1_per_m → renamed)
```

#### Isometric-only

```
index_step_target_angle_degree          [n_steps] float Target hold angle (source: target_deg)
index_step_ramp_from_angle_degree       [n_steps] float Starting angle for ramp (source: ramp_from_deg)
```

#### Isovelocity-only

```
index_step_velocity_degree_per_second   [n_steps] float Commanded velocity (source: velocity_deg_s)
index_step_theta_start_degree           [n_steps] float Start angle for isovelocity ramp (source: theta_start_deg)
index_step_pre_stim_t0_second           [n_steps] float Pre-stim window start (source: pre_stim_t0)
index_step_pre_stim_t1_second           [n_steps] float Pre-stim window end   (source: pre_stim_t1)
index_step_iso_t0_second                [n_steps] float Isovelocity ramp onset  (source: iso_t0)
index_step_iso_t1_second                [n_steps] float Isovelocity ramp offset (source: iso_t1)
index_step_mean_xforce_stim_newton      [n_steps] float Mean x-axis force during stim window (source: mean_xforce_stim)
index_step_guard_triggered              [n_steps] bool  Whether the isovelocity angle guard fired (source: guard_triggered)
index_step_guard_angle_degree           [n_steps] float Angle at which guard fired; NaN if not triggered (source: guard_angle_deg)
index_step_velocity_seg1_degree_per_second [n_steps] float Velocity of segment 1 (source: velocity_seg1_deg_s)
index_step_velocity_seg2_degree_per_second [n_steps] float Velocity of segment 2, mirror only (source: velocity_seg2_deg_s)
```

#### Block metadata (present when `use_block_stim` is active)

```
index_step_block_number                 [n_steps] int   Block number (source: block_index; avoids index_step_block_index doubling)
index_step_block_direction              [n_steps] str   Bending direction for the block (source: block_direction)
index_step_block_stim_sides             [n_steps] str   Stim sides for the block (source: block_stim_sides)
index_step_stim_voltage_left_volt       [n_steps] float Left-side stim voltage (source: left_stim_voltage)
index_step_stim_voltage_right_volt      [n_steps] float Right-side stim voltage (source: right_stim_voltage)
```

#### Motor-position reference fields

```
index_step_cumulative_commanded_motor_microstep [n_steps] float  Cumulative Teknic motor-shaft microsteps commanded
                                                          (round(self._motor_continuous_step_pos)); NOT re-zeroed
                                                          between steps. Units: motor-shaft microsteps (1600/rev
                                                          on the ClearPath shaft; gear ratio 5:1 to specimen frame
                                                          = 8000 microsteps per output-shaft revolution). NOT
                                                          output-shaft steps or degrees.
                                                          (source: cumulative_commanded_steps; renamed D9)
index_step_encoder_cumulative_degree    [n_steps] float  Running sum of per-step net encoder displacement
                                                          (angle_measured[-1] - angle_measured[0] per step);
                                                          accumulated without re-zeroing across step boundaries.
                                                          Encoder is US Digital E6 on the gearbox output
                                                          (specimen frame) — values are in specimen-frame degrees.
                                                          (source: encoder_cumulative_deg)
```

> **Dropped (D4).** `trial_index`, `cycle_index`, and `step_order` per-step scalars are not written
> as `index_step_`* arrays. `trial_index` and `cycle_index` are redundant with `step_manifest` +
> per-sample timeseries streams. `step_order` is reconstructable from per-sample `sequence_index` +
> `block_index`.

### `index_cycle_`* — per-cycle design grid (`single_finite` only; D10)

Parallel arrays, **one element per bending cycle** (length C = total cycles in trial). `**single_finite` protocols only** — `segmented_finite` protocols have no per-cycle design grid. Join to per-sample data via the `cycle_index` stream already in `timeseries` (§3b-i). Key names locked at M0 (D10); **implemented in M2c** (exporter writes these from the `organize_cycles`-built per-cycle arrays; `index_cycle_operating_point`/`_units` follow the PI-approved native-unit mapping — `angle`→degree, `strain`→strain fraction, `strain_pct`→percent, `curvature`→per_meter; `index_cycle_motor_amplitude_degree` is always written in degrees).

> **D10 (PI sign-off 2026-06-22).** `Lonoff`/`Ronoff` are omitted (redundant with per-sample `stim_side`/`stim_state`). Does not add `index_step_`* columns (D4 intact). `protocol_activation_start_cycle` is 0-based (see FLAG below); all other `index_cycle_index` values are 1-based.

```
index_cycle_index                       [C] int     1-based cycle number
index_cycle_frequency_hertz             [C] float   Commanded cycle frequency for this cycle
index_cycle_operating_point             [C] float   Drive amplitude in protocol-native unit
                                                      (matches index_step_operating_point pattern)
index_cycle_operating_point_units       [C] str     Unit string for operating_point (e.g. "degree", "millimeter")
index_cycle_motor_amplitude_degree      [C] float   Motor commanded amplitude (degree)
index_cycle_active                      [C] bool    Whether stimulation was commanded this cycle
index_cycle_activation_duty             [C] float   Stim duty cycle, dimensionless (0–1)
index_cycle_activation_phase            [C] float   Stim activation phase, dimensionless (0–1)
```

### `inertial_` — **RETIRED (D12)**

> **D12 (PI sign-off 2026-06-22; finalized with Points 1-3).** The `inertial_` prefix section is retired. All computed MOI and inertia-correction parameters are now written under `calibration_inertia_*` (see `calibration_` section above). Field mapping:
>
>
> | Old source                                 | New `calibration_inertia_*` key                                               |
> | ------------------------------------------ | ----------------------------------------------------------------------------- |
> | `inertial_specimen_moi`                    | `calibration_inertia_specimen_moi_gram_millimeter_squared`                    |
> | `inertial_total_moi`                       | `calibration_inertia_total_moi_gram_millimeter_squared`                       |
> | `inertial_calibration_profile/I_est`       | `calibration_inertia_apparatus_moi_gram_millimeter_squared` (`*(180/pi)*1e9`) |
> | `inertial_calibration_profile/bias_est`    | `calibration_inertia_bias_newton_meter`                                       |
> | `inertial_calibration_profile/axis_sensor` | `calibration_inertia_axis_sensor`                                             |
> | `inertial_specimen_from_geometry`          | `calibration_inertia_specimen_from_geometry`                                  |
> | `inertial_system_moi`                      | **DROPPED** (phantom — no source attr; "apparatus" now = I_est)               |
> | `inertial_moi_provenance`                  | **DROPPED** (geometry/mass/density are first-class `measurement_specimen_*`)  |
> | `inertial_system_from_profile`             | **DROPPED**                                                                   |
>
>
> "inertial" → "inertia" normalization and the full flattening of `inertial_calibration_profile` are applied by M1; the `inertial_` prefix does not appear in any new code after M1.

### `measurement_` — ALL spatial geometry (animal AND apparatus)

**Any distance or length → `measurement_`**, regardless of whether it describes the specimen or the fixture.

```
measurement_specimen_frustum_heights_millimeter     [3]   (Current: specimen_geometry_heights_mm)
measurement_specimen_frustum_depths_millimeter      [3]   (Current: specimen_geometry_depths_mm; the "width" stations)
measurement_specimen_frustum_positions_millimeter   [3]   (Current: specimen_geometry_positions_mm)
measurement_specimen_density_gram_per_cubic_millimeter     (Current: specimen_geometry_density_g_per_mm3; MOVED here from the
                                                            retired inertial_specimen_density_* — it is a measured input, not a
                                                            computed output, Point 1)

measurement_specimen_inertia_frustum_height_millimeter      float  Back-end (maximum) cross-section height of the elliptical
                                                                frustum MOI model, in mm. (Current: frustum_height_mm;
                                                                v2.7, renamed with "inertia" segment 2026-07-01)
measurement_specimen_inertia_frustum_width_millimeter       float  Back-end (maximum) cross-section width of the elliptical
                                                                frustum MOI model, in mm. (Current: frustum_width_mm;
                                                                v2.7, renamed 2026-07-01)
measurement_specimen_inertia_frustum_length_millimeter      float  Axial span of the frustum (head-to-tail span modeled), in mm.
                                                                (Current: frustum_length_mm; v2.7, renamed 2026-07-01)
measurement_specimen_inertia_frustum_density_gram_per_cubic_millimeter  float  Material density for the frustum MOI model (g/mm^3).
                                                                Distinct from measurement_specimen_density_gram_per_cubic_millimeter
                                                                (which comes from specimen_geometry_density_g_per_mm3 / the station
                                                                model). (Current: frustum_density_g_per_mm3; v2.7, renamed 2026-07-01)
measurement_specimen_bodylength_millimeter                 (Current: fishlen_TL — total length)
measurement_specimen_standardlength_millimeter            (Current: fishlen_SL — standard length)
measurement_specimen_body_mass_gram                       (Current: fishmass)
measurement_clamp_separation_millimeter                   (Current: dclamp)
measurement_clamp_offset_vertical_millimeter              (Current: dvert)
measurement_clamp_offset_horizontal_millimeter            (Current: dhoriz)
measurement_target_muscle_depth_millimeter                (Current: target_muscle_depth_mm)
measurement_specimen_local_body_width_millimeter          Local body width at the test-section site (mm). Distinct from whole-body summary morphometrics. (Current: xsec_width; renamed 2026-07-02, retired measurement_xsec_width_millimeter)
measurement_specimen_local_body_height_millimeter         Local body height at the test-section site (mm). Distinct from whole-body summary morphometrics. (Current: xsec_height; renamed 2026-07-02, retired measurement_xsec_height_millimeter)
```

> **Body length — LOCKED:** two keys, compound words run together per naming convention: `measurement_specimen_bodylength_millimeter` (TL) and `measurement_specimen_standardlength_millimeter` (SL).

### `note_` — bench and post-hoc notes

Raw file holds bench-side notes only. Layer-3 exclusion verdicts live in the hub, not here.

```
note_bench                              apparatus/bench-side observations (replaces note_stim_clamp / stim_clamp_notices)
note_posthoc                            post-hoc analysis notes (new; populated during analysis, not at acquisition)
```

*(Previous `note_stim_clamp` is retired; its content moves to `note_bench`.)*

### `protocol_` — global protocol settings (defaults)

Global defaults. Per-step *realized* values live in `index_step_`*.

```
protocol_type                           file-level protocol selector (was: test_type)
protocol_starting_angle_degree          (was: starting_angle_deg)
protocol_duration_second                (was: duration)
protocol_isometric_initial / _final / _num_steps / _target_unit
protocol_isovelocity_min_velocity / _max_velocity / _starting_strain / _*_unit
protocol_step_frequency_hertz / _step_amplitude / _step_amplitude_unit / _step_curvature_per_meter
protocol_amplitude_frequency_exponent / _velocity_exponent
protocol_sampling_mode                  'single_finite' | 'segmented_finite'
protocol_acquisition_mode               'continuous' for the stitched isometric engine (run-computed)
protocol_guard_triggered                bool — isovelocity angle guard fired (run-computed)
protocol_guard_angle_degree             angle (deg) the isovelocity guard fired at; NaN if not (run-computed)
... (per protocol type — see bender_routing_spec.py for the full list)

# config-sourced timing defaults (routed under protocol_)
protocol_wait_before_second             (was: waitbefore)
protocol_wait_after_second              (was: waitafter)
protocol_ramp_duration_second           (was: rampdur)
protocol_prepoststim_duration_second    (was: prepoststim_dur)
protocol_prepoststim_separation_second  (was: prepoststim_sep)
protocol_prestim_time_second            (was: prestim_time)
protocol_poststim_time_second           (was: poststim_time)
protocol_ramp_mode                      linear/exponential (was: ramp_mode_default)
protocol_amplitude_step_velocity_degree_per_second  (was: amp_step_vel)

# Single_finite per-cycle design scalars (D10 — single_finite only; locked at M0)
protocol_activation_start_cycle         int     0-based index of first activated cycle
                                                  (FLAG: 0-based unlike index_cycle_index which is 1-based;
                                                  PI-confirmed asymmetry, do not normalize)
protocol_cycles_per_step                int     Nominal bending cycles per step (single_finite; 1 for non-sweep)
protocol_cycle_count                    int     Total bending cycles in trial
protocol_end_cycle_count                int     Number of post-stim passive cycles at end of trial
```

### `session_` — recording-session logistics

```
session_date                            ISO-8601 timestamp of run start (YYYY-MM-DDTHH:MM:SS)
session_apparatus_id                    rig identity label (was: apparatus_id; current Python attr: apparatus_id)
session_analyst                         operator initials/name (was: analyst)
session_simulated                       boolean — True if acquisition used synthetic DAQ data; False if real hardware
                                          (Current: 01_Metadata.attrs['simulated'] — one-line rename approved, pending code change)
session_protocol_template               template config file name
```

### `specimen_` — non-spatial animal facts

Identity and biology only — **no geometry** (geometry → `measurement_`). All six fields below are now confirmed logged.

```
specimen_id                             (conformant; currently in protocol_metadata/attrs)
specimen_genusspecies                   binomial name, e.g. "Micropterus salmoides" (Current: genus_species — rename pending)
specimen_sex                            biological sex, e.g. "M" / "F" / "unknown" (now written; Current: specimen_sex in bender_settings)
specimen_muscle_type                    preparation region, e.g. "epaxial" (now written; Current: specimen_muscle_type)
specimen_prep_condition                 prep state, e.g. "in_vivo" / "in_vitro" (Current: prep_condition)
```

> **FLAG E — DONE (D9).** `genus_species` → `specimen_genusspecies` implemented: Python attr `bender.specimen_genusspecies` in `bender_functions.py`; `BENDER_ROUTING` source key `specimen_genusspecies`; GUI session key `gui_genus_species` intentionally kept for JSON template compatibility (see §8).

### `step_manifest` — per-step timing and operating-point index

JSON dataset in `metadata/`. Present for **all** protocols.


| Field                | Type   | Notes                                                                                                                   |
| -------------------- | ------ | ----------------------------------------------------------------------------------------------------------------------- |
| `step_index`         | int    | 1-based; matches `step_NNN` subgroup index in `segmented_finite` (`step_001` → `step_index = 1`); 1 for `single_finite` |
| `wall_clock_start`   | string | ISO-8601 real-world start of this step                                                                                  |
| `duration_second`    | float  | Recorded duration of this step                                                                                          |
| `rest_before_second` | float  | Unrecorded rest gap before this step; 0 for `step_index = 1` (first step)                                               |


Additional fields for `**segmented_finite`** only:


| Field                   | Type   | Notes                                                                          |
| ----------------------- | ------ | ------------------------------------------------------------------------------ |
| `operating_point`       | float  | Independent variable for this step                                             |
| `operating_point_units` | string | e.g. `millimeter` for FL/isometric, `bodylength_per_second` for FV/isovelocity |


`single_finite` manifests have exactly 1 row (`step_index = 1`, `rest_before_second = 0`, no `operating_point` fields).

### `99_Unrouted` — catch-all subgroup for unclassified fields

**Location:** `metadata/99_Unrouted` (subgroup; D7 confirmed).

Any config-sourced field that the routing ledger does not yet classify into a prefix group is written
here rather than at the `metadata` root or the file root. This prevents unclassified fields from
polluting the flat namespace while the routing ledger is being completed.

- Present only when there are unrouted fields; absent from a fully-routed file.
- **After M1** this is the only subgroup permitted under `metadata`; all other subgroups
(`calibration_link`, `inertial_calibration_profile`, `protocol_metadata`, `bender_settings`,
`trial_index`) are removed or flattened as part of M1/M2.
- Keys inside `99_Unrouted` follow the same `<description>_<unit>` convention; they are candidates
for routing in a follow-up pass and should not be read as canonical.
- Root holds only `schema_version` and the two groups `metadata` / `timeseries`; `99_Unrouted`
moves from root to `metadata/99_Unrouted` in M1 (D7).

---

## 5. Raw + calibration-values principle (D11 — supersedes prior §5)

> **D11 (PI sign-off 2026-06-22) SUPERSEDES the prior "raw + decoded, both in file" rule.** The new architecture is: raw voltages + embedded calibration VALUES in the raw file; calibrated output is NOT ground truth in the raw file.

> **D11 amendment (Point 4, finalized 2026-06-22).** The Python exporter ALSO writes a small `derived/` group into the raw file for **live RA inspection during a run** — NOT as ground truth. Contents:
>
> - `derived/forcetorque_calibrated` — ATI 6×6 matrix applied to the raw F/T rows (reuses the already-computed `self.forcetorque` / `apply_calibration_forcetorque`). Always written when the calibration matrix is present.
> - `derived/torque_inertia_corrected` — calibrated primary torque minus `I_sum * alpha`, where `I_sum` = the SUM of whichever inertias are present: apparatus (`I_est`) and/or specimen MOI. Written whenever at least one is present; silently skipped (no error) when neither is.
> - **Unit guard (mandatory):** keep torque in N*m. Use `I_est` in its native `N*m/(deg/s^2)` with `alpha` in deg/s^2, OR convert fully to SI (`I[kg*m^2] * alpha[rad/s^2]`). The stored g*mm^2 values are for cross-comparability only and must NOT be multiplied by deg/s^2. Do not reuse the dead, unit-inconsistent `get_corrected_torque`.
> - The R pipeline ignores `derived/` entirely and re-derives everything from raw + embedded calibration VALUES. This group exists so the RA can distinguish a clean signal from a dirty one and see the muscle signal with apparatus inertia removed, live.

1. **Raw is immutable.** `aidata` (raw ADC voltages, all channels) is written as-acquired and never overwritten, modified, or deleted. It is the permanent audit trail.
2. **Calibration VALUES are embedded at collection.** `calibration_forcetorque_matrix` (6×6 ATI matrix), `calibration_sono_left`/`right` (breakpoint tables), and inertia parameters (`calibration_inertia_`*) are written to `metadata/` at write time so the file is self-contained for re-processing without the original config.
3. **No calibrated decoded array in timeseries.** The Python exporter does **not** write a calibrated F/T array to `timeseries`; `forcetorque_raw` is never written (M2b). The inspection-only `derived/forcetorque_calibrated` is written ONLY when a REAL `calibration_forcetorque_matrix` is embedded (the `np.eye(6)` identity fallback is refused, so an absent matrix means "not calibrated" and no derived output).
4. **Downstream re-derives from raw + embedded VALUES.** R reads `aidata` rows 0–5 (identified via `daq_ai_channel_map`) + `calibration_forcetorque_matrix` to produce calibrated F/T. All correction and calibration happens in R and lands in hub `derived/` only.
5. **Rationale:** clean separation of immutable raw archive from analysis outputs; calibration values travel with the file; no decoded F/T to go stale if recalibration is needed.

In `timeseries`, this means:

- `aidata` (raw ADC voltages, `[8 × N]`) — immutable, always present; F/T = rows 0–5 via `daq_ai_channel_map`
- `forcetorque_raw` — **never written** (M2b); preserved only in legacy pre-M2b files where it already exists
- Calibrated F/T, inertial-corrected torque — authoritative copies live in hub `derived/` (R-computed); the raw file additionally carries an exporter-written `derived/` group for live RA inspection only (D11 amendment, Point 4), written ONLY when a REAL calibration matrix is embedded, which R ignores

---

## 6. Current-file → target mapping (audited keys)

From the audit of `2026-06-04_bass13_bender_24_isometric.h5` and the sim-validation session of `2026-06-12`. Resolves the known `metadata` payload:


| Current key (location)                                             | Status                                           | Target key                                                               | Action                                                                                                                                     |
| ------------------------------------------------------------------ | ------------------------------------------------ | ------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------ |
| `calibration` (bender_settings)                                    | present                                          | `calibration_forcetorque_matrix`                                         | rename                                                                                                                                     |
| `sono_cal_left`                                                    | present                                          | `calibration_sono_left_volt_millimeter_breakpoints`                      | rename (M2b — breakpoint-table suffix, was draft `calibration_sono_left`)                                                                  |
| `sono_cal_right`                                                   | present                                          | `calibration_sono_right_volt_millimeter_breakpoints`                     | rename (M2b — breakpoint-table suffix, was draft `calibration_sono_right`)                                                                 |
| `genus_species` (bender_settings)                                  | present                                          | `specimen_genusspecies`                                                  | **done (D9)** — Python attr + routing spec use `specimen_genusspecies`; GUI key `gui_genus_species` kept for JSON compat                   |
| `prep_condition` (bender_settings)                                 | present                                          | `specimen_prep_condition`                                                | rename                                                                                                                                     |
| `specimen_sex` (bender_settings)                                   | present                                          | `specimen_sex`                                                           | conformant                                                                                                                                 |
| `specimen_muscle_type` (bender_settings)                           | present                                          | `specimen_muscle_type`                                                   | conformant                                                                                                                                 |
| `session_analyst` (bender_settings)                                | present                                          | `session_analyst`                                                        | conformant                                                                                                                                 |
| `session_date` (01_Metadata top attrs)                             | present                                          | `session_date`                                                           | conformant; schema rename to `session_date` pending (already correct)                                                                      |
| `apparatus_id` (01_Metadata top attrs)                             | present                                          | `session_apparatus_id`                                                   | rename                                                                                                                                     |
| `simulated` (01_Metadata top attrs)                                | present                                          | `session_simulated`                                                      | **done (D9)** — Python attr `bender.session_simulated` + routing spec already use `session_simulated`                                      |
| `specimen_geometry_heights_mm`                                     | present                                          | `measurement_specimen_frustum_heights_millimeter`                        | rename                                                                                                                                     |
| `specimen_geometry_depths_mm`                                      | present                                          | `measurement_specimen_frustum_depths_millimeter`                         | rename                                                                                                                                     |
| `specimen_geometry_positions_mm`                                   | present                                          | `measurement_specimen_frustum_positions_millimeter`                      | rename                                                                                                                                     |
| `fishlen_TL`                                                       | present                                          | `measurement_specimen_bodylength_millimeter`                             | **done (D9)** — routing spec maps `fishlen_TL` to canonical name; Python attr kept as `fishlen_TL` for JSON template compat                |
| `fishlen_SL`                                                       | present                                          | `measurement_specimen_standardlength_millimeter`                         | **done (D9)** — routing spec maps `fishlen_SL` to canonical name; Python attr kept as `fishlen_SL` for JSON template compat                |
| `fishmass`                                                         | present                                          | `measurement_specimen_body_mass_gram`                                    | rename                                                                                                                                     |
| `dclamp`                                                           | present                                          | `measurement_clamp_separation_millimeter`                                | rename                                                                                                                                     |
| `dvert`                                                            | present                                          | `measurement_clamp_offset_vertical_millimeter`                           | rename                                                                                                                                     |
| `dhoriz`                                                           | present                                          | `measurement_clamp_offset_horizontal_millimeter`                         | rename                                                                                                                                     |
| `clamp_plate_extension_mm`                                         | present                                          | `calibration_inertia_apparatus_plate_to_plate_millimeter`                | reclassified to calibration_inertia_apparatus_* (v2.7)                                                                                     |
| `xsec_width`                                                       | present                                          | `measurement_xsec_width_millimeter`                                      | rename                                                                                                                                     |
| `xsec_height`                                                      | present                                          | `measurement_xsec_height_millimeter`                                     | rename                                                                                                                                     |
| `target_muscle_depth_mm`                                           | present                                          | `measurement_target_muscle_depth_millimeter`                             | rename                                                                                                                                     |
| `inertial_torque_specimen_primary`                                 | empty `(0,)`                                     | `calibration_inertia_specimen_moi_gram_millimeter_squared`               | rename + reclassify (D12)                                                                                                                  |
| `inertial_torque_system_primary`                                   | empty `(0,)`                                     | —                                                                        | **DROP** (phantom; "apparatus" now = I_est, Point 3, D12)                                                                                  |
| `inertial_torque_total_primary`                                    | empty `(0,)`                                     | `calibration_inertia_total_moi_gram_millimeter_squared`                  | rename + reclassify (D12)                                                                                                                  |
| `inertial_calibration_profile/I_est` (subgroup)                    | present                                          | `calibration_inertia_apparatus_moi_gram_millimeter_squared`              | flatten + unit-convert `*(180/pi)*1e9` (Point 2, D12)                                                                                      |
| `inertial_calibration_profile/bias_est` (subgroup)                 | present                                          | `calibration_inertia_bias_newton_meter`                                  | flatten (Point 2, D12)                                                                                                                     |
| `inertial_calibration_profile/axis_sensor` (subgroup)              | present                                          | `calibration_inertia_axis_sensor`                                        | flatten (Point 2, D12)                                                                                                                     |
| `specimen_geometry_density_g_per_mm3`                              | present                                          | `measurement_specimen_density_gram_per_cubic_millimeter`                 | move to measurement_ (measured input, Point 1)                                                                                             |
| `frustum_height_mm` (attr on Bender / frustum_inputs dict)         | present when frustum path used                   | `measurement_specimen_inertia_frustum_height_millimeter`                 | first-class key; read from frustum_inputs dict by exporter (v2.7, renamed 2026-07-01)                                                      |
| `frustum_width_mm` (attr on Bender / frustum_inputs dict)          | present when frustum path used                   | `measurement_specimen_inertia_frustum_width_millimeter`                  | first-class key; read from frustum_inputs dict by exporter (v2.7, renamed 2026-07-01)                                                      |
| `frustum_length_mm` (attr on Bender / frustum_inputs dict)         | present when frustum path used                   | `measurement_specimen_inertia_frustum_length_millimeter`                 | first-class key; read from frustum_inputs dict by exporter (v2.7, renamed 2026-07-01)                                                      |
| `frustum_clamp_offset_mm` (frustum_inputs dict)                    | present when frustum path used                   | `calibration_inertia_apparatus_aor_to_clamp_millimeter`                  | apparatus geometry (not specimen); many-to-one with specimen_profile_clamp_offset_mm; omit if not positive (v2.7 amended 2026-07-01)       |
| `specimen_profile_clamp_offset_mm`                                 | present when live geometry path used             | `calibration_inertia_apparatus_aor_to_clamp_millimeter`                  | primary source for AOR-to-clamp key; many-to-one (v2.7 amended 2026-07-01); replaces measurement_profile_clamp_offset_millimeter (retired) |
| `frustum_density_g_per_mm3` (attr on Bender / frustum_inputs dict) | present when frustum path used                   | `measurement_specimen_inertia_frustum_density_gram_per_cubic_millimeter` | first-class key; read from frustum_inputs dict by exporter (v2.7, renamed 2026-07-01)                                                      |
| `frustum_tip_scale` (attr on Bender)                               | present when frustum path used                   | —                                                                        | **DROP** (redundant when N-station geometry captures each cross-section; v2.7)                                                             |
| `primary_torque_raw`                                               | empty `(0,)`                                     | —                                                                        | delete from metadata (derived series → hub)                                                                                                |
| `primary_torque_corrected`                                         | empty `(0,)`                                     | —                                                                        | delete from metadata (derived series → hub)                                                                                                |
| `stim_clamp_notices`                                               | empty `(0,)`                                     | `note_bench`                                                             | rename                                                                                                                                     |
| `trial_index/`* (subgroup)                                         | present (minus `trial_index`/`cycle_index` cols) | `index_step_`* parallel arrays                                           | move + flatten (later)                                                                                                                     |
| `protocol_metadata/step_order`                                     | **dropped**                                      | —                                                                        | done (reconstructable from per-sample `sequence_index` + `block_index`)                                                                    |
| `protocol_metadata/block_sequence`                                 | `**protocol_block_sequence` (JSON)**             | —                                                                        | done (routed, `json.dumps`)                                                                                                                |
| `protocol_metadata` (subgroup)                                     | **removed**                                      | —                                                                        | done (attr-mirrors now canonical; redundant dict-only keys dropped)                                                                        |
| `calibration_link/use_inertial_calibration`                        | present                                          | —                                                                        | **DROP** (no correction at export; D3, D12)                                                                                                |
| `calibration_link/calibration_file`                                | present                                          | `calibration_inertia_file`                                               | flatten from subgroup (M1); "inertial" -> "inertia" (D3, D12)                                                                              |
| `calibration_link/calibration_available`                           | present                                          | —                                                                        | **DROP** (no correction at export; D3, D12)                                                                                                |
| `daq_ai_sample_rate_hz`                                            | present                                          | `daq_ai_sample_rate_hertz`                                               | rename                                                                                                                                     |
| `t` (timeseries)                                                   | present                                          | `time_second`                                                            | rename                                                                                                                                     |
| `tnorm` (timeseries)                                               | present                                          | `time_normalized`                                                        | rename                                                                                                                                     |
| `angle_cmd` (timeseries)                                           | present                                          | `angle_commanded_degree`                                                 | rename                                                                                                                                     |
| `angle_measured` (timeseries)                                      | present                                          | `angle_measured_degree`                                                  | rename                                                                                                                                     |
| `anglevel_cmd` (timeseries)                                        | present                                          | `angular_velocity_commanded_degree_per_second`                           | rename                                                                                                                                     |
| `S1stimcmd` (timeseries)                                           | present                                          | `stim_channel1_command_volt`                                             | rename                                                                                                                                     |
| `S2stimcmd` (timeseries)                                           | present                                          | `stim_channel2_command_volt`                                             | rename                                                                                                                                     |
| `forcetorque` (timeseries)                                         | **dropped**                                      | —                                                                        | done (identical copy of `forcetorque_raw`; not split, PI decision)                                                                         |
| `forcetorque_raw` (timeseries)                                     | present                                          | —                                                                        | **NEVER WRITTEN** (D11, M2b — unconditional; calibrated F/T → `derived/forcetorque_calibrated`, gated on a real embedded matrix); preserved only in legacy pre-M2b files |


---

## 7. Governing rules

1. **Immutability.** `aidata` (and all raw streams) are written as-acquired and never modified post-hoc. All correction happens in R and lands in `derived` in the hub copy. R is the authoritative correction path.
2. **Raw + calibration VALUES in file; decoded output in hub only (D11 — supersedes prior Rule 2).** The raw file stores raw ADC voltages (`aidata`) and embedded calibration VALUES (`calibration_forcetorque_matrix`, `calibration_sono_`*, `calibration_inertia_`*). Calibrated F/T, corrected torque, and all derived quantities are computed in R and written only to hub `derived/`. The exporter does not compute or write a decoded F/T array. See §5 for the full rule.
3. **Decode via metadata, not per-channel.** Channel identity lives in `daq_ai_channel_map`. The exporter does not hard-slice `aidata` into named per-channel datasets; R reads `daq_ai_channel_map` to identify F/T rows and apply the calibration matrix (D11; see §3d).
4. **Distance rule.** Any length/distance → `measurement_`, animal or apparatus. `specimen_` holds only non-spatial facts (id, species, sex, muscle type — no geometry, no mass).
5. **Step parameters appear twice by design.** Global defaults in `protocol_`*; per-step realized values in `index_step_`*. Intentional, not duplication.
6. **Two timeseries structures, determined by `protocol_sampling_mode`.** `single_finite`: flat channel arrays directly in `timeseries/`, one continuous monotonic time axis. `segmented_finite`: `step_NNN/` subgroups in `timeseries/`, each with a local time axis starting at 0; `step_manifest` in `metadata/` is the only cross-step index. Subgroup prefix is always `step_` — never `trial_` or any other prefix. The locked protocol→mode mapping is in §10.
7. **Two files, shared group names.** Raw file = `metadata` + `timeseries` (+ an inspection-only exporter-written `derived/`, D11 amendment / Point 4). Hub copy = the two groups + the authoritative R-computed `derived/`. Group names identical across both; the raw-file `derived/` is RA-inspection only and is not ground truth.
8. **Don't reuse one variable for two physical quantities.** Curvature, strain, angle, width, frequency are distinct even though all floats. Verify each assignment by physical role.
9. **Naming convention is locked.** `<description>_<unit>`, fully spelled-out units, no abbreviations. See §1.
10. **Every user-facing GUI input field must have a corresponding exported HDF5 key.** No GUI input may be silently dropped from export — silence means a reader cannot distinguish "value not recorded" from "field was zero". Applies to all new fields going forward. A retroactive audit of existing GUI fields is a separate follow-up, not part of this rule's effective date (2026-07-01).

---

## 8. Resolved flags (D9 — all done)

All schema-level decisions are locked. All code-side renames are complete as of D9 (commit following M2):


| Flag                    | Decision                                                                                                  | Status (D9)                                                                                                                                                               |
| ----------------------- | --------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| A — ATI units           | newton / newton_meter (SI)                                                                                | **DONE** — `ts_rename` in exporter already uses `_newton` / `_newton_meter` in channel names                                                                             |
| B — stim channel naming | `stim_channel1_command_volt` / `stim_channel2_command_volt`                                               | **DONE** — `bender_h5_export.py` `ts_rename` maps `S1stimcmd` → `stim_channel1_command_volt`, `S2stimcmd` → `stim_channel2_command_volt`                                  |
| C — body length keys    | `measurement_specimen_bodylength_millimeter` (TL) + `measurement_specimen_standardlength_millimeter` (SL) | **DONE** — `BENDER_ROUTING` maps `fishlen_TL` / `fishlen_SL` to canonical names; Python attrs and GUI session keys intentionally keep old names for JSON template compat  |
| D — simulated key       | `session_simulated`                                                                                       | **DONE** — Python attr `bender.session_simulated` and `BENDER_ROUTING` source key already use `session_simulated`                                                        |
| E — species key         | `specimen_genusspecies`                                                                                   | **DONE** — Python attr `bender.specimen_genusspecies` renamed; `BENDER_ROUTING` source key `specimen_genusspecies`; GUI session key `gui_genus_species` intentionally kept for JSON template compat |


- **Decision 2 — body-width timing.** RESOLVED 2026-07-02: `measurement_specimen_body_width_millimeter` retired; superseded by `measurement_specimen_local_body_width_millimeter` (routed from `xsec_width`, renamed same date).

---

## 9. Backlog (logged, not blocking)

- Update `tier2_data_curation_sop` + `jlab_folder_architecture` to the new key names; update any provenance assertions (`apparatus_id` → `session_apparatus_id`, etc.).

---

## 10. Locked protocol → sampling mode mapping

PI decision. This mapping is immutable; changing a protocol's mode requires a schema version bump.


| Protocol                | `protocol_sampling_mode` | Notes                                                       |
| ----------------------- | ------------------------ | ----------------------------------------------------------- |
| dynamic                 | `single_finite`          | One finite DAQmx task per trial; flat continuous timeseries |
| sweep (frequency_sweep) | `single_finite`          | One finite DAQmx task per trial; flat continuous timeseries |
| isometric               | `segmented_finite`       | One finite DAQmx task per step; `step_NNN/` subgroups       |
| isovelocity             | `segmented_finite`       | One finite DAQmx task per step; `step_NNN/` subgroups       |
| FL (force-length)       | `segmented_finite`       | One finite DAQmx task per step; `step_NNN/` subgroups       |
| FV (force-velocity)     | `segmented_finite`       | One finite DAQmx task per step; `step_NNN/` subgroups       |


> **Current implementation note.** Isometric and isovelocity currently run as a stitched continuous timeline (one NI task per trial, not per step). The `segmented_finite` entry above is the **target architecture** — this schema defines the goal, not the current code state. Do not align future code to `single_finite` for these protocols; align to `segmented_finite` when the acquisition engine is updated.

---

## 11. Fixed-schema null-sentinel rule (2026-07-02)

The schema is **fixed**: every canonical metadata **scalar** is **always emitted**. There is no
omit-on-zero and no skip-on-None. When a value is unrecorded/absent, the writer emits a **typed NA
sentinel** instead of dropping the key. A reader can therefore assume the full canonical key set is
present in every post-amendment file and never has to distinguish "key absent" from "value missing."

### 11.1 Per-dtype sentinel table

| dtype       | present value            | absent (NA) sentinel                    |
| ----------- | ------------------------ | --------------------------------------- |
| float       | the float                | `NaN`                                   |
| int / count | stored as float          | float `NaN`                             |
| string      | the string               | `"NA"`                                  |
| bool        | `"true"` / `"false"`     | `"NA"` (3-state string)                 |
| array       | the array                | empty array (shape `(0,)`)              |

Booleans are stored as **3-state strings** (`"true"`, `"false"`, `"NA"`) — never native HDF5 bool —
so absence is representable in the same field. Readers must parse the string, not coerce with
`bool()` (the Python string `"false"` is truthy).

Fields representing counts or indices are integer-valued floats; a present value will always be
whole-numbered. Readers may cast to int after confirming non-NaN.

The writer infers the sentinel dtype for an absent key from its **name suffix** (D0-B):
numeric suffixes (`_millimeter`, `_second`, `_hertz`, `_degree`, `_volt`, `_count`, `_index`,
`_squared`, `_celsius`, `_newton`, `_gram`, `_num_steps`) -> `NaN`; bool suffixes (`_enabled`,
`_simulated`, `_triggered`, `_bilateral`, `_randomize`, `_from_geometry`) and any no-match ->
`"NA"`. A *present* value always keeps its native type; suffix inference is used only to pick the
sentinel for an absent field.

### 11.2 Guardrail 1 — required fields reject NA

A field marked `required: True` must be **present AND non-NA**. Validation treats a required field
that is missing, `""`, `"NA"`, or `NaN` as a **CRITICAL** failure — not a warning. NA is legitimate
only for optional / unrecorded fields. (The one narrow exception: `specimen_id` may be NA for a
simulated run, where there is no physical specimen.)

### 11.3 Guardrail 2 — canonical schema only; dropped keys stay dropped

"Always present" means every field in the **canonical** schema is emitted. It must **not** resurrect
deliberately-dropped legacy keys. The NA-fill sweep iterates the routing table only; it never
invents keys. The following stay **dropped** and must never be re-added by the always-present rule:

- `calibration_inertia_system_moi` (and `_used`, `_available` variants)
- the `protocol_metadata` subgroup (flattened to `protocol_*` scalars)
- `calibration_forcetorque_matrix_used` / `_available`
- `measurement_xsec_width_millimeter` (retired 2026-07-02; replaced by `measurement_specimen_local_body_width_millimeter`)
- `measurement_xsec_height_millimeter` (retired 2026-07-02; replaced by `measurement_specimen_local_body_height_millimeter`)
- `measurement_specimen_body_width_millimeter` (retired -- superseded by `measurement_specimen_local_body_width_millimeter`)

### 11.4 Boundary — this rule covers metadata SCALARS ONLY

The null-sentinel rule applies to **metadata scalars** (attributes and scalar datasets under
`metadata/`). It does **not** apply to **timeseries channels**. A timeseries channel that was not
recorded is an **absent or empty dataset** — never a NaN-filled array of run length. A NaN-filled
timeseries would falsely imply the channel was acquired; an absent/empty dataset correctly signals
it was not. Matrix datasets (e.g. `calibration_forcetorque_matrix`) likewise remain present-only
when real: an absent calibration matrix means uncalibrated, not a NaN matrix.

