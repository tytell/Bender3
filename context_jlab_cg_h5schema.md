# Bender / CritterGripper — HDF5 Schema

## Authority metadata
- **Purpose:** Defines the on-disk structure, group layout, key-naming taxonomy, and unit conventions for Bender HDF5 files.
- **Authority level:** Tier 2 (binding for all H5 write/read tasks — export routine, R pipeline, GUI→backend assignment).
- **Scope:** `bender_h5_export.py` (write), `01_calibration.R` / `02_correct.R` / `03_analyze.R` (read), GUI variable assignment in `bender_streamlit_gui.py`.
- **Owner:** PI (schema decisions are PI-only).
- **Version:** v2.3 (Phase 0 + index/stim schema: ledger-driven canonical exporter; `forcetorque_raw` kept `[6×N]` not split; index fields `cycle_index`/`step_index`/`sequence_index`/`block_index`/`block_direction`/`block_stim_sides` converted to per-sample timeseries arrays; `trial_index` + `step_order` dropped; `protocol_metadata` subgroup removed, `block_sequence` → `protocol_block_sequence` JSON; `99_Unrouted` fallback)
- **Last reviewed:** 2026-06-15

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

| Physical quantity | Spelled-out unit | Example field |
|---|---|---|
| Length / distance | `millimeter` | `clamp_separation_millimeter` |
| Mass | `gram` | `body_mass_gram` |
| Temperature | `celsius` | `temperature_celsius` |
| Angle | `degree` | `angle_commanded_degree` |
| Angular velocity | `degree_per_second` | `velocity_degree_per_second` |
| Strain | `percent` (positive = lengthening, negative = shortening) | `strain_percent` |
| Time | `second` | `time_second` |
| Frequency | `hertz` | `sample_rate_hertz` |
| Force | `newton` | `force_x_newton` |
| Torque | `newton_meter` | `torque_z_newton_meter` |
| Stiffness | `newton_per_millimeter` (or context-appropriate spelled-out derived unit) | — |
| Voltage | `volt` | `stim_command_volt` |
| Dimensionless ratio | no suffix | `time_normalized` |

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

---

## 3. `timeseries`

Continuous streams, sample-indexed at the acquisition rate, sharing one `time_second` axis.

> **Phase 0 status.** The canonical **dataset names** below are implemented in `bender_h5_export.py`, and the index/stim fields are now **per-sample timeseries arrays** (§3b-i). The on-disk group skeleton is still `02_TimeSeries/trial_XXXX` (one subgroup per step). The full §2 restructure — collapse to a single continuous `timeseries` group + `index_step_*` ranges, drop the numeric `01_/02_` prefixes — is a later, separate step. The per-sample index arrays are designed to concatenate cleanly into that flat timeline; the `index_step_*` range form (§4) is not yet written.

### 3a. Motor and time axes

| Dataset | Shape | Notes |
|---|---|---|
| `time_second` | `[N]` | Time axis in seconds. *(Was `t`; renamed in Phase 0.)* |
| `time_normalized` | `[N]` | Normalized time axis, dimensionless. *(Was `tnorm`.)* |
| `angle_commanded_degree` | `[N]` | Commanded motor angle (degrees). *(Was `angle_cmd`.)* |
| `angle_measured_degree` | `[N]` | Encoder angle, gearbox output / specimen frame (degrees). *(Was `angle_measured`.)* |
| `angular_velocity_commanded_degree_per_second` | `[N]` | Commanded angular velocity. *(Was `anglevel_cmd`.)* |

### 3b. Stimulator command streams

| Dataset | Shape | Notes |
|---|---|---|
| `stim_channel1_command_volt` | `[N]` | Stim command waveform, channel 1 (volts). *(Was `S1stimcmd`.)* |
| `stim_channel2_command_volt` | `[N]` | Stim command waveform, channel 2 (volts). *(Was `S2stimcmd`.)* |
| `stim_side` | `[N]` | Per-sample categorical: `none` / `left` / `right` / `both` (handles co-contraction). |
| `stim_state` | `[N]` | Per-sample stim phase: `passive` / `off` / `on` (3-value string). |
| `stim_type` | `[N]` | Per-sample activity: `active` / `passive`. |
| `instantaneous_frequency_hertz` | `[N]` | Chirp instantaneous frequency (frequency_sweep only). *(Was `sweep_instantaneous_freq`.)* |

> **Stim fields — RESOLVED (PI).** `stim_state` stays a 3-value string (`passive`/`off`/`on`). `stim_type` stays per-sample (`active`/`passive`) — `active_passive` mixes both within one experiment. `stim_side` already exists as the per-sample categorical above (`none`/`left`/`right`/`both`). No code change.

### 3b-i. Per-sample index streams (LOCKED — converted from per-trial scalar attrs)

Each is a length-`[N]` array; within the current `trial_XXXX` skeleton the step-protocol fields are constant across a trial's samples (broadcast) and concatenate cleanly when the timeline is later flattened.

| Dataset | Shape | Notes |
|---|---|---|
| `cycle_index` | `[N]` | Cycle number per sample; `-1` = not a numbered cycle (all step-protocol samples). *(Was `cycle_index_by_sample`.)* |
| `step_index` | `[N]` | Ramp/step number at each sample (step protocols only). |
| `sequence_index` | `[N]` | Pre-shuffle parameter index at each sample (step protocols only). |
| `block_index` | `[N]` | Block number at each sample (step protocols only). |
| `block_direction` | `[N]` | Bending direction for the block (string; step protocols only). |
| `block_stim_sides` | `[N]` | Stim sides for the block (string; step protocols only). |

**Dropped:** `trial_index` (redundant with the file/trial identity + `sequence_index`) and `step_order` (reconstructable from per-sample `sequence_index` + `block_index`). Dynamic / frequency_sweep carry only `cycle_index` (no discrete shuffled steps/blocks).

### 3c. Raw analog input buffer

| Dataset | Shape | Notes |
|---|---|---|
| `aidata` | `[8 × N]` | Raw analog input buffer. **IMMUTABLE** — written as-acquired, never decoded or modified at write time. Channel identities live in `metadata/daq_ai_channel_map`. |

### 3d. ATI 6-axis force/torque (decoded, in `timeseries`)

Raw + calibration principle (see Governing Rule 2): both the immutable raw voltages **and** the calibrated physical-unit data live in the raw file. The raw buffer `aidata` rows 0–5 are the immutable original; `forcetorque_raw` is the calibrated `[6 × N]` array computed at write time from `aidata` + `metadata/calibration_forcetorque_matrix`.

| Dataset | Shape | Notes |
|---|---|---|
| `forcetorque_raw` | `[6 × N]` | Calibrated F/T in physical units, **unmodified post-write** (immutable decoded original). **Kept as a single `[6 × N]` array — NOT split into per-channel datasets (PI decision).** |

> **Force/torque channel order — LOCKED (PI decision).** `forcetorque_raw` rows are, in order: `x_force`, `y_force`, `z_force` (newton), then `x_torque`, `y_torque`, `z_torque` (newton_meter). The Python writer does **not** split or rename these into per-channel datasets; the R pipeline slices `forcetorque_raw` and names/decodes the channels downstream. Row index → channel: `0=x_force, 1=y_force, 2=z_force, 3=x_torque, 4=y_torque, 5=z_torque`.

**Not in `timeseries`:** inertial-corrected torque series, `primary_torque_*`, and any computed derived quantities. These are regenerated in R from `forcetorque_raw` + correction parameters and written only to `derived` in the hub copy.

---

## 4. `metadata`

Flat. **No subgroups.** Every key carries a descriptive prefix; alphabetical sort groups related fields.

### `calibration_` — sensor → physical-unit conversions

```
calibration_forcetorque_matrix          [6×6]   ATI calibration matrix (raw voltages → newton / newton_meter)
                                                  (Current: calibration; rename pending)
calibration_sono_left                   [4 or 8] Sono calibration breakpoints, left channel
calibration_sono_right                  [4 or 8] Sono calibration breakpoints, right channel
```

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
```

### `index_` — structural indices into the continuous `timeseries`

Parallel arrays, **one element per step**. Replaces old per-trial subgroups (`trial_index`).

```
index_step_sample_start                 [n_steps]
index_step_sample_end                   [n_steps]
index_step_type                         [n_steps]
index_step_target_angle_degree          [n_steps]   (Current: index_step_target_angle_deg)
index_step_target_value_native          [n_steps]
index_step_curvature_per_meter          [n_steps]   (Current: index_step_curvature_1_per_m)
index_step_stim_voltage_left_volt       [n_steps]
index_step_stim_voltage_right_volt      [n_steps]
index_step_stim_t0_second               [n_steps]
index_step_stim_t1_second               [n_steps]
index_step_t_active_start_second        [n_steps]
index_step_t_active_end_second          [n_steps]
... (one column per per-step realized parameter currently in trial_index)
```

### `inertial_` — inertial-correction parameters + provenance

Parameters only. Inertial-torque **time series** are derived (hub `derived`), not here.

```
inertial_specimen_moi                   specimen moment of inertia (parameter)
inertial_system_moi                     full system MOI (parameter)
inertial_total_moi                      total MOI including clamp (parameter)
inertial_moi_provenance                 which frustum dims / calibration fed each MOI
```

### `measurement_` — ALL spatial geometry (animal AND apparatus)

**Any distance or length → `measurement_`**, regardless of whether it describes the specimen or the fixture.

```
measurement_specimen_frustum_heights_millimeter     [3]   (Current: specimen_geometry_heights_mm)
measurement_specimen_frustum_depths_millimeter      [3]   (Current: specimen_geometry_depths_mm; the "width" stations)
measurement_specimen_frustum_positions_millimeter   [3]   (Current: specimen_geometry_positions_mm)
measurement_specimen_bodylength_millimeter                 (Current: fishlen_TL — total length)
measurement_specimen_standardlength_millimeter            (Current: fishlen_SL — standard length)
measurement_specimen_body_width_millimeter                ← TEST-SECTION BODY WIDTH (was missing — now written)
measurement_specimen_body_mass_gram                       (Current: fishmass)
measurement_clamp_separation_millimeter                   (Current: dclamp)
measurement_clamp_offset_vertical_millimeter              (Current: dvert)
measurement_clamp_offset_horizontal_millimeter            (Current: dhoriz)
measurement_clamp_plate_extension_millimeter              (Current: clamp_plate_extension_mm)
measurement_target_muscle_depth_millimeter                (Current: target_muscle_depth_mm)
measurement_xsec_width_millimeter                         (Current: xsec_width)
measurement_xsec_height_millimeter                        (Current: xsec_height)
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

Global defaults. Per-step *realized* values live in `index_step_*`.

```
protocol_type                           file-level protocol selector (was: test_type)
protocol_starting_angle_degree          (was: starting_angle_deg)
protocol_duration_second                (was: duration)
protocol_isometric_initial / _final / _num_steps / _target_unit
protocol_isovelocity_min_velocity / _max_velocity / _starting_strain / _*_unit
protocol_step_frequency_hertz / _step_amplitude / _step_amplitude_unit / _step_curvature_per_meter
protocol_amplitude_frequency_exponent / _velocity_exponent
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

> **FLAG E — RESOLVED.** Rename `genus_species` → `specimen_genusspecies` approved (compound word, no underscore). Pending rename in `bender_functions.py` Python attr, GUI widget key, and `bender_h5_export.py`.

---

## 5. Raw + decoded principle (LOCKED)

1. **Always store both.** Whenever raw sensor data and a calibration/correction exist, store **both** the raw (immutable original) and the decoded (calibrated) form. `forcetorque_raw` is the canonical example.
2. **Raw is immutable.** Once written, a raw array is never overwritten, modified, or deleted. It is the permanent audit trail.
3. **Decoded is computed at write time** from the live calibration and written alongside the raw.
4. **Downstream always re-derives from raw + parameters.** The R pipeline reads `forcetorque_raw` + `calibration_forcetorque_matrix` and recomputes corrected/inertial series. It never modifies the raw array.
5. **Rationale:** signal integrity, recalibration without re-acquisition, audit trail, reproducibility.

In `timeseries`, this means:
- `aidata` (raw ADC voltages) — immutable, always present
- `forcetorque_raw` (calibrated F/T, decoded at write time) — immutable post-write
- Inertial-corrected torque — hub `derived` only, never in the raw file

---

## 6. Current-file → target mapping (audited keys)

From the audit of `2026-06-04_bass13_bender_24_isometric.h5` and the sim-validation session of `2026-06-12`. Resolves the known `metadata` payload:

| Current key (location) | Status | Target key | Action |
|---|---|---|---|
| `calibration` (bender_settings) | present | `calibration_forcetorque_matrix` | rename |
| `sono_cal_left` | present | `calibration_sono_left` | rename |
| `sono_cal_right` | present | `calibration_sono_right` | rename |
| `genus_species` (bender_settings) | present | `specimen_genusspecies` | rename |
| `prep_condition` (bender_settings) | present | `specimen_prep_condition` | rename |
| `specimen_sex` (bender_settings) | present | `specimen_sex` | conformant |
| `specimen_muscle_type` (bender_settings) | present | `specimen_muscle_type` | conformant |
| `session_analyst` (bender_settings) | present | `session_analyst` | conformant |
| `session_date` (01_Metadata top attrs) | present | `session_date` | conformant; schema rename to `session_date` pending (already correct) |
| `apparatus_id` (01_Metadata top attrs) | present | `session_apparatus_id` | rename |
| `simulated` (01_Metadata top attrs) | present | `session_simulated` | rename |
| `specimen_geometry_heights_mm` | present | `measurement_specimen_frustum_heights_millimeter` | rename |
| `specimen_geometry_depths_mm` | present | `measurement_specimen_frustum_depths_millimeter` | rename |
| `specimen_geometry_positions_mm` | present | `measurement_specimen_frustum_positions_millimeter` | rename |
| `fishlen_TL` | present | `measurement_specimen_bodylength_millimeter` | rename |
| `fishlen_SL` | present | `measurement_specimen_standardlength_millimeter` | rename |
| `fishmass` | present | `measurement_specimen_body_mass_gram` | rename |
| `dclamp` | present | `measurement_clamp_separation_millimeter` | rename |
| `dvert` | present | `measurement_clamp_offset_vertical_millimeter` | rename |
| `dhoriz` | present | `measurement_clamp_offset_horizontal_millimeter` | rename |
| `clamp_plate_extension_mm` | present | `measurement_clamp_plate_extension_millimeter` | rename |
| `xsec_width` | present | `measurement_xsec_width_millimeter` | rename |
| `xsec_height` | present | `measurement_xsec_height_millimeter` | rename |
| `target_muscle_depth_mm` | present | `measurement_target_muscle_depth_millimeter` | rename |
| `inertial_torque_specimen_primary` | empty `(0,)` | `inertial_specimen_moi` (param) | rename + reclassify |
| `inertial_torque_system_primary` | empty `(0,)` | `inertial_system_moi` | rename + reclassify |
| `inertial_torque_total_primary` | empty `(0,)` | `inertial_total_moi` | rename + reclassify |
| `primary_torque_raw` | empty `(0,)` | — | delete from metadata (derived series → hub) |
| `primary_torque_corrected` | empty `(0,)` | — | delete from metadata (derived series → hub) |
| `stim_clamp_notices` | empty `(0,)` | `note_bench` | rename |
| `trial_index/*` (subgroup) | present (minus `trial_index`/`cycle_index` cols) | `index_step_*` parallel arrays | move + flatten (later) |
| `protocol_metadata/step_order` | **dropped** | — | done (reconstructable from per-sample `sequence_index` + `block_index`) |
| `protocol_metadata/block_sequence` | **`protocol_block_sequence` (JSON)** | — | done (routed, `json.dumps`) |
| `protocol_metadata` (subgroup) | **removed** | — | done (attr-mirrors now canonical; redundant dict-only keys dropped) |
| `calibration_link` (empty group) | empty | — | delete |
| `daq_ai_sample_rate_hz` | present | `daq_ai_sample_rate_hertz` | rename |
| `t` (timeseries) | present | `time_second` | rename |
| `tnorm` (timeseries) | present | `time_normalized` | rename |
| `angle_cmd` (timeseries) | present | `angle_commanded_degree` | rename |
| `angle_measured` (timeseries) | present | `angle_measured_degree` | rename |
| `anglevel_cmd` (timeseries) | present | `angular_velocity_commanded_degree_per_second` | rename |
| `S1stimcmd` (timeseries) | present | `stim_channel1_command_volt` | rename |
| `S2stimcmd` (timeseries) | present | `stim_channel2_command_volt` | rename |
| `forcetorque` (timeseries) | **dropped** | — | done (identical copy of `forcetorque_raw`; not split, PI decision) |
| `forcetorque_raw` (timeseries) | present | `forcetorque_raw` `[6×N]` | conformant (keep, immutable, not split) |

---

## 7. Governing rules

1. **Immutability.** `aidata` (and all raw streams) are written as-acquired and never modified post-hoc. All correction happens in R and lands in `derived` in the hub copy. R is the authoritative correction path.
2. **Raw + decoded, both in file.** Whenever raw data + calibration parameters exist, store both: raw (immutable original) and decoded (computed at write). Never overwrite the raw. See §5 for the full rule.
3. **Decode via metadata, not per-channel.** Channel identity lives in `daq_ai_channel_map`. Decoded F/T channels in `timeseries` (§3d) are the single decoded form; the export does not additionally hard-slice `aidata` into named channels.
4. **Distance rule.** Any length/distance → `measurement_`, animal or apparatus. `specimen_` holds only non-spatial facts (id, species, sex, muscle type — no geometry, no mass).
5. **Step parameters appear twice by design.** Global defaults in `protocol_*`; per-step realized values in `index_step_*`. Intentional, not duplication.
6. **No per-step subgroups.** One continuous time series + `index_step_*` ranges. `trial_*` group naming is retired.
7. **Two files, shared group names.** Raw file = `metadata` + `timeseries`. Hub copy = same two groups + `derived`. Group names identical across both.
8. **Don't reuse one variable for two physical quantities.** Curvature, strain, angle, width, frequency are distinct even though all floats. Verify each assignment by physical role.
9. **Naming convention is locked.** `<description>_<unit>`, fully spelled-out units, no abbreviations. See §1.

---

## 8. Resolved flags / pending code changes

All schema-level decisions are locked. Remaining work is code-side renames only:

| Flag | Decision | Code change pending |
|---|---|---|
| A — ATI units | newton / newton_meter (SI) | channel names already use `newton` / `newton_meter` in this doc |
| B — stim channel naming | `stim_channel1_command_volt` / `stim_channel2_command_volt` | rename `S1stimcmd` / `S2stimcmd` in export |
| C — body length keys | `measurement_specimen_bodylength_millimeter` (TL) + `measurement_specimen_standardlength_millimeter` (SL) | rename `fishlen_TL` / `fishlen_SL` in export |
| D — simulated key | `session_simulated` | one-line rename in `bender_h5_export.py` |
| E — species key | `specimen_genusspecies` | rename `genus_species` in `bender_functions.py`, GUI, export |

- **Decision 2 — body-width timing.** Is the rig collecting fish *before* the schema rewrite ships? If yes, fix `measurement_specimen_body_width_millimeter` as a standalone commit on the current structure first.

---

## 9. Backlog (logged, not blocking)

- Update `tier2_data_curation_sop` + `jlab_folder_architecture` to the new key names; update any provenance assertions (`apparatus_id` → `session_apparatus_id`, etc.).
