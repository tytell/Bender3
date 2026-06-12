# Bender / CritterGripper — HDF5 Schema

## Authority metadata
- **Purpose:** Defines the on-disk structure, group layout, key-naming taxonomy, and unit conventions for Bender HDF5 files.
- **Authority level:** Tier 2 (binding for all H5 write/read tasks — export routine, R pipeline, GUI→backend assignment).
- **Scope:** `bender_h5_export.py` (write), `01_calibration.R` / `02_correct.R` / `03_analyze.R` (read), GUI variable assignment in `bender_streamlit_gui.py`.
- **Owner:** PI (schema decisions are PI-only).
- **Version:** v2.1
- **Last reviewed:** 2026-06-12

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

Continuous streams, sample-indexed at the acquisition rate, sharing one `time_second` axis. Step boundaries are encoded in `index_step_*` in `metadata` (§4). **No per-trial or per-step subgroups.**

### 3a. Motor and time axes

| Dataset | Shape | Notes |
|---|---|---|
| `time_second` | `[N]` | Time axis in seconds. *(Current name: `t` — rename pending.)* |
| `time_normalized` | `[N]` | Normalized time axis, dimensionless. *(Current: `tnorm`.)* |
| `angle_commanded_degree` | `[N]` | Commanded motor angle (degrees). *(Current: `angle_cmd`.)* |
| `angle_measured_degree` | `[N]` | Encoder angle, gearbox output / specimen frame (degrees). *(Current: `angle_measured`.)* |
| `angular_velocity_commanded_degree_per_second` | `[N]` | Commanded angular velocity. *(Current: `anglevel_cmd`.)* |

### 3b. Stimulator command streams

| Dataset | Shape | Notes |
|---|---|---|
| `stim_channel1_command_volt` | `[N]` | Stim command waveform, channel 1 (volts). *(Current: `S1stimcmd`.)* |
| `stim_channel2_command_volt` | `[N]` | Stim command waveform, channel 2 (volts). *(Current: `S2stimcmd`.)* |
| `stim_side` | `[N]` | Per-sample categorical label. *(Candidate to fold into `index_step_*` — see §7.)* |
| `stim_state` | `[N]` | Per-sample categorical label. |
| `stim_type` | `[N]` | Per-sample categorical label. |

### 3c. Raw analog input buffer

| Dataset | Shape | Notes |
|---|---|---|
| `aidata` | `[8 × N]` | Raw analog input buffer. **IMMUTABLE** — written as-acquired, never decoded or modified at write time. Channel identities live in `metadata/daq_ai_channel_map`. |

### 3d. ATI 6-axis force/torque (decoded, in `timeseries`)

Raw + calibration principle (see Governing Rule 2): both the immutable raw voltages **and** the calibrated physical-unit channels live in the raw file. The raw buffer `aidata` rows 0–5 are the immutable original; the decoded channels below are computed at write time from `aidata` + `metadata/calibration_forcetorque_matrix`.

| Dataset | Shape | Notes |
|---|---|---|
| `forcetorque_raw` | `[6 × N]` | Calibrated F/T in physical units (newton / newton_meter), **unmodified post-write** (immutable decoded original). *(Current: `forcetorque_raw`.)* |
| `force_x_newton` | `[N]` | x-axis force (newton). *(Decoded slice; current: `forcetorque[0,:]`.)* |
| `force_y_newton` | `[N]` | y-axis force (newton). |
| `force_z_newton` | `[N]` | z-axis force (newton). |
| `torque_x_newton_meter` | `[N]` | x-axis torque (newton_meter). |
| `torque_y_newton_meter` | `[N]` | y-axis torque (newton_meter). |
| `torque_z_newton_meter` | `[N]` | z-axis torque (newton_meter). Primary bending axis in most protocols. |

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

### `daq_` — acquisition / DAQ configuration

```
daq_ai_channel_map                      channel index → identity (e.g. ai0=Fx … ai6=sono)
daq_ai_sample_rate_hertz                AI + encoder sample clock (hertz) (Current: daq_ai_sample_rate_hz)
daq_ao_do_sample_rate_hertz             AO stim + DO motor sample clock (hertz)
daq_instrumentation                     list of active sensors this session (was: instrumentation)
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
protocol_name                           (was: protocol)
protocol_isometric_ramp_duration_second
protocol_isometric_target_angle_degree
protocol_isovelocity_velocity_degree_per_second
protocol_dynamic_frequency_hertz
protocol_dynamic_amplitude_degree
... (per protocol type)
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
| `trial_index/*` (subgroup) | present | `index_step_*` parallel arrays | move + flatten |
| `protocol_metadata/step_order` | present | resolve to `index_step_*` or `protocol_*` | decision pending |
| `calibration_link` (empty group) | empty | — | delete |
| `daq_ai_sample_rate_hz` | present | `daq_ai_sample_rate_hertz` | rename |
| `t` (timeseries) | present | `time_second` | rename |
| `tnorm` (timeseries) | present | `time_normalized` | rename |
| `angle_cmd` (timeseries) | present | `angle_commanded_degree` | rename |
| `angle_measured` (timeseries) | present | `angle_measured_degree` | rename |
| `anglevel_cmd` (timeseries) | present | `angular_velocity_commanded_degree_per_second` | rename |
| `S1stimcmd` (timeseries) | present | `stim_channel1_command_volt` | rename |
| `S2stimcmd` (timeseries) | present | `stim_channel2_command_volt` | rename |
| `forcetorque` (timeseries) | present | `force_x/y/z_newton` + `torque_x/y/z_newton_meter` | split + rename |
| `forcetorque_raw` (timeseries) | present | `forcetorque_raw` | conformant (keep, immutable) |

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
- Resolve whether per-sample `stim_side` / `stim_state` / `stim_type` streams stay in `timeseries` or fold into `index_step_*` (currently redundant with step ranges + stim params).
- Confirm `protocol_metadata/step_order` target (`index_step_*` vs `protocol_*`).
