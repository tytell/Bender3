# Bender → HDF5 Routing — Decision Log

- **Spec module:** `bender_routing_spec.py` (`BENDER_ROUTING`, `EXCLUDED`, `MISSING_REQUIRED`)
- **Schema:** `context_jlab_cg_h5schema.md` v2.1
- **Status:** DRAFT — built from a full trace of `bender_functions.py` (`__init__`, `record_motor_signal`, `run_experiment`, `run`, inertial-model builders, `update_metadata`) and the GUI assignment paths. **No exporter code changed.**
- **Contract:** every public (non-`_`) Bender attribute must be in `BENDER_ROUTING` **or** `EXCLUDED`; anything else makes the future exporter raise (`assert_full_coverage`).

---

## Coverage check (verified)

| Instance state | routed | excluded | unmapped |
|---|---|---|---|
| Fresh `Bender(config_A)` | 70 | 44 | **0** |
| Configured (protocol params + runtime arrays set) | 90+ | 42+ | **0** |
| `BENDER_ROUTING` total entries | 132 | — | — |
| `EXCLUDED` total entries | — | 72 | — |
| ROUTING ∩ EXCLUDED overlap | — | — | **∅** |

`BENDER_ROUTING` (130) and `EXCLUDED` (74) are a deliberate superset of any single
instance: they cover attributes that appear only after a run, only after specific
protocol selection, or only on certain inertial models. Coverage on a live instance
is always a subset, and both tested states report **0 unmapped**.

---

## MC questions — answers recorded (CONFIRMED by user)

| # | Field(s) | Decision (confirmed) | Rationale |
|---|---|---|---|
| 1 | `post_trial_notes` | → `note_posthoc` | QC/post-experiment text. Avoids collision with `stim_clamp_notices → note_bench` (schema §6). |
| 2 | inertial MOI / mass | `specimen_moi_{specimen,profile,frustum} → inertial_specimen_moi` (first present); `i_total_system → inertial_total_moi`; `total_mass → EXCLUDED` | Matches schema §4 inertial params; mass is not a schema field and is recomputed in R. |
| 3 | `starting_angle_deg` | → `protocol_starting_angle_degree` | A protocol setting, not session logistics. |
| 4 | `temp_C_room` / `temp_C_tank` | → `session_temperature_room_celsius` / `session_temperature_tank_celsius` | Session environmental facts; not spatial → not `measurement_`. |
| 5 | `segment` | → `specimen_segment` | Anatomical identity field. |
| 6 | per-sample `stim_side/stim_state/stim_type` | kept as `timeseries` datasets | Derived per-sample categorical labels from the **commanded** S1/S2 (NOT `stim_monitor`/`ai7`). Confirmed to stay full timeseries. |
| 7 | protocol param families | **spelled-out renames** under `protocol_*`, value + `_unit` (+ `_input_mode`) pattern; family prefix kept | See "Protocol value/unit convention" below. |
| 8 | geometry aliases | `dbend` + `test_segment_position_mm → measurement_bend_position_millimeter`; `test_segment_length_mm → EXCLUDED` | Bend/segment position is spatial → `measurement_`. `test_segment_length_mm` is the **same quantity as `dclamp`** (kept in sync by `_sync_dclamp_test_segment_aliases`), already routed to `measurement_clamp_separation_millimeter`; routing the alias would duplicate the distance and break Rule 4. |

> **To change any decision:** edit the single matching row in `BENDER_ROUTING`
> (or move it to/from `EXCLUDED`) and rerun `assert_full_coverage`.

---

## Protocol value/unit convention (MC #7, confirmed)

**File-level selector.** The whole file is one protocol (`run_experiment` sets one
`test_type`). It is written once as **`protocol_type`** (e.g. `"isometric"`). The
`protocol_` prefix remains the umbrella for all recipe fields; `protocol_type` is the
category selector (more accurate than `protocol_name`; cleaner umbrella than `test_mode`).

**Value + unit pair.** Parameters whose unit is mode-dependent are split:

| Role | Field | Example |
|---|---|---|
| value (numbers) | `protocol_isometric_initial` | `-5.0` |
| resolved unit | `protocol_isometric_target_unit` | `"percent"` |
| raw convention (provenance, option B) | `protocol_isometric_target_input_mode` | `"strain"` |

The `*_mode` Bender attr is resolved at write time to a spelled-out SI unit
(`degree` / `percent` / `per_meter` / `degree_per_second`), and the raw input convention
(`strain` vs `strain_pct` vs `angle`…) is preserved in `_input_mode`. One `_unit`/`_input_mode`
pair covers every value sharing the same mode.

Mode→pair groups applied:
- `all_amps` (value) ↔ `all_amps_mode` → `protocol_step_amplitude` + `protocol_step_amplitude_unit`
- `curve_input_values` ↔ `curve_input_mode` → `protocol_curve_input` + `protocol_curve_input_unit`
- `isometric_initial/final/mirror_target_*` ↔ `isometric_mode` → `protocol_isometric_target_unit`
- `isovelocity_min_vel/max_vel` ↔ `isovelocity_velocity_mode` → `protocol_isovelocity_velocity_unit`
- `isovelocity_starting_strain` ↔ `isovelocity_starting_strain_mode` → `protocol_isovelocity_starting_strain_unit`

Fixed-unit params keep a direct spelled key (`all_freqs → protocol_step_frequency_hertz`,
`all_curves → protocol_step_curvature_per_meter`, `stim_pulse_rate → protocol_stim_pulse_rate_hertz`, …).

**Analysis-safety (confirmed):** these are protocol *metadata*, not measured signals. Raw
arrays (`aidata`, `forcetorque_raw`) and canonical-unit timeseries are untouched. Explicit
number + spelled unit is strictly more self-describing for R/AI readers, and resolving to a
canonical unit makes units consistent across files. The `index_step_*` realized values should
follow the same value+unit pattern for consistency (exporter work, noted below).

---

## Routing rationale by group

### `metadata / specimen_` — non-spatial identity
`fishcode→specimen_id`, `genus_species→specimen_genusspecies`, `prep_condition→specimen_prep_condition`,
`specimen_sex`, `specimen_muscle_type` (conformant), `segment→specimen_segment`.
Geometry deliberately excluded from this group per Rule 4.

### `metadata / measurement_` — ALL spatial geometry (animal + apparatus)
Direct from schema §4/§6 mapping table: body/standard length, mass, xsec width/height, clamp
separation/offsets/plate extension, target muscle depth, the three frustum arrays, and the
two inertial-profile spatial inputs (`profile_length`, `profile_clamp_offset`). Rule 4: any
distance → `measurement_`, regardless of animal vs fixture.

### `metadata / session_` — recording-session logistics
`session_date` (computed at run start), `apparatus_id→session_apparatus_id`, `session_analyst`,
`simulation_mode→session_simulated` (bool on every file), `simulation_material`,
`config_name→session_protocol_template`, `acquisition_start`, and the two temperatures
(`session_temperature_room/tank_celsius`).

### `metadata / calibration_`
`calibration→calibration_forcetorque_matrix` (6×6 → newton/newton_meter), `sono_cal_left/right`.

### `metadata / daq_`
Sample rates (spelled-out hertz); `input_channels` + `input_channel_names` → `daq_ai_channel_map`
(computed index:identity) and `daq_instrumentation` (active sensor list).

### `metadata / inertial_` — parameters + provenance only
MOI params (decision #2), calibration linkage (`use_inertial_calibration`, `inertial_calibration_file`,
`inertial_calibration_profile` dict→subgroup), `frustum_inputs` dict, and the profile model
inputs/densities. **No inertial-torque time series** — those are regenerated in R (`derived`, hub only).

### `metadata / note_`
`stim_clamp_notices→note_bench` (schema §6); `post_trial_notes→note_posthoc` (decision #1).

### `metadata / protocol_` — global defaults
`test_type→protocol_type` (file-level selector), motion/timing/stim defaults, recruitment/lateral
settings, and the dynamic / isometric / isovelocity / frequency-sweep param families spelled out
under the value + `_unit` (+ `_input_mode`) convention (decision #7; see section above).
Per-step *realized* values are separately destined for `index_step_*` (see MISSING_REQUIRED).

### `timeseries`
`t→time_second`, `tnorm→time_normalized`, `angle→angle_commanded_degree`,
`angle_measured→angle_measured_degree`, `anglevel→angular_velocity_commanded_degree_per_second`,
`S1stimcmd/S2stimcmd→stim_channel1/2_command_volt`, `aidata` (immutable 8×N),
`forcetorque_raw` (immutable decoded 6×N), `forcetorque` **split** into the six named
force/torque channels, decoded `sono_left/right_mm→sono_left/right_millimeter`,
`cycle_index_history→cycle_index_by_sample`, `sweep_instantaneous_freq→instantaneous_frequency_hertz`,
and the three per-sample stim categoricals (decision #6).

---

## EXCLUDED — why each class is not in the raw file

| Class | Members (examples) | Reason |
|---|---|---|
| Internal caches / motor state | all `_`-prefixed; `_motor_continuous_step_pos`, `_encoder_cumulative_deg`, … | Observability/control state, not data. |
| Objects / transport dicts | `master_logger`, `h5_protocol_metadata`, `trial_records` | Object or container; contents routed individually (records → timeseries + index). |
| Machine-specific | `outputfile`, `timestamp` | Absolute path / compact stamp; portable `filename` + `start_time_iso` written instead. |
| DAQ scratch buffers | `stimcmdhi`, `dig` | Placeholder NI buffers, not acquired data. |
| Derived series → hub | `primary_torque_raw`, `angledata`, `is_stim_cycle` | Recomputed in R; `angledata` duplicates `angle_measured`. |
| Legacy aliases | `initial/final/num_steps/mode`, `min_vel/max_vel/starting_strain`, `test_segment_*`, `dbend` | Mirrored onto canonical names; routing would duplicate. |
| Mass | `total_mass` | Not a schema field (decision #2). |
| Legacy frustum scalars | `frustum_height_mm`, … | Captured inside `frustum_inputs` → `inertial_frustum_inputs`. |
| Hardware-config provenance | `device_name`, `motor_*`, `encoder_*`, channel lists, `forcetorque_calibration_file`, timing defaults | Already written verbatim to `metadata` by the config-provenance writer; re-routing would double-write. **If you want any explicitly routed instead, move it into `BENDER_ROUTING`.** |

---

## MISSING_REQUIRED — schema fields with no Bender source

| Key | Intended source | Action needed |
|---|---|---|
| `measurement_specimen_body_width_millimeter` | GUI | Add a widget; distinct from `xsec_width` (test-section cut width). |
| `note_posthoc` | analysis pipeline | Populated downstream; or fed by `post_trial_notes` (decision #1). |
| `inertial_moi_provenance` | computed | Assemble from active inertial model + `frustum_inputs`. |
| `index_step_*` (parallel arrays) | computed | Derive `sample_start/end/type/target/...` from `trial_records` step boundaries. |

---

## Open items before locking

All 8 MC decisions + the protocol value/unit naming are **confirmed**. Remaining work is
exporter/GUI implementation, out of scope for this spec:

1. `MISSING_REQUIRED` wire-ups: `measurement_specimen_body_width_millimeter` (GUI widget), `index_step_*` derivation, `inertial_moi_provenance` assembly.
2. Exporter: resolve each `*_mode` attr to a spelled-out `_unit` and emit the `_input_mode` provenance sibling at write time.
3. Exporter: apply the same value+`_unit` pattern to `index_step_*` realized values for consistency.
4. Confirm the hardware-config provenance set should stay `EXCLUDED` (relying on the config writer) rather than being explicitly routed under `daq_`/`session_`.
5. Code-side schema renames still pending (`genus_species → specimen_genusspecies`, `simulated → session_simulated`) per schema §8.
