---
name: ""
overview: ""
todos: []
isProject: false
---

# HDF5 schema Phase-0 -> Phase-1 migration: flat groups + index_step_*

Deferred migration referenced in the segmented-finite plan's "Out of scope" section
(session E) and in `.cursorrules` HDF5 cleanup. Steps 1-5 of the segmented-finite
plan are committed (`840155d`..`4503a06`): `daq_collection_type`, `step_manifest`,
one-based `step_NNN` subgroups, `step_buffer_s` bookends, and the 2s reset floor are
live. This plan finishes the structural migration: lowercase flat `metadata`/`timeseries`
groups, drop the numeric `01`_/`02_` prefixes, flatten metadata subgroups, and add the
`index_step_*` parallel arrays.

Authority: Tier-2 schema doc `context_jlab_cg_h5schema.md` is ground truth (PI-owned).
Every step is self-contained: a fresh session can execute one step from this file alone
(files touched, the single commit it produces, validation, recommended model, and a
do-not-bundle boundary).

CAVEAT (carried from segmented-finite plan): Step 3 (`5d26659`) and Step 4 (`6a6fa11`)
are tagged [RIG-VALIDATION NEEDED] and not yet hardware-confirmed. This migration is
schema STRUCTURE, independent of per-step acquisition timing, so it proceeds regardless;
any step that re-touches the acquisition loop inherits that rig-validation debt (none of
the steps below do).

---

## Locked decisions (PI sign-off, 2026-06-19)

- **D1 - Group rename: GO.** `01_Metadata` -> `metadata`, `02_TimeSeries` -> `timeseries`.
Flat, lowercase, no numeric prefixes. This is the point of the migration.
- **D2 - Hard cut, no shim.** No dual-name reader. Old `01_/02_` files are read by
checking out the pre-migration commit if ever needed. Single operator + RA; no
compatibility guarantees. **Record the cutover commit hash in this file** (see
"Cutover record" below) once Step M1 lands so before/after is unambiguous.
- **D3 - calibration_link: relocate, don't delete.** The subgroup holds REAL captured
values (`use_inertial_calibration`, `calibration_file`, `calibration_available`,
`bender_h5_export.py:307-310`). Schema doc 6 implies deletion but conflicts with code
writing data -> flatten to flat `calibration_*` keys, matching the prefix taxonomy.
Keep only `calibration_inertia_file` (provenance); drop `calibration_inertia_used`
and `calibration_inertia_available` — no correction runs at export, so availability
flags are meaningless. Normalize "inertial" -> "inertia" throughout. Confirm name at
M0 (spec freeze), not M1.
- *D4 - index_step_ column set = code, not invention.** Enumerate exactly the per-step
scalar fields the exporter already collects (today written under
`01_Metadata/trial_index/<col>`, `bender_h5_export.py:501-531`). The schema doc's "..."
is non-exhaustive by design; Steps 1-5 already fixed the real fields. See "index_step_*
column inventory" below. Do NOT add new columns.
- **D5 - daq_collection_type: keep in-file.** Cheap, self-describing, human/AI-readable.
Keep alongside `protocol_sampling_mode` (the key the schema doc keys structure off).
Redundant but not in conflict; drop neither.
- **D6 - Root-attr mapping: approved.** `post-trial notes` -> `note_bench`/`note_posthoc`,
`start_time_iso` -> `session_date`, `test_type` -> `protocol_type`. Only genuinely
file-level attrs (not duplicating a flat metadata key) stay at root. M1 must print the
final surviving root-attr list for PI eyeball before finalizing. See "Root-attr
disposition" below.
- **D7 - 99_Unrouted: under metadata.** Move from root to `metadata/99_Unrouted` (or a
flat `unrouted_*` prefix - confirm at M1) so root holds only the two groups + true
file-level attrs.
- **D8 - R pipeline: 01/02 fresh, 03 PORT (amended 2026-06-22).** `01_calibration.R`
and `02_correct.R` are fresh (no existing R applies the FT matrix or does a
time-series inertial correction). `03_analyze.R` is a PORT of the validated
muscle/mechanics codebase (`bender_scup_muscle/R` + `ScupMechanics/R`) through one
new flat-schema loader — not a fresh author. `bender_data.R` reads Bender2-era groups
(`/Calibrated`, `/NominalStimulus`, `/RawInput`) that do not exist in this schema —
leave it untouched as a historical reader. Not a refactor target.
**M3 must follow D9** (field renames A-E) so R scripts are written once against
stable final field names; Session F moves before Session D.
- **D9 - Field renames: SEPARATE.** Do NOT fold schema-doc 8 Flags A-E
(`genus_species` -> `specimen_genusspecies`, `simulated` -> `session_simulated`,
`fishlen_TL/SL` -> `measurement_specimen_*length_millimeter`, etc.) into this migration.
One-logical-change-per-commit. Land the structural migration clean, then a follow-up
Phase-0 session does A-E once the new structure is stable. **D9 must execute before
M3** so `03_analyze.R` is written once against stable final field names.
- **D10 - single_finite per-cycle metadata block: ADD (PI sign-off 2026-06-22).**
The muscle work-loop / passive-stiffness analysis needs per-cycle nominal design
for `single_finite` combo/constant-freq trials. Per-sample `stim_side`, `stim_state`,
and `cycle_index` already exist (§3b/§3b-i) so stim labeling and cycle membership
need no new keys; only the per-cycle DESIGN GRID does. New parallel arrays
(`index_cycle_*`, length C = total bending cycles; join to samples via `cycle_index`)
  - protocol scalars:
  - `index_cycle_index` (1-based), `index_cycle_frequency_hertz`,
  `index_cycle_operating_point`, `index_cycle_operating_point_units`
  (user's drive amplitude in its native unit — matches `index_step_operating_point`
  pattern), `index_cycle_motor_amplitude_degree`, `index_cycle_active` (bool),
  `index_cycle_activation_duty`, `index_cycle_activation_phase`
  - `protocol_activation_start_cycle` (int, 0-based — FLAG: 0-based unlike
  `index_cycle_index`), `protocol_cycles_per_step`, `protocol_cycle_count`,
  `protocol_end_cycle_count`
  - `Lonoff`/`Ronoff` omitted (redundant with per-sample `stim_side`/`stim_state`).
  - Does NOT add `index_step_*` columns (D4 intact).
  - Exact final key names confirmed at M0. Implementation step deferred (see Deferred).
- **D11 amendment — derived/ group (PI sign-off 2026-06-22, Point 4 FINAL; IMPLEMENTED).**
  The exporter writes a `derived/` group into the HDF5 file at export time for live RA
  inspection. NOT ground truth; R ignores it and re-derives from raw. Contents:
  - `derived/forcetorque_calibrated` — ATI 6×6 matrix applied to the raw F/T rows
    (reuse the already-computed `self.forcetorque` / `apply_calibration_forcetorque`,
    `bender_functions.py:630/1113`). Always written when the matrix is present.
  - `derived/torque_inertia_corrected` — calibrated primary torque minus `I_sum·α`,
    where `I_sum` = the SUM of whichever inertias are present: apparatus (`I_est`)
    and/or specimen MOI. Written whenever at least one is present; silently skipped
    (no error) when neither is. (Gate is `I_est`/`specimen_moi` presence — NOT the
    dropped `system_moi`.)
  - The two datasets use SEPARATE try/except blocks so a failure in one does not
    suppress the other (RA gets partial output, not silence).
  - Unit guard (IMPLEMENTED): torque stays in N·m. `torque_inertia_corrected` reads
    `I_est` from the in-memory `inertial_calibration_profile` dict in its native
    `N·m/(deg/s²)` (NOT the converted stored g·mm² value) and converts specimen MOI
    g·mm² -> `N·m/(deg/s²)` via `*1e-9*(pi/180)` before summing; `α` is in deg/s². The
    dead, unit-inconsistent `get_corrected_torque` is NOT reused.
- **D12 - calibration_inertia_ consolidation (PI sign-off 2026-06-22; FINALIZED with
  Points 1-3; IMPLEMENTED).** All computed MOI lives under `calibration_inertia_*`; raw
  measured inputs (geometry/mass/density) stay under `measurement_specimen_*`. Units live
  in the KEY SUFFIX (no `attrs['units']`):
  - `calibration_inertia_specimen_moi_gram_millimeter_squared` (geometry-derived inline
    at export; replaces `inertial_specimen_moi`). Live path is the analytic N-station
    integral (`set_specimen_geometry_inertial_model`, `bender_functions.py:2225`) — exact
    per-segment, NO `num_samples` knob; the station lists are the parameters.
  - `calibration_inertia_specimen_from_geometry` (bool — analytic MOI available).
  - `calibration_inertia_total_moi_gram_millimeter_squared` (replaces `inertial_total_moi`;
    = THEORETICAL apparatus baseline + specimen, sourced from `i_total_system` — NOT
    apparatus_moi + specimen_moi; see Deferred flag 5).
  - Flatten `inertial_calibration_profile` subgroup (Point 2):
    `I_est` -> `calibration_inertia_apparatus_moi_gram_millimeter_squared` (unit-convert
    `*(180/pi)*1e9` at write time — I_est is `N·m/(deg/s²)`, a per-degree MOI);
    `bias_est` -> `calibration_inertia_bias_newton_meter`;
    `axis_sensor` -> `calibration_inertia_axis_sensor`.
  - DROP `calibration_inertia_system_moi` (Point 3 — phantom, no source attr; "apparatus"
    now means rig-alone = I_est).
  - DROP `calibration_inertia_moi_provenance` (Point 1 — geometry/mass/density are now
    first-class `measurement_specimen_*` fields).
  - MOVE specimen density `specimen_geometry_density_g_per_mm3` ->
    `measurement_specimen_density_gram_per_cubic_millimeter` (measured input, Point 1).
  - DROP from schema/routing the flags `use_theoretical_inertial_correction` (gates
    nothing) and `use_frustum_inertial_model` (notebook auto-build switch only; removing
    the attribute itself is a separate optional code cleanup) (Point 1).
  - `calibration_inertia_matrix` (apparatus-MOI lookup table; empty until apparatus-
    calibration workflow implemented — see Deferred flag 4).
  - `calibration_inertia_file` (source file; replaces `calibration_inertial_file`).
  - DROP `calibration_inertia_used`, `calibration_inertia_available` (no correction at export).
  - Normalize "inertial" -> "inertia" in ALL field names and plan text.

---

## index_step_* column inventory (D4)

Source of truth: the per-step `entry` dicts in `_run_force_length_steps`
(`bender_functions.py:4052-4118`, isometric) and `_run_isovelocity_steps`
(`bender_functions.py:4913-4987`, isovelocity), as currently flattened into the
`trial_index` subgroup by the exporter (`bender_h5_export.py:501-531`). The migration
re-emits these as flat `metadata/index_step_<name>` parallel arrays (one element per step).

Scalar per-step fields present today (drop the timeseries arrays `t`/`angle_cmd`/
`anglevel_cmd`/`tnorm`/`S1stimcmd`/`S2stimcmd`/`aidata`/`angle_measured`/`forcetorque*`/
`primary_torque_raw`, which are timeseries, not index columns):

- Shared: `step_index`, `recruitment`, `unilateral_posture_lateral_index`,
`motor_positive_bend_toward_lateral_index`, `bilateral_mirror_motor`,
`stim_t0`, `stim_t1`, `t_pre_baseline_start`, `t_pre_baseline_end`,
`t_active_start`, `t_active_end`, `t_post_baseline_start`, `t_post_baseline_end`,
`wall_clock_start`, `rest_before_second`, `duration_second`,
`operating_point`, `operating_point_units`
- Isometric-only: `target_deg`, `ramp_from_deg`
- Isovelocity-only: `velocity_deg_s`, `theta_start_deg`, `pre_stim_t0`, `pre_stim_t1`,
`iso_t0`, `iso_t1`, `mean_xforce_stim`, `guard_triggered`, `guard_angle_deg`,
`velocity_seg1_deg_s`, `velocity_seg2_deg_s` (mirror only)
- Block metadata (when `use_block_stim`): `block_index`, `block_direction`,
`block_stim_sides`, `left_stim_voltage`, `right_stim_voltage`
- Motor-position reference fields added by `_record_motor_position_reference(...)`
- DROPPED (schema decision, already removed from manifest): `trial_index`, `cycle_index`
per-step scalars, `step_order`.

Naming: prefix each with `index_step`_ and apply schema-doc spelled-out units where the
doc already names a target (e.g. `target_deg` -> `index_step_target_angle_degree`,
`duration_second` stays, `*_t0/_t1` -> `*_second`). M2 must produce the exact final
name map for PI confirmation - this inventory is the input set, not the final key strings.

---

## Root-attr disposition (D6)

Current root `f.attrs` (`bender_h5_export.py:294-298`):


| Root attr today    | Disposition                                                             | Target                                    |
| ------------------ | ----------------------------------------------------------------------- | ----------------------------------------- |
| `schema_version`   | KEEP at root (file-level provenance)                                    | `schema_version` (root)                   |
| `test_type`        | MOVE (duplicates protocol selector)                                     | `metadata/protocol_type` (already routed) |
| `post-trial notes` | MOVE                                                                    | `metadata/note_bench` (+ `note_posthoc`)  |
| `filename`         | DROP from root (already mirrored at `g_meta.attrs['filename']`, `:305`) | `metadata/filename`                       |
| `start_time_iso`   | MOVE                                                                    | `metadata/session_date` (verify routing)  |


After M1 the only root-level attr expected to survive is `schema_version`, plus the two
groups `metadata` and `timeseries`. M1 prints the actual surviving root set for eyeball.

---

## Steps (each = one commit, one model, do-not-bundle)

### M0 - Spec freeze (no code)

**Model:** Sonnet | **Validation:** PI review | **Do not bundle.**

- Update `context_jlab_cg_h5schema.md` so it literally matches D1-D12: confirm the
`index_step_`* final name map (input from the inventory above), the `calibration_*`
flat key names (D3), `metadata/99_Unrouted` placement (D7), root-attr survivors (D6),
D10 `index_cycle_*` + `protocol_*` final key names, D11 architecture (embedded matrix
values, `calibration_forcetorque_file`, removal of `forcetorque_raw` with fallback
guard, supersession of §3d/Governing Rule 2/F/T channel-order lock), and D12
`calibration_inertia_*` consolidation (drop `calibration_inertia_used`/`_available`).
Note: "inertial" -> "inertia" rename is resolved in M1 (structural flattening), not here.
- Bump schema version v2.4 -> v2.5 with a changelog entry enumerating D10, D11, D12,
D8 amendment, and the "inertia" normalization.
- Output: the schema doc is the build target for M1-M3. No `.py` changes.

**Commit:** `docs(schema): freeze Phase-1 flat-group target (metadata/timeseries, index_step_*, calibration_* flat, D10/D11/D12)`

---

### M1 - Group structure: rename + flatten + remove dead branch (writer + ALL readers)

**Model:** Opus | **Validation:** Mac (unit tests + sim export) | **Do not bundle with M2.**

- `bender_h5_export.py`: `create_group('01_Metadata')` -> `'metadata'`,
`create_group('02_TimeSeries')` -> `'timeseries'`; fix the flat-write `_tname`
literal (`:363`).
- Flatten `calibration_link` subgroup -> flat `calibration_inertia_file` attr only;
drop `use_inertial_calibration` and `calibration_available` (D3, D12).
- Flatten `inertial_calibration_profile` subgroup -> flat `calibration_inertia_`*
scalars per D12 (`:562-566`). Apply "inertial" -> "inertia" normalization to ALL
`inertial_*` keys M1 touches when flattening — this resolves Deferred flag 2 in
place rather than deferring it to D9 or M2b. D9 Flags A-E renames are still SEPARATE
(D9); only the inertia-prefix normalization rides here.
- Remove the dead `trial_{i:04d}` branch + stale comments (`:368-370`) - isometric/FL
are `segmented` since `bender_functions.py:4126`, so only `step_NNN` (segmented) and
flat (single_finite) remain.
- Reconcile root attrs per D6; move `99_Unrouted` under `metadata` per D7.
- Update EVERY H5 reader (hard cut, D2): `bender_h5_plot_helpers.py`,
`bender_h5_explore.py`, `validate_plot_h5_batch.py`, `plot_sono_length_vs_time.py`,
the `h5_explorer` route in `bender_streamlit_gui.py`, `bender_tk_gui.py`.
- Update tests: `tests/test_canonical_export_schema.py`, `tests/test_continuous_isometric.py`.
- Print the surviving root-attr set during a sim export for PI eyeball (D6).
- M1 is PURELY STRUCTURAL: LEAVE the `trial_index` subgroup intact and unchanged. After
M1 `metadata` is flat EXCEPT for the interim `trial_index` subgroup, which M2 removes as
part of the full `trial_index` -> `index_step_*` swap. Do not touch the `trial_index`
writer (`bender_h5_export.py:501-531`) in M1; readers of `trial_index` keep working until
M2. (Boundary confirmed by PI 2026-06-19.)

**Commit:** `refactor(h5): flat metadata/timeseries groups, drop 01_/02_ prefixes, flatten subgroups (hard cut)`
**Cutover record:** write the resulting commit hash into the "Cutover record" section below.

---

### M2 - full trial_index -> index_step_* swap (writer)

**Model:** Opus | **Validation:** Mac (unit tests + sim segmented export) | **Do not bundle.**

- M2 OWNS the entire `trial_index` -> `index_step_*` swap (PI decision 2026-06-19):
remove the `trial_index` subgroup writer (`bender_h5_export.py:501-531`) AND write the
replacement flat `metadata/index_step_<name>` parallel arrays in the same commit. After
M2 there is no `trial_index` subgroup and `metadata` is fully flat.
- Build `index_step_*` from the per-step `entry` scalars (inventory above). Apply the
final name map confirmed at M0.
- Update any reader that still reads `trial_index` (none should remain after this commit)
to read `index_step_*` instead.
- One element per step; string columns as `S` arrays, numeric as float/int, matching the
existing manifest-column logic.
- `step_manifest` (already JSON in metadata) is unchanged - it stays the cross-step
timing index; `index_step_*` carries the realized per-step parameters.
- Update `tests/test_canonical_export_schema.py` to assert flat `index_step_*` arrays and
the absence of the `trial_index` subgroup.

**Commit:** `feat(h5): write flat index_step_* per-step arrays; remove trial_index subgroup`

---

### M3 - R pipeline: 01/02 fresh, 03 PORT against the flat schema

**Model:** Opus | **Validation:** R reads a sim `single_finite` AND a sim `segmented_finite` file | **Do not bundle.**
**Depends on:** M0, M1, M2, M2b (D11 embedding; deferred), and D9 field renames landing first.

Existing R analysis corpus (Desktop/bender_projects: `bender_scup_muscle/R` +
`ScupMechanics/R`) is coupled to Bender2-era `/Calibrated`, `/NominalStimulus`,
`/RawInput` via one loader (`load_bender_data.R`) — NOT to `01_/02_/trial_XXXX`.
That loader is the only schema-coupled file; all downstream analysis math is reusable.

- **Shared seam:** ONE new flat-schema loader (replacing `load_bender_data.R`), reading
`timeseries/*` + flat `metadata/*` (including D10 `index_cycle_*` / `protocol_*` +
per-sample `stim_side`/`stim_state`/`cycle_index`). Port bias subtraction + 9th-order
Butterworth low-pass. Emits the same analysis-tibble column names the ported math
expects.
- **New `01_calibration.R` (mostly fresh):** reads the EMBEDDED
`calibration_forcetorque_matrix` VALUES (D11 — never an external `.cal` path) +
raw F/T rows of `aidata` identified via `daq_ai_channel_map`; reads
`calibration_sono_left`/`right` VALUE tables; reads MOI from
`calibration_inertia_*` (D12). Writes calibrated output to hub `derived/` only.
No cal-matrix math exists in the current R corpus — fresh.
- **New `02_correct.R` (fresh, NOT a port):** time-series inertial correction
`tau_corrected = tau_raw - I * theta_ddot` (`theta_ddot` from `timeseries` angle,
`I` from `calibration_inertia_*`); writes to hub `derived/` only, never the raw
file. Port `estimate_body_torque.R` axis-selection + torque QC as a pre-step.
No existing R does this series correction — budget as new numerical work.
- `**03_analyze.R` (PORT + ADAPT — bulk of reusable IP):** port `calc_work`,
`calc_muscle_torque`, `calc_muscle_metrics`, `calc_muscle_work_time`,
`set_cycle_types`, `stiffness_damping`, `get_mechanics_by_half_cycle` (finish
`work.Nm` stub), `interpolate_even_phase`, `estimate_muscle_mass`,
`detect_trial_type`, `filter_*` QC stack, Fig 3/4/5 plot stack. Must handle:
`single_finite` protocol branching (constant-freq / combo / sweep), freq x curvature
(+ duty x phase) grids, half-cycle ramp-edge handling, both passive-reference modes,
Coughlin + sinusoid + symmetry QC. Reads D10 `index_cycle_*` + per-sample
`stim_side`/`stim_state`/`cycle_index` for stim labeling and the passive template.
`segmented_finite` (isometric/isovelocity/FL/FV) analysis is fresh — parse
`step_manifest` + `index_step_*`.
- Leave `bender_data.R` and the Desktop/bender_projects corpus untouched (D8).

**Commit:** `feat(R): flat-schema loader + 01_calibration (embedded cal), 02_correct (inertial series), 03_analyze (ported muscle/mechanics)`

---

### M4 - End-to-end validation

**Model:** Opus | **Validation:** full round-trip | **Do not bundle.**

- Sim run -> `export_primary_h5` -> open in `bender_h5_explore` -> R `01/02/03` read ->
`derived/` round-trip, for one `segmented_finite` (isometric or isovelocity) and one
`single_finite` (dynamic or sweep) file.
- Confirm: flat groups, no `01_/02`_ prefixes, `step_NNN` one-based, `index_step_*`
present, root holds only `schema_version` + the two groups, R reads without path errors.
- Document any residual gaps as a follow-up (e.g. the deferred Flags A-E session, D9).

**Commit:** `test(h5): end-to-end Phase-1 schema validation (sim segmented + single_finite)`

---

## Execution order / session map


| Session | Step                                                | Model  | Hardware    | Commits       |
| ------- | --------------------------------------------------- | ------ | ----------- | ------------- |
| A       | M0                                                  | Sonnet | none (docs) | 1             |
| B       | M1                                                  | Opus   | Mac only    | 1             |
| C       | M2                                                  | Opus   | Mac only    | 1             |
| C2      | M2b — D11 embedding (deferred; no session assigned) | Opus   | Mac         | 1             |
| D       | Flags A-E field renames (D9)                        | Opus   | Mac         | separate plan |
| E       | M3                                                  | Opus   | Mac (R)     | 1             |
| F       | M4                                                  | Opus   | Mac         | 1             |


Note: Session D (D9 field renames) moved before Session E (M3) so `03_analyze.R`
is written once against stable final field names. M2b is deferred pending scope
assignment (see Deferred flag 1); it must land before M3 executes.

**Cross-session invariants:**

- Flat lowercase `metadata` / `timeseries`; no numeric prefixes anywhere in new code.
- Hard cut (D2): no dual-name reader; cutover hash recorded below.
- M1 is purely structural and leaves `trial_index` intact; M2 owns the FULL
`trial_index` -> `index_step_`* swap (write + subgroup removal) in one commit.
- `index_step_*` columns = code-derived (D4); no invented fields.
- `daq_collection_type` AND `protocol_sampling_mode` both kept (D5).
- `99_Unrouted` under `metadata` (D7); root holds only `schema_version` + two groups.
- Field renames (A-E) stay OUT of this migration (D9); D9 executes before M3.
- `index_cycle_*` / `protocol_*` per-cycle block = single_finite only (D10); does not
touch `index_step_*` (D4 intact).
- Raw archive file = raw-voltage `aidata` + embedded calibration VALUES + all
parameters; no calibrated output from exporter (D11; implementation deferred M2b).
- All MOI fields under `calibration_inertia_*`; "inertia" (not "inertial") throughout
(D12).

---

## Deferred

Items locked in principle but without an assigned implementation session. Do not
start M3 until flag 1 is resolved. Flag 2 is closed (resolved into M1).


| Flag | Description                                                                                                                                                                                                                                                                                                                    | Unblocks |
| ---- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | -------- |
| 1    | **D10/D11 implementation steps.** M2b (D11: embed calibration values + remove calibrated `forcetorque_raw`) and M2c (D10: write `index_cycle_*`/`protocol_*` per-cycle block) each need their own commit and session. Scope and session not yet assigned.                                                                      | M3       |
| 2    | **D9 / inertia rename conflict — RESOLVED into M1.** `calibration_inertia_*` normalization (D12) rides with M1's structural flattening of `inertial_calibration_profile`; the "inertial" -> "inertia" rename applies to all `inertial_*` keys M1 touches. D9 Flags A-E remain separate. No M0 decision needed.                 | CLOSED   |
| 3    | **Apparatus-MOI lookup key.** `calibration_inertia_matrix` lookup must be keyed by (clamp span, distance-from-AOR). Confirm that a distance-from-AOR measurement field is written per run before M2b is implemented.                                                                                                           | M2b      |
| 4    | **Apparatus-inertia calibration workflow.** Running empty-apparatus trials, deriving MOI from measured accelerations, and populating `calibration_inertia_matrix` is a separate acquisition + analysis procedure. `calibration_inertia_matrix` ships as an empty field until this workflow is built. Requires a separate plan. | post-M4  |


---

## Cutover record

- Pre-migration HEAD at planning: `4503a06` (segmented-finite Step 5 + checkbox no-op).
- Migration cutover commit (fill in after M1 lands): `_________`_.
- Files written with `01_Metadata`/`02_TimeSeries` exist only at or before the cutover
commit's parent; read them by checking out that commit.

