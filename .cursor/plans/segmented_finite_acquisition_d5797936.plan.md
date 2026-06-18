# Segmented-finite acquisition: final implementation plan

Every decision is LOCKED. Each step is self-contained: a fresh session can execute
one step from this file alone (files touched, exact insertion points, the single
commit it produces, validation, recommended model, and a do-not-bundle boundary).

---

## Locked decisions (D1–D4 + reset/hold)

- **D1 — Naming.** `02_TimeSeries` subgroups are `step_001`, `step_002`, … :
  ONE-based, 3-digit, `step_` prefix. `step_index` is ONE-based too (PI updating
  the Tier-2 schema doc to match). Invariant: `step_NNN` suffix == `step_index`,
  both one-based. No zero-based references in new code. Filename convention
  separate and unchanged.

- **D2 — Recording buffer.** Hardcoded `step_buffer_s = 1.0` (no GUI field).
  Every segmented step's recorded window:
  `[1s neutral hold @0°, energized] → [move/ramp] → [active] → [return to neutral] → [1s neutral hold @0°, energized]`.
  Both bookends INSIDE the step subgroup, separate from and additional to
  `pre_baseline_s`/`hold_duration_s`/`post_baseline_s`/`pre_hold_s` (unchanged,
  at target/start). `rest_before_second` = UNRECORDED gap between step N
  post-buffer and step N+1 pre-buffer.

- **D3 — daq_collection_type.** Protocol-derived only, in `01_Metadata`, never in
  GUI: dynamic/sweep → `"continuous"`; isometric/isovelocity/FL/FV → `"segmented"`.
  Gate: do NOT label isometric `"segmented"` until Step 4 lands.

- **D4 — Reset floor.** Keep name `rest_between_steps_s`, operator-settable,
  hardcoded 2s minimum enforced at Run (not keystroke). Exact block message:
  `"Reset duration must be at least 2s to ensure clean recording buffers between
  steps. Physiological rest requirements will typically exceed this minimum."`

- **Reset/hold always-on.** `reset_between_steps` and `hold_motor_between_steps`
  forced True internally in the run path for segmented protocols; checkboxes become
  no-ops (removal is later cleanup). Steps 3/4 do not depend on checkbox state.

---

## Steps (each = one commit, one model, do-not-bundle)

### Step 1 — Mode maps + step_manifest scaffold + one-based step_index
**Recommended model:** Sonnet  
**Validation:** Mac only (no hardware needed)  
**Do not bundle with Step 2.**

**What it does:**
- Add `protocol_sampling_mode` derivation: `"single_finite"` for dynamic/sweep,
  `"segmented_finite"` for isometric/isovelocity/FL/FV.
- Add `daq_collection_type` derivation: dynamic/sweep → `"continuous"`;
  isovelocity/FV → `"segmented"`; isometric/FL → `"continuous"` (gate until Step 4).
- Write both fields to `01_Metadata` via `export_primary_h5`.
- Replace the hardcoded `protocol_acquisition_mode='continuous'` in
  `bender_routing_spec.py:242` with protocol-derived `daq_collection_type`.
- Write `step_manifest` JSON string to `01_Metadata` for segmented protocols
  (initially populated from existing trial_records; later enriched by Steps 3/4).
- Change `step_index` from zero-based to ONE-based in
  `bender_functions.py` (`_build_trial_record`, `_run_isovelocity_steps`, etc.)
  and `bender_h5_export.py`. Invariant: `step_NNN` suffix == `step_index`.
- Update test `tests/test_continuous_isometric.py:134` to expect new field name/value.

**Files touched:**
- `bender_functions.py` — step_index one-based in `_build_trial_record` (~:4034),
  isovelocity step loop (~:4729–4906)
- `bender_routing_spec.py` — replace hardcoded `protocol_acquisition_mode` (~:242)
- `bender_h5_export.py` — write `protocol_sampling_mode`, `daq_collection_type`,
  `step_manifest` to `01_Metadata`; step_index one-based in exporter loop (~:348)
- `tests/test_continuous_isometric.py` — update assertion at ~:134

**Commit message:** `feat: add protocol_sampling_mode, daq_collection_type, step_manifest scaffold; one-based step_index`

---

### Step 2 — export_primary_h5 branch: continuous flat vs segmented subgroups
**Recommended model:** Sonnet  
**Validation:** Mac only  
**Do not bundle with Step 1 or Step 3.**

**What it does:**
- Branch `export_primary_h5` on `daq_collection_type`:
  - `"continuous"` (dynamic/sweep): write channel arrays as flat datasets
    directly under `02_TimeSeries` (no `step_NNN` subgroup).
  - `"segmented"` (isovelocity/FV; isometric/FL still `"continuous"` via gate):
    keep existing subgroup structure (renamed `step_NNN` once Step 3 lands).
    For now, segmented path keeps `trial_{i:04d}` subgroups — renaming is
    coupled to Step 3 (isovelocity) and Step 4 (isometric).
- For the `"continuous"` branch: single flat write of all channels; no loop
  over trial records needed (only one entry for dynamic/sweep).
- Update export tests to verify flat layout for continuous protocols.

**Files touched:**
- `bender_h5_export.py` — branch logic in `export_primary_h5`

**Commit message:** `feat: branch export_primary_h5 on daq_collection_type; flat write for continuous protocols`

---

### Step 3 — isovelocity/FV: step_NNN subgroups + 1s bookends + per-step timing
**Recommended model:** Opus  
**Validation:** Mac (unit tests) then user-confirmed rig run  
**Do not bundle with Step 2 or Step 4.**

**What it does:**
- In `_run_isovelocity_steps` loop (~:4729), rename subgroup from `trial_{i:04d}`
  to `step_{step_index:03d}` (one-based).
- Add 1s neutral pre-hold (`step_buffer_s = 1.0`) at t=0 before the velocity
  ramp in each step timeline; reuse existing trailing neutral hold (`:4567–4578`)
  as post-buffer — do NOT stack a second post-buffer on top of it.
- Capture `wall_clock_start = time.time()` before each `run()` call (~:4830).
- Compute `rest_before_second` = wall clock delta from previous step's post-buffer
  end to this step's pre-buffer start (inter-`run()` gap).
- Capture `duration_second = len(self.t) / daq_ai_sample_rate_hz`.
- Populate `operating_point` (velocity in deg/s) and `operating_point_units`
  (`"degree_per_second"`).
- Enrich `step_manifest` with per-step timing from the captured values.
- Force `reset_between_steps = True` and `hold_motor_between_steps = True`
  in the run path (not via checkbox).
- Flip `daq_collection_type` for isovelocity/FV to `"segmented"` (already set in
  Step 1; no change needed here unless gate requires adjustment).

**Files touched:**
- `bender_functions.py` — `_run_isovelocity_steps`, `bender_gui_preview.py`
  (OFF run path — separate edit in same commit).
- `bender_h5_export.py` — subgroup naming, step_manifest population.

**Safety checklist:** motor stays energized, reset-in-finally correct,
return-to-home verified, no mid-step device reset added.

**Commit message:** `feat: isovelocity step_NNN subgroups + 1s bookends + per-step timing (segmented_finite)`

---

### Step 4 — isometric/FL: per-step run() loop + 1s bookends + step_NNN
**Recommended model:** Opus  
**Validation:** Mac (unit tests) then MANDATORY user-confirmed rig run + safety checklist  
**Do not bundle with Step 3 or Step 5.**

**What it does:**
- Restructure `_run_force_length_steps` (~:3881) from one stitched FINITE task
  to a per-step `run()` loop mirroring `_run_isovelocity_steps`.
- Add 1s neutral pre-hold (`step_buffer_s = 1.0`) as net-new neutral hold at t=0
  before the ramp (currently the timeline opens on the ramp at ~:3788).
- Remove in-waveform gap hold (~:3724) — rest gaps are now unrecorded between
  `run()` calls, not stitched into the waveform.
- Add 1s neutral post-hold after return-to-neutral.
- Capture `wall_clock_start`, `rest_before_second`, `duration_second`,
  `operating_point` (angle in degrees), `operating_point_units` (`"degree"`).
- Flip the D3 gate: change isometric/FL `daq_collection_type` to `"segmented"`.
- Rename subgroups from `trial_{i:04d}` to `step_{step_index:03d}`.
- Force `reset_between_steps = True` and `hold_motor_between_steps = True`.
- Port equivalent timeline changes to `bender_gui_preview.py` (OFF run path —
  separate edit in same commit).

**Files touched:**
- `bender_functions.py` — `_run_force_length_steps`, `_build_isometric_continuous_timeline`
- `bender_h5_export.py` — subgroup naming for isometric/FL
- `bender_gui_preview.py` — timeline preview (OFF run path)
- `tests/test_continuous_isometric.py` — update assertions for per-step structure

**Safety checklist:**
- Motor stays ENERGIZED across steps (no de-energize between).
- One `reset_device()` per `run()` call in `finally` (no mid-step resets).
- Return-to-home before and after each step.
- `hold_motor_between_steps=True` asserted before loop.
- Idle hold word reasserted after each reset.
- Start and end at home (angle 0).

**Commit message:** `feat: isometric per-step run() loop + 1s bookends + step_NNN (segmented_finite)`

---

### Step 5 — rest_between_steps_s 2s floor + checkbox no-ops
**Recommended model:** Opus  
**Validation:** Mac then rig confirm  
**Do not bundle with Step 4.**

**What it does:**
- Enforce 2s minimum for `rest_between_steps_s` at Run (not keystroke).
  Exact block message: `"Reset duration must be at least 2s to ensure clean
  recording buffers between steps. Physiological rest requirements will
  typically exceed this minimum."`
- Make `reset_between_steps` and `hold_motor_between_steps` checkboxes no-ops
  in `bender_streamlit_gui.py` (~:2495, ~:2515) — values are forced True in the
  run path. Checkboxes remain visible (removal is later cleanup).
- Optional backstop in `bender_functions.py`: assert both flags True at run entry
  for segmented protocols.

**Files touched:**
- `bender_streamlit_gui.py`
- `bender_functions.py` (optional backstop)

**Commit message:** `feat: enforce 2s rest_between_steps_s floor; reset/hold flags always-on for segmented`

---

## Execution order / session map

| Session | Steps | Model | Hardware | Commits |
|---------|-------|-------|----------|---------|
| A | 1 + 2 | Sonnet | Mac only | 2 |
| B | 3 | Opus | Mac then user-confirmed rig | 1 |
| C | 4 | Opus | Mac then MANDATORY rig + safety checklist | 1 |
| D | 5 + checkbox no-op | Opus | Mac then rig confirm | 1 |
| E | h5 migration (Phase 0) | Opus | Separate planning pass | TBD |

**Cross-session invariants:**
- One-based `step_index == step_NNN` throughout.
- `step_buffer_s = 1.0` bookends inside each segmented step, separate from at-target holds.
- `daq_collection_type` protocol-derived with the isometric gate until Step 4.
- `rest_between_steps_s` name kept; 2s Run floor from Step 5.
- `reset_between_steps` and `hold_motor_between_steps` forced True in run path
  (not via checkbox) from Step 3 onward.
