# CritterGripper EVOLUTION_PLAN

## Authority metadata

- Purpose: Prioritized roadmap/checklist for hardening and release readiness.
- Authority level: Tier 3 (hybrid).
- Scope: Planning and execution sequencing.
- Owner: Project maintainer.
- Last reviewed: 2026-04-23.

## How to use this file

- This file is advisory by default.
- It becomes binding only when a task/user explicitly references specific items or sprints.
- It does not override Tier 1 (`.cursorrules`) or Tier 2 (`.cursor/specs/ui_style.md`) constraints.
- Use the canonical user-facing triad consistently: **Hardware configuration**, **Morphometrics**, **Protocol**.

This document proposes a structured set of 100 critical “finalization” questions across:
- Hardware Safety
- UI/UX
- Data Integrity
- Edge Cases

Each question is formatted as a checkbox and ranked from `1` (Critical/Safety) to `100` (Polishing/Aesthetic). Questions are then grouped into 10 sprints of 10 questions each.

## Sprint 1 Pre-Analysis (first 10 questions)
Based on a code read-through of the current UI + DAQ + export pipeline, the files most likely involved in addressing Sprint 1 are:
- `bender_streamlit_gui.py`: integrate/upgrade the “Abort” flow, ensure “Stop DAQ & reset NI device” prevents any downstream export/QC, block “Run experiment”/save after stop/timeout until re-acknowledged, and gate export/QC so it never runs with stale/partial buffers.
- `bender_functions.py`: strengthen cancellation/safe-shutdown semantics around `Bender.run` (task ordering), ensure timeout/exception ends in a known safe electrical state without downstream export, and add explicit flags/hooks for GUI-driven interruption.
- `bender_daq_kill.py`: confirm `daq_emergency_stop` is sufficient for immediate stimulation/motor safety and consider enhancements to verify “tasks stopped + outputs cleared” with reliable operator messaging.
- `bender_h5_export.py`: add/enforce checks so export is refused (or clearly marked invalid) when acquisition did not complete successfully and `trial_records` is empty/partial.

## Sprint 1: Preventing Hardware Damage (1-10)
Goal: Make emergency stop, abort/cancel, timeouts, and motion/stimulation bounds safe and deterministic so hardware never transitions into an unsafe or confusing state.

- [ ] 1. Does the emergency stop (sidebar “Stop DAQ & reset NI device”) immediately halt all NI tasks and clear outputs (analog/digital) even if acquisition is currently running?
- [ ] 2. If the user clicks “Abort” during a pending-confirmation warning, does it fully prevent DAQ from starting and never transition into saving/export steps?
- [ ] 3. If a run times out or throws inside `Bender.run`, does the code always stop all NI tasks and leave stimulation/motor outputs in a known safe state (for example, 0V and motor disabled)?
- [ ] 4. During `wait_until_done` timeouts, does the system detect the timeout reason and show a clear operator message without attempting export/QC using stale buffers?
- [ ] 5. Are there software-enforced motion/stimulus bounds (max velocity, sample rates > 0, stim pulse rate/duty, voltage limits) validated before hardware tasks start?
- [ ] 6. Is every “motion is too fast / invalid” constraint guaranteed to run before any AO/DO tasks are started?
- [ ] 7. Does `daq_emergency_stop` reset the correct device(s) for the currently configured `device_name`, and if `device_name` is missing, does it fail safely with no ambiguity about what was reset?
- [ ] 8. After an emergency stop, does the GUI block further “Proceed/Run experiment” actions until the operator re-applies setup and acknowledges the stop?
- [ ] 9. If the NI device is reserved/busy, does the app prevent partial task creation and clearly differentiate “DAQ busy” vs “DAQ config invalid”?
- [ ] 10. Is there protection against overlapping runs (double-clicks, refresh/re-entry, re-entrant UI code) so new `Task()` instances never start while a previous run is still active?

## Sprint 2: Making Operator Workflow Unambiguous (11-20)
Goal: Ensure the operator always knows what is selected vs applied, what will happen next, and that UI controls cannot trigger unsafe or unintended transitions.

- [ ] 11. Does the “Run experiment” button become disabled while acquisition is in progress to prevent concurrent `b.run_experiment` executions?
- [ ] 12. When users change Step 2 selections (folder/filename), does the UI explicitly indicate whether that choice is “selected” vs “applied/committed” to the experiment object?
- [ ] 13. If users change folder/filename after applying setup, does the app automatically mark setup as unconfirmed (or warn) before allowing run/save?
- [ ] 14. Are browse/selection failures (tk dialogs) surfaced with actionable guidance and a fallback manual entry option?
- [ ] 15. Do stepwise navigation controls prevent skipping required safety-critical selections (hardware config, destination, morphometrics basics) even when using “Back/Next”?
- [ ] 16. When soft warnings are generated (missing calibration, incomplete morphometrics, destination incomplete), do “Proceed/Abort” messages clearly state consequences for hardware and saving?
- [ ] 17. Is there a clear “what will be written where” summary immediately before the operator proceeds to DAQ acquisition?
- [ ] 18. If “Stop DAQ & reset NI device” is pressed during/around a run, does the UI prevent further export/QC steps until the operator restarts the workflow?
- [ ] 19. Are UI states consistent across full workflow and template mode so the operator cannot hit a hidden “different semantics” trap?
- [ ] 20. Are help texts/tooltips accurate for every critical widget that affects hardware motion/stimulation and saving destinations?

## Sprint 3: Ensuring Data Files Are Trustworthy (21-30)
Goal: Prevent corrupted or misleading outputs; ensure exports and QC are only written when acquisition has produced coherent trial data.

- [ ] 21. Does `export_primary_h5` always avoid overwriting prior experimental data (unique naming) unless explicitly configured otherwise?
- [ ] 22. Does export detect empty/partial `trial_records` (for example, aborted/failed acquisition) and refuse to write or clearly mark file invalid?
- [ ] 23. Are array lengths in the HDF5 datasets consistent and trimmed/padded in a deterministic way to avoid corrupted time series?
- [ ] 24. Are calibration-related metadata fields (for example, `use_inertial_calibration`, `calibration_file`) stored correctly and validated for existence?
- [ ] 25. Are schema version and `test_type` recorded consistently across both root attributes and group metadata?
- [ ] 26. Is post-trial notes appending safe and idempotent so repeated exports do not duplicate the same note block?
- [ ] 27. Do QC plots use a deterministic naming scheme tied to the chosen `.h5` and run metadata, so operators can always match QC to data?
- [ ] 28. If QC save fails, does the GUI report “acquisition succeeded but QC failed” (not “all succeeded”)?
- [ ] 29. Do save operations validate that target directories exist and are writable, with clear error messaging and no silent partial writes?
- [ ] 30. When HDF5 export throws mid-write, is there a strategy to prevent leaving behind misleading partial files (rename/cleanup/atomicity)?

## Sprint 4: Hardening Path & Session Edge Cases (31-40)
Goal: Make path composition, normalization, and session restore robust so the wrong directory never gets used.

- [ ] 31. Does the app robustly handle Windows long paths, special characters, and UNC/network paths during output composition?
- [ ] 32. When users select a folder but leave filename blank, does the composed output resolve to empty and block run/save rather than using stale `bender.outputfile`?
- [ ] 33. If a user pastes a full path into the filename field, does the app normalize it consistently to basename to avoid mixing folder+filename incorrectly?
- [ ] 34. Are path comparisons correct on Windows (case-insensitive, normalization) so setup confirmation does not flicker unexpectedly?
- [ ] 35. If the selected directory is deleted or permissions change after selection but before run, does the app preflight and block instead of failing mid-acquisition?
- [ ] 36. Do session restore/autosave operations never restore stale `bender.outputfile` that no longer matches the currently selected folder/filename?
- [ ] 37. Are multiple tabs/sessions isolated so one operator cannot affect another operator’s session_state through shared Streamlit server state?
- [ ] 38. Do Streamlit reruns (widget interactions) mid-flow preserve invariants so in-memory experiment state never diverges from UI?
- [ ] 39. When users switch hardware config modules, does the app correctly invalidate dependent state so the next run does not reuse mismatched calibration channels?
- [ ] 40. Are missing morphometrics/template files handled gracefully so the operator cannot apply partial measurements and proceed unknowingly?

## Sprint 5: Reliability for Templates & Procedure Loading (41-50)
Goal: Ensure template and procedure loading are deterministic, validated, and cannot accidentally mismatch data destinations or experiment parameters.

- [ ] 41. When loading a hardware config module, does the app clear stale experiment buffers so new runs never reuse old data arrays?
- [ ] 42. When applying setup (data path), does the dirty/applied state tracking reliably reflect whether the experiment object actually changed?
- [ ] 43. When applying morphometrics templates, does the app always apply all required blocks (intrinsic, clamp geometry, mounted/inertial profile) and mark measurement state correctly?
- [ ] 44. Are there guardrails preventing morphometrics application without a valid data folder when templates assume those paths?
- [ ] 45. Do template workflow loaders correctly handle “reload module” without losing user-selected data destination selections?
- [ ] 46. If JSON parsing fails for morphometrics/template files, does the error message name the offending file/key and provide next steps?
- [ ] 47. After morphometrics application, are derived flags (for example, inertial usage and computed derived parameters) updated consistently in the `Bender` object?
- [ ] 48. After template loading, does the app enforce that required measurement values exist before enabling procedure/run actions?
- [ ] 49. Does protocol preview generation match exactly what DAQ acquisition will use (no parameter drift between preview and run)?
- [ ] 50. If protocol parameters change after refresh preview, does the UI force the operator to refresh preview before allowing run?

## Sprint 6: Safety Guards for Protocol Parameters (51-60)
Goal: Prevent invalid protocol parameters from generating hardware-level faults, NaNs, or unsafe motor/stimulus outputs.

- [ ] 51. Is there an explicit upper bound (or enforced constraint) on stimulation amplitude/voltages, duty cycle, and pulse timing before DAQ starts?
- [ ] 52. Are computed duty-cycle/pulse-width parameters validated to prevent overflow, negative timing, or NaN propagation into DAQ scheduling?
- [ ] 53. Are timelines validated for monotonicity and finite spacing (for example, `t[1] - t[0] > 0` and no NaN/Inf)?
- [ ] 54. In unit conversions (angle/curvature/strain), do conversion functions guard against division by zero and invalid geometry inputs?
- [ ] 55. Does the motor pulse generator enforce “Motion is too fast!” before any AO/DO tasks begin output scheduling?
- [ ] 56. For isometric/isovelocity protocols, are produced timelines validated so they are not empty, one-sample, or inconsistent in length?
- [ ] 57. Are the “stim cycles per step” and index-range assumptions validated so invalid user settings cannot create out-of-range indexing?
- [ ] 58. If the operator selects an unsupported test type, does the app fail fast before any hardware communication?
- [ ] 59. Are config-provided sample rates sanitized and validated across AI and AO/DO channels (finite, > 0, consistent)?
- [ ] 60. Does the app provide a sufficient “preflight preview” that operators can inspect for motor and stimulation characteristics before enabling hardware?

## Sprint 7: Robust Error Reporting & Recovery (61-70)
Goal: Make failures understandable and recoverable so operators can return to a safe state quickly.

- [ ] 61. When DAQ fails, does the friendly error path include the exception type plus key context (device name, sample rates, outputfile)?
- [ ] 62. Are error messages consistent and actionable across all GUI entry points (Streamlit vs other interfaces)?
- [ ] 63. Does the app log enough metadata so failures can be reproduced later (timestamp, test type, configuration module, data destination)?
- [ ] 64. If an error occurs after acquisition but before export, does the app clearly label which phase failed?
- [ ] 65. If QC fails due to missing/invalid HDF5 content, does the GUI report that clearly and point to likely causes?
- [ ] 66. Does the app prevent Streamlit crashes on file I/O exceptions (permissions, missing directories, invalid HDF5 reads)?
- [ ] 67. Are advanced diagnostics features safe and do they avoid leaking secrets or irrelevant internal paths?
- [ ] 68. For template/morphometrics loaders, does the UI indicate schema/version compatibility to prevent silent misinterpretation?
- [ ] 69. Are common filesystem errors (permission denied, wrong directory, locked file) handled gracefully with retries or next steps?
- [ ] 70. Are there recovery pathways to retry failed operations without needing a full browser refresh?

## Sprint 8: Responsiveness & Scalability (71-80)
Goal: Keep the UI responsive and prevent resource spikes from extremely large inputs.

- [ ] 71. Does the app avoid expensive recomputation on every rerun by caching preview computations where feasible?
- [ ] 72. Are large settings/diagnostic tables rendered in a way that doesn’t freeze the UI?
- [ ] 73. Are charts truncated to a safe max point count (for example, `pv_pts`) so Plotly rendering remains responsive?
- [ ] 74. Do file discovery functions (modules/templates) run efficiently and avoid repeated costly filesystem walks?
- [ ] 75. Does QC figure generation reuse cached intermediate results when appropriate?
- [ ] 76. Are there guardrails preventing extremely large user-provided timelines that could blow memory (duration/frequency lists)?
- [ ] 77. Do long operations show progress indicators consistently (spinners/status updates) and do they not time out silently?
- [ ] 78. When reading HDF5 for preview/explorer, does it avoid loading full datasets when only metadata is needed?
- [ ] 79. Are resources cleaned up after long operations so repeated runs do not leak memory/handles?
- [ ] 80. Are autosave/restore snapshot operations size-limited so session snapshots do not become excessive?

## Sprint 9: Quality Gates via Tests & CI (81-90)
Goal: Ensure changes do not regress safety/data behaviors and keep coverage for critical logic.

- [ ] 81. Do automated tests cover every motion test type including boundary constraints (sweep/step/dynamic/step_change/isometric/isovelocity/calibration)?
- [ ] 82. Do tests include failure injection for DAQ-missing and export failure modes (permission denied, invalid output dir)?
- [ ] 83. Are there tests for output path composition and normalization edge cases (filename-only vs folder+filename, full-path input pasted)?
- [ ] 84. Do tests validate QC naming/linkage to the exact `.h5` file saved?
- [ ] 85. Do tests verify that dirty/applied tracking behaves correctly across reruns, refresh preview, and apply setup?
- [ ] 86. Do tests ensure `unique_filepath` never overwrites existing data?
- [ ] 87. Do tests validate that `export_primary_h5` writes correct schema groups/datasets for representative `trial_records`?
- [ ] 88. Are tests CI-friendly so they run without NI-DAQmx installed (simulation mode only)?
- [ ] 89. Are style/lint checks included in CI to catch runtime hazards (missing imports, syntax errors, unused keys)?
- [ ] 90. Are there performance regression checks that detect significant slowdowns in preview/export logic?

## Sprint 10: UI Polish & Operator Confidence (91-100)
Goal: Improve clarity, consistency, and usability so operators trust the app and onboarding time drops.

- [ ] 91. Are step labels and navigation consistent across all routes (stepwise, scratch, full workflow)?
- [ ] 92. Is the “selected vs applied” UX consistent across the setup tab, sidebar status check, and any template workflow routes?
- [ ] 93. Are all critical buttons styled consistently (colors, emphasis, width) and do they match operator expectations?
- [ ] 94. Does CSS injection remain robust across Streamlit version changes (no broken selectors/layout drift)?
- [ ] 95. Are accessibility basics covered (contrast, focus order, readable font sizes, keyboard navigability)?
- [ ] 96. Is wording standardized for “Apply setup”, “Apply section”, “Apply morphometrics”, “Proceed”, and “Abort” to avoid ambiguity?
- [ ] 97. Are tooltips or “operator cheat sheets” available for the safety-critical steps so new operators can follow them?
- [ ] 98. Does the app run startup preflight checks (nidaqmx availability, device name present, output directories) and communicate results?
- [ ] 99. Are session restore and autosave behaviors transparent enough that operators understand what state was reused?
- [ ] 100. Does the UI provide a concise end-of-run summary (what ran, where it was saved, timestamps, QC outputs) in one place?


---

# Backlog item: JLAB_ROOT path anchor — 2026-05-28

**What:** Add a project-root anchor so the app resolves all file/folder paths relative to a single base, rather than requiring full absolute paths everywhere.

**Why:**
- Matches existing JLab convention from jlab_folder_architecture SOP ("Set JLAB_ROOT once per machine").
- Satisfies .cursorrules rule "Paths read from manifest, never hardcoded."
- Addresses GUI audit items 7, 59-63 (file selection ergonomics).
- Cross-machine portable — same relative path works on Mac dev and Windows rig.

**Approach (preferred): environment variable**
- Set once per machine: `export JLAB_ROOT="/Users/yjimenez/.../01_JimenezLab"` (Mac) or equivalent on Windows.
- App reads via `os.environ.get('JLAB_ROOT', <fallback>)`.
- All template / data / config paths in the UI become short relative paths.
- Full path resolved internally: `os.path.join(JLAB_ROOT, uselative_path)`.

**Touches:**
- bender_streamlit_gui.py — every file/folder text_input and selectbox.
- bender_functions.py — any path resolution that bypasses the UI.
- Template loaders for morphometrics, protocols, hardware config.

**Risk:**
- Touches every file path in the app. Requires the validation-gate discipline from jlab_project_consolidation_sop §5: run end-to-end against a known-good output before declaring done.
- Build as its own focused session, not bundled with other fixes.

**Status:** Logged. Not started.

---

# GUI Spec Audit — 2026-05-28

Master violation list of bender_streamlit_gui.py against ux_spec.md and ui_style.md.
**101 items total.** Full table lives in `.cursor/plans/gui_spec_audit_list_75d7dbea.plan.md`.

Disposition each item as: **FIXED** / **DEFERRED** (with reason) / **SPEC-AMENDED**. None dropped silently.

## Triage tiers

### Tier 1 — Crash or correctness risk (this session)
- Items **52–58** (widget-write rule): same crash class hit 3× during 2026-05-27 session.
- Items **9, 10, 25** (hardware/run actions firing without Apply): scientific-correctness risk.

### Tier 2 — Root of original GUI bugs (next session)
- Items **1–8, 11–17** (Apply-commit model violations): why status icons misfire and values feel unstable.
- Items **30–35** (missing guarded init): the field-clearing class.

### Tier 3 — Cleanup
- Items **36, 78** (theme/large-text widgets): remov
- Items **64–66, 69–73** (layout: dead stepwise route, duplicate section numbers, three-mode launcher).

### Tier 4 — Decide, do not drop
- Items **26–29** (`gui_` prefix violations on `fld_*` / `morpho_*`): **DEFERRED** for safety — renames must change reads, writes, and HDF5 saves together or data breaks. Revisit as dedicated session.
- Items **59–63** (filesystem selectboxes vs paste-only §7): **DECISION NEEDED** — either fix to match spec, or amend §7 to allow short known-list pickers. Likely resolves naturally with JLAB_ROOT path anchor (see backlog above).
- Items **84–101** (missing button `type=`): cosmetic. Fix during a styling pass.

**Status:** Tier 1 in progress. Tiers 2–4 pending.

## Conditional field validation
**Bug:** Stimulation fields (`Stim cycles per step`, etc.) are flagged as required even when `Enable stimulation` is unchecked. Same class as the `prep_condition` bug — validator does not respect conditional gates.
**Spec ref:** ux_spec §3.
**Action:** Audit validator for fields that should be required only when an enabling checkbox is true. Common pattern, likely affects other conditional groups (sono, stim monitor, profile inertial).

## Measurements green check fails after Apply
**Bug:** Measurements section status icon does not turn green even after fields are filled and "Apply section" / "Apply all morphometrics" is clicked.
**Spec ref:** ux_spec §3 (Dirty state).
**Audit cross-ref:** Related to items 11, 13, 49 (status icons reading unapplied state).
**Pre-existing:** Yes — present before Batch B fixes.
**Action:** Audit `_workflow_ready_state` / `_measurements_fields_ok` / `_morpho_apply_dirty` chain. Confirm Apply actually clears dirty state AND that the dirty check reads the post-Apply baseline correctly.
