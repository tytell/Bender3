# CritterGripper / Bender — GUI UX Specification

## Authority metadata
- Purpose: Defines the user-experience contract for the Streamlit operator console.
- Authority level: Tier 2 (binding for all UI tasks).
- Scope: `bender_streamlit_gui.py` and any file rendering UI.
- Read by: Cursor on every session via `.cursor/specs/`.
- Last reviewed: 2026-06-03.

This document states the intended user experience as a contract. When fixing
bugs or adding features, conform to this spec. If a request conflicts with
this spec, stop and ask before proceeding.

---

## 1. Core design philosophy

The app is a **lab instrument control panel**, not a web form or dashboard.
Scientists work through a fixed sequence to configure and run a bending
experiment. The interface should feel like a physical instrument panel:
deliberate, predictable, and safe.

- Values are **staged** by the user, then **committed** with an explicit Apply.
- Nothing the user types affects the experiment until they Apply it.
- The page is a **single scrollable sequence**, not tabs. All sections always
  rendered, so all fields always live in session_state.

---

## 2. Layout

### Home (launcher)
The app opens to a launcher with two primary entry points:
- **Run Experiment** → the single-page experiment workflow (Section 2.1).
- **Review Data** → the H5 explorer for all post-hoc review: open recorded
  `.h5` files, plot series, and append notes. This is the single home for
  reviewing data, whether from the run you just finished or from older files.

An offline **Simulate Bending Mechanics** tool may remain accessible from
home for hardware-free exploration.

### 2.1 Run Experiment (single scrolling page)
The page is organized as **four vertical sections**, each with exactly **one
Apply**. Sections stack vertically (never side-by-side); columns are allowed
only for fields *within* a single section. Sections are ordered by how often
their values change between runs (most stable first):

1. **Config / Filepath** — hardware configuration + data file path. A single
   Apply ("Apply setup") commits both. This section carries the **run gate**:
   the experiment sections below stay hidden until Config + Path are applied.
   The gate is a safety feature and must be preserved.
2. **Specimen** — identity (genus/species, specimen ID, prep/segment label),
   universal morphometrics (lengths, mass), session temperature (editable),
   and prep condition. One Apply commits only this section's fields. There is
   a single "prep condition" field — no duplicate.
3. **Clamp geometry + inertial correction** — clamp spacing/offsets,
   cross-section, mounted body profile, density, and the inertial-correction
   flag. One Apply commits only this section's fields (the former separate
   clamp-geometry and profile/inertial Applies are merged into one).
4. **Experimental protocol** — protocol/run fields, preview, run, save.

- Each Apply commits **only** its own section's fields.
- No stepwise tabs. All sections always rendered, so all fields always live
  in session_state.
- Review of recorded data does **not** live here. All post-hoc review,
  plotting, and note-appending belong to the Review Data module. Run
  Experiment ends at Preview/Run; reviewing what you collected happens in
  Review Data.
- Template loading happens inline within the relevant sections — it is not a
  separate "mode." There is one experiment engine; loading a template is just
  populating fields on that engine.
- Sidebar carries navigation icons (one per section) for quick scroll-to.

### 2.2 Deferred — Training modes (NSF Broader Impacts)
The earlier "Build from scratch / Templates / Step-by-step" entry points are
retired as separate engines. Their pedagogical intent is preserved as a
**deferred** feature: a training overlay built on top of the stable single-page
engine, added only after the core engine is stable. Planned overlays:
- Guided/step-by-step presentation (one section at a time, progress, hints).
- Template-first flow for returning users.
Build order: stabilize the engine first, then layer training UX on top and
document it for the grant.

---

## 3. The Apply-commit model (non-negotiable)

This is the central interaction contract.

- The user fills fields freely. Typing has **no downstream consequence**.
- Each of the four sections (2.1) has **exactly one** explicit **Apply**
  button, and that Apply commits only that section's fields.
- Values are committed to `st.session_state` and validated **only** when Apply
  is clicked.
- Status icons (green check / yellow caution) update **only** after Apply.
- Downstream actions (hardware config load, protocol build, run) trigger
  **only** after the relevant Apply.

### Dirty state
- After Apply, a section is "clean" until the user edits a field again.
- Editing a field after Apply marks the section "dirty" (yellow icon returns).
- Re-applying clears dirty state.

---

## 4. Session state persistence rules

What must survive each event:

| Event | What survives | What resets |
|-------|--------------|-------------|
| Scrolling between sections | Everything | Nothing |
| Loading a template | Template-populated fields update; rest survives | Nothing |
| Snapshot load | All scientific data fields | UI preferences |
| Workflow reset to home | UI preferences only | Scientific data fields |

- Tab/section navigation never clears scientific data.
- UI preferences are not restored by snapshot or autosave — they belong to
  the user's current session.

### Theme and accessibility removal
Color/theme adjustability and the large-text toggle are **removed**. They
added UX overhead and caused repeated widget-write crashes for no scientific
benefit. The app ships with a single fixed theme. Users who need larger text
use the browser zoom (Cmd/Ctrl +). No theme or text-size widgets remain in
session_state.

---

## 5. The Streamlit widget-write constraint (hard rule)

Streamlit forbids writing to a widget's key in `st.session_state` **after**
that widget is instantiated. This single rule has caused the recurring
`StreamlitAPIException` / `StreamlitValueAssignmentNotAllowedError` crashes.
The rules below are mandatory and must be followed without exception.

### 5.1 Initialization — DO
Initialize every non-button widget key **once**, before the widget renders,
with a guard. Never use an unconditional assignment.
```python
# Correct — guarded, runs once
if "gui_specimen_id" not in st.session_state:
    st.session_state["gui_specimen_id"] = ""
```
```python
# WRONG — unconditional, overwrites user input every rerun
st.session_state["gui_specimen_id"] = ""
```

### 5.2 Buttons — DO NOT
Never initialize or write a button key in session_state. Buttons are
momentary and self-managed. Any line that assigns to a `gui_btn_*` or
`*_btn_*` key must be removed.

### 5.3 Validation — DO
Read the widget's return value directly, not session_state, to avoid the
one-rerun lag on status icons.
```python
# Correct
specimen_id = st.text_input("Specimen ID", key="gui_specimen_id")
if specimen_id:
    show_green_icon()
```

### 5.4 Restore operations — REQUIRED exclusions
Any restore (`snapshot`, `autosave`) MUST filter keys through
`_is_restore_safe_key()` and exclude all of the following. Do not write any
of these during restore:
- All button keys: `gui_btn_*`, `*_btn_*`
- Any widget key that has already rendered earlier in the same script run
- Any leftover UI preference keys from removed theme/text-size widgets
  (these should no longer exist, but exclude defensively if encountered)

### 5.5 When a widget-write error appears — fix procedure
1. Read the traceback to find the exact key being written (e.g.
   `gui_ui_large_text`).
2. Locate the write. If it is a restore operation, add the key to the
   `_is_restore_safe_key()` exclusion list. If it is a button key, remove
   the write entirely. If it is an init, convert it to the guarded form in 5.1.
3. Do not move the widget or change its render order to "make room" for the
   write — exclude or guard the write instead.
4. Re-run and confirm the crash is gone before moving on.

---

## 6. Templates

- The user can load templates for: hardware config, specimen identity,
  and protocols.
- Loading a template populates the relevant fields, then the user reviews
  and clicks Apply to commit.
- Template loading must work identically in the single-page layout.
- Saving a template captures the current applied values.

---

## 7. File and folder selection

- All file/folder selection uses **pasted text input**, not dropdowns.
- No filesystem-browsing dropdowns (they cause rerun conflicts).
- Folder paths validate with `os.path.isdir()`; file paths with
  `os.path.isfile()`, after the input, with clear success/error feedback.

---

## 8. Snapshots

- Save snapshot writes session state to JSON in `SessionSnapshots/`.
- Load snapshot restores scientific data fields, excluding the keys named
  in Section 5.
- Snapshot restore reuses the existing autosave restore pipeline and its
  safe-key filtering.

---

## 9. Safety (instrument control)

- Never trigger hardware-affecting actions on keystroke — only on Apply / Run.
- Warn before run if force/torque calibration is missing.
- Hardware import failures (`nidaqmx`) on a non-rig machine are expected and
  must not be treated as bugs.
- The app must remain fully GUI-testable on a Mac without hardware.

---

## 10. Known deferred items (not bugs)

- Bilateral pre/post stimulus in `make_stimuli()` — approved design, validate
  on the rig.
- `bio_*` / "biometrics" key rename to "morphometrics" — deferred refactor.
- `prepoststim_sep` wiring — part of the bilateral stim extension.

---

*bender_ux_spec_v1_0 | Jimenez Lab | 2026-05-28*
