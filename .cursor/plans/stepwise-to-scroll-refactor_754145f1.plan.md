---
name: stepwise-to-scroll-refactor
overview: Restructure the stepwise workflow rendering in `bender_streamlit_gui.py` so all experiment sections render sequentially on a single scrollable page while preserving all template-loading behavior, Apply-button commit semantics, and existing session_state keys.
todos:
  - id: map-step-gates
    content: Isolate all stepwise-only visibility/chrome gates and keep template gates untouched.
    status: completed
  - id: flatten-render-order
    content: Update workflow render flow to always show sections sequentially in target order for the single-page route.
    status: completed
  - id: compat-preservation
    content: Preserve all Apply handlers, template loaders, and session_state keys unchanged.
    status: completed
  - id: regression-checks
    content: Run targeted UI regression checks for template permutations, apply workflows, preview/run, and review sections.
    status: completed
isProject: false
---

# Single-Page Scroll Refactor Plan

## What currently drives stepwise rendering
- Route gate: `main()` branches on `_nav_route()`; `stepwise` path uniquely renders `_render_stepwise_rail()` and enables step-gated section helpers.
- Step index state: `gui_stepwise_step` (via `_stepwise_step()`) controls visibility.
- Section gates:
  - `_show_hw_config_section()` and `_show_data_path_section()` only show at step `0` in stepwise.
  - `_show_full_sec2()` (Measurements block) only shows at step `1`.
  - `_show_sec3_through_6()` (Experiment type/parameters + preview + run + save) only shows at step `2`.
  - `_show_sec7_and_8()` (Visualize + Add note) only shows at step `>=3`.
- Stepwise-only chrome: `_render_app_chrome()` injects `.bnd-stepwise-active`, compacts layout CSS, and shows stepwise Home UI.

## Minimum change to render all sections sequentially (no tabs/step rail)
- Keep section bodies and Apply handlers as-is; only remove stepwise gating around visibility/chrome.
- In `main()`:
  - Stop calling `_render_stepwise_rail()` for the workflow page.
  - Force sequential rendering by making all section gates evaluate true for the workflow route you use for this page.
- Adjust section visibility helpers to preserve template-only logic while removing `gui_stepwise_step` dependence:
  - Keep `_tpl_only_procedure()` checks unchanged.
  - For the single-page route, return `True` for: hardware, data path, measurements, protocol/run/preview/save, visualize/review.
- Keep route-level non-workflow pages unchanged (`landing`, `sim_compare`, `h5_explorer`).
- Keep ordering explicit in the render flow to match your target sequence:
  1) Hardware Configuration
  2) Data File Path
  3) Specimen Identity
  4) Morphometrics & Conditions
  5) Clamp Geometry
  6) Mounted Body Profile (inertial model)
  7) Protocol / Run
  8) Experiment Preview
  9) Review Data

## Exact behaviors to preserve
- Template loading behavior and order (must remain byte-for-byte equivalent in control flow):
  - Hardware config load/build actions in section 1.
  - Biometrics template load/save block in section 3.
  - Protocol template load/save in section 4.
- Apply semantics:
  - Setup `Apply setup` still commits config/path only when clicked.
  - Biometrics/form submits (`Apply specimen identity`, `Apply section`, `Apply all biometrics`) remain unchanged.
  - Procedure `Apply procedure` / `Refresh experiment preview` remain unchanged.
- All existing session keys, including stepwise keys (keep for backward compatibility even if unused in new layout).
- Existing helper side effects and dirty-baseline tracking (`_touch_*_baseline_if_clean`, `_soft_apply_reminder`, confirmation flags).

## Hidden/conditional rendering that needs special handling
- Template mode checkbox gates (`gui_tpl_chk_config`, `gui_tpl_chk_biometrics`, `gui_tpl_have_protocol_template`) currently hide large chunks; this logic must remain.
- Section hide toggles (`gui_bio_hide`, `gui_exp_hide`, `gui_sec6_hide`, `gui_sec7_hide`, `gui_sec8_hide`) should remain functional in scroll mode.
- Stepwise-only nuance: `_stepwise_on_data_file_path_step()` affects setup-button help/toast phrasing; this should be neutralized safely so behavior is deterministic in single-page mode.
- Existing stop condition (`if 'bender' not in st.session_state: st.stop()`) must still trigger exactly as now when setup hasn’t instantiated/loaded hardware.

## Risks and mitigations
- Risk: accidental behavior drift in template workflow due to over-broad gate edits.
  - Mitigation: only remove stepwise-step checks; do not alter template predicates or load/apply functions.
- Risk: stale CSS/chrome assumptions from `.bnd-stepwise-active` affecting spacing.
  - Mitigation: ensure stepwise CSS injection path is not active for single-page workflow; verify standard workflow chrome remains.
- Risk: subtle Apply-path differences from `_sw_dp` branching.
  - Mitigation: keep `_apply_setup_action()` intact and set deterministic non-stepwise path for `_sw_dp`.
- Risk: hidden reliance on `gui_stepwise_step` in diagnostics logs/sidebar.
  - Mitigation: retain key + helpers, but treat as inert; verify logs/checklist still render without route-step coupling.

## Verification checklist (post-implementation)
- Single-page workflow displays all 9 target sections in order with no step rail.
- Template mode permutations still hide/show the same sections as before.
- Apply buttons in setup, biometrics, and procedure still gate writes exactly on click.
- Protocol preview and Run behavior unchanged.
- Review/visualization and note append still operate on the selected file.
- No session key removals or renames; legacy sessions with `gui_stepwise_step` load without errors.