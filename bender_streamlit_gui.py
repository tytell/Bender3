"""
CritterGripper — Streamlit UI for Bender experiment dispatch.

Run from the project directory:
    streamlit run bender_streamlit_gui.py

Select ``test_type`` first; edit fields, use **Apply** to copy them onto the
``Bender`` instance, optionally **save / load protocol templates** (procedure only) and
**save / load biometrics templates** in **section 3**, **Check required fields** to validate,
**Refresh preview** for a table/plot of commanded motion (no DAQ), then **Run experiment**
when hardware is ready.
"""
from __future__ import annotations

import base64
import importlib
import json
import math
import os
import sys
from typing import Optional

import h5py
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

# Project root on path for `bender_functions` and config modules
_ROOT = os.path.dirname(os.path.abspath(__file__))
_LOGO_PATH = os.path.join(_ROOT, 'assets', 'jimenez_biomechanics_logo.png')
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)


def _nsf_logo_path() -> Optional[str]:
    for name in ('nsf_logo.png', 'nsf_logo.svg'):
        p = os.path.join(_ROOT, 'assets', name)
        if os.path.isfile(p):
            return p
    return None


def _img_data_uri(path: str) -> str:
    ext = os.path.splitext(path)[1].lower()
    with open(path, 'rb') as f:
        raw = f.read()
    b64 = base64.b64encode(raw).decode('ascii')
    mime = 'image/svg+xml' if ext == '.svg' else 'image/png'
    return f'data:{mime};base64,{b64}'


from bender_daq_kill import daq_emergency_stop  # noqa: E402
from bender_config_builder import (  # noqa: E402
    discover_config_modules,
    effective_load_module_name,
    parse_comma_list,
    parse_n_floats,
    read_base_defaults,
    render_generated_config,
    sanitize_config_module_stem,
)
from bender_functions import Bender  # noqa: E402
from bender_gui_preview import build_protocol_preview  # noqa: E402
from bender_biometrics_templates import (  # noqa: E402
    apply_biometrics_template_to_session,
    biometrics_template_display_label,
    default_biometrics_templates_dir,
    list_biometrics_template_files,
    load_biometrics_template,
    save_biometrics_template,
    sanitize_biometrics_filename_stem,
)
from bender_protocol_templates import (  # noqa: E402
    apply_template_to_session_state,
    build_procedure_dict_from_updates,
    default_templates_dir,
    list_template_files,
    load_protocol_template,
    save_protocol_template,
    sanitize_template_filename_stem,
    snapshot_bender_procedure,
    template_display_label,
)

_BND_LS_ACTION_MARK = '<div class="bnd-ls-action" aria-hidden="true"></div>'


def _inject_load_save_button_theme() -> None:
    """Inject CSS for persistent action buttons.

    Visual tiers (main panel):
    - **Load/Save (red, full width):** `_load_save_button` — commits that read/write disk or
      apply hardware + paths (e.g. Load hardware configuration and data path, Write config file and load).
    - **Mode / branch (default Streamlit primary + secondary, two columns):**
      `_hardware_configuration_mode_toggle` — mutually exclusive setup paths only; not file I/O.
    """
    st.markdown(
        """
<style>
section[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action) button[kind="primary"],
section[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action) button[kind="secondary"] {
    background-color: #b91c1c !important;
    background-image: none !important;
    color: #ffffff !important;
    border: 1px solid #7f1d1d !important;
}
section[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action) button[kind="primary"]:hover,
section[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action) button[kind="secondary"]:hover {
    background-color: #991b1b !important;
    border-color: #450a0a !important;
    color: #ffffff !important;
}
section[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action) button[kind="primary"]:focus-visible,
section[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action) button[kind="secondary"]:focus-visible {
    box-shadow: 0 0 0 0.2rem rgba(185, 28, 28, 0.45) !important;
}
/* Bordered layout panels (Streamlit container border=True) — stronger frame in main panel */
section[data-testid="stMain"] [data-testid="stVerticalBlockBorderWrapper"] {
    border-width: 2px !important;
    border-color: #94a3b8 !important;
    border-radius: 10px !important;
    padding: 0.65rem 0.85rem 0.85rem 0.85rem !important;
    background: #f8fafc !important;
    margin-bottom: 0.75rem !important;
}
/* Alert blocks: clearer frame */
section[data-testid="stMain"] div[data-testid="stAlert"] {
    border: 1px solid #cbd5e1 !important;
    border-radius: 8px !important;
}
</style>
""",
        unsafe_allow_html=True,
    )


def _inject_stepwise_compact_layout_css() -> None:
    """Tighten main padding and step rail when ``.bnd-stepwise-active`` is present (see chrome)."""
    st.markdown(
        """
<style>
body:has(.bnd-stepwise-active) .main .block-container {
    padding-top: 0.35rem !important;
    padding-bottom: 0.65rem !important;
}
body:has(.bnd-stepwise-active) section[data-testid="stMain"] [data-testid="stVerticalBlockBorderWrapper"] {
    padding: 0.4rem 0.55rem 0.5rem 0.55rem !important;
    margin-bottom: 0.35rem !important;
}
body:has(.bnd-stepwise-active) section[data-testid="stMain"] [data-testid="stProgress"] {
    margin-top: 0 !important;
    margin-bottom: 0.3rem !important;
}
body:has(.bnd-stepwise-active) section[data-testid="stMain"] [data-testid="stImage"] {
    margin-top: 0 !important;
    margin-bottom: 0 !important;
}
body:has(.bnd-stepwise-active) section[data-testid="stMain"] h1,
body:has(.bnd-stepwise-active) section[data-testid="stMain"] h2,
body:has(.bnd-stepwise-active) section[data-testid="stMain"] h3 {
    margin-top: 0.4rem !important;
    margin-bottom: 0.25rem !important;
}
body:has(.bnd-stepwise-active) section[data-testid="stMain"] [data-testid="stHeader"] {
    margin-top: 0.25rem !important;
}
</style>
""",
        unsafe_allow_html=True,
    )


def _load_save_button(
    label: str, *, key: str, help: Optional[str] = None, button_type: str = 'primary', **kwargs
) -> bool:
    """Full-width Load/Save-style action; pair with `_inject_load_save_button_theme` for red styling."""
    kwargs.pop('use_container_width', None)
    with st.container():
        st.markdown(_BND_LS_ACTION_MARK, unsafe_allow_html=True)
        return bool(
            st.button(label, key=key, help=help, type=button_type, use_container_width=True, **kwargs)
        )


def _hardware_configuration_mode_toggle() -> str:
    """Load existing vs Build new: paired wide buttons (selected = primary), not red load/save styling."""
    st.session_state.setdefault('gui_config_setup_mode', 'Load existing')
    mode = str(st.session_state.get('gui_config_setup_mode', 'Load existing'))
    if mode not in ('Load existing', 'Build new'):
        mode = 'Load existing'
        st.session_state['gui_config_setup_mode'] = mode
    st.caption('**Build new** writes a new `.py` from the fields below, then **Write config file and load**.')
    c1, c2 = st.columns(2)
    with c1:
        sel_load = mode == 'Load existing'
        if st.button(
            'Load existing',
            key='gui_hw_cfg_mode_load',
            use_container_width=True,
            type='primary' if sel_load else 'secondary',
            help='Use a module already on disk (dropdown or typed name).',
        ):
            if not sel_load:
                st.session_state['gui_config_setup_mode'] = 'Load existing'
                st.rerun()
    with c2:
        sel_build = mode == 'Build new'
        if st.button(
            'Build new',
            key='gui_hw_cfg_mode_build',
            use_container_width=True,
            type='primary' if sel_build else 'secondary',
            help='Fill channel and DAQ fields, then write a new `.py` and load it.',
        ):
            if not sel_build:
                st.session_state['gui_config_setup_mode'] = 'Build new'
                st.rerun()
    return str(st.session_state.get('gui_config_setup_mode', 'Load existing'))


def _shared_experiment_dir() -> str:
    """Protocol & biometrics JSON files use the same folder as **Data folder** (section 2) when it exists."""
    d = str(st.session_state.get('gui_data_folder') or '').strip()
    if d:
        norm = os.path.normpath(d)
        if os.path.isdir(norm):
            return norm
    return os.path.normpath(default_templates_dir(_ROOT))


def _protocol_template_option_label(p: Optional[str]):
    if p is None:
        return '— Choose a template —'
    return f'{template_display_label(p)}  ·  {os.path.basename(p)}'


def _biometrics_template_option_label(p: Optional[str]):
    if p is None:
        return '— Choose a biometrics file —'
    return biometrics_template_display_label(p)
from bender_h5_export import export_primary_h5, save_universal_qc_figure  # noqa: E402
from bender_h5_plot_helpers import (  # noqa: E402
    align_xy,
    h5_custom_plot_summary,
    list_h5_plot_variables,
    read_h5_series,
)

MOTION_TYPES = frozenset(
    {'dynamic', 'frequency_sweep', 'frequency_step', 'curvature_step', 'step_change'}
)

# Extra fields for motion-series protocols (not fully enumerated in get_dispatch_schema).
MOTION_GUI_FIELDS = [
    ('duration', 'float', 'Duration (s)'),
    ('all_freqs', 'list_float', 'Frequencies Hz (comma-separated, e.g. 1, 5)'),
    (
        'all_amps',
        'list_float',
        'Amplitudes (comma-separated; see all_amps_mode: strain=decimal ε, strain_pct=percent)',
    ),
    ('all_amps_mode', 'select', 'Amplitude interpretation (strain vs strain_pct, etc.)'),
    ('cycles_per_step', 'int', 'Cycles per step'),
    ('n_end_cycles', 'int', 'End cycles'),
    ('randomize', 'bool', 'Randomize order'),
    ('random_seed', 'optional_int', 'Random seed (empty = None)'),
    ('stim_cycles_in_step', 'list_int', 'Stim cycle indices per step (comma-separated, e.g. 2, 3)'),
    ('is_stim', 'bool', 'Stimulation enabled'),
    ('stim_pulse_rate', 'float', 'Stim pulse rate (Hz)'),
    ('S1volts', 'float', 'S1 stim voltage'),
    ('S2volts', 'float', 'S2 stim voltage'),
]

STEP_CHANGE_EXTRA = [
    ('step_change_frequencies', 'list_float', 'step_change frequencies (comma-separated)'),
    ('step_change_curves', 'list_float', 'step_change curves / amplitudes (comma-separated)'),
    ('step_change_cycles_per_step', 'list_int', 'cycles per step (comma-separated)'),
]

# Required by :meth:`Bender.make_cycles_frequency_sweep` but not used by other motion modes.
FREQUENCY_SWEEP_ONLY_FIELDS = [
    (
        'amplitude_frequency_exponent',
        'float',
        'Amplitude–frequency exponent α (command θ ∝ f^α vs sweep start)',
    ),
]

_MOTION_ROW_BY_NAME = {row[0]: row for row in MOTION_GUI_FIELDS}


def _motion_parameter_rows(test_type: str):
    """
    Rows (name, kind, label) shown for each motion ``test_type``, aligned with
    :meth:`Bender.run_experiment` branches (dynamic vs sweep/step vs step_change).
    """
    # Shared blocks (see ``organize_cycles`` + stim wiring in ``run_experiment``).
    freq_amp_mode = [
        _MOTION_ROW_BY_NAME['all_freqs'],
        _MOTION_ROW_BY_NAME['all_amps'],
        _MOTION_ROW_BY_NAME['all_amps_mode'],
    ]
    organize_block = [
        _MOTION_ROW_BY_NAME['cycles_per_step'],
        _MOTION_ROW_BY_NAME['n_end_cycles'],
        _MOTION_ROW_BY_NAME['randomize'],
        _MOTION_ROW_BY_NAME['random_seed'],
        _MOTION_ROW_BY_NAME['stim_cycles_in_step'],
    ]
    stim_block = [
        _MOTION_ROW_BY_NAME['is_stim'],
        _MOTION_ROW_BY_NAME['stim_pulse_rate'],
        _MOTION_ROW_BY_NAME['S1volts'],
        _MOTION_ROW_BY_NAME['S2volts'],
    ]
    if test_type == 'dynamic':
        # Duration comes from ``sum(period_by_cycle)`` after ``organize_cycles``, not ``self.duration``.
        return [*freq_amp_mode, *organize_block, *stim_block]
    if test_type == 'frequency_sweep':
        return [
            _MOTION_ROW_BY_NAME['duration'],
            *freq_amp_mode,
            *FREQUENCY_SWEEP_ONLY_FIELDS,
            *stim_block,
        ]
    if test_type in ('frequency_step', 'curvature_step'):
        return [_MOTION_ROW_BY_NAME['duration'], *freq_amp_mode, *stim_block]
    if test_type == 'step_change':
        return [*list(STEP_CHANGE_EXTRA), *stim_block]
    return list(MOTION_GUI_FIELDS)


ALL_AMPS_MODE_OPTIONS = ('strain', 'strain_pct', 'curvature', 'angle')

DATA_FOLDER_HELP = (
    'Choose the folder where experiment files live. Enter the folder path only — do not put the file name here. '
    'Runs, exports, protocol templates, and biometrics files can all use this folder.'
)
DATA_FILE_NAME_HELP = (
    'This is the name of your saved measurements file (HDF5). Enter only the file name, not the full path. '
    'The app joins it with the data folder field above to build where data is saved. You may type .h5 or leave it off. '
    'If that exact file already exists, the app uses a new name like my_run_001.h5 so nothing is overwritten.'
)

ISOMETRIC_STIM_JSON_HELP = (
    'Optional JSON for **each isometric step** (merged onto defaults). Common keys: **ramp_duration_s**, '
    '**hold_duration_s**, **settle_before_stim_s** (quiet time at target angle before stim), '
    '**stim_duration_s** (use null to stim through the rest of the hold), **inter_step_interval_s** '
    '(idle time between steps; 0 = back-to-back), **is_stim**, **stim_pulse_rate** (Hz), **stim_voltage** (V), '
    '**device_name** (null = use your NI config). Example: {"ramp_duration_s": 2, "hold_duration_s": 5, "stim_pulse_rate": 75}'
)
ISOMETRIC_STIM_OVERRIDES_HELP = (
    'Rare. JSON merged into stim **routing** (not timing): e.g. **recruitment**, **lateral_mode**, '
    '**bilateral_mirror_motor**, **bilateral_sequential_left_frac**. Leave `{}` unless you need overrides.'
)
ISOVELOCITY_STIM_JSON_HELP = (
    'Optional JSON merged onto isovelocity defaults: **settle_before_stim_s**, **stim_duration_s** (null = rest of iso), '
    '**pre_iso_stim_duration_s**, **is_stim**, **stim_pulse_rate**, **stim_voltage**, **device_name**, '
    'and optionally **iso_duration_s** / **pre_hold_s** to override the main fields.'
)
ISOVELOCITY_STIM_OVERRIDES_HELP = ISOMETRIC_STIM_OVERRIDES_HELP

RANDOM_SEED_HELP = (
    'Only used when **randomize order** is checked. Enter an integer for a **reproducible** shuffle (same order '
    'every run); leave empty for a **different** random order each time.'
)

RECRUITMENT_FIELD_HELP = (
    'How stimulation and/or motor commands are routed across left vs right: simultaneous, sequential halves of a '
    'step, or unilateral (one side).'
)

ISOVELOCITY_WIDGET_LABEL = {
    'isovelocity_min_vel': 'Minimum angular velocity (deg/s)',
    'isovelocity_max_vel': 'Maximum angular velocity (deg/s)',
    'isovelocity_starting_strain': 'Starting posture (strain / curvature / angle)',
    'isovelocity_num_steps': 'Number of velocity steps',
    'isovelocity_starting_strain_mode': 'Unit for starting posture value',
    'isovelocity_randomize': 'Randomize order of velocity steps',
    'isovelocity_random_seed': 'Random seed (optional)',
    'isovelocity_iso_duration_s': 'Constant-velocity bend duration (s)',
    'isovelocity_pre_hold_s': 'Pre-hold at starting angle (s)',
}

ISOVELOCITY_FIELD_HELP = {
    'isovelocity_min_vel': 'Lower end of the angular velocity sweep (deg/s). Steps are spaced linearly between min and max.',
    'isovelocity_max_vel': 'Upper end of the angular velocity sweep (deg/s).',
    'isovelocity_starting_strain': 'Initial posture before each velocity step; interpreted with **starting strain mode** (e.g. decimal ε vs motor angle).',
    'isovelocity_num_steps': 'How many commanded angular velocities between min and max (inclusive). Each is a separate trial from the same start posture.',
    'isovelocity_starting_strain_mode': 'Whether **starting posture** is decimal ε, percent strain, curvature κ, or motor angle (deg).',
    'isovelocity_randomize': 'Shuffle the order of velocity steps. Use **random seed** if you want the same order every time.',
    'isovelocity_random_seed': RANDOM_SEED_HELP,
    'isovelocity_iso_duration_s': 'After the quiet **pre-hold**, how long (seconds) the motor holds **constant angular velocity** before the step ends. This is not the pre-hold time—that is **pre-hold**.',
    'isovelocity_pre_hold_s': 'Time (s) at the starting posture before each trial’s constant-velocity segment begins.',
    'recruitment': RECRUITMENT_FIELD_HELP,
    'lateral_mode': (
        'Expert only. Leave **blank** unless you need a custom stim-router label. Normal left/right/bilateral behavior '
        'is set with **recruitment** above; this field overrides that name inside merged stim JSON.'
    ),
    'bilateral_mirror_motor': (
        'When **recruitment** is **bilateral sequential**, mirror the commanded bend between the first and second '
        'half of each step (left vs right). Does **not** turn on bilateral testing by itself—choose **recruitment** for '
        'which sides are used.'
    ),
    'bilateral_sequential_left_frac': 'Share of each step spent on the **left** side before the right (0–1).',
}

ISOMETRIC_FIELD_HELP = {
    'isometric_initial': 'First target in the sweep; units follow **isometric mode** (strain, curvature, or angle).',
    'isometric_final': 'Last target in the sweep (same units as initial).',
    'isometric_num_steps': 'Number of steps between initial and final (endpoints included).',
    'isometric_mode': 'How **initial** / **final** are interpreted: decimal ε, percent strain, κ, or motor angle.',
    'isometric_randomize': 'Shuffle step order. Use **random seed** for a fixed order across runs.',
    'isometric_random_seed': RANDOM_SEED_HELP,
    'isometric_inter_step_interval_s': 'Idle time after each step’s acquisition before the next ramp (0 = back-to-back).',
    'recruitment': RECRUITMENT_FIELD_HELP,
    'lateral_mode': (
        'Expert only. Leave **blank** unless you need a custom stim-router label. **Recruitment** sets routine '
        'left/right/bilateral behavior; this overrides the router name in merged stim params.'
    ),
    'bilateral_mirror_motor': (
        'When **recruitment** is **bilateral sequential**, mirror the hold/ramp between left and right halves of each '
        'step. Which sides are active is still set by **recruitment**, not this checkbox alone.'
    ),
    'bilateral_sequential_left_frac': 'Share of each step on the **left** side before the right (0–1).',
}

MOTION_FIELD_HELP = {
    'random_seed': RANDOM_SEED_HELP,
    'duration': 'Total motion timeline length (s) for protocols that use a single continuous run.',
    'all_freqs': 'Comma-separated list of frequencies (Hz) for each step or segment.',
    'all_amps': 'Comma-separated amplitudes; meaning depends on **amplitude interpretation** (strain vs angle, etc.).',
    'all_amps_mode': 'How each entry in **amplitudes** is converted to curvature / motor angle.',
    'cycles_per_step': 'Oscillation cycles at each frequency step before advancing.',
    'n_end_cycles': 'Extra cycles at the last frequency at the end of the run.',
    'randomize': 'Shuffle the order of frequency/amplitude steps when applicable.',
    'stim_cycles_in_step': 'Which cycle indices (1-based) within each step receive stimulation.',
    'is_stim': 'Enable patterned stimulation during the motion.',
    'stim_pulse_rate': 'Carrier pulse rate for stimulation (Hz).',
    'S1volts': 'Stimulus channel 1 amplitude (V).',
    'S2volts': 'Stimulus channel 2 amplitude (V).',
    'amplitude_frequency_exponent': 'For frequency sweep: exponent α so amplitude ∝ f^α relative to the sweep start.',
    'step_change_frequencies': 'Comma-separated frequencies for each step-change segment.',
    'step_change_curves': 'Comma-separated curvature / amplitude targets per segment.',
    'step_change_cycles_per_step': 'Comma-separated integers: cycles per step-change segment.',
}


def _format_strain_or_amp_mode(opt: str) -> str:
    """User-facing labels: make decimal vs percent strain explicit."""
    o = str(opt)
    if o == 'strain':
        return 'strain — decimal ε (0.05 = 5%)'
    if o == 'strain_pct':
        return 'strain_pct — percent scale (5 = 5%)'
    if o == 'curvature':
        return 'curvature — κ (1/m)'
    if o == 'angle':
        return 'angle — motor (deg)'
    return o


RECRUITMENT_OPTIONS = (
    'bilateral_simultaneous',
    'bilateral_sequential',
    'left',
    'right',
)

BILATERAL_MIRROR_LABEL = 'Perform test on both sides (bilateral)'
LATERAL_MODE_LABEL = 'Stim routing override (optional; experts only)'


def _parse_float_list(s: str):
    s = (s or '').strip()
    if not s:
        return None
    out = []
    for part in s.split(','):
        part = part.strip()
        if part:
            out.append(float(part))
    return out if out else None


def _parse_int_list(s: str):
    s = (s or '').strip()
    if not s:
        return None
    out = []
    for part in s.split(','):
        part = part.strip()
        if part:
            out.append(int(round(float(part))))
    return out if out else None


def _widget_key(name: str) -> str:
    return f'fld_{name}'


def _get_session_value(b: Bender, name: str, default=None):
    if hasattr(b, name):
        v = getattr(b, name)
        if v is not None:
            return v
    return default


def _clip_action_words(text: str, max_words: int = 5) -> str:
    words = str(text).strip().split()
    if len(words) <= max_words:
        return str(text).strip()
    return ' '.join(words[:max_words])


def _st_error_actions(headline: str, actions: list[str], *, max_words: int = 5) -> None:
    """Show ``st.error`` with a headline and bulleted next-step phrases (≤ ``max_words`` each)."""
    parts = [f'**{headline.strip()}**']
    for a in actions:
        s = str(a).strip()
        if s:
            parts.append(f'- {_clip_action_words(s, max_words)}')
    st.error('\n\n'.join(parts))


def _st_error_detail(headline: str, actions: list[str], detail: str) -> None:
    """Like `_st_error_actions` plus an expander with full text (e.g. server or Python message)."""
    _st_error_actions(headline, actions)
    with st.expander('Details'):
        st.code(str(detail).strip() or '(empty)', language=None)


def _missing_fields_to_actions(missing: list[str]) -> list[str]:
    """Turn validator missing names into short imperative lines (capped for UI size)."""
    names = [str(n).strip() for n in missing if str(n).strip()]
    lines = [_clip_action_words(f'Set {n}', 5) for n in names[:12]]
    if len(names) > 12:
        lines.append('Review other required fields')
    return lines


def _render_field(b: Bender, name: str, kind: str, label: str, *, help_text: Optional[str] = None):
    """
    Render one widget; return value to setattr, or None to skip.

    Streamlit keeps keyed widget values in ``st.session_state``. Passing ``value=`` / ``index=``
    every rerun is ignored once the key exists, which can leave the UI stale (e.g. empty text)
    while ``Bender`` still holds good data. We only seed session state when the key is missing.
    """
    sk = _widget_key(name)
    cur = _get_session_value(b, name)
    h = help_text if help_text else None

    if kind == 'float':
        if sk not in st.session_state:
            st.session_state[sk] = float(cur) if cur is not None else 0.0
        return float(st.number_input(label, key=sk, format='%.6g', help=h))

    if kind == 'int':
        if sk not in st.session_state:
            st.session_state[sk] = int(cur) if cur is not None else 0
        return int(st.number_input(label, key=sk, step=1, help=h))

    if kind == 'optional_int':
        use_key = f'{sk}_use'
        if use_key not in st.session_state:
            st.session_state[use_key] = cur is not None
        use = st.checkbox(f'{label} (use custom)', key=use_key, help=h)
        if not use:
            return None
        if sk not in st.session_state:
            st.session_state[sk] = int(cur) if cur is not None else 0
        return int(st.number_input(label, key=sk, step=1))

    if kind == 'bool':
        if sk not in st.session_state:
            st.session_state[sk] = bool(cur) if cur is not None else False
        return bool(st.checkbox(label, key=sk, help=h))

    if kind == 'str':
        if sk not in st.session_state:
            st.session_state[sk] = str(cur or '')
        return str(st.text_input(label, key=sk, help=h))

    if kind == 'select':
        opts = list(ALL_AMPS_MODE_OPTIONS)
        dv = str(cur or 'strain')
        if dv not in opts:
            dv = 'strain'
        if sk not in st.session_state:
            st.session_state[sk] = dv
        return str(st.selectbox(label, opts, key=sk, format_func=_format_strain_or_amp_mode, help=h))

    if kind == 'list_float':
        if sk not in st.session_state:
            if cur is not None:
                try:
                    st.session_state[sk] = ', '.join(str(float(x)) for x in list(cur))
                except (TypeError, ValueError):
                    st.session_state[sk] = ''
            else:
                st.session_state[sk] = ''
        s = str(st.text_input(label, key=sk, help=h))
        return _parse_float_list(s)

    if kind == 'list_int':
        if sk not in st.session_state:
            if cur is not None:
                try:
                    st.session_state[sk] = ', '.join(str(int(x)) for x in list(cur))
                except (TypeError, ValueError):
                    st.session_state[sk] = ''
            else:
                st.session_state[sk] = ''
        s = str(st.text_input(label, key=sk, help=h))
        return _parse_int_list(s)

    if kind == 'json_dict':
        if sk not in st.session_state:
            st.session_state[sk] = json.dumps(cur, indent=2) if isinstance(cur, dict) else '{}'
        jh = help_text if help_text is not None else f'JSON object for `{name}` (see protocol docs).'
        s = str(st.text_area(label, height=120, key=sk, help=jh))
        try:
            return json.loads(s)
        except json.JSONDecodeError:
            _st_error_actions(
                f'Invalid JSON: {name}',
                ['Fix commas and brackets', 'Match expected JSON shape', 'Compare protocol docs'],
            )
            return None

    return None


def _clear_fld_session_keys():
    """Drop procedure widget keys so the next render re-seeds from ``Bender`` (e.g. new config)."""
    for k in list(st.session_state.keys()):
        if k.startswith('fld_'):
            del st.session_state[k]


def _friendly_error_actions(err: Exception, *, action: str) -> tuple[str, list[str]]:
    """Headline and bulleted next steps (≤5 words each) for common failures."""
    msg = str(err or '').strip()
    low = msg.lower()
    if 'export_primary_h5 requires bender.outputfile or outputfile' in low:
        return (
            'HDF5 path not set.',
            ['Open section 2', 'Set data folder', 'Set file name', 'Apply data path'],
        )
    if 'please enter a note first' in low:
        return ('No note text.', ['Type your note', 'Click append again'])
    if 'file not found' in low and 'h5' in low:
        return ('File missing.', ['Pick file from list', 'Check folder path'])
    if 'selected file is not .h5' in low:
        return ('Not an HDF5 file.', ['Select a .h5 file', 'Use section 8 list'])
    if 'kaleido' in low and action == 'save_qc':
        return (
            'QC PNG export failed.',
            ['Run pip install kaleido', 'Or use HTML output'],
        )
    if action == 'save_h5':
        return (
            'HDF5 save failed.',
            ['Open Technical details', 'Fix error shown there'],
        )
    if action == 'save_qc':
        return (
            'QC figure save failed.',
            ['Open Technical details', 'Fix error shown there'],
        )
    if action == 'run_experiment':
        if '-200077' in msg or 'sampclk_rate' in low or 'daqmx_sampclk_rate' in low:
            return (
                'DAQ rejected sample rate.',
                ['Check DAQ rates in config', 'Check motion timeline dt', 'Open Technical details'],
            )
        if '-200560' in msg or 'wait until done' in low:
            return (
                'DAQ acquisition timed out.',
                ['Check USB connection', 'Close other DAQ software', 'Retry run experiment'],
            )
        if '-50103' in msg or 'resource reserved' in low or 'specified resource is reserved' in low:
            return (
                'DAQ still reserved elsewhere.',
                ['Close NI MAX panels', 'Close extra Streamlit tabs', 'Wait then retry'],
            )
        return (
            'Experiment did not run.',
            ['Open Technical details', 'Fix error shown there'],
        )
    return (
        'Something went wrong.',
        ['Open Technical details', 'Fix error shown there'],
    )


def _show_friendly_error(err: Exception, *, action: str):
    h, bullets = _friendly_error_actions(err, action=action)
    _st_error_actions(h, bullets)
    with st.expander('Technical details'):
        st.exception(err)


def _needs_missing_calibration_confirmation(b: Bender) -> bool:
    """True when inertial-calibration use is requested but file path is absent or missing."""
    use_cal = bool(getattr(b, 'use_inertial_calibration', False))
    cal_file = str(getattr(b, 'inertial_calibration_file', '') or '').strip()
    if not use_cal:
        return False
    return (not cal_file) or (not os.path.isfile(cal_file))


def _append_note_to_h5_file(h5_path: str, note_text: str):
    """Append a dated note block to an existing .h5 file's post-trial notes attrs."""
    note = str(note_text or '').strip()
    if not note:
        raise ValueError('Please enter a note first.')
    if not os.path.isfile(h5_path):
        raise ValueError(f'File not found: {h5_path}')
    ts = pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')
    with h5py.File(h5_path, 'r+') as f:
        prev = str(f.attrs.get('post-trial notes', '') or '').strip()
        merged = note if not prev else f"{prev}\n\n--- Post-Experiment QC note ({ts}) ---\n{note}"
        f.attrs['post-trial notes'] = merged
        if '01_Metadata' in f:
            f['01_Metadata'].attrs['post-trial notes'] = merged


def _candidate_review_files(data_file_path: str):
    """Files in the data directory likely useful for review (h5/images/html)."""
    p = str(data_file_path or '').strip()
    if not p:
        return []
    d = os.path.dirname(p) or '.'
    if not os.path.isdir(d):
        return []
    exts = {'.h5', '.png', '.jpg', '.jpeg', '.webp', '.html'}
    out = []
    for n in sorted(os.listdir(d)):
        fp = os.path.join(d, n)
        if os.path.isfile(fp) and os.path.splitext(n.lower())[1] in exts:
            out.append(fp)
    return out


def _ensure_review_file_selection(files: list) -> None:
    if not files:
        return
    sel = st.session_state.get('gui_review_selected')
    if sel not in files:
        st.session_state['gui_review_selected'] = files[0]


def _read_qc_trial_index(b: Bender) -> int:
    """Which in-memory `trial_records` entry the QC figure uses (after Run / export)."""
    tr = list(getattr(b, 'trial_records', []) or [])
    n = len(tr)
    if n <= 1:
        return 0
    try:
        ix = int(st.session_state.get('gui_qc_trial_index'))
    except (TypeError, ValueError):
        ix = n - 1
    return int(max(0, min(ix, n - 1)))


def _qc_figure_base_path(b: Bender, selected_path: str, trial_idx: int) -> Optional[str]:
    """Output stem for `save_universal_qc_figure` when tying the QC file name to a chosen `.h5`."""
    p = str(selected_path or '').strip()
    if not p.lower().endswith('.h5'):
        return None
    proc = str(getattr(b, 'test_type', 'unknown'))
    stem = p[:-3]
    return stem + f'_{proc}_trial{int(trial_idx):03d}_qc'


def _compose_output_h5_path() -> str:
    """Full `.h5` path from **section 2** **Data folder** + **Data file name**."""
    folder = str(st.session_state.get('gui_data_folder') or '').strip()
    fn = str(st.session_state.get('gui_data_filename') or '').strip()
    if not fn:
        return ''
    if not fn.lower().endswith('.h5'):
        fn = fn + '.h5'
    if folder:
        return os.path.normpath(os.path.join(folder, fn))
    return os.path.normpath(fn)


def _output_path_anchor_for_review(b: Optional[Bender] = None) -> str:
    """Path used to locate the data directory (composed section 2 path, else ``b.outputfile``)."""
    p = _compose_output_h5_path().strip()
    if p:
        return p
    inst = b if b is not None else st.session_state.get('bender')
    if inst is not None:
        return str(getattr(inst, 'outputfile', '') or '').strip()
    return ''


def _normalize_config_module_name(raw: str) -> str:
    """Strip whitespace and optional ``.py`` suffix for ``importlib.import_module``."""
    s = str(raw or '').strip()
    if s.lower().endswith('.py'):
        s = s[:-3].strip()
    return s


def _selected_config_matches_bender(b: Bender, eff_raw: str) -> bool:
    """True if ``b`` was loaded from the same module name as ``eff_raw`` (dropdown + override)."""
    eff = _normalize_config_module_name(eff_raw)
    if not eff:
        return False
    loaded = _normalize_config_module_name(getattr(b, 'config_name', '') or '')
    return bool(loaded) and loaded == eff


def _apply_loaded_config_module(raw_mod: str) -> Optional[str]:
    """Instantiate ``Bender`` from a config module and refresh session. Returns an error message or ``None``."""
    _cm = _normalize_config_module_name(raw_mod)
    if not _cm:
        return 'Enter a config module name.'
    if _ROOT not in sys.path:
        sys.path.insert(0, _ROOT)
    try:
        st.session_state['bender'] = Bender(_cm)
        st.session_state['cfg_mod'] = _cm
        b0 = st.session_state['bender']
        outp0 = str(getattr(b0, 'outputfile', '') or '').strip()
        if outp0:
            n0 = os.path.normpath(outp0)
            st.session_state['gui_pending_data_folder'] = os.path.dirname(n0) or ''
            st.session_state['gui_pending_data_filename'] = os.path.basename(n0)
        else:
            st.session_state['gui_pending_data_folder'] = ''
            st.session_state['gui_pending_data_filename'] = ''
        _meta0 = getattr(b0, 'h5_protocol_metadata', {}) or {}
        st.session_state['gui_pending_genus_species'] = str(_meta0.get('genus_species', '') or '')
        st.session_state['gui_pending_specimen_id'] = str(_meta0.get('specimen_id', '') or '')
        st.session_state['gui_pending_post_notes'] = str(getattr(b0, 'post_trial_notes', '') or '')
        _init_biometrics_session_state(b0, force=True)
        _clear_fld_session_keys()
        st.session_state.pop('gui_tpl_bio_done', None)
        return None
    except ImportError as e:
        st.session_state.pop('bender', None)
        return (
            f'Could not import config module `{_cm}`.\n\n'
            f'- Enter the **module name only** (no `.py`), e.g. `jimenez_bender_config_A`.\n'
            f'- The file must be importable from the app folder: `{_ROOT}`\n'
            f'- Current working directory: `{os.getcwd()}`\n\n'
            f'Detail: {e}'
        )
    except Exception as e:
        st.session_state.pop('bender', None)
        return f'{type(e).__name__}: {e}'


def _sec1_apply_composed_path_to_bender() -> Optional[str]:
    """Copy **Data folder** + **Data file name** onto ``bender.outputfile``. Returns an error message or ``None``."""
    b1 = st.session_state.get('bender')
    if b1 is None:
        return 'Load hardware configuration first.'
    outp = _compose_output_h5_path().strip()
    if not outp:
        return 'Set **Data file name** first.'
    b1.outputfile = outp
    return None


def _paths_equal_norm(a: str, b: str) -> bool:
    if not str(a).strip() or not str(b).strip():
        return False
    try:
        return os.path.normcase(os.path.normpath(a)) == os.path.normcase(os.path.normpath(b))
    except Exception:
        return str(a).strip() == str(b).strip()


def _session_float(key: str) -> Optional[float]:
    try:
        v = st.session_state.get(key)
        if v is None:
            return None
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def _collect_input_sanity_warnings() -> list[str]:
    """Heuristic checks: blanks, zeros, and sub‑1 mm values that often mean a units slip."""
    out: list[str] = []

    df = str(st.session_state.get('gui_data_folder') or '').strip()
    fn = str(st.session_state.get('gui_data_filename') or '').strip()
    if not df:
        out.append('Data folder is empty.')
    else:
        try:
            if not os.path.isdir(os.path.normpath(df)):
                out.append('Data folder is not a valid directory.')
        except Exception:
            out.append('Data folder path is invalid.')
    if not fn:
        out.append('Data file name is empty.')

    if not str(st.session_state.get('gui_specimen_id') or '').strip():
        out.append('Specimen ID is blank.')
    if not str(st.session_state.get('gui_genus_species') or '').strip():
        out.append('Genus-species is blank.')

    for label, key in (
        ('Clamp spacing (dclamp)', 'bio_dclamp'),
        ('Cross-section width', 'bio_xsec'),
        ('Cross-section height', 'bio_xsec_height'),
        ('Profile length', 'bio_prof_L'),
    ):
        v = _session_float(key)
        if v is not None and 0 < v < 1.0:
            out.append(f'{label} is {v:g} mm (< 1 mm — check units).')

    for label, key in (('Clamp spacing', 'bio_dclamp'), ('TL', 'bio_fishlen_TL'), ('SL', 'bio_fishlen_SL')):
        v = _session_float(key)
        if v is not None and v <= 0:
            out.append(f'{label} is zero or negative.')

    m = _session_float('bio_fishmass')
    if m is not None and m <= 0:
        out.append('Mass is zero or negative.')

    for label, key in (
        ('Profile proximal H', 'bio_prof_ph'),
        ('Profile proximal W', 'bio_prof_pw'),
        ('Profile distal H', 'bio_prof_dh'),
        ('Profile distal W', 'bio_prof_dw'),
    ):
        v = _session_float(key)
        if v is not None and 0 < v < 1.0:
            out.append(f'{label} is {v:g} mm (< 1 mm).')

    # Dedupe while keeping order
    seen: set[str] = set()
    uniq: list[str] = []
    for w in out:
        if w not in seen:
            seen.add(w)
            uniq.append(w)
    return uniq


def _render_sidebar_stepwise_progress() -> None:
    """Stepwise HW + path summary (sidebar)."""
    b = st.session_state.get('bender')
    st.markdown('**Progress**')
    if b is None:
        st.caption('Step 1: load hardware.')
        return
    st.caption(f"Module `{getattr(b, 'config_name', '?')}`")
    composed = _compose_output_h5_path().strip()
    applied = str(getattr(b, 'outputfile', '') or '').strip()
    if not applied:
        st.caption('Step 2: apply data path.')
    elif not composed:
        st.caption('Step 2: set a file name.')
    elif _paths_equal_norm(applied, composed):
        st.caption(f"HDF5 → `{os.path.basename(applied)}`")
    else:
        st.caption('Step 2: form path ≠ experiment — use **Set data file path**.')
    st.divider()


def _render_sidebar_input_checks() -> None:
    """Suspicious / missing inputs when hardware is loaded."""
    if st.session_state.get('bender') is None:
        return
    warns = _collect_input_sanity_warnings()
    if not warns:
        return
    st.markdown('**Checks**')
    st.warning('\n\n'.join(f'- {w}' for w in warns[:10]))
    if len(warns) > 10:
        st.caption(f'… and {len(warns) - 10} more.')
    st.divider()


def _seed_cfg_build_from_base(base: str) -> None:
    """Fill ``gui_cfg_bld_*`` widget defaults from an existing config module."""
    try:
        d = read_base_defaults(base)
    except Exception:
        return
    st.session_state['gui_cfg_bld_forcetorque_calibration_file'] = str(d['forcetorque_calibration_file'])
    st.session_state['gui_cfg_bld_positive_motor_direction'] = str(d['positive_motor_direction'])
    st.session_state['gui_cfg_bld_specimen_lateral_index'] = int(d['specimen_lateral_index_on_positive_motor_side'])
    st.session_state['gui_cfg_bld_primary_bending_axis'] = str(d['primary_bending_axis'])
    st.session_state['gui_cfg_bld_device_name'] = str(d['device_name'])
    st.session_state['gui_cfg_bld_daq_ai_sr'] = float(d['daq_ai_sample_rate_hz'])
    st.session_state['gui_cfg_bld_daq_ao_sr'] = float(d['daq_ao_do_sample_rate_hz'])
    st.session_state['gui_cfg_bld_motor_steps'] = int(d['motor_full_steps_per_rev'])
    st.session_state['gui_cfg_bld_motor_gear'] = int(d['motor_gear_ratio'])
    st.session_state['gui_cfg_bld_stim_channels'] = ', '.join(str(x) for x in d['stim_channels'])
    st.session_state['gui_cfg_bld_motor_port'] = str(d['motor_port'])
    st.session_state['gui_cfg_bld_encoder_chan'] = str(d['encoder_chan'])
    st.session_state['gui_cfg_bld_SG_chan'] = ', '.join(str(x) for x in d['SG_chan'])
    st.session_state['gui_cfg_bld_SG_name'] = ', '.join(str(x) for x in d['SG_name'])
    st.session_state['gui_cfg_bld_S1side'] = str(d['S1side'])
    st.session_state['gui_cfg_bld_S2side'] = str(d['S2side'])
    st.session_state['gui_cfg_bld_use_sono'] = bool(d['use_sono'])
    st.session_state['gui_cfg_bld_sono_channel'] = ', '.join(str(x) for x in d['sono_channel'])
    st.session_state['gui_cfg_bld_sono_name'] = ', '.join(str(x) for x in d['sono_name'])
    st.session_state['gui_cfg_bld_encoder_ppr'] = int(d['encoder_pulses_per_rev'])
    st.session_state['gui_cfg_bld_motor_axis'] = str(d['motor_axis'])
    st.session_state['gui_cfg_bld_bending_axis_sensor'] = str(d['bending_axis_sensor'])
    st.session_state['gui_cfg_bld_bending_axis_specimen'] = str(d['bending_axis_specimen'])
    st.session_state['gui_cfg_bld_stim_monitor_chan'] = ', '.join(str(x) for x in d['stim_monitor_chan'])
    st.session_state['gui_cfg_bld_stim_monitor_name'] = ', '.join(str(x) for x in d['stim_monitor_name'])
    st.session_state['gui_cfg_bld_sono_internal_samplefreq'] = int(d['sono_internal_samplefreq'])
    st.session_state['gui_cfg_bld_sono_cal_left'] = ', '.join(str(x) for x in d['sono_cal_left'])
    st.session_state['gui_cfg_bld_sono_cal_right'] = ', '.join(str(x) for x in d['sono_cal_right'])
    st.session_state['gui_cfg_bld_amp_step_vel'] = int(d['amp_step_vel'])
    rm = str(d['ramp_mode_default'])
    st.session_state['gui_cfg_bld_ramp_mode_default'] = rm if rm in ('linear', 'exponential') else 'linear'
    st.session_state['gui_cfg_bld_waitbefore'] = float(d['waitbefore'])
    st.session_state['gui_cfg_bld_waitafter'] = float(d['waitafter'])
    st.session_state['gui_cfg_bld_rampdur'] = float(d['rampdur'])
    st.session_state['gui_cfg_bld_prepoststim_dur'] = float(d['prepoststim_dur'])
    st.session_state['gui_cfg_bld_prepoststim_sep'] = float(d['prepoststim_sep'])
    st.session_state['gui_cfg_bld_prestim_time'] = float(d['prestim_time'])
    st.session_state['gui_cfg_bld_poststim_time'] = float(d['poststim_time'])


def _maybe_seed_cfg_build_fields() -> None:
    base = str(st.session_state.get('gui_cfg_build_base') or 'jimenez_bender_config_A')
    if st.session_state.get('gui_cfg_build_seeded_for') == base:
        return
    _seed_cfg_build_from_base(base)
    st.session_state['gui_cfg_build_seeded_for'] = base


def _flush_pending_load_config_session():
    """Apply deferred session updates before widgets bind to these keys (Streamlit restriction)."""
    if 'gui_pending_genus_species' in st.session_state:
        st.session_state['gui_genus_species'] = st.session_state.pop('gui_pending_genus_species')
    if 'gui_pending_specimen_id' in st.session_state:
        st.session_state['gui_specimen_id'] = st.session_state.pop('gui_pending_specimen_id')
    if 'gui_pending_data_folder' in st.session_state:
        st.session_state['gui_data_folder'] = st.session_state.pop('gui_pending_data_folder')
    if 'gui_pending_data_filename' in st.session_state:
        st.session_state['gui_data_filename'] = st.session_state.pop('gui_pending_data_filename')
    if 'gui_pending_post_notes' in st.session_state:
        st.session_state['gui_post_notes'] = st.session_state.pop('gui_pending_post_notes')


def _ensure_gui_data_path_session_keys():
    """Seed folder/filename from legacy single-field ``gui_outputfile`` if needed."""
    if 'gui_data_folder' in st.session_state and 'gui_data_filename' in st.session_state:
        return
    leg = str(st.session_state.get('gui_outputfile', '') or '').strip()
    if leg:
        norm = os.path.normpath(leg)
        st.session_state['gui_data_folder'] = os.path.dirname(norm) or ''
        st.session_state['gui_data_filename'] = os.path.basename(norm)
    else:
        st.session_state['gui_data_folder'] = ''
        st.session_state['gui_data_filename'] = ''


def _sync_genus_species_to_bender(b: Bender) -> None:
    """Store identity, notebook-style specimen metadata, and protocol IDs for HDF5 export."""
    meta = dict(getattr(b, 'h5_protocol_metadata', {}) or {})
    meta['genus_species'] = str(st.session_state.get('gui_genus_species') or '').strip()
    sid = str(st.session_state.get('gui_specimen_id') or '').strip()
    meta['specimen_id'] = sid
    setattr(b, 'fishcode', sid)

    def _str_attr(sess_key: str, bender_attr: str) -> None:
        if sess_key not in st.session_state:
            return
        s = str(st.session_state[sess_key] or '').strip()
        meta[bender_attr] = s
        setattr(b, bender_attr, s)

    def _float_attr(sess_key: str, bender_attr: str) -> None:
        if sess_key not in st.session_state:
            return
        try:
            v = float(st.session_state[sess_key])
        except (TypeError, ValueError):
            return
        meta[bender_attr] = v
        setattr(b, bender_attr, v)

    _str_attr('bio_segment', 'segment')
    _float_attr('bio_fishmass', 'fishmass')
    _float_attr('bio_fishlen_TL', 'fishlen_TL')
    _float_attr('bio_fishlen_SL', 'fishlen_SL')
    _float_attr('bio_xsec_height', 'xsec_height')
    _float_attr('bio_dvert', 'dvert')
    _float_attr('bio_dhoriz', 'dhoriz')

    b.h5_protocol_metadata = meta


def _apply_pair(b: Bender, name: str, value):
    if value is None and name == 'random_seed':
        setattr(b, name, None)
        return
    if value is None:
        return
    if name == 'lateral_mode' and isinstance(value, str) and not value.strip():
        setattr(b, 'lateral_mode', None)
        return
    setattr(b, name, value)


def _apply_form_updates(b: Bender, updates: dict, tt: str):
    if tt == 'calibration':
        emb = st.session_state.get('gui_calibration_embedded_base')
        if isinstance(emb, dict):
            proc = emb.get('procedure')
            if isinstance(proc, dict):
                for k, v in proc.items():
                    _apply_pair(b, k, v)
    for k, v in updates.items():
        _apply_pair(b, k, v)
    b.test_type = tt


def _apply_procedure_form_to_bender(b: Bender, updates: dict, tt: str) -> None:
    """Sync biometrics flags, copy procedure fields onto ``b``, and mirror any QC note text from session."""
    _sync_biometric_flags_from_session(b)
    _apply_form_updates(b, updates, tt)
    _sync_genus_species_to_bender(b)
    _pn = str(st.session_state.get('gui_post_notes') or '').strip()
    if _pn:
        b.post_trial_notes = _pn


def _consume_pending_protocol_template(valid_test_types: list) -> None:
    """Apply a template queued by **Load template** (runs before procedure widgets render)."""
    path = st.session_state.pop('gui_pending_protocol_template_path', None)
    if not path:
        return
    try:
        data = load_protocol_template(path)
        ok, msg = apply_template_to_session_state(
            st.session_state,
            data,
            widget_key=_widget_key,
            valid_test_types=list(valid_test_types),
        )
        st.session_state['gui_protocol_load_feedback'] = (ok, msg)
    except OSError as e:
        st.session_state['gui_protocol_load_feedback'] = (False, f'Could not read template: {e}')
    except json.JSONDecodeError as e:
        st.session_state['gui_protocol_load_feedback'] = (False, f'Invalid JSON: {e}')
    except Exception as e:
        st.session_state['gui_protocol_load_feedback'] = (False, f'{type(e).__name__}: {e}')


def _consume_pending_biometrics_template() -> None:
    path = st.session_state.pop('gui_pending_biometrics_path', None)
    if not path:
        return
    try:
        data = load_biometrics_template(path)
        ok, msg = apply_biometrics_template_to_session(st.session_state, data)
        st.session_state['gui_biometrics_load_feedback'] = (ok, msg)
        if ok:
            st.session_state.pop('gui_tpl_bio_done', None)
    except OSError as e:
        st.session_state['gui_biometrics_load_feedback'] = (False, f'Could not read file: {e}')
    except json.JSONDecodeError as e:
        st.session_state['gui_biometrics_load_feedback'] = (False, f'Invalid JSON: {e}')
    except Exception as e:
        st.session_state['gui_biometrics_load_feedback'] = (False, f'{type(e).__name__}: {e}')


def _apply_intrinsic_biometrics_to_bender(b: Bender) -> None:
    """TL/SL, mass, density, environment, offsets, identity → ``b``."""
    b.fishlen_TL = float(st.session_state['bio_fishlen_TL'])
    b.fishlen_SL = float(st.session_state['bio_fishlen_SL'])
    b.fishmass = float(st.session_state['bio_fishmass'])
    b.specimen_profile_density_g_per_mm3 = float(st.session_state['bio_prof_rho'])
    b.temp_C_room = float(st.session_state['bio_temp_room'])
    b.temp_C_tank = float(st.session_state['bio_temp_tank'])
    b.dvert = float(st.session_state['bio_dvert'])
    b.dhoriz = float(st.session_state['bio_dhoriz'])
    b.fishcode = str(st.session_state.get('gui_specimen_id') or '')
    b.segment = str(st.session_state.get('bio_segment') or '')
    _sync_genus_species_to_bender(b)


def _apply_clamp_geometry_to_bender(b: Bender) -> None:
    """Clamp spacing, bend position, cross-section at the test segment → ``b``."""
    b.dclamp = float(st.session_state['bio_dclamp'])
    b.test_segment_length_mm = float(st.session_state['bio_dclamp'])
    b.dbend = float(st.session_state['bio_dbend'])
    b.test_segment_position_mm = float(st.session_state['bio_dbend'])
    b.xsec_width = float(st.session_state['bio_xsec'])
    b.xsec_height = float(st.session_state['bio_xsec_height'])


def _apply_mounted_profile_inertial_to_bender(b: Bender) -> None:
    """Tapered outline + outline length + density for rotational inertia model → ``b``."""
    stations = b.make_profile_stations(
        st.session_state['bio_prof_ph'],
        st.session_state['bio_prof_pw'],
        st.session_state['bio_prof_dh'],
        st.session_state['bio_prof_dw'],
    )
    b.set_profiled_specimen_inertial_model(
        stations,
        st.session_state['bio_prof_L'],
        st.session_state['bio_prof_rho'],
        clamp_offset_mm=float(st.session_state['bio_prof_clamp']),
        num_samples=int(st.session_state['bio_prof_samples']),
    )


def _apply_all_biometrics_to_bender(b: Bender) -> None:
    """Copy intrinsic, clamp, mounted-profile, and inertial flags from section 3 biometrics session state onto ``b``."""
    _sync_biometric_flags_from_session(b)
    _apply_intrinsic_biometrics_to_bender(b)
    _apply_clamp_geometry_to_bender(b)
    _apply_mounted_profile_inertial_to_bender(b)


def _sync_biometric_flags_from_session(b: Bender):
    """Copy biometrics panel session state onto ``b`` (flags + specimen geometry aliases)."""
    if 'bio_use_theoretical_inertial' in st.session_state:
        b.use_theoretical_inertial_correction = bool(st.session_state['bio_use_theoretical_inertial'])
    if 'bio_dclamp' in st.session_state:
        v = float(st.session_state['bio_dclamp'])
        b.dclamp = v
        b.test_segment_length_mm = v
    if 'bio_dbend' in st.session_state:
        v = float(st.session_state['bio_dbend'])
        b.dbend = v
        b.test_segment_position_mm = v
    if 'bio_xsec' in st.session_state:
        b.xsec_width = float(st.session_state['bio_xsec'])
    if 'bio_temp_room' in st.session_state:
        b.temp_C_room = float(st.session_state['bio_temp_room'])
    if 'bio_temp_tank' in st.session_state:
        b.temp_C_tank = float(st.session_state['bio_temp_tank'])
    if 'bio_xsec_height' in st.session_state:
        b.xsec_height = float(st.session_state['bio_xsec_height'])
    if 'bio_dvert' in st.session_state:
        b.dvert = float(st.session_state['bio_dvert'])
    if 'bio_dhoriz' in st.session_state:
        b.dhoriz = float(st.session_state['bio_dhoriz'])
    if 'gui_specimen_id' in st.session_state:
        b.fishcode = str(st.session_state.get('gui_specimen_id') or '')
    if 'bio_segment' in st.session_state:
        b.segment = str(st.session_state['bio_segment'] or '')
    if 'bio_fishmass' in st.session_state:
        b.fishmass = float(st.session_state['bio_fishmass'])
    if 'bio_fishlen_TL' in st.session_state:
        b.fishlen_TL = float(st.session_state['bio_fishlen_TL'])
    if 'bio_fishlen_SL' in st.session_state:
        b.fishlen_SL = float(st.session_state['bio_fishlen_SL'])


def _init_biometrics_session_state(b: Bender, *, force: bool = False):
    """Seed Streamlit widget keys from ``b`` (``force`` overwrites after config reload)."""

    def _put(key, val):
        if force or key not in st.session_state:
            st.session_state[key] = val

    dc = getattr(b, 'dclamp', None)
    xw = getattr(b, 'xsec_width', None)
    _meta_b = getattr(b, 'h5_protocol_metadata', {}) or {}
    _put(
        'gui_specimen_id',
        str(_meta_b.get('specimen_id') or getattr(b, 'fishcode', '') or '').strip(),
    )
    _put('bio_segment', str(getattr(b, 'segment', '') or ''))
    _fm = getattr(b, 'fishmass', None)
    _put('bio_fishmass', float(_fm) if _fm is not None and math.isfinite(float(_fm)) else 0.0)
    _ftl = getattr(b, 'fishlen_TL', None)
    _put('bio_fishlen_TL', float(_ftl) if _ftl is not None and math.isfinite(float(_ftl)) else 0.0)
    _fsl = getattr(b, 'fishlen_SL', None)
    _put('bio_fishlen_SL', float(_fsl) if _fsl is not None and math.isfinite(float(_fsl)) else 0.0)
    _xh = getattr(b, 'xsec_height', None)
    _put(
        'bio_xsec_height',
        float(_xh) if _xh is not None and math.isfinite(float(_xh)) else (float(xw) if xw is not None else 8.0),
    )
    _put('bio_dvert', float(getattr(b, 'dvert', 0.0) or 0.0))
    _put('bio_dhoriz', float(getattr(b, 'dhoriz', 0.0) or 0.0))
    _put('bio_dclamp', float(dc) if dc is not None else 10.0)
    _put('bio_xsec', float(xw) if xw is not None else 8.0)
    _put('bio_dbend', float(getattr(b, 'dbend', 0.0) or 0.0))
    _put('bio_temp_room', float(getattr(b, 'temp_C_room', 22.0) or 22.0))
    _put('bio_temp_tank', float(getattr(b, 'temp_C_tank', 22.0) or 22.0))
    _put(
        'bio_use_theoretical_inertial',
        bool(getattr(b, 'use_theoretical_inertial_correction', False)),
    )
    # Simple 2-station profile defaults (mm)
    _put('bio_prof_L', float(getattr(b, 'specimen_profile_length_mm', 25.0) or 25.0))
    _put('bio_prof_rho', float(getattr(b, 'specimen_profile_density_g_per_mm3', 1.03e-3) or 1.03e-3))
    _put('bio_prof_ph', float(4.0))
    _put('bio_prof_pw', float(6.0))
    _put('bio_prof_dh', float(3.5))
    _put('bio_prof_dw', float(5.0))
    _put('bio_prof_clamp', float(getattr(b, 'specimen_profile_clamp_offset_mm', 20.0) or 20.0))
    _put('bio_prof_samples', int(getattr(b, 'specimen_profile_num_samples', 120) or 120))


def _nav_route() -> str:
    return str(st.session_state.get('gui_app_route') or 'landing')


def _stepwise_step() -> int:
    v = st.session_state.get('gui_stepwise_step', 0)
    try:
        return int(v)
    except (TypeError, ValueError):
        return 0


def _stepwise_on_data_file_path_step() -> bool:
    return _nav_route() == 'stepwise' and _stepwise_step() == 1


def _template_hide_config_build_new() -> bool:
    """In template mode, user said they have a saved config — skip the *Build new* path in section 1."""
    return _nav_route() == 'templates' and bool(st.session_state.get('gui_tpl_chk_config', True))


def _tpl_only_procedure() -> bool:
    """True when user has config + biometrics files but still needs to set the procedure in the app (no protocol JSON yet)."""
    if _nav_route() != 'templates':
        return False
    if not bool(st.session_state.get('gui_tpl_chk_config', True)):
        return False
    if not bool(st.session_state.get('gui_tpl_chk_biometrics', True)):
        return False
    # Third box: "I already have a protocol template" — when checked, use full sections (load protocol from file).
    return not bool(st.session_state.get('gui_tpl_have_protocol_template', False))


def _show_hw_config_section() -> bool:
    """Section 1 · hardware configuration (no data-folder block)."""
    if _nav_route() == 'templates' and _tpl_only_procedure():
        return False
    if _nav_route() == 'stepwise':
        return _stepwise_step() == 0
    return True


def _show_data_path_section() -> bool:
    """Section 2 · data folder & file name (+ **Load hardware configuration and data path** for load-existing)."""
    if _nav_route() == 'templates' and _tpl_only_procedure():
        return False
    if _nav_route() == 'stepwise':
        return _stepwise_step() == 1
    return True


def _show_full_sec2() -> bool:
    """Section 3 · biometrics."""
    if _nav_route() == 'templates' and _tpl_only_procedure():
        return False
    if _nav_route() == 'stepwise' and _stepwise_step() != 2:
        return False
    return True


def _show_sec3_through_6() -> bool:
    if _nav_route() != 'stepwise':
        return True
    return _stepwise_step() == 3


def _show_sec7_and_8() -> bool:
    if _nav_route() != 'stepwise':
        return True
    return _stepwise_step() == 4


def _render_landing_page() -> None:
    st.markdown(
        '<div class="bnd-landing-page" aria-hidden="true"></div>',
        unsafe_allow_html=True,
    )
    st.markdown(
        """
<style>
body:has(.bnd-landing-page) section[data-testid="stMain"] .block-container {
    padding-top: 1rem !important;
}
/* Landing: equal-tall action buttons (row is only buttons, so tops align). */
body:has(.bnd-landing-page) section[data-testid="stMain"] button[kind="primary"],
body:has(.bnd-landing-page) section[data-testid="stMain"] button[kind="secondary"] {
    min-height: 3.15rem !important;
}
</style>
""",
        unsafe_allow_html=True,
    )
    try:
        c0, c1, c2 = st.columns([2.4, 1, 1], vertical_alignment='center')
    except TypeError:
        c0, c1, c2 = st.columns([2.4, 1, 1])
    with c0:
        st.title('The CritterGripper App')
        st.caption('Non-destructive multiaxial platform for the mechanical testing of bending organisms.')
    with c1:
        if os.path.isfile(_LOGO_PATH):
            try:
                st.image(_LOGO_PATH, use_container_width=True)
            except TypeError:
                try:
                    st.image(_LOGO_PATH, use_column_width=True)
                except TypeError:
                    st.image(_LOGO_PATH, width=160)
    with c2:
        nsf = _nsf_logo_path()
        if nsf:
            st.markdown(
                f'<img src="{_img_data_uri(nsf)}" style="max-width:100%;height:auto;max-height:80px;" alt="NSF"/>',
                unsafe_allow_html=True,
            )
            st.caption('NSF support acknowledgment.')

    st.markdown('### Choose a workflow')
    st.caption('Same capabilities — pick the layout that fits how you like to work.')

    a, b, c = st.columns(3)
    with a:
        st.markdown('**Build from scratch**')
        st.caption(
            'All sections visible: hardware, data path, biometrics, experiment, preview, run, save, plots, notes.'
        )
    with b:
        st.markdown('**Templates**')
        st.caption('Load hardware `.py` and biometrics JSON from disk; optional protocol template checklist.')
    with c:
        st.markdown('**Step-by-step**')
        st.caption('One focus per step with a progress bar; good for training or checklist-style sessions.')

    st.markdown('<div style="height:0.85rem"></div>', unsafe_allow_html=True)
    ba, bb, bc = st.columns(3)
    with ba:
        if st.button('Start full workflow', key='land_scratch', use_container_width=True, type='primary'):
            st.session_state['gui_app_route'] = 'scratch'
            st.session_state.pop('gui_stepwise_step', None)
            st.rerun()
    with bb:
        if st.button('Template workflow', key='land_templates', use_container_width=True):
            st.session_state['gui_app_route'] = 'templates'
            st.session_state.pop('gui_stepwise_step', None)
            st.session_state.pop('gui_tpl_bio_done', None)
            st.session_state.pop('gui_tpl_cfg_autoloaded', None)
            st.rerun()
    with bc:
        if st.button('Step-by-step', key='land_stepwise', use_container_width=True):
            st.session_state['gui_app_route'] = 'stepwise'
            st.session_state['gui_stepwise_step'] = 0
            st.rerun()


def _render_app_chrome() -> None:
    if _nav_route() == 'stepwise':
        st.markdown('<div class="bnd-stepwise-active" aria-hidden="true"></div>', unsafe_allow_html=True)
        _inject_stepwise_compact_layout_css()
        # Wide-enough columns so "Home" stays one line; align row vertically.
        try:
            c_home, c_spacer, c_logo = st.columns([2.2, 7.6, 2.2], vertical_alignment='center')
        except TypeError:
            c_home, c_spacer, c_logo = st.columns([2.2, 7.6, 2.2])
        with c_home:
            if st.button('Home', key='gui_nav_home', type='secondary'):
                st.session_state['gui_app_route'] = 'landing'
                st.rerun()
        with c_spacer:
            pass
        with c_logo:
            lg0, lg1 = st.columns([1, 1])
            with lg1:
                if os.path.isfile(_LOGO_PATH):
                    try:
                        st.image(_LOGO_PATH, width=40)
                    except Exception:
                        pass
        _inject_load_save_button_theme()
        return

    h0, h1, h2, h3 = st.columns([1, 2.2, 1, 1.2])
    with h0:
        if st.button('← Home', key='gui_nav_home'):
            st.session_state['gui_app_route'] = 'landing'
            st.rerun()
    with h1:
        st.markdown('### The CritterGripper App')
        _mode = _nav_route()
        if _mode == 'scratch':
            st.caption('Full workflow — all sections below are visible.')
        elif _mode == 'templates':
            st.caption('Template mode — checklist at top controls which sections appear.')
    with h2:
        if os.path.isfile(_LOGO_PATH):
            try:
                st.image(_LOGO_PATH, width=120)
            except Exception:
                pass
    with h3:
        nsf = _nsf_logo_path()
        if nsf:
            st.markdown(
                f'<img src="{_img_data_uri(nsf)}" style="max-height:56px;width:auto;" alt="NSF"/>',
                unsafe_allow_html=True,
            )
    st.divider()
    _inject_load_save_button_theme()


def _render_sidebar() -> None:
    """Sidebar: stepwise progress, input checks, emergency stop."""
    with st.sidebar:
        if _nav_route() == 'stepwise':
            _render_sidebar_stepwise_progress()
        _render_sidebar_input_checks()
        st.markdown('### Emergency stop')
        st.caption('Stops DAQ tasks and resets the NI device. Safe to use if something slips.')
        if st.button(
            'Stop DAQ & reset NI device',
            key='gui_kill_daq',
            type='primary',
            use_container_width=True,
            help='Stops outputs and resets tasks.',
        ):
            dev = None
            if st.session_state.get('bender') is not None:
                dev = getattr(st.session_state['bender'], 'device_name', None)
            ok, msg = daq_emergency_stop(dev)
            if ok:
                st.success(msg)
            else:
                st.warning(msg)
        st.divider()


_STEPWISE_TAB_LABELS = (
    '1 · Hardware configuration',
    '2 · Data file path',
    '3 · Biometrics',
    '4 · Experiment & run',
    '5 · Visualize & notes',
)

# Short labels for tab buttons (keep narrow columns readable)
_STEPWISE_TAB_SHORT = (
    'Hardware',
    'Data',
    'Biometrics',
    'Experiment',
    'Visualize',
)


def _render_stepwise_rail() -> None:
    st.session_state.setdefault('gui_stepwise_step', 0)
    cur = int(st.session_state.get('gui_stepwise_step', 0) or 0)
    cur = max(0, min(4, cur))
    if cur != st.session_state.get('gui_stepwise_step'):
        st.session_state['gui_stepwise_step'] = cur

    with st.container(border=True):
        st.progress(min(1.0, (cur + 1) / 5.0), text=f'Step {cur + 1} / 5 · {_STEPWISE_TAB_SHORT[cur]}')

        tab_cols = st.columns(5, gap='small')
        for i in range(5):
            with tab_cols[i]:
                _lbl = f'{i + 1}. {_STEPWISE_TAB_SHORT[i]}'
                if st.button(
                    _lbl,
                    key=f'gui_sw_tab_{i}',
                    use_container_width=True,
                    type='primary' if i == cur else 'secondary',
                    disabled=i == cur,
                    help=_STEPWISE_TAB_LABELS[i],
                ):
                    st.session_state['gui_stepwise_step'] = i
                    st.rerun()

        n1, n2, _ = st.columns([1, 1, 4])
        with n1:
            if st.button(
                '← Previous',
                key='gui_sw_back',
                disabled=cur <= 0,
                use_container_width=True,
                help='Move to the previous step',
            ):
                st.session_state['gui_stepwise_step'] = cur - 1
                st.rerun()
        with n2:
            if st.button(
                'Next →',
                key='gui_sw_next',
                disabled=cur >= 4,
                use_container_width=True,
                help='Move to the next step',
            ):
                st.session_state['gui_stepwise_step'] = cur + 1
                st.rerun()

    if st.session_state.get('bender') is None and _stepwise_step() >= 2:
        with st.container(border=True):
            w1, w2 = st.columns([3, 1])
            with w1:
                st.warning(
                    'No hardware configuration is loaded yet. Go back to **step 1** to choose a module, then **step 2** '
                    'to set folder, file name, and **Set data file path** before continuing to biometrics.'
                )
            with w2:
                if st.button('Go to step 1', key='gui_stepwise_go_a', use_container_width=True):
                    st.session_state['gui_stepwise_step'] = 0
                    st.rerun()


def _cb_tpl_config_module_changed() -> None:
    st.session_state['gui_tpl_reload_config'] = True


def _render_template_procedure_strip() -> None:
    """Compact loaders: config + biometrics from disk, then experiment sections (procedure edited in the app)."""
    st.subheader('Load saved files')
    st.caption('Order: hardware `.py` → data path → biometrics JSON. Reload module after changing the list or typed name.')

    if st.session_state.pop('gui_tpl_reload_config', False):
        eff = effective_load_module_name(
            typed=str(st.session_state.get('gui_cfg_mod') or ''),
            selected=str(st.session_state.get('gui_load_cfg_select') or ''),
        )
        if not eff:
            st.session_state['gui_tpl_config_load_err'] = 'Choose a config module from the list or type a name to override.'
        else:
            err = _apply_loaded_config_module(eff)
            if err:
                st.session_state['gui_tpl_config_load_err'] = err
            else:
                st.rerun()
    if err_msg := st.session_state.pop('gui_tpl_config_load_err', None):
        _st_error_detail(
            'Hardware config load failed.',
            ['Check module name', 'Fix errors in Details'],
            err_msg,
        )

    st.subheader('1 · Hardware configuration')
    with st.container(border=True):
        _cfg_mods = discover_config_modules(_ROOT)
        st.session_state.setdefault('gui_cfg_mod', str(st.session_state.get('cfg_mod') or 'jimenez_bender_config_A'))
        _default_pick = st.session_state['gui_cfg_mod'] if st.session_state['gui_cfg_mod'] in _cfg_mods else _cfg_mods[0]
        st.session_state.setdefault('gui_load_cfg_select', _default_pick)
        if st.session_state.get('bender') is None and not st.session_state.get('gui_tpl_cfg_autoloaded'):
            st.session_state['gui_tpl_cfg_autoloaded'] = True
            st.session_state['gui_tpl_reload_config'] = True
            st.rerun()
        st.selectbox(
            'Hardware configuration module',
            options=_cfg_mods,
            key='gui_load_cfg_select',
            on_change=_cb_tpl_config_module_changed,
        )
        st.text_input('Or type module name (overrides list)', key='gui_cfg_mod', placeholder='e.g. my_lab_config')
        if _load_save_button(
            'Load hardware configuration',
            key='gui_btn_load_config_tpl',
            help='Reload the selected hardware `.py` module from disk (use after typing an override name).',
        ):
            st.session_state['gui_tpl_reload_config'] = True
            st.rerun()
        if st.session_state.get('bender') is not None:
            st.success(f"Hardware configuration active: `{getattr(st.session_state['bender'], 'config_name', '?')}`")

    st.divider()
    st.subheader('2 · Data file path')
    with st.container(border=True):
        _dfc, _fnc = st.columns(2)
        with _dfc:
            st.text_input(
                'Data folder',
                key='gui_data_folder',
                placeholder=r'Example: C:\Users\me\Data\Experiments',
                help=DATA_FOLDER_HELP,
            )
        with _fnc:
            st.text_input(
                'Data file name',
                key='gui_data_filename',
                placeholder='my_experiment.h5',
                help=DATA_FILE_NAME_HELP,
            )

    st.divider()
    _bio_tpl_dir = _shared_experiment_dir()
    _bio_tpl_list = list_biometrics_template_files(_bio_tpl_dir)
    _bio_opts: list = [None] + _bio_tpl_list
    st.selectbox(
        'Biometrics JSON',
        _bio_opts,
        format_func=_biometrics_template_option_label,
        key='gui_biometrics_template_select',
    )
    if _load_save_button('Load biometrics into form', key='gui_biometrics_btn_load_tpl'):
        _bp = st.session_state.get('gui_biometrics_template_select')
        if not _bp:
            st.session_state['gui_biometrics_load_feedback'] = (False, 'Choose a biometrics file first.')
        else:
            st.session_state['gui_pending_biometrics_path'] = _bp
        st.rerun()
    if bf := st.session_state.pop('gui_biometrics_load_feedback', None):
        ok_bf, txt_bf = bf
        if ok_bf:
            st.success(txt_bf)
        else:
            _st_error_detail(
                'Biometrics load failed.',
                ['Check JSON file', 'Read Details'],
                txt_bf,
            )
    if st.session_state.get('bender') is not None:
        if _load_save_button(
            'Apply loaded biometrics to experiment',
            key='gui_tpl_apply_bio_to_bender',
        ):
            _apply_all_biometrics_to_bender(st.session_state['bender'])
            st.session_state['gui_tpl_bio_done'] = True
            st.success('Biometrics applied. Procedure sections are below.')
            st.rerun()


def _template_procedure_gate() -> None:
    if _nav_route() != 'templates' or not _tpl_only_procedure():
        return
    _render_template_procedure_strip()
    if not st.session_state.get('bender'):
        st.info('Load **hardware configuration** above to continue.')
        st.stop()
    if not st.session_state.get('gui_tpl_bio_done'):
        st.info('Load a **biometrics** file and click **Apply loaded biometrics to experiment**.')
        st.stop()


def main():
    _ensure_gui_data_path_session_keys()
    if 'gui_post_notes' not in st.session_state:
        st.session_state['gui_post_notes'] = ''
    if 'gui_genus_species' not in st.session_state:
        st.session_state['gui_genus_species'] = ''
    if 'gui_specimen_id' not in st.session_state:
        st.session_state['gui_specimen_id'] = ''
    st.session_state.setdefault('gui_app_route', 'landing')
    if _nav_route() == 'landing':
        _render_landing_page()
        return

    _flush_pending_load_config_session()
    _consume_pending_biometrics_template()

    _render_app_chrome()
    _render_sidebar()

    if _nav_route() == 'templates':
        st.session_state.pop('gui_tpl_need_procedure', None)
        st.markdown('**Template mode** — check what you already have on disk.')
        st.checkbox('I have a saved hardware **config**', value=True, key='gui_tpl_chk_config')
        st.checkbox('I have a **biometrics** JSON file', value=True, key='gui_tpl_chk_biometrics')
        st.checkbox('I already have a **protocol** template (JSON)', value=False, key='gui_tpl_have_protocol_template')
        if not _tpl_only_procedure():
            st.info(
                'All sections below match **Build from scratch**, or use **Load template into form** in **section 4** when you '
                'already have a protocol file. Uncheck the third box above if you only need to configure the procedure after '
                'loading config and biometrics.'
            )

    if _nav_route() == 'stepwise':
        _render_stepwise_rail()
    else:
        st.caption('Sections follow the numbered order on the page.')

    _cfg_mods = discover_config_modules(_ROOT)
    _template_procedure_gate()

    if _show_hw_config_section():
        st.subheader('1 · Hardware configuration')
        with st.container(border=True):
            st.session_state.setdefault('gui_cfg_mod', str(st.session_state.get('cfg_mod') or 'jimenez_bender_config_A'))
            _default_pick = st.session_state['gui_cfg_mod'] if st.session_state['gui_cfg_mod'] in _cfg_mods else _cfg_mods[0]
            st.session_state.setdefault('gui_load_cfg_select', _default_pick)
            st.session_state.setdefault('gui_cfg_build_base', _cfg_mods[0])
            st.session_state.setdefault('gui_cfg_build_out', '')
            st.session_state.setdefault('gui_cfg_build_overwrite', False)

            if _template_hide_config_build_new():
                st.caption('Template mode: load-only here. Uncheck "I have a saved config" above to use **Build new**.')
                st.session_state['gui_config_setup_mode'] = 'Load existing'
                mode = 'Load existing'
            else:
                mode = _hardware_configuration_mode_toggle()

            if mode == 'Load existing':
                st.selectbox(
                    'Hardware configuration module',
                    options=_cfg_mods,
                    key='gui_load_cfg_select',
                    help='Python modules in this project folder that define hardware (DAQ, motor, channels).',
                )
                st.text_input(
                    'Override: type module name',
                    key='gui_cfg_mod',
                    help='If you fill this, it replaces the dropdown. Leave empty to use the list.',
                    placeholder='e.g. my_lab_config',
                )
            else:
                st.selectbox(
                    'Start from template (base module for import *)',
                    options=_cfg_mods,
                    key='gui_cfg_build_base',
                    help='Other settings from the template stay unless you override them below.',
                )
                _maybe_seed_cfg_build_fields()
                st.text_input(
                    'Save new config as (module name, no `.py`)',
                    key='gui_cfg_build_out',
                    placeholder='e.g. lab_setup_2026',
                    help='Writes a new `.py` file in this folder and loads it.',
                )
                with st.expander('Calibration, direction & axis labels', expanded=False):
                    st.text_input('Force/torque calibration file', key='gui_cfg_bld_forcetorque_calibration_file')
                    st.text_input('Positive motor direction (`left` / `right`)', key='gui_cfg_bld_positive_motor_direction')
                    st.number_input(
                        'Specimen lateral index on positive motor side',
                        key='gui_cfg_bld_specimen_lateral_index',
                        step=1,
                        format='%d',
                    )
                    st.text_input('Motor axis (config `motor_axis`)', key='gui_cfg_bld_motor_axis')
                    st.text_input('Bending axis — sensor (`bending_axis_sensor`)', key='gui_cfg_bld_bending_axis_sensor')
                    st.text_input('Primary bending axis (e.g. zTorque)', key='gui_cfg_bld_primary_bending_axis')
                    st.text_input(
                        'Bending axis — specimen (`bending_axis_specimen`: dorsoventral / lateral / anteroposterior)',
                        key='gui_cfg_bld_bending_axis_specimen',
                    )
                    st.text_input('S1 side label', key='gui_cfg_bld_S1side')
                    st.text_input('S2 side label', key='gui_cfg_bld_S2side')
                with st.expander('DAQ rates & device', expanded=True):
                    st.text_input('NI-DAQ device name', key='gui_cfg_bld_device_name')
                    st.number_input('AI + encoder sample rate (Hz)', key='gui_cfg_bld_daq_ai_sr', min_value=1.0, format='%.f')
                    st.number_input('AO + DO sample rate (Hz)', key='gui_cfg_bld_daq_ao_sr', min_value=1.0, format='%.f')
                with st.expander('Motor, encoder & stim channels', expanded=False):
                    st.number_input('Motor full steps per rev', key='gui_cfg_bld_motor_steps', min_value=1, step=1, format='%d')
                    st.number_input('Motor gear ratio', key='gui_cfg_bld_motor_gear', min_value=1, step=1, format='%d')
                    st.number_input('Encoder pulses per rev', key='gui_cfg_bld_encoder_ppr', min_value=1, step=1, format='%d')
                    st.text_input('Motor port', key='gui_cfg_bld_motor_port')
                    st.text_input('Encoder counter channel', key='gui_cfg_bld_encoder_chan')
                    st.text_input('Stim AO channels (comma-separated)', key='gui_cfg_bld_stim_channels')
                with st.expander('Strain / ATI channels', expanded=False):
                    st.text_input('SG AI channels (comma-separated)', key='gui_cfg_bld_SG_chan')
                    st.text_input('SG channel names (comma-separated)', key='gui_cfg_bld_SG_name')
                with st.expander('Sonomicrometry', expanded=False):
                    st.checkbox('Use sonomicrometry', key='gui_cfg_bld_use_sono')
                    st.text_input('Sono AI channels (comma-separated)', key='gui_cfg_bld_sono_channel')
                    st.text_input('Sono names (comma-separated)', key='gui_cfg_bld_sono_name')
                    st.number_input(
                        'Sono internal sample frequency (`sono_internal_samplefreq`)',
                        key='gui_cfg_bld_sono_internal_samplefreq',
                        min_value=1,
                        step=1,
                        format='%d',
                    )
                    st.text_input(
                        'Sono cal left [V_lo, V_hi, mm_lo, mm_hi] comma-separated',
                        key='gui_cfg_bld_sono_cal_left',
                        help='Four numbers, same order as the config file.',
                    )
                    st.text_input(
                        'Sono cal right [V_lo, V_hi, mm_lo, mm_hi] comma-separated',
                        key='gui_cfg_bld_sono_cal_right',
                    )
                with st.expander('Stim monitor AI (optional)', expanded=False):
                    st.text_input(
                        'Stim monitor AI channels (comma-separated; empty = none)',
                        key='gui_cfg_bld_stim_monitor_chan',
                    )
                    st.text_input('Stim monitor channel names (comma-separated)', key='gui_cfg_bld_stim_monitor_name')
                with st.expander('Default timing, ramp & motion', expanded=False):
                    st.number_input('waitbefore (s)', key='gui_cfg_bld_waitbefore', min_value=0.0, format='%.6g')
                    st.number_input('waitafter (s)', key='gui_cfg_bld_waitafter', min_value=0.0, format='%.6g')
                    st.number_input('rampdur (s)', key='gui_cfg_bld_rampdur', min_value=0.0, format='%.6g')
                    st.number_input('prepoststim_dur (s)', key='gui_cfg_bld_prepoststim_dur', min_value=0.0, format='%.6g')
                    st.number_input('prepoststim_sep (s)', key='gui_cfg_bld_prepoststim_sep', format='%.6g')
                    st.number_input('prestim_time', key='gui_cfg_bld_prestim_time', format='%.6g')
                    st.number_input('poststim_time', key='gui_cfg_bld_poststim_time', format='%.6g')
                    st.number_input('amp_step_vel', key='gui_cfg_bld_amp_step_vel', min_value=1, step=1, format='%d')
                    st.selectbox(
                        'ramp_mode_default',
                        options=['linear', 'exponential'],
                        key='gui_cfg_bld_ramp_mode_default',
                    )
                    st.caption('`units` / `unit_rules` dicts stay from the template unless you edit the generated `.py` file.')
                st.checkbox(
                    'Overwrite if a `.py` file with that name already exists',
                    key='gui_cfg_build_overwrite',
                )
                if _load_save_button('Write config file and load', key='gui_btn_write_load_config'):
                    base = str(st.session_state.get('gui_cfg_build_base') or '').strip()
                    out_raw = str(st.session_state.get('gui_cfg_build_out') or '').strip()
                    out_stem = sanitize_config_module_stem(out_raw)
                    if not base:
                        _st_error_actions('Build new blocked.', ['Pick template module above'])
                    elif not out_raw:
                        _st_error_actions('Build new blocked.', ['Enter new module file name'])
                    else:
                        stim_ch = parse_comma_list(str(st.session_state.get('gui_cfg_bld_stim_channels') or ''))
                        sg_ch = parse_comma_list(str(st.session_state.get('gui_cfg_bld_SG_chan') or ''))
                        sg_nm = parse_comma_list(str(st.session_state.get('gui_cfg_bld_SG_name') or ''))
                        sm_ch = parse_comma_list(str(st.session_state.get('gui_cfg_bld_stim_monitor_chan') or ''))
                        sm_nm = parse_comma_list(str(st.session_state.get('gui_cfg_bld_stim_monitor_name') or ''))
                        sono_cal_ok = True
                        sono_lf: list[float] = []
                        sono_rf: list[float] = []
                        try:
                            sono_lf = parse_n_floats(str(st.session_state.get('gui_cfg_bld_sono_cal_left') or ''), 4)
                            sono_rf = parse_n_floats(str(st.session_state.get('gui_cfg_bld_sono_cal_right') or ''), 4)
                        except ValueError as e:
                            sono_cal_ok = False
                            _st_error_detail(
                                'Sono calibration invalid.',
                                ['Enter four comma-separated numbers', 'Match V and mm pairs'],
                                str(e),
                            )
                        if not stim_ch:
                            _st_error_actions('Stim channels missing.', ['Add comma-separated port lines'])
                        elif not sg_ch or not sg_nm or len(sg_ch) != len(sg_nm):
                            _st_error_actions('SG lists mismatch.', ['Match SG channel name counts'])
                        elif (bool(sm_ch) ^ bool(sm_nm)) or (sm_ch and sm_nm and len(sm_ch) != len(sm_nm)):
                            _st_error_actions('Stim monitor mismatch.', ['Fill both lists or clear both'])
                        elif not sono_cal_ok:
                            pass
                        else:
                            sono_ch = parse_comma_list(str(st.session_state.get('gui_cfg_bld_sono_channel') or ''))
                            sono_nm = parse_comma_list(str(st.session_state.get('gui_cfg_bld_sono_name') or ''))
                            use_sono = bool(st.session_state.get('gui_cfg_bld_use_sono'))
                            if use_sono and (not sono_ch or not sono_nm or len(sono_ch) != len(sono_nm)):
                                _st_error_actions('Sono lists mismatch.', ['Match sono name count', 'Or disable sonomicrometry'])
                            else:
                                path = os.path.join(_ROOT, out_stem + '.py')
                                if os.path.isfile(path) and not st.session_state.get('gui_cfg_build_overwrite'):
                                    _st_error_actions('Config file exists.', ['Enable overwrite checkbox', 'Or pick new name'])
                                else:
                                    rm = str(st.session_state.get('gui_cfg_bld_ramp_mode_default') or 'linear')
                                    if rm not in ('linear', 'exponential'):
                                        rm = 'linear'
                                    assignments = {
                                        'forcetorque_calibration_file': str(
                                            st.session_state.get('gui_cfg_bld_forcetorque_calibration_file') or 'FT56491.cal'
                                        ),
                                        'positive_motor_direction': str(
                                            st.session_state.get('gui_cfg_bld_positive_motor_direction') or 'left'
                                        ),
                                        'specimen_lateral_index_on_positive_motor_side': int(
                                            st.session_state.get('gui_cfg_bld_specimen_lateral_index') or -1
                                        ),
                                        'motor_axis': str(st.session_state.get('gui_cfg_bld_motor_axis') or 'z'),
                                        'bending_axis_sensor': str(
                                            st.session_state.get('gui_cfg_bld_bending_axis_sensor') or 'z'
                                        ),
                                        'primary_bending_axis': str(
                                            st.session_state.get('gui_cfg_bld_primary_bending_axis') or 'zTorque'
                                        ),
                                        'bending_axis_specimen': str(
                                            st.session_state.get('gui_cfg_bld_bending_axis_specimen') or 'dorsoventral'
                                        ),
                                        'device_name': str(st.session_state.get('gui_cfg_bld_device_name') or 'Dev1'),
                                        'daq_ai_sample_rate_hz': float(st.session_state.get('gui_cfg_bld_daq_ai_sr') or 1000.0),
                                        'daq_ao_do_sample_rate_hz': float(
                                            st.session_state.get('gui_cfg_bld_daq_ao_sr') or 60000.0
                                        ),
                                        'motor_full_steps_per_rev': int(st.session_state.get('gui_cfg_bld_motor_steps') or 1600),
                                        'motor_gear_ratio': int(st.session_state.get('gui_cfg_bld_motor_gear') or 5),
                                        'encoder_pulses_per_rev': int(st.session_state.get('gui_cfg_bld_encoder_ppr') or 10000),
                                        'stim_channels': stim_ch,
                                        'motor_port': str(st.session_state.get('gui_cfg_bld_motor_port') or 'port0'),
                                        'encoder_chan': str(st.session_state.get('gui_cfg_bld_encoder_chan') or 'ctr0'),
                                        'SG_chan': sg_ch,
                                        'SG_name': sg_nm,
                                        'stim_monitor_chan': sm_ch,
                                        'stim_monitor_name': sm_nm,
                                        'S1side': str(st.session_state.get('gui_cfg_bld_S1side') or 'left'),
                                        'S2side': str(st.session_state.get('gui_cfg_bld_S2side') or 'right'),
                                        'use_sono': use_sono,
                                        'sono_channel': sono_ch if use_sono else [],
                                        'sono_name': sono_nm if use_sono else [],
                                        'sono_internal_samplefreq': int(
                                            st.session_state.get('gui_cfg_bld_sono_internal_samplefreq') or 241
                                        ),
                                        'sono_cal_left': sono_lf,
                                        'sono_cal_right': sono_rf,
                                        'amp_step_vel': int(st.session_state.get('gui_cfg_bld_amp_step_vel') or 10),
                                        'ramp_mode_default': rm,
                                        'waitbefore': float(st.session_state.get('gui_cfg_bld_waitbefore') or 3.0),
                                        'waitafter': float(st.session_state.get('gui_cfg_bld_waitafter') or 4.0),
                                        'rampdur': float(st.session_state.get('gui_cfg_bld_rampdur') or 0.25),
                                        'prepoststim_dur': float(st.session_state.get('gui_cfg_bld_prepoststim_dur') or 0.06),
                                        'prepoststim_sep': float(st.session_state.get('gui_cfg_bld_prepoststim_sep') or 1.0),
                                        'prestim_time': float(st.session_state.get('gui_cfg_bld_prestim_time') or -2.0),
                                        'poststim_time': float(st.session_state.get('gui_cfg_bld_poststim_time') or 2.0),
                                    }
                                    src = render_generated_config(base, assignments)
                                    with open(path, 'w', encoding='utf-8') as f:
                                        f.write(src)
                                    importlib.invalidate_caches()
                                    err = _apply_loaded_config_module(out_stem)
                                    if err:
                                        _st_error_detail(
                                            'Generated config did not load.',
                                            ['Fix assignments above', 'Read Details message'],
                                            err,
                                        )
                                    else:
                                        st.session_state['gui_cfg_mod'] = out_stem
                                        st.success(f'Wrote and loaded `{out_stem}`')
                                        st.rerun()

    if _show_hw_config_section() and _show_data_path_section():
        st.divider()
        st.caption(
            '**Data file path** (section 2, next) is separate from **hardware configuration** (section 1): the `.py` '
            'module does not set your HDF5 save location.'
        )

    if _show_data_path_section():
        st.subheader('2 · Data file path')
        st.caption('HDF5 save location (separate from the hardware `.py` module).')
        _sw_dp = _stepwise_on_data_file_path_step()
        with st.container(border=True):
            if _sw_dp:
                st.caption('Red button sets the path; it also reloads the module if you changed it in step 1.')
            st.session_state.setdefault('gui_sec1_hide', False)
            st.checkbox(
                'Collapse data path fields (values are kept)',
                key='gui_sec1_hide',
                help='Hides folder and file name inputs to reduce scrolling. Uncheck to edit again.',
            )
            if st.session_state.get('gui_sec1_hide'):
                st.caption('Expand the checkbox above to edit **Data folder** and **Data file name**.')
            else:
                df_col, fn_col = st.columns(2)
                with df_col:
                    st.text_input(
                        'Data folder',
                        key='gui_data_folder',
                        placeholder=r'Example: C:\Users\me\Data\Experiments',
                        help=DATA_FOLDER_HELP,
                    )
                with fn_col:
                    st.text_input(
                        'Data file name',
                        key='gui_data_filename',
                        placeholder='my_experiment.h5',
                        help=DATA_FILE_NAME_HELP,
                    )
                full_out = _compose_output_h5_path()
                if full_out:
                    st.caption(f'**Save path:** `{full_out}`')
                else:
                    st.caption('Enter a file name to preview the full path.')
                anchor = _output_path_anchor_for_review()
                qc_files = _candidate_review_files(anchor) if anchor else []
                if qc_files:
                    st.caption(
                        f'{len(qc_files)} data-related file(s) in this folder — use **section 8** for plots and **section 9** for notes.'
                    )
                elif full_out:
                    st.caption('No matching files in the folder yet; run or save to create some.')
            _cfg_mode = str(st.session_state.get('gui_config_setup_mode', 'Load existing'))
            if _cfg_mode == 'Load existing':
                _load_lbl = 'Set data file path' if _sw_dp else 'Load hardware configuration and data path'
                _load_hlp = (
                    'Sets save path; reloads the hardware module if needed. Does not start DAQ.'
                    if _sw_dp
                    else 'Loads the module from section 1 and sets the save path. Does not start DAQ; re-click after changes.'
                )
                if _load_save_button(_load_lbl, key='gui_btn_load_config', help=_load_hlp):
                    eff = effective_load_module_name(
                        typed=str(st.session_state.get('gui_cfg_mod') or ''),
                        selected=str(st.session_state.get('gui_load_cfg_select') or ''),
                    )
                    if not eff:
                        _st_error_actions('No config selected.', ['Pick list or type name'])
                    else:
                        b0 = st.session_state.get('bender')
                        need_hw_reload = b0 is None or not _selected_config_matches_bender(b0, eff)
                        err = _apply_loaded_config_module(eff) if need_hw_reload else None
                        if err:
                            _st_error_detail(
                                'Hardware config load failed.',
                                ['Check module name', 'Fix errors in Details'],
                                err,
                            )
                        else:
                            perr = _sec1_apply_composed_path_to_bender()
                            if perr:
                                if need_hw_reload:
                                    _st_error_detail(
                                        'Hardware OK; path failed.',
                                        ['Fix folder name', 'Fix file name'],
                                        perr,
                                    )
                                else:
                                    _st_error_detail(
                                        'Data path not applied.',
                                        ['Fix folder name', 'Fix file name'],
                                        perr,
                                    )
                            else:
                                if _sw_dp and not need_hw_reload:
                                    st.toast('Data file path set.')
                                elif _sw_dp and need_hw_reload:
                                    st.toast('Hardware configuration loaded and data file path set.')
                                else:
                                    st.toast('Hardware configuration and data file path applied.')
                                if need_hw_reload:
                                    st.success(f'Loaded `{_normalize_config_module_name(eff)}`')
                                else:
                                    st.success('Data file path set on the experiment object.')
                                st.rerun()
            if st.session_state.get('bender') is not None and _cfg_mode == 'Build new':
                _apply_lbl = 'Set data file path' if _sw_dp else 'Apply data file path'
                _apply_hlp = (
                    'Copy **Data folder** and **Data file name** to the experiment save path (after **Write config file and load**).'
                    if _sw_dp
                    else (
                        'Copy **Data folder** and **Data file name** to the experiment save path. Use after **Write config file and load** '
                        'when you change where files should be saved.'
                    )
                )
                if _load_save_button(_apply_lbl, key='gui_sec1_apply_paths', help=_apply_hlp):
                    perr = _sec1_apply_composed_path_to_bender()
                    if perr:
                        _st_error_detail(
                            'Data path not applied.',
                            ['Fix folder name', 'Fix file name'],
                            perr,
                        )
                    else:
                        st.toast('Data file path set.' if _sw_dp else 'Data file path applied.')
                        st.rerun()

    if 'bender' not in st.session_state:
        st.stop()

    b: Bender = st.session_state['bender']
    _init_biometrics_session_state(b, force=False)
    _sync_biometric_flags_from_session(b)
    _ensure_review_file_selection(
        _candidate_review_files(_output_path_anchor_for_review(b)) if _output_path_anchor_for_review(b) else []
    )
    schema = b.get_dispatch_schema()
    test_types = list(schema['test_types'])
    _consume_pending_protocol_template(test_types)

    if _show_full_sec2():
        st.subheader('3 · Biometrics')

        st.markdown('**Biometrics templates**')
        st.caption(
            'Save or reload this section as JSON in your **Data folder** (**section 2**). After **Load biometrics**, use **Apply** '
            'in each block below (or **Apply all biometrics**) to update the experiment object.'
        )
        _df_check = str(st.session_state.get('gui_data_folder') or '').strip()
        _bio_tpl_dir = _shared_experiment_dir()
        if not (_df_check and os.path.isdir(os.path.normpath(_df_check))):
            st.caption(
                f'**Data folder** is not set or not found on disk—listing template JSON files from `{_bio_tpl_dir}` until '
                '**section 2** points to a valid folder.'
            )
        if bf := st.session_state.pop('gui_biometrics_load_feedback', None):
            ok_bf, txt_bf = bf
            if ok_bf:
                st.success(txt_bf)
            else:
                _st_error_detail(
                    'Biometrics load failed.',
                    ['Check JSON file', 'Read Details'],
                    txt_bf,
                )
        if bfs := st.session_state.pop('gui_biometrics_save_feedback', None):
            ok_s, txt_s = bfs
            if ok_s:
                st.success(txt_s)
            else:
                _st_error_detail(
                    'Biometrics save failed.',
                    ['Check name and folder', 'Read Details'],
                    txt_s,
                )
        _bio_tpl_list = list_biometrics_template_files(_bio_tpl_dir)
        _bio_opts: list = [None] + _bio_tpl_list
        _bio_pick = st.selectbox(
            'Biometrics file to load',
            _bio_opts,
            format_func=_biometrics_template_option_label,
            key='gui_biometrics_template_select',
        )
        c_bl, c_bs = st.columns(2)
        with c_bl:
            if _load_save_button(
                'Load biometrics into form',
                key='gui_biometrics_btn_load',
                help='Fills biometrics widgets from the file. Use **Apply** in each block below (or **Apply all**) to update `Bender`.',
            ):
                if not _bio_pick:
                    st.session_state['gui_biometrics_load_feedback'] = (False, 'Choose a biometrics file first.')
                else:
                    st.session_state['gui_pending_biometrics_path'] = _bio_pick
                st.rerun()
        with c_bs:
            st.caption('Independent of **protocol templates** in **section 4**.')

        st.text_input('Save biometrics as (name)', key='gui_biometrics_new_name', placeholder='e.g. Zebrafish adult default')
        st.text_area('Description (optional)', key='gui_biometrics_new_desc', height=50, placeholder='Optional note stored in the JSON file.')
        st.checkbox('Overwrite if same file name exists', key='gui_biometrics_overwrite')
        if _load_save_button('Save biometrics', key='gui_biometrics_btn_save'):
            _bn = str(st.session_state.get('gui_biometrics_new_name') or '').strip()
            _bd = str(st.session_state.get('gui_biometrics_new_desc') or '').strip()
            _bst = sanitize_biometrics_filename_stem(_bn or 'biometrics')
            _bout = os.path.normpath(os.path.join(_shared_experiment_dir(), f'{_bst}.json'))
            try:
                if os.path.isfile(_bout) and not bool(st.session_state.get('gui_biometrics_overwrite')):
                    st.session_state['gui_biometrics_save_feedback'] = (
                        False,
                        f'File already exists: `{_bout}`. Enable **Overwrite** or change the name.',
                    )
                else:
                    os.makedirs(os.path.dirname(_bout) or '.', exist_ok=True)
                    save_biometrics_template(
                        _bout,
                        name=_bn or _bst,
                        description=_bd,
                        session_state=st.session_state,
                    )
                    st.session_state['gui_biometrics_save_feedback'] = (True, f'Saved `{_bout}`')
            except Exception as e:
                st.session_state['gui_biometrics_save_feedback'] = (False, f'{type(e).__name__}: {e}')
            st.rerun()

        st.divider()
        st.markdown('**Specimen identity**')
        st.caption('Export metadata (`genus_species`, `specimen_id`) and notebook-style `fishcode` (mirrors specimen ID).')
        id1, id2 = st.columns(2)
        with id1:
            st.text_input(
                'Genus-species',
                key='gui_genus_species',
                placeholder='e.g. Danio rerio',
                help='Stored in the exported `.h5` under protocol metadata (`genus_species`) when you run or export.',
            )
        with id2:
            st.text_input(
                'Specimen ID',
                key='gui_specimen_id',
                placeholder='e.g. fish-042 or prep code',
                help='Primary specimen label; also written to `fishcode` on the experiment object for notebook compatibility.',
            )
        st.text_input('Segment / preparation label (`segment`)', key='bio_segment', placeholder='e.g. whole body, hemi')

        st.divider()
        st.session_state.setdefault('gui_bio_hide', False)
        if st.session_state.get('gui_bio_hide'):
            st.caption(
                'Section body hidden. Uncheck **Hide section** at the bottom to edit biometrics inputs.'
            )
        else:
            st.markdown('### Intrinsic biometrics')
            st.caption(
                'Whole-specimen lengths, mass, material density, environment, and offsets. Identity fields above are applied '
                'with this block.'
            )
            L1, L2 = st.columns(2)
            with L1:
                st.number_input('Total length TL (`fishlen_TL`)', min_value=0.0, format='%.6g', key='bio_fishlen_TL')
                st.number_input('Total length SL (`fishlen_SL`)', min_value=0.0, format='%.6g', key='bio_fishlen_SL')
            with L2:
                st.number_input('Vertical offset `dvert` (mm)', min_value=0.0, format='%.6g', key='bio_dvert')
                st.number_input('Horizontal offset `dhoriz` (mm)', min_value=0.0, format='%.6g', key='bio_dhoriz')
            m1, m2 = st.columns(2)
            with m1:
                st.number_input('Mass `fishmass` (g)', min_value=0.0, format='%.6g', key='bio_fishmass')
            with m2:
                st.number_input('Specimen density (g / mm³)', min_value=1e-9, format='%.6g', key='bio_prof_rho')
            t1, t2 = st.columns(2)
            with t1:
                st.number_input('temp_C_room', min_value=-5.0, max_value=60.0, format='%.3f', key='bio_temp_room')
            with t2:
                st.number_input('temp_C_tank', min_value=-5.0, max_value=60.0, format='%.3f', key='bio_temp_tank')
            if st.button('Apply intrinsic biometrics', key='bio_btn_intrinsic'):
                _apply_intrinsic_biometrics_to_bender(b)
                st.toast('Intrinsic biometrics applied.')

            st.divider()
            st.markdown('### Clamp geometry')
            st.caption('Spacing between clamps and where the bend is measured along the body; cross-section at the test segment for strain ↔ motor mapping.')
            st.number_input(
                'Test segment length = clamp spacing (`dclamp` / `test_segment_length_mm`)',
                min_value=0.001,
                format='%.6g',
                key='bio_dclamp',
            )
            st.number_input(
                'Test segment position (`dbend` / `test_segment_position_mm`)',
                min_value=0.0,
                format='%.6g',
                key='bio_dbend',
            )
            x1, x2 = st.columns(2)
            with x1:
                st.number_input('Width `xsec_width` (mm)', min_value=0.001, format='%.6g', key='bio_xsec')
            with x2:
                st.number_input('Height `xsec_height` (mm)', min_value=0.001, format='%.6g', key='bio_xsec_height')
            if st.button('Apply clamp geometry', key='bio_btn_clamp'):
                _apply_clamp_geometry_to_bender(b)
                st.toast('Clamp geometry applied.')

            st.divider()
            st.markdown('### Mounted body profile (inertial model)')
            st.caption(
                'Tapered outline between proximal and distal cross-sections, outline length, and integration settings. '
                'Uses the same **density** as intrinsic biometrics. Drives rotational inertia for corrections.'
            )
            st.number_input(
                'Specimen outline length for profile model (mm)',
                min_value=0.001,
                format='%.6g',
                key='bio_prof_L',
                help='Length along the specimen used with the tapered outline below (`specimen_profile_length_mm`).',
            )
            p1, p2 = st.columns(2)
            with p1:
                st.number_input('Proximal height (mm)', min_value=0.001, format='%.6g', key='bio_prof_ph')
                st.number_input('Proximal width (mm)', min_value=0.001, format='%.6g', key='bio_prof_pw')
            with p2:
                st.number_input('Distal height (mm)', min_value=0.001, format='%.6g', key='bio_prof_dh')
                st.number_input('Distal width (mm)', min_value=0.001, format='%.6g', key='bio_prof_dw')
            p3, p4 = st.columns(2)
            with p3:
                st.number_input('Clamp offset (mm)', min_value=0.0, format='%.6g', key='bio_prof_clamp')
            with p4:
                st.number_input('Profile integration samples', min_value=20, max_value=400, step=10, key='bio_prof_samples')
            if st.button('Apply mounted profile & inertia model', key='bio_btn_profile'):
                _apply_mounted_profile_inertial_to_bender(b)
                st.toast('Profile / inertia model applied.')

            st.divider()
            st.checkbox(
                'Check here to perform inertial correction',
                key='bio_use_theoretical_inertial',
                help=(
                    'Subtracts model **system** inertia (from calibration, if loaded) and **specimen** inertia from the '
                    'profile above when correcting measured torque.'
                ),
            )

        if st.button(
            'Apply all biometrics',
            key='bio_btn_apply_all',
            help='Copies intrinsic biometrics, clamp geometry, mounted profile / inertia model, and the inertial-correction flag onto the experiment object.',
        ):
            _apply_all_biometrics_to_bender(b)
            st.toast('Biometrics applied.')

        st.checkbox(
            'Hide section (values stay; unhide to edit)',
            key='gui_bio_hide',
            help='Collapse biometrics fields after you are done editing or applying.',
        )

    if _show_sec3_through_6():

        st.divider()
        st.subheader('4 · Experiment type & parameters')

        st.markdown('**Protocol templates (load)**')
        st.caption(
            'Lists `.json` files in your **Data folder** (**section 2**). **Load template into form** sets **experiment type** '
            'and procedure fields; click **Apply** (below or in **section 6**) to copy onto the Bender object. Does not '
            'change **section 3** biometrics—use **Load biometrics** there or enter fields manually.'
        )
        if fb := st.session_state.pop('gui_protocol_load_feedback', None):
            ok_fb, txt_fb = fb
            if ok_fb:
                st.success(txt_fb)
            else:
                _st_error_detail(
                    'Protocol load failed.',
                    ['Check JSON template', 'Read Details'],
                    txt_fb,
                )
        _tpl_folder_top = _shared_experiment_dir()
        _tpl_files_top = list_template_files(_tpl_folder_top)
        _tpl_options_top: list = [None] + _tpl_files_top
        _tpl_pick_top = st.selectbox(
            'Template to load',
            _tpl_options_top,
            format_func=_protocol_template_option_label,
            key='gui_protocol_template_select',
        )
        c_tl_top, c_ts_top = st.columns(2)
        with c_tl_top:
            if _load_save_button(
                'Load template into form',
                key='gui_protocol_btn_load',
                help='Sets **Experiment type** and procedure widgets from the file. Then click **Apply** to copy onto the Bender object.',
            ):
                if not _tpl_pick_top:
                    st.session_state['gui_protocol_load_feedback'] = (False, 'Choose a template file first.')
                else:
                    st.session_state['gui_pending_protocol_template_path'] = _tpl_pick_top
                st.rerun()
        with c_ts_top:
            st.caption('Same directory as your **Data folder** and **Save** in **Procedure fields**.')

        st.divider()

        tt = st.selectbox('Experiment type (test_type)', test_types, key='test_type_select')
        b.test_type = tt

        st.session_state.setdefault('gui_exp_hide', False)
        st.caption(
            'Fields are stored in the browser session. Use **Apply** below (or in **section 6**) to copy them onto the '
            'experiment object before **Refresh preview** or **Run**. Collapsing only hides this panel; widgets keep '
            'your entries. If anything looks out of sync, use **Load hardware configuration and data path** again '
            '(**sections 1–2**, load existing).'
        )

        updates = {}

        with st.expander('Procedure fields', expanded=not bool(st.session_state.get('gui_exp_hide'))):
            if tt == 'isometric':
                st.caption(
                    '**Isometric** turns strain or curvature targets into motor angles using **test segment length** '
                    'and **cross-section width** from **section 3** (same as clamp spacing `dclamp`). '
                    'Those values are copied to the experiment object when you **Run** or **Apply**.'
                )
                st.markdown('**Required**')
                for key in schema['isometric_required']:
                    label = key.replace('_', ' ')
                    updates[key] = _render_field(
                        b,
                        key,
                        'float' if 'steps' not in key else 'int',
                        label,
                        help_text=ISOMETRIC_FIELD_HELP.get(key),
                    )
                if 'isometric_num_steps' in updates and updates['isometric_num_steps'] is not None:
                    updates['isometric_num_steps'] = int(updates['isometric_num_steps'])
                st.markdown('**Optional**')
                for key in schema['isometric_optional']:
                    if key == 'isometric_stim_params':
                        updates[key] = _render_field(
                            b,
                            key,
                            'json_dict',
                            'Isometric step timing & stimulation (JSON, optional)',
                            help_text=ISOMETRIC_STIM_JSON_HELP,
                        )
                    elif key == 'isometric_stim_overrides':
                        updates[key] = _render_field(
                            b,
                            key,
                            'json_dict',
                            'Isometric stim routing overrides (JSON, advanced)',
                            help_text=ISOMETRIC_STIM_OVERRIDES_HELP,
                        )
                    elif key == 'recruitment':
                        skr = _widget_key('recruitment')
                        cur_r = _get_session_value(b, 'recruitment', 'bilateral_simultaneous')
                        if skr not in st.session_state:
                            st.session_state[skr] = cur_r if cur_r in RECRUITMENT_OPTIONS else RECRUITMENT_OPTIONS[0]
                        updates[key] = st.selectbox(
                            'recruitment',
                            list(RECRUITMENT_OPTIONS),
                            key=skr,
                            help=RECRUITMENT_FIELD_HELP,
                        )
                    elif key == 'lateral_mode':
                        skl = _widget_key('lateral_mode')
                        if skl not in st.session_state:
                            st.session_state[skl] = str(_get_session_value(b, key) or '')
                        updates[key] = st.text_input(LATERAL_MODE_LABEL, key=skl, help=ISOMETRIC_FIELD_HELP.get('lateral_mode'))
                    elif key in ('bilateral_mirror_motor',):
                        updates[key] = _render_field(
                            b, key, 'bool', BILATERAL_MIRROR_LABEL, help_text=ISOMETRIC_FIELD_HELP.get(key)
                        )
                    elif key == 'bilateral_sequential_left_frac':
                        updates[key] = _render_field(
                            b, key, 'float', key, help_text=ISOMETRIC_FIELD_HELP.get(key)
                        )
                    elif key == 'isometric_mode':
                        modes = list(ALL_AMPS_MODE_OPTIONS)
                        skm = _widget_key('isometric_mode')
                        cur_m = str(_get_session_value(b, key, 'strain'))
                        if skm not in st.session_state:
                            st.session_state[skm] = cur_m if cur_m in modes else 'strain'
                        updates[key] = st.selectbox(
                            'isometric_mode (how to interpret isometric_initial / isometric_final)',
                            modes,
                            key=skm,
                            format_func=_format_strain_or_amp_mode,
                            help=ISOMETRIC_FIELD_HELP.get(key),
                        )
                    elif key == 'isometric_inter_step_interval_s':
                        skg = _widget_key('isometric_inter_step_interval_s')
                        if skg not in st.session_state:
                            v0 = _get_session_value(b, key, 0.0)
                            st.session_state[skg] = float(v0) if v0 is not None else 0.0
                        updates[key] = float(
                            st.number_input(
                                'Time between steps (s)',
                                min_value=0.0,
                                format='%.6g',
                                key=skg,
                                help=(
                                    'Seconds to wait after each step finishes (after acquisition) before the next '
                                    'ramp/stimulus period. Use **0** for back-to-back steps.'
                                ),
                            )
                        )
                    elif 'random_seed' in key:
                        sks = _widget_key(key)
                        if sks not in st.session_state:
                            v0 = _get_session_value(b, key)
                            st.session_state[sks] = '' if v0 is None else str(v0)
                        s = st.text_input('Random seed (optional)', key=sks, help=RANDOM_SEED_HELP)
                        if not str(s).strip():
                            updates[key] = None
                        else:
                            try:
                                updates[key] = int(s)
                            except ValueError:
                                _st_error_actions(
                                    'Random seed invalid.',
                                    ['Use whole number only', 'Or leave field blank'],
                                )
                                updates[key] = None
                    else:
                        kind = 'bool' if 'randomize' in key else 'str'
                        updates[key] = _render_field(
                            b, key, kind, key, help_text=ISOMETRIC_FIELD_HELP.get(key)
                        )

            elif tt == 'isovelocity':
                st.markdown('**Required**')
                for key in schema['isovelocity_required']:
                    kind = 'int' if 'num_steps' in key else 'float'
                    lbl = ISOVELOCITY_WIDGET_LABEL.get(key, key.replace('_', ' '))
                    updates[key] = _render_field(
                        b, key, kind, lbl, help_text=ISOVELOCITY_FIELD_HELP.get(key)
                    )
                if 'isovelocity_num_steps' in updates and updates['isovelocity_num_steps'] is not None:
                    updates['isovelocity_num_steps'] = int(updates['isovelocity_num_steps'])
                st.markdown('**Optional**')
                for key in schema['isovelocity_optional']:
                    if key == 'isovelocity_stim_params':
                        updates[key] = _render_field(
                            b,
                            key,
                            'json_dict',
                            'Isovelocity segment timing & stimulation (JSON, optional)',
                            help_text=ISOVELOCITY_STIM_JSON_HELP,
                        )
                    elif key == 'isovelocity_stim_overrides':
                        updates[key] = _render_field(
                            b,
                            key,
                            'json_dict',
                            'Isovelocity stim routing overrides (JSON, advanced)',
                            help_text=ISOVELOCITY_STIM_OVERRIDES_HELP,
                        )
                    elif key == 'recruitment':
                        skr = _widget_key('recruitment')
                        cur_r = _get_session_value(b, 'recruitment', 'bilateral_simultaneous')
                        if skr not in st.session_state:
                            st.session_state[skr] = cur_r if cur_r in RECRUITMENT_OPTIONS else RECRUITMENT_OPTIONS[0]
                        updates[key] = st.selectbox(
                            'recruitment',
                            list(RECRUITMENT_OPTIONS),
                            key=skr,
                            help=RECRUITMENT_FIELD_HELP,
                        )
                    elif key == 'lateral_mode':
                        skl = _widget_key('lateral_mode')
                        if skl not in st.session_state:
                            st.session_state[skl] = str(_get_session_value(b, key) or '')
                        updates[key] = st.text_input(
                            LATERAL_MODE_LABEL, key=skl, help=ISOVELOCITY_FIELD_HELP.get('lateral_mode')
                        )
                    elif key in ('bilateral_mirror_motor',):
                        updates[key] = _render_field(
                            b, key, 'bool', BILATERAL_MIRROR_LABEL, help_text=ISOVELOCITY_FIELD_HELP.get(key)
                        )
                    elif key == 'bilateral_sequential_left_frac':
                        updates[key] = _render_field(
                            b, key, 'float', key, help_text=ISOVELOCITY_FIELD_HELP.get(key)
                        )
                    elif key == 'isovelocity_starting_strain_mode':
                        modes = list(ALL_AMPS_MODE_OPTIONS)
                        skm = _widget_key('isovelocity_starting_strain_mode')
                        cur_m = str(_get_session_value(b, key, 'strain'))
                        if skm not in st.session_state:
                            st.session_state[skm] = cur_m if cur_m in modes else 'strain'
                        updates[key] = st.selectbox(
                            ISOVELOCITY_WIDGET_LABEL.get(
                                key, 'isovelocity_starting_strain_mode (unit for isovelocity_starting_strain)'
                            ),
                            modes,
                            key=skm,
                            format_func=_format_strain_or_amp_mode,
                            help=ISOVELOCITY_FIELD_HELP.get(key),
                        )
                    elif 'random_seed' in key:
                        sks = _widget_key(key)
                        if sks not in st.session_state:
                            v0 = _get_session_value(b, key)
                            st.session_state[sks] = '' if v0 is None else str(v0)
                        s = st.text_input('Random seed (optional)', key=sks, help=RANDOM_SEED_HELP)
                        if not str(s).strip():
                            updates[key] = None
                        else:
                            try:
                                updates[key] = int(s)
                            except ValueError:
                                _st_error_actions(
                                    'Random seed invalid.',
                                    ['Use whole number only', 'Or leave field blank'],
                                )
                                updates[key] = None
                    else:
                        kind = 'bool' if 'randomize' in key else 'float'
                        lbl = ISOVELOCITY_WIDGET_LABEL.get(key, key.replace('_', ' '))
                        updates[key] = _render_field(
                            b, key, kind, lbl, help_text=ISOVELOCITY_FIELD_HELP.get(key)
                        )

            elif tt == 'calibration':
                st.markdown('**Required**')
                bases = [x for x in test_types if x != 'calibration']
                cur = _get_session_value(b, 'calibration_base_test_type', 'dynamic')
                if cur not in bases:
                    cur = 'dynamic'
                sk_cal = _widget_key('calibration_base_test_type')
                if sk_cal not in st.session_state:
                    st.session_state[sk_cal] = cur
                if st.session_state[sk_cal] not in bases:
                    st.session_state[sk_cal] = bases[0]
                updates['calibration_base_test_type'] = st.selectbox(
                    'calibration_base_test_type', bases, key=sk_cal
                )
                st.info(
                    'Calibration runs the **base** motion protocol. Set **test_type** to that base '
                    '(e.g. dynamic), click **Apply** in **section 4** or **6** (and **Refresh preview** if you use it), '
                    'then switch back to **calibration** before running.'
                )
                st.markdown('**Optional**')
                for key in schema['calibration_optional']:
                    sko = _widget_key(key)
                    if sko not in st.session_state:
                        st.session_state[sko] = str(_get_session_value(b, key) or '')
                    updates[key] = st.text_input(key, key=sko)

            elif tt in MOTION_TYPES:
                st.markdown('**Motion-series parameters** (procedure-specific)')
                if tt == 'dynamic':
                    st.caption(
                        'Timeline length follows **cycles** (freq × cycles_per_step + end cycles), not **Duration**. '
                        'If `period_by_cycle` is not set yet, call **organize_cycles** from a notebook/script after **Apply**.'
                    )
                elif tt == 'step_change':
                    st.caption('Motion length is computed inside **make_cycles_step_change**; **Duration** does not apply.')
                fields = _motion_parameter_rows(tt)
                for name, kind, label in fields:
                    updates[name] = _render_field(
                        b, name, kind, label, help_text=MOTION_FIELD_HELP.get(name)
                    )

            else:
                st.warning(f'No dedicated field panel for {tt!r} yet; use notebook or extend this script.')

            st.divider()
            st.markdown('**Save current procedure as template**')
            st.caption(
                '**Save** stores the current **experiment type** and all **procedure fields** for that type '
                '(dynamic / sweeps / steps / isometric / isovelocity / calibration). **Section 2** (data path) and '
                '**section 3** biometrics are not included—use **Save biometrics** there. For **calibration**, the template '
                'can also embed the **base protocol** from the Bender object (e.g. frequencies & strains) if you **Apply** '
                'that base type before saving. Use **Protocol templates (load)** above, then **Apply** and **Run**.'
            )
            if sf := st.session_state.pop('gui_protocol_save_feedback', None):
                ok_sf, txt_sf = sf
                if ok_sf:
                    st.success(txt_sf)
                else:
                    _st_error_detail(
                        'Protocol save failed.',
                        ['Check name and folder', 'Read Details'],
                        txt_sf,
                    )
            st.text_input(
                'Template name',
                key='gui_protocol_new_name',
                placeholder='e.g. Protocol A (any test_type)',
            )
            st.text_area(
                'Description (optional)',
                key='gui_protocol_new_desc',
                height=70,
                placeholder='e.g. Isometric 5 steps; or dynamic 1/3/5 Hz × strains (strain_pct); or calibration + base',
            )
            st.checkbox('Overwrite if a file with the same name already exists', key='gui_protocol_overwrite')
            if _load_save_button('Save template', key='gui_protocol_btn_save'):
                _name = str(st.session_state.get('gui_protocol_new_name') or '').strip()
                _desc = str(st.session_state.get('gui_protocol_new_desc') or '').strip()
                _stem = sanitize_template_filename_stem(_name or 'protocol')
                _out = os.path.normpath(os.path.join(_shared_experiment_dir(), f'{_stem}.json'))
                try:
                    if os.path.isfile(_out) and not bool(st.session_state.get('gui_protocol_overwrite')):
                        st.session_state['gui_protocol_save_feedback'] = (
                            False,
                            f'File already exists: `{_out}`. Enable **Overwrite** or pick a different name.',
                        )
                    else:
                        os.makedirs(os.path.dirname(_out) or '.', exist_ok=True)
                        _proc = build_procedure_dict_from_updates(updates)
                        _base = None
                        if tt == 'calibration':
                            _btt = _proc.get('calibration_base_test_type')
                            if _btt:
                                _snap = snapshot_bender_procedure(b, schema, str(_btt))
                                if _snap:
                                    _base = {'test_type': str(_btt), 'procedure': _snap}
                        save_protocol_template(
                            _out,
                            name=_name or _stem,
                            description=_desc,
                            test_type=tt,
                            procedure=_proc,
                            base_protocol=_base,
                        )
                        st.session_state['gui_protocol_save_feedback'] = (True, f'Saved `{_out}`')
                except Exception as e:
                    st.session_state['gui_protocol_save_feedback'] = (False, f'{type(e).__name__}: {e}')
                st.rerun()

        ap1, ap2 = st.columns([1, 4])
        with ap1:
            if st.button(
                'Apply procedure',
                key='gui_exp_apply',
                help='Copy procedure fields onto the experiment object (not **Run experiment**).',
            ):
                _apply_procedure_form_to_bender(b, updates, tt)
                st.toast('Settings applied.')
        with ap2:
            st.caption('Same as **Apply procedure** in section 6.')

        st.checkbox(
            'Hide section (values stay; unhide to edit)',
            key='gui_exp_hide',
            help='Collapse **Procedure fields** after you finish editing, saving a template, or clicking **Apply**.',
        )

        st.divider()
        st.session_state.setdefault('gui_sec4_hide', False)
        st.subheader('5 · Experiment preview (table & plot, no DAQ)')
        st.caption(
            'Uses the same motion math as a real run for **commanded** angle / velocity. '
            '**Refresh preview** copies the form onto the experiment object, then recomputes commands. '
            'For **dynamic**, preview calls **organize_cycles** and updates `period_by_cycle`, so a following **Run** '
            'matches the preview if you do not overwrite those arrays elsewhere. '
            'Set **test_segment_length_mm** and **xsec_width** in **section 3** (or **Apply** there) so preview matches strain geometry.'
        )
        if st.session_state.get('gui_sec4_hide'):
            st.caption('Preview panel hidden. Uncheck **Hide section** below to show controls and plots.')
        else:
            c_ap4, _ = st.columns([1, 4])
            with c_ap4:
                if st.button(
                    'Apply (procedure + biometrics)',
                    key='gui_preview_apply',
                    help='Copy procedure and biometrics onto `Bender` (like section 6 **Apply procedure** plus biometrics).',
                ):
                    _sync_biometric_flags_from_session(b)
                    _sync_genus_species_to_bender(b)
                    _apply_procedure_form_to_bender(b, updates, tt)
                    st.toast('Settings applied.')
            pv_on = st.checkbox('Show last preview', value=True, key='gui_show_preview')
            pv_pts = st.slider('Preview plot resolution (max points)', 400, 12000, 6000, step=200, key='gui_preview_pts')
            if st.button('Refresh preview'):
                _sync_biometric_flags_from_session(b)
                _sync_genus_species_to_bender(b)
                _apply_form_updates(b, updates, tt)
                st.session_state['gui_last_preview'] = build_protocol_preview(
                    b, requested_test_type=tt, max_plot_points=int(pv_pts)
                )
                st.session_state['gui_last_preview_tt'] = tt

            if pv_on and st.session_state.get('gui_last_preview') is not None:
                prev = st.session_state['gui_last_preview']
                if st.session_state.get('gui_last_preview_tt') != tt:
                    st.warning('Test type changed since the last preview — click **Refresh preview** to update.')
                if prev.get('error'):
                    _st_error_actions(
                        'Preview failed.',
                        ['Fix highlighted parameters', 'Click Refresh preview'],
                    )
                    with st.expander('Preview error detail'):
                        st.code(str(prev['error']))
                elif prev.get('ok'):
                    if tt == 'calibration':
                        st.info(
                            f"Calibration uses base protocol **{prev.get('motion_test_type')}** for motion "
                            '(same as **Run experiment**).'
                        )
                    if prev.get('table'):
                        st.markdown('**Summary table**')
                        st.dataframe(pd.DataFrame(prev['table']), use_container_width=True, hide_index=True)
                    tp = prev.get('t_plot')
                    ap = prev.get('angle_plot')
                    vp = prev.get('anglevel_plot')
                    if tp is not None and ap is not None and len(tp) > 0:
                        if tt == 'isometric' and prev.get('preview_isometric'):
                            st.markdown('**Command preview** (isometric: ramp–hold per step, same timing as run)')
                        elif tt == 'isovelocity' and prev.get('preview_isovelocity'):
                            st.markdown(
                                '**Command preview** (isovelocity: pre-hold at start angle + constant ω segment per step, '
                                'same timing as run)'
                            )
                        else:
                            st.markdown('**Command preview** (motor angle and angular velocity)')
                        fig = go.Figure()
                        fig.add_trace(go.Scatter(x=tp, y=ap, mode='lines', name='Commanded angle (deg)'))
                        if vp is not None and len(vp) > 0:
                            fig.add_trace(
                                go.Scatter(
                                    x=tp,
                                    y=vp,
                                    mode='lines',
                                    name='Commanded anglevel (deg/s)',
                                    yaxis='y2',
                                )
                            )
                        if vp is not None and len(vp) > 0:
                            fig.update_layout(
                                height=420,
                                margin=dict(l=48, r=48, t=40, b=40),
                                xaxis_title='Time (s)',
                                yaxis=dict(title='Angle (deg)'),
                                yaxis2=dict(title='Anglevel (deg/s)', overlaying='y', side='right'),
                                legend=dict(yanchor='top', y=0.99, xanchor='left', x=0.01),
                            )
                        else:
                            fig.update_layout(
                                height=420,
                                margin=dict(l=48, r=48, t=40, b=40),
                                xaxis_title='Time (s)',
                                yaxis=dict(title='Angle (deg)'),
                                legend=dict(yanchor='top', y=0.99, xanchor='left', x=0.01),
                            )
                        st.plotly_chart(fig, use_container_width=True)
                    elif prev.get('table') and tt in ('isometric', 'isovelocity'):
                        st.caption(
                            'Step protocols: table lists setpoints; refresh preview after fixing errors to see the plot.'
                        )
                else:
                    st.warning('Preview incomplete; click **Refresh preview** again.')

        st.checkbox(
            'Hide section (values stay; unhide to edit)',
            key='gui_sec4_hide',
            help='Collapse preview controls and plots when you are done reviewing.',
        )

        st.divider()
        st.session_state.setdefault('gui_sec5_hide', False)
        st.subheader('6 · Save, validate, and run')
        if st.session_state.get('gui_sec5_hide'):
            st.caption('Run controls hidden. Uncheck **Hide section** below.')
        else:
            if st.button('View current settings'):
                _sync_biometric_flags_from_session(b)
                _sync_genus_species_to_bender(b)
                _apply_form_updates(b, updates, tt)
                settings_rows = [
                    {'group': 'experiment', 'name': 'test_type', 'value': tt},
                    {
                        'group': 'export',
                        'name': 'data_file_target_h5',
                        'value': _compose_output_h5_path() or getattr(b, 'outputfile', None),
                    },
                    {
                        'group': 'specimen',
                        'name': 'genus_species',
                        'value': (getattr(b, 'h5_protocol_metadata', {}) or {}).get('genus_species', ''),
                    },
                    {
                        'group': 'specimen',
                        'name': 'specimen_id',
                        'value': (getattr(b, 'h5_protocol_metadata', {}) or {}).get('specimen_id', ''),
                    },
                    {'group': 'biometric', 'name': 'test_segment_length_mm', 'value': getattr(b, 'dclamp', None)},
                    {'group': 'biometric', 'name': 'test_segment_position_mm', 'value': getattr(b, 'dbend', None)},
                    {'group': 'biometric', 'name': 'xsec_width', 'value': getattr(b, 'xsec_width', None)},
                    {'group': 'biometric', 'name': 'temp_C_room', 'value': getattr(b, 'temp_C_room', None)},
                    {'group': 'biometric', 'name': 'temp_C_tank', 'value': getattr(b, 'temp_C_tank', None)},
                ]
                for k, v in sorted(updates.items(), key=lambda kv: kv[0]):
                    settings_rows.append({'group': 'parameter', 'name': k, 'value': str(v)})
                st.dataframe(pd.DataFrame(settings_rows), use_container_width=True, hide_index=True)
            c1, c2 = st.columns(2)
            with c1:
                if st.button(
                    'Apply procedure',
                    help='Copy procedure fields to the experiment object (no DAQ).',
                ):
                    _apply_procedure_form_to_bender(b, updates, tt)
                    st.toast('Settings applied.')
            with c2:
                if st.button(
                    'Check required fields',
                    help='Runs validation for the current test type (does not talk to hardware).',
                ):
                    _sync_biometric_flags_from_session(b)
                    _sync_genus_species_to_bender(b)
                    _apply_form_updates(b, updates, tt)
                    rep = b.validate_dispatch_setup(test_type=tt)
                    if rep['ok']:
                        st.success('All required fields for this procedure are set.')
                    else:
                        _st_error_actions(
                            'Required fields missing.',
                            _missing_fields_to_actions(list(rep['missing'])),
                        )
            daq_ok = st.checkbox('Hardware: I intend to run DAQ', value=False)
            bio_confirm = st.checkbox(
                'Biometrics applied (section 3)',
                value=False,
                key='gui_run_biometrics_confirm',
                help='Confirm you loaded JSON or used Apply / Apply all in section 3.',
            )
            needs_cal_confirm = _needs_missing_calibration_confirmation(b)
            if needs_cal_confirm:
                st.warning('No calibration file detected. Are you sure you wish to proceed?')
            ok_wo_cal = st.checkbox(
                'Yes, proceed without calibration file',
                key='gui_confirm_run_without_calibration',
                disabled=not needs_cal_confirm,
            )
            _, _run_big, _ = st.columns([1, 2, 1])
            with _run_big:
                if st.button(
                    'Run experiment',
                    type='primary',
                    use_container_width=True,
                    disabled=not daq_ok or not bio_confirm,
                    help='Starts DAQ acquisition.',
                ):
                    _sync_biometric_flags_from_session(b)
                    _apply_form_updates(b, updates, tt)
                    _sync_genus_species_to_bender(b)
                    outp = _compose_output_h5_path().strip()
                    if outp:
                        b.outputfile = outp
                    notes_in = str(st.session_state.get('gui_post_notes') or '').strip()
                    if needs_cal_confirm and not ok_wo_cal:
                        st.info('Run canceled. Check "Yes, proceed without calibration file" to continue.')
                        return
                    try:
                        _status_factory = getattr(st, 'status', None)
                        if callable(_status_factory):
                            with _status_factory('Run in progress…', expanded=True) as run_status:
                                run_status.write('DAQ…')
                                with st.spinner('Acquiring (DAQ)…'):
                                    b.run_experiment(test_type=tt)
                                st.success('Acquisition finished.')
                                run_status.write('HDF5…')
                                with st.spinner('Writing data file (.h5)…'):
                                    rep = export_primary_h5(
                                        b,
                                        post_trial_notes=notes_in if notes_in else None,
                                        outputfile=outp or None,
                                        append_post_trial_notes=bool(st.session_state.get('gui_qc_notes_append', True)),
                                    )
                                qix = _read_qc_trial_index(b)
                                sel_h5 = str(st.session_state.get('gui_review_selected') or '').strip()
                                qc_base = _qc_figure_base_path(b, sel_h5, qix)
                                run_status.write('QC plot…')
                                with st.spinner('Saving QC plot…'):
                                    qc_path, _ = save_universal_qc_figure(b, qc_trial_index=qix, base_path=qc_base)
                                _st_done = getattr(run_status, 'update', None)
                                if callable(_st_done):
                                    _st_done(label='Run finished', state='complete')
                                st.success('Data has been saved! Check data folder to confirm before proceeding.')
                                st.info(f"Data file: `{rep['outputfile']}`  |  QC plot: `{qc_path}`")
                                if bool(st.session_state.get('gui_qc_notes_append', True)):
                                    st.session_state['gui_post_notes'] = ''
                                else:
                                    st.session_state['gui_post_notes'] = str(rep.get('post_trial_notes') or '')
                        else:
                            with st.spinner('Acquiring (DAQ)…'):
                                b.run_experiment(test_type=tt)
                            st.success('Acquisition finished.')
                            with st.spinner('Writing data file (.h5)…'):
                                rep = export_primary_h5(
                                    b,
                                    post_trial_notes=notes_in if notes_in else None,
                                    outputfile=outp or None,
                                    append_post_trial_notes=bool(st.session_state.get('gui_qc_notes_append', True)),
                                )
                            qix = _read_qc_trial_index(b)
                            sel_h5 = str(st.session_state.get('gui_review_selected') or '').strip()
                            qc_base = _qc_figure_base_path(b, sel_h5, qix)
                            with st.spinner('Saving QC plot…'):
                                qc_path, _ = save_universal_qc_figure(b, qc_trial_index=qix, base_path=qc_base)
                            st.success('Data has been saved! Check data folder to confirm before proceeding.')
                            st.info(f"Data file: `{rep['outputfile']}`  |  QC plot: `{qc_path}`")
                            if bool(st.session_state.get('gui_qc_notes_append', True)):
                                st.session_state['gui_post_notes'] = ''
                            else:
                                st.session_state['gui_post_notes'] = str(rep.get('post_trial_notes') or '')
                    except Exception as e:
                        _show_friendly_error(e, action='run_experiment')

        st.checkbox(
            'Hide section (values stay; unhide to edit)',
            key='gui_sec5_hide',
            help='Collapse validate/run controls after you finish.',
        )

        st.divider()
        st.session_state.setdefault('gui_sec6_hide', False)
        st.subheader('7 · Save data here')
        st.caption('Writes from current in-memory **trial_records** (after **Run** or a prior save). No new DAQ.')
        if st.session_state.get('gui_sec6_hide'):
            st.caption('Save controls hidden. Uncheck **Hide section** below.')
        else:

            def _export_h5_from_session():
                _sync_genus_species_to_bender(b)
                outp = _compose_output_h5_path().strip()
                if not outp and not getattr(b, 'outputfile', None):
                    _st_error_actions(
                        'No HDF5 path set.',
                        ['Open section 2', 'Set data folder', 'Set file name'],
                    )
                    return None
                notes_in = str(st.session_state.get('gui_post_notes') or '').strip()
                if outp:
                    b.outputfile = outp
                return export_primary_h5(
                    b,
                    post_trial_notes=notes_in if notes_in else None,
                    outputfile=outp or None,
                    append_post_trial_notes=bool(st.session_state.get('gui_qc_notes_append', True)),
                )

            def _save_qc_plot_only():
                qix = _read_qc_trial_index(b)
                sel_h5 = str(st.session_state.get('gui_review_selected') or '').strip()
                qc_base = _qc_figure_base_path(b, sel_h5, qix)
                return save_universal_qc_figure(b, qc_trial_index=qix, base_path=qc_base)

            _, c_big, _ = st.columns([1, 2, 1])
            with c_big:
                if _load_save_button('Save Data File (.h5) and QC Plot', key='gui_save_h5_and_qc'):
                    try:
                        with st.spinner('Writing data file (.h5)…'):
                            rep = _export_h5_from_session()
                        if rep is None:
                            pass
                        else:
                            qix = _read_qc_trial_index(b)
                            sel_h5 = str(st.session_state.get('gui_review_selected') or '').strip()
                            qc_base = _qc_figure_base_path(b, sel_h5, qix)
                            try:
                                with st.spinner('Saving QC plot…'):
                                    qc_path, _ = save_universal_qc_figure(b, qc_trial_index=qix, base_path=qc_base)
                            except Exception as e:
                                _show_friendly_error(e, action='save_qc')
                                st.warning(f"The data file was saved: `{rep['outputfile']}`")
                            else:
                                st.success('Data file and QC plot saved.')
                                st.info(f"Data file: `{rep['outputfile']}`  |  QC plot: `{qc_path}`")
                            if bool(st.session_state.get('gui_qc_notes_append', True)):
                                st.session_state['gui_post_notes'] = ''
                            else:
                                st.session_state['gui_post_notes'] = str(rep.get('post_trial_notes') or '')
                    except Exception as e:
                        _show_friendly_error(e, action='save_h5')

            e1, e2 = st.columns(2)
            with e1:
                if _load_save_button('Only save Data File (.h5)', key='gui_save_h5_only'):
                    try:
                        rep = _export_h5_from_session()
                        if rep is not None:
                            st.success(f"{rep['message']}  →  `{rep['outputfile']}`")
                            if bool(st.session_state.get('gui_qc_notes_append', True)):
                                st.session_state['gui_post_notes'] = ''
                            else:
                                st.session_state['gui_post_notes'] = str(rep.get('post_trial_notes') or '')
                    except Exception as e:
                        _show_friendly_error(e, action='save_h5')
            with e2:
                if _load_save_button('Only Save QC Plot', key='gui_save_qc_only'):
                    try:
                        qc_path, _ = _save_qc_plot_only()
                        st.success(f'QC plot saved: `{qc_path}`')
                    except Exception as e:
                        _show_friendly_error(e, action='save_qc')

        st.checkbox(
            'Hide section (values stay; unhide to edit)',
            key='gui_sec6_hide',
            help='Collapse save buttons when finished.',
        )

    if _show_sec7_and_8():

        st.divider()
        st.session_state.setdefault('gui_sec7_hide', False)
        st.subheader('8 · Visualize experimental data')
        if st.session_state.get('gui_sec7_hide'):
            st.caption('Visualization panel hidden. Uncheck **Hide section** below.')
        if not st.session_state.get('gui_sec7_hide'):
            data_path = _output_path_anchor_for_review(b)
            review_files = _candidate_review_files(data_path) if data_path else []
            if not review_files:
                st.info(
                    'No `.h5` / image / HTML files found in the current data folder yet. '
                    'Set **section 2** **Data folder** and **Data file name**, or run/export to create files.'
                )
            else:
                selected_file = st.session_state.get('gui_review_selected')
                if selected_file not in review_files:
                    st.session_state['gui_review_selected'] = review_files[0]
                    selected_file = review_files[0]
                st.caption(
                    f'**Selected file:** `{os.path.basename(selected_file)}` — choose or confirm the same file in **section 9 · Add note**.'
                )
                ext = os.path.splitext(str(selected_file).lower())[1]
                if ext in ('.png', '.jpg', '.jpeg', '.webp'):
                    st.image(selected_file, caption=os.path.basename(selected_file))
                else:
                    st.caption(f"Selected file: `{selected_file}`")
    
                st.markdown('**Custom plots from data file**')
                if ext != '.h5':
                    st.info('Choose a **`.h5`** data file in **section 9 · Add note** to plot saved time series (torque, angle, stim, etc.).')
                else:
                    summ = h5_custom_plot_summary(selected_file)
                    if not summ['ok']:
                        st.warning(summ.get('error') or 'Could not read this HDF5 file.')
                    else:
                        m1, m2, m3, m4 = st.columns(4)
                        with m1:
                            st.metric('test_type (file)', summ['test_type'] or '—')
                        with m2:
                            st.metric('Trials', str(summ['n_trials']))
                        with m3:
                            st.metric('schema', summ['schema_version'] or '—')
                        with m4:
                            ax = summ.get('primary_bending_axis') or ''
                            st.metric('primary axis (saved)', ax or '—')
    
                        trial_names = summ['trial_names']
                        if not trial_names:
                            st.info('No `trial_*` groups found under `02_TimeSeries`.')
                        else:
                            def_trial = trial_names[0]
                            if st.session_state.get('gui_h5_trial_select') not in trial_names:
                                st.session_state['gui_h5_trial_select'] = def_trial
                            trial_sel = st.selectbox(
                                'Trial',
                                trial_names,
                                key='gui_h5_trial_select',
                                help='Each trial is one acquisition segment stored in the file.',
                            )
                            plot_vars = list_h5_plot_variables(
                                selected_file,
                                trial_sel,
                                channel_names=summ.get('channel_names'),
                            )
                            if not plot_vars:
                                st.info('No plottable numeric series in this trial (need length ≥ 2).')
                            else:
                                ids = [v['id'] for v in plot_vars]
                                id_to_label = {v['id']: v['label'] for v in plot_vars}
                                tt0 = plot_vars[0].get('trial_test_type') or ''
                                if tt0:
                                    st.caption(f'**test_type (this trial):** `{tt0}`')
    
                                n_panel = st.number_input(
                                    'Number of figure panels',
                                    min_value=1,
                                    max_value=4,
                                    value=1,
                                    step=1,
                                    key='gui_h5_n_panels',
                                    help='Each panel is a separate Plotly figure. Use multiple Y traces in one panel for several lines on the same axes.',
                                )
                                for p in range(int(n_panel)):
                                    st.markdown(f'**Panel {p + 1}**')
                                    cxa, cya = st.columns(2)
                                    with cxa:
                                        st.selectbox(
                                            'X (time, angle, stim, …)',
                                            ids,
                                            key=f'gui_h5_x_{p}',
                                            format_func=lambda i: id_to_label.get(i, i),
                                        )
                                    with cya:
                                        st.multiselect(
                                            'Y — one or more traces (same plot)',
                                            ids,
                                            key=f'gui_h5_y_{p}',
                                            format_func=lambda i: id_to_label.get(i, i),
                                            help='Pick several series to overlay on the same axes (e.g. raw vs corrected torque).',
                                        )
    
                                if st.button('Generate plots from file', key='gui_h5_gen_plots'):
                                    for p in range(int(n_panel)):
                                        x_id = st.session_state.get(f'gui_h5_x_{p}')
                                        y_ids = st.session_state.get(f'gui_h5_y_{p}') or []
                                        if not x_id or not y_ids:
                                            _st_error_actions(
                                                f'Panel {p + 1} incomplete.',
                                                ['Choose X axis', 'Choose one or more Y'],
                                            )
                                            continue
                                        try:
                                            x_raw = read_h5_series(selected_file, trial_sel, x_id)
                                            fig_p = go.Figure()
                                            x_label = id_to_label.get(x_id, x_id)
                                            for y_id in y_ids:
                                                y_raw = read_h5_series(selected_file, trial_sel, y_id)
                                                x, y = align_xy(x_raw, y_raw)
                                                y_label = id_to_label.get(y_id, y_id)
                                                fig_p.add_trace(
                                                    go.Scatter(x=x, y=y, mode='lines', name=y_label),
                                                )
                                            fig_p.update_layout(
                                                title=f'Panel {p + 1}: {os.path.basename(selected_file)} · {trial_sel}',
                                                xaxis_title=x_label,
                                                yaxis_title='Y (multiple traces)',
                                                height=420,
                                                legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='right', x=1),
                                            )
                                            st.plotly_chart(fig_p, use_container_width=True)
                                        except Exception as e:
                                            _st_error_actions(
                                                f'Panel {p + 1} plot failed.',
                                                ['Check trial selection', 'Open Technical details'],
                                            )
                                            with st.expander('Technical details'):
                                                st.exception(e)
    
                st.markdown('**Add post-experiment note to selected file**')
                note_file = st.text_area(
                    'New note text (appended for .h5 files)',
                    key='gui_selected_file_note',
                    height=90,
                )
                if _load_save_button('Append note to selected file', key='gui_append_note_to_h5'):
                    try:
                        if ext != '.h5':
                            raise ValueError('Selected file is not .h5. Choose a data file to append notes.')
                        _append_note_to_h5_file(selected_file, note_file)
                        st.success('Note appended to selected data file.')
                    except Exception as e:
                        _show_friendly_error(e, action='save_h5')

        st.checkbox(
            'Hide section (values stay; unhide to edit)',
            key='gui_sec7_hide',
            help='Collapse plots and file preview when finished.',
        )

        st.divider()
        st.session_state.setdefault('gui_sec8_hide', False)
        st.subheader('9 · Add note')
        st.caption(
            'Write down anything you noticed about the specimen, apparatus, data quality, etc. that may impact analysis '
            'downstream. Pick a data-folder file below for QC plot naming and **section 8** plots. '
            'Notes here are saved into the `.h5` on export. After **Run** or **Save Data File (.h5) and QC Plot**, the app '
            'writes data and a QC plot; PNG needs **kaleido** (`pip install kaleido`), otherwise HTML.'
        )
        if st.session_state.get('gui_sec8_hide'):
            st.caption('Note controls hidden. Uncheck **Hide section** below.')
        if not st.session_state.get('gui_sec8_hide'):
            data_path_qc = _output_path_anchor_for_review(b)
            review_files_qc = _candidate_review_files(data_path_qc) if data_path_qc else []
            if not review_files_qc:
                st.info(
                    'Set **section 2** **Data folder** and **Data file name** (or save once) so files appear here.'
                )
            else:
                _ensure_review_file_selection(review_files_qc)
                st.selectbox(
                    'Data folder file (QC plot name, section 8 plots)',
                    review_files_qc,
                    key='gui_review_selected',
                    format_func=lambda fp: os.path.basename(fp),
                    help=(
                        'Files in your **data folder**. Choosing an `.h5` sets the stem for the next QC plot export '
                        '(**section 7**).'
                    ),
                )

            tr_qc = list(getattr(b, 'trial_records', []) or [])
            n_tr = len(tr_qc)
            if n_tr > 1:
                opts_ix = list(range(n_tr))
                _cur_ix = st.session_state.get('gui_qc_trial_index')
                if _cur_ix not in opts_ix:
                    st.session_state['gui_qc_trial_index'] = opts_ix[-1]

                def _qc_trial_label(i: int) -> str:
                    r = tr_qc[i]
                    stp = r.get('step_index', r.get('trial_index', r.get('cycle_index')))
                    ttp = str(r.get('test_type', '') or '')
                    parts = [f'Trial {i}']
                    if ttp:
                        parts.append(ttp)
                    if stp is not None:
                        parts.append(f'step {stp}')
                    return ' — '.join(parts) if len(parts) > 1 else parts[0]

                st.selectbox(
                    'In-memory trial for QC plot (this session)',
                    opts_ix,
                    format_func=_qc_trial_label,
                    key='gui_qc_trial_index',
                    help=(
                        'QC plot uses **trial_records** in memory after Run or save. '
                        'Pick which segment when the protocol produced several (e.g. multiple isometric steps).'
                    ),
                )
            else:
                st.caption(
                    '**QC trial:** only one in-memory segment — no trial choice needed. '
                    '(After multi-step runs, a trial list appears here.)'
                )

            st.markdown('**Note text for next save (stored on `.h5`)**')
            st.text_area(
                'Comment for export (optional)',
                key='gui_post_notes',
                height=100,
                help=(
                    'On save/export, this text is **appended** to any existing note with a timestamp. Leave empty to skip. '
                    'Uncheck **Append…** below to replace the stored note entirely.'
                ),
            )
            st.checkbox(
                'Append new text to existing note in this data file path (when the file already exists)',
                value=True,
                key='gui_qc_notes_append',
            )

        st.checkbox(
            'Hide section (values stay; unhide to edit)',
            key='gui_sec8_hide',
            help='Collapse file picker and note fields when finished.',
        )


if __name__ == '__main__':
    main()
