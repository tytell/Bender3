"""
Browser-based prototype GUI for Bender dispatch parameters.

Run from the project directory:
    streamlit run bender_streamlit_gui.py

Select ``test_type`` first; edit fields, use **Save parameters** to copy them onto the
``Bender`` instance, **Check required fields** to validate, **Refresh preview** for a
table/plot of commanded motion (no DAQ), then **Run experiment** when hardware is ready.
"""
from __future__ import annotations

import json
import os
import sys

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

# Project root on path for `bender_functions` and config modules
_ROOT = os.path.dirname(os.path.abspath(__file__))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from bender_functions import Bender  # noqa: E402
from bender_gui_preview import build_protocol_preview  # noqa: E402
from bender_h5_export import export_primary_h5, save_universal_qc_figure  # noqa: E402

MOTION_TYPES = frozenset(
    {'dynamic', 'frequency_sweep', 'frequency_step', 'curvature_step', 'step_change'}
)

# Extra fields for motion-series protocols (not fully enumerated in get_dispatch_schema).
MOTION_GUI_FIELDS = [
    ('duration', 'float', 'Duration (s)'),
    ('all_freqs', 'list_float', 'Frequencies Hz (comma-separated, e.g. 1, 5)'),
    ('all_amps', 'list_float', 'Amplitudes (comma-separated; meaning from all_amps_mode)'),
    ('all_amps_mode', 'select', 'Amplitude interpretation'),
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

RECRUITMENT_OPTIONS = (
    'bilateral_simultaneous',
    'bilateral_sequential',
    'left',
    'right',
)


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


def _render_field(b: Bender, name: str, kind: str, label: str):
    """
    Render one widget; return value to setattr, or None to skip.

    Streamlit keeps keyed widget values in ``st.session_state``. Passing ``value=`` / ``index=``
    every rerun is ignored once the key exists, which can leave the UI stale (e.g. empty text)
    while ``Bender`` still holds good data. We only seed session state when the key is missing.
    """
    sk = _widget_key(name)
    cur = _get_session_value(b, name)

    if kind == 'float':
        if sk not in st.session_state:
            st.session_state[sk] = float(cur) if cur is not None else 0.0
        return float(st.number_input(label, key=sk, format='%.6g'))

    if kind == 'int':
        if sk not in st.session_state:
            st.session_state[sk] = int(cur) if cur is not None else 0
        return int(st.number_input(label, key=sk, step=1))

    if kind == 'optional_int':
        use_key = f'{sk}_use'
        if use_key not in st.session_state:
            st.session_state[use_key] = cur is not None
        use = st.checkbox(f'{label} (use custom)', key=use_key)
        if not use:
            return None
        if sk not in st.session_state:
            st.session_state[sk] = int(cur) if cur is not None else 0
        return int(st.number_input(label, key=sk, step=1))

    if kind == 'bool':
        if sk not in st.session_state:
            st.session_state[sk] = bool(cur) if cur is not None else False
        return bool(st.checkbox(label, key=sk))

    if kind == 'str':
        if sk not in st.session_state:
            st.session_state[sk] = str(cur or '')
        return str(st.text_input(label, key=sk))

    if kind == 'select':
        opts = list(ALL_AMPS_MODE_OPTIONS)
        dv = str(cur or 'strain')
        if dv not in opts:
            dv = 'strain'
        if sk not in st.session_state:
            st.session_state[sk] = dv
        return str(st.selectbox(label, opts, key=sk))

    if kind == 'list_float':
        if sk not in st.session_state:
            if cur is not None:
                try:
                    st.session_state[sk] = ', '.join(str(float(x)) for x in list(cur))
                except (TypeError, ValueError):
                    st.session_state[sk] = ''
            else:
                st.session_state[sk] = ''
        s = str(st.text_input(label, key=sk))
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
        s = str(st.text_input(label, key=sk))
        return _parse_int_list(s)

    if kind == 'json_dict':
        if sk not in st.session_state:
            st.session_state[sk] = json.dumps(cur, indent=2) if isinstance(cur, dict) else '{}'
        s = str(st.text_area(label, height=120, key=sk))
        try:
            return json.loads(s)
        except json.JSONDecodeError:
            st.error(f'Invalid JSON for {name}')
            return None

    return None


def _clear_fld_session_keys():
    """Drop procedure widget keys so the next render re-seeds from ``Bender`` (e.g. new config)."""
    for k in list(st.session_state.keys()):
        if k.startswith('fld_'):
            del st.session_state[k]


def _apply_pair(b: Bender, name: str, value):
    if value is None and name == 'random_seed':
        setattr(b, name, None)
        return
    if value is None:
        return
    setattr(b, name, value)


def _apply_form_updates(b: Bender, updates: dict, tt: str):
    for k, v in updates.items():
        _apply_pair(b, k, v)
    b.test_type = tt


def _sync_biometric_flags_from_session(b: Bender):
    """Keep inertia flags aligned when the biometrics panel is collapsed."""
    if 'bio_use_theoretical_inertial' in st.session_state:
        b.use_theoretical_inertial_correction = bool(st.session_state['bio_use_theoretical_inertial'])
    if 'bio_use_frustum' in st.session_state:
        b.use_frustum_inertial_model = bool(st.session_state['bio_use_frustum'])


def _init_biometrics_session_state(b: Bender, *, force: bool = False):
    """Seed Streamlit widget keys from ``b`` (``force`` overwrites after config reload)."""

    def _put(key, val):
        if force or key not in st.session_state:
            st.session_state[key] = val

    dc = getattr(b, 'dclamp', None)
    xw = getattr(b, 'xsec_width', None)
    _put('bio_dclamp', float(dc) if dc is not None else 10.0)
    _put('bio_xsec', float(xw) if xw is not None else 8.0)
    fi = getattr(b, 'frustum_inputs', None) or {}
    _put('bio_f_height', float(fi.get('height_mm', 4.0)))
    _put('bio_f_width', float(fi.get('width_mm', 6.0)))
    _put('bio_f_length', float(fi.get('length_mm', 25.0)))
    _put('bio_f_rho', float(fi.get('density_g_per_mm3', 1.03e-3)))
    _put('bio_f_tip', float(fi.get('tip_scale', 0.85)))
    _put('bio_f_clamp_off', float(fi.get('clamp_offset_mm', 20.0)))
    _put('bio_f_samples', 100)
    _put('bio_use_frustum', bool(getattr(b, 'use_frustum_inertial_model', True)))
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


def main():
    if 'gui_outputfile' not in st.session_state:
        st.session_state['gui_outputfile'] = ''
    if 'gui_post_notes' not in st.session_state:
        st.session_state['gui_post_notes'] = ''

    st.title('Bender dispatch (prototype)')
    st.caption(
        'Work top to bottom: **1 specimen geometry** → **2 experiment** → **3 preview** → **QC note** (optional) '
        '→ **4 run** → **5 export**. Use **Save parameters** before preview or run when you change procedure fields.'
    )

    with st.sidebar:
        cfg_mod = st.text_input('Config module (no .py)', value='jimenez_bender_config_A')
        if st.button('Load config & build experiment object', type='primary'):
            try:
                st.session_state['bender'] = Bender(cfg_mod)
                st.session_state['cfg_mod'] = cfg_mod
                b0 = st.session_state['bender']
                st.session_state['gui_outputfile'] = str(getattr(b0, 'outputfile', '') or '')
                st.session_state['gui_post_notes'] = str(getattr(b0, 'post_trial_notes', '') or '')
                _init_biometrics_session_state(b0, force=True)
                _clear_fld_session_keys()
                st.success(f'Loaded {cfg_mod}')
            except Exception as e:
                st.error(str(e))
                st.session_state.pop('bender', None)

        st.divider()
        st.subheader('Output / export')
        st.text_input(
            'Data File (.h5)',
            key='gui_outputfile',
            help='Full path for the saved experiment file. Example: C:\\Data\\run001.h5',
        )
        st.number_input(
            'QC plot: trial index',
            min_value=0,
            value=0,
            step=1,
            key='gui_qc_trial_index',
            help=(
                'After acquisition, each run stores one or more entries in `trial_records`. '
                'This index picks **which trial** the universal QC figure uses (0 = first trial). '
                'Typical motion runs have a single trial → leave at **0**. '
                'Multi-step protocols (e.g. many isometric steps) may append several trials → increase the index '
                'to QC a specific step.'
            ),
        )
        st.checkbox('After run: save data file (.h5)', value=True, key='gui_save_h5')
        st.checkbox('After run: save QC figure (PNG or HTML)', value=True, key='gui_save_qc')
        st.caption('PNG needs **kaleido** (`pip install kaleido`). Otherwise an HTML QC file is written.')

    if 'bender' not in st.session_state:
        st.info('Use the sidebar to load a Bender config module.')
        return

    b: Bender = st.session_state['bender']
    _init_biometrics_session_state(b, force=False)
    _sync_biometric_flags_from_session(b)
    schema = b.get_dispatch_schema()
    test_types = list(schema['test_types'])

    st.subheader('1 · Specimen geometry & biometrics')
    st.session_state.setdefault('gui_bio_hide', False)
    hide_bio = st.checkbox(
        'Hide biometrics & geometry (values stay on the experiment object; uncheck to edit again)',
        key='gui_bio_hide',
    )
    if hide_bio:
        st.caption(
            'Section hidden. Uncheck **Hide biometrics & geometry** to edit **d_clamp**, **xsec_width**, '
            'frustum, or profiled specimen inputs.'
        )
    else:
        st.markdown('**Bending geometry** (strain ↔ curvature ↔ motor angle)')
        st.caption(
            '**d_clamp**: effective bending lever / clamp spacing (mm). **xsec_width**: specimen width for strain '
            '(mm). Both are required for strain-based amplitudes and for preview/run math.'
        )
        g1, g2 = st.columns(2)
        with g1:
            st.number_input('d_clamp (mm)', min_value=0.001, format='%.6g', key='bio_dclamp')
        with g2:
            st.number_input('Cross-section width xsec_width (mm)', min_value=0.001, format='%.6g', key='bio_xsec')
        if st.button('Apply bending geometry to experiment object', key='bio_btn_geom'):
            b.dclamp = float(st.session_state['bio_dclamp'])
            b.xsec_width = float(st.session_state['bio_xsec'])
            st.success('Updated `dclamp` and `xsec_width`.')

        st.divider()
        st.markdown('**Elliptical frustum (specimen inertia)**')
        st.caption(
            'Feeds `set_frustum_inertial_model` — used with theoretical inertial correction when enabled.'
        )
        f1, f2, f3 = st.columns(3)
        with f1:
            st.number_input('Height (mm)', min_value=0.001, format='%.6g', key='bio_f_height')
            st.number_input('Width (mm)', min_value=0.001, format='%.6g', key='bio_f_width')
        with f2:
            st.number_input('Length (mm)', min_value=0.001, format='%.6g', key='bio_f_length')
            st.number_input('Density (g / mm³)', min_value=1e-9, format='%.6g', key='bio_f_rho')
        with f3:
            st.number_input('Tip scale (0–1)', min_value=0.0, max_value=1.0, format='%.6g', key='bio_f_tip')
            st.number_input('Clamp offset (mm)', min_value=0.0, format='%.6g', key='bio_f_clamp_off')
        st.number_input('Frustum integration samples', min_value=10, max_value=500, step=10, key='bio_f_samples')
        st.checkbox('Use frustum-based specimen inertia flag', key='bio_use_frustum')
        if st.button('Apply frustum inertial model', key='bio_btn_frustum'):
            b.use_frustum_inertial_model = bool(st.session_state['bio_use_frustum'])
            b.set_frustum_inertial_model(
                st.session_state['bio_f_height'],
                st.session_state['bio_f_width'],
                st.session_state['bio_f_length'],
                st.session_state['bio_f_rho'],
                tip_scale=float(st.session_state['bio_f_tip']),
                clamp_offset_mm=float(st.session_state['bio_f_clamp_off']),
                num_samples=int(st.session_state['bio_f_samples']),
            )
            st.success('Frustum inertia model updated.')

        st.divider()
        st.markdown('**Optional: profiled specimen (proximal → distal, two stations)**')
        st.caption('Builds a simple two-segment outline via `make_profile_stations` + `set_profiled_specimen_inertial_model`.')
        p1, p2 = st.columns(2)
        with p1:
            st.number_input('Proximal height (mm)', min_value=0.001, format='%.6g', key='bio_prof_ph')
            st.number_input('Proximal width (mm)', min_value=0.001, format='%.6g', key='bio_prof_pw')
        with p2:
            st.number_input('Distal height (mm)', min_value=0.001, format='%.6g', key='bio_prof_dh')
            st.number_input('Distal width (mm)', min_value=0.001, format='%.6g', key='bio_prof_dw')
        p3, p4 = st.columns(2)
        with p3:
            st.number_input('Specimen length (mm)', min_value=0.001, format='%.6g', key='bio_prof_L')
            st.number_input('Density (g / mm³)', min_value=1e-9, format='%.6g', key='bio_prof_rho')
        with p4:
            st.number_input('Clamp offset (mm)', min_value=0.0, format='%.6g', key='bio_prof_clamp')
            st.number_input('Profile integration samples', min_value=20, max_value=400, step=10, key='bio_prof_samples')
        if st.button('Apply profiled inertial model', key='bio_btn_profile'):
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
            st.success('Profiled specimen inertia model updated.')

        st.divider()
        st.checkbox(
            'Use theoretical inertial correction (specimen + hardware model)',
            key='bio_use_theoretical_inertial',
        )

    st.divider()
    st.subheader('2 · Experiment type & parameters')

    tt = st.selectbox('Experiment type (test_type)', test_types, key='test_type_select')
    b.test_type = tt

    st.session_state.setdefault('gui_exp_hide', False)
    st.checkbox(
        'Collapse procedure fields (values stay for preview & run — open the bar below to edit)',
        key='gui_exp_hide',
    )
    st.caption(
        'Fields are stored in the browser session. **Refresh preview** copies them onto the experiment object. '
        'Collapsing only hides this panel; widgets still keep your entries. If anything looks out of sync, use '
        '**Save parameters** or reload the config module.'
    )

    updates = {}

    with st.expander('Procedure fields', expanded=not bool(st.session_state.get('gui_exp_hide'))):
        if tt == 'isometric':
            st.markdown('**Required**')
            for key in schema['isometric_required']:
                label = key.replace('_', ' ')
                updates[key] = _render_field(b, key, 'float' if 'steps' not in key else 'int', label)
            if 'isometric_num_steps' in updates and updates['isometric_num_steps'] is not None:
                updates['isometric_num_steps'] = int(updates['isometric_num_steps'])
            st.markdown('**Optional**')
            for key in schema['isometric_optional']:
                if key in ('isometric_stim_params', 'isometric_stim_overrides'):
                    updates[key] = _render_field(b, key, 'json_dict', key)
                elif key == 'recruitment':
                    skr = _widget_key('recruitment')
                    cur_r = _get_session_value(b, 'recruitment', 'bilateral_simultaneous')
                    if skr not in st.session_state:
                        st.session_state[skr] = cur_r if cur_r in RECRUITMENT_OPTIONS else RECRUITMENT_OPTIONS[0]
                    updates[key] = st.selectbox('recruitment', list(RECRUITMENT_OPTIONS), key=skr)
                elif key == 'lateral_mode':
                    skl = _widget_key('lateral_mode')
                    if skl not in st.session_state:
                        st.session_state[skl] = str(_get_session_value(b, key) or '')
                    updates[key] = st.text_input('lateral_mode', key=skl)
                elif key in ('bilateral_mirror_motor',):
                    updates[key] = _render_field(b, key, 'bool', key)
                elif key == 'bilateral_sequential_left_frac':
                    updates[key] = _render_field(b, key, 'float', key)
                elif key == 'isometric_mode':
                    modes = ('strain', 'strain_pct', 'curvature', 'angle')
                    skm = _widget_key('isometric_mode')
                    cur_m = str(_get_session_value(b, key, 'strain'))
                    if skm not in st.session_state:
                        st.session_state[skm] = cur_m if cur_m in modes else 'strain'
                    updates[key] = st.selectbox('isometric_mode', modes, key=skm)
                elif 'random_seed' in key:
                    sks = _widget_key(key)
                    if sks not in st.session_state:
                        v0 = _get_session_value(b, key)
                        st.session_state[sks] = '' if v0 is None else str(v0)
                    s = st.text_input(key.replace('_', ' '), key=sks)
                    if not str(s).strip():
                        updates[key] = None
                    else:
                        try:
                            updates[key] = int(s)
                        except ValueError:
                            st.error(f'{key}: need integer or empty')
                            updates[key] = None
                else:
                    kind = 'bool' if 'randomize' in key else 'str'
                    updates[key] = _render_field(b, key, kind, key)

        elif tt == 'isovelocity':
            st.markdown('**Required**')
            for key in schema['isovelocity_required']:
                kind = 'int' if 'num_steps' in key else 'float'
                updates[key] = _render_field(b, key, kind, key.replace('_', ' '))
            if 'isovelocity_num_steps' in updates and updates['isovelocity_num_steps'] is not None:
                updates['isovelocity_num_steps'] = int(updates['isovelocity_num_steps'])
            st.markdown('**Optional**')
            for key in schema['isovelocity_optional']:
                if key in ('isovelocity_stim_params', 'isovelocity_stim_overrides'):
                    updates[key] = _render_field(b, key, 'json_dict', key)
                elif key == 'recruitment':
                    skr = _widget_key('recruitment')
                    cur_r = _get_session_value(b, 'recruitment', 'bilateral_simultaneous')
                    if skr not in st.session_state:
                        st.session_state[skr] = cur_r if cur_r in RECRUITMENT_OPTIONS else RECRUITMENT_OPTIONS[0]
                    updates[key] = st.selectbox('recruitment', list(RECRUITMENT_OPTIONS), key=skr)
                elif key == 'lateral_mode':
                    skl = _widget_key('lateral_mode')
                    if skl not in st.session_state:
                        st.session_state[skl] = str(_get_session_value(b, key) or '')
                    updates[key] = st.text_input('lateral_mode', key=skl)
                elif key in ('bilateral_mirror_motor',):
                    updates[key] = _render_field(b, key, 'bool', key)
                elif key == 'bilateral_sequential_left_frac':
                    updates[key] = _render_field(b, key, 'float', key)
                elif key == 'isovelocity_starting_strain_mode':
                    modes = ('strain', 'strain_pct', 'curvature', 'angle')
                    skm = _widget_key('isovelocity_starting_strain_mode')
                    cur_m = str(_get_session_value(b, key, 'strain'))
                    if skm not in st.session_state:
                        st.session_state[skm] = cur_m if cur_m in modes else 'strain'
                    updates[key] = st.selectbox('isovelocity_starting_strain_mode', modes, key=skm)
                elif 'random_seed' in key:
                    sks = _widget_key(key)
                    if sks not in st.session_state:
                        v0 = _get_session_value(b, key)
                        st.session_state[sks] = '' if v0 is None else str(v0)
                    s = st.text_input(key.replace('_', ' '), key=sks)
                    if not str(s).strip():
                        updates[key] = None
                    else:
                        try:
                            updates[key] = int(s)
                        except ValueError:
                            st.error(f'{key}: need integer or empty')
                            updates[key] = None
                else:
                    kind = 'bool' if 'randomize' in key else 'float'
                    updates[key] = _render_field(b, key, kind, key)

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
                '(e.g. dynamic), **Save parameters** (and **Refresh preview** if you use it), then switch '
                'back to **calibration** before running.'
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
                    'If `period_by_cycle` is not set yet, call **organize_cycles** from a notebook/script after **Save parameters**.'
                )
            elif tt == 'step_change':
                st.caption('Motion length is computed inside **make_cycles_step_change**; **Duration** does not apply.')
            fields = _motion_parameter_rows(tt)
            for name, kind, label in fields:
                updates[name] = _render_field(b, name, kind, label)

        else:
            st.warning(f'No dedicated field panel for {tt!r} yet; use notebook or extend this script.')

    st.divider()
    st.subheader('3 · Experiment preview (table & plot, no DAQ)')
    st.caption(
        'Uses the same motion math as a real run for **commanded** angle / velocity. '
        '**Refresh preview** copies the form onto the experiment object, then recomputes commands. '
        'For **dynamic**, preview calls **organize_cycles** and updates `period_by_cycle`, so a following **Run** '
        'matches the preview if you do not overwrite those arrays elsewhere. '
        'You still need **dclamp** and **xsec_width** on the object (specimen metadata).'
    )
    pv_on = st.checkbox('Show last preview', value=True, key='gui_show_preview')
    pv_pts = st.slider('Preview plot resolution (max points)', 400, 12000, 6000, step=200, key='gui_preview_pts')
    if st.button('Refresh preview'):
        _sync_biometric_flags_from_session(b)
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
            st.error(prev['error'])
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
            if tp is not None and ap is not None and vp is not None and len(tp) > 0:
                st.markdown('**Command preview** (motor angle and angular velocity)')
                fig = go.Figure()
                fig.add_trace(go.Scatter(x=tp, y=ap, mode='lines', name='Commanded angle (deg)'))
                fig.add_trace(
                    go.Scatter(
                        x=tp,
                        y=vp,
                        mode='lines',
                        name='Commanded anglevel (deg/s)',
                        yaxis='y2',
                    )
                )
                fig.update_layout(
                    height=420,
                    margin=dict(l=48, r=48, t=40, b=40),
                    xaxis_title='Time (s)',
                    yaxis=dict(title='Angle (deg)'),
                    yaxis2=dict(title='Anglevel (deg/s)', overlaying='y', side='right'),
                    legend=dict(yanchor='top', y=0.99, xanchor='left', x=0.01),
                )
                st.plotly_chart(fig, use_container_width=True)
            elif prev.get('table') and tt in ('isometric', 'isovelocity'):
                st.caption('Step protocols: table lists setpoints; full continuous traces are built during the run.')
        else:
            st.warning('Preview incomplete; click **Refresh preview** again.')

    with st.expander('Post-experiment / QC assessment (written into the data file)', expanded=False):
        st.caption(
            'Use this **after** you review **QC plots** and the **rig / setup**. Your text is stored in the '
            '**Data File (.h5)** as HDF5 attributes (`post-trial notes` on the root and in `01_Metadata`), so it '
            'travels with the experiment record for analysis or archival.'
        )
        st.text_area(
            'New QC / rig / run-quality comment (optional)',
            key='gui_post_notes',
            height=100,
            help=(
                'On export, this text is **appended** to the note already stored in that file (if any), with a '
                'timestamp. Leave empty to re-export data without adding a new paragraph. Uncheck **Append…** '
                'below to replace the stored note entirely.'
            ),
        )
        st.checkbox(
            'Append new text to existing note in this data file path (when the file already exists)',
            value=True,
            key='gui_qc_notes_append',
        )

    st.divider()
    st.subheader('4 · Save, validate, and run')
    c1, c2, c3 = st.columns(3)
    with c1:
        if st.button(
            'Save parameters',
            help='Copies the form into the in-memory experiment object (same as what preview and run use).',
        ):
            _sync_biometric_flags_from_session(b)
            for k, v in updates.items():
                _apply_pair(b, k, v)
            b.test_type = tt
            _pn = str(st.session_state.get('gui_post_notes') or '').strip()
            if _pn:
                b.post_trial_notes = _pn
            st.success('Form values stored on the experiment object.')
    with c2:
        if st.button(
            'Check required fields',
            help='Runs validation for the current test type (does not talk to hardware).',
        ):
            _sync_biometric_flags_from_session(b)
            _apply_form_updates(b, updates, tt)
            rep = b.validate_dispatch_setup(test_type=tt)
            if rep['ok']:
                st.success('All required fields for this procedure are set.')
            else:
                st.error('Still needed: ' + ', '.join(rep['missing']))
    with c3:
        daq_ok = st.checkbox('Hardware: I intend to run DAQ', value=False)
        if st.button('Run experiment', type='primary', disabled=not daq_ok):
            _sync_biometric_flags_from_session(b)
            _apply_form_updates(b, updates, tt)
            outp = (st.session_state.get('gui_outputfile') or '').strip()
            if outp:
                b.outputfile = outp
            notes_in = str(st.session_state.get('gui_post_notes') or '').strip()
            try:
                with st.spinner('Acquiring (DAQ)…'):
                    b.run_experiment(test_type=tt)
                st.success('Acquisition finished.')
                if st.session_state.get('gui_save_h5'):
                    with st.spinner('Writing data file (.h5)…'):
                        rep = export_primary_h5(
                            b,
                            post_trial_notes=notes_in if notes_in else None,
                            outputfile=outp or None,
                            append_post_trial_notes=bool(st.session_state.get('gui_qc_notes_append', True)),
                        )
                    st.success(f"{rep['message']}  →  `{rep['outputfile']}`")
                    if bool(st.session_state.get('gui_qc_notes_append', True)):
                        st.session_state['gui_post_notes'] = ''
                    else:
                        st.session_state['gui_post_notes'] = str(rep.get('post_trial_notes') or '')
                if st.session_state.get('gui_save_qc'):
                    qix = int(st.session_state.get('gui_qc_trial_index') or 0)
                    with st.spinner('Saving QC figure…'):
                        qc_path, _ = save_universal_qc_figure(b, qc_trial_index=qix)
                    st.success(f'QC saved: `{qc_path}`')
            except Exception as e:
                st.exception(e)

    st.divider()
    st.subheader('5 · Export without re-running')
    e1, e2 = st.columns(2)
    with e1:
        if st.button('Export data file only (.h5, from current trial_records)'):
            outp = (st.session_state.get('gui_outputfile') or '').strip()
            if not outp and not getattr(b, 'outputfile', None):
                st.error('Set **Data File (.h5)** in the sidebar first.')
            else:
                notes_in = str(st.session_state.get('gui_post_notes') or '').strip()
                if outp:
                    b.outputfile = outp
                try:
                    rep = export_primary_h5(
                        b,
                        post_trial_notes=notes_in if notes_in else None,
                        outputfile=outp or None,
                        append_post_trial_notes=bool(st.session_state.get('gui_qc_notes_append', True)),
                    )
                    st.success(f"{rep['message']}  →  `{rep['outputfile']}`")
                    if bool(st.session_state.get('gui_qc_notes_append', True)):
                        st.session_state['gui_post_notes'] = ''
                    else:
                        st.session_state['gui_post_notes'] = str(rep.get('post_trial_notes') or '')
                except Exception as e:
                    st.exception(e)
    with e2:
        if st.button('Save QC figure only'):
            try:
                qix = int(st.session_state.get('gui_qc_trial_index') or 0)
                qc_path, _ = save_universal_qc_figure(b, qc_trial_index=qix)
                st.success(f'QC saved: `{qc_path}`')
            except Exception as e:
                st.exception(e)

    with st.expander('Raw dispatch schema (JSON)'):
        st.json(schema)


if __name__ == '__main__':
    main()
