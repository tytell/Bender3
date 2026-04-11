"""Full CritterGripper NiceGUI workflow (sections 1–6 + H5 explorer + landing)."""
from __future__ import annotations

import importlib
import os
import sys
from dataclasses import dataclass, field
from typing import Any, Dict, List, Literal, Tuple

import pandas as pd
import plotly.graph_objects as go
from nicegui import app, ui

from bender_biometrics_templates import (
    apply_biometrics_template_to_session,
    biometrics_template_display_label,
    list_biometrics_template_files,
    load_biometrics_template,
    save_biometrics_template,
    sanitize_biometrics_filename_stem,
)
from bender_config_builder import (
    discover_config_modules,
    effective_load_module_name,
    parse_comma_list,
    parse_n_floats,
    read_base_defaults,
    render_generated_config,
    sanitize_config_module_stem,
)
from bender_daq_kill import daq_emergency_stop
from bender_functions import Bender
from bender_gui_preview import build_protocol_preview
from bender_h5_export import export_primary_h5, save_universal_qc_figure
from bender_h5_explore import (
    align_xy,
    build_series_catalog_legacy,
    build_series_catalog_v2,
    detect_h5_schema,
    list_v2_trials,
)
from bender_protocol_templates import (
    apply_template_to_session_state,
    build_procedure_dict_from_updates,
    list_template_files,
    load_protocol_template,
    save_protocol_template,
    sanitize_template_filename_stem,
    snapshot_bender_procedure,
    template_display_label,
)

from .constants import (
    BIO_DBEND_FIELD_HELP,
    BIO_DENSITY_PRESET_LABELS,
    BIO_PROF_CLAMP_FIELD_HELP,
    DATA_FILE_NAME_HELP,
    DATA_FOLDER_HELP,
)
from .procedure_ops import gather_procedure_updates, widget_key as proc_widget_key
from .procedure_ui import render_procedure_fields
from .session_logic import (
    apply_all_biometrics_to_bender,
    apply_form_updates,
    apply_procedure_form_to_bender,
    apply_specimen_identity_to_bender,
    clear_fld_session_keys,
    compose_output_h5_path,
    ensure_gui_defaults,
    init_biometrics_session_state,
    needs_missing_calibration_confirmation,
    qc_figure_base_path,
    read_qc_trial_index,
    section2_destination_incomplete,
    shared_experiment_dir,
    sync_biometric_flags_from_session,
    sync_genus_species_to_bender,
)
from .styles import (
    ALERT,
    BTN_COMMIT,
    BTN_MODE_OFF,
    BTN_MODE_ON,
    BTN_SIDE,
    CAPTION,
    LABEL,
    MAIN_COLUMN,
    PANEL,
    SECTION_TITLE,
    SHELL_BODY,
    SKIP_LINK,
    SUCCESS,
)

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)


def normalize_config_module_name(raw: str) -> str:
    s = str(raw or '').strip()
    if s.lower().endswith('.py'):
        s = s[:-3].strip()
    return s


def selected_config_matches_bender(b: Bender, eff_raw: str) -> bool:
    eff = normalize_config_module_name(eff_raw)
    if not eff:
        return False
    loaded = normalize_config_module_name(getattr(b, 'config_name', '') or '')
    return bool(loaded) and loaded == eff


def apply_data_path_to_bender(b: Bender, sess: Dict[str, Any]) -> str | None:
    outp = compose_output_h5_path(sess).strip()
    if not outp:
        return 'Set **Data file name** first.'
    b.outputfile = outp
    return None


def protocol_template_option_label(p: str) -> str:
    return f'{template_display_label(p)}  ·  {os.path.basename(p)}'


def _client_route() -> str:
    """Per-browser-tab route (not shared across users). Survives SPA navigation; resets on full reload."""
    r = app.storage.client.setdefault('gui_route', 'landing')
    return r if r in ('landing', 'workflow', 'h5') else 'landing'


def _set_client_route(route: Literal['landing', 'workflow', 'h5']) -> None:
    app.storage.client['gui_route'] = route


@dataclass
class AppState:
    bender: Bender | None = None
    sess: Dict[str, Any] = field(default_factory=dict)
    config_mode: Literal['Load existing', 'Build new'] = 'Load existing'
    cfg_mod: str = ''
    load_cfg_select: str = ''
    gui_cfg_build_base: str = ''
    gui_cfg_build_out: str = ''
    gui_cfg_build_overwrite: bool = False
    build: Dict[str, Any] = field(default_factory=dict)
    build_seeded_for: str = ''
    hardware_error: str = ''
    protocol_feedback: Tuple[bool, str] | None = None
    biometrics_feedback: Tuple[bool, str] | None = None


state = AppState()


def _init_hw_defaults() -> None:
    mods = discover_config_modules(_ROOT)
    if not mods:
        mods = ['jimenez_bender_config_A']
    state.load_cfg_select = state.load_cfg_select if state.load_cfg_select in mods else mods[0]
    state.gui_cfg_build_base = state.gui_cfg_build_base if state.gui_cfg_build_base in mods else mods[0]
    if not state.cfg_mod.strip():
        state.cfg_mod = 'jimenez_bender_config_A'


def apply_loaded_config_module(raw_mod: str) -> str | None:
    cm = normalize_config_module_name(raw_mod)
    if not cm:
        return 'Enter a config module name.'
    if _ROOT not in sys.path:
        sys.path.insert(0, _ROOT)
    try:
        b = Bender(cm)
        state.bender = b
        outp0 = str(getattr(b, 'outputfile', '') or '').strip()
        if outp0:
            n0 = os.path.normpath(outp0)
            state.sess['gui_data_folder'] = os.path.dirname(n0) or ''
            state.sess['gui_data_filename'] = os.path.basename(n0)
        else:
            state.sess['gui_data_folder'] = state.sess.get('gui_data_folder', '')
            state.sess['gui_data_filename'] = state.sess.get('gui_data_filename', '')
        ensure_gui_defaults(state.sess)
        init_biometrics_session_state(b, state.sess, force=True)
        clear_fld_session_keys(state.sess)
        state.sess['test_type_select'] = getattr(b, 'test_type', '') or ''
        meta0 = getattr(b, 'h5_protocol_metadata', {}) or {}
        state.sess['gui_genus_species'] = str(meta0.get('genus_species', '') or '')
        state.sess['gui_specimen_id'] = str(meta0.get('specimen_id', '') or getattr(b, 'fishcode', '') or '')
        state.sess['gui_post_notes'] = str(getattr(b, 'post_trial_notes', '') or '')
        return None
    except ImportError as e:
        state.bender = None
        return f'Could not import `{cm}`: {e}'
    except Exception as e:
        state.bender = None
        return f'{type(e).__name__}: {e}'


def seed_build_from_base() -> None:
    base = str(state.gui_cfg_build_base or '').strip()
    if state.build_seeded_for == base:
        return
    try:
        d = read_base_defaults(base)
    except Exception:
        return
    state.build = {
        'forcetorque_calibration_file': str(d['forcetorque_calibration_file']),
        'positive_motor_direction': str(d['positive_motor_direction']),
        'specimen_lateral_index_on_positive_motor_side': int(d['specimen_lateral_index_on_positive_motor_side']),
        'motor_axis': str(d['motor_axis']),
        'bending_axis_sensor': str(d['bending_axis_sensor']),
        'primary_bending_axis': str(d['primary_bending_axis']),
        'bending_axis_specimen': str(d['bending_axis_specimen']),
        'device_name': str(d['device_name']),
        'daq_ai_sample_rate_hz': float(d['daq_ai_sample_rate_hz']),
        'daq_ao_do_sample_rate_hz': float(d['daq_ao_do_sample_rate_hz']),
        'motor_full_steps_per_rev': int(d['motor_full_steps_per_rev']),
        'motor_gear_ratio': int(d['motor_gear_ratio']),
        'encoder_pulses_per_rev': int(d['encoder_pulses_per_rev']),
        'stim_channels': ', '.join(str(x) for x in d['stim_channels']),
        'motor_port': str(d['motor_port']),
        'encoder_chan': str(d['encoder_chan']),
        'SG_chan': ', '.join(str(x) for x in d['SG_chan']),
        'SG_name': ', '.join(str(x) for x in d['SG_name']),
        'stim_monitor_chan': ', '.join(str(x) for x in d['stim_monitor_chan']),
        'stim_monitor_name': ', '.join(str(x) for x in d['stim_monitor_name']),
        'S1side': str(d['S1side']),
        'S2side': str(d['S2side']),
        'use_sono': bool(d['use_sono']),
        'sono_channel': ', '.join(str(x) for x in d['sono_channel']),
        'sono_name': ', '.join(str(x) for x in d['sono_name']),
        'sono_internal_samplefreq': int(d['sono_internal_samplefreq']),
        'sono_cal_left': ', '.join(str(x) for x in d['sono_cal_left']),
        'sono_cal_right': ', '.join(str(x) for x in d['sono_cal_right']),
        'amp_step_vel': int(d['amp_step_vel']),
        'ramp_mode_default': str(d['ramp_mode_default'])
        if str(d['ramp_mode_default']) in ('linear', 'exponential')
        else 'linear',
        'waitbefore': float(d['waitbefore']),
        'waitafter': float(d['waitafter']),
        'rampdur': float(d['rampdur']),
        'prepoststim_dur': float(d['prepoststim_dur']),
        'prepoststim_sep': float(d['prepoststim_sep']),
        'prestim_time': float(d['prestim_time']),
        'poststim_time': float(d['poststim_time']),
    }
    state.build_seeded_for = base


def maybe_build_write() -> None:
    base = str(state.gui_cfg_build_base or '').strip()
    out_raw = str(state.gui_cfg_build_out or '').strip()
    out_stem = sanitize_config_module_stem(out_raw)
    if not base or not out_raw:
        state.hardware_error = 'Pick base module and enter new file name.'
        return
    stim_ch = parse_comma_list(str(state.build.get('stim_channels') or ''))
    sg_ch = parse_comma_list(str(state.build.get('SG_chan') or ''))
    sg_nm = parse_comma_list(str(state.build.get('SG_name') or ''))
    sm_ch = parse_comma_list(str(state.build.get('stim_monitor_chan') or ''))
    sm_nm = parse_comma_list(str(state.build.get('stim_monitor_name') or ''))
    try:
        sono_lf = parse_n_floats(str(state.build.get('sono_cal_left') or ''), 4)
        sono_rf = parse_n_floats(str(state.build.get('sono_cal_right') or ''), 4)
    except ValueError as e:
        state.hardware_error = f'Sono calibration invalid: {e}'
        return
    if not stim_ch:
        state.hardware_error = 'Add stim AO channels.'
        return
    if not sg_ch or not sg_nm or len(sg_ch) != len(sg_nm):
        state.hardware_error = 'SG channels/names mismatch.'
        return
    if (bool(sm_ch) ^ bool(sm_nm)) or (sm_ch and sm_nm and len(sm_ch) != len(sm_nm)):
        state.hardware_error = 'Stim monitor lists mismatch.'
        return
    sono_ch = parse_comma_list(str(state.build.get('sono_channel') or ''))
    sono_nm = parse_comma_list(str(state.build.get('sono_name') or ''))
    use_sono = bool(state.build.get('use_sono'))
    if use_sono and (not sono_ch or not sono_nm or len(sono_ch) != len(sono_nm)):
        state.hardware_error = 'Sono mismatch or disable sono.'
        return
    path = os.path.join(_ROOT, out_stem + '.py')
    if os.path.isfile(path) and not state.gui_cfg_build_overwrite:
        state.hardware_error = 'File exists — enable overwrite.'
        return
    rm = str(state.build.get('ramp_mode_default') or 'linear')
    if rm not in ('linear', 'exponential'):
        rm = 'linear'
    assignments = {
        'forcetorque_calibration_file': str(state.build.get('forcetorque_calibration_file') or 'FT56491.cal'),
        'positive_motor_direction': str(state.build.get('positive_motor_direction') or 'left'),
        'specimen_lateral_index_on_positive_motor_side': int(
            state.build.get('specimen_lateral_index_on_positive_motor_side') or -1
        ),
        'motor_axis': str(state.build.get('motor_axis') or 'z'),
        'bending_axis_sensor': str(state.build.get('bending_axis_sensor') or 'z'),
        'primary_bending_axis': str(state.build.get('primary_bending_axis') or 'zTorque'),
        'bending_axis_specimen': str(state.build.get('bending_axis_specimen') or 'dorsoventral'),
        'device_name': str(state.build.get('device_name') or 'Dev1'),
        'daq_ai_sample_rate_hz': float(state.build.get('daq_ai_sample_rate_hz') or 1000.0),
        'daq_ao_do_sample_rate_hz': float(state.build.get('daq_ao_do_sample_rate_hz') or 60000.0),
        'motor_full_steps_per_rev': int(state.build.get('motor_full_steps_per_rev') or 1600),
        'motor_gear_ratio': int(state.build.get('motor_gear_ratio') or 5),
        'encoder_pulses_per_rev': int(state.build.get('encoder_pulses_per_rev') or 10000),
        'stim_channels': stim_ch,
        'motor_port': str(state.build.get('motor_port') or 'port0'),
        'encoder_chan': str(state.build.get('encoder_chan') or 'ctr0'),
        'SG_chan': sg_ch,
        'SG_name': sg_nm,
        'stim_monitor_chan': sm_ch,
        'stim_monitor_name': sm_nm,
        'S1side': str(state.build.get('S1side') or 'left'),
        'S2side': str(state.build.get('S2side') or 'right'),
        'use_sono': use_sono,
        'sono_channel': sono_ch if use_sono else [],
        'sono_name': sono_nm if use_sono else [],
        'sono_internal_samplefreq': int(state.build.get('sono_internal_samplefreq') or 241),
        'sono_cal_left': sono_lf,
        'sono_cal_right': sono_rf,
        'amp_step_vel': int(state.build.get('amp_step_vel') or 10),
        'ramp_mode_default': rm,
        'waitbefore': float(state.build.get('waitbefore') or 3.0),
        'waitafter': float(state.build.get('waitafter') or 4.0),
        'rampdur': float(state.build.get('rampdur') or 0.25),
        'prepoststim_dur': float(state.build.get('prepoststim_dur') or 0.06),
        'prepoststim_sep': float(state.build.get('prepoststim_sep') or 1.0),
        'prestim_time': float(state.build.get('prestim_time') or -2.0),
        'poststim_time': float(state.build.get('poststim_time') or 2.0),
    }
    src = render_generated_config(base, assignments)
    with open(path, 'w', encoding='utf-8') as f:
        f.write(src)
    importlib.invalidate_caches()
    err = apply_loaded_config_module(out_stem)
    if err:
        state.hardware_error = err
    else:
        state.cfg_mod = out_stem
        state.hardware_error = ''
        ui.notify(f'Wrote and loaded `{out_stem}`', type='positive')


@ui.refreshable
def hardware_section() -> None:
    ui.label('1 · Hardware configuration').classes(SECTION_TITLE)
    ui.markdown('**Build new** writes a new `.py`, then **Write config file and load**.').classes(CAPTION + ' mb-2')
    mods = discover_config_modules(_ROOT)
    with ui.row().classes('w-full gap-3'):
        ui.button('Load existing', on_click=lambda: _set_hw_mode('Load existing')).classes(
            BTN_MODE_ON if state.config_mode == 'Load existing' else BTN_MODE_OFF
        ).props('unelevated no-caps')
        ui.button('Build new', on_click=lambda: _set_hw_mode('Build new')).classes(
            BTN_MODE_ON if state.config_mode == 'Build new' else BTN_MODE_OFF
        ).props('unelevated no-caps')
    with ui.column().classes(PANEL + ' mt-3'):
        if state.config_mode == 'Load existing':
            ui.label('Hardware configuration module').classes(LABEL)
            ui.select(mods, value=state.load_cfg_select).classes('w-full').bind_value(state, 'load_cfg_select')
            ui.label('Override: type module name').classes(LABEL)
            ui.input(placeholder='e.g. my_lab_config').classes('w-full').bind_value(state, 'cfg_mod')
            ui.button('Load hardware configuration', on_click=_on_load_hardware).classes(BTN_COMMIT + ' mt-3').props(
                'unelevated no-caps'
            )
        else:
            seed_build_from_base()
            ui.label('Start from template (base module)').classes(LABEL)
            ui.select(mods, value=state.gui_cfg_build_base, on_change=_on_build_base_changed).classes('w-full')
            ui.label('Save new config as (module name, no `.py`)').classes(LABEL)
            ui.input(placeholder='e.g. lab_setup_2026').classes('w-full').bind_value(state, 'gui_cfg_build_out')
            with ui.expansion('Calibration, direction & axis labels', icon='tune').classes('w-full mt-2'):
                for k in (
                    'forcetorque_calibration_file',
                    'positive_motor_direction',
                    'specimen_lateral_index_on_positive_motor_side',
                    'motor_axis',
                    'bending_axis_sensor',
                    'primary_bending_axis',
                    'bending_axis_specimen',
                    'S1side',
                    'S2side',
                ):
                    if k == 'specimen_lateral_index_on_positive_motor_side':
                        ui.number(k.replace('_', ' '), format='%.0f').classes('w-full').bind_value(state.build, k)
                    else:
                        ui.input(k.replace('_', ' ')).classes('w-full').bind_value(state.build, k)
            with ui.expansion('DAQ rates & device', icon='memory', value=True).classes('w-full'):
                ui.input('NI-DAQ device name').classes('w-full').bind_value(state.build, 'device_name')
                ui.number('AI + encoder sample rate (Hz)', min=1.0, format='%.0f').classes('w-full').bind_value(
                    state.build, 'daq_ai_sample_rate_hz'
                )
                ui.number('AO + DO sample rate (Hz)', min=1.0, format='%.0f').classes('w-full').bind_value(
                    state.build, 'daq_ao_do_sample_rate_hz'
                )
            with ui.expansion('Motor, encoder & stim', icon='settings').classes('w-full'):
                ui.number('Motor full steps per rev', min=1, format='%.0f').classes('w-full').bind_value(
                    state.build, 'motor_full_steps_per_rev'
                )
                ui.number('Motor gear ratio', min=1, format='%.0f').classes('w-full').bind_value(
                    state.build, 'motor_gear_ratio'
                )
                ui.number('Encoder pulses per rev', min=1, format='%.0f').classes('w-full').bind_value(
                    state.build, 'encoder_pulses_per_rev'
                )
                ui.input('Motor port').classes('w-full').bind_value(state.build, 'motor_port')
                ui.input('Encoder channel').classes('w-full').bind_value(state.build, 'encoder_chan')
                ui.input('Stim AO channels (comma-separated)').classes('w-full').bind_value(state.build, 'stim_channels')
            with ui.expansion('Strain / ATI', icon='straighten').classes('w-full'):
                ui.input('SG AI channels').classes('w-full').bind_value(state.build, 'SG_chan')
                ui.input('SG names').classes('w-full').bind_value(state.build, 'SG_name')
            with ui.expansion('Sonomicrometry', icon='graphic_eq').classes('w-full'):
                ui.checkbox('Use sonomicrometry').bind_value(state.build, 'use_sono')
                ui.input('Sono channels').classes('w-full').bind_value(state.build, 'sono_channel')
                ui.input('Sono names').classes('w-full').bind_value(state.build, 'sono_name')
                ui.number('Sono internal sample frequency', min=1, format='%.0f').classes('w-full').bind_value(
                    state.build, 'sono_internal_samplefreq'
                )
                ui.input('Sono cal left').classes('w-full').bind_value(state.build, 'sono_cal_left')
                ui.input('Sono cal right').classes('w-full').bind_value(state.build, 'sono_cal_right')
            with ui.expansion('Stim monitor AI (optional)', icon='monitor_heart').classes('w-full'):
                ui.input('Stim monitor AI channels').classes('w-full').bind_value(state.build, 'stim_monitor_chan')
                ui.input('Stim monitor names').classes('w-full').bind_value(state.build, 'stim_monitor_name')
            with ui.expansion('Default timing, ramp & motion', icon='schedule').classes('w-full'):
                for k, mn in (
                    ('waitbefore', 0),
                    ('waitafter', 0),
                    ('rampdur', 0),
                    ('prepoststim_dur', 0),
                    ('prepoststim_sep', None),
                    ('prestim_time', None),
                    ('poststim_time', None),
                ):
                    if mn is not None:
                        ui.number(k, min=mn, format='%.6g').classes('w-full').bind_value(state.build, k)
                    else:
                        ui.number(k, format='%.6g').classes('w-full').bind_value(state.build, k)
                ui.number('amp_step_vel', min=1, format='%.0f').classes('w-full').bind_value(state.build, 'amp_step_vel')
                ui.select({'linear': 'linear', 'exponential': 'exponential'}).classes('w-full').bind_value(
                    state.build, 'ramp_mode_default'
                )
            ui.checkbox('Overwrite if file exists').bind_value(state, 'gui_cfg_build_overwrite').classes('mt-2')
            ui.button('Write config file and load', on_click=_on_write_hw).classes(BTN_COMMIT + ' mt-2').props(
                'unelevated no-caps'
            )
    if state.hardware_error:
        ui.markdown(state.hardware_error).classes(ALERT + ' mt-2 whitespace-pre-wrap')
    if state.bender:
        ui.markdown(f"Hardware active: `{getattr(state.bender, 'config_name', '?')}`").classes(SUCCESS + ' mt-2')


def _set_hw_mode(m: Literal['Load existing', 'Build new']) -> None:
    state.config_mode = m
    state.hardware_error = ''
    if m == 'Build new':
        state.build_seeded_for = ''
        seed_build_from_base()
    hardware_section.refresh()


def _on_build_base_changed(e: Any) -> None:
    state.gui_cfg_build_base = e.value
    state.build_seeded_for = ''
    seed_build_from_base()
    hardware_section.refresh()


def _on_load_hardware() -> None:
    state.hardware_error = ''
    eff = effective_load_module_name(typed=state.cfg_mod, selected=state.load_cfg_select)
    if not eff:
        state.hardware_error = 'Pick a module or type override.'
    else:
        err = apply_loaded_config_module(eff)
        state.hardware_error = err or ''
        if not err:
            ui.notify(f'Loaded `{normalize_config_module_name(eff)}`', type='positive')
    hardware_section.refresh()
    _refresh_workflow()


def _on_write_hw() -> None:
    state.hardware_error = ''
    maybe_build_write()
    hardware_section.refresh()
    if not state.hardware_error:
        _refresh_workflow()


def _refresh_workflow() -> None:
    data_path_section.refresh()
    biometrics_section.refresh()
    experiment_section.refresh()


@ui.refreshable
def data_path_section() -> None:
    ui.label('2 · Data file path').classes(SECTION_TITLE)
    if not state.bender:
        ui.markdown('Load **hardware** in section 1 first.').classes(ALERT)
        return
    ui.caption('HDF5 save location (separate from the hardware `.py` module).').classes(CAPTION)
    with ui.column().classes(PANEL):
        with ui.row().classes('w-full gap-3 flex-wrap'):
            ui.input('Data folder', placeholder=r'C:\Data\Experiments').classes('flex-1 min-w-[12rem]').bind_value(
                state.sess, 'gui_data_folder'
            ).props(f'title="{DATA_FOLDER_HELP}"')
            ui.input('Data file name', placeholder='my_experiment.h5').classes('flex-1 min-w-[12rem]').bind_value(
                state.sess, 'gui_data_filename'
            ).props(f'title="{DATA_FILE_NAME_HELP}"')
        p = compose_output_h5_path(state.sess)
        if p:
            ui.markdown(f'**Save path:** `{p}`').classes(CAPTION)
        ui.button('Load hardware configuration and data path', on_click=_on_load_hw_and_path).classes(BTN_COMMIT + ' mt-2').props(
            'unelevated no-caps'
        )


def _on_load_hw_and_path() -> None:
    if not state.bender:
        ui.notify('Load hardware first.', type='warning')
        return
    eff = effective_load_module_name(typed=state.cfg_mod, selected=state.load_cfg_select)
    if not eff:
        ui.notify('Resolve config module in section 1.', type='warning')
        return
    need_reload = not selected_config_matches_bender(state.bender, eff)
    if need_reload:
        err = apply_loaded_config_module(eff)
        if err:
            ui.notify(err, type='negative')
            _refresh_workflow()
            return
    perr = apply_data_path_to_bender(state.bender, state.sess)
    if perr:
        ui.notify(perr, type='warning')
    else:
        ui.notify('Hardware path and data file path applied.', type='positive')
    _refresh_workflow()


@ui.refreshable
def biometrics_section() -> None:
    ui.label('3 · Biometrics').classes(SECTION_TITLE)
    if not state.bender:
        ui.markdown('Load hardware first.').classes(ALERT)
        return
    init_biometrics_session_state(state.bender, state.sess, force=False)
    tpl_dir = shared_experiment_dir(state.sess, _ROOT)
    if not (str(state.sess.get('gui_data_folder') or '').strip() and os.path.isdir(os.path.normpath(str(state.sess.get('gui_data_folder'))))):
        ui.caption(f'Listing biometrics from `{tpl_dir}` until section 2 folder is valid.').classes(CAPTION)
    bio_opts: Dict[str, str] = {'': '— Choose a biometrics file —'}
    for p in list_biometrics_template_files(tpl_dir):
        bio_opts[p] = biometrics_template_display_label(p)
    if state.biometrics_feedback:
        ok, msg = state.biometrics_feedback
        state.biometrics_feedback = None
        ui.markdown(msg).classes(SUCCESS if ok else ALERT)
    ui.markdown('**Biometrics templates**').classes(LABEL)
    ui.select(bio_opts, value=str(state.sess.get('gui_biometrics_template_pick') or '')).classes('w-full').bind_value(
        state.sess, 'gui_biometrics_template_pick'
    )
    ui.button('Load biometrics into form', on_click=_on_load_biometrics).classes(BTN_COMMIT).props('unelevated no-caps')
    ui.input('Save biometrics as (name)').classes('w-full').bind_value(state.sess, 'gui_biometrics_new_name')
    ui.textarea('Description (optional)').classes('w-full').bind_value(state.sess, 'gui_biometrics_new_desc').props('rows=2')
    ui.checkbox('Overwrite if exists').bind_value(state.sess, 'gui_biometrics_overwrite')
    ui.button('Save biometrics', on_click=_on_save_biometrics).classes(BTN_COMMIT + ' mt-1').props('unelevated no-caps')

    with ui.column().classes(PANEL + ' mt-3'):
        ui.markdown('**Specimen identity**').classes(LABEL)
        with ui.row().classes('w-full gap-2 flex-wrap'):
            ui.input('Genus-species').classes('flex-1').bind_value(state.sess, 'gui_genus_species')
            ui.input('Specimen ID').classes('flex-1').bind_value(state.sess, 'gui_specimen_id')
        ui.input('Segment / preparation').classes('w-full').bind_value(state.sess, 'bio_segment')
        ui.button('Apply specimen identity', on_click=_on_apply_identity).classes(BTN_SIDE).props('unelevated no-caps')

        ui.markdown('**Intrinsic**').classes(LABEL + ' mt-3')
        with ui.row().classes('gap-2 w-full flex-wrap'):
            ui.number('TL', format='%.6g').classes('flex-1').bind_value(state.sess, 'bio_fishlen_TL')
            ui.number('SL', format='%.6g').classes('flex-1').bind_value(state.sess, 'bio_fishlen_SL')
        ui.number('Mass (g)', format='%.6g').classes('w-full').bind_value(state.sess, 'bio_fishmass')

        ui.markdown('**Experimental conditions**').classes(LABEL + ' mt-3')
        with ui.row().classes('gap-2 w-full'):
            ui.number('Room °C', format='%.3f').classes('flex-1').bind_value(state.sess, 'bio_temp_room')
            ui.number('Tank °C', format='%.3f').classes('flex-1').bind_value(state.sess, 'bio_temp_tank')
        ui.input('Prep condition').classes('w-full').bind_value(state.sess, 'bio_prep_condition')

        ui.markdown('**Clamp geometry**').classes(LABEL + ' mt-3')
        ui.number('Clamp spacing / segment length', format='%.6g').classes('w-full').bind_value(state.sess, 'bio_dclamp').props(
            f'title="{BIO_DBEND_FIELD_HELP}"'
        )
        ui.number('Along-body distance to segment center (mm)', format='%.6g').classes('w-full').bind_value(
            state.sess, 'bio_dbend'
        )
        with ui.row().classes('gap-2 w-full'):
            ui.number('Width', format='%.6g').classes('flex-1').bind_value(state.sess, 'bio_xsec')
            ui.number('Height', format='%.6g').classes('flex-1').bind_value(state.sess, 'bio_xsec_height')
        with ui.row().classes('gap-2 w-full'):
            ui.number('dvert', format='%.6g').classes('flex-1').bind_value(state.sess, 'bio_dvert')
            ui.number('dhoriz', format='%.6g').classes('flex-1').bind_value(state.sess, 'bio_dhoriz')

        ui.markdown('**Mounted profile**').classes(LABEL + ' mt-3')
        ui.select({x: x for x in BIO_DENSITY_PRESET_LABELS}, label='Density preset').classes('w-full').bind_value(
            state.sess, 'bio_prof_rho_preset'
        )
        ui.number('Specimen density (g/mm³)', format='%.6g', min=1e-12).classes('w-full').bind_value(
            state.sess, 'bio_prof_rho'
        )
        ui.number('Outline length (mm)', format='%.6g').classes('w-full').bind_value(state.sess, 'bio_prof_L')
        with ui.row().classes('gap-2 w-full flex-wrap'):
            ui.number('Prox H', format='%.6g').classes('flex-1').bind_value(state.sess, 'bio_prof_ph')
            ui.number('Prox W', format='%.6g').classes('flex-1').bind_value(state.sess, 'bio_prof_pw')
            ui.number('Dist H', format='%.6g').classes('flex-1').bind_value(state.sess, 'bio_prof_dh')
            ui.number('Dist W', format='%.6g').classes('flex-1').bind_value(state.sess, 'bio_prof_dw')
        ui.number('Clamp offset (mm)', format='%.6g').classes('w-full').bind_value(state.sess, 'bio_prof_clamp').props(
            f'title="{BIO_PROF_CLAMP_FIELD_HELP}"'
        )
        ui.number('Profile samples', min=20, max=400, step=10, format='%.0f').classes('w-full').bind_value(
            state.sess, 'bio_prof_samples'
        )
        ui.checkbox('Perform inertial correction').bind_value(state.sess, 'bio_use_theoretical_inertial')
        ui.button('Apply all biometrics', on_click=_on_apply_all_bio).classes(BTN_COMMIT + ' mt-3').props('unelevated no-caps')


def _on_load_biometrics() -> None:
    path = str(state.sess.get('gui_biometrics_template_pick') or '').strip()
    if not path:
        state.biometrics_feedback = (False, 'Choose a biometrics file first.')
    else:
        try:
            data = load_biometrics_template(path)
            ok, msg = apply_biometrics_template_to_session(state.sess, data)
            state.biometrics_feedback = (ok, msg)
        except Exception as e:
            state.biometrics_feedback = (False, f'{type(e).__name__}: {e}')
    biometrics_section.refresh()


def _on_save_biometrics() -> None:
    if not state.bender:
        return
    bn = str(state.sess.get('gui_biometrics_new_name') or '').strip()
    bd = str(state.sess.get('gui_biometrics_new_desc') or '').strip()
    bst = sanitize_biometrics_filename_stem(bn or 'biometrics')
    bout = os.path.normpath(os.path.join(shared_experiment_dir(state.sess, _ROOT), f'{bst}.json'))
    try:
        if os.path.isfile(bout) and not bool(state.sess.get('gui_biometrics_overwrite')):
            state.biometrics_feedback = (False, f'File exists: `{bout}`')
        else:
            os.makedirs(os.path.dirname(bout) or '.', exist_ok=True)
            save_biometrics_template(bout, name=bn or bst, description=bd, session_state=state.sess)
            state.biometrics_feedback = (True, f'Saved `{bout}`')
    except Exception as e:
        state.biometrics_feedback = (False, f'{type(e).__name__}: {e}')
    biometrics_section.refresh()


def _on_apply_identity() -> None:
    if state.bender:
        apply_specimen_identity_to_bender(state.bender, state.sess)
        ui.notify('Specimen identity applied.', type='positive')


def _on_apply_all_bio() -> None:
    if state.bender:
        from .session_logic import sync_bio_prof_rho_from_density_preset

        sync_bio_prof_rho_from_density_preset(state.sess)
        apply_all_biometrics_to_bender(state.bender, state.sess)
        ui.notify('Biometrics applied.', type='positive')


@ui.refreshable
def procedure_panel() -> None:
    b = state.bender
    if not b:
        return
    schema = b.get_dispatch_schema()
    test_types = list(schema['test_types'])
    tt = str(state.sess.get('test_type_select') or b.test_type or '')
    if tt not in test_types:
        tt = test_types[0]
        state.sess['test_type_select'] = tt
    with ui.expansion('Procedure fields', value=not bool(state.sess.get('gui_exp_hide'))).classes('w-full'):
        render_procedure_fields(
            b, tt, schema, test_types, state.sess, caption_class=CAPTION, label_class=LABEL
        )
        ui.separator().classes('my-3')
        ui.markdown('**Save procedure as template**').classes(LABEL)
        ui.input('Template name').classes('w-full').bind_value(state.sess, 'gui_protocol_new_name')
        ui.textarea('Description').classes('w-full').bind_value(state.sess, 'gui_protocol_new_desc').props('rows=2')
        ui.checkbox('Overwrite if exists').bind_value(state.sess, 'gui_protocol_overwrite')
        with ui.row().classes('w-full gap-2 mt-2'):
            ui.button('Apply procedure', on_click=_on_apply_procedure).classes(BTN_SIDE + ' flex-1').props('unelevated no-caps')
            ui.button('Save template', on_click=_on_save_protocol).classes(BTN_SIDE + ' flex-1').props('unelevated no-caps')
    ui.checkbox('Hide procedure section').bind_value(state.sess, 'gui_exp_hide')


def _current_updates() -> Tuple[Dict[str, Any], str | None]:
    b = state.bender
    assert b is not None
    schema = b.get_dispatch_schema()
    tt = str(state.sess.get('test_type_select') or b.test_type)
    return gather_procedure_updates(b, tt, schema, list(schema['test_types']), state.sess)


def _on_apply_procedure() -> None:
    if not state.bender:
        return
    updates, err = _current_updates()
    if err:
        ui.notify(f'Procedure: {err}', type='negative')
        return
    tt = str(state.sess.get('test_type_select') or state.bender.test_type)
    apply_procedure_form_to_bender(state.bender, updates, tt, state.sess)
    ui.notify('Procedure applied.', type='positive')


def _on_save_protocol() -> None:
    if not state.bender:
        return
    b = state.bender
    updates, err = _current_updates()
    if err:
        ui.notify(err, type='negative')
        return
    tt = str(state.sess.get('test_type_select') or b.test_type)
    name = str(state.sess.get('gui_protocol_new_name') or '').strip()
    desc = str(state.sess.get('gui_protocol_new_desc') or '').strip()
    stem = sanitize_template_filename_stem(name or 'protocol')
    out = os.path.normpath(os.path.join(shared_experiment_dir(state.sess, _ROOT), f'{stem}.json'))
    schema = b.get_dispatch_schema()
    try:
        if os.path.isfile(out) and not bool(state.sess.get('gui_protocol_overwrite')):
            ui.notify(f'File exists: {out}', type='warning')
            return
        os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
        proc = build_procedure_dict_from_updates(updates)
        base = None
        if tt == 'calibration':
            btt = proc.get('calibration_base_test_type')
            if btt:
                snap = snapshot_bender_procedure(b, schema, str(btt))
                if snap:
                    base = {'test_type': str(btt), 'procedure': snap}
        save_protocol_template(out, name=name or stem, description=desc, test_type=tt, procedure=proc, base_protocol=base)
        ui.notify(f'Saved `{out}`', type='positive')
    except Exception as e:
        ui.notify(f'{type(e).__name__}: {e}', type='negative')


@ui.refreshable
def experiment_section() -> None:
    if not state.bender:
        ui.markdown('Load hardware to configure experiment.').classes(ALERT)
        return
    b = state.bender
    schema = b.get_dispatch_schema()
    test_types = list(schema['test_types'])
    cur_tt = str(state.sess.get('test_type_select') or b.test_type or '')
    if cur_tt not in test_types:
        cur_tt = test_types[0]
    state.sess['test_type_select'] = cur_tt
    b.test_type = cur_tt

    ui.label('4 · Experiment type & parameters').classes(SECTION_TITLE)
    if state.protocol_feedback:
        ok, msg = state.protocol_feedback
        state.protocol_feedback = None
        ui.markdown(msg).classes(SUCCESS if ok else ALERT)

    tpl_dir = shared_experiment_dir(state.sess, _ROOT)
    opts: Dict[str, str] = {'': '— Choose a template —'}
    for p in list_template_files(tpl_dir):
        opts[p] = protocol_template_option_label(p)
    pick = str(state.sess.get('gui_protocol_template_pick') or '')
    if pick not in opts:
        pick = ''
    ui.markdown('**Protocol templates (load)**').classes(LABEL)
    ui.select(opts, value=pick).classes('w-full').bind_value(state.sess, 'gui_protocol_template_pick')
    ui.button('Load template into form', on_click=_on_load_protocol_template).classes(BTN_COMMIT + ' mt-1').props(
        'unelevated no-caps'
    )

    def _tt_change(e: Any) -> None:
        state.sess['test_type_select'] = e.value
        b.test_type = e.value
        procedure_panel.refresh()

    ui.label('Experiment type (`test_type`)').classes(LABEL + ' mt-2')
    ui.select({t: t for t in test_types}, value=cur_tt, on_change=_tt_change).classes('w-full')

    procedure_panel()

    ui.label('5 · Experiment preview (no DAQ)').classes(SECTION_TITLE + ' mt-6')
    with ui.expansion('Preview', value=not bool(state.sess.get('gui_sec4_hide'))).classes('w-full'):
        ui.button('Apply (procedure + biometrics flags)', on_click=_on_preview_apply).classes(BTN_COMMIT + ' mb-2').props(
            'unelevated no-caps'
        )
        ui.checkbox('Show last preview').bind_value(state.sess, 'gui_show_preview')
        ui.slider(min=400, max=12000, step=200, value=state.sess.get('gui_preview_pts', 6000)).classes('w-full').bind_value(
            state.sess, 'gui_preview_pts'
        )
        ui.button('Refresh preview', on_click=_on_refresh_preview).classes(BTN_COMMIT + ' mt-2').props('unelevated no-caps')
        prev = state.sess.get('gui_last_preview')
        if state.sess.get('gui_show_preview') and isinstance(prev, dict):
            if prev.get('error'):
                ui.markdown(str(prev['error'])).classes(ALERT)
            elif prev.get('ok') and prev.get('table'):
                ui.markdown('**Summary**')
                ui.html(pd.DataFrame(prev['table']).to_html(classes='text-sm', index=False), sanitize=False)
            tp, ap, vp = prev.get('t_plot'), prev.get('angle_plot'), prev.get('anglevel_plot')
            if tp is not None and ap is not None and len(tp) > 0:
                fig = go.Figure()
                fig.add_trace(go.Scatter(x=tp, y=ap, mode='lines', name='angle (deg)'))
                if vp is not None and len(vp) > 0:
                    fig.add_trace(go.Scatter(x=tp, y=vp, mode='lines', name='anglevel', yaxis='y2'))
                    fig.update_layout(yaxis2=dict(title='deg/s', overlaying='y', side='right'), height=420)
                else:
                    fig.update_layout(height=420)
                ui.plotly(fig).classes('w-full')
    ui.checkbox('Hide preview section').bind_value(state.sess, 'gui_sec4_hide')

    ui.label('6 · Save, validate, and run').classes(SECTION_TITLE + ' mt-6')
    with ui.expansion('Run', value=not bool(state.sess.get('gui_sec5_hide'))).classes('w-full'):
        ui.button('View current settings', on_click=_on_view_settings).classes(BTN_COMMIT + ' mb-2').props('unelevated no-caps')
        with ui.row().classes('gap-2 w-full'):
            ui.button('Apply procedure', on_click=_on_apply_procedure).classes(BTN_SIDE + ' flex-1').props('unelevated no-caps')
            ui.button('Check required fields', on_click=_on_validate).classes(BTN_SIDE + ' flex-1').props('unelevated no-caps')
        ui.checkbox('Hardware: I intend to run DAQ').bind_value(state.sess, 'gui_run_daq_ok')
        ui.checkbox('Biometrics applied (section 3)').bind_value(state.sess, 'gui_run_biometrics_confirm')
        needs_cal = needs_missing_calibration_confirmation(b)
        if needs_cal:
            ui.markdown('No calibration file detected.').classes(ALERT)
        cb_cal = ui.checkbox('Proceed without calibration file').bind_value(
            state.sess, 'gui_confirm_run_without_calibration'
        )
        if not needs_cal:
            cb_cal.disable()
        needs_dest = section2_destination_incomplete(state.sess)
        if needs_dest:
            ui.markdown('Section 2 path incomplete.').classes(ALERT)
        cb_dest = ui.checkbox('Proceed without section 2 save path').bind_value(
            state.sess, 'gui_confirm_run_without_destination'
        )
        if not needs_dest:
            cb_dest.disable()
        ui.checkbox('Append QC notes to H5').bind_value(state.sess, 'gui_qc_notes_append')
        ui.button('Run experiment', on_click=_on_run_experiment).classes(BTN_COMMIT + ' mt-3').props('unelevated no-caps')
    ui.checkbox('Hide run section').bind_value(state.sess, 'gui_sec5_hide')
    ui.separator().classes('my-4')
    ui.input('Post-trial / QC notes').classes('w-full').bind_value(state.sess, 'gui_post_notes')


def _on_load_protocol_template() -> None:
    if not state.bender:
        return
    path = str(state.sess.get('gui_protocol_template_pick') or '').strip()
    if not path:
        state.protocol_feedback = (False, 'Choose a template first.')
    else:
        try:
            data = load_protocol_template(path)
            schema = state.bender.get_dispatch_schema()
            ok, msg = apply_template_to_session_state(
                state.sess, data, widget_key=proc_widget_key, valid_test_types=list(schema['test_types'])
            )
            state.protocol_feedback = (ok, msg)
            if ok:
                tt = str(state.sess.get('test_type_select') or '')
                if tt:
                    state.bender.test_type = tt
                ui.notify(msg, type='positive')
        except Exception as e:
            state.protocol_feedback = (False, f'{type(e).__name__}: {e}')
    experiment_section.refresh()


def _on_preview_apply() -> None:
    if not state.bender:
        return
    updates, err = _current_updates()
    if err:
        ui.notify(err, type='negative')
        return
    tt = str(state.sess.get('test_type_select') or state.bender.test_type)
    sync_biometric_flags_from_session(state.bender, state.sess)
    sync_genus_species_to_bender(state.bender, state.sess)
    apply_procedure_form_to_bender(state.bender, updates, tt, state.sess)
    ui.notify('Applied.', type='positive')


def _on_refresh_preview() -> None:
    if not state.bender:
        return
    updates, err = _current_updates()
    if err:
        ui.notify(err, type='negative')
        return
    tt = str(state.sess.get('test_type_select') or state.bender.test_type)
    sync_biometric_flags_from_session(state.bender, state.sess)
    sync_genus_species_to_bender(state.bender, state.sess)
    apply_form_updates(state.bender, updates, tt, state.sess)
    pts = int(state.sess.get('gui_preview_pts') or 6000)
    state.sess['gui_last_preview'] = build_protocol_preview(state.bender, requested_test_type=tt, max_plot_points=pts)
    state.sess['gui_last_preview_tt'] = tt
    experiment_section.refresh()
    ui.notify('Preview updated.', type='positive')


def _on_view_settings() -> None:
    if not state.bender:
        return
    updates, err = _current_updates()
    if err:
        ui.notify(err, type='negative')
        return
    tt = str(state.sess.get('test_type_select') or state.bender.test_type)
    sync_biometric_flags_from_session(state.bender, state.sess)
    sync_genus_species_to_bender(state.bender, state.sess)
    apply_form_updates(state.bender, updates, tt, state.sess)
    rows = [
        {'group': 'experiment', 'name': 'test_type', 'value': tt},
        {'group': 'export', 'name': 'h5', 'value': compose_output_h5_path(state.sess) or getattr(state.bender, 'outputfile', '')},
    ]
    for k, v in sorted(updates.items(), key=lambda kv: kv[0]):
        rows.append({'group': 'parameter', 'name': k, 'value': str(v)})
    with ui.dialog() as d, ui.card():
        ui.html(pd.DataFrame(rows).to_html(classes='text-xs', index=False), sanitize=False)
        ui.button('Close', on_click=d.close)
    d.open()


def _on_validate() -> None:
    if not state.bender:
        return
    updates, err = _current_updates()
    if err:
        ui.notify(err, type='negative')
        return
    tt = str(state.sess.get('test_type_select') or state.bender.test_type)
    sync_biometric_flags_from_session(state.bender, state.sess)
    sync_genus_species_to_bender(state.bender, state.sess)
    apply_form_updates(state.bender, updates, tt, state.sess)
    rep = state.bender.validate_dispatch_setup(test_type=tt)
    if rep['ok']:
        ui.notify('All required fields set.', type='positive')
    else:
        ui.notify(f"Missing: {', '.join(rep['missing'][:8])}", type='warning')


def _on_run_experiment() -> None:
    if not state.bender:
        return
    if not state.sess.get('gui_run_daq_ok'):
        ui.notify('Confirm DAQ intent.', type='warning')
        return
    if not state.sess.get('gui_run_biometrics_confirm'):
        ui.notify('Confirm biometrics applied.', type='warning')
        return
    needs_cal = needs_missing_calibration_confirmation(state.bender)
    if needs_cal and not state.sess.get('gui_confirm_run_without_calibration'):
        ui.notify('Calibration confirmation required.', type='warning')
        return
    needs_dest = section2_destination_incomplete(state.sess)
    if needs_dest and not state.sess.get('gui_confirm_run_without_destination'):
        ui.notify('Set section 2 or confirm proceed without path.', type='warning')
        return
    updates, err = _current_updates()
    if err:
        ui.notify(err, type='negative')
        return
    tt = str(state.sess.get('test_type_select') or state.bender.test_type)
    b = state.bender
    sync_biometric_flags_from_session(b, state.sess)
    apply_form_updates(b, updates, tt, state.sess)
    sync_genus_species_to_bender(b, state.sess)
    outp = compose_output_h5_path(state.sess).strip()
    if outp:
        b.outputfile = outp
    notes_in = str(state.sess.get('gui_post_notes') or '').strip()
    try:
        b.run_experiment(test_type=tt)
        rep = export_primary_h5(
            b,
            post_trial_notes=notes_in if notes_in else None,
            outputfile=outp or None,
            append_post_trial_notes=bool(state.sess.get('gui_qc_notes_append', True)),
        )
        ntr = len(getattr(b, 'trial_records', []) or [])
        qix = read_qc_trial_index(state.sess, ntr)
        sel_h5 = str(rep.get('outputfile') or outp or '')
        qc_base = qc_figure_base_path(b, sel_h5, qix)
        qc_path, _ = save_universal_qc_figure(b, qc_trial_index=qix, base_path=qc_base)
        ui.notify(f"Saved `{rep['outputfile']}` | QC `{qc_path}`", type='positive')
    except Exception as e:
        ui.notify(f'Run failed: {type(e).__name__}: {e}', type='negative')


@ui.refreshable
def h5_explorer_section() -> None:
    ui.label('H5 data explorer').classes(SECTION_TITLE)
    ui.caption('Paths are read on this machine.').classes(CAPTION)
    ui.input('Path to `.h5`').classes('w-full').bind_value(state.sess, 'gui_h5_explore_path')
    ui.button('Load file & refresh series', on_click=_on_h5_load).classes(BTN_COMMIT + ' mt-2').props('unelevated no-caps')

    loaded = str(state.sess.get('gui_h5_explore_loaded_path') or '').strip()
    if not loaded or not os.path.isfile(loaded):
        if loaded:
            ui.markdown('File not found.').classes(ALERT)
        return
    if str(state.sess.get('gui_h5_explore_schema') or '') == 'v2':
        trials = list_v2_trials(loaded)

        def _trial_changed(e: Any) -> None:
            state.sess['gui_h5_explore_trial'] = e.value
            try:
                cat, n2 = build_series_catalog_v2(loaded, str(e.value))
                state.sess['gui_h5_explore_catalog'] = cat
                state.sess['gui_h5_explore_notes'] = n2
            except Exception as ex:
                ui.notify(str(ex), type='negative')
            h5_explorer_section.refresh()

        if trials:
            cur = str(state.sess.get('gui_h5_explore_trial') or trials[0])
            if cur not in trials:
                cur = trials[0]
            ui.select({t: t for t in trials}, label='Trial', value=cur, on_change=_trial_changed).classes('w-full mb-2')
    cat = state.sess.get('gui_h5_explore_catalog') or {}
    if not cat:
        return
    keys = sorted(cat.keys())
    ui.select({k: k for k in keys}, label='X').classes('w-full').bind_value(state.sess, 'gui_h5_explore_x')
    ui.select({k: k for k in keys}, label='Y').classes('w-full').bind_value(state.sess, 'gui_h5_explore_y')
    try:
        xk = str(state.sess.get('gui_h5_explore_x') or keys[0])
        yk = str(state.sess.get('gui_h5_explore_y') or keys[min(1, len(keys) - 1)])
        xd, yd, n = align_xy(cat, xk, yk)
        if n > 0:
            fig = go.Figure(go.Scatter(x=xd, y=yd, mode='lines', name=yk))
            fig.update_layout(height=480, xaxis_title=xk, yaxis_title=yk)
            ui.plotly(fig).classes('w-full')
    except Exception as e:
        ui.markdown(str(e)).classes(ALERT)


def _on_h5_load() -> None:
    pin = str(state.sess.get('gui_h5_explore_path') or '').strip()
    upl = str(state.sess.get('gui_h5_explore_upload_path') or '').strip()
    chosen = ''
    if pin and os.path.isfile(pin):
        chosen = os.path.normpath(pin)
    elif upl and os.path.isfile(upl):
        chosen = upl
    if not chosen:
        ui.notify('No valid file.', type='warning')
        return
    try:
        schema = detect_h5_schema(chosen)
        notes: List[str] = []
        if schema == 'v2':
            trials = list_v2_trials(chosen)
            if not trials:
                ui.notify('No trials in v2 file.', type='negative')
                return
            tid = str(state.sess.get('gui_h5_explore_trial') or trials[0])
            if tid not in trials:
                tid = trials[0]
            state.sess['gui_h5_explore_trial'] = tid
            cat, n2 = build_series_catalog_v2(chosen, tid)
            notes.extend(n2)
        elif schema == 'legacy':
            state.sess['gui_h5_explore_trial'] = ''
            cat, n2 = build_series_catalog_legacy(chosen)
            notes.extend(n2)
        else:
            ui.notify('Unknown H5 schema.', type='negative')
            return
        state.sess['gui_h5_explore_loaded_path'] = chosen
        state.sess['gui_h5_explore_schema'] = schema
        state.sess['gui_h5_explore_catalog'] = cat
        state.sess['gui_h5_explore_notes'] = notes
        ui.notify(f'Loaded {os.path.basename(chosen)}', type='positive')
    except Exception as e:
        ui.notify(str(e), type='negative')
    h5_explorer_section.refresh()


def _on_daq_kill() -> None:
    if not state.bender:
        ui.notify('No hardware loaded.', type='warning')
        return
    dev = str(getattr(state.bender, 'device_name', '') or '').strip()
    if not dev:
        ui.notify('No device_name on config.', type='warning')
        return
    ok, msg = daq_emergency_stop(dev)
    ui.notify(msg, type='positive' if ok else 'negative')


def _go(route: Literal['landing', 'workflow', 'h5']) -> None:
    _set_client_route(route)
    shell_body.refresh()


@ui.refreshable
def shell_body() -> None:
    ui.query('body').classes(SHELL_BODY)
    with ui.header(elevated=True).classes('items-center justify-between bg-white border-b border-slate-200'):
        ui.label('CritterGripper').classes('text-lg font-semibold text-slate-700')
        with ui.row().classes('gap-2'):
            ui.button('Workflow', on_click=lambda: _go('workflow')).props('flat no-caps')
            ui.button('H5 explorer', on_click=lambda: _go('h5')).props('flat no-caps')
            ui.button('Landing', on_click=lambda: _go('landing')).props('flat no-caps')
            ui.button('DAQ stop', on_click=_on_daq_kill, color='red').props('flat no-caps')
    with ui.column().classes(MAIN_COLUMN):
        ui.link('Skip to main content', '#main-content').classes(SKIP_LINK)
        route = _client_route()
        if route == 'landing':
            with ui.element('div').props('id=main-content').classes('w-full flex flex-col items-center py-8 md:py-12'):
                ui.label('CritterGripper').classes(
                    'text-3xl md:text-4xl font-semibold text-slate-700 tracking-tight text-center mt-2'
                )
                ui.markdown(
                    'Industrial / scientific console for **hardware**, **biometrics**, and **protocol** — '
                    'preview commanded motion without DAQ, validate, then run when you are ready.'
                ).classes('text-slate-500 text-center max-w-xl text-base mt-3 mb-2')
                ui.markdown('Choose where to start:').classes('text-sm text-slate-500 mt-6 mb-2')
                with ui.row().classes('w-full max-w-2xl gap-4 justify-center flex-wrap'):
                    ui.button('Enter workflow', on_click=lambda: _go('workflow')).classes(
                        BTN_COMMIT + ' flex-1 min-w-[10rem] min-h-[3rem]'
                    ).props('unelevated no-caps')
                    ui.button('H5 data explorer', on_click=lambda: _go('h5')).classes(
                        BTN_COMMIT + ' flex-1 min-w-[10rem] min-h-[3rem]'
                    ).props('unelevated no-caps')
                ui.markdown('Use the header anytime to switch views. **DAQ stop** halts NI outputs on the loaded device.').classes(
                    'text-slate-500 text-center max-w-xl text-sm mt-8'
                )
        elif route == 'h5':
            with ui.element('div').props('id=main-content').classes('w-full'):
                h5_explorer_section()
        else:
            with ui.element('div').props('id=main-content').classes('w-full space-y-6'):
                ensure_gui_defaults(state.sess)
                _init_hw_defaults()
                hardware_section()
                data_path_section()
                biometrics_section()
                experiment_section()


@ui.page('/')
def page_main() -> None:
    # Each browser tab gets its own route; default is landing (fixes shared module state across clients).
    app.storage.client.setdefault('gui_route', 'landing')
    shell_body()


def run_app(*, port: int = 8765) -> None:
    ui.run(title='CritterGripper', port=port, favicon='🦴', reload=True, show=True, tailwind=True)
