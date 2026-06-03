import os
from pathlib import Path
import sys
import tempfile

import h5py
import numpy as np
import streamlit as st
from streamlit.testing.v1 import AppTest

sys.path.append(str(Path(__file__).resolve().parents[1]))

import bender_h5_explore as h5x
import bender_json_persistent as jp
import bender_streamlit_gui as gui


def test_h5_backend_handles_invalid_paths_and_types():
    assert h5x.detect_h5_schema('') == 'unknown'
    assert h5x.detect_h5_schema('not_a_real_file.h5') == 'unknown'
    assert h5x.list_h5_group_children('not_a_real_file.h5', '') == []
    assert h5x.list_h5_group_children('not_a_real_file.h5', '../bad') == []


def test_h5_backend_roundtrip_generic_file_and_attribute_writes():
    fd, path = tempfile.mkstemp(suffix='.h5')
    os.close(fd)
    try:
        with h5py.File(path, 'w') as f:
            g = f.create_group('g')
            g.create_dataset('x', data=np.linspace(0.0, 1.0, 50))
            g.create_dataset('ft', data=np.zeros((6, 50), dtype=np.float64))
            g.attrs['species'] = 'test'

        assert h5x.detect_h5_schema(path) in ('generic', 'browse')
        rows = h5x.list_h5_group_children(path, 'g')
        assert any(r['name'] == 'x' and r['plot'] == '1d' for r in rows)
        assert any(r['name'] == 'ft' and r['plot'] == 'six' for r in rows)

        x = h5x.read_h5_series_1d(path, 'g/x')
        tz = h5x.read_h5_series_1d(path, 'g/ft', channel=5)
        assert x.size == 50
        assert tz.size == 50

        notes = h5x.write_h5_user_attributes(
            path,
            'g',
            updates={'species': ('str', 'trout')},
            delete_names=[],
            additions=[('valid_bool', 'bool', 'true')],
        )
        assert notes
    finally:
        os.remove(path)


def test_h5_backend_rejects_invalid_paths_and_channels():
    fd, path = tempfile.mkstemp(suffix='.h5')
    os.close(fd)
    try:
        with h5py.File(path, 'w') as f:
            f.create_dataset('ft', data=np.zeros((6, 10), dtype=np.float64))

        try:
            h5x.read_h5_series_1d(path, 'ft')
            assert False, 'expected channel-required failure'
        except ValueError:
            pass

        try:
            h5x.read_h5_series_1d(path, 'ft', channel=8)
            assert False, 'expected invalid channel failure'
        except ValueError:
            pass

        try:
            h5x.read_h5_series_1d(path, '../ft', channel=0)
            assert False, 'expected invalid path failure'
        except ValueError:
            pass
    finally:
        os.remove(path)


def test_streamlit_routes_render_with_stressed_session_state():
    routes = ('scratch', 'templates', 'stepwise', 'h5_explorer')
    for route in routes:
        at = AppTest.from_file('bender_streamlit_gui.py')
        at.session_state['gui_app_route'] = route

        # Stress key areas: config/path/morphometrics/protocol-related session values.
        at.session_state['gui_cfg_mod'] = 'nonexistent_config'
        at.session_state['gui_data_folder'] = r'Z:\definitely\not\real\folder'
        at.session_state['gui_data_filename'] = ''
        at.session_state['morpho_fishmass'] = -1.0
        at.session_state['morpho_dclamp'] = -50.0
        at.session_state['morpho_fishlen_TL'] = np.nan
        # Stress the specimen-geometry parser inputs (mismatched / invalid lists).
        at.session_state['morpho_geom_x'] = '1, 5, 7'
        at.session_state['morpho_geom_y'] = '3, 2'
        at.session_state['morpho_geom_pos'] = 'abc, 0, 10'
        at.session_state['gui_run_morphometrics_confirm'] = False
        at.session_state['gui_confirm_run_without_calibration'] = False
        at.session_state['gui_confirm_run_without_destination'] = False
        at.session_state['gui_h5_explore_path'] = r'C:\not_real\fake.h5'

        at.run(timeout=120)
        assert not at.exception


def _clear_streamlit_session_state():
    for k in list(st.session_state.keys()):
        st.session_state.pop(k, None)


def test_autosave_roundtrip_and_start_fresh_cleanup(tmp_path, monkeypatch):
    _clear_streamlit_session_state()
    monkeypatch.setattr(gui, '_ROOT', str(tmp_path))
    st.session_state['gui_app_route'] = 'scratch'
    st.session_state['gui_data_folder'] = str(tmp_path)
    st.session_state['gui_data_filename'] = 'trial.h5'
    st.session_state['morpho_fishmass'] = 10.0
    st.session_state['gui_setup_confirmed'] = True
    st.session_state['gui_measurements_confirmed'] = True
    st.session_state['gui_protocol_confirmed'] = True
    st.session_state['gui_default_state_baseline'] = {}
    st.session_state['gui_recovered_state_baseline'] = {}

    ok, saved_path = gui._write_snapshot_payload(source='autosave', update_latest=True)
    assert ok
    assert os.path.isfile(saved_path)
    assert os.path.isfile(gui._autosave_latest_path())

    payload, err = gui._load_latest_autosave()
    assert err is None
    assert payload is not None
    assert payload['state']['gui_data_filename'] == 'trial.h5'
    assert 'gui_setup_confirmed' not in payload['state']
    assert 'gui_measurements_confirmed' not in payload['state']
    assert 'gui_protocol_confirmed' not in payload['state']

    _clear_streamlit_session_state()
    gui._restore_autosave_payload(payload)
    assert st.session_state['gui_data_filename'] == 'trial.h5'
    assert st.session_state['morpho_fishmass'] == 10.0
    assert st.session_state['gui_session_source'] == 'restored'

    gui._reset_workflow_session_to_home(clear_autosave=True)
    assert st.session_state.get('gui_app_route') == 'landing'
    assert st.session_state.get('gui_session_source') == 'fresh'
    assert not os.path.isfile(gui._autosave_latest_path())


def test_autosave_excludes_ndarray_caches(tmp_path, monkeypatch):
    _clear_streamlit_session_state()
    monkeypatch.setattr(gui, '_ROOT', str(tmp_path))
    st.session_state['gui_app_route'] = 'scratch'
    st.session_state['gui_data_folder'] = str(tmp_path)
    st.session_state['gui_last_preview'] = {
        'ok': True,
        't': np.array([0.0, 1.0]),
        'angle': np.array([0.0, 5.0]),
        't_plot': np.array([0.0, 1.0]),
    }
    st.session_state['gui_h5_explore_catalog'] = {'Time (s)': np.linspace(0.0, 1.0, 4)}

    ok, _ = gui._write_snapshot_payload(source='autosave', update_latest=True)
    assert ok
    payload, err = gui._load_latest_autosave()
    assert err is None
    assert payload is not None
    state = payload['state']
    assert 'gui_last_preview' not in state
    assert 'gui_h5_explore_catalog' not in state


def test_autosave_rejects_non_persistable_gui_value(tmp_path, monkeypatch):
    _clear_streamlit_session_state()
    monkeypatch.setattr(gui, '_ROOT', str(tmp_path))
    st.session_state['gui_app_route'] = 'scratch'
    st.session_state['gui_data_folder'] = str(tmp_path)

    class _Opaque:
        pass

    st.session_state['gui_opaque_test_holder'] = _Opaque()

    ok, msg = gui._write_snapshot_payload(source='autosave', update_latest=True)
    assert not ok
    assert 'JsonPersistTypeError' in msg or 'Non-JSON-persistable' in msg
    assert st.session_state.get('gui_autosave_last_error')


def test_to_json_persistent_never_stringifies_ndarray():
    arr = np.array([-3.0, -2.998])
    out = jp.to_json_persistent(arr)
    assert isinstance(out, list)
    assert out == [-3.0, -2.998]
    assert not isinstance(out, str)


def test_autosave_corrupt_and_schema_mismatch_fallback(tmp_path, monkeypatch):
    _clear_streamlit_session_state()
    monkeypatch.setattr(gui, '_ROOT', str(tmp_path))
    os.makedirs(gui._session_snapshots_dir(), exist_ok=True)

    with open(gui._autosave_latest_path(), 'w', encoding='utf-8') as f:
        f.write('{bad json')
    gui._bootstrap_autosave_recovery()
    assert 'Autosave ignored' in str(st.session_state.get('gui_recovery_banner_message') or '')
    assert st.session_state.get('gui_recovery_banner_level') == 'warning'

    _clear_streamlit_session_state()
    with open(gui._autosave_latest_path(), 'w', encoding='utf-8') as f:
        f.write('{"schema_version":999,"state":{}}')
    gui._bootstrap_autosave_recovery()
    assert 'Autosave ignored' in str(st.session_state.get('gui_recovery_banner_message') or '')


def test_state_origin_map_classifies_default_recovered_and_user():
    _clear_streamlit_session_state()
    st.session_state['gui_default_state_baseline'] = {'gui_a': 1, 'gui_b': 2}
    st.session_state['gui_recovered_state_baseline'] = {'gui_b': 9}
    origin = gui._build_state_origin_map({'gui_a': 1, 'gui_b': 9, 'gui_c': 3})
    assert origin['gui_a'] == 'default'
    assert origin['gui_b'] == 'recovered'
    assert origin['gui_c'] == 'user'


class _DummyBender:
    def __init__(self, config_name: str = 'dummy_cfg'):
        self.config_name = config_name
        self.outputfile = ''

    def validate_dispatch_setup(self, test_type=None):
        return {'ok': True, 'missing': [], 'test_type': str(test_type or 'dynamic')}


def test_checklist_requires_confirmation_plus_valid_state():
    _clear_streamlit_session_state()
    b = _DummyBender('dummy_cfg')
    st.session_state['bender'] = b
    st.session_state['gui_load_cfg_select'] = 'dummy_cfg'
    st.session_state['gui_config_setup_mode'] = 'Load existing'
    st.session_state['gui_data_folder'] = r'C:\tmp'
    st.session_state['gui_data_filename'] = 'x.h5'
    b.outputfile = gui._compose_output_h5_path()
    gui._mark_data_path_applied()

    ready = gui._workflow_ready_state(b, 'dynamic')
    assert ready['setup_ok']
    assert not ready['measurements_ok']
    assert ready['run_disabled'] is False

    st.session_state['gui_genus_species'] = 'Danio rerio'
    st.session_state['gui_specimen_id'] = 'fish-1'
    st.session_state['morpho_fishmass'] = 2.0
    st.session_state['morpho_dclamp'] = 10.0
    st.session_state['morpho_xsec'] = 2.0
    gui._mark_morpho_applied()
    ready2 = gui._workflow_ready_state(b, 'dynamic')
    assert ready2['measurements_ok']


def test_run_button_state_when_no_bender():
    _clear_streamlit_session_state()
    disabled, reason = gui._run_button_state()
    assert disabled
    assert 'Setup' in reason


def test_protocol_checklist_ignores_form_widget_status_section():
    _clear_streamlit_session_state()
    b = _DummyBender('dummy_cfg')
    st.session_state['bender'] = b
    st.session_state['test_type_select'] = 'dynamic'
    gui._mark_procedure_applied()
    checks_by_sec = {gui._CHK_SEC_EXP: ['Frequencies field is blank.']}
    assert gui._protocol_confirmed_for_checklist(b, checks_by_sec)


def test_refresh_confirmation_flags_clear_when_dirty():
    _clear_streamlit_session_state()
    b = _DummyBender('dummy_cfg')
    st.session_state['bender'] = b
    st.session_state['gui_setup_confirmed'] = True
    st.session_state['gui_measurements_confirmed'] = True
    st.session_state['gui_protocol_confirmed'] = True
    st.session_state['gui_data_folder'] = r'C:\tmp'
    st.session_state['gui_data_filename'] = 'x.h5'
    gui._mark_data_path_applied()
    st.session_state['gui_morpho_applied_sig'] = gui._morpho_fingerprint()
    st.session_state['gui_proc_applied_sig'] = gui._procedure_fingerprint()

    st.session_state['gui_data_filename'] = 'changed.h5'
    st.session_state['gui_morpho_apply_invalidated'] = True
    st.session_state['gui_proc_apply_invalidated'] = True
    gui._refresh_confirmation_flags()

    assert not st.session_state['gui_setup_confirmed']
    assert not st.session_state['gui_measurements_confirmed']
    assert not st.session_state['gui_protocol_confirmed']


def test_repair_data_path_from_applied_signature():
    _clear_streamlit_session_state()
    st.session_state['gui_data_folder'] = ''
    st.session_state['gui_data_filename'] = ''
    st.session_state['gui_data_path_applied_sig'] = (
        '',
        r'C:\Users\jimen\Desktop\BenderCode\Bender_PC2026\TestData\demo.h5',
    )
    gui._repair_data_path_fields_from_session()
    assert st.session_state['gui_data_folder'].lower().endswith('testdata')
    assert st.session_state['gui_data_filename'] == 'demo.h5'


def test_fix_lines_and_review_usage_helpers():
    _clear_streamlit_session_state()
    lines = gui._build_checklist_fix_lines(
        b=None,
        setup_ok=False,
        measurements_ok=False,
        protocol_ok=False,
    )
    assert any('Setup:' in row for row in lines)
    assert any('Measurements:' in row for row in lines)
    assert any('Protocol/Run:' in row for row in lines)
    assert any('Hardware:' in row for row in lines)

    gui._mark_review_data_used()
    assert st.session_state.get('gui_review_data_used') is True


def test_stepwise_does_not_show_redundant_hardware_not_loaded_banner():
    at = AppTest.from_file('bender_streamlit_gui.py')
    at.session_state['gui_app_route'] = 'stepwise'
    at.session_state['gui_stepwise_step'] = 4
    at.run(timeout=120)
    assert not at.exception
    rendered = '\n'.join(getattr(n, 'value', '') for n in list(at.markdown) + list(at.caption))
    assert 'Hardware not loaded.' not in rendered
