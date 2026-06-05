"""Hardware config path resolution and selection matching (no DAQ)."""
from __future__ import annotations

import os
import sys
from pathlib import Path

import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import bender_streamlit_gui as gui


class _DummyBender:
    def __init__(self, config_name: str = 'lab_cfg'):
        self.config_name = config_name


def _clear_session() -> None:
    for key in list(st.session_state.keys()):
        del st.session_state[key]


def test_resolve_hardware_config_import_target_from_absolute_py(tmp_path):
    cfg = tmp_path / 'lab_cfg.py'
    cfg.write_text('device_name = "Dev1"\n', encoding='utf-8')
    stem, abs_path = gui._resolve_hardware_config_import_target(str(cfg))
    assert stem == 'lab_cfg'
    assert os.path.normpath(abs_path) == os.path.normpath(str(cfg))


def test_resolve_hardware_config_import_target_from_module_stem(monkeypatch, tmp_path):
    monkeypatch.setattr(gui, '_ROOT', str(tmp_path))
    cfg = tmp_path / 'jimenez_bender_config_A.py'
    cfg.write_text('device_name = "Dev1"\n', encoding='utf-8')
    stem, abs_path = gui._resolve_hardware_config_import_target('jimenez_bender_config_A')
    assert stem == 'jimenez_bender_config_A'
    assert abs_path == os.path.normpath(str(cfg))


def test_selected_config_matches_bender_same_stem_different_paths(tmp_path, monkeypatch):
    monkeypatch.setattr(gui, '_ROOT', str(tmp_path))
    a = tmp_path / 'a' / 'lab_cfg.py'
    b = tmp_path / 'b' / 'lab_cfg.py'
    a.parent.mkdir()
    b.parent.mkdir()
    a.write_text('x = 1\n', encoding='utf-8')
    b.write_text('x = 2\n', encoding='utf-8')

    _clear_session()
    bender = _DummyBender('lab_cfg')
    st.session_state['gui_loaded_cfg_abs_path'] = str(a)
    st.session_state['gui_load_cfg_file_path'] = str(b)

    assert not gui._selected_config_matches_bender(bender, 'lab_cfg')
    st.session_state['gui_load_cfg_file_path'] = str(a)
    assert gui._selected_config_matches_bender(bender, 'lab_cfg')


def test_raw_mod_for_hardware_config_load_prefers_file_path(tmp_path, monkeypatch):
    monkeypatch.setattr(gui, '_ROOT', str(tmp_path))
    cfg = tmp_path / 'external_cfg.py'
    cfg.write_text('device_name = "Dev1"\n', encoding='utf-8')

    _clear_session()
    st.session_state['gui_load_cfg_file_path'] = str(cfg)

    assert gui._raw_mod_for_hardware_config_load(module_stem='external_cfg') == os.path.normpath(str(cfg))


def test_edit_loaded_config_and_save_new_persists_change(tmp_path, monkeypatch):
    """Unified config editor flow: load a config, edit one field, save to a NEW module, reload it,
    and confirm the edited value persisted (and inherited fields are unchanged)."""
    import importlib

    import bender_config_builder as cb

    # A standalone base config the editor "loads" and seeds from.
    base = tmp_path / 'base_cfg.py'
    base.write_text(
        'device_name = "Dev1"\n'
        'waitbefore = 3.0\n'
        'waitafter = 4.0\n'
        'sono_distance = ""\n'
        # Channel lists referenced by the generated config tail (input_channels rebuild).
        'SG_chan = ["ai0", "ai1"]\n'
        'SG_name = ["xForce", "yForce"]\n'
        'use_sono = False\n'
        'sono_channel = []\n'
        'sono_name = []\n'
        'stim_monitor_chan = []\n'
        'stim_monitor_name = []\n',
        encoding='utf-8',
    )
    monkeypatch.syspath_prepend(str(tmp_path))
    importlib.invalidate_caches()

    # Seed (load) defaults from the base, then edit one field as the GUI would.
    d = cb.read_base_defaults('base_cfg')
    assert d['waitbefore'] == 3.0
    edited = {'waitbefore': 9.5, 'sono_distance': '12.5,14.2', 'waitafter': d['waitafter']}

    # Save to a NEW module that inherits from the base via ``import *`` and overrides edits.
    src = cb.render_generated_config('base_cfg', edited)
    new_mod = tmp_path / 'edited_cfg.py'
    new_mod.write_text(src, encoding='utf-8')
    importlib.invalidate_caches()

    reloaded = importlib.import_module('edited_cfg')
    assert reloaded.waitbefore == 9.5  # edited value persisted
    assert reloaded.sono_distance == '12.5,14.2'  # edited value persisted
    assert reloaded.waitafter == 4.0  # inherited from base, unchanged
    assert reloaded.device_name == 'Dev1'  # inherited from base via import *
