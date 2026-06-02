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
