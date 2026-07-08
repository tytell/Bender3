"""Regression test: aor-to-clamp distance must survive Apply with no specimen geometry.

Verifies the fix in _apply_clamp_geometry_to_bender that unconditionally writes
specimen_profile_clamp_offset_mm onto the Bender object.  Before the fix, the
value was only forwarded inside set_specimen_geometry_inertial_model -- skipped
on every empty-apparatus calibration run -- so the exporter emitted NaN.

Test must FAIL on pre-fix code and PASS after.
"""
import math
from pathlib import Path
import sys

import h5py
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))
sys.path.append(str(Path(__file__).resolve().parents[1] / 'templates' / 'configs'))

import streamlit as st  # noqa: E402
from bender_functions import Bender  # noqa: E402
from bender_h5_export import export_primary_h5  # noqa: E402
import bender_streamlit_gui as gui  # noqa: E402

_AOR_MM = 20.0  # positive so the exporter's > 0 gate passes


def _clear_session():
    for k in list(st.session_state.keys()):
        st.session_state.pop(k, None)


def _set_section4_keys(aor_mm: float = _AOR_MM):
    """Populate the Section 4 session-state keys that _apply_clamp_geometry_to_bender reads."""
    st.session_state['morpho_dclamp'] = 10.0
    st.session_state['morpho_dbend'] = 0.0
    st.session_state['morpho_xsec'] = 5.0        # must be > 0 (min_value=0.001 in widget)
    st.session_state['morpho_xsec_height'] = 5.0  # must be > 0
    st.session_state['morpho_dvert'] = 0.0
    st.session_state['morpho_dhoriz'] = 0.0
    st.session_state['morpho_muscle_depth'] = 0.0
    st.session_state['morpho_clamp_plate_extension'] = 15.0
    st.session_state['morpho_prof_clamp'] = aor_mm
    # Leave morpho_geom_x/y/pos blank so the geometry model is NOT invoked.
    # This isolates the fix: pre-fix code would have left specimen_profile_clamp_offset_mm = 0.0.
    st.session_state['morpho_geom_x'] = ''
    st.session_state['morpho_geom_y'] = ''
    st.session_state['morpho_geom_pos'] = ''


def _make_bender():
    b = Bender('jimenez_bender_config_A')
    b.session_simulated = True
    b.is_stim = False
    b.stim_pulse_rate = 75.0
    b.isometric_initial = -5.0
    b.isometric_final = 5.0
    b.isometric_num_steps = 2
    b.isometric_mode = 'angle'
    return b


def test_aor_survives_apply_without_specimen_geometry(tmp_path):
    """aor-to-clamp key must be finite and equal to the input on an empty-apparatus run."""
    _clear_session()
    b = _make_bender()
    _set_section4_keys()

    ok = gui._apply_clamp_geometry_to_bender(b)
    assert ok, '_apply_clamp_geometry_to_bender returned False -- check morpho_xsec validity'

    assert b.specimen_profile_clamp_offset_mm == _AOR_MM, (
        f'specimen_profile_clamp_offset_mm not written to Bender: '
        f'expected {_AOR_MM}, got {b.specimen_profile_clamp_offset_mm}'
    )

    out = str(tmp_path / 'aor_cal.h5')
    b.outputfile = out
    b.run_experiment(test_type='isometric')
    res = export_primary_h5(b, outputfile=out)
    actual_path = res['outputfile']

    with h5py.File(actual_path, 'r') as f:
        aor_val = f['metadata'].attrs['calibration_inertia_apparatus_aor_to_clamp_millimeter']

    assert math.isfinite(aor_val), (
        f'calibration_inertia_apparatus_aor_to_clamp_millimeter is NaN -- '
        f'specimen_profile_clamp_offset_mm was not written unconditionally'
    )
    assert aor_val == pytest.approx(_AOR_MM), (
        f'Expected {_AOR_MM}, got {aor_val}'
    )
