import warnings
from pathlib import Path
import sys

from streamlit.testing.v1 import AppTest

sys.path.append(str(Path(__file__).resolve().parents[1]))

import bender_simulation as sim


def test_oscillatory_handles_malformed_numeric_inputs_without_crash():
    cases = [
        dict(
            length_mm=30,
            width_mm=25.4,
            material_key='polyurethane',
            freq_hz=1.0,
            amplitude_mm=2.0,
            n_cycles_shown=3.0,
            points_per_cycle=200,
            noise_std_N=0.01,
        ),
        dict(
            length_mm=0,
            width_mm=0,
            material_key='unknown',
            freq_hz=0,
            amplitude_mm=0,
            n_cycles_shown=0,
            points_per_cycle=0,
            noise_std_N=-1,
        ),
        dict(
            length_mm=-5,
            width_mm=-3,
            material_key=None,
            freq_hz=-1,
            amplitude_mm=-2,
            n_cycles_shown=-3,
            points_per_cycle=-10,
            noise_std_N='0.2',
        ),
        dict(
            length_mm=float('nan'),
            width_mm=float('nan'),
            material_key='silicone',
            freq_hz=float('nan'),
            amplitude_mm=float('nan'),
            n_cycles_shown=float('nan'),
            points_per_cycle='abc',
            noise_std_N='bad',
        ),
    ]

    for case in cases:
        with warnings.catch_warnings(record=True) as seen:
            warnings.simplefilter('always')
            out = sim.oscillatory_viscoelastic_timeseries(**case)
        assert out['t'].size >= 80
        assert out['delta_mm'].size == out['F'].size
        runtime_warnings = [w for w in seen if issubclass(w.category, RuntimeWarning)]
        assert not runtime_warnings


def test_static_and_quasistatic_curve_accept_unexpected_types():
    for val in [120, 0, -10, float('nan'), '3', 'x', None]:
        d, f1, f2 = sim.static_stiffness_comparison_delta_grid(val, delta_max_mm='bad', n_points='bad')
        assert d.size >= 2
        assert d.size == f1.size == f2.size

    d, f = sim.force_displacement_comparison_curve(
        length_mm='oops',
        width_mm='oops',
        material_key='silicone',
        delta_max_mm='bad',
        n_points='bad',
        noise_std_N='bad',
    )
    assert d.size >= 2
    assert d.size == f.size


def test_streamlit_apptest_sim_compare_route_runs_with_extreme_but_valid_values():
    at = AppTest.from_file('bender_streamlit_gui.py')
    at.session_state['gui_app_route'] = 'sim_compare'
    at.run(timeout=120)
    assert not at.exception

    at.session_state['gui_sim_cmp_oscillatory'] = True
    at.session_state['gui_sim_cmp_L_a'] = 5.0
    at.session_state['gui_sim_cmp_L_b'] = 150.0
    at.session_state['gui_sim_cmp_w_a'] = 2.0
    at.session_state['gui_sim_cmp_w_b'] = 50.0
    at.session_state['gui_sim_cmp_freq_hz'] = 8.0
    at.session_state['gui_sim_cmp_amp_mm'] = 12.0
    at.session_state['gui_sim_cmp_cycles'] = 10.0
    at.session_state['gui_sim_cmp_ppc'] = 400
    at.session_state['gui_sim_cmp_noise'] = 0.2
    at.run(timeout=120)
    assert not at.exception
