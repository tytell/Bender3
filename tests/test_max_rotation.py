"""Pre-run maximum commanded motor rotation (all protocol types)."""
from pathlib import Path
import sys
import numpy as np
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

from bender_functions import Bender
from shared.utilities import compute_max_rotation_deg


def _minimal_dynamic_bender():
    b = Bender('jimenez_bender_config_A')
    b.dclamp = 10.0
    b.xsec_width = 2.0
    b.all_freqs = [1.0]
    b.all_amps = [0.05]
    b.all_amps_mode = 'strain'
    b.cycles_per_step = 2
    b.n_end_cycles = 0
    b.randomize = False
    b.waitbefore = 0.0
    b.waitafter = 0.0
    b.rampdur = 0.01
    b.amp_step_vel = 5.0
    return b


@pytest.mark.parametrize(
    'test_type,setup',
    [
        ('dynamic', _minimal_dynamic_bender),
        (
            'isometric',
            lambda: (
                (b := _minimal_dynamic_bender()),
                setattr(b, 'test_type', 'isometric'),
                setattr(b, 'isometric_initial', 0.02),
                setattr(b, 'isometric_final', 0.06),
                setattr(b, 'isometric_num_steps', 3),
                setattr(b, 'isometric_mode', 'strain'),
                b,
            )[-1],
        ),
        (
            'isovelocity',
            lambda: (
                (b := _minimal_dynamic_bender()),
                setattr(b, 'test_type', 'isovelocity'),
                setattr(b, 'isovelocity_min_vel', 1.0),
                setattr(b, 'isovelocity_max_vel', 3.0),
                setattr(b, 'isovelocity_starting_strain', 0.03),
                setattr(b, 'isovelocity_starting_strain_mode', 'strain'),
                setattr(b, 'isovelocity_velocity_mode', 'angle_vel'),
                setattr(b, 'isovelocity_num_steps', 2),
                setattr(b, 'isovelocity_iso_duration_s', 0.2),
                setattr(b, 'isovelocity_pre_hold_s', 0.1),
                b,
            )[-1],
        ),
        (
            'frequency_sweep',
            lambda: (
                (b := _minimal_dynamic_bender()),
                setattr(b, 'test_type', 'frequency_sweep'),
                setattr(b, 'all_freqs', [0.5, 2.0]),
                setattr(b, 'duration', 1.0),
                setattr(b, 'amplitude_frequency_exponent', 0.0),
                b,
            )[-1],
        ),
        (
            'step_change',
            lambda: (
                (b := _minimal_dynamic_bender()),
                setattr(b, 'test_type', 'step_change'),
                setattr(b, 'step_change_frequencies', [1.0]),
                setattr(b, 'step_change_curves', [0.1]),
                setattr(b, 'step_change_cycles_per_step', [2]),
                b,
            )[-1],
        ),
    ],
)
def test_compute_max_rotation_deg_positive_finite(test_type, setup):
    b = setup()
    rep = b.validate_dispatch_setup(test_type=test_type)
    assert rep['ok'], rep.get('missing')
    deg = compute_max_rotation_deg({'bender': b, 'test_type': test_type})
    assert np.isfinite(deg) and deg > 0
    assert rep['max_rotation_deg'] == pytest.approx(deg)


def test_isometric_single_step_fallback():
    b = _minimal_dynamic_bender()
    b.test_type = 'isometric'
    b.isometric_initial = 0.02
    b.isometric_final = 0.04
    b.isometric_num_steps = 1
    b.isometric_mode = 'strain'
    deg = compute_max_rotation_deg({'bender': b, 'test_type': 'isometric'})
    assert deg > 0


def test_announce_pre_run_max_rotation_prints_line(capsys):
    b = _minimal_dynamic_bender()
    b._announce_pre_run_max_rotation('dynamic')
    out = capsys.readouterr().out
    assert 'maximum of' in out
    assert getattr(b, 'max_commanded_rotation_deg', None) is not None
