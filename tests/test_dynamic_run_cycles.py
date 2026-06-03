"""Dynamic protocol: organize_cycles at run/preview without a prior preview."""
from pathlib import Path
import sys
from unittest.mock import patch

import numpy as np
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

from bender_functions import Bender


def _minimal_dynamic_bender():
    b = Bender('jimenez_bender_config_A')
    b.dclamp = 10.0
    b.xsec_width = 2.0
    b.all_freqs = [1.0, 2.0]
    b.all_amps = [0.05]
    b.all_amps_mode = 'strain'
    b.cycles_per_step = 2
    b.n_end_cycles = 1
    b.randomize = False
    b.stim_cycles_in_step = []
    b.stim_pulse_rate = 75.0
    b.is_stim = False
    return b


def test_organize_cycles_for_dynamic_run_builds_period_by_cycle():
    b = _minimal_dynamic_bender()
    assert not hasattr(b, 'period_by_cycle')
    b._organize_cycles_for_dynamic_run()
    assert hasattr(b, 'period_by_cycle')
    pbc = np.asarray(b.period_by_cycle, dtype=float)
    assert pbc.size > 0
    assert np.all(np.isfinite(pbc))
    assert np.all(pbc > 0)


def test_organize_cycles_for_dynamic_run_rebuilds_after_freq_change():
    b = _minimal_dynamic_bender()
    b._organize_cycles_for_dynamic_run()
    dur1 = float(np.sum(b.period_by_cycle))
    b.all_freqs = [0.5]
    b._organize_cycles_for_dynamic_run()
    dur2 = float(np.sum(b.period_by_cycle))
    assert dur2 > dur1


def test_organize_cycles_for_dynamic_run_uses_existing_curves_without_mode():
    """Avoid re-interpreting all_amps as strain when GUI already stored all_curves."""
    b = _minimal_dynamic_bender()
    b.all_amps = [10.0]
    b.all_amps_mode = None
    b.all_curves = [0.1]
    b._organize_cycles_for_dynamic_run()
    assert np.isclose(float(np.max(np.asarray(b.amp_by_cycle))), np.rad2deg(0.1 * (10.0 / 1000.0)))


def test_organize_cycles_respects_muscle_depth_for_strain_metadata():
    b = _minimal_dynamic_bender()
    b.muscle_depth_mm = 0.5
    b.all_curves = [0.2]
    b.all_amps_mode = None
    b._organize_cycles_for_dynamic_run()
    lever = (2.0 / 2.0 - 0.5) / 1000.0
    expected = lever * 0.2
    assert float(np.asarray(b.organized_strains, dtype=float).reshape(-1)[0]) == pytest.approx(expected)


def test_organize_cycles_for_dynamic_run_requires_dclamp():
    b = _minimal_dynamic_bender()
    b.dclamp = None
    with pytest.raises(ValueError, match='dclamp'):
        b._organize_cycles_for_dynamic_run()


def test_run_experiment_dynamic_invokes_organize_at_start():
    b = _minimal_dynamic_bender()
    with patch.object(
        b,
        '_organize_cycles_for_dynamic_run',
        side_effect=RuntimeError('stop-after-organize'),
    ) as oc:
        with pytest.raises(RuntimeError, match='stop-after-organize'):
            b.run_experiment(test_type='dynamic')
    oc.assert_called_once()
