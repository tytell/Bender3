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
    b.target_muscle_depth_mm = 0.5
    b.all_curves = [0.2]
    b.all_amps_mode = None
    b._organize_cycles_for_dynamic_run()
    lever = (2.0 / 2.0 - 0.5) / 1000.0
    expected = lever * 0.2
    assert float(np.asarray(b.organized_strains, dtype=float).reshape(-1)[0]) == pytest.approx(expected)


def test_cycle_tags_follow_generated_motion_timeline():
    b = _minimal_dynamic_bender()
    b.waitbefore = 3.0
    b.waitafter = 2.0
    b.prestim_time = -2.0  # Deliberately differs from waitbefore.
    b._organize_cycles_for_dynamic_run()
    angle, anglevel, tnorm, t = b.make_cycles_dynamic(
        b.period_by_cycle, b.freq_by_cycle, b.amp_by_cycle
    )
    b.t = t
    b.angle = angle
    b.anglevel = anglevel
    b.tnorm = tnorm
    b.aidata = np.zeros((1, t.size), dtype=float)

    b.make_cycle_tags()
    tags = np.asarray(b.cycle_index_history)

    assert np.all(tags[t < 0.0] == -1)
    assert float(t[np.flatnonzero(tags >= 0)[0]]) == pytest.approx(0.0)
    cycle_edges = np.concatenate(([0.0], np.cumsum(np.asarray(b.period_by_cycle))))
    for i, (start_s, end_s) in enumerate(zip(cycle_edges[:-1], cycle_edges[1:])):
        midpoint = (start_s + end_s) / 2.0
        sample = int(np.argmin(np.abs(t - midpoint)))
        assert tags[sample] == i
    assert np.all(tags[t >= cycle_edges[-1]] == -1)


def test_organize_cycles_for_dynamic_run_requires_dclamp():
    b = _minimal_dynamic_bender()
    b.dclamp = None
    with pytest.raises(ValueError, match='dclamp'):
        b._organize_cycles_for_dynamic_run()


def test_run_experiment_dynamic_invokes_organize_at_start():
    b = _minimal_dynamic_bender()
    with patch.object(b, '_announce_pre_run_max_rotation'), patch.object(
        b,
        '_organize_cycles_for_dynamic_run',
        side_effect=RuntimeError('stop-after-organize'),
    ) as oc:
        with pytest.raises(RuntimeError, match='stop-after-organize'):
            b.run_experiment(test_type='dynamic')
    oc.assert_called_once()


def _minimal_sweep_bender():
    b = Bender('jimenez_bender_config_A')
    b.dclamp = 10.0
    b.xsec_width = 2.0
    b.all_freqs = [1.0, 4.0]
    b.all_amps = [0.05]
    b.all_amps_mode = 'strain'
    b.duration = 2.0
    b.waitbefore = 0.1
    b.waitafter = 0.1
    b.amplitude_frequency_exponent = 0.0
    b.is_stim = False
    b.stim_pulse_rate = 75.0
    b.session_simulated = True
    return b


def test_ensure_cycle_stim_arrays_for_sweep_sets_empty_arrays():
    # Issue 13: a fresh sweep Bender has no per-cycle stim arrays (organize_cycles never ran).
    b = _minimal_sweep_bender()
    for name in ('freq_by_cycle', 'phase_by_cycle', 'duty_by_cycle', 'stimburstdur', 'period_by_cycle'):
        assert not hasattr(b, name)
    b._ensure_cycle_stim_arrays_for_sweep()
    for name in ('freq_by_cycle', 'phase_by_cycle', 'duty_by_cycle', 'stimburstdur', 'period_by_cycle'):
        arr = np.asarray(getattr(b, name), dtype=float)
        assert arr.size == 0


def test_frequency_sweep_run_does_not_raise_on_missing_cycle_arrays():
    # Regression: building stimulation_params for the sweep used to raise AttributeError on
    # self.stimburstdur / duty_by_cycle / freq_by_cycle / phase_by_cycle when organize_cycles
    # had never populated them. The sweep run path must now populate them itself.
    b = _minimal_sweep_bender()
    b.run_experiment(test_type='frequency_sweep')
    for name in ('freq_by_cycle', 'phase_by_cycle', 'duty_by_cycle', 'stimburstdur', 'period_by_cycle'):
        assert hasattr(b, name)
    assert np.asarray(b.freq_by_cycle, dtype=float).size == 0


def test_frequency_sweep_run_clears_stale_dynamic_cycle_arrays():
    # A dynamic run first leaves per-cycle arrays populated; a subsequent sweep must clear them
    # so the cycle-indexed stim loop (meaningless for a continuous chirp) does not fire.
    b = _minimal_sweep_bender()
    b.cycles_per_step = 2
    b.n_end_cycles = 0
    b.randomize = False
    b.stim_cycles_in_step = []
    b._organize_cycles_for_dynamic_run()
    assert np.asarray(b.freq_by_cycle, dtype=float).size > 0
    b.run_experiment(test_type='frequency_sweep')
    assert np.asarray(b.freq_by_cycle, dtype=float).size == 0
