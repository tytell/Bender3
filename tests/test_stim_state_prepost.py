"""Regression: dynamic pre/post (-1) stim tagging vs numbered-cycle passive/off/on."""
from pathlib import Path
import sys

import numpy as np

sys.path.append(str(Path(__file__).resolve().parents[1]))

from bender_h5_export import _compute_stim_phase_and_activity  # noqa: E402


def test_numbered_cycles_keep_passive_off_on():
    # Samples: quiet cycle0 | pulse cycle0 | quiet cycle1 (no stim) | pulse cycle2 | off cycle2
    cyc = np.array([0, 0, 1, 2, 2], dtype=float)
    on = np.array([False, True, False, True, False])
    phase, activity = _compute_stim_phase_and_activity(cyc, on, stim_enabled=True)

    assert list(phase) == ['off', 'on', 'passive', 'on', 'off']
    assert list(activity) == ['active', 'active', 'passive', 'active', 'active']


def test_prepost_minus1_is_off_on_not_passive():
    # Shared -1 for pre and post on a stim trial: quiet=off, pulse=on (historical categories).
    cyc = np.array([-1, -1, -1, -1, -1], dtype=float)
    on = np.array([False, True, False, True, False])
    phase, activity = _compute_stim_phase_and_activity(cyc, on, stim_enabled=True)

    assert list(phase) == ['off', 'on', 'off', 'on', 'off']
    assert list(activity) == ['active', 'active', 'active', 'active', 'active']
    assert 'passive' not in phase


def test_mixed_timeline_cycles_and_minus1_both_off_on():
    # pre(-1) quiet, pre pulse, cycle0 off, cycle0 on, post(-1) quiet, post pulse
    cyc = np.array([-1, -1, 0, 0, -1, -1], dtype=float)
    on = np.array([False, True, False, True, False, True])
    phase, activity = _compute_stim_phase_and_activity(cyc, on, stim_enabled=True)

    assert list(phase) == ['off', 'on', 'off', 'on', 'off', 'on']
    assert list(activity) == ['active', 'active', 'active', 'active', 'active', 'active']


def test_stim_disabled_stays_passive():
    cyc = np.array([-1, 0, 0], dtype=float)
    on = np.array([True, True, False])
    phase, activity = _compute_stim_phase_and_activity(cyc, on, stim_enabled=False)
    assert list(phase) == ['passive', 'passive', 'passive']
    assert list(activity) == ['passive', 'passive', 'passive']
