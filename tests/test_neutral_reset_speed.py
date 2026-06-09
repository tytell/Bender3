"""Unit tests for the speed-capped neutral-reset ramp.

The return-to-neutral ramp is open-loop stepper motion: a fixed duration regardless of amplitude
slews large resets too fast and loses steps (residual scales with commanded speed). The ramp
duration is therefore sized so the constant slew rate of the linear ramp never exceeds
``reset_max_speed_deg_per_s`` (duration = |from_deg| / max_speed, floored at two AI samples).
"""
from __future__ import annotations

import numpy as np
import pytest

from bender_functions import Bender


class _ResetBender(Bender):
    """Minimal instance exposing only what the reset-duration helper reads."""

    def __init__(self, *, daq_hz=1000.0, max_speed=15.0):
        self.daq_ai_sample_rate_hz = float(daq_hz)
        self.reset_max_speed_deg_per_s = float(max_speed)


def test_duration_is_amplitude_over_max_speed():
    b = _ResetBender(max_speed=15.0)
    # 45 deg at 15 deg/s -> 3 s (well above the AI-sample floor).
    assert b._neutral_reset_ramp_duration_s(45.0) == pytest.approx(45.0 / 15.0)
    # Sign does not matter; only the magnitude sets the duration.
    assert b._neutral_reset_ramp_duration_s(-45.0) == pytest.approx(45.0 / 15.0)


def test_duration_scales_linearly_with_amplitude():
    """Constant slew rate: doubling the amplitude doubles the ramp duration."""
    b = _ResetBender(max_speed=15.0)
    d1 = b._neutral_reset_ramp_duration_s(20.0)
    d2 = b._neutral_reset_ramp_duration_s(40.0)
    assert d2 == pytest.approx(2.0 * d1)


def test_peak_commanded_speed_does_not_exceed_cap():
    """The realized ramp's peak |velocity| must not exceed the speed cap."""
    b = _ResetBender(max_speed=15.0)
    from_deg = 82.0
    ramp_s = b._neutral_reset_ramp_duration_s(from_deg)
    _t, _angle, anglevel = b._timeline_ramp_hold(from_deg, 0.0, ramp_s, 0.0, b.daq_ai_sample_rate_hz)
    peak = float(np.max(np.abs(np.asarray(anglevel, dtype=float))))
    # Allow a small numerical margin from sample-grid quantization of the ramp.
    assert peak <= 15.0 * 1.05


def test_lower_max_speed_gives_longer_ramp():
    fast = _ResetBender(max_speed=30.0)
    slow = _ResetBender(max_speed=10.0)
    assert slow._neutral_reset_ramp_duration_s(60.0) > fast._neutral_reset_ramp_duration_s(60.0)


def test_tiny_amplitude_is_floored_to_two_ai_samples():
    """A near-zero reset still spans >= 2 AI samples so the ramp cannot collapse."""
    b = _ResetBender(daq_hz=1000.0, max_speed=15.0)
    min_ramp_s = 2.0 / b.daq_ai_sample_rate_hz
    # 0.001 deg / 15 deg/s = 6.7e-5 s, far below the floor.
    assert b._neutral_reset_ramp_duration_s(0.001) == pytest.approx(min_ramp_s)


@pytest.mark.parametrize('bad', [0.0, -5.0, float('nan'), float('inf')])
def test_nonpositive_or_nonfinite_max_speed_raises(bad):
    b = _ResetBender(max_speed=15.0)
    b.reset_max_speed_deg_per_s = bad
    with pytest.raises(ValueError, match='reset_max_speed_deg_per_s'):
        b._neutral_reset_ramp_duration_s(30.0)


def test_preview_helper_matches_backend_duration():
    """Preview must size the reset ramp exactly like the run (single source of truth)."""
    from bender_gui_preview import _preview_neutral_reset_ramp_s

    b = _ResetBender(daq_hz=1000.0, max_speed=15.0)
    for amp in (5.0, 45.0, 82.0):
        assert _preview_neutral_reset_ramp_s(b, amp, b.daq_ai_sample_rate_hz) == pytest.approx(
            b._neutral_reset_ramp_duration_s(amp)
        )
