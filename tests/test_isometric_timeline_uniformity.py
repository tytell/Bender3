"""Regression: isometric ramp-hold timelines must stay uniformly sampled.

A zero-angle-change step with a finite ramp (the first step of a non-block isometric
run, or any step that repeats the previous target) used to collapse the ramp to a
single sample while still offsetting the hold by the full ramp duration. That opened a
``ramp_s``-wide gap as ``t[1]-t[0]``, which:

  * made ``record_motor_signal`` infer a ~1/ramp_s AI rate (snapped back to the config
    rate by the setter guard), and
  * under-sized the FINITE AI buffer (``len(t)``) relative to the AO/DO buffer (built
    across the full ``t[0]..t[-1]`` span),

so the motor/stim generation ran longer than the acquisition and was stopped before
completing -> NI warning 200010 on every isometric trial. These tests pin the timeline
uniformity and the AI-vs-AO/DO duration invariant.
"""
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

from bender_functions import Bender


def _bender():
    b = Bender('jimenez_bender_config_A')  # ai = 1000 Hz, ao/do = 60000 Hz
    return b


def test_ramp_hold_zero_angle_change_is_uniform():
    """ramp_s > 0 with no angle change must sample the ramp span flat, not collapse it."""
    b = _bender()
    ramp_s, hold_s, daq_hz = 2.0, 7.0, 1000.0
    t, angle, anglevel = b._timeline_ramp_hold(10.0, 10.0, ramp_s, hold_s, daq_hz)

    dt = t[1] - t[0]
    # No gap: the first spacing is the AI period, not the whole ramp.
    assert dt == pytest.approx(1.0 / daq_hz, rel=1e-9)
    # Uniform spacing across the entire ramp+hold timeline (the gap bug broke this).
    assert np.allclose(np.diff(t), dt, rtol=0, atol=1e-9)
    # Span matches the sample count: t[-1]-t[0] == (len-1)*dt.
    assert (t[-1] - t[0]) == pytest.approx((t.size - 1) * dt, abs=1e-9)
    # Same length as a real ramp of the same durations (no missing ramp samples).
    t_real, _, _ = b._timeline_ramp_hold(0.0, 10.0, ramp_s, hold_s, daq_hz)
    assert t.size == t_real.size
    # Commanded angle is flat and velocity is ~0 over the whole (no-op) move.
    assert np.allclose(angle, 10.0)
    assert np.allclose(anglevel, 0.0, atol=1e-9)


def test_ramp_hold_zero_ramp_still_collapses():
    """ramp_s <= 0 (no time for a ramp) still collapses to a single ramp sample."""
    b = _bender()
    t, angle, _ = b._timeline_ramp_hold(0.0, 5.0, 0.0, 3.0, 1000.0)
    # First sample is the (collapsed) ramp point at t=0; the hold follows immediately at dt.
    assert t[0] == 0.0
    assert (t[1] - t[0]) == pytest.approx(0.001, rel=1e-9)
    assert angle[0] == 5.0


def test_zero_angle_step_ai_not_shorter_than_aodo():
    """End to end: AI acquisition must not finish before AO/DO generation (avoids NI 200010)."""
    b = _bender()
    t, angle, anglevel = b._timeline_ramp_hold(8.0, 8.0, 2.0, 7.0, 1000.0)
    b.record_motor_signal(t, angle, anglevel, tnorm=np.zeros_like(t))
    b.record_stim_signal(np.zeros_like(t), np.zeros_like(t))
    # AI rate is the clean config rate, not the corrupted ~0.5 Hz the gap produced.
    assert b.daq_ai_sample_rate_hz == 1000.0
    b.make_motor_stepper_pulses(
        daq_ao_do_sample_rate_hz=b.daq_ao_do_sample_rate_hz,
        motor_gear_ratio=b.motor_gear_ratio,
        motor_full_steps_per_rev=b.motor_full_steps_per_rev,
    )
    ai_duration_s = len(b.t) / float(b.daq_ai_sample_rate_hz)
    aodo_duration_s = len(b.tout) / float(b.daq_ao_do_sample_rate_hz)
    # The AO/DO finite generation is triggered by AI start; it must finish first.
    assert ai_duration_s >= aodo_duration_s


def test_record_motor_signal_rejects_gapped_timeline():
    """A gapped timeline must fail fast before it can mis-size the DAQ buffers."""
    b = _bender()
    # Normal dt, but a 2 s gap in the middle (isolates the span check from rate inference).
    t = np.concatenate([np.arange(0.0, 2.0, 0.001), 4.0 + np.arange(0.0, 5.0, 0.001)])
    ang = np.zeros_like(t)
    with pytest.raises(ValueError, match='not uniformly sampled'):
        b.record_motor_signal(t, ang, ang)


def test_record_motor_signal_accepts_uniform_timeline():
    """A uniform timeline (the fixed builder output) passes the guard unchanged."""
    b = _bender()
    t, angle, anglevel = b._timeline_ramp_hold(0.0, 0.0, 2.0, 5.0, 1000.0)
    b.record_motor_signal(t, angle, anglevel)
    assert b.daq_ai_sample_rate_hz == 1000.0
    assert len(b.t) == t.size
