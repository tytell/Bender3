"""Regression: daq_ai_sample_rate_hz must never hold a corrupt/implausible value.

A stale 0.5 Hz (equal to a typical drive frequency) under-samples dynamic sine
commands into a ramp, so the motor steps one direction and "walks". The validated
property snaps any invalid/implausibly low rate back to the hardware config rate,
breaking the self-perpetuating loop where record_motor_signal re-saved 1/dt each run.
"""
from pathlib import Path
import sys

import numpy as np

sys.path.append(str(Path(__file__).resolve().parents[1]))

from bender_functions import Bender


def _dynamic_bender():
    b = Bender('jimenez_bender_config_A')
    b.dclamp = 10.0
    b.xsec_width = 2.0
    b.all_freqs = [0.5]
    b.all_amps = [5.0]
    b.all_amps_mode = 'angle'
    b.cycles_per_step = 5
    b.n_end_cycles = 0
    b.randomize = False
    b.stim_cycles_in_step = []
    b.stim_pulse_rate = 75.0
    b.is_stim = False
    return b


def test_config_rate_loaded():
    b = Bender('jimenez_bender_config_A')
    assert b.daq_ai_sample_rate_hz == 1000.0


def test_corrupt_low_rate_is_rejected():
    b = Bender('jimenez_bender_config_A')
    b.daq_ai_sample_rate_hz = 0.5  # the corruption we observed
    assert b.daq_ai_sample_rate_hz == 1000.0


def test_various_invalid_rates_reset_to_config():
    b = Bender('jimenez_bender_config_A')
    for bad in (0.0, -10.0, 0.5, 49.9, float('nan'), float('inf'), None, 'oops'):
        b.daq_ai_sample_rate_hz = bad
        assert b.daq_ai_sample_rate_hz == 1000.0, f'rate {bad!r} should reset to config'


def test_valid_rate_is_accepted():
    b = Bender('jimenez_bender_config_A')
    b.daq_ai_sample_rate_hz = 2000.0
    assert b.daq_ai_sample_rate_hz == 2000.0


def test_record_motor_signal_cannot_drive_rate_to_half_hz():
    """A 2 s-spaced timeline must not persist a 0.5 Hz rate (the old failure loop)."""
    b = Bender('jimenez_bender_config_A')
    t = np.arange(0.0, 10.0, 2.0)  # dt = 2 s -> 1/dt = 0.5 Hz
    ang = np.zeros_like(t)
    av = np.zeros_like(t)
    b.record_motor_signal(t, ang, av)
    assert b.daq_ai_sample_rate_hz == 1000.0


def test_dynamic_command_is_a_sine_not_a_ramp_after_attempted_corruption():
    """Even after a corruption attempt, the dynamic command must be a well-sampled,
    reversing sine (the symptom was a single monotonic ramp with ~9 samples)."""
    b = _dynamic_bender()
    b.daq_ai_sample_rate_hz = 0.5  # attempt the corruption; guard resets to 1000
    b._organize_cycles_for_dynamic_run()
    angle, anglevel, tnorm, t = b.make_cycles_dynamic(
        b.period_by_cycle, b.freq_by_cycle, b.amp_by_cycle, record_protocol=False
    )
    # Many samples (not the ~9 the 0.5 Hz bug produced).
    assert t.size > 1000
    # The command actually oscillates: velocity changes sign (reverses) many times.
    sign_changes = int(np.sum(np.diff(np.sign(anglevel[anglevel != 0])) != 0))
    assert sign_changes >= 5, 'dynamic command should reverse direction, not ramp one way'
