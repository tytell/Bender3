"""Fix 1: motor stays energized/braked across segment boundaries (hold_motor_between_steps).

Two independent behaviors:
- Motor side: with the flag ON (default), the ENABLE line stays high through the last DO
  sample and run() reasserts ENABLE-high in the inter-segment gap via a short-lived DO task.
  Regardless of the flag, run() also PRE-energizes a de-energized driver (idle word + dwell)
  before starting a waveform, so STEP pulses are never sent during the drive's enable sequence.
- DAQ side: acquisition pause/flush/device-reset is unchanged regardless of the flag.

The make_motor_stepper_pulses checks are hardware-free (pure signal math, runs on Mac); the
run() reassert is exercised with NI symbols mocked (no hardware), mirroring
test_run_twice_device_reset.py.
"""
from pathlib import Path
import sys
from unittest.mock import patch, MagicMock

import numpy as np
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

import bender_daq_kill
import bender_functions
from bender_functions import Bender

AO_HZ = 60000.0
GEAR = 5
STEPS_PER_REV = 1600
DAQ_HZ = 1000.0

# Enable bit is P0.2 in the packbits column order [0,0,0,0,0, ENABLE, STEP, DIR] (MSB first),
# so it lands at bit 2 of each uint32 DO word.
ENABLE_BIT = 2


def _bender():
    return Bender('jimenez_bender_config_A')


def _last_enable_bit(b):
    """Build one ramp+hold segment and return the ENABLE bit of the final DO sample."""
    t, a, v = b._timeline_ramp_hold(0.0, 10.0, 0.5, 0.5, DAQ_HZ)
    b.record_motor_signal(t, a, v)
    b.record_stim_signal(None, None)
    b.make_motor_stepper_pulses(
        daq_ao_do_sample_rate_hz=AO_HZ,
        motor_gear_ratio=GEAR,
        motor_full_steps_per_rev=STEPS_PER_REV,
    )
    return int((int(b.dig[-1]) >> ENABLE_BIT) & 1)


def test_flag_defaults_true():
    assert _bender().hold_motor_between_steps is True


def test_hold_on_keeps_enable_high_through_last_sample():
    b = _bender()
    assert b.hold_motor_between_steps is True
    assert _last_enable_bit(b) == 1


def test_hold_off_drops_enable_on_last_sample_legacy():
    b = _bender()
    b.hold_motor_between_steps = False
    assert _last_enable_bit(b) == 0


def test_pack_motor_do_word_matches_column_order():
    b = _bender()
    # [0,0,0,0,0, ENABLE, STEP, DIR] MSB-first -> ENABLE=bit2(4), STEP=bit1(2), DIR=bit0(1).
    assert b._pack_motor_do_word(enable=0, step=0, direction=0) == 0
    assert b._pack_motor_do_word(enable=1, step=0, direction=0) == 4
    assert b._pack_motor_do_word(enable=1, step=1, direction=1) == 7
    assert b._pack_motor_do_word(enable=0, step=1, direction=0) == 2


def test_idle_word_enable_bit_consistent_with_waveform_packing():
    """The idle word's ENABLE bit decodes the same way as a real DO sample's ENABLE bit."""
    b = _bender()
    idle = b._pack_motor_do_word(enable=1, step=0, direction=0)
    assert (idle >> ENABLE_BIT) & 1 == 1
    # And it carries no STEP/DIR set.
    assert idle == 4


def test_idle_word_reflects_last_waveform_direction():
    """Task B: the held idle word reuses the waveform's final DIR bit so an energized hold does
    not toggle DIR across the gap. A reverse-ending waveform -> idle word DIR=reverse (word 5)."""
    b = _bender()
    t, a, v = b._timeline_ramp_hold(0.0, -10.0, 0.5, 0.5, DAQ_HZ)
    b.record_motor_signal(t, a, v)
    b.record_stim_signal(None, None)
    b.make_motor_stepper_pulses(
        daq_ao_do_sample_rate_hz=AO_HZ, motor_gear_ratio=GEAR, motor_full_steps_per_rev=STEPS_PER_REV,
    )
    assert b._last_motor_direction_bit == 1
    idle = b._pack_motor_do_word(enable=1, step=0, direction=b._last_motor_direction_bit)
    assert idle == 5                                  # ENABLE(4) | DIR(1), STEP clear


# --- run() reassert across the gap (NI mocked) ----------------------------------------

class _CtxTask:
    """nidaqmx.Task stand-in supporting the ``with`` protocol used by Bender.run()."""

    def __init__(self, *_a, **_k):
        self.ai_channels = MagicMock()
        self.ci_channels = MagicMock()
        self.ao_channels = MagicMock()
        self.do_channels = MagicMock()
        self.timing = MagicMock()
        self.triggers = MagicMock()
        self.in_stream = MagicMock()
        self.out_stream = MagicMock()
        self.stop = MagicMock()
        self.start = MagicMock()
        self.wait_until_done = MagicMock(return_value=None)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _minimal_dynamic_bender(hold):
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
    b.simulation_mode = False
    b.hold_motor_between_steps = hold
    # Keep the mocked acquisition fast: the pre-energize dwell is real wall-clock sleep.
    b.motor_enable_dwell_s = 0.01
    return b


def _run_dynamic_capturing_do_writes(hold):
    """Run one mocked dynamic acquisition; return the single-sample DO writes (idle words)."""
    dig_writer_instance = MagicMock()
    dig_writer_instance.write_many_sample_port_uint32.side_effect = (
        lambda data, *a, **k: int(np.asarray(data).reshape(-1).size)
    )
    one_sample_calls: list = []
    dig_writer_instance.write_one_sample_port_uint32.side_effect = (
        lambda word, *a, **k: one_sample_calls.append(int(word))
    )

    def _fake_emergency_stop(device_name=None, release_motor_enable_line=None):
        return True, 'mock reset'

    patches = [
        patch.object(bender_functions, 'Task', _CtxTask),
        patch.object(bender_functions, 'AnalogMultiChannelReader', MagicMock()),
        patch.object(bender_functions, 'CounterReader', MagicMock()),
        patch.object(bender_functions, 'AnalogMultiChannelWriter', MagicMock()),
        patch.object(
            bender_functions, 'DigitalSingleChannelWriter',
            MagicMock(return_value=dig_writer_instance),
        ),
        patch.object(bender_functions, 'daq', MagicMock()),
        patch.object(bender_functions, 'TerminalConfiguration', MagicMock()),
        patch.object(bender_daq_kill, 'daq_emergency_stop', _fake_emergency_stop),
    ]
    b = _minimal_dynamic_bender(hold)
    for p in patches:
        p.start()
    try:
        b.run_experiment(test_type='dynamic')
    finally:
        for p in patches:
            p.stop()
    return one_sample_calls


# STEP is P0.1 -> bit 1; DIR is P0.0 -> bit 0 in the same column order.
STEP_BIT = 1
DIR_BIT = 0


def test_run_reasserts_enable_high_in_gap_when_holding():
    """With the flag ON, run() reasserts an enable-high, step-low idle word after the device
    reset. The DIR bit now mirrors the waveform's final direction (Task B: no spurious DIR flip
    across the gap), so it may be 0 or 1 -- assert ENABLE high and STEP low, not the literal 4."""
    calls = _run_dynamic_capturing_do_writes(hold=True)
    assert calls, 'expected at least one enable-high idle write when holding the motor'
    for w in calls:
        assert (w >> ENABLE_BIT) & 1 == 1, f'idle word must keep ENABLE high, saw {w}'
        assert (w >> STEP_BIT) & 1 == 0, f'idle word must keep STEP low (no pulse in gap), saw {w}'
        assert (w >> DIR_BIT) & 1 in (0, 1), f'idle word DIR bit malformed, saw {w}'


def test_run_legacy_mode_writes_are_pre_energize_only():
    """With the flag OFF, run() still leaves the motor de-energized BETWEEN segments (the
    waveform's last sample drops ENABLE -- covered above), but it now PRE-energizes the driver
    before each waveform starts (lost-step guard: the drive ignores STEP input while its enable
    sequence runs). So single-sample DO writes exist, every one is the enable-high/step-low idle
    form, and none marks the driver as held across the gap."""
    calls = _run_dynamic_capturing_do_writes(hold=False)
    assert calls, 'expected pre-energize idle writes in legacy mode'
    for w in calls:
        assert (w >> ENABLE_BIT) & 1 == 1, f'pre-energize word must set ENABLE, saw {w}'
        assert (w >> STEP_BIT) & 1 == 0, f'pre-energize word must keep STEP low, saw {w}'
