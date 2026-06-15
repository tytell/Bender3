"""Regression: the dynamic run path must release/reset the NI device after EVERY run.

Background: commit 595b762 added ``finally: _stop_run_tasks()`` in ``Bender.run()`` but it
only stopped the NI task objects; it never reset the device. On Windows/NI a second
back-to-back FINITE acquisition then wedges in ``wait_until_done()``. The fix resets the
device (``bender_daq_kill.daq_emergency_stop``) in the ``finally`` after every run, so two
consecutive runs both tear down cleanly.

Also guards the ``stimburstdur`` 0/0 divide that warned when stim was disabled.
"""
from pathlib import Path
import sys
import warnings
from unittest.mock import patch, MagicMock

import numpy as np
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

import bender_daq_kill
import bender_functions
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
    b.session_simulated = False
    return b


class _CtxTask:
    """Stand-in for nidaqmx.Task supporting the ``with`` protocol used by Bender.run()."""

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


def _patched_daq_run(reset_counter):
    """Context: replace every NI symbol used by Bender.run() with safe fakes.

    ``reset_counter`` is a list whose length is incremented each time the device-reset
    teardown fires, so callers can assert it ran once per run.
    """

    def _fake_emergency_stop(device_name=None, release_motor_enable_line=None):
        reset_counter.append({'device': device_name, 'release': release_motor_enable_line})
        return True, 'mock reset'

    dig_writer_instance = MagicMock()
    # ``write_many_sample_port_uint32`` must report it wrote every sample (== len(dig)).
    dig_writer_instance.write_many_sample_port_uint32.side_effect = (
        lambda data, *a, **k: int(np.asarray(data).reshape(-1).size)
    )

    return [
        patch.object(bender_functions, 'Task', _CtxTask),
        patch.object(bender_functions, 'AnalogMultiChannelReader', MagicMock()),
        patch.object(bender_functions, 'CounterReader', MagicMock()),
        patch.object(bender_functions, 'AnalogMultiChannelWriter', MagicMock()),
        patch.object(
            bender_functions,
            'DigitalSingleChannelWriter',
            MagicMock(return_value=dig_writer_instance),
        ),
        patch.object(bender_functions, 'daq', MagicMock()),
        patch.object(bender_functions, 'TerminalConfiguration', MagicMock()),
        patch.object(bender_daq_kill, 'daq_emergency_stop', _fake_emergency_stop),
    ]


def test_two_consecutive_dynamic_runs_each_reset_device():
    b = _minimal_dynamic_bender()
    reset_calls: list = []
    patches = _patched_daq_run(reset_calls)
    for p in patches:
        p.start()
    try:
        b.run_experiment(test_type='dynamic')
        b.run_experiment(test_type='dynamic')
    finally:
        for p in patches:
            p.stop()

    # The device-reset teardown must fire on EVERY run, not just the first.
    assert len(reset_calls) >= 2, (
        f'expected device reset on each run, saw {len(reset_calls)} reset(s)'
    )
    # It should target the configured device each time.
    assert all(c['device'] == b.device_name for c in reset_calls[:2])


def test_run_finally_resets_even_when_acquisition_raises():
    """A failing acquisition must still reset the device in the finally -- AND release the
    motor (whole-port TRISTATE power-up). With the energized power-up HIGH now effective on
    the device, a failure would otherwise leave the motor powered/holding; the documented
    contract is that failure leaves the motor de-energized and safe."""
    b = _minimal_dynamic_bender()
    reset_calls: list = []
    patches = _patched_daq_run(reset_calls)

    class _ExplodingTask(_CtxTask):
        def __init__(self, *a, **k):
            super().__init__(*a, **k)
            self.wait_until_done = MagicMock(side_effect=RuntimeError('boom'))

    # Swap the Task patch for one whose acquisition raises.
    patches[0] = patch.object(bender_functions, 'Task', _ExplodingTask)
    for p in patches:
        p.start()
    try:
        with pytest.raises(RuntimeError, match='boom'):
            b.run_experiment(test_type='dynamic')
    finally:
        for p in patches:
            p.stop()

    assert len(reset_calls) >= 1, 'device reset must run even when acquisition fails'
    # The finally's reset on the failed run must carry the motor-release option so the
    # driver ends de-energized (failure-safe), not held by the energized power-up state.
    assert reset_calls[-1]['release'], (
        f"failed run must release the motor in its final reset, saw {reset_calls[-1]}"
    )


def test_stimburstdur_no_divide_warning_when_stim_disabled():
    b = _minimal_dynamic_bender()
    b.is_stim = False
    b.stim_pulse_rate = 0.0
    with warnings.catch_warnings():
        warnings.simplefilter('error', RuntimeWarning)
        with np.errstate(divide='raise', invalid='raise'):
            b._organize_cycles_for_dynamic_run()
    sbd = np.asarray(b.stimburstdur, dtype=float)
    assert sbd.size > 0
    assert np.all(np.isfinite(sbd))
    assert np.all(sbd == 0.0)
