"""Direction bit in make_motor_stepper_pulses must track the actual step (dstep),
not commanded velocity, so forward+return motion nets zero steps and no step is
clocked out against its own motion (the open-loop drift fixed here).

Hardware-free: make_motor_stepper_pulses is pure signal math and runs on Mac.
"""
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

from bender_functions import Bender

AO_HZ = 60000.0
GEAR = 5
STEPS_PER_REV = 1600
STEPSIZE = 360.0 / STEPS_PER_REV
DAQ_HZ = 1000.0


def _bender():
    return Bender('jimenez_bender_config_A')


def _pulses(b, t, angle, anglevel):
    """Run the converter and return (tout, motorstep, motordirection) as int arrays."""
    b.record_motor_signal(np.asarray(t, float), np.asarray(angle, float), np.asarray(anglevel, float))
    b.record_stim_signal(None, None)
    tout, dig, motorstep, motordirection = b.make_motor_stepper_pulses(
        daq_ao_do_sample_rate_hz=AO_HZ,
        motor_gear_ratio=GEAR,
        motor_full_steps_per_rev=STEPS_PER_REV,
    )
    return tout, motorstep.astype(int), motordirection.astype(int)


def _net_steps(motorstep, motordirection):
    """Signed step total: bit 0 = forward (+1), bit 1 = reverse (-1)."""
    signed = np.where(motordirection == 0, 1, -1) * motorstep
    return int(signed.sum())


def _ramp_hold(b, a0, a1, ramp_s=0.5, hold_s=1.0):
    return b._timeline_ramp_hold(a0, a1, ramp_s, hold_s, DAQ_HZ)


def _concat(*segs):
    t_parts, a_parts, v_parts = [], [], []
    t_off = 0.0
    for i, (t, a, v) in enumerate(segs):
        t = np.asarray(t, float)
        if i == 0:
            t_parts.append(t)
            a_parts.append(np.asarray(a, float))
            v_parts.append(np.asarray(v, float))
        else:
            t_parts.append(t[1:] + t_off)
            a_parts.append(np.asarray(a, float)[1:])
            v_parts.append(np.asarray(v, float)[1:])
        t_off = t_parts[-1][-1]
    return np.concatenate(t_parts), np.concatenate(a_parts), np.concatenate(v_parts)


def _assert_every_pulse_matches_dstep(t, angle, motorstep, motordirection):
    """At every emitted step, the direction bit equals the sign of that step's dstep."""
    tout = np.arange(t[0], t[-1], 1.0 / AO_HZ)
    poshi = np.interp(tout, t, angle) * GEAR
    stepnum = np.round(poshi / STEPSIZE)
    dstep = np.diff(stepnum)
    step_idx = np.flatnonzero(motorstep == 1)
    assert step_idx.size > 0, "expected at least one step pulse"
    for k in step_idx:
        expected = 1 if dstep[k - 1] < 0 else 0
        assert motordirection[k] == expected, (
            f"step at sample {k} tagged dir={motordirection[k]} but dstep={dstep[k-1]}"
        )


def test_forward_return_cycle_nets_zero():
    b = _bender()
    t, angle, anglevel = _concat(_ramp_hold(b, 0.0, 30.0), _ramp_hold(b, 30.0, 0.0))
    tout, motorstep, motordirection = _pulses(b, t, angle, anglevel)
    assert _net_steps(motorstep, motordirection) == 0
    _assert_every_pulse_matches_dstep(t, angle, motorstep, motordirection)


def test_dir_settled_before_each_step():
    """DIR leads every STEP edge: with a hold at the apex, the idle sample before a
    step already carries that step's direction (guards the ea54d45 setup-timing bug)."""
    b = _bender()
    t, angle, anglevel = _concat(
        _ramp_hold(b, 0.0, 20.0, ramp_s=0.4, hold_s=1.5),
        _ramp_hold(b, 20.0, 0.0, ramp_s=0.4, hold_s=1.5),
    )
    tout, motorstep, motordirection = _pulses(b, t, angle, anglevel)
    step_idx = np.flatnonzero(motorstep == 1)
    for k in step_idx:
        if k == 0:
            continue
        assert motordirection[k - 1] == motordirection[k], (
            f"DIR not settled before step at sample {k}: "
            f"prev={motordirection[k-1]}, here={motordirection[k]}"
        )
    # The first forward step must be preceded by forward (0), never the stale REVERSE
    # that velhi<=0 produced during the zero-velocity start.
    first_step = int(step_idx[0])
    assert motordirection[first_step] == 0
    assert motordirection[max(0, first_step - 1)] == 0


def test_isovelocity_like_ramp_and_return_balance():
    """Rightward constant-velocity ramp then linear return to 0 nets zero (no 'left of
    zero' residual). Symmetric for the leftward sign."""
    b = _bender()
    for sign in (+1.0, -1.0):
        peak = 15.0 * sign
        # constant-velocity ramp 0 -> peak (no hold), then linear return peak -> 0
        ramp = _ramp_hold(b, 0.0, peak, ramp_s=0.5, hold_s=0.0)
        ret = _ramp_hold(b, peak, 0.0, ramp_s=0.5, hold_s=0.0)
        t, angle, anglevel = _concat(ramp, ret)
        tout, motorstep, motordirection = _pulses(b, t, angle, anglevel)
        assert _net_steps(motorstep, motordirection) == 0, f"sign={sign}"
        _assert_every_pulse_matches_dstep(t, angle, motorstep, motordirection)


def test_multi_step_isometric_accumulation_nets_zero():
    """Many forward-ramp/hold segments each followed by a separate reset-to-0 segment
    (mirrors the per-run() isometric path) must sum to zero net steps across all steps."""
    b = _bender()
    targets = [5.0, 12.0, 19.0, 26.0, 8.0]
    total = 0
    for tgt in targets:
        # forward ramp + hold (one acquisition segment), ending at tgt
        t, angle, anglevel = _ramp_hold(b, 0.0, tgt, ramp_s=0.4, hold_s=1.0)
        _, ms, md = _pulses(b, t, angle, anglevel)
        total += _net_steps(ms, md)
        # separate neutral-reset segment tgt -> 0
        t, angle, anglevel = _ramp_hold(b, tgt, 0.0, ramp_s=0.4, hold_s=0.0)
        _, ms, md = _pulses(b, t, angle, anglevel)
        total += _net_steps(ms, md)
    assert total == 0


def test_direction_convention_unchanged():
    """A pure forward ramp emits only forward (bit 0) steps; a pure reverse ramp emits
    only reverse (bit 1) steps."""
    b = _bender()
    t, angle, anglevel = _ramp_hold(b, 0.0, 10.0, ramp_s=0.5, hold_s=0.0)
    _, ms, md = _pulses(b, t, angle, anglevel)
    assert np.all(md[ms == 1] == 0)
    assert _net_steps(ms, md) == round(10.0 * GEAR / STEPSIZE)

    t, angle, anglevel = _ramp_hold(b, 0.0, -10.0, ramp_s=0.5, hold_s=0.0)
    _, ms, md = _pulses(b, t, angle, anglevel)
    assert np.all(md[ms == 1] == 1)
    assert _net_steps(ms, md) == -round(10.0 * GEAR / STEPSIZE)
