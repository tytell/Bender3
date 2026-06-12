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

from bender_functions import Bender, DIR_HOLD_AFTER_STEP_SAMPLES

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


def _iso_block(b, th0, vel, approach_from_deg, post_baseline_s=1.0, mirror_stim_side=None):
    return b._isovelocity_one_block(
        th0, vel,
        pre_hold_s=0.3, iso_duration_s=0.2,
        stim_onset_s=None, stim_duration_s=None, is_stim=False,
        spr=75.0, stim_voltage=5.0, daq_hz=DAQ_HZ,
        recruitment='bilateral_simultaneous', sequential_left_frac=0.5,
        mirror_stim_side=mirror_stim_side, post_baseline_s=post_baseline_s,
        approach_from_deg=approach_from_deg,
    )


@pytest.mark.parametrize('th0', [0.0, 5.0, 10.0, -8.0])
def test_isovelocity_from_rest_returns_to_center(th0):
    """An isovelocity block that starts from rest (approach_from_deg=0) must begin and
    end at 0 and net zero steps -- no leftward 'kickback' by the starting strain."""
    b = _bender()
    b.daq_ai_sample_rate_hz = DAQ_HZ
    blk = _iso_block(b, th0, 20.0, approach_from_deg=0.0)
    assert blk['angle'][0] == pytest.approx(0.0)
    assert blk['angle'][-1] == pytest.approx(0.0)
    _, ms, md = _pulses(b, blk['t'], blk['angle'], blk['anglevel'])
    assert _net_steps(ms, md) == 0


def test_isovelocity_continuation_segment_has_no_approach():
    """The bilateral-mirror second segment continues from a mid-trajectory angle, so
    approach_from_deg=None must leave the timeline starting at th0 (no spurious 0->th0)."""
    b = _bender()
    b.daq_ai_sample_rate_hz = DAQ_HZ
    blk = _iso_block(b, 7.5, -20.0, approach_from_deg=None, mirror_stim_side='right')
    assert blk['angle'][0] == pytest.approx(7.5)


def test_isovelocity_multi_step_no_accumulated_drift():
    """Several from-rest isovelocity steps in a row each net zero, so the prep does not
    walk off center across steps."""
    b = _bender()
    b.daq_ai_sample_rate_hz = DAQ_HZ
    total = 0
    for vel in (8.0, 16.0, 24.0):
        blk = _iso_block(b, 6.0, vel, approach_from_deg=0.0)
        _, ms, md = _pulses(b, blk['t'], blk['angle'], blk['anglevel'])
        total += _net_steps(ms, md)
    assert total == 0


def test_isovelocity_nonzero_start_requires_positive_pre_hold():
    """A non-zero starting angle with no pre-hold time cannot ramp into position; fail loud."""
    b = _bender()
    b.daq_ai_sample_rate_hz = DAQ_HZ
    with pytest.raises(ValueError, match='pre_hold_s must be > 0'):
        b._timeline_prehold_isovelocity(5.0, 20.0, 0.0, 0.2, DAQ_HZ, approach_from_deg=0.0)


def test_direction_convention_unchanged():
    """A pure forward ramp emits only forward (bit 0) steps; a pure reverse ramp emits
    only reverse (bit 1) steps."""
    b = _bender()
    t, angle, anglevel = _ramp_hold(b, 0.0, 10.0, ramp_s=0.5, hold_s=0.0)
    _, ms, md = _pulses(b, t, angle, anglevel)
    assert np.all(md[ms == 1] == 0)
    assert _net_steps(ms, md) == round(10.0 * GEAR / STEPSIZE)

    # Second ramp is an independent trajectory starting from neutral, not a continuation
    # of the first. Re-anchor the cross-segment step-quantization phase the same way a real
    # protocol does at its start (command_start_position_zero), so the carry does not treat
    # the +10 deg endpoint as the starting point of the reverse demo.
    b._motor_continuous_step_pos = None
    t, angle, anglevel = _ramp_hold(b, 0.0, -10.0, ramp_s=0.5, hold_s=0.0)
    _, ms, md = _pulses(b, t, angle, anglevel)
    assert np.all(md[ms == 1] == 1)
    assert _net_steps(ms, md) == -round(10.0 * GEAR / STEPSIZE)


# --- Cross-segment sub-step carry (Fix 2: no fractional-step accumulation) -------------

def _flat_angle_for_microsteps(microsteps: float) -> float:
    """Specimen angle (deg) whose commanded motor position equals ``microsteps`` microsteps."""
    return float(microsteps) * STEPSIZE / GEAR


def test_first_segment_has_no_boundary_step():
    """A fresh phase (None) reproduces the old ``motorstep[0]=0``: the protocol's first
    segment never emits a spurious boundary step (no carry to apply yet)."""
    b = _bender()
    b._motor_continuous_step_pos = None
    angle = np.full(3, _flat_angle_for_microsteps(0.6))   # commanded = 0.6 microsteps
    t = np.array([0.0, 0.001, 0.002])
    _, ms, md = _pulses(b, t, angle, np.zeros(3))
    assert ms[0] == 0
    assert _net_steps(ms, md) == 0


def test_carry_recovers_boundary_step_deferred_one_sample():
    """With a carried sub-step phase, the inter-segment catch-up step the old forced
    ``motorstep[0]=0`` silently dropped is still emitted -- but DEFERRED one AO sample so the
    held idle DIR gets setup lead before the first STEP edge (Task B). The step is moved, not
    dropped: net steps stay 1 and the carried phase still advances.

    Previous segment ended at 0.4 microsteps (emitted round(0.4)=0). This segment commands a
    flat 0.6 microsteps; crossing the 0.5 boundary clocks exactly one forward microstep, now at
    sample 1 (sample 0 holds for DIR lead).
    """
    b = _bender()
    b._motor_continuous_step_pos = 0.4
    angle = np.full(3, _flat_angle_for_microsteps(0.6))
    t = np.array([0.0, 0.001, 0.002])
    _, ms, md = _pulses(b, t, angle, np.zeros(3))
    assert ms[0] == 0                      # boundary sample emits no step (DIR setup lead)
    assert ms[1] == 1                      # catch-up step deferred to the next sample
    assert md[1] == 0                      # forward (bit 0)
    assert _net_steps(ms, md) == 1
    # Phase advanced to the new commanded continuous position for the next segment.
    assert b._motor_continuous_step_pos == pytest.approx(0.6)


def test_continuous_chain_emits_round_of_commanded_no_accumulation():
    """Across a chain of continuous hold-ended segments on one Bender (mirrors the per-run()
    protocol path), cumulative emitted steps equal round(commanded absolute position) at every
    segment boundary. This is the driftless invariant the carry guarantees."""
    b = _bender()
    b._motor_continuous_step_pos = None
    chain = [0.0, 3.3, 7.7, 1.1, 9.9, 0.4, 5.5]
    cumulative = 0
    for a0, a1 in zip(chain[:-1], chain[1:]):
        t, angle, anglevel = _ramp_hold(b, a0, a1, ramp_s=0.4, hold_s=0.3)
        _, ms, md = _pulses(b, t, angle, anglevel)
        cumulative += _net_steps(ms, md)
        expected = round(a1 * GEAR / STEPSIZE) - round(chain[0] * GEAR / STEPSIZE)
        assert cumulative == expected, f"after segment ending {a1} deg: {cumulative} != {expected}"


# --- Task B: energized between-ramp hold must not creep ---------------------------------

def test_trailing_hold_keeps_last_step_direction_reverse():
    """After the final STEP of a reverse ramp, the trailing hold samples keep DIR=reverse (1)
    instead of snapping to forward (0): an energized hold must not toggle DIR."""
    b = _bender()
    t, angle, anglevel = _ramp_hold(b, 0.0, -10.0, ramp_s=0.5, hold_s=1.0)
    _, ms, md = _pulses(b, t, angle, anglevel)
    last_step = int(np.flatnonzero(ms == 1)[-1])
    assert md[last_step] == 1                       # last step was reverse
    assert np.all(md[last_step + 1:] == 1)          # trailing hold holds reverse, no snap to 0


def test_trailing_hold_keeps_last_step_direction_forward():
    """Forward ramp: the trailing hold keeps DIR=forward (0)."""
    b = _bender()
    t, angle, anglevel = _ramp_hold(b, 0.0, 10.0, ramp_s=0.5, hold_s=1.0)
    _, ms, md = _pulses(b, t, angle, anglevel)
    last_step = int(np.flatnonzero(ms == 1)[-1])
    assert md[last_step] == 0
    assert np.all(md[last_step + 1:] == 0)


def test_last_motor_direction_bit_tracks_final_dir():
    """make_motor_stepper_pulses records the final DO DIR bit so run() can reassert it."""
    b = _bender()
    _pulses(b, *_ramp_hold(b, 0.0, 10.0, ramp_s=0.5, hold_s=0.5))
    assert b._last_motor_direction_bit == 0
    b._motor_continuous_step_pos = None
    _pulses(b, *_ramp_hold(b, 0.0, -10.0, ramp_s=0.5, hold_s=0.5))
    assert b._last_motor_direction_bit == 1


def test_hold_only_segment_keeps_previous_direction():
    """A segment with no steps (pure hold) keeps the previous segment's DIR -- it must never
    snap the energized hold to forward when the last motion was reverse."""
    b = _bender()
    # Reverse ramp leaves the last direction bit at reverse (1).
    _pulses(b, *_ramp_hold(b, 0.0, -10.0, ramp_s=0.5, hold_s=0.2))
    assert b._last_motor_direction_bit == 1
    # A flat hold (no commanded change -> no steps) must hold DIR=reverse, not snap to forward.
    t = np.array([0.0, 0.001, 0.002, 0.003])
    angle = np.full(t.size, -10.0)
    _, ms, md = _pulses(b, t, angle, np.zeros(t.size))
    assert np.all(ms == 0)                          # no step pulses during the hold
    assert np.all(md == 1)                          # DIR held at reverse


# --- DIR hold-after-step: a reversal's DIR flip must not land on a STEP falling edge ----

def test_dir_holds_after_last_step_at_reversal():
    """At a turnaround, DIR keeps the previous step's direction for at least
    DIR_HOLD_AFTER_STEP_SAMPLES samples after that step's pulse: flipping DIR on the very
    next sample coincides with the STEP falling edge and mis-latches the last step before
    a reversal on edge-sensitive drive inputs (~2-microsteps-per-cycle dynamic drift)."""
    b = _bender()
    t, angle, anglevel = _concat(
        _ramp_hold(b, 0.0, 20.0, ramp_s=0.4, hold_s=1.0),
        _ramp_hold(b, 20.0, 0.0, ramp_s=0.4, hold_s=1.0),
    )
    _, ms, md = _pulses(b, t, angle, anglevel)
    step_idx = np.flatnonzero(ms == 1)
    rev_steps = step_idx[md[step_idx] == 1]
    fwd_steps = step_idx[md[step_idx] == 0]
    assert rev_steps.size > 0 and fwd_steps.size > 0
    first_rev = int(rev_steps[0])
    last_fwd_before = int(fwd_steps[fwd_steps < first_rev][-1])
    for off in range(1, DIR_HOLD_AFTER_STEP_SAMPLES + 1):
        assert md[last_fwd_before + off] == md[last_fwd_before], (
            f"DIR flipped {off} sample(s) after the last forward step at {last_fwd_before} "
            f"(falling-edge race window)"
        )
    # LEAD is preserved: DIR is already reverse on the sample before the first reverse step.
    assert md[first_rev - 1] == md[first_rev] == 1


def test_dir_never_flips_within_hold_window_sine():
    """Across a multi-cycle sine (the dynamic waveform shape), DIR never changes on the
    DIR_HOLD_AFTER_STEP_SAMPLES samples following any STEP pulse."""
    b = _bender()
    tt = np.arange(0.0, 4.0, 1.0 / DAQ_HZ)
    angle = 6.0 * np.sin(2 * np.pi * 1.0 * tt)
    anglevel = np.gradient(angle, tt)
    _, ms, md = _pulses(b, tt, angle, anglevel)
    step_idx = np.flatnonzero(ms == 1)
    assert step_idx.size > 0
    n = md.size
    for k in step_idx:
        for off in range(1, DIR_HOLD_AFTER_STEP_SAMPLES + 1):
            if k + off < n:
                assert md[k + off] == md[k], (
                    f"DIR flipped {off} sample(s) after step at {k}"
                )


def test_sine_cycles_net_zero_steps():
    """Each full sine cycle nets zero signed steps (the per-cycle drift invariant the
    DIR hold protects at the drive input is already true at the emission level)."""
    b = _bender()
    # End exactly on the cycle boundary (inclusive endpoint) so the commanded final
    # position is exactly 0; an arange-exclusive end lands mid-cycle and legitimately
    # rounds to a nonzero final step.
    tt = np.arange(int(4.0 * DAQ_HZ) + 1) / DAQ_HZ
    angle = 6.0 * np.sin(2 * np.pi * 1.0 * tt)
    anglevel = np.gradient(angle, tt)
    _, ms, md = _pulses(b, tt, angle, anglevel)
    assert _net_steps(ms, md) == 0


def test_reversed_continuation_first_step_has_dir_lead():
    """When a new segment reverses relative to the held idle DIR, the first STEP is deferred so
    DIR is settled (>=1 sample) before the step edge -- the boundary missed/extra-step source."""
    b = _bender()
    # First segment ends forward, leaving the carried phase just shy of a boundary so the next
    # (reverse) segment would otherwise clock its catch-up step on sample 0.
    b._motor_continuous_step_pos = 0.6                # emitted round(0.6)=1
    angle = np.full(3, _flat_angle_for_microsteps(0.4))  # commanded 0.4 -> round 0 (reverse step)
    t = np.array([0.0, 0.001, 0.002])
    _, ms, md = _pulses(b, t, angle, np.zeros(3))
    assert ms[0] == 0                               # sample 0 holds: DIR setup lead
    assert ms[1] == 1                               # reverse catch-up step deferred to sample 1
    assert md[1] == 1                               # reverse (bit 1)
    assert _net_steps(ms, md) == -1                 # step moved, not dropped
