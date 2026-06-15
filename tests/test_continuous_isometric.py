"""Continuous isometric acquisition: one finite task for all steps/blocks.

Pins the contract of the rewritten ``_run_force_length_steps`` /
``_build_isometric_continuous_timeline``:

  * the stitched timeline is uniformly sampled (passes ``record_motor_signal``'s 8*dt guard)
    for many asymmetric steps,
  * per-step ``slice_bounds`` overlap by exactly one sample at each junction and keep every
    per-step array at its full local length ``Li`` (identical shapes to the old per-step flow),
  * a block run produces ``num_steps * n_blocks`` trials (NOT ``num_steps``), each tagged with
    its block direction and signed target,
  * the run starts and ends at home,
  * ``reset_between_steps`` and the inter-step gap are embedded as dead time (not in entries),
  * ``mean_xforce_stim`` is gone and ``acquisition_mode='continuous'`` is recorded,
  * stim stays sample-aligned to motion in every per-step entry.
"""
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

from bender_functions import Bender


DAQ_HZ = 1000.0


def _bender():
    b = Bender('jimenez_bender_config_A')  # ai = 1000 Hz
    b.dclamp = 10.0
    b.xsec_width = 2.0
    return b


def _sim_bender():
    b = _bender()
    b.session_simulated = True
    b.device_name = 'SimDev'
    return b


def _build(b, targets, **kw):
    """Call the timeline helper with sane non-stim defaults."""
    params = dict(
        ramp_duration_s=0.5,
        hold_duration_s=2.0,
        pre_baseline_s=0.0,
        post_baseline_s=0.0,
        stim_onset_s=None,
        settle_before_stim_s=0.2,
        stim_duration_s=None,
        is_stim=False,
        stim_pulse_rate_hz=75.0,
        stim_voltage=5.0,
        daq_hz=DAQ_HZ,
    )
    params.update(kw)
    return b._build_isometric_continuous_timeline(np.asarray(targets, dtype=float), **params)


def test_helper_slice_bounds_overlap_contiguous():
    b = _bender()
    t_all, angle_all, anglevel_all, s1, s2, bounds, locs = _build(b, [10.0, 20.0, 30.0])
    # Every per-step slice keeps the FULL local length Li (overlapping slices, not Li-1).
    for (s, e), loc in zip(bounds, locs):
        assert e - s == loc['t'].size
    # Junction overlap: each step (i>0) shares exactly one sample with the previous slice.
    for i in range(1, len(bounds)):
        assert bounds[i][0] == bounds[i - 1][1] - 1
    # Concrete shape for equal-length steps.
    L = locs[0]['t'].size
    assert bounds == [(0, L), (L - 1, 2 * L - 1), (2 * L - 2, 3 * L - 2)]
    # Stitched timeline is one uniform grid; stim arrays match motion length.
    assert np.allclose(np.diff(t_all), 1.0 / DAQ_HZ, atol=1e-9)
    assert s1.size == t_all.size and s2.size == t_all.size
    assert angle_all.size == t_all.size and anglevel_all.size == t_all.size


def test_helper_timeline_uniform_for_many_asymmetric_steps():
    b = _bender()
    targets = [12.0, -37.0, 5.0, -22.0, 44.0, -3.0, 28.0]
    t_all, angle_all, anglevel_all, _s1, _s2, _b, _l = _build(b, targets)
    # The uniformity guard in record_motor_signal must accept the stitched grid unchanged.
    b.record_motor_signal(t_all, angle_all, anglevel_all)
    assert b.daq_ai_sample_rate_hz == DAQ_HZ
    assert b.t.size == t_all.size
    # span == (n-1)*dt exactly (grid-snap), not merely within tolerance.
    assert (t_all[-1] - t_all[0]) == pytest.approx((t_all.size - 1) / DAQ_HZ, abs=1e-9)


def test_block_run_six_entries_shapes_and_tags():
    b = _sim_bender()
    out = b.isometric(
        5.0, 45.0, 3, mode='angle',
        stim_params={
            'block_sequence': [
                {'direction': 'right', 'stim_sides': 'left'},
                {'direction': 'left', 'stim_sides': 'right'},
            ],
            'ramp_duration_s': 0.5,
            'hold_duration_s': 2.0,
            'pre_baseline_s': 0.0,
            'post_baseline_s': 0.0,
            'is_stim': False,
        },
    )
    # 3 steps x 2 blocks -> 6 trials, NOT 3.
    assert len(out) == 6
    assert [e['step_index'] for e in out] == list(range(6))
    assert [e['trial_index'] for e in out] == list(range(6))
    for e in out:
        n = e['t'].size
        assert e['aidata'].shape[1] == n
        assert e['angle_measured'].size == n
        assert e['S1stimcmd'].size == n and e['S2stimcmd'].size == n
        assert e['anglevel_cmd'].size == n
        assert 'mean_xforce_stim' not in e
    # Each block tagged; the two blocks bend in opposite directions (signs flip).
    assert all(e['block_direction'] == 'right' for e in out[:3])
    assert all(e['block_direction'] == 'left' for e in out[3:])
    assert all(e['block_index'] == 0 for e in out[:3])
    assert all(e['block_index'] == 1 for e in out[3:])
    assert np.sign(out[0]['target_deg']) == -np.sign(out[3]['target_deg'])
    assert all(np.sign(e['target_deg']) == np.sign(out[0]['target_deg']) for e in out[:3])
    # Start at home (first step ramps from 0) and end at home (final return ramp).
    assert abs(out[0]['ramp_from_deg']) < 1e-9
    assert abs(float(b._last_commanded_angle_deg)) < 1e-6
    # Continuous-mode metadata, user-selectable keys preserved.
    md = b.h5_protocol_metadata
    assert md.get('acquisition_mode') == 'continuous'
    for k in ('rest_between_steps_s', 'reset_between_steps', 'isometric_inter_step_interval_s'):
        assert k in md


def test_reset_between_steps_embeds_return_ramps_without_changing_entry_shapes():
    b = _bender()
    targets = [20.0, 30.0, 40.0]
    t0, *_rest0, bounds0, locs0 = _build(b, targets, reset_steps=False)
    t1, *_rest1, bounds1, locs1 = _build(b, targets, reset_steps=True)
    # Per-step entry lengths are identical (reset ramps are dead time, not in any entry).
    assert [l['t'].size for l in locs0] == [l['t'].size for l in locs1]
    for (s, e), loc in zip(bounds1, locs1):
        assert e - s == loc['t'].size
    # With reset_between_steps, each step (i>0) ramps from 0.
    assert all(abs(locs1[i]['prev']) < 1e-9 for i in range(1, len(locs1)))
    # The embedded return-to-zero ramps make the stitched waveform longer.
    assert t1.size > t0.size


def test_inter_step_gap_is_dead_time_not_in_entries():
    b = _bender()
    targets = [20.0, 30.0, 40.0]
    t0, *_r0, _b0, locs0 = _build(b, targets, inter_step_hold_s=0.0)
    t1, *_r1, _b1, locs1 = _build(b, targets, inter_step_hold_s=1.0)
    # Gap holds do not change per-step entry shapes...
    assert [l['t'].size for l in locs0] == [l['t'].size for l in locs1]
    # ...but they lengthen the single acquisition (settle time preserved, recorded-then-discarded).
    assert t1.size > t0.size


def test_stim_stays_sample_aligned_to_motion():
    b = _sim_bender()
    out = b.isometric(
        10.0, 40.0, 2, mode='angle',
        stim_params={
            'block_sequence': [{'direction': 'left', 'stim_sides': 'left'}],
            'ramp_duration_s': 0.5,
            'hold_duration_s': 3.0,
            'pre_baseline_s': 0.0,
            'post_baseline_s': 0.0,
            'settle_before_stim_s': 0.5,
            'stim_duration_s': 1.0,
            'is_stim': True,
            'left_stim_voltage': 5.0,
            'right_stim_voltage': 5.0,
        },
    )
    saw_stim = False
    for e in out:
        t = e['t']
        on = (np.abs(e['S1stimcmd']) > 1e-9) | (np.abs(e['S2stimcmd']) > 1e-9)
        assert e['S1stimcmd'].size == t.size and e['S2stimcmd'].size == t.size
        if on.any():
            saw_stim = True
            assert t[on].min() >= e['stim_t0'] - 1e-6
            assert t[on].max() < e['stim_t1'] + 1e-3
    assert saw_stim, 'expected stim pulses within the active window of at least one step'
