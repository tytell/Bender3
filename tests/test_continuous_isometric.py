"""Segmented isometric acquisition: one FINITE run() per step (segmented_finite).

Pins the contract of the per-step ``_run_force_length_steps`` / ``_build_isometric_one_step``
(the conversion from the old stitched single-finite task + post-hoc slicing):

  * each step is its own recorded window, bracketed by 1 s neutral bookends at home (0 deg),
  * the step starts at home and ends at home (the return-to-neutral is recorded INSIDE the step),
  * ``step_index`` is ONE-based and contiguous; a block run yields ``num_steps * n_blocks`` trials,
    each tagged with its block direction and signed target,
  * ``daq_collection_type == 'segmented'`` and ``protocol_sampling_mode == 'segmented_finite'``,
  * per-step timing fields (``wall_clock_start``, ``duration_second``, ``rest_before_second``) and
    ``operating_point`` (= ``target_deg``) / ``operating_point_units`` (``'degree'``) are recorded,
  * ``mean_xforce_stim`` is gone and stim stays sample-aligned to motion within each step.
"""
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

from bender_functions import Bender, SEGMENTED_STEP_BUFFER_S


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


# --- single-step builder (no DAQ; pure timeline) ---------------------------------------------

def _build_step(b, target, **kw):
    params = dict(
        prev_deg=0.0,
        ramp_duration_s=0.5,
        hold_duration_s=2.0,
        pre_baseline_s=0.0,
        post_baseline_s=0.0,
        stim_onset_s=None,
        settle_before_stim_s=0.2,
        stim_duration_s=None,
        is_stim=False,
        spr=75.0,
        stim_voltage=5.0,
        daq_hz=DAQ_HZ,
        post_buffer_s=SEGMENTED_STEP_BUFFER_S,
    )
    params.update(kw)
    return b._build_isometric_one_step(float(target), **params)


def test_one_step_starts_and_ends_at_home():
    b = _bender()
    step = _build_step(b, 30.0)
    t = step['t']
    angle = step['angle']
    assert t.size == angle.size
    # Ramp opens from home (prev_deg=0) and the recorded window ends at home after the return ramp.
    assert abs(float(angle[0])) < 1e-9
    assert abs(float(angle[-1])) < 1e-9
    # The target is actually reached somewhere in the middle (the hold).
    assert np.max(np.abs(angle)) == pytest.approx(30.0, abs=1e-6)


def test_one_step_post_bookend_is_flat_at_home():
    b = _bender()
    step = _build_step(b, 25.0, post_buffer_s=SEGMENTED_STEP_BUFFER_S)
    angle = step['angle']
    n_buf = int(round(SEGMENTED_STEP_BUFFER_S * DAQ_HZ))
    # The trailing post-bookend is a flat neutral hold at 0 deg of (approximately) step_buffer_s.
    tail = angle[-n_buf:]
    assert np.allclose(tail, 0.0, atol=1e-9)
    # Velocity is zero throughout the post-bookend (held, not moving).
    assert np.allclose(step['anglevel'][-n_buf:], 0.0, atol=1e-9)


def test_one_step_return_to_neutral_is_speed_capped():
    b = _bender()
    b.reset_max_speed_deg_per_s = 15.0
    target = 45.0
    step = _build_step(b, target, post_buffer_s=0.0)  # isolate the return ramp (no trailing hold)
    anglevel = step['anglevel']
    # The forward ramp (0 -> target) honors the operator ramp_duration_s and may be fast; only the
    # RETURN-to-neutral ramp (the negative-velocity tail, target -> 0) must be speed-capped at
    # reset_max_speed_deg_per_s, so a large target does not slew too fast and lose open-loop steps.
    return_vel = anglevel[anglevel < -1e-9]
    assert return_vel.size > 0, 'expected a recorded return-to-neutral ramp'
    assert np.max(np.abs(return_vel)) <= 15.0 + 1e-6
    # Forward ramp here is intentionally faster (90 deg/s over 0.5 s) -- it is operator-sized.
    assert np.max(anglevel) > 15.0


# --- full per-step run (simulated DAQ) -------------------------------------------------------

def _block_run(b, **overrides):
    sp = {
        'block_sequence': [
            {'direction': 'right', 'stim_sides': 'left'},
            {'direction': 'left', 'stim_sides': 'right'},
        ],
        'ramp_duration_s': 0.5,
        'hold_duration_s': 2.0,
        'pre_baseline_s': 0.0,
        'post_baseline_s': 0.0,
        'is_stim': False,
        # Keep the inter-step gap at 0 so the test does not actually sleep between run() calls.
        'inter_step_interval_s': 0.0,
    }
    sp.update(overrides)
    return b.isometric(5.0, 45.0, 3, mode='angle', stim_params=sp)


def test_block_run_six_entries_shapes_and_tags():
    b = _sim_bender()
    out = _block_run(b)
    # 3 steps x 2 blocks -> 6 trials, NOT 3.
    assert len(out) == 6
    assert [e['step_index'] for e in out] == list(range(1, 7))
    assert [e['trial_index'] for e in out] == list(range(6))
    for e in out:
        n = e['t'].size
        assert e['aidata'].shape[1] == n
        assert e['angle_measured'].size == n
        assert e['S1stimcmd'].size == n and e['S2stimcmd'].size == n
        assert e['anglevel_cmd'].size == n
        assert 'mean_xforce_stim' not in e
        # Every recorded step starts and ends at home.
        assert abs(float(e['angle_cmd'][0])) < 1e-9
        assert abs(float(e['angle_cmd'][-1])) < 1e-9
    # Each block tagged; the two blocks bend in opposite directions (signs flip).
    assert all(e['block_direction'] == 'right' for e in out[:3])
    assert all(e['block_direction'] == 'left' for e in out[3:])
    assert all(e['block_index'] == 0 for e in out[:3])
    assert all(e['block_index'] == 1 for e in out[3:])
    assert np.sign(out[0]['target_deg']) == -np.sign(out[3]['target_deg'])
    assert all(np.sign(e['target_deg']) == np.sign(out[0]['target_deg']) for e in out[:3])
    # Start at home (first step ramps from 0) and end at home (final return).
    assert abs(out[0]['ramp_from_deg']) < 1e-9
    assert abs(float(b._last_commanded_angle_deg)) < 1e-6


def test_block_run_is_segmented_finite_with_step_manifest_fields():
    b = _sim_bender()
    out = _block_run(b)
    # Gate D3 flipped: isometric is now a per-step run() loop.
    assert b.daq_collection_type == 'segmented'
    assert b.protocol_sampling_mode == 'segmented_finite'
    for e in out:
        # operating_point / units route into the step_manifest (schema 4).
        assert e['operating_point'] == pytest.approx(e['target_deg'])
        assert e['operating_point_units'] == 'degree'
        assert 'wall_clock_start' in e and e['wall_clock_start']
        assert e['duration_second'] == pytest.approx(e['t'].size / DAQ_HZ, abs=1e-9)
        assert 'rest_before_second' in e and e['rest_before_second'] >= 0.0
    # The first step has no preceding recorded step.
    assert out[0]['rest_before_second'] == pytest.approx(0.0, abs=1e-9)
    # Preserved user-selectable metadata.
    md = b.h5_protocol_metadata
    for k in ('rest_between_steps_s', 'reset_between_steps', 'isometric_inter_step_interval_s'):
        assert k in md
    # reset/hold forced on for segmented protocols.
    assert md['reset_between_steps'] is True
    assert md['hold_motor_between_steps'] is True


def test_each_step_has_one_second_neutral_bookends():
    b = _sim_bender()
    out = _block_run(b)
    n_buf = int(round(SEGMENTED_STEP_BUFFER_S * DAQ_HZ))
    for e in out:
        a = e['angle_cmd']
        # Pre-bookend: the recorded window opens with a 1 s flat neutral hold at home.
        assert np.allclose(a[:n_buf], 0.0, atol=1e-9)
        # Post-bookend: it closes with a 1 s flat neutral hold at home.
        assert np.allclose(a[-n_buf:], 0.0, atol=1e-9)
        # Event markers sit after the pre-bookend (shifted by step_buffer_s).
        assert e['t_active_start'] >= SEGMENTED_STEP_BUFFER_S - 1e-9


def test_stim_stays_sample_aligned_to_motion():
    b = _sim_bender()
    out = _block_run(
        b,
        block_sequence=[{'direction': 'left', 'stim_sides': 'left'}],
        hold_duration_s=3.0,
        settle_before_stim_s=0.5,
        stim_duration_s=1.0,
        is_stim=True,
        left_stim_voltage=5.0,
        right_stim_voltage=5.0,
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
