"""Unit tests for isometric / isovelocity block-sequence helpers."""
from __future__ import annotations

import numpy as np
import pytest

from bender_functions import Bender


class _BlockHelperBender(Bender):
    """Minimal instance for block routing helpers (no hardware config import)."""

    def __init__(self):
        self.S1side = 'left'
        self.S2side = 'right'
        self.positive_motor_direction = 'right'
        self.specimen_lateral_index_on_positive_motor_side = -1
        self.specimen_side_index_left = -1
        self.specimen_side_index_right = 1
        self.stim_pulse_rate = 75.0


@pytest.fixture
def b():
    return _BlockHelperBender()


def test_normalize_block_direction(b):
    assert b._normalize_block_direction('left') == 'left'
    assert b._normalize_block_direction('RIGHT') == 'right'
    with pytest.raises(ValueError, match='left.*right'):
        b._normalize_block_direction('center')


def test_normalize_stim_sides(b):
    assert b._normalize_stim_sides(None) == 'left'
    assert b._normalize_stim_sides('both') == 'both'
    assert b._normalize_stim_sides('Bilateral') == 'both'
    # OFF aliases (shared by isometric and isovelocity per-block routing).
    assert b._normalize_stim_sides('off') == 'off'
    assert b._normalize_stim_sides('OFF') == 'off'
    assert b._normalize_stim_sides('none') == 'off'
    assert b._normalize_stim_sides('no_stim') == 'off'
    with pytest.raises(ValueError, match='left.*right.*both'):
        b._normalize_stim_sides('sequential')


def test_normalize_block_sequence_empty_is_legacy(b):
    assert b._normalize_block_sequence(None) is None
    assert b._normalize_block_sequence([]) is None


def test_normalize_block_sequence_valid(b):
    seq = b._normalize_block_sequence([
        {'direction': 'left', 'stim_sides': 'right'},
        {'direction': 'right', 'stim_sides': 'both'},
    ])
    assert seq == [
        {'direction': 'left', 'stim_sides': 'right'},
        {'direction': 'right', 'stim_sides': 'both'},
    ]


def test_stim_sides_to_recruitment(b):
    assert b._stim_sides_to_recruitment('left') == 'left_unilateral'
    assert b._stim_sides_to_recruitment('right') == 'right_unilateral'
    assert b._stim_sides_to_recruitment('both') == 'bilateral_simultaneous'
    # OFF deposits no stim, so the recruitment label is moot; falls through to bilateral.
    assert b._stim_sides_to_recruitment('off') == 'bilateral_simultaneous'


def test_validate_block_sequence_voltages(b):
    seq = [{'direction': 'left', 'stim_sides': 'left'}]
    b._validate_block_sequence_voltages(seq, 5.0, 5.0)
    with pytest.raises(ValueError, match='left_stim_voltage'):
        b._validate_block_sequence_voltages(seq, 0.0, 5.0)
    seq_both = [{'direction': 'left', 'stim_sides': 'both'}]
    with pytest.raises(ValueError, match='right_stim_voltage'):
        b._validate_block_sequence_voltages(seq_both, 5.0, -1.0)


def test_validate_block_sequence_voltages_off_needs_no_voltage(b):
    """An OFF block requires no stim voltage (zero/blank voltages are fine)."""
    seq_off = [{'direction': 'left', 'stim_sides': 'off'}]
    b._validate_block_sequence_voltages(seq_off, 0.0, 0.0)
    # OFF mixed with LEFT still requires the LEFT voltage only.
    seq_mixed = [
        {'direction': 'left', 'stim_sides': 'off'},
        {'direction': 'right', 'stim_sides': 'left'},
    ]
    b._validate_block_sequence_voltages(seq_mixed, 5.0, 0.0)
    with pytest.raises(ValueError, match='left_stim_voltage'):
        b._validate_block_sequence_voltages(seq_mixed, 0.0, 0.0)


def test_normalize_block_sequence_with_off(b):
    """OFF is a valid per-block stim_sides for the shared isometric/isovelocity block sequence."""
    seq = b._normalize_block_sequence([
        {'direction': 'left', 'stim_sides': 'off'},
        {'direction': 'right', 'stim_sides': 'both'},
    ])
    assert seq == [
        {'direction': 'left', 'stim_sides': 'off'},
        {'direction': 'right', 'stim_sides': 'both'},
    ]


def test_route_stim_sides_volts_left_only(b):
    t = np.linspace(0, 0.1, 200)
    active = (t >= 0.02) & (t < 0.08)
    s1, s2 = b._route_stim_sides_volts(t, active, 75.0, 'left', 4.0, 6.0)
    assert np.max(np.abs(s1)) == pytest.approx(4.0)
    assert np.max(np.abs(s2)) == 0.0


def test_route_stim_sides_volts_both_different_voltages(b):
    t = np.linspace(0, 0.1, 200)
    active = (t >= 0.02) & (t < 0.08)
    s1, s2 = b._route_stim_sides_volts(t, active, 75.0, 'both', 3.5, 7.0)
    assert np.max(np.abs(s1)) == pytest.approx(3.5)
    assert np.max(np.abs(s2)) == pytest.approx(7.0)


def test_route_stim_sides_volts_off_is_zero(b):
    """OFF yields all-zero S1/S2 even with an active window and nonzero voltages.

    This is the single routing function both isometric (``_run_force_length_steps``) and
    isovelocity (``_isovelocity_one_block``) call for per-block stim, so OFF behaves identically
    in both protocols.
    """
    t = np.linspace(0, 0.1, 200)
    active = (t >= 0.02) & (t < 0.08)
    s1, s2 = b._route_stim_sides_volts(t, active, 75.0, 'off', 5.0, 5.0)
    assert np.count_nonzero(s1) == 0
    assert np.count_nonzero(s2) == 0
    assert s1.shape == t.shape and s2.shape == t.shape


def test_preview_append_neutral_reset_skips_noop():
    """At neutral with zero ramp, preview must not append a spurious reset segment."""
    from bender_gui_preview import _preview_append_neutral_reset

    b = _BlockHelperBender()
    t_chunks, a_chunks, w_chunks, s1_chunks, s2_chunks = [], [], [], [], []
    toff, last_deg = _preview_append_neutral_reset(
        b, 0.0, 0.0, 1000.0, t_chunks, a_chunks, w_chunks, s1_chunks, s2_chunks, 0.0,
    )
    assert toff == 0.0
    assert last_deg == 0.0
    assert t_chunks == []


def test_preview_append_neutral_reset_skips_zero_displacement_with_ramp():
    """Already-neutral reset with a positive ramp is still degenerate (0°→0°); must skip."""
    from bender_gui_preview import _preview_append_neutral_reset

    b = _BlockHelperBender()
    t_chunks, a_chunks, w_chunks, s1_chunks, s2_chunks = [], [], [], [], []
    toff, last_deg = _preview_append_neutral_reset(
        b, 0.0, 2.0, 1000.0, t_chunks, a_chunks, w_chunks, s1_chunks, s2_chunks, 0.0,
    )
    assert toff == 0.0
    assert last_deg == 0.0
    assert t_chunks == []


def test_run_neutral_reset_segment_skips_noop(b):
    """At neutral with zero ramp, backend must not run DAQ for a reset segment."""
    b.daq_ai_sample_rate_hz = 1000.0
    b.daq_ao_do_sample_rate_hz = 1000.0
    result = b._run_neutral_reset_segment(0.0, 0.0, 'mock_device')
    assert result == 0.0


def test_run_neutral_reset_segment_skips_zero_displacement_with_ramp(b):
    """Already-neutral reset with a positive ramp is degenerate; backend must skip DAQ."""
    b.daq_ai_sample_rate_hz = 1000.0
    b.daq_ao_do_sample_rate_hz = 1000.0
    result = b._run_neutral_reset_segment(0.0, 2.0, 'mock_device')
    assert result == 0.0


def test_run_neutral_reset_segment_requires_positive_ramp_for_real_reset(b):
    """A needed reset (motor off neutral) cannot take zero time; must raise, not silently skip."""
    b.daq_ai_sample_rate_hz = 1000.0
    b.daq_ao_do_sample_rate_hz = 1000.0
    with pytest.raises(ValueError, match='must be > 0 s'):
        b._run_neutral_reset_segment(12.5, 0.0, 'mock_device')


def test_preview_append_neutral_reset_requires_positive_ramp_for_real_reset():
    """Preview mirrors the backend: a real reset with non-positive ramp raises."""
    from bender_gui_preview import _preview_append_neutral_reset

    b = _BlockHelperBender()
    t_chunks, a_chunks, w_chunks, s1_chunks, s2_chunks = [], [], [], [], []
    with pytest.raises(ValueError, match='must be > 0 s'):
        _preview_append_neutral_reset(
            b, 12.5, 0.0, 1000.0, t_chunks, a_chunks, w_chunks, s1_chunks, s2_chunks, 0.0,
        )


def test_resolve_stim_onset_duration_migrates_legacy_settle(b):
    onset, dur = b._resolve_stim_onset_duration_s(
        {'settle_before_stim_s': 0.5, 'stim_duration_s': 2.0},
        segment_duration_s=5.0,
    )
    assert onset == pytest.approx(0.5)
    assert dur == pytest.approx(2.0)


def test_resolve_stim_onset_duration_migrates_pre_iso(b):
    onset, dur = b._resolve_stim_onset_duration_s(
        {'pre_iso_stim_duration_s': 0.1, 'stim_duration_s': 0.15},
        segment_duration_s=0.2,
    )
    assert onset == pytest.approx(-0.1)
    assert dur == pytest.approx(0.15)


def test_isometric_negative_onset_budget_is_ramp_duration(b):
    """Isometric negative onset within the ramp passes unchanged; beyond it is auto-clamped."""
    ramp_s = 2.0
    hold_s = 5.0
    # Onset -1.5 s starts 1.5 s into a 2 s ramp: within budget, no clamp.
    sp_ok = {'is_stim': True, 'stim_onset_s': -1.5, 'stim_duration_s': 1.0}
    notices = b._validate_stim_timing_for_steps(
        sp_ok,
        test_type='isometric',
        num_steps=3,
        pre_hold_at_start_s=ramp_s,
        segment_duration_s=hold_s,
    )
    assert notices == []
    assert sp_ok['stim_onset_s'] == pytest.approx(-1.5)
    # Onset -2.5 s would start before the ramp begins (budget exceeded): auto-clamp to -ramp.
    sp_clamp = {'is_stim': True, 'stim_onset_s': -2.5, 'stim_duration_s': 1.0}
    notices = b._validate_stim_timing_for_steps(
        sp_clamp,
        test_type='isometric',
        num_steps=3,
        pre_hold_at_start_s=ramp_s,
        segment_duration_s=hold_s,
    )
    assert notices  # clamping reported
    assert sp_clamp['stim_onset_s'] == pytest.approx(-ramp_s)


def test_isometric_segment_duration_is_hold_duration_s(b):
    """The isometric active segment used for clamping is hold_duration_s from the stim params."""
    ramp_s = 2.0
    # Stim window onset 0 + duration 4 s fits inside a 5 s hold: no clamp.
    sp_ok = {'is_stim': True, 'stim_onset_s': 0.0, 'stim_duration_s': 4.0, 'hold_duration_s': 5.0}
    notices = b._validate_stim_timing_for_steps(
        sp_ok,
        test_type='isometric',
        num_steps=2,
        pre_hold_at_start_s=ramp_s,
        segment_duration_s=float(sp_ok.get('hold_duration_s', 5.0)),
    )
    assert notices == []
    assert sp_ok['stim_duration_s'] == pytest.approx(4.0)
    # Same window with a shorter 3 s hold now bleeds past the segment end: auto-clamp duration.
    sp_bad = {'is_stim': True, 'stim_onset_s': 0.0, 'stim_duration_s': 4.0, 'hold_duration_s': 3.0}
    notices = b._validate_stim_timing_for_steps(
        sp_bad,
        test_type='isometric',
        num_steps=2,
        pre_hold_at_start_s=ramp_s,
        segment_duration_s=float(sp_bad.get('hold_duration_s', 3.0)),
    )
    assert any('active segment' in m for m in notices)
    assert sp_bad['stim_duration_s'] == pytest.approx(3.0)


def test_validate_stim_timing_bounds_rejects_pre_hold_bleed(b):
    with pytest.raises(ValueError, match='before the allowed pre-hold'):
        b._validate_stim_timing_bounds(
            step_index=1,
            stim_onset_s=-0.5,
            stim_duration_s=0.1,
            pre_hold_at_start_s=0.3,
            segment_duration_s=0.2,
            protocol_label='isovelocity',
        )


def test_validate_stim_timing_bounds_rejects_segment_bleed(b):
    with pytest.raises(ValueError, match='past the active segment'):
        b._validate_stim_timing_bounds(
            step_index=2,
            stim_onset_s=0.1,
            stim_duration_s=0.5,
            pre_hold_at_start_s=0.3,
            segment_duration_s=0.4,
            protocol_label='isovelocity',
        )


def test_validate_stim_timing_bounds_accepts_valid_window(b):
    b._validate_stim_timing_bounds(
        step_index=1,
        stim_onset_s=-0.2,
        stim_duration_s=0.15,
        pre_hold_at_start_s=0.3,
        segment_duration_s=0.2,
        protocol_label='isovelocity',
    )


def test_lateral_index_for_block_direction(b):
    assert b._lateral_index_for_block_direction('left') == b.specimen_side_index_left
    assert b._lateral_index_for_block_direction('right') == b.specimen_side_index_right
