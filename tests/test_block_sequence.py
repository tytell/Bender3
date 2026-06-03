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


def test_validate_block_sequence_voltages(b):
    seq = [{'direction': 'left', 'stim_sides': 'left'}]
    b._validate_block_sequence_voltages(seq, 5.0, 5.0)
    with pytest.raises(ValueError, match='left_stim_voltage'):
        b._validate_block_sequence_voltages(seq, 0.0, 5.0)
    seq_both = [{'direction': 'left', 'stim_sides': 'both'}]
    with pytest.raises(ValueError, match='right_stim_voltage'):
        b._validate_block_sequence_voltages(seq_both, 5.0, -1.0)


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


def test_run_neutral_reset_segment_skips_noop(b):
    """At neutral with zero ramp, backend must not run DAQ for a reset segment."""
    b.daq_ai_sample_rate_hz = 1000.0
    b.daq_ao_do_sample_rate_hz = 1000.0
    result = b._run_neutral_reset_segment(0.0, 0.0, 'mock_device')
    assert result == 0.0


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


def test_lateral_index_for_block_direction(b):
    assert b._lateral_index_for_block_direction('left') == b.specimen_side_index_left
    assert b._lateral_index_for_block_direction('right') == b.specimen_side_index_right
