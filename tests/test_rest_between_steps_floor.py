"""Step 5 (D4): 2 s minimum on rest_between_steps_s enforced at Run for stepped protocols.

Tests the pure helper _segmented_rest_floor_error directly (no DAQ, no hardware, Mac-only).
The floor is a GUI Run-layer policy; the engine deliberately accepts 0 s gaps so tests
can call _run_force_length_steps / _run_isovelocity_steps without sleeping.
"""
from pathlib import Path
import sys

sys.path.append(str(Path(__file__).resolve().parents[1]))

from bender_streamlit_gui import _segmented_rest_floor_error, SEGMENTED_REST_FLOOR_S

_EXACT_MSG = (
    'Reset duration must be at least 2s to ensure clean recording buffers between '
    'steps. Physiological rest requirements will typically exceed this minimum.'
)


def test_floor_constant_is_two():
    assert SEGMENTED_REST_FLOOR_S == 2.0


# --- isometric blocks below floor -------------------------------------------------------

def test_isometric_below_floor_returns_exact_message():
    assert _segmented_rest_floor_error('isometric', 1.0) == _EXACT_MSG


def test_isometric_zero_returns_exact_message():
    assert _segmented_rest_floor_error('isometric', 0.0) == _EXACT_MSG


def test_isometric_just_below_floor_returns_message():
    assert _segmented_rest_floor_error('isometric', 1.9999) == _EXACT_MSG


# --- isovelocity blocks below floor -----------------------------------------------------

def test_isovelocity_below_floor_returns_exact_message():
    assert _segmented_rest_floor_error('isovelocity', 0.0) == _EXACT_MSG


def test_isovelocity_one_second_returns_message():
    assert _segmented_rest_floor_error('isovelocity', 1.5) == _EXACT_MSG


# --- at or above floor returns None -----------------------------------------------------

def test_isometric_at_floor_returns_none():
    assert _segmented_rest_floor_error('isometric', 2.0) is None


def test_isometric_above_floor_returns_none():
    assert _segmented_rest_floor_error('isometric', 5.0) is None


def test_isovelocity_at_floor_returns_none():
    assert _segmented_rest_floor_error('isovelocity', 2.0) is None


def test_isovelocity_above_floor_returns_none():
    assert _segmented_rest_floor_error('isovelocity', 10.0) is None


# --- non-stepped protocols always return None -------------------------------------------

def test_dynamic_zero_returns_none():
    assert _segmented_rest_floor_error('dynamic', 0.0) is None


def test_frequency_sweep_zero_returns_none():
    assert _segmented_rest_floor_error('frequency_sweep', 0.0) is None


def test_unknown_type_returns_none():
    assert _segmented_rest_floor_error('unknown', 0.0) is None


# --- bad value types handled gracefully (no raise, returns None) ------------------------

def test_none_value_returns_none():
    assert _segmented_rest_floor_error('isometric', None) is None


def test_string_value_returns_none():
    assert _segmented_rest_floor_error('isometric', 'not_a_number') is None
