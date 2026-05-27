"""Stress tests for :func:`bender_functions.convert_to_curvature`."""
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

from bender_functions import convert_to_curvature


@pytest.mark.parametrize(
    'mode,value,dc,xw,expect_raises',
    [
        ('strain', 0.05, 10.0, 2.0, False),
        ('strain_pct', 5.0, 10.0, 2.0, False),
        ('curvature', 0.5, 10.0, None, False),
        ('angle', 3.0, 10.0, None, False),
        ('strain_rate', 0.1, 10.0, 2.0, False),
        ('strain_pct_rate', 10.0, 10.0, 2.0, False),
        ('curvature_rate', 0.2, 10.0, None, False),
        ('angle_vel', 5.0, 10.0, None, False),
        ('strain', -0.02, 10.0, 2.0, False),
        ('strain', 0.0, 10.0, 2.0, False),
        ('strain', np.nan, 10.0, 2.0, False),
        ('strain', 0.05, 10.0, None, True),
        ('angle', 1.0, None, None, True),
        ('bogus', 1.0, 10.0, 2.0, True),
    ],
)
def test_convert_to_curvature_modes(mode, value, dc, xw, expect_raises):
    if expect_raises:
        with pytest.raises(ValueError):
            convert_to_curvature(value, mode, dclamp_mm=dc, xsec_width_mm=xw)
        return
    out = convert_to_curvature(value, mode, dclamp_mm=dc, xsec_width_mm=xw)
    arr = np.asarray(out, dtype=float).reshape(-1)
    assert arr.size >= 1
    assert np.all(np.isfinite(arr) | np.isnan(arr))


def test_angle_vel_matches_identity_for_deg_per_s():
    dc = 12.0
    v = 7.5
    kdot = float(convert_to_curvature(v, 'angle_vel', dclamp_mm=dc))
    back = float(np.rad2deg(kdot * (dc / 1000.0)))
    assert back == pytest.approx(v)
