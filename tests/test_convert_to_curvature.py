"""Stress tests for :func:`bender_functions.convert_to_curvature`."""
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

from bender_functions import (
    Bender,
    _strain_lever_arm_m,
    convert_to_curvature,
    curvature_to_strain_fraction,
)


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


def test_muscle_depth_zero_matches_surface_half_width():
    xw = 10.0
    eps = 0.05
    k_surface = float(convert_to_curvature(eps, 'strain', xsec_width_mm=xw, target_muscle_depth_mm=0.0))
    k_default = float(convert_to_curvature(eps, 'strain', xsec_width_mm=xw))
    assert k_surface == pytest.approx(k_default)


def test_muscle_depth_increases_curvature_for_same_strain():
    eps = 0.05
    xw = 10.0
    k0 = float(convert_to_curvature(eps, 'strain', xsec_width_mm=xw, target_muscle_depth_mm=0.0))
    k3 = float(convert_to_curvature(eps, 'strain', xsec_width_mm=xw, target_muscle_depth_mm=3.0))
    assert k3 > k0


@pytest.mark.parametrize('md', [-0.1, 5.0, 5.0001])
def test_muscle_depth_rejected_at_or_beyond_half_width(md):
    with pytest.raises(ValueError, match='target_muscle_depth_mm'):
        _strain_lever_arm_m(10.0, md)


def test_strain_lever_arm_round_trip_three_sites():
    """Forward convert_to_curvature; inverse via preview path and organize_cycles formula."""
    eps0, xw, md = 0.05, 10.0, 3.0
    kappa = float(convert_to_curvature(eps0, 'strain', xsec_width_mm=xw, target_muscle_depth_mm=md))
    lever_m = _strain_lever_arm_m(xw, md)
    strain_preview = float(np.asarray(kappa, dtype=float) * lever_m)
    strain_org = float(curvature_to_strain_fraction(kappa, xsec_width_mm=xw, target_muscle_depth_mm=md))

    b = Bender('jimenez_bender_config_A')
    b.organize_cycles(
        [kappa],
        [1.0],
        False,
        1,
        0,
        52.0,
        xw,
        0,
        [0.0],
        [0.0],
        0.0,
        target_muscle_depth_mm=md,
    )
    strain_from_organize = float(np.asarray(b.organized_strains, dtype=float).reshape(-1)[0])

    assert strain_preview == pytest.approx(eps0, rel=1e-9, abs=1e-12)
    assert strain_org == pytest.approx(eps0, rel=1e-9, abs=1e-12)
    assert strain_from_organize == pytest.approx(eps0, rel=1e-9, abs=1e-12)
    assert strain_preview == pytest.approx(strain_org)
    assert strain_preview == pytest.approx(strain_from_organize)
