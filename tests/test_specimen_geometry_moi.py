"""Ground-truth tests for the user-defined specimen geometry -> volume -> center-axis I_spec model.

Reference (uniform solid cylinder, radius r, length L, density rho, about a CENTER
transverse axis):
    volume   = pi * r**2 * L
    mass     = volume * rho
    I_center = m * (3*r**2 + L**2) / 12
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from bender_functions import Bender  # noqa: E402


def _fresh_bender():
    b = Bender.__new__(Bender)  # avoid hardware/config __init__ on a non-rig machine
    return b


def test_uniform_cylinder_volume_mass_and_center_axis_moi():
    r = 3.0          # mm radius -> height = depth = diameter = 2r
    L = 40.0         # mm length
    rho = 1.06e-3    # g/mm^3 (muscle-like)
    diameter = 2.0 * r

    b = _fresh_bender()
    # Two stations spanning -L/2 .. +L/2 about the AoR center; constant cross-section.
    b.set_specimen_geometry_inertial_model(
        heights_mm=[diameter, diameter],
        depths_mm=[diameter, diameter],
        positions_mm=[-L / 2.0, L / 2.0],
        density_g_per_mm3=rho,
    )

    vol_truth = np.pi * r**2 * L
    mass_truth = vol_truth * rho
    moi_truth = mass_truth * (3.0 * r**2 + L**2) / 12.0

    assert b.specimen_volume_mm3 == pytest.approx(vol_truth, rel=1e-9)
    assert b.specimen_mass_specimen == pytest.approx(mass_truth, rel=1e-9)
    assert b.specimen_moi_specimen == pytest.approx(moi_truth, rel=1e-9)


def test_multi_station_matches_uniform_cylinder():
    # Extra interior stations with the same dimensions must not change the result.
    r = 2.5
    L = 30.0
    rho = 1.0e-3
    d = 2.0 * r
    b = _fresh_bender()
    b.set_specimen_geometry_inertial_model(
        heights_mm=[d, d, d, d],
        depths_mm=[d, d, d, d],
        positions_mm=[-L / 2.0, -L / 6.0, L / 6.0, L / 2.0],
        density_g_per_mm3=rho,
    )
    vol_truth = np.pi * r**2 * L
    moi_truth = (vol_truth * rho) * (3.0 * r**2 + L**2) / 12.0
    assert b.specimen_volume_mm3 == pytest.approx(vol_truth, rel=1e-9)
    assert b.specimen_moi_specimen == pytest.approx(moi_truth, rel=1e-9)


def test_unsorted_positions_are_handled():
    r = 2.0
    L = 20.0
    rho = 1.0e-3
    d = 2.0 * r
    b = _fresh_bender()
    b.set_specimen_geometry_inertial_model(
        heights_mm=[d, d],
        depths_mm=[d, d],
        positions_mm=[L / 2.0, -L / 2.0],  # reversed
        density_g_per_mm3=rho,
    )
    vol_truth = np.pi * r**2 * L
    assert b.specimen_volume_mm3 == pytest.approx(vol_truth, rel=1e-9)


def test_mismatched_list_lengths_raise():
    b = _fresh_bender()
    with pytest.raises(ValueError):
        b.set_specimen_geometry_inertial_model(
            heights_mm=[2.0, 3.0, 4.0],
            depths_mm=[2.0, 3.0],
            positions_mm=[0.0, 1.0],
            density_g_per_mm3=1.0e-3,
        )


@pytest.mark.parametrize(
    "heights,depths,positions",
    [
        ([2.0, -1.0], [2.0, 2.0], [0.0, 1.0]),      # negative dimension
        ([2.0, 0.0], [2.0, 2.0], [0.0, 1.0]),        # zero dimension
        ([2.0, np.nan], [2.0, 2.0], [0.0, 1.0]),     # NaN dimension
        ([2.0, 2.0], [2.0, 2.0], [0.0, np.nan]),     # NaN position
    ],
)
def test_invalid_numeric_inputs_raise(heights, depths, positions):
    b = _fresh_bender()
    with pytest.raises(ValueError):
        b.set_specimen_geometry_inertial_model(
            heights_mm=heights,
            depths_mm=depths,
            positions_mm=positions,
            density_g_per_mm3=1.0e-3,
        )


def test_single_station_rejected():
    b = _fresh_bender()
    with pytest.raises(ValueError):
        b.set_specimen_geometry_inertial_model(
            heights_mm=[2.0],
            depths_mm=[2.0],
            positions_mm=[0.0],
            density_g_per_mm3=1.0e-3,
        )
