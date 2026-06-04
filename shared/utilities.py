"""
Protocol-wide helpers (motor angle demand, etc.).

``compute_max_rotation_deg`` reuses the same timeline / conversion paths as
``bender_gui_preview.build_protocol_preview`` and ``Bender`` step protocols.
"""
from __future__ import annotations

from typing import Any, Dict, Optional, Union

import numpy as np

Number = Union[int, float]


def format_max_rotation_message(max_rotation_deg: Number) -> str:
    """Single informational line shown in the GUI and printed at run start."""
    return f"This run will move a maximum of {float(max_rotation_deg):.1f}°"


def _effective_test_type(protocol_params: Dict[str, Any]) -> str:
    tt = protocol_params.get('test_type')
    if tt is not None and str(tt).strip():
        return str(tt).strip()
    b = protocol_params.get('bender')
    if b is not None:
        return str(getattr(b, 'test_type', None) or 'dynamic')
    raise ValueError("protocol_params must include 'test_type' or a 'bender' with test_type set.")


def _max_rotation_from_preview(b: Any, test_type: str) -> float:
    from bender_gui_preview import build_protocol_preview

    prev = build_protocol_preview(b, requested_test_type=test_type)
    if prev.get('error'):
        raise ValueError(str(prev['error']))
    ang = prev.get('angle')
    if ang is None:
        raise ValueError(f"No commanded angle timeline for test_type={test_type!r}.")
    a = np.asarray(ang, dtype=float).reshape(-1)
    if a.size == 0 or not np.all(np.isfinite(a)):
        raise ValueError("Commanded angle timeline is empty or non-finite.")
    return float(np.max(np.abs(a)))


def _max_rotation_isometric_fallback(b: Any) -> float:
    """When preview refuses num_steps < 2, use the same κ→deg conversion as :meth:`Bender.isometric`."""
    from bender_functions import convert_to_curvature

    initial = getattr(b, 'isometric_initial', None)
    final = getattr(b, 'isometric_final', None)
    num_steps = getattr(b, 'isometric_num_steps', None)
    if initial is None or final is None or num_steps is None:
        raise ValueError("isometric max rotation needs isometric_initial, isometric_final, isometric_num_steps.")
    mode = getattr(b, 'isometric_mode', None) or 'strain'
    dc = b._effective_dclamp_mm() if hasattr(b, '_effective_dclamp_mm') else getattr(b, 'dclamp', None)
    if dc is None:
        raise ValueError("isometric max rotation needs clamp spacing (dclamp / test_segment_length_mm).")
    xw = getattr(b, 'xsec_width', None)
    md = float(getattr(b, 'target_muscle_depth_mm', 0.0) or 0.0)
    vals = np.linspace(float(initial), float(final), int(num_steps))
    kappa = convert_to_curvature(vals, mode, dclamp_mm=float(dc), xsec_width_mm=xw, target_muscle_depth_mm=md)
    targets_deg_raw = np.rad2deg(np.asarray(kappa, dtype=float) * (float(dc) / 1000.0))
    peak = float(np.max(np.abs(targets_deg_raw)))
    sp = getattr(b, 'isometric_stim_params', None) or {}
    ml_n = sp.get('isometric_mirror_target_left', getattr(b, 'isometric_mirror_target_left', None))
    mr_n = sp.get('isometric_mirror_target_right', getattr(b, 'isometric_mirror_target_right', None))
    mirror_bm = bool(sp.get('bilateral_mirror_motor', getattr(b, 'bilateral_mirror_motor', False)))
    if mirror_bm and ml_n is not None and mr_n is not None:
        kL = convert_to_curvature(
            np.asarray([float(ml_n)]), mode, dclamp_mm=float(dc), xsec_width_mm=xw, target_muscle_depth_mm=md
        )
        kR = convert_to_curvature(
            np.asarray([float(mr_n)]), mode, dclamp_mm=float(dc), xsec_width_mm=xw, target_muscle_depth_mm=md
        )
        peak = max(
            peak,
            float(np.abs(np.rad2deg(float(np.asarray(kL).reshape(-1)[0]) * (float(dc) / 1000.0)))),
            float(np.abs(np.rad2deg(float(np.asarray(kR).reshape(-1)[0]) * (float(dc) / 1000.0)))),
        )
    rec = getattr(b, 'recruitment', 'bilateral_simultaneous')
    if hasattr(b, '_normalize_recruitment'):
        rec = b._normalize_recruitment(rec)
    if hasattr(b, 'recruitment_unilateral_lateral_index') and hasattr(b, 'motor_command_sign_for_bend_toward_index'):
        uidx = b.recruitment_unilateral_lateral_index(rec)
        if uidx is not None:
            peak *= abs(float(b.motor_command_sign_for_bend_toward_index(uidx)))
    return peak


def compute_max_rotation_deg(protocol_params: Dict[str, Any]) -> float:
    """
    Maximum absolute commanded motor angle (degrees) for a protocol.

    Parameters
    ----------
    protocol_params : dict
        Must include ``bender`` (a configured :class:`~bender_functions.Bender`).
        Optional ``test_type`` overrides ``bender.test_type``.

    Returns
    -------
    float
        Peak |commanded angle| over the full protocol timeline (preview path), or
        peak hold target for isometric when preview requires ≥2 steps.
    """
    b = protocol_params.get('bender')
    if b is None:
        raise ValueError("protocol_params must include 'bender'.")
    tt = _effective_test_type(protocol_params)
    motion_tt = tt
    if tt == 'calibration':
        motion_tt = str(getattr(b, 'calibration_base_test_type', 'dynamic') or 'dynamic')

    if tt == 'isometric':
        n = getattr(b, 'isometric_num_steps', None)
        try:
            if n is not None and int(n) < 2:
                return _max_rotation_isometric_fallback(b)
        except (TypeError, ValueError):
            pass

    try:
        return _max_rotation_from_preview(b, tt)
    except ValueError as preview_err:
        if tt == 'isometric':
            return _max_rotation_isometric_fallback(b)
        raise preview_err
