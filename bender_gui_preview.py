"""
Protocol preview for the Streamlit GUI: motor command timeline without DAQ.

Mirrors motion-generation branches in :meth:`bender_functions.Bender.run_experiment`
(isometric / isovelocity / calibration base / dynamic / sweep / step_change).
"""
from __future__ import annotations

import math
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

PreviewResult = Dict[str, Any]


def _strain_mode_caption(mode: Any) -> str:
    """Short label for preview tables: decimal ε vs percent, etc."""
    m = str(mode or 'strain')
    return {
        'strain': 'decimal ε (0.05 ≙ 5%)',
        'strain_pct': 'percent / strain_pct (5 ≙ 5%)',
        'curvature': 'curvature κ (1/m)',
        'angle': 'motor angle (deg)',
    }.get(m, m)


def _as_float_list(x: Any) -> Optional[List[float]]:
    if x is None:
        return None
    try:
        arr = np.asarray(x, dtype=float).reshape(-1)
        if arr.size == 0:
            return None
        return [float(v) for v in arr.tolist()]
    except (TypeError, ValueError):
        return None


def _downsample(t, y, max_points: int):
    t = np.asarray(t, dtype=float).reshape(-1)
    y = np.asarray(y, dtype=float).reshape(-1)
    n = min(t.size, y.size)
    if n == 0:
        return t[:0], y[:0]
    t = t[:n]
    y = y[:n]
    if n <= max_points:
        return t, y
    idx = np.linspace(0, n - 1, num=max_points, dtype=int)
    return t[idx], y[idx]


def _motor_angle_to_kappa_strain(
    angle_deg: np.ndarray, *, dclamp_mm: float, xsec_width_mm: Optional[float]
) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """Inverse of motor command: κ (1/m) and surface engineering strain ε from angle (deg)."""
    ang = np.asarray(angle_deg, dtype=float).reshape(-1)
    kappa = np.deg2rad(ang) * 1000.0 / float(dclamp_mm)
    if xsec_width_mm is None or not math.isfinite(float(xsec_width_mm)) or float(xsec_width_mm) <= 0:
        return kappa, None
    half_m = (float(xsec_width_mm) / 2.0) / 1000.0
    return kappa, kappa * half_m


def _attach_native_unit_series(out: PreviewResult, b: Any, t: np.ndarray, angle: np.ndarray) -> None:
    dc = getattr(b, '_effective_dclamp_mm', lambda: None)()
    if dc is None:
        dc = getattr(b, 'dclamp', None)
    xw = getattr(b, 'xsec_width', None)
    if dc is None:
        return
    kappa, strain = _motor_angle_to_kappa_strain(angle, dclamp_mm=float(dc), xsec_width_mm=xw)
    out['curvature'] = kappa
    if strain is not None:
        out['strain'] = strain


def _stim_table_rows(b: Any, sp: dict, *, recruitment: str) -> List[dict]:
    spr_raw = sp.get('stim_pulse_rate', None)
    spr = float(spr_raw) if spr_raw is not None else float(getattr(b, 'stim_pulse_rate', 75.0) or 75.0)
    return [
        {'metric': 'stimulation enabled', 'value': bool(sp.get('is_stim', False))},
        {'metric': 'pulse rate (Hz)', 'value': spr},
        {'metric': 'stim voltage (V)', 'value': float(sp.get('stim_voltage', 5.0))},
        {'metric': 'recruitment', 'value': recruitment},
    ]


def _stim_for_ramp_hold(
    b: Any,
    t: np.ndarray,
    ramp_s: float,
    hold_s: float,
    sp: dict,
    rec: str,
    *,
    mirror: bool = False,
    hold_windows: Optional[Tuple[Tuple[float, float], Tuple[float, float]]] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    spr_raw = sp.get('stim_pulse_rate', None)
    spr = float(spr_raw) if spr_raw is not None else float(getattr(b, 'stim_pulse_rate', 75.0) or 75.0)
    stim_voltage = float(sp.get('stim_voltage', 5.0))
    settle = float(sp.get('settle_before_stim_s', 0.5))
    stim_duration_s = sp.get('stim_duration_s', None)
    is_stim = bool(sp.get('is_stim', False))
    seq_frac = float(sp.get('bilateral_sequential_left_frac', getattr(b, 'bilateral_sequential_left_frac', 0.5)))
    t = np.asarray(t, dtype=float).reshape(-1)
    if mirror and hold_windows is not None:
        h1, h2 = hold_windows
        active_l = (t >= h1[0] + settle) & (t < h1[1])
        active_r = (t >= h2[0] + settle) & (t < h2[1])
        if stim_duration_s is not None:
            active_l &= t < (h1[0] + settle + float(stim_duration_s))
            active_r &= t < (h2[0] + settle + float(stim_duration_s))
        if is_stim and (np.any(active_l) or np.any(active_r)):
            p_l = b._pulse_carrier_volts(t, active_l, spr, stim_voltage)
            p_r = b._pulse_carrier_volts(t, active_r, spr, stim_voltage)
            s1 = np.zeros_like(t)
            s2 = np.zeros_like(t)
            b._deposit_stim_on_side(p_l, 'left', s1, s2)
            b._deposit_stim_on_side(p_r, 'right', s1, s2)
            return s1, s2
        return np.zeros_like(t), np.zeros_like(t)
    t_stim0 = float(ramp_s) + settle
    if stim_duration_s is None:
        t_stim1 = float(ramp_s) + float(hold_s)
    else:
        t_stim1 = t_stim0 + float(stim_duration_s)
    t_stim1 = min(t_stim1, float(t[-1]) + 1e-9)
    active = (t >= t_stim0) & (t < t_stim1)
    if is_stim and np.any(active):
        pulse = b._pulse_carrier_volts(t, active, spr, stim_voltage)
        return b._route_recruitment_stim(pulse, rec, sequential_left_frac=seq_frac)
    return np.zeros_like(t), np.zeros_like(t)


def _default_stim_duties_phases(b: Any) -> Tuple[List[float], List[float]]:
    d = getattr(b, 'all_stimduties', None)
    p = getattr(b, 'all_stimphases', None)
    if d is None or (isinstance(d, (list, tuple, np.ndarray)) and len(d) == 0):
        d = [0.3]
    if p is None or (isinstance(p, (list, tuple, np.ndarray)) and len(p) == 0):
        p = [0.5]
    d = _as_float_list(d) or [0.3]
    p = _as_float_list(p) or [0.5]
    return d, p


def build_protocol_preview(
    b: Any,
    *,
    requested_test_type: str,
    max_plot_points: int = 8000,
) -> PreviewResult:
    """
    Return a dict with keys:
      - ok: bool
      - error: optional str
      - requested_test_type, motion_test_type: str
      - table: list of dict rows for st.dataframe
      - t, angle, anglevel: optional 1d arrays (motor preview)
      - t_plot, angle_plot, anglevel_plot: downsampled for plotting
    """
    req = str(requested_test_type)
    motion_tt = req
    if req == 'calibration':
        motion_tt = str(getattr(b, 'calibration_base_test_type', 'dynamic') or 'dynamic')

    out: PreviewResult = {
        'ok': False,
        'error': None,
        'requested_test_type': req,
        'motion_test_type': motion_tt,
        'table': [],
        't': None,
        'angle': None,
        'anglevel': None,
        'stim_s1': None,
        'stim_s2': None,
        'stim_total': None,
        'strain': None,
        'curvature': None,
        't_plot': None,
        'angle_plot': None,
        'anglevel_plot': None,
        'stim_plot': None,
        'strain_plot': None,
        'curvature_plot': None,
    }

    try:
        if req in ('isometric', 'isovelocity'):
            out.update(_preview_step_protocols(b, req))
        elif motion_tt == 'dynamic':
            out.update(_preview_dynamic(b, max_plot_points))
        elif motion_tt == 'frequency_sweep':
            out.update(_preview_frequency_sweep(b, max_plot_points))
        elif motion_tt == 'frequency_step':
            out.update(_preview_frequency_step(b, max_plot_points))
        elif motion_tt == 'curvature_step':
            out.update(_preview_curvature_step(b, max_plot_points))
        elif motion_tt == 'step_change':
            out.update(_preview_step_change(b, max_plot_points))
        else:
            out['error'] = f'No preview implemented for {req!r}.'
            return out

        if out.get('error'):
            return out
        out['ok'] = True
        t = out.get('t')
        ang = out.get('angle')
        av = out.get('anglevel')
        if t is not None and ang is not None:
            tp, ap = _downsample(t, ang, max_plot_points)
            out['t_plot'] = tp
            out['angle_plot'] = ap
            if av is not None:
                _, vp = _downsample(t, av, max_plot_points)
                out['anglevel_plot'] = vp
            else:
                out['anglevel_plot'] = None
            _attach_native_unit_series(out, b, np.asarray(t), np.asarray(ang))
            if out.get('strain') is not None:
                _, out['strain_plot'] = _downsample(t, out['strain'], max_plot_points)
            if out.get('curvature') is not None:
                _, out['curvature_plot'] = _downsample(t, out['curvature'], max_plot_points)
            s1 = out.get('stim_s1')
            s2 = out.get('stim_s2')
            if s1 is not None:
                s1a = np.asarray(s1, dtype=float).reshape(-1)
                s2a = np.zeros_like(s1a) if s2 is None else np.asarray(s2, dtype=float).reshape(-1)
                n_st = min(s1a.size, s2a.size, np.asarray(t).size)
                stot = s1a[:n_st] + s2a[:n_st]
                out['stim_total'] = stot
                _, out['stim_plot'] = _downsample(t, stot, max_plot_points)
        return out
    except Exception as e:
        out['error'] = f'{type(e).__name__}: {e}'
        return out


def _isovelocity_stim_params_from_b(b: Any) -> dict:
    """Defaults + ``isovelocity_stim_params`` merged like :meth:`bender_functions.Bender.isovelocity`."""
    sp = {
        'iso_duration_s': float(getattr(b, 'isovelocity_iso_duration_s', 0.2) or 0.2),
        'pre_hold_s': float(getattr(b, 'isovelocity_pre_hold_s', 0.3) or 0.3),
        'settle_before_stim_s': 0.02,
        'pre_iso_stim_duration_s': 0.0,
        'stim_duration_s': None,
        'is_stim': False,
        'stim_pulse_rate': None,
        'stim_voltage': 5.0,
        'device_name': None,
    }
    user = getattr(b, 'isovelocity_stim_params', None)
    if isinstance(user, dict):
        sp.update(user)
    return b._stim_params_with_lateral(sp)


def _preview_concat_isovelocity_timeline(
    b: Any, theta0_fixed: float, velocities_deg_per_s: np.ndarray, sp: dict
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Concatenate per-step isovelocity command timelines (same geometry as
    :meth:`bender_functions.Bender._run_isovelocity_steps`, without DAQ).
    Each step starts again from ``theta0_fixed``.
    """
    vels = np.atleast_1d(np.asarray(velocities_deg_per_s, dtype=float)).reshape(-1)
    n = int(vels.size)
    if n < 1:
        raise ValueError('isovelocity preview needs at least one velocity step.')
    spr_raw = sp.get('stim_pulse_rate', None)
    spr = float(spr_raw) if spr_raw is not None else float(getattr(b, 'stim_pulse_rate', 75.0) or 75.0)
    daq_hz = float(getattr(b, 'daq_ai_sample_rate_hz', 0.0) or 0.0)
    if not np.isfinite(daq_hz) or daq_hz <= 0:
        daq_hz = 1000.0
    pre_hold_s = float(sp['pre_hold_s'])
    iso_duration_s = float(sp['iso_duration_s'])
    settle_before_stim_s = float(sp.get('settle_before_stim_s', 0.02))
    pre_iso_stim_duration_s = float(sp.get('pre_iso_stim_duration_s', 0.0))
    stim_duration_s = sp.get('stim_duration_s', None)
    is_stim = bool(sp.get('is_stim', False))
    stim_voltage = float(sp.get('stim_voltage', 5.0))
    rec = b._normalize_recruitment(
        sp.get('recruitment', sp.get('lateral_mode', getattr(b, 'recruitment', 'bilateral_simultaneous')))
    )
    bm = bool(sp.get('bilateral_mirror_motor', getattr(b, 'bilateral_mirror_motor', False)))
    rec = b._recruitment_with_bilateral_mirror_motor(rec, bm)
    mirror = bm and rec == 'bilateral_sequential'
    seq_frac = float(
        sp.get('bilateral_sequential_left_frac', getattr(b, 'bilateral_sequential_left_frac', 0.5))
    )
    th0_fixed = float(theta0_fixed)
    t_chunks: List[np.ndarray] = []
    a_chunks: List[np.ndarray] = []
    w_chunks: List[np.ndarray] = []
    s1_chunks: List[np.ndarray] = []
    s2_chunks: List[np.ndarray] = []
    toff = 0.0
    for i in range(n):
        v_mag = abs(float(vels[i]))
        if mirror:
            v1 = v_mag * b.motor_command_sign_for_bend_toward_index(b.specimen_side_index_left)
            v2 = v_mag * b.motor_command_sign_for_bend_toward_index(b.specimen_side_index_right)
            d1 = b._isovelocity_one_block(
                th0_fixed,
                v1,
                pre_hold_s=pre_hold_s,
                iso_duration_s=iso_duration_s,
                settle_before_stim_s=settle_before_stim_s,
                pre_iso_stim_duration_s=pre_iso_stim_duration_s,
                stim_duration_s=stim_duration_s,
                is_stim=is_stim,
                spr=spr,
                stim_voltage=stim_voltage,
                daq_hz=daq_hz,
                recruitment=rec,
                sequential_left_frac=seq_frac,
                mirror_stim_side='left',
            )
            th_mid = float(d1['angle'][-1])
            d2 = b._isovelocity_one_block(
                th_mid,
                v2,
                pre_hold_s=pre_hold_s,
                iso_duration_s=iso_duration_s,
                settle_before_stim_s=settle_before_stim_s,
                pre_iso_stim_duration_s=pre_iso_stim_duration_s,
                stim_duration_s=stim_duration_s,
                is_stim=is_stim,
                spr=spr,
                stim_voltage=stim_voltage,
                daq_hz=daq_hz,
                recruitment=rec,
                sequential_left_frac=seq_frac,
                mirror_stim_side='right',
            )
            off = float(d1['t'][-1])
            t_seg = np.concatenate([d1['t'], d2['t'][1:] + off])
            a_seg = np.concatenate([d1['angle'], d2['angle'][1:]])
            w_seg = np.concatenate([d1['anglevel'], d2['anglevel'][1:]])
            s1_seg = np.concatenate([d1['s1'], d2['s1'][1:]])
            s2_seg = np.concatenate([d1['s2'], d2['s2'][1:]])
        else:
            v_sign = float(vels[i])
            uidx = b.recruitment_unilateral_lateral_index(rec)
            if uidx is not None:
                v_sign = v_mag * b.motor_command_sign_for_bend_toward_index(uidx)
            d0 = b._isovelocity_one_block(
                th0_fixed,
                v_sign,
                pre_hold_s=pre_hold_s,
                iso_duration_s=iso_duration_s,
                settle_before_stim_s=settle_before_stim_s,
                pre_iso_stim_duration_s=pre_iso_stim_duration_s,
                stim_duration_s=stim_duration_s,
                is_stim=is_stim,
                spr=spr,
                stim_voltage=stim_voltage,
                daq_hz=daq_hz,
                recruitment=rec,
                sequential_left_frac=seq_frac,
                mirror_stim_side=None,
            )
            t_seg, a_seg, w_seg = d0['t'], d0['angle'], d0['anglevel']
            s1_seg, s2_seg = d0['s1'], d0['s2']
        t_chunks.append(np.asarray(t_seg, dtype=float) + toff)
        a_chunks.append(np.asarray(a_seg, dtype=float))
        w_chunks.append(np.asarray(w_seg, dtype=float))
        s1_chunks.append(np.asarray(s1_seg, dtype=float))
        s2_chunks.append(np.asarray(s2_seg, dtype=float))
        toff = float(t_chunks[-1][-1])
    t = np.concatenate(t_chunks)
    angle = np.concatenate(a_chunks)
    anglevel = np.concatenate(w_chunks)
    s1 = np.concatenate(s1_chunks)
    s2 = np.concatenate(s2_chunks)
    return t, angle, anglevel, s1, s2


def _isometric_stim_params_from_b(b: Any) -> dict:
    sp = {
        'ramp_duration_s': 2.0,
        'hold_duration_s': 5.0,
        'settle_before_stim_s': 0.5,
        'stim_duration_s': None,
        'inter_step_interval_s': None,
        'is_stim': False,
        'stim_pulse_rate': None,
        'stim_voltage': 5.0,
        'device_name': None,
    }
    user = getattr(b, 'isometric_stim_params', None)
    if isinstance(user, dict):
        sp.update(user)
    if sp.get('inter_step_interval_s', None) is None:
        sp['inter_step_interval_s'] = float(getattr(b, 'isometric_inter_step_interval_s', 0.0) or 0.0)
    return sp


def _preview_concat_isometric_timeline(
    b: Any, targets_deg: np.ndarray, sp: dict, *, mode: str = 'strain'
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Match :meth:`Bender._run_force_length_steps` ramp/hold segments (and optional mirror path)."""
    targets_deg = np.atleast_1d(np.asarray(targets_deg, dtype=float)).reshape(-1)
    num_steps = int(targets_deg.size)
    daq_hz = float(getattr(b, 'daq_ai_sample_rate_hz', 0.0) or 0.0)
    if not np.isfinite(daq_hz) or daq_hz <= 0:
        daq_hz = 1000.0
    ramp = float(sp.get('ramp_duration_s', 2.0))
    hold = float(sp.get('hold_duration_s', 5.0))
    gap_s = float(sp.get('inter_step_interval_s', 0.0) or 0.0)
    rec = b._normalize_recruitment(
        sp.get('recruitment', sp.get('lateral_mode', getattr(b, 'recruitment', 'bilateral_simultaneous')))
    )
    bm = bool(sp.get('bilateral_mirror_motor', getattr(b, 'bilateral_mirror_motor', False)))
    rec = b._recruitment_with_bilateral_mirror_motor(rec, bm)
    mirror = bm and rec == 'bilateral_sequential'
    dc = getattr(b, '_effective_dclamp_mm', lambda: None)()
    xw = getattr(b, 'xsec_width', None)
    mhl = mhr = None
    if mirror and dc is not None:
        ml_n = sp.get('isometric_mirror_target_left', getattr(b, 'isometric_mirror_target_left', None))
        mr_n = sp.get('isometric_mirror_target_right', getattr(b, 'isometric_mirror_target_right', None))
        if ml_n is not None and mr_n is not None:
            try:
                from bender_functions import convert_to_curvature

                fL, fR = float(ml_n), float(mr_n)
                if np.isfinite(fL) and np.isfinite(fR):
                    kL = convert_to_curvature(np.asarray([fL]), mode, dclamp_mm=float(dc), xsec_width_mm=xw)
                    kR = convert_to_curvature(np.asarray([fR]), mode, dclamp_mm=float(dc), xsec_width_mm=xw)
                    mhl = float(
                        np.abs(np.rad2deg(float(np.asarray(kL, dtype=float).reshape(-1)[0]) * (float(dc) / 1000.0)))
                    )
                    mhr = float(
                        np.abs(np.rad2deg(float(np.asarray(kR, dtype=float).reshape(-1)[0]) * (float(dc) / 1000.0)))
                    )
            except (TypeError, ValueError):
                mhl = mhr = None

    t_chunks: List[np.ndarray] = []
    a_chunks: List[np.ndarray] = []
    s1_chunks: List[np.ndarray] = []
    s2_chunks: List[np.ndarray] = []
    toff = 0.0

    for i in range(num_steps):
        if i > 0 and gap_s > 0:
            last_ang = float(a_chunks[-1][-1]) if a_chunks else float(targets_deg[0])
            dt = 1.0 / daq_hz
            n_g = max(2, int(round(gap_s / dt)) + 1)
            tg = np.linspace(0.0, gap_s, n_g, dtype=float)
            ag = np.full(n_g, last_ang, dtype=float)
            zg = np.zeros(n_g, dtype=float)
            t_chunks.append(tg + toff)
            a_chunks.append(ag)
            s1_chunks.append(zg)
            s2_chunks.append(zg)
            toff = float(t_chunks[-1][-1])

        target = float(targets_deg[i])
        prev = float(targets_deg[i - 1]) if i > 0 else float(targets_deg[0])
        if mirror:
            kw = {}
            if mhl is not None:
                kw['mirror_abs_deg_left'] = mhl
            if mhr is not None:
                kw['mirror_abs_deg_right'] = mhr
            tloc, aloc, _wloc, h1, h2 = b._timeline_mirror_two_holds(prev, target, ramp, hold, daq_hz, **kw)
            s1loc, s2loc = _stim_for_ramp_hold(b, tloc, ramp, hold, sp, rec, mirror=True, hold_windows=(h1, h2))
        else:
            tloc, aloc, _wloc = b._timeline_ramp_hold(prev, target, ramp, hold, daq_hz)
            s1loc, s2loc = _stim_for_ramp_hold(b, tloc, ramp, hold, sp, rec)
        t_chunks.append(tloc + toff)
        a_chunks.append(aloc)
        s1_chunks.append(s1loc)
        s2_chunks.append(s2loc)
        toff = float(t_chunks[-1][-1])

    t = np.concatenate(t_chunks)
    angle = np.concatenate(a_chunks)
    anglevel = np.gradient(angle, t, edge_order=1) if t.size >= 2 else np.zeros_like(angle)
    s1 = np.concatenate(s1_chunks)
    s2 = np.concatenate(s2_chunks)
    return t, angle, anglevel, s1, s2


def _preview_step_protocols(b: Any, req: str) -> PreviewResult:
    r: PreviewResult = {
        'table': [],
        't': None,
        'angle': None,
        'anglevel': None,
        'error': None,
    }
    if req == 'isometric':
        initial = getattr(b, 'isometric_initial', None)
        final = getattr(b, 'isometric_final', None)
        n = getattr(b, 'isometric_num_steps', None)
        mode = getattr(b, 'isometric_mode', None) or 'strain'
        if initial is None or final is None or n is None:
            r['error'] = 'Isometric preview needs isometric_initial, isometric_final, isometric_num_steps.'
            return r
        n = int(n)
        if n < 2:
            r['error'] = 'isometric_num_steps should be at least 2 for a useful preview.'
            return r
        vals = np.linspace(float(initial), float(final), n)
        seq_idx = np.arange(vals.size, dtype=int)
        if bool(getattr(b, 'isometric_randomize', False)) and vals.size > 1:
            rs = getattr(b, 'isometric_random_seed', None)
            rng = np.random.default_rng(None if rs is None else int(rs))
            rng.shuffle(seq_idx)
            vals = vals[seq_idx]
        dc = getattr(b, '_effective_dclamp_mm', lambda: None)()
        if dc is None:
            r['error'] = (
                'Isometric preview needs clamp spacing (mm): set **test_segment_length_mm** (or `dclamp`) '
                'in biometrics / section 2.'
            )
            return r
        xw = getattr(b, 'xsec_width', None)
        try:
            from bender_functions import convert_to_curvature
        except ImportError:
            r['error'] = 'Could not import convert_to_curvature for isometric preview.'
            return r
        kappa = convert_to_curvature(vals, mode, dclamp_mm=float(dc), xsec_width_mm=xw)
        targets_deg = np.rad2deg(kappa * (float(dc) / 1000.0))
        sp = _isometric_stim_params_from_b(b)
        for k in ('isometric_mirror_target_left', 'isometric_mirror_target_right'):
            if k not in sp and getattr(b, k, None) is not None:
                sp[k] = getattr(b, k)
        rec = b._normalize_recruitment(
            sp.get('recruitment', sp.get('lateral_mode', getattr(b, 'recruitment', 'bilateral_simultaneous')))
        )
        mirror_bm = bool(sp.get('bilateral_mirror_motor', getattr(b, 'bilateral_mirror_motor', False)))
        rec = b._recruitment_with_bilateral_mirror_motor(rec, mirror_bm)
        uidx = b.recruitment_unilateral_lateral_index(rec)
        if uidx is not None:
            targets_deg = targets_deg * b.motor_command_sign_for_bend_toward_index(uidx)

        for i, v in enumerate(vals):
            r['table'].append(
                {
                    'step_index': i,
                    'sequence_index': int(seq_idx[i]),
                    'target_value': float(v),
                    'target_deg': float(targets_deg[i]),
                    'curvature_1_per_m': float(kappa[i]),
                    'input_mode': _strain_mode_caption(mode),
                }
            )
        try:
            t, angle, anglevel, s1, s2 = _preview_concat_isometric_timeline(b, targets_deg, sp, mode=str(mode))
        except Exception as e:
            r['error'] = f'Isometric timeline preview failed: {type(e).__name__}: {e}'
            return r
        r['t'] = t
        r['angle'] = angle
        r['anglevel'] = anglevel
        r['stim_s1'] = s1
        r['stim_s2'] = s2
        r['table'].extend(_stim_table_rows(b, sp, recruitment=str(rec)))
        r['preview_isometric'] = True
        return r

    min_v = getattr(b, 'isovelocity_min_vel', None)
    max_v = getattr(b, 'isovelocity_max_vel', None)
    ss = getattr(b, 'isovelocity_starting_strain', None)
    n = getattr(b, 'isovelocity_num_steps', None)
    mode = getattr(b, 'isovelocity_starting_strain_mode', None) or 'strain'
    if min_v is None or max_v is None or ss is None or n is None:
        r['error'] = (
            'Isovelocity preview needs isovelocity_min_vel, isovelocity_max_vel, '
            'isovelocity_starting_strain, isovelocity_num_steps.'
        )
        return r
    n = int(n)
    vels = np.linspace(float(min_v), float(max_v), max(1, n))
    seq_idx = np.arange(vels.size, dtype=int)
    if bool(getattr(b, 'isovelocity_randomize', False)) and vels.size > 1:
        rs = getattr(b, 'isovelocity_random_seed', None)
        rng = np.random.default_rng(None if rs is None else int(rs))
        rng.shuffle(seq_idx)
        vels = vels[seq_idx]
    dc = getattr(b, '_effective_dclamp_mm', lambda: None)()
    if dc is None:
        r['error'] = (
            'Isovelocity preview needs clamp spacing (mm): set **test_segment_length_mm** (or `dclamp`) '
            'in biometrics / section 2.'
        )
        return r
    xw = getattr(b, 'xsec_width', None)
    try:
        from bender_functions import convert_to_curvature
    except ImportError:
        r['error'] = 'Could not import convert_to_curvature for isovelocity preview.'
        return r
    k0 = convert_to_curvature(float(ss), mode, dclamp_mm=float(dc), xsec_width_mm=xw)
    k0 = float(np.asarray(k0).reshape(-1)[0])
    theta0 = float(np.rad2deg(k0 * (float(dc) / 1000.0)))
    velocity_mode = str(getattr(b, 'isovelocity_velocity_mode', None) or 'angle_vel').lower()
    kdot = convert_to_curvature(vels, velocity_mode, dclamp_mm=float(dc), xsec_width_mm=xw)
    vels_deg = np.rad2deg(np.asarray(kdot, dtype=float) * (float(dc) / 1000.0))
    sp = _isovelocity_stim_params_from_b(b)
    rec_iso = b._normalize_recruitment(
        sp.get('recruitment', sp.get('lateral_mode', getattr(b, 'recruitment', 'bilateral_simultaneous')))
    )
    uidx0 = b.recruitment_unilateral_lateral_index(rec_iso)
    if uidx0 is not None:
        theta0 = theta0 * b.motor_command_sign_for_bend_toward_index(uidx0)

    r['table'].append(
        {'row': 'starting strain', 'value': float(ss), 'unit': _strain_mode_caption(mode)}
    )
    vel_caption = {
        'strain_rate': 'dε/dt (1/s)',
        'strain_pct_rate': 'd(% strain)/dt (%/s)',
        'curvature_rate': 'dκ/dt (1/m/s)',
        'angle_vel': 'deg/s',
    }.get(velocity_mode, velocity_mode)
    for i, v in enumerate(vels):
        r['table'].append(
            {
                'row': f'velocity step {i}',
                'value': float(v),
                'unit': vel_caption,
                'motor_deg_s': float(vels_deg[i]),
                'sequence_index': int(seq_idx[i]),
            }
        )
    try:
        t, angle, anglevel, s1, s2 = _preview_concat_isovelocity_timeline(b, theta0, vels_deg, sp)
    except Exception as e:
        r['error'] = f'Isovelocity timeline preview failed: {type(e).__name__}: {e}'
        return r
    r['t'] = t
    r['angle'] = angle
    r['anglevel'] = anglevel
    r['stim_s1'] = s1
    r['stim_s2'] = s2
    r['table'].extend(_stim_table_rows(b, sp, recruitment=str(rec_iso)))
    r['preview_isovelocity'] = True
    return r


def _preview_dynamic(b: Any, max_plot_points: int) -> PreviewResult:
    r: PreviewResult = {'table': [], 't': None, 'angle': None, 'anglevel': None, 'error': None}
    organize = getattr(b, '_organize_cycles_for_dynamic_run', None)
    if not callable(organize):
        r['error'] = 'Dynamic preview requires a Bender with _organize_cycles_for_dynamic_run.'
        return r
    try:
        organize()
    except ValueError as e:
        r['error'] = str(e)
        return r
    spr = float(getattr(b, 'stim_pulse_rate', 0.0) or 0.0)
    angle, anglevel, tnorm, t = b.make_cycles_dynamic(
        b.period_by_cycle, b.freq_by_cycle, b.amp_by_cycle, record_protocol=False
    )
    r['t'] = np.asarray(t, dtype=float).reshape(-1)
    r['angle'] = np.asarray(angle, dtype=float).reshape(-1)
    r['anglevel'] = np.asarray(anglevel, dtype=float).reshape(-1)
    r['table'] = [
        {'metric': 'cycles (with end padding)', 'value': int(np.asarray(b.freq_by_cycle).size)},
        {'metric': 'approx. motion duration (s)', 'value': float(np.sum(b.period_by_cycle))},
        {'metric': 'time samples', 'value': int(r['t'].size)},
    ]
    if bool(getattr(b, 'is_stim', False)):
        s1, s2 = b.make_stimuli(is_stim=True, t_basis=r['t'], tnorm_basis=tnorm, stim_pulse_rate=spr)
        r['stim_s1'] = np.asarray(s1, dtype=float).reshape(-1)
        r['stim_s2'] = np.asarray(s2, dtype=float).reshape(-1)
        r['table'].extend(
            [
                {'metric': 'stimulation enabled', 'value': True},
                {'metric': 'pulse rate (Hz)', 'value': spr},
            ]
        )
    else:
        r['table'].append({'metric': 'stimulation enabled', 'value': False})
    return r


def _curves_from_amps(b: Any, mode: str) -> np.ndarray:
    aa = getattr(b, 'all_amps', None)
    if aa is None:
        raise ValueError('Set all_amps for this protocol.')
    conv = b.get_all_amps(aa, mode=mode)
    return np.asarray(conv['curvature_1_per_m'], dtype=float).reshape(-1)


def _preview_frequency_sweep(b: Any, max_plot_points: int) -> PreviewResult:
    r: PreviewResult = {'table': [], 't': None, 'angle': None, 'anglevel': None, 'error': None}
    dur = getattr(b, 'duration', None)
    exp = getattr(b, 'amplitude_frequency_exponent', None)
    if dur is None:
        r['error'] = 'frequency_sweep preview needs duration (s).'
        return r
    if exp is None:
        r['error'] = 'frequency_sweep preview needs amplitude_frequency_exponent.'
        return r
    af = getattr(b, 'all_freqs', None)
    mode = getattr(b, 'all_amps_mode', None) or 'strain'
    curves = _curves_from_amps(b, mode)
    angle, anglevel, tnorm, sweep_freq, t = b.make_cycles_frequency_sweep(
        af,
        curves,
        float(exp),
        float(dur),
        b.waitbefore,
        nominal_frequency=getattr(b, 'sweep_nominal_frequency', None),
        nominal_curvature=getattr(b, 'sweep_nominal_curvature', None),
    )
    r['t'] = np.asarray(t, dtype=float).reshape(-1)
    r['angle'] = np.asarray(angle, dtype=float).reshape(-1)
    r['anglevel'] = np.asarray(anglevel, dtype=float).reshape(-1)
    r['table'] = [
        {'metric': 'motion duration (s)', 'value': float(dur)},
        {'metric': 'samples', 'value': int(r['t'].size)},
    ]
    return r


def _preview_frequency_step(b: Any, max_plot_points: int) -> PreviewResult:
    r: PreviewResult = {'table': [], 't': None, 'angle': None, 'anglevel': None, 'error': None}
    dur = getattr(b, 'duration', None)
    if dur is None:
        r['error'] = 'frequency_step preview needs duration (s).'
        return r
    af = getattr(b, 'all_freqs', None)
    mode = getattr(b, 'all_amps_mode', None) or 'strain'
    curves = _curves_from_amps(b, mode)
    angle, anglevel, tnorm, _, t = b.make_cycles_frequency_step(
        af,
        curves,
        float(dur),
        b.waitbefore,
        nominal_frequency=getattr(b, 'step_nominal_frequency', None),
        nominal_curvature=getattr(b, 'step_nominal_curvature', None),
    )
    r['t'] = np.asarray(t, dtype=float).reshape(-1)
    r['angle'] = np.asarray(angle, dtype=float).reshape(-1)
    r['anglevel'] = np.asarray(anglevel, dtype=float).reshape(-1)
    r['table'] = [{'metric': 'motion duration (s)', 'value': float(dur)}, {'metric': 'samples', 'value': int(r['t'].size)}]
    return r


def _preview_curvature_step(b: Any, max_plot_points: int) -> PreviewResult:
    r: PreviewResult = {'table': [], 't': None, 'angle': None, 'anglevel': None, 'error': None}
    dur = getattr(b, 'duration', None)
    if dur is None:
        r['error'] = 'curvature_step preview needs duration (s).'
        return r
    af = getattr(b, 'all_freqs', None)
    mode = getattr(b, 'all_amps_mode', None) or 'strain'
    curves = _curves_from_amps(b, mode)
    angle, anglevel, tnorm, t = b.make_cycles_curvature_step(
        af,
        curves,
        float(dur),
        b.waitbefore,
        nominal_frequency=getattr(b, 'step_nominal_frequency', None),
        nominal_curvature=getattr(b, 'step_nominal_curvature', None),
    )
    r['t'] = np.asarray(t, dtype=float).reshape(-1)
    r['angle'] = np.asarray(angle, dtype=float).reshape(-1)
    r['anglevel'] = np.asarray(anglevel, dtype=float).reshape(-1)
    r['table'] = [{'metric': 'motion duration (s)', 'value': float(dur)}, {'metric': 'samples', 'value': int(r['t'].size)}]
    return r


def _preview_step_change(b: Any, max_plot_points: int) -> PreviewResult:
    r: PreviewResult = {'table': [], 't': None, 'angle': None, 'anglevel': None, 'error': None}
    freqs = getattr(b, 'step_change_frequencies', None)
    curves = getattr(b, 'step_change_curves', None)
    cps = getattr(b, 'step_change_cycles_per_step', None)
    if freqs is None or curves is None or cps is None:
        r['error'] = 'step_change preview needs step_change_frequencies, step_change_curves, step_change_cycles_per_step.'
        return r
    angle, anglevel, tnorm, t, movedur = b.make_cycles_step_change(
        freqs,
        curves,
        cps,
        dclamp=getattr(b, 'dclamp', None),
        amp_step_vel=getattr(b, 'step_change_amp_step_vel', None),
    )
    r['t'] = np.asarray(t, dtype=float).reshape(-1)
    r['angle'] = np.asarray(angle, dtype=float).reshape(-1)
    r['anglevel'] = np.asarray(anglevel, dtype=float).reshape(-1)
    r['table'] = [
        {'metric': 'computed motion duration (s)', 'value': float(movedur)},
        {'metric': 'samples', 'value': int(r['t'].size)},
    ]
    return r
