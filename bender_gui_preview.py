"""
Protocol preview for the Streamlit GUI: motor command timeline without DAQ.

Mirrors motion-generation branches in :meth:`bender_functions.Bender.run_experiment`
(isometric / isovelocity / calibration base / dynamic / sweep / step_change).
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np

PreviewResult = Dict[str, Any]


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
        't_plot': None,
        'angle_plot': None,
        'anglevel_plot': None,
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
        if t is not None and ang is not None and av is not None:
            tp, ap = _downsample(t, ang, max_plot_points)
            _, vp = _downsample(t, av, max_plot_points)
            out['t_plot'] = tp
            out['angle_plot'] = ap
            out['anglevel_plot'] = vp
        return out
    except Exception as e:
        out['error'] = f'{type(e).__name__}: {e}'
        return out


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
        xs = np.linspace(float(initial), float(final), n)
        for i, v in enumerate(xs):
            r['table'].append({'step_index': i, 'target_value': float(v), 'input_mode': str(mode)})
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
    r['table'].append({'row': 'starting strain', 'value': float(ss), 'unit': str(mode)})
    for i, v in enumerate(vels):
        r['table'].append({'row': f'velocity step {i}', 'value': float(v), 'unit': 'deg/s'})
    return r


def _preview_dynamic(b: Any, max_plot_points: int) -> PreviewResult:
    r: PreviewResult = {'table': [], 't': None, 'angle': None, 'anglevel': None, 'error': None}
    dc = getattr(b, 'dclamp', None)
    xw = getattr(b, 'xsec_width', None)
    if dc is None:
        r['error'] = 'Dynamic preview needs dclamp (mm) on the Bender — usually set with specimen metadata.'
        return r
    af = _as_float_list(getattr(b, 'all_freqs', None))
    aa = getattr(b, 'all_amps', None)
    mode = getattr(b, 'all_amps_mode', None) or 'strain'
    if xw is None:
        r['error'] = (
            'Dynamic preview needs xsec_width (mm) on the Bender (organize_cycles uses it for strain metadata).'
        )
        return r
    if not af:
        r['error'] = 'Set all_freqs (Hz list) for dynamic preview.'
        return r
    if aa is None:
        r['error'] = 'Set all_amps for dynamic preview.'
        return r
    conv = b.get_all_amps(aa, mode=mode)
    all_curves = np.asarray(conv['curvature_1_per_m'], dtype=float).reshape(-1)
    randomize = bool(getattr(b, 'randomize', False))
    rs = getattr(b, 'random_seed', None)
    if randomize and rs is not None:
        np.random.seed(int(rs))
    cps = int(getattr(b, 'cycles_per_step', 0) or 0)
    nec = int(getattr(b, 'n_end_cycles', 0) or 0)
    if cps <= 0:
        r['error'] = 'cycles_per_step must be a positive integer.'
        return r
    stim_ix = getattr(b, 'stim_cycles_in_step', None)
    if stim_ix is None:
        stim_ix = []
    stim_ix = np.asarray(stim_ix, dtype=int).reshape(-1).tolist()
    duties, phases = _default_stim_duties_phases(b)
    spr = float(getattr(b, 'stim_pulse_rate', 0.0) or 0.0)
    b.organize_cycles(
        list(all_curves),
        af,
        randomize,
        cps,
        nec,
        float(dc),
        float(xw),
        stim_ix,
        duties,
        phases,
        spr,
    )
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
