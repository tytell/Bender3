"""
Validate, schema-audit, and plot Bender schema-v2 HDF5 exports.

Usage:
  python validate_plot_h5_batch.py FILE.h5 [FILE2.h5 ...]
  python validate_plot_h5_batch.py --acquired-date 2026-06-15
  python validate_plot_h5_batch.py --prefix 2026-06-15
  python validate_plot_h5_batch.py --today   # by file mtime (not acquisition date)
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from datetime import date
from typing import Any, Dict, List, Optional, Tuple

import h5py
import matplotlib.pyplot as plt
import numpy as np

from bender_h5_plot_helpers import _read_input_channel_names, h5_custom_plot_summary
from bender_routing_spec import BENDER_ROUTING, MISSING_REQUIRED

DEFAULT_FOLDER = r'C:\Users\jimen\Desktop\bender_summer_2026'

# Metadata keys marked required=True in bender_routing_spec (ledger canonical names).
REQUIRED_METADATA_KEYS = sorted(
    {route['key'] for route in BENDER_ROUTING.values() if route.get('required') and route.get('tier') == 'metadata'}
)

REQUIRED_ROOT_ATTRS = ('schema_version', 'test_type', 'start_time_iso', 'filename')

REQUIRED_META_DATASETS = (
    'calibration_forcetorque_matrix',
    'daq_instrumentation',
    'trial_index',
)

REQUIRED_TRIAL_DATASETS = (
    'time_second',
    'time_normalized',
    'angle_commanded_degree',
    'angle_measured_degree',
    'angular_velocity_commanded_degree_per_second',
    'forcetorque_raw',
    'aidata',
    'stim_channel1_command_volt',
    'stim_channel2_command_volt',
    'stim_side',
    'stim_state',
    'stim_type',
    'cycle_index',
)

STEP_PROTOCOL_DATASETS = (
    'step_index',
    'sequence_index',
    'block_index',
    'block_direction',
    'block_stim_sides',
)

TORQUE_ROW_LABELS = ('xTorque', 'yTorque', 'zTorque')
TORQUE_ROW_IDX = {'xTorque': 3, 'yTorque': 4, 'zTorque': 5}

FILENAME_PATTERN = re.compile(
    r'^(sim_)?(\d{4}-\d{2}-\d{2})_([A-Za-z0-9]+)_bender_(\d+)_([a-z_]+)\.h5$'
)


def _meta_lookup(f: h5py.File, key: str) -> Any:
    """Return a metadata value from ``01_Metadata`` attrs or top-level datasets."""
    meta = f['01_Metadata']
    if key in meta.attrs:
        return meta.attrs[key]
    if key in meta:
        obj = meta[key]
        if isinstance(obj, h5py.Dataset):
            return obj[()]
    return None


def _is_empty_meta(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, (str, bytes)):
        return len(str(value).strip()) == 0
    arr = np.asarray(value)
    if arr.size == 0:
        return True
    if arr.dtype.kind in ('U', 'S', 'O'):
        flat = [str(x).strip() for x in arr.reshape(-1)]
        return all(len(x) == 0 for x in flat)
    return False


def audit_h5_schema(path: str) -> Dict[str, Any]:
    """
    Check a schema-v2 export against critical contract fields from
    ``context_jlab_cg_h5schema.md`` / ``bender_routing_spec.py``.
    """
    critical: List[str] = []
    warnings: List[str] = []
    info: List[str] = []
    stats: Dict[str, Any] = {'path': path, 'basename': os.path.basename(path)}

    try:
        with h5py.File(path, 'r') as f:
            # --- root attrs ---
            for key in REQUIRED_ROOT_ATTRS:
                if key not in f.attrs:
                    critical.append(f'Root attr missing: {key}')
            stats['schema_version'] = str(f.attrs.get('schema_version', ''))
            stats['test_type'] = str(f.attrs.get('test_type', ''))
            stats['start_time_iso'] = str(f.attrs.get('start_time_iso', ''))

            if '01_Metadata' not in f:
                critical.append('Missing group 01_Metadata')
                return {'ok': False, 'critical': critical, 'warnings': warnings, 'info': info, 'stats': stats}
            if '02_TimeSeries' not in f:
                critical.append('Missing group 02_TimeSeries')
                return {'ok': False, 'critical': critical, 'warnings': warnings, 'info': info, 'stats': stats}

            meta = f['01_Metadata']
            simulated = bool(_meta_lookup(f, 'session_simulated'))
            stats['session_simulated'] = simulated

            # --- required metadata (ledger) ---
            for key in REQUIRED_METADATA_KEYS:
                val = _meta_lookup(f, key)
                if val is None and key not in meta.attrs and key not in meta:
                    critical.append(f'Required metadata missing: {key}')
                elif key == 'specimen_id' and _is_empty_meta(val) and not simulated:
                    critical.append('Required metadata empty: specimen_id')
                elif _is_empty_meta(val) and key not in ('specimen_id',):
                    warnings.append(f'Required metadata empty: {key}')

            for ds in REQUIRED_META_DATASETS:
                if ds not in meta:
                    critical.append(f'Required metadata dataset/group missing: {ds}')

            if 'calibration_link' not in meta:
                warnings.append('Missing subgroup calibration_link (legacy path still usable via attrs)')
            else:
                cl = meta['calibration_link']
                for attr in ('calibration_file', 'calibration_available', 'use_inertial_calibration'):
                    if attr not in cl.attrs:
                        warnings.append(f'calibration_link missing attr: {attr}')

            if 'trial_index' in meta and 'trial_names' not in meta['trial_index']:
                critical.append('trial_index missing trial_names dataset')

            # --- filename convention ---
            basename = os.path.basename(path)
            m = FILENAME_PATTERN.match(basename)
            if not m:
                warnings.append(
                    f'Filename does not match YYYY-MM-DD_<specimen>_bender_<NN>_<protocol>.h5: {basename}'
                )
            else:
                sim_prefix, file_date, file_specimen, _run_num, file_protocol = m.groups()
                stats['filename_date'] = file_date
                stats['filename_specimen'] = file_specimen
                stats['filename_protocol'] = file_protocol
                sid = str(_meta_lookup(f, 'specimen_id') or '').strip()
                if sid and file_specimen != sid:
                    warnings.append(f'Filename specimen {file_specimen!r} != metadata specimen_id {sid!r}')
                if stats['test_type'] and file_protocol != stats['test_type']:
                    warnings.append(
                        f'Filename protocol {file_protocol!r} != test_type {stats["test_type"]!r}'
                    )
                session_date = str(_meta_lookup(f, 'session_date') or '')[:10]
                if session_date and file_date != session_date:
                    warnings.append(
                        f'Filename date {file_date} != session_date {session_date} '
                        '(specimen prep date vs acquisition date?)'
                    )
                if simulated and not sim_prefix:
                    warnings.append('session_simulated=True but filename lacks sim_ prefix')
                if not simulated and sim_prefix:
                    warnings.append('Filename has sim_ prefix but session_simulated=False')

            # --- per-trial timeseries ---
            trials = sorted(k for k in f['02_TimeSeries'].keys() if str(k).startswith('trial_'))
            stats['n_trials'] = len(trials)
            if not trials:
                critical.append('No trial_* groups in 02_TimeSeries')

            step_protocols = {'isometric', 'isovelocity'}
            file_tt = str(stats.get('test_type') or '').lower()

            for tn in trials:
                tg = f['02_TimeSeries'][tn]
                try:
                    t = _trial_time(tg)
                except KeyError:
                    critical.append(f'{tn}: missing time_second')
                    continue
                n = t.size
                if n < 2:
                    critical.append(f'{tn}: fewer than 2 time samples')

                for ds in REQUIRED_TRIAL_DATASETS:
                    if ds not in tg:
                        critical.append(f'{tn}: missing timeseries dataset {ds}')

                if file_tt in step_protocols:
                    for ds in STEP_PROTOCOL_DATASETS:
                        if ds not in tg:
                            warnings.append(f'{tn}: step protocol missing {ds}')

                if 'forcetorque_raw' in tg:
                    ft = np.asarray(tg['forcetorque_raw'][()])
                    if ft.shape != (6, n):
                        critical.append(f'{tn}: forcetorque_raw shape {ft.shape} != (6, {n})')
                    elif not np.any(np.isfinite(ft[5])):
                        critical.append(f'{tn}: zTorque row is all non-finite')

                if 'aidata' in tg:
                    ad = np.asarray(tg['aidata'][()])
                    ch = _read_input_channel_names(f)
                    n_ch = len(ch) if ch else ad.shape[0]
                    if ad.shape[0] != n_ch or ad.shape[1] != n:
                        critical.append(f'{tn}: aidata shape {ad.shape} != ({n_ch}, {n})')

            # --- documented schema gaps (informational, not export bugs) ---
            for gap in MISSING_REQUIRED:
                info.append(f'Known schema gap (not in raw file): {gap["key"]} — {gap.get("note", "")}')
            info.append(
                'Phase-0 on-disk layout uses 01_Metadata / 02_TimeSeries groups '
                '(flat metadata/timeseries is a later restructure per schema doc §2).'
            )
            if _meta_lookup(f, 'measurement_specimen_body_width_millimeter') is None:
                info.append('measurement_specimen_body_width_millimeter absent (no GUI source yet; expected).')

            # --- run-specific flags ---
            if bool(_meta_lookup(f, 'protocol_guard_triggered')):
                guard_deg = _meta_lookup(f, 'protocol_guard_angle_degree')
                warnings.append(f'Isovelocity guard triggered at {guard_deg} deg')

    except OSError as exc:
        critical.append(f'Cannot open file: {exc}')

    ok = len(critical) == 0
    return {'ok': ok, 'critical': critical, 'warnings': warnings, 'info': info, 'stats': stats}


def _stim_monitor_index(channel_names: List[str]) -> Optional[int]:
    for i, name in enumerate(channel_names):
        if str(name).lower() == 'stim_monitor':
            return i
    return None


def _plot_stim_panel(
    ax,
    t: np.ndarray,
    data: Dict[str, np.ndarray],
    *,
    include_command: bool,
) -> None:
    has_cmd = include_command and (
        np.any(np.isfinite(data.get('S1', np.array([]))))
        or np.any(np.isfinite(data.get('S2', np.array([]))))
    )
    has_mon = 'stim_monitor' in data and np.any(np.isfinite(data['stim_monitor']))
    if include_command and np.any(np.isfinite(data.get('S1', np.array([])))):
        ax.plot(t, data['S1'], color='teal', lw=0.8, label='S1 command (V)')
    if include_command and np.any(np.isfinite(data.get('S2', np.array([])))):
        ax.plot(t, data['S2'], '--', color='gray', lw=0.8, label='S2 command (V)')
    if has_mon:
        ax.plot(
            t, data['stim_monitor'], color='crimson', lw=0.9,
            label='stim_monitor (AI, V)',
        )
    if not has_cmd and not has_mon:
        ax.text(
            0.5, 0.5, 'No stim traces in file',
            transform=ax.transAxes, ha='center', va='center', fontsize=9,
        )
    ax.set_ylabel('Stim (V)' if include_command else 'Stim monitor (V)')
    ax.grid(True, alpha=0.3)
    if has_cmd or has_mon:
        ax.legend(loc='upper right', fontsize=8)


def plot_torque_vs_time(
    path: str,
    out_path: Optional[str] = None,
    *,
    include_stim_command: bool = False,
) -> str:
    """Save raw x/y/z torque and stim traces vs stitched time for all trials."""
    with h5py.File(path, 'r') as f:
        data = _stitch_trials(f)
        tt = str(f.attrs.get('test_type', 'unknown'))
        simulated = bool(_meta_lookup(f, 'session_simulated'))
        primary = str(_meta_lookup(f, 'daq_primary_bending_axis') or 'zTorque')
        if primary not in TORQUE_ROW_IDX:
            primary = 'zTorque'

    t = data['t']
    if t.size == 0:
        raise ValueError(f'{path}: no time samples to plot')
    ft = data['forcetorque_raw']
    if ft.size == 0 or ft.shape[0] < 6:
        raise ValueError(f'{path}: missing forcetorque_raw')

    if out_path is None:
        out_path = os.path.splitext(path)[0] + '_torque_vs_time.png'

    series = [
        ('xTorque', ft[3], 'darkorange'),
        ('yTorque', ft[4], 'purple'),
        ('zTorque', ft[5], 'firebrick'),
    ]
    fig, axes = plt.subplots(4, 1, figsize=(12, 9), sharex=True, layout='constrained')
    for ax, (name, y, color) in zip(axes[:3], series):
        lw = 1.2 if name == primary else 0.9
        ax.plot(t, y, color=color, lw=lw, label=f'{name} raw')
        ax.set_ylabel(f'{name} (N·m)')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right', fontsize=8)

    ax_stim = axes[3]
    _plot_stim_panel(ax_stim, t, data, include_command=include_stim_command)
    ax_stim.set_xlabel('Time (s)')

    sim_tag = ' [SIM]' if simulated else ''
    fig.suptitle(
        f'{os.path.basename(path)}{sim_tag}\nTorque + stim — {tt} (all steps stitched)',
        fontsize=11,
    )
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def _pick_dataset(group: h5py.Group, *names: str) -> Optional[str]:
    for nm in names:
        if nm in group:
            return nm
    return None


def _load_sono_cal(f: h5py.File, sono_name: str) -> np.ndarray:
    side = 'right' if 'right' in sono_name.lower() else 'left'
    meta = f['01_Metadata']
    canon = f'calibration_sono_{side}_millimeter_per_volt'
    if canon in meta:
        return np.asarray(meta[canon][()], dtype=float)
    if 'bender_settings' in meta:
        legacy = f'sono_cal_{side}'
        if legacy in meta['bender_settings']:
            return np.asarray(meta['bender_settings'][legacy][()], dtype=float)
    raise KeyError(f'No sono calibration for {sono_name!r}')


def _sono_voltage_span(cal: np.ndarray) -> Tuple[float, float]:
    cal = np.asarray(cal, dtype=float).reshape(-1)
    if cal.size >= 4 and cal.size % 2 == 0:
        v = cal[: cal.size // 2]
        return float(np.min(v)), float(np.max(v))
    return float(cal[0]), float(cal[1])


def _sono_mm_from_volts(volts: np.ndarray, cal: np.ndarray) -> np.ndarray:
    cal = np.asarray(cal, dtype=float).reshape(-1)
    if cal.size >= 4 and cal.size % 2 == 0:
        n = cal.size // 2
        v_pts, mm_pts = cal[:n], cal[n:]
        slope, intercept = np.polyfit(v_pts, mm_pts, 1)
        return np.asarray(volts, dtype=float) * slope + intercept
    v_lo, v_hi = float(cal[0]), float(cal[1])
    mm_lo, mm_hi = float(cal[2]), float(cal[3])
    slope = (mm_hi - mm_lo) / (v_hi - v_lo)
    intercept = mm_lo - slope * v_lo
    return np.asarray(volts, dtype=float) * slope + intercept


def _inter_trial_gap_s(f: h5py.File) -> float:
    meta = f['01_Metadata']
    for group_name in ('bender_settings', 'protocol_metadata'):
        if group_name not in meta:
            continue
        gs = meta[group_name]
        for key in (
            'isometric_inter_step_interval_s',
            'isovelocity_inter_step_interval_s',
            'inter_step_interval_s',
        ):
            if key in gs.attrs:
                try:
                    return float(gs.attrs[key] or 0.0)
                except (TypeError, ValueError):
                    pass
    return 0.0


def _trial_time(group: h5py.Group) -> np.ndarray:
    key = _pick_dataset(group, 'time_second', 't')
    if not key:
        raise KeyError('missing time_second / t')
    return np.asarray(group[key][()], dtype=float).reshape(-1)


def validate_h5(path: str) -> Dict:
    issues: List[str] = []
    warnings: List[str] = []
    stats: Dict = {'path': path, 'basename': os.path.basename(path)}

    try:
        with h5py.File(path, 'r') as f:
            summ = h5_custom_plot_summary(path)
            if not summ.get('ok'):
                issues.append(str(summ.get('error') or 'h5_custom_plot_summary failed'))
                return {'ok': False, 'issues': issues, 'warnings': warnings, 'stats': stats}

            stats.update(
                {
                    'test_type': summ['test_type'],
                    'schema_version': summ['schema_version'],
                    'n_trials': summ['n_trials'],
                    'channels': summ['channel_names'],
                    'primary_axis': summ['primary_bending_axis'],
                    'start_time_iso': str(f.attrs.get('start_time_iso', '')),
                }
            )

            ch = _read_input_channel_names(f)
            if len(ch) != summ['n_trials'] and len(ch) < 6:
                warnings.append(f'Unexpected channel count: {len(ch)}')

            trials = sorted(k for k in f['02_TimeSeries'].keys() if str(k).startswith('trial_'))
            if not trials:
                issues.append('No trial_* groups in 02_TimeSeries')

            sono_idx = next((i for i, n in enumerate(ch) if str(n).lower().startswith('sono_')), None)
            sono_name = ch[sono_idx] if sono_idx is not None else None
            v_lo = v_hi = None
            if sono_name:
                try:
                    cal = _load_sono_cal(f, sono_name)
                    v_lo, v_hi = _sono_voltage_span(cal)
                    stats['sono_cal_points'] = int(cal.size // 2) if cal.size % 2 == 0 else 2
                    stats['sono_v_span'] = [v_lo, v_hi]
                except KeyError as e:
                    warnings.append(str(e))

            axis_idx = {'xTorque': 3, 'yTorque': 4, 'zTorque': 5}.get(
                str(stats.get('primary_axis') or 'zTorque'), 5
            )

            trial_durs: List[float] = []
            z_finite: List[float] = []
            sono_in_cal: List[float] = []
            angle_ranges: List[float] = []

            for tn in trials:
                tg = f['02_TimeSeries'][tn]
                try:
                    t = _trial_time(tg)
                except KeyError:
                    issues.append(f'{tn}: missing time axis')
                    continue
                n = t.size
                if n < 2:
                    issues.append(f'{tn}: fewer than 2 samples')
                    continue
                trial_durs.append(float(t[-1] - t[0]))

                for ds, expected_ndim in (
                    ('forcetorque_raw', 2),
                    ('aidata', 2),
                ):
                    if ds not in tg:
                        issues.append(f'{tn}: missing {ds}')
                    else:
                        arr = np.asarray(tg[ds][()])
                        if expected_ndim == 2 and arr.ndim == 2:
                            width = arr.shape[1] if arr.shape[0] in (6, 8) else arr.shape[0]
                            if width != n:
                                issues.append(f'{tn}: {ds} width {width} != time n={n}')

                ang_key = _pick_dataset(tg, 'angle_measured_degree', 'angle_measured')
                if ang_key:
                    ang = np.asarray(tg[ang_key][()], dtype=float)
                    if ang.size == n and np.any(np.isfinite(ang)):
                        angle_ranges.append(float(np.nanmax(ang) - np.nanmin(ang)))
                else:
                    warnings.append(f'{tn}: no measured angle dataset')

                if 'forcetorque_raw' in tg:
                    ft = np.asarray(tg['forcetorque_raw'][()], dtype=float)
                    if ft.ndim == 2 and ft.shape[0] >= 6:
                        z_finite.append(float(np.mean(np.isfinite(ft[axis_idx]))))

                if sono_idx is not None and v_lo is not None and 'aidata' in tg:
                    v = np.asarray(tg['aidata'][()][sono_idx], dtype=float)
                    m = np.isfinite(v) & (v >= v_lo) & (v <= v_hi)
                    sono_in_cal.append(float(np.mean(m)))

            stats['trial_duration_s'] = trial_durs
            stats['total_trial_duration_s'] = round(float(sum(trial_durs)), 3)
            stats['ztorque_finite_frac'] = round(float(np.mean(z_finite)), 4) if z_finite else None
            stats['sono_in_cal_frac'] = round(float(np.mean(sono_in_cal)), 4) if sono_in_cal else None
            stats['angle_range_deg'] = angle_ranges

            if stats.get('ztorque_finite_frac') is not None and stats['ztorque_finite_frac'] < 0.9:
                warnings.append(f'Low primary torque finite fraction: {stats["ztorque_finite_frac"]}')
            if stats.get('sono_in_cal_frac') is not None and stats['sono_in_cal_frac'] < 0.5:
                warnings.append(
                    f'Sono mostly outside cal span [{v_lo:.2g}, {v_hi:.2g}]: '
                    f'{stats["sono_in_cal_frac"] * 100:.1f}% in-cal'
                )

            tt = str(stats.get('test_type') or '')
            if tt and tt not in os.path.basename(path).lower():
                warnings.append(f'Filename does not contain test_type={tt!r}')

    except OSError as e:
        issues.append(f'Cannot open file: {e}')

    ok = len(issues) == 0
    return {'ok': ok, 'issues': issues, 'warnings': warnings, 'stats': stats}


def _stitch_trials(f: h5py.File) -> Dict[str, np.ndarray]:
    gap = _inter_trial_gap_s(f)
    g_ts = f['02_TimeSeries']
    trials = sorted(k for k in g_ts.keys() if str(k).startswith('trial_'))
    t_parts: List[np.ndarray] = []
    ang_cmd: List[np.ndarray] = []
    ang_meas: List[np.ndarray] = []
    ang_vel: List[np.ndarray] = []
    ft_parts: List[np.ndarray] = []
    s1_parts: List[np.ndarray] = []
    s2_parts: List[np.ndarray] = []
    stim_mon_parts: List[np.ndarray] = []
    sono_parts: List[np.ndarray] = []
    ch = _read_input_channel_names(f)
    stim_mon_idx = _stim_monitor_index(ch)
    sono_idx = next((i for i, n in enumerate(ch) if str(n).lower().startswith('sono_')), None)
    cal = None
    if sono_idx is not None:
        sono_name = ch[sono_idx]
        cal = _load_sono_cal(f, sono_name)

    t_off = 0.0
    for tn in trials:
        tg = g_ts[tn]
        t = _trial_time(tg)
        n = t.size
        if n == 0:
            continue
        t_parts.append(t - t[0] + t_off)

        def _read_1d(names: Tuple[str, ...]) -> np.ndarray:
            key = _pick_dataset(tg, *names)
            if not key:
                return np.full(n, np.nan)
            a = np.asarray(tg[key][()], dtype=float).reshape(-1)
            return a if a.size == n else np.full(n, np.nan)

        ang_cmd.append(_read_1d(('angle_commanded_degree', 'angle_cmd')))
        ang_meas.append(_read_1d(('angle_measured_degree', 'angle_measured')))
        ang_vel.append(_read_1d(('angular_velocity_commanded_degree_per_second', 'anglevel_cmd')))
        ft = np.asarray(tg['forcetorque_raw'][()], dtype=float) if 'forcetorque_raw' in tg else None
        ft_parts.append(ft if (ft is not None and ft.shape == (6, n)) else np.full((6, n), np.nan))
        s1k = _pick_dataset(tg, 'stim_channel1_command_volt', 'S1stimcmd')
        s2k = _pick_dataset(tg, 'stim_channel2_command_volt', 'S2stimcmd')
        s1_parts.append(np.asarray(tg[s1k][()], float).reshape(-1) if s1k else np.full(n, np.nan))
        s2_parts.append(np.asarray(tg[s2k][()], float).reshape(-1) if s2k else np.full(n, np.nan))
        if stim_mon_idx is not None and 'aidata' in tg:
            ad = np.asarray(tg['aidata'][()], dtype=float)
            sm = ad[stim_mon_idx].reshape(-1) if ad.shape[0] > stim_mon_idx else np.array([])
            stim_mon_parts.append(sm if sm.size == n else np.full(n, np.nan))
        else:
            stim_mon_parts.append(np.full(n, np.nan))
        if sono_idx is not None and cal is not None:
            v = np.asarray(tg['aidata'][()][sono_idx], dtype=float)
            L = _sono_mm_from_volts(v, cal)
            v_lo, v_hi = _sono_voltage_span(cal)
            L = np.where((v >= v_lo) & (v <= v_hi), L, np.nan)
            sono_parts.append(L)
        dt = (t[-1] - t[0]) / max(1, n - 1) if n > 1 else 0.0
        t_off = float(t_parts[-1][-1] + dt + gap)

    out = {
        't': np.concatenate(t_parts) if t_parts else np.array([]),
        'angle_cmd': np.concatenate(ang_cmd) if ang_cmd else np.array([]),
        'angle_meas': np.concatenate(ang_meas) if ang_meas else np.array([]),
        'ang_vel': np.concatenate(ang_vel) if ang_vel else np.array([]),
        'forcetorque_raw': np.concatenate(ft_parts, axis=1) if ft_parts else np.array([]),
        'S1': np.concatenate(s1_parts) if s1_parts else np.array([]),
        'S2': np.concatenate(s2_parts) if s2_parts else np.array([]),
        'stim_monitor': np.concatenate(stim_mon_parts) if stim_mon_parts else np.array([]),
    }
    if sono_parts:
        out['sono_mm'] = np.concatenate(sono_parts)
    return out


def plot_qc_from_h5(path: str, out_path: Optional[str] = None) -> str:
    with h5py.File(path, 'r') as f:
        tt = str(f.attrs.get('test_type', 'unknown'))
        data = _stitch_trials(f)
        primary = str(h5_custom_plot_summary(path).get('primary_bending_axis') or 'zTorque')
        axis_idx = {'xTorque': 3, 'yTorque': 4, 'zTorque': 5}.get(primary, 5)
        off = [a for a in ('xTorque', 'yTorque', 'zTorque') if a != primary]

    t = data['t']
    has_sono = 'sono_mm' in data and np.any(np.isfinite(data['sono_mm']))
    n_rows = 6 if has_sono else 5
    fig, axes = plt.subplots(n_rows, 1, figsize=(12, 2.2 * n_rows), sharex=True, layout='constrained')
    if n_rows == 1:
        axes = [axes]

    axes[0].plot(t, data['angle_cmd'], '--', color='black', lw=0.8, label='angle_cmd')
    axes[0].plot(t, data['angle_meas'], color='royalblue', lw=0.8, label='angle_measured')
    ax0r = axes[0].twinx()
    ax0r.plot(t, data['ang_vel'], color='orange', lw=0.6, label='anglevel_cmd')
    axes[0].set_ylabel('Angle (deg)')
    ax0r.set_ylabel('Velocity (deg/s)')
    axes[0].legend(loc='upper right', fontsize=8)
    axes[0].grid(True, alpha=0.3)

    ft = data['forcetorque_raw']
    if ft.size:
        axes[1].plot(t, ft[axis_idx], color='firebrick', lw=0.8, label=f'{primary} raw')
        axes[1].set_ylabel(f'{primary} (N·m)')
        axes[1].legend(loc='upper right', fontsize=8)
        axes[1].grid(True, alpha=0.3)
        off_idx = {'xTorque': 3, 'yTorque': 4, 'zTorque': 5}
        for i, ax_name in enumerate(off[:2]):
            axes[2 + i].plot(t, ft[off_idx[ax_name]], lw=0.8, label=f'{ax_name} raw')
            axes[2 + i].set_ylabel(f'{ax_name} (N·m)')
            axes[2 + i].legend(loc='upper right', fontsize=8)
            axes[2 + i].grid(True, alpha=0.3)

    stim_row = 4
    if 'stim_monitor' in data and np.any(np.isfinite(data['stim_monitor'])):
        axes[stim_row].plot(
            t, data['stim_monitor'], color='crimson', lw=0.9,
            label='stim_monitor (AI, V)',
        )
    else:
        axes[stim_row].text(
            0.5, 0.5, 'No stim_monitor channel in file',
            transform=axes[stim_row].transAxes, ha='center', va='center', fontsize=9,
        )
    axes[stim_row].set_ylabel('Stim monitor (V)')
    axes[stim_row].legend(loc='upper right', fontsize=8)
    axes[stim_row].grid(True, alpha=0.3)

    if has_sono:
        axes[5].plot(t, data['sono_mm'], color='darkcyan', lw=0.8, label='sono length')
        axes[5].set_ylabel('Sono (mm)')
        axes[5].legend(loc='upper right', fontsize=8)
        axes[5].grid(True, alpha=0.3)
        axes[5].set_xlabel('Time (s)')
    else:
        axes[stim_row].set_xlabel('Time (s)')

    fig.suptitle(f'QC (from H5): {tt} | all steps\n{os.path.basename(path)}', fontsize=11)
    if out_path is None:
        out_path = os.path.splitext(path)[0] + '_qc_from_h5.png'
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def plot_sono_from_h5(path: str, out_path: Optional[str] = None, *, smooth_s: float = 0.0) -> str:
    with h5py.File(path, 'r') as f:
        data = _stitch_trials(f)
    t = data['t']
    L = data.get('sono_mm')
    if L is None or not np.any(np.isfinite(L)):
        raise ValueError(f'{path}: no in-cal sono length to plot')
    if smooth_s > 0 and t.size > 2:
        dt = float(np.median(np.diff(t[np.isfinite(t)])))
        if np.isfinite(dt) and dt > 0:
            w = int(max(3, round(smooth_s / dt)))
            kernel = np.ones(w) / w
            filled = np.where(np.isfinite(L), L, np.nan)
            sm = np.convolve(np.nan_to_num(filled, nan=np.nanmean(L)), kernel, mode='same')
            L = np.where(np.isfinite(L), sm, np.nan)
    if out_path is None:
        out_path = os.path.splitext(path)[0] + '_sono_length_vs_time.png'
    fig, ax = plt.subplots(figsize=(12, 4.5), layout='constrained')
    ax.plot(t, L, color='#0ea5e9', lw=1.0)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Muscle length (mm)')
    ax.set_title(os.path.basename(path))
    ax.grid(True, alpha=0.3)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def files_acquired_on(folder: str, day_prefix: str) -> List[str]:
    """Return ``.h5`` paths whose root ``start_time_iso`` begins with ``day_prefix`` (YYYY-MM-DD)."""
    out: List[str] = []
    for name in os.listdir(folder):
        if not name.lower().endswith('.h5'):
            continue
        full = os.path.join(folder, name)
        if not os.path.isfile(full):
            continue
        try:
            with h5py.File(full, 'r') as f:
                st = str(f.attrs.get('start_time_iso', ''))
                if st.startswith(day_prefix):
                    out.append(full)
        except OSError:
            continue
    return sorted(out)


def files_modified_on(folder: str, day: date) -> List[str]:
    out: List[str] = []
    for name in os.listdir(folder):
        if not name.lower().endswith('.h5'):
            continue
        full = os.path.join(folder, name)
        if os.path.isfile(full) and date.fromtimestamp(os.path.getmtime(full)) == day:
            out.append(full)
    return sorted(out)


def files_with_prefix(folder: str, prefix: str, *, recursive: bool = True) -> List[str]:
    """Return ``.h5`` paths whose basename starts with ``prefix``."""
    out: List[str] = []
    if not prefix:
        return out
    if recursive:
        for dirpath, _, filenames in os.walk(folder):
            for name in filenames:
                if name.lower().endswith('.h5') and name.startswith(prefix):
                    out.append(os.path.join(dirpath, name))
    else:
        for name in os.listdir(folder):
            if name.lower().endswith('.h5') and name.startswith(prefix):
                full = os.path.join(folder, name)
                if os.path.isfile(full):
                    out.append(full)
    return sorted(out)


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('files', nargs='*', help='HDF5 paths')
    p.add_argument(
        '--prefix',
        default='',
        help='Plot all .h5 files in --folder whose names start with this string (e.g. 2026-06-15)',
    )
    p.add_argument(
        '--acquired-date',
        default='',
        metavar='YYYY-MM-DD',
        help='All .h5 in --folder whose start_time_iso begins with this date',
    )
    p.add_argument('--today', action='store_true', help='Shortcut for --acquired-date=<today>')
    p.add_argument('--folder', default=DEFAULT_FOLDER, help='Folder for --prefix / --acquired-date')
    p.add_argument('--smooth-sono', type=float, default=0.0, metavar='SEC', help='Sono moving avg (s)')
    p.add_argument('--schema-audit', action='store_true', help='Run full schema contract audit')
    p.add_argument('--torque', action='store_true', help='Save torque vs time PNG for each file')
    p.add_argument(
        '--stim-command',
        action='store_true',
        help='On torque plots, include S1/S2 command voltages alongside stim_monitor',
    )
    p.add_argument('--no-qc', action='store_true')
    p.add_argument('--no-sono', action='store_true')
    p.add_argument('--report', default='', help='Write combined audit report to this text path')
    args = p.parse_args(argv)

    acquired = str(args.acquired_date or '').strip()
    if args.today and not acquired:
        acquired = date.today().isoformat()

    paths = list(args.files)
    if args.prefix:
        paths.extend(files_with_prefix(args.folder, args.prefix))
    if acquired:
        paths.extend(files_acquired_on(args.folder, acquired))
    paths = sorted(set(os.path.normpath(p) for p in paths if p))
    if not paths:
        if args.prefix:
            print(
                f'No .h5 files starting with {args.prefix!r} under {args.folder}',
                file=sys.stderr,
            )
        elif acquired:
            print(
                f'No .h5 files with start_time_iso beginning {acquired!r} under {args.folder}',
                file=sys.stderr,
            )
        p.error('Provide FILE.h5 paths, --prefix, or --acquired-date')

    report_lines: List[str] = []
    n_ok = 0
    n_schema_ok = 0
    for path in paths:
        print('=' * 72)
        print(os.path.basename(path))
        report_lines.append('=' * 72)
        report_lines.append(os.path.basename(path))
        if not os.path.isfile(path):
            print('  SKIP — file not found')
            report_lines.append('  SKIP — file not found')
            continue

        if args.schema_audit:
            audit = audit_h5_schema(path)
            schema_status = 'PASS' if audit['ok'] else 'FAIL'
            print(f'  Schema audit: {schema_status}')
            report_lines.append(f'  Schema audit: {schema_status}')
            st = audit['stats']
            print(f'  schema={st.get("schema_version")}  test_type={st.get("test_type")}  '
                  f'trials={st.get("n_trials")}  start={st.get("start_time_iso")}')
            report_lines.append(
                f'  schema={st.get("schema_version")}  test_type={st.get("test_type")}  '
                f'trials={st.get("n_trials")}  start={st.get("start_time_iso")}'
            )
            for msg in audit.get('critical', []):
                print(f'  CRITICAL: {msg}')
                report_lines.append(f'  CRITICAL: {msg}')
            for msg in audit.get('warnings', []):
                print(f'  WARN: {msg}')
                report_lines.append(f'  WARN: {msg}')
            for msg in audit.get('info', [])[:3]:
                print(f'  INFO: {msg}')
            if audit['ok']:
                n_schema_ok += 1

        rep = validate_h5(path)
        status = 'PASS' if rep['ok'] else 'FAIL'
        print(f'  Data validation: {status}')
        report_lines.append(f'  Data validation: {status}')
        st = rep['stats']
        print(f'  duration={st.get("total_trial_duration_s")} s')
        if st.get('sono_in_cal_frac') is not None:
            print(f'  sono in-cal: {st["sono_in_cal_frac"] * 100:.1f}%')
        if st.get('ztorque_finite_frac') is not None:
            print(f'  primary torque finite: {st["ztorque_finite_frac"] * 100:.1f}%')
        if rep['issues']:
            for msg in rep['issues']:
                print(f'  ISSUE: {msg}')
                report_lines.append(f'  ISSUE: {msg}')
        if rep['warnings']:
            for msg in rep['warnings']:
                print(f'  WARN: {msg}')
                report_lines.append(f'  WARN: {msg}')
        if rep['ok']:
            n_ok += 1

        if args.torque:
            try:
                tp = plot_torque_vs_time(path, include_stim_command=args.stim_command)
                print(f'  Saved torque plot: {tp}')
                report_lines.append(f'  Saved torque plot: {tp}')
            except Exception as e:
                print(f'  Torque plot failed: {e}', file=sys.stderr)
                report_lines.append(f'  Torque plot failed: {e}')

        if not args.no_qc:
            try:
                qc = plot_qc_from_h5(path)
                print(f'  Saved QC plot: {qc}')
            except Exception as e:
                print(f'  QC plot failed: {e}', file=sys.stderr)
        if not args.no_sono:
            tt = str(st.get('test_type') or '').lower()
            smooth = args.smooth_sono or (0.2 if 'dynamic' in tt or 'isovelocity' in tt else 0.0)
            try:
                sono = plot_sono_from_h5(path, smooth_s=smooth)
                print(f'  Saved sono plot: {sono}')
            except Exception as e:
                print(f'  Sono plot skipped: {e}')

    print('=' * 72)
    if args.schema_audit:
        print(f'Schema audit: {n_schema_ok}/{len(paths)} passed')
    print(f'Data validation: {n_ok}/{len(paths)} passed')
    if args.report:
        with open(args.report, 'w', encoding='utf-8') as fh:
            fh.write('\n'.join(report_lines))
            fh.write('\n')
        print(f'Wrote report: {args.report}')
    return 0 if n_ok == len(paths) else 1


if __name__ == '__main__':
    raise SystemExit(main())
