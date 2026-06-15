"""
Plot measured sonomicrometry muscle length (mm) vs time from Bender v2 HDF5 exports.

QC options (recommended for noisy / extrapolated traces):
  - Full two-point calibration: L = slope*V + intercept from [V_lo, V_hi, mm_lo, mm_hi]
  - Mask samples outside the calibrated voltage span (gaps in the line)
  - Optional hold-only window when trial attrs define pre/active/post baselines

Usage:
  python plot_sono_length_vs_time.py FILE.h5 [FILE2.h5 ...]
  python plot_sono_length_vs_time.py --fuzzy-batch   # re-plot known noisy summer-2026 files
"""
from __future__ import annotations

import argparse
import os
import sys
from typing import List, Optional, Tuple

import h5py
import matplotlib.pyplot as plt
import numpy as np

from bender_h5_plot_helpers import _read_input_channel_names

FUZZY_BATCH = [
    r'C:\Users\jimen\Desktop\bender_summer_2026\2026-06-04_snapper01_bender_15_isometric.h5',
    r'C:\Users\jimen\Desktop\bender_summer_2026\2026-06-04_snapper01_bender_18_isovelocity.h5',
    r'C:\Users\jimen\Desktop\bender_summer_2026\2026-06-04_snapper01_bender_20_dynamic.h5',
]


def sono_cal_params(cal: np.ndarray) -> Tuple[float, float, float, float]:
    """Return (v_lo, v_hi, slope, intercept) for 2- or multi-point cal tables."""
    cal = np.asarray(cal, dtype=float).reshape(-1)
    if cal.size >= 4 and cal.size % 2 == 0:
        n = cal.size // 2
        volts, mm = cal[:n], cal[n:]
        v_lo, v_hi = float(np.min(volts)), float(np.max(volts))
        slope, intercept = np.polyfit(volts, mm, 1)
        return v_lo, v_hi, float(slope), float(intercept)
    v_lo, v_hi, mm_lo, mm_hi = [float(x) for x in cal[:4]]
    slope = (mm_hi - mm_lo) / (v_hi - v_lo)
    intercept = mm_lo - slope * v_lo
    return v_lo, v_hi, slope, intercept


def sono_mm_from_volts(volts: np.ndarray, cal: np.ndarray, *, mode: str = 'full') -> np.ndarray:
    v_lo, v_hi, slope, intercept = sono_cal_params(cal)
    v = np.asarray(volts, dtype=float)
    if mode == 'slope':
        return v * slope
    if mode == 'full':
        return v * slope + intercept
    raise ValueError(f"Unknown cal mode {mode!r}; use 'full' or 'slope'.")


def in_cal_mask(volts: np.ndarray, cal: np.ndarray) -> np.ndarray:
    v_lo, v_hi, _, _ = sono_cal_params(cal)
    v = np.asarray(volts, dtype=float)
    return (v >= v_lo) & (v <= v_hi)


def trial_hold_mask(t_local: np.ndarray, tg: h5py.Group) -> Optional[np.ndarray]:
    """True during pre + active + post baseline hold if metadata present."""
    if 't_pre_baseline_start' not in tg.attrs:
        return None
    t = np.asarray(t_local, dtype=float)
    t0 = float(tg.attrs['t_pre_baseline_start'])
    t1 = float(tg.attrs.get('t_post_baseline_end', t[-1]))
    return (t >= t0) & (t <= t1)


def moving_average(y: np.ndarray, window: int) -> np.ndarray:
    """Centered moving average; leaves NaNs as NaN."""
    y = np.asarray(y, dtype=float)
    w = int(max(1, window))
    if w <= 1:
        return y
    out = y.copy()
    fin = np.isfinite(y)
    if np.count_nonzero(fin) < w:
        return y
    filled = np.where(fin, y, np.interp(np.arange(y.size), np.flatnonzero(fin), y[fin]))
    kernel = np.ones(w, dtype=float) / w
    sm = np.convolve(filled, kernel, mode='same')
    out[fin] = sm[fin]
    return out


def _trial_time_seconds(tg: h5py.Group) -> np.ndarray:
    for key in ('time_second', 't'):
        if key in tg:
            return np.asarray(tg[key][()], dtype=float)
    raise KeyError('trial group has no time_second or t dataset')


def _load_sono_cal_table(f: h5py.File, sono_name: str) -> Tuple[np.ndarray, str]:
    side = 'right' if 'right' in sono_name.lower() else 'left'
    meta = f['01_Metadata']
    canon = f'calibration_sono_{side}_millimeter_per_volt'
    if canon in meta:
        return np.asarray(meta[canon][()], dtype=float), canon
    gs = meta.get('bender_settings')
    legacy = f'sono_cal_{side}'
    if gs is not None and legacy in gs:
        return np.asarray(gs[legacy][()], dtype=float), legacy
    raise KeyError(f'no sono calibration for {sono_name!r}')


def inter_trial_gap_s(f: h5py.File) -> float:
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


def load_stitched_series(
    path: str,
    *,
    cal_mode: str = 'full',
    hold_only: bool = False,
    in_cal_only: bool = False,
) -> Tuple[np.ndarray, np.ndarray, str, dict]:
    meta: dict = {}
    with h5py.File(path, 'r') as f:
        if '02_TimeSeries' not in f:
            raise ValueError(f'{path}: no 02_TimeSeries group')
        ch = _read_input_channel_names(f)
        sono = next((n for n in ch if str(n).lower().startswith('sono_')), None)
        if not sono:
            raise ValueError(f'{path}: no sono_* channel in {ch}')
        idx = ch.index(sono)
        cal, cal_key = _load_sono_cal_table(f, sono)
        v_lo, v_hi, slope, intercept = sono_cal_params(cal)
        meta.update(
            {
                'sono_name': sono,
                'cal_key': cal_key,
                'v_lo': v_lo,
                'v_hi': v_hi,
                'mm_lo': float(cal[2]),
                'mm_hi': float(cal[3]),
                'slope': slope,
                'intercept': intercept,
                'cal_mode': cal_mode,
                'hold_only': hold_only,
                'in_cal_only': in_cal_only,
            }
        )
        gap = inter_trial_gap_s(f)
        g_ts = f['02_TimeSeries']
        trials = sorted(k for k in g_ts.keys() if str(k).startswith('trial_'))
        meta['n_trials'] = len(trials)
        t_parts: List[np.ndarray] = []
        L_parts: List[np.ndarray] = []
        t_off = 0.0
        for tn in trials:
            tg = g_ts[tn]
            t_local = _trial_time_seconds(tg)
            v = np.asarray(tg['aidata'][()][idx, :], dtype=float)
            L = sono_mm_from_volts(v, cal, mode=cal_mode)
            mask = np.ones(t_local.shape, dtype=bool)
            if hold_only:
                hm = trial_hold_mask(t_local, tg)
                if hm is not None:
                    mask &= hm
            if in_cal_only:
                mask &= in_cal_mask(v, cal)
            L_plot = np.where(mask, L, np.nan)
            t_parts.append(t_local + t_off)
            L_parts.append(L_plot)
            t_off = float(t_local[-1]) + t_off + gap
    return np.concatenate(t_parts), np.concatenate(L_parts), sono, meta


def plot_and_save(
    path: str,
    *,
    cal_mode: str = 'full',
    hold_only: bool = True,
    in_cal_only: bool = True,
    min_finite_frac: float = 0.08,
    smooth_window_s: float = 0.0,
    out_path: Optional[str] = None,
) -> str:
    t, L, sono, meta = load_stitched_series(
        path, cal_mode=cal_mode, hold_only=hold_only, in_cal_only=in_cal_only
    )
    finite_frac = float(np.mean(np.isfinite(L))) if L.size else 0.0
    sparse_in_cal = bool(in_cal_only and finite_frac < min_finite_frac)

    if out_path is None:
        out_path = os.path.splitext(path)[0] + '_sono_length_vs_time.png'

    parts = [os.path.basename(path)]
    if cal_mode == 'full':
        parts.append(
            f"L = {meta['slope']:.4g}·V + {meta['intercept']:.2g} mm"
        )
    else:
        parts.append(f"L = V × {meta['slope']:.4g} mm/V")
    notes = []
    if hold_only:
        notes.append('hold windows only')
    if in_cal_only:
        notes.append(f"V in [{meta['v_lo']:.2g}, {meta['v_hi']:.2g}]")
    if sparse_in_cal:
        notes.append(f"sparse in-cal ({finite_frac * 100:.1f}% of samples)")
    if notes:
        parts.append('; '.join(notes))

    fig, ax = plt.subplots(figsize=(12, 4.5), layout='constrained')
    if sparse_in_cal:
        m = np.isfinite(L)
        ax.plot(
            t[m], L[m], color='#0ea5e9', marker='.', ms=3,
            linestyle='none', label=f'{sono} length (in-cal)',
        )
    else:
        ax.plot(t, L, color='#0ea5e9', lw=1.0, label=f'{sono} length')
    if smooth_window_s > 0 and t.size > 2 and not sparse_in_cal:
        dt = float(np.median(np.diff(t[np.isfinite(t)])))
        if np.isfinite(dt) and dt > 0:
            w = int(max(3, round(smooth_window_s / dt)))
            L = moving_average(L, w)
            parts.append(f'{smooth_window_s * 1000:.0f} ms mov. avg.')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Muscle length (mm)')
    ax.set_title('\n'.join(parts))
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right')
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('files', nargs='*', help='HDF5 file paths')
    p.add_argument(
        '--fuzzy-batch',
        action='store_true',
        help='Re-plot known noisy summer-2026 files with QC defaults',
    )
    p.add_argument(
        '--cal',
        choices=('full', 'slope'),
        default='full',
        help='full = slope*V+intercept; slope = V*mm/V only',
    )
    p.add_argument('--no-hold-only', action='store_true', help='Include motor ramps')
    p.add_argument('--no-in-cal-only', action='store_true', help='Do not mask outside cal V')
    p.add_argument(
        '--smooth',
        type=float,
        default=0.0,
        metavar='SEC',
        help='Moving-average window in seconds (e.g. 0.05 for noisy dynamic/isovelocity)',
    )
    args = p.parse_args(argv)

    paths = list(args.files)
    if args.fuzzy_batch:
        paths.extend(FUZZY_BATCH)
    if not paths:
        p.error('Provide FILE.h5 paths or --fuzzy-batch')

    hold_only = not args.no_hold_only
    in_cal_only = not args.no_in_cal_only
    smooth_s = float(args.smooth or 0.0)
    if args.fuzzy_batch and smooth_s <= 0:
        smooth_s = 0.2
    for path in paths:
        if not os.path.isfile(path):
            print(f'SKIP (missing): {path}', file=sys.stderr)
            continue
        use_hold = hold_only
        if use_hold and 'dynamic' in os.path.basename(path).lower():
            use_hold = False
        out = plot_and_save(
            path,
            cal_mode=args.cal,
            hold_only=use_hold,
            in_cal_only=in_cal_only,
            smooth_window_s=smooth_s,
        )
        print('Saved', out)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
