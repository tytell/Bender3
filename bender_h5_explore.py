"""
Read experiment `.h5` files for the **H5 Data Explorer** (Streamlit): schema v2 (``02_TimeSeries``)
and legacy layouts (``Calibrated`` / ``NominalStimulus``).
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Tuple

import h5py
import numpy as np

from bender_functions import convert_to_curvature

# Matches :attr:`Bender.forcetorque` row order (Fx, Fy, Fz, Tx, Ty, Tz).
FT_ROW_LABELS: Tuple[str, ...] = (
    'Fx (N)',
    'Fy (N)',
    'Fz (N)',
    'Tx (N·m)',
    'Ty (N·m)',
    'Tz (N·m)',
)


def detect_h5_schema(path: str) -> str:
    """Return ``v2``, ``legacy``, or ``unknown``."""
    if not path or not os.path.isfile(path):
        return 'unknown'
    try:
        with h5py.File(path, 'r') as f:
            if '02_TimeSeries' in f:
                return 'v2'
            if 'Calibrated' in f and 'NominalStimulus' in f:
                return 'legacy'
    except OSError:
        return 'unknown'
    return 'unknown'


def list_v2_trials(path: str) -> List[str]:
    with h5py.File(path, 'r') as f:
        g = f.get('02_TimeSeries')
        if g is None:
            return []
        return sorted(k for k in g.keys() if str(k).startswith('trial_'))


def _as_float_1d(ds: h5py.Dataset) -> np.ndarray:
    return np.asarray(ds[...], dtype=np.float64).ravel()


def _trial_dclamp_mm(trial_group: h5py.Group) -> Optional[float]:
    for key in ('dclamp', 'test_segment_length_mm', 'dvert'):
        if key in trial_group.attrs:
            try:
                v = float(trial_group.attrs[key])
                if np.isfinite(v) and v > 0:
                    return v
            except (TypeError, ValueError):
                continue
    return None


def build_series_catalog_v2(path: str, trial: str) -> Tuple[Dict[str, np.ndarray], List[str]]:
    """
    Build human-readable series name → 1D ``float64`` arrays (original lengths; align in UI).

    Returns
    -------
    catalog, notes
    """
    notes: List[str] = []
    out: Dict[str, np.ndarray] = {}
    with h5py.File(path, 'r') as f:
        tg = f['02_TimeSeries'][trial]
        if 't' in tg:
            out['Time (s)'] = _as_float_1d(tg['t'])
        else:
            notes.append('No `t` dataset — time axis unavailable.')
        n_ref = len(out.get('Time (s)', np.array([])))

        if 'tnorm' in tg:
            out['Normalized time (tnorm)'] = _as_float_1d(tg['tnorm'])

        if 'angle_measured' in tg:
            out['Angle measured (deg)'] = _as_float_1d(tg['angle_measured'])
        if 'angle_cmd' in tg:
            out['Angle command (deg)'] = _as_float_1d(tg['angle_cmd'])

        dc = _trial_dclamp_mm(tg)
        if dc is not None and 'angle_measured' in tg:
            try:
                ang = _as_float_1d(tg['angle_measured'])
                kappa = convert_to_curvature(ang, 'angle', dclamp_mm=dc)
                out['Curvature κ (1/m) from angle'] = np.asarray(kappa, dtype=np.float64).ravel()
            except Exception as e:
                notes.append(f'Curvature from angle skipped: {e}')
        elif 'angle_measured' in tg:
            notes.append('Curvature κ not computed (no positive `dclamp` / `test_segment_length_mm` on trial).')

        for key, label in (
            ('primary_torque_corrected', 'Primary torque corrected (N·m)'),
            ('primary_torque_raw', 'Primary torque raw (N·m)'),
        ):
            if key in tg:
                out[label] = _as_float_1d(tg[key])

        if 'forcetorque' in tg:
            ft = np.asarray(tg['forcetorque'][...], dtype=np.float64)
            if ft.ndim == 2:
                if ft.shape[0] == 6:
                    for i, name in enumerate(FT_ROW_LABELS):
                        out[name] = ft[i, :].ravel()
                elif ft.shape[1] == 6:
                    for i, name in enumerate(FT_ROW_LABELS):
                        out[name] = ft[:, i].ravel()
                else:
                    notes.append(f'`forcetorque` shape {ft.shape} not 6×N or N×6 — skipped.')
            else:
                notes.append('`forcetorque` is not 2D — skipped.')

        if 'forcetorque_corrected' in tg:
            ftc = np.asarray(tg['forcetorque_corrected'][...], dtype=np.float64)
            if ftc.ndim == 2 and (ftc.shape[0] == 6 or ftc.shape[1] == 6):
                prefix = 'FT corrected — '
                if ftc.shape[0] == 6:
                    for i, name in enumerate(FT_ROW_LABELS):
                        out[prefix + name] = ftc[i, :].ravel()
                else:
                    for i, name in enumerate(FT_ROW_LABELS):
                        out[prefix + name] = ftc[:, i].ravel()

        # Drop empty
        out = {k: v for k, v in out.items() if v.size > 0}
        if n_ref == 0 and out:
            n_ref = min(len(v) for v in out.values())
            notes.append('Using shortest series length as reference (no `t`).')

    return out, notes


def build_series_catalog_legacy(path: str) -> Tuple[Dict[str, np.ndarray], List[str]]:
    notes: List[str] = []
    out: Dict[str, np.ndarray] = {}
    with h5py.File(path, 'r') as f:
        ns = f['NominalStimulus']
        cal = f['Calibrated']
        out['Time (s)'] = _as_float_1d(ns['t'])
        if 'tnorm' in ns:
            out['Normalized time (tnorm)'] = _as_float_1d(ns['tnorm'])

        torque_paths = (
            ('z Torque (N·m)', 'zTorque'),
            ('x Torque (N·m)', 'xTorque'),
            ('y Torque (N·m)', 'yTorque'),
        )
        for label, leaf in torque_paths:
            p = f'Calibrated/{leaf}'
            if p in f:
                out[label] = _as_float_1d(f[p])

        enc_candidates = ('Encoder', 'encoder')
        for leaf in enc_candidates:
            p = f'Calibrated/{leaf}'
            if p in f:
                out['Encoder / angle (deg)'] = _as_float_1d(f[p])
                break

        dclamp: Optional[float] = None
        for obj in (ns, cal, f):
            for key in ('dclamp', 'Dclamp', 'test_segment_length_mm'):
                if key in obj.attrs:
                    try:
                        v = float(obj.attrs[key])
                        if np.isfinite(v) and v > 0:
                            dclamp = v
                            break
                    except (TypeError, ValueError):
                        continue
            if dclamp is not None:
                break

        if dclamp and 'Encoder / angle (deg)' in out:
            try:
                ang = out['Encoder / angle (deg)']
                out['Curvature κ (1/m) from angle'] = np.asarray(
                    convert_to_curvature(ang, 'angle', dclamp_mm=dclamp), dtype=np.float64
                ).ravel()
            except Exception as e:
                notes.append(f'Legacy curvature from angle skipped: {e}')
        elif 'Encoder / angle (deg)' in out:
            notes.append('Curvature κ not computed (no dclamp in file attrs).')

        out = {k: v for k, v in out.items() if v.size > 0}

    return out, notes


def align_xy(
    catalog: Dict[str, np.ndarray], x_key: str, y_key: str
) -> Tuple[np.ndarray, np.ndarray, int]:
    """Trim X and Y to common length; returns ``x, y, n``."""
    if x_key not in catalog or y_key not in catalog:
        raise KeyError('Unknown X or Y series.')
    x = np.asarray(catalog[x_key], dtype=np.float64).ravel()
    y = np.asarray(catalog[y_key], dtype=np.float64).ravel()
    n = int(min(x.size, y.size))
    if n <= 0:
        return x[:0], y[:0], 0
    return x[:n], y[:n], n
