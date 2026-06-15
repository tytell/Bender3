"""
Read Bender-exported HDF5 files for Streamlit custom plots (time series from ``02_TimeSeries``).
"""
from __future__ import annotations

import ast
from typing import Any, Dict, List, Optional, Tuple

import h5py
import numpy as np


def _coerce_channel_name_list(val: Any) -> List[str]:
    """Turn stored HDF5 attr/dataset value into a list of channel name strings."""
    if val is None:
        return []
    if isinstance(val, bytes):
        val = val.decode('utf-8', errors='replace')
    if isinstance(val, str):
        s = val.strip()
        if s.startswith('[') and s.endswith(']'):
            try:
                parsed = ast.literal_eval(s)
                if isinstance(parsed, (list, tuple)):
                    return [str(x) for x in parsed]
            except (SyntaxError, ValueError, TypeError):
                pass
        return [s] if s else []
    arr = np.asarray(val)
    if arr.ndim == 0:
        return _coerce_channel_name_list(arr.item())
    if arr.dtype.kind in 'SU':
        out = []
        for x in arr.ravel().tolist():
            if isinstance(x, bytes):
                out.append(x.decode('utf-8', errors='replace'))
            else:
                out.extend(_coerce_channel_name_list(str(x)))
        return out
    return []


def _decode_attr(val: Any) -> str:
    if val is None:
        return ''
    if isinstance(val, bytes):
        return val.decode('utf-8', errors='replace')
    return str(val)


def _read_input_channel_names(f: h5py.File) -> List[str]:
    # Canonical (Phase 0): instrumentation list at 01_Metadata/daq_instrumentation.
    try:
        meta = f['01_Metadata']
        if 'daq_instrumentation' in meta:
            return _coerce_channel_name_list(meta['daq_instrumentation'][()])
        v = meta.attrs.get('daq_instrumentation', None)
        if v is not None:
            return _coerce_channel_name_list(v)
    except Exception:
        pass
    # Legacy fallback: pre-Phase-0 files wrote 01_Metadata/bender_settings/input_channel_names.
    try:
        gs = f['01_Metadata/bender_settings']
        if 'input_channel_names' in gs:
            return _coerce_channel_name_list(gs['input_channel_names'][()])
        v = gs.attrs.get('input_channel_names', None)
        if v is not None:
            return _coerce_channel_name_list(v)
    except Exception:
        pass
    return []


def h5_custom_plot_summary(path: str) -> Dict[str, Any]:
    """
    Return file-level metadata for the custom plot UI.

    Keys: ok, error (if not ok), test_type, schema_version, trial_names, channel_names,
    primary_bending_axis, n_trials.
    """
    try:
        with h5py.File(path, 'r') as f:
            tt = _decode_attr(f.attrs.get('test_type', ''))
            schema = _decode_attr(f.attrs.get('schema_version', ''))
            if '02_TimeSeries' not in f:
                return {
                    'ok': False,
                    'error': 'This file has no `02_TimeSeries` group (not a Bender experiment export?).',
                    'test_type': tt,
                    'schema_version': schema,
                    'trial_names': [],
                    'channel_names': [],
                    'primary_bending_axis': '',
                    'n_trials': 0,
                }
            g_ts = f['02_TimeSeries']
            trial_names = sorted(k for k in g_ts.keys() if str(k).startswith('trial_'))
            ch = _read_input_channel_names(f)
            paxis = ''
            try:
                meta = f['01_Metadata']
                paxis = _decode_attr(meta.attrs.get('daq_primary_bending_axis', ''))
                if not paxis and 'bender_settings' in meta:
                    # Legacy fallback for pre-Phase-0 files.
                    paxis = _decode_attr(meta['bender_settings'].attrs.get('primary_bending_axis', ''))
            except Exception:
                pass
            return {
                'ok': True,
                'error': None,
                'test_type': tt,
                'schema_version': schema,
                'trial_names': trial_names,
                'channel_names': ch,
                'primary_bending_axis': paxis,
                'n_trials': len(trial_names),
            }
    except OSError as e:
        return {
            'ok': False,
            'error': str(e),
            'test_type': '',
            'schema_version': '',
            'trial_names': [],
            'channel_names': [],
            'primary_bending_axis': '',
            'n_trials': 0,
        }


def _numeric_1d_plottable(arr: np.ndarray) -> bool:
    if arr.ndim != 1 or arr.size < 2:
        return False
    return np.issubdtype(arr.dtype, np.number) or np.issubdtype(arr.dtype, np.floating)


def _numeric_row_plottable(row: np.ndarray) -> bool:
    v = np.asarray(row).reshape(-1)
    if v.size < 2:
        return False
    return np.issubdtype(v.dtype, np.number) or np.issubdtype(v.dtype, np.floating)


def list_h5_plot_variables(path: str, trial_name: str, channel_names: Optional[List[str]] = None) -> List[Dict[str, Any]]:
    """
    Build selectable series for X/Y.

    Each item: ``id`` (stable string), ``label`` (UI), ``trial_test_type`` (optional).
    Ids for matrix rows: ``{dataset}@row{i}`` when time runs along columns, else ``{dataset}@col{j}``.
    """
    channel_names = list(channel_names or [])
    out: List[Dict[str, Any]] = []
    with h5py.File(path, 'r') as f:
        tg = f['02_TimeSeries'][trial_name]
        trial_tt = _decode_attr(tg.attrs.get('test_type', ''))
        for key in sorted(tg.keys()):
            ds = tg[key]
            if not isinstance(ds, h5py.Dataset):
                continue
            try:
                arr = np.asarray(ds[()])
            except Exception:
                continue
            if arr.ndim == 1 and _numeric_1d_plottable(arr):
                out.append(
                    {
                        'id': str(key),
                        'label': f'{key}  (n={arr.size})',
                        'trial_test_type': trial_tt,
                    }
                )
            elif arr.ndim == 2:
                r, c = arr.shape
                if c >= r:
                    for i in range(r):
                        row = arr[i, :]
                        if not _numeric_row_plottable(row):
                            continue
                        nm = channel_names[i] if i < len(channel_names) else f'row {i}'
                        out.append(
                            {
                                'id': f'{key}@row{i}',
                                'label': f'{key} / {nm}  (n={row.size})',
                                'trial_test_type': trial_tt,
                            }
                        )
                else:
                    for j in range(c):
                        col = arr[:, j]
                        if not _numeric_row_plottable(col):
                            continue
                        nm = channel_names[j] if j < len(channel_names) else f'col {j}'
                        out.append(
                            {
                                'id': f'{key}@col{j}',
                                'label': f'{key} / {nm}  (n={col.size})',
                                'trial_test_type': trial_tt,
                            }
                        )
    return out


def read_h5_series(path: str, trial_name: str, var_id: str) -> np.ndarray:
    """Load a 1-D float vector for plotting."""
    with h5py.File(path, 'r') as f:
        tg = f['02_TimeSeries'][trial_name]
        if '@row' in var_id:
            name, _, idx_s = var_id.partition('@row')
            i = int(idx_s)
            arr = np.asarray(tg[name][()])
            if arr.ndim != 2:
                raise ValueError(f'{name} is not 2-D')
            return np.asarray(arr[i, :], dtype=float).reshape(-1)
        if '@col' in var_id:
            name, _, idx_s = var_id.partition('@col')
            j = int(idx_s)
            arr = np.asarray(tg[name][()])
            if arr.ndim != 2:
                raise ValueError(f'{name} is not 2-D')
            return np.asarray(arr[:, j], dtype=float).reshape(-1)
        arr = np.asarray(tg[var_id][()])
        return np.asarray(arr, dtype=float).reshape(-1)


def align_xy(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    n = int(min(x.size, y.size))
    if n <= 0:
        return x[:0], y[:0]
    return x[:n].astype(float), y[:n].astype(float)
