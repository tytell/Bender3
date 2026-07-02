"""
Read experiment `.h5` files for the **H5 Data Explorer** (Streamlit): schema v2 (``timeseries``),
legacy layouts (``Calibrated`` / optional ``NominalStimulus`` / ``RawInput``), and a **generic**
fallback that discovers 1D numeric series (and 6-channel FT blocks) anywhere in the file.
"""
from __future__ import annotations

import os
from typing import Any, Dict, List, Optional, Tuple, Union

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


def _top_level_key_ci(f: h5py.Group, wanted: Tuple[str, ...]) -> Optional[str]:
    """Resolve a top-level group name case-insensitively; return actual key or ``None``."""
    lower = {str(k).lower(): k for k in f.keys()}
    for w in wanted:
        if w.lower() in lower:
            return str(lower[w.lower()])
    return None


def _dataset_in_group_ci(g: h5py.Group, wanted: Tuple[str, ...]) -> Optional[str]:
    lower = {str(k).lower(): k for k in g.keys()}
    for w in wanted:
        if w.lower() in lower:
            return str(lower[w.lower()])
    return None


def _attrs_float(keys: Tuple[str, ...], *objs: Optional[Union[h5py.Group, h5py.File]]) -> Optional[float]:
    for obj in objs:
        if obj is None:
            continue
        for key in keys:
            if key in obj.attrs:
                try:
                    v = float(obj.attrs[key])
                    if np.isfinite(v) and v > 0:
                        return v
                except (TypeError, ValueError):
                    continue
    return None


def _guess_sample_rate_hz(
    f: h5py.File,
    cal_name: Optional[str],
    ns_name: Optional[str],
    ri_name: Optional[str],
) -> float:
    fs = _attrs_float(
        (
            'SampleFrequency',
            'sample_frequency',
            'sample_rate_hz',
            'SampleRate',
            'fs',
            'Fs',
            'DAQRate',
        ),
        f[ri_name] if ri_name else None,
        f[ns_name] if ns_name else None,
        f[cal_name] if cal_name else None,
        f,
    )
    return float(fs) if fs is not None else 1000.0


def _has_6_channel_ft(ds: h5py.Dataset) -> bool:
    if ds.ndim != 2:
        return False
    a, b = int(ds.shape[0]), int(ds.shape[1])
    return a == 6 or b == 6


def _legacy_group_has_torque_like(cal: h5py.Group) -> bool:
    for k in cal.keys():
        item = cal[k]
        if not isinstance(item, h5py.Dataset):
            continue
        kl = str(k).lower()
        if 'torque' in kl or kl in ('zt', 'tx', 'ty', 'tz', 'zt_torque'):
            return True
    return False


def _legacy_heuristic(f: h5py.File) -> bool:
    cal = _top_level_key_ci(f, ('Calibrated', 'CALIBRATED', 'calibrated'))
    ns = _top_level_key_ci(
        f, ('NominalStimulus', 'NOMINALSTIMULUS', 'Nominal_Stimulus', 'nominal_stimulus')
    )
    if cal and ns:
        return True
    if cal and _legacy_group_has_torque_like(f[cal]):
        return True
    ri = _top_level_key_ci(f, ('RawInput', 'RAWINPUT', 'rawinput'))
    if ri:
        g = f[ri]
        for ft_key in ('forcetransducer', 'Forcetransducer', 'force_torque', 'FT'):
            if ft_key in g:
                ds = g[ft_key]
                if isinstance(ds, h5py.Dataset) and _has_6_channel_ft(ds):
                    return True
    return False


def _h5_dataset_dtype_plottable(dt: np.dtype) -> bool:
    """True if elements can be coerced to float for plotting (excludes compound / opaque / strings)."""
    if getattr(dt, 'names', None):
        return False
    try:
        k = dt.kind
    except Exception:
        return False
    return k in 'fiub'


def classify_h5_dataset_plot_mode(ds: h5py.Dataset) -> str:
    """Return ``none``, ``1d``, or ``six`` for explorer / UI logic."""
    if ds.size < 2 or not _h5_dataset_dtype_plottable(ds.dtype):
        return 'none'
    if ds.ndim == 1:
        return '1d'
    if ds.ndim == 2:
        a, b = int(ds.shape[0]), int(ds.shape[1])
        if a == 1 or b == 1:
            return '1d'
        if a == 6 and b == 6:
            return 'none'
        if a == 6 or b == 6:
            return 'six'
    return 'none'


def _generic_count_plottable(f: h5py.File) -> int:
    """Rough count of datasets the generic catalog would expose (cap during walk)."""
    n = 0

    def visitor(_name: str, obj: h5py.Group | h5py.Dataset) -> None:
        nonlocal n
        if n >= 200:
            return
        if not isinstance(obj, h5py.Dataset):
            return
        if obj.size < 2:
            return
        if not _h5_dataset_dtype_plottable(obj.dtype):
            return
        mode = classify_h5_dataset_plot_mode(obj)
        if mode == '1d':
            n += 1
        elif mode == 'six':
            n += 6

    try:
        f.visititems(visitor)
    except Exception:
        pass
    return n


def detect_h5_schema(path: str) -> str:
    """Return ``v2``, ``legacy``, ``generic``, ``browse``, or ``unknown``."""
    if not path or not os.path.isfile(path):
        return 'unknown'
    try:
        with h5py.File(path, 'r') as f:
            if 'timeseries' in f:
                return 'v2'
            if _legacy_heuristic(f):
                return 'legacy'
            if _generic_count_plottable(f) >= 1:
                return 'generic'
            return 'browse'
    except OSError:
        return 'unknown'


def list_v2_trials(path: str) -> List[str]:
    with h5py.File(path, 'r') as f:
        g = f.get('timeseries')
        if g is None:
            return []
        return sorted(k for k in g.keys() if str(k).startswith('step_'))


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

    def _pick(group, *names):
        """Return the first dataset name present in ``group`` (canonical first, legacy fallback)."""
        for nm in names:
            if nm in group:
                return nm
        return None

    with h5py.File(path, 'r') as f:
        tg = f['timeseries'][trial]
        t_key = _pick(tg, 'time_second', 't')
        if t_key:
            out['Time (s)'] = _as_float_1d(tg[t_key])
        else:
            notes.append('No `time_second` dataset — time axis unavailable.')
        n_ref = len(out.get('Time (s)', np.array([])))

        tnorm_key = _pick(tg, 'time_normalized', 'tnorm')
        if tnorm_key:
            out['Normalized time'] = _as_float_1d(tg[tnorm_key])

        ang_meas_key = _pick(tg, 'angle_measured_degree', 'angle_measured')
        ang_cmd_key = _pick(tg, 'angle_commanded_degree', 'angle_cmd')
        if ang_meas_key:
            out['Angle measured (deg)'] = _as_float_1d(tg[ang_meas_key])
        if ang_cmd_key:
            out['Angle command (deg)'] = _as_float_1d(tg[ang_cmd_key])

        dc = _trial_dclamp_mm(tg)
        if dc is not None and ang_meas_key:
            try:
                ang = _as_float_1d(tg[ang_meas_key])
                kappa = convert_to_curvature(ang, 'angle', dclamp_mm=dc)
                out['Curvature κ (1/m) from angle'] = np.asarray(kappa, dtype=np.float64).ravel()
            except Exception as e:
                notes.append(f'Curvature from angle skipped: {e}')
        elif ang_meas_key:
            notes.append('Curvature κ not computed (no positive `dclamp` / `test_segment_length_mm` on trial).')

        # primary_torque_* is no longer written (derived → R/hub); read only for legacy files.
        for key, label in (
            ('primary_torque_corrected', 'Primary torque corrected (N·m)'),
            ('primary_torque_raw', 'Primary torque raw (N·m)'),
        ):
            if key in tg:
                out[label] = _as_float_1d(tg[key])

        ft_key = _pick(tg, 'forcetorque_raw', 'forcetorque')
        if ft_key:
            ft = np.asarray(tg[ft_key][...], dtype=np.float64)
            if ft.ndim == 2:
                if ft.shape[0] == 6:
                    for i, name in enumerate(FT_ROW_LABELS):
                        out[name] = ft[i, :].ravel()
                elif ft.shape[1] == 6:
                    for i, name in enumerate(FT_ROW_LABELS):
                        out[name] = ft[:, i].ravel()
                else:
                    notes.append(f'`{ft_key}` shape {ft.shape} not 6×N or N×6 — skipped.')
            else:
                notes.append(f'`{ft_key}` is not 2D — skipped.')

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
        cal_name = _top_level_key_ci(f, ('Calibrated', 'CALIBRATED', 'calibrated'))
        ns_name = _top_level_key_ci(
            f, ('NominalStimulus', 'NOMINALSTIMULUS', 'Nominal_Stimulus', 'nominal_stimulus')
        )
        ri_name = _top_level_key_ci(f, ('RawInput', 'RAWINPUT', 'rawinput'))

        cal = f[cal_name] if cal_name else None
        ns = f[ns_name] if ns_name else None
        ri = f[ri_name] if ri_name else None

        t_ds: Optional[h5py.Dataset] = None
        if ns is not None:
            tk = _dataset_in_group_ci(ns, ('t', 'time', 'Time', 'T'))
            if tk is not None:
                t_ds = ns[tk]
        if t_ds is None and cal is not None:
            tk = _dataset_in_group_ci(cal, ('t', 'time', 'Time'))
            if tk is not None:
                t_ds = cal[tk]
        if t_ds is not None:
            out['Time (s)'] = _as_float_1d(t_ds)

        if ns is not None:
            tnk = _dataset_in_group_ci(ns, ('tnorm', 'Tnorm', 't_norm'))
            if tnk is not None:
                out['Normalized time (tnorm)'] = _as_float_1d(ns[tnk])

        torque_specs = (
            ('z Torque (N·m)', ('zTorque', 'ZTorque', 'z_torque', 'zt_torque', 'TorqueZ', 'torque_z', 'Tz')),
            ('x Torque (N·m)', ('xTorque', 'XTorque', 'x_torque', 'TorqueX', 'torque_x', 'Tx')),
            ('y Torque (N·m)', ('yTorque', 'YTorque', 'y_torque', 'TorqueY', 'torque_y', 'Ty')),
        )
        if cal is not None:
            for label, leaves in torque_specs:
                leaf = _dataset_in_group_ci(cal, leaves)
                if leaf is not None:
                    out[label] = _as_float_1d(cal[leaf])

        enc_leaves = (
            'Encoder',
            'encoder',
            'Angle',
            'angle',
            'angle_measured',
            'Angle_measured',
        )
        if cal is not None:
            elk = _dataset_in_group_ci(cal, enc_leaves)
            if elk is not None:
                out['Encoder / angle (deg)'] = _as_float_1d(cal[elk])

        sono_specs = (
            (
                'Left sono / length',
                ('sono_left', 'Sono_left', 'left_sono', 'Left_Muscle_Length', 'left_muscle_length'),
            ),
            (
                'Right sono / length',
                ('sono_right', 'Sono_right', 'right_sono', 'Right_Muscle_Length', 'right_muscle_length'),
            ),
        )
        if cal is not None:
            for label, leaves in sono_specs:
                if label in out:
                    continue
                sk = _dataset_in_group_ci(cal, leaves)
                if sk is not None:
                    out[label] = _as_float_1d(cal[sk])
            for sub in ('Sonomicrometry', 'sonomicrometry', 'Sono'):
                if sub not in cal:
                    continue
                sg = cal[sub]
                if not isinstance(sg, h5py.Group):
                    continue
                for label, leaves in sono_specs:
                    if label in out:
                        continue
                    sk = _dataset_in_group_ci(sg, leaves)
                    if sk is not None:
                        out[label] = _as_float_1d(sg[sk])

        if ri is not None:
            ft_key: Optional[str] = None
            for fk in ('forcetransducer', 'Forcetransducer', 'force_torque', 'FT'):
                if fk in ri:
                    ft_key = fk
                    break
            if ft_key is not None:
                ft_ds = ri[ft_key]
                if isinstance(ft_ds, h5py.Dataset) and _has_6_channel_ft(ft_ds):
                    ft = np.asarray(ft_ds[...], dtype=np.float64)
                    if ft.shape[0] == 6:
                        for i, name in enumerate(FT_ROW_LABELS):
                            out['RawInput — ' + name] = ft[i, :].ravel()
                    else:
                        for i, name in enumerate(FT_ROW_LABELS):
                            out['RawInput — ' + name] = ft[:, i].ravel()

        dclamp: Optional[float] = None
        for obj in (ns, cal, f):
            if obj is None:
                continue
            v = _attrs_float(('dclamp', 'Dclamp', 'test_segment_length_mm'), obj)
            if v is not None:
                dclamp = v
                break

        if dclamp is not None and 'Encoder / angle (deg)' in out:
            try:
                ang = out['Encoder / angle (deg)']
                out['Curvature κ (1/m) from angle'] = np.asarray(
                    convert_to_curvature(ang, 'angle', dclamp_mm=dclamp), dtype=np.float64
                ).ravel()
            except Exception as e:
                notes.append(f'Legacy curvature from angle skipped: {e}')
        elif 'Encoder / angle (deg)' in out:
            notes.append('Curvature κ not computed (no dclamp in file attrs).')

        if 'Time (s)' not in out:
            lengths = [int(v.size) for v in out.values()]
            nmax = max(lengths) if lengths else 0
            if nmax > 0:
                fs = _guess_sample_rate_hz(f, cal_name, ns_name, ri_name)
                out['Time (s)'] = np.arange(nmax, dtype=np.float64) / fs
                notes.append(
                    'No time vector found — using t = arange(n)/'
                    f'{fs:g} Hz (from attrs if present, else default 1000).'
                )

        out = {k: v for k, v in out.items() if v.size > 0}

    if not out:
        gcat, gnotes = build_series_catalog_generic(path)
        return gcat, notes + gnotes

    return out, notes


def build_series_catalog_generic(path: str, max_series: int = 120) -> Tuple[Dict[str, np.ndarray], List[str]]:
    """Discover 1D numeric datasets (and 6×N force/torque blocks) anywhere in the file."""
    notes: List[str] = ['Generic layout: numeric series discovered by walking the HDF5 tree.']
    out: Dict[str, np.ndarray] = {}

    def visitor(name: str, obj: h5py.Group | h5py.Dataset) -> None:
        if len(out) >= max_series:
            return
        if not isinstance(obj, h5py.Dataset) or obj.size < 2:
            return
        if not _h5_dataset_dtype_plottable(obj.dtype):
            return
        label_base = name[1:] if name.startswith('/') else name
        pretty = label_base.replace('/', ' › ')
        mode = classify_h5_dataset_plot_mode(obj)
        if mode == '1d':
            out[pretty] = _as_float_1d(obj)
            return
        if mode != 'six':
            return
        arr = np.asarray(obj[...], dtype=np.float64)
        sh = arr.shape
        if sh[0] == 6:
            for i, lbl in enumerate(FT_ROW_LABELS):
                if len(out) >= max_series:
                    return
                out[f'{pretty} — {lbl}'] = arr[i, :].ravel()
        else:
            for i, lbl in enumerate(FT_ROW_LABELS):
                if len(out) >= max_series:
                    return
                out[f'{pretty} — {lbl}'] = arr[:, i].ravel()

    with h5py.File(path, 'r') as f:
        try:
            f.visititems(visitor)
        except Exception as e:
            notes.append(f'Walk truncated: {e}')

    out = {k: v for k, v in out.items() if v.size > 0}
    if not out:
        notes.append('No plottable 1D numeric series found.')
    return out, notes


def _validate_h5_rel_parts(rel: str) -> Tuple[str, ...]:
    rel = (rel or '').strip().strip('/')
    if not rel:
        return tuple()
    parts = tuple(p for p in rel.split('/') if p)
    if any(p in ('.', '..', '') for p in parts):
        raise ValueError('Invalid HDF5 path')
    return parts


def h5_join_internal_path(parent: str, name: str) -> str:
    p = (parent or '').strip().strip('/')
    n = str(name).strip('/')
    if not n:
        raise ValueError('Empty name')
    if p:
        return f'{p}/{n}'
    return n


def list_h5_group_children(path: str, group_path: str) -> List[dict]:
    """
    List direct children of ``group_path`` (use ``''`` for the file root).

    Each entry: ``name``, ``kind`` (``group`` | ``dataset``), ``shape`` (str),
    ``plot`` (``none`` | ``1d`` | ``six``).
    """
    try:
        parts = _validate_h5_rel_parts(group_path)
    except ValueError:
        return []
    out: List[dict] = []
    try:
        with h5py.File(path, 'r') as f:
            g: h5py.Group = f
            for p in parts:
                if p not in g:
                    return []
                nx = g[p]
                if not isinstance(nx, h5py.Group):
                    return []
                g = nx
            names = sorted(g.keys(), key=lambda x: str(x).lower())
            for name in names:
                item = g[name]
                if isinstance(item, h5py.Group):
                    out.append({'name': str(name), 'kind': 'group', 'shape': '', 'plot': 'none'})
                elif isinstance(item, h5py.Dataset):
                    shp = '×'.join(str(x) for x in item.shape)
                    plt = classify_h5_dataset_plot_mode(item)
                    out.append({'name': str(name), 'kind': 'dataset', 'shape': shp, 'plot': plt})
    except OSError:
        return []
    return out


def read_h5_series_1d(
    path: str,
    dataset_internal_path: str,
    channel: Optional[int] = None,
) -> np.ndarray:
    """
    Load a numeric series as 1D ``float64``. For 6×N or N×6 datasets, pass ``channel`` in ``0..5``
    (same ordering as :data:`FT_ROW_LABELS`).
    """
    parts = _validate_h5_rel_parts(dataset_internal_path)
    if not parts:
        raise ValueError('Empty dataset path')
    with h5py.File(path, 'r') as f:
        ds: Union[h5py.Group, h5py.Dataset, h5py.File] = f
        for i, p in enumerate(parts):
            if p not in ds:
                raise KeyError('/'.join(parts[: i + 1]))
            ds = ds[p]
        if not isinstance(ds, h5py.Dataset):
            raise TypeError('Path does not resolve to a dataset')
        mode = classify_h5_dataset_plot_mode(ds)
        if mode == 'none':
            raise ValueError(f'Unsupported dataset shape {ds.shape} / dtype {ds.dtype}')
        arr = np.asarray(ds[...], dtype=np.float64)
        if mode == '1d':
            return arr.ravel()
        if channel is None:
            raise ValueError('This dataset is 6-channel — choose channel 0–5 (e.g. 5 = Tz).')
        c = int(channel)
        if c < 0 or c > 5:
            raise ValueError('Channel must be 0–5')
        if arr.shape[0] == 6:
            return arr[c, :].ravel()
        return arr[:, c].ravel()


def _attr_value_display_and_kind(value: Any) -> Tuple[str, str]:
    """Return (text_for_text_input, kind) where kind is float|int|bool|str|raw."""
    if isinstance(value, bytes):
        try:
            return (value.decode('utf-8'), 'str')
        except Exception:
            return (repr(value), 'raw')
    if isinstance(value, str):
        return (value, 'str')
    if isinstance(value, (bool, np.bool_)):
        return ('true' if bool(value) else 'false', 'bool')
    if isinstance(value, (np.floating, float)):
        t = float(value)
        if np.isfinite(t) and abs(t - round(t)) < 1e-12 and abs(t) < 1e15:
            return (repr(t), 'float')
        return (repr(t), 'float')
    if isinstance(value, (np.integer, int)):
        return (str(int(value)), 'int')
    if isinstance(value, np.ndarray):
        if value.size == 1:
            return _attr_value_display_and_kind(value.flat[0])
        return (np.array2string(value, separator=', ', max_line_width=120, threshold=12), 'raw')
    return (str(value), 'raw')


def list_h5_attribute_rows(path: str, internal_path: str) -> List[dict]:
    """
    Attributes on the object at ``internal_path`` (``''`` = file root).

    Each row: ``name``, ``value_text``, ``kind``, ``editable`` (bool).
    """
    rows: List[dict] = []
    with h5py.File(path, 'r') as f:
        obj: Union[h5py.File, h5py.Group, h5py.Dataset] = f
        parts = _validate_h5_rel_parts(internal_path)
        for p in parts:
            if p not in obj:
                raise KeyError(f'No such path segment: {p!r}')
            obj = obj[p]
        for k in sorted(obj.attrs.keys(), key=lambda x: str(x).lower()):
            v = obj.attrs[k]
            text, kind = _attr_value_display_and_kind(v)
            rows.append(
                {
                    'name': str(k),
                    'value_text': text,
                    'kind': kind,
                    'editable': kind in ('float', 'int', 'bool', 'str'),
                }
            )
    return rows


def _parse_attr_text(text: str, kind: str) -> Any:
    s = (text or '').strip()
    if kind == 'str':
        return s
    if kind == 'float':
        return float(s)
    if kind == 'int':
        return int(float(s))
    if kind == 'bool':
        sl = s.lower()
        if sl in ('true', '1', 'yes', 'on'):
            return True
        if sl in ('false', '0', 'no', 'off'):
            return False
        raise ValueError('use true or false')
    raise ValueError(f'cannot parse kind {kind!r}')


def write_h5_user_attributes(
    path: str,
    internal_path: str,
    *,
    updates: Dict[str, Tuple[str, str]],
    delete_names: List[str],
    additions: Optional[List[Tuple[str, str, str]]] = None,
) -> List[str]:
    """
    Open ``path`` read/write and apply attribute changes.

    ``updates``: attr name -> (kind, new_text) for editable kinds only.
    ``additions``: optional list of (name, kind, text).
    """
    notes: List[str] = []
    additions = additions or []
    with h5py.File(path, 'r+') as f:
        obj: Union[h5py.File, h5py.Group, h5py.Dataset] = f
        parts = _validate_h5_rel_parts(internal_path)
        for p in parts:
            if p not in obj:
                raise KeyError(f'No such path segment: {p!r}')
            obj = obj[p]
        for dn in delete_names:
            if dn in obj.attrs:
                del obj.attrs[dn]
                notes.append(f'Deleted attribute {dn!r}')
        for name, (kind, text) in updates.items():
            if kind == 'raw':
                notes.append(f'Skipped {name!r} (unsupported type in UI)')
                continue
            try:
                val = _parse_attr_text(text, kind)
            except Exception as e:
                raise ValueError(f'{name!r}: {e}') from e
            obj.attrs[name] = val
            notes.append(f'Updated {name!r}')
        for name, kind, text in additions:
            name = str(name).strip()
            if not name:
                continue
            if name in obj.attrs:
                raise ValueError(f'Attribute {name!r} already exists')
            val = _parse_attr_text(text, kind)
            obj.attrs[name] = val
            notes.append(f'Added {name!r}')
    return notes


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
