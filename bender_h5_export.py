"""
Experiment **data file** (HDF5 / `.h5`) export and universal QC figure — same behavior as ``bender_run.ipynb`` cells.

Used by the Streamlit GUI and callable from scripts/notebooks to avoid duplicating logic.
"""
from __future__ import annotations

import importlib
import inspect
import json
import os
import sys
from datetime import datetime
from typing import Any, Dict, Optional, Tuple

import h5py
import numpy as np


def _is_non_data_object(value: Any) -> bool:
    """True for values that are not experiment data and must never be serialized into the file.

    Catches callables, file handles, and instances of non-builtin classes (e.g. the
    ``MasterLogger``). These were previously dumped into the file via ``str(value)``,
    storing opaque reprs instead of failing loudly / skipping cleanly.
    """
    if callable(value):
        return True
    if hasattr(value, 'read') and callable(getattr(value, 'read', None)):
        return True
    if isinstance(value, (str, bytes, bool, int, float, complex, list, tuple, dict, set)):
        return False
    if isinstance(value, (np.generic, np.ndarray)):
        return False
    return (getattr(type(value), '__module__', '') or '') not in ('builtins',)


def _store_h5_metadata_value(group: Any, key: str, value: Any, *, max_attr_array_size: int = 512) -> None:
    """Store one metadata value, fail-loud: never stringify arrays or opaque objects.

    * ``None`` and non-data objects (loggers, callables, file handles) are skipped, leaving a
      short ``<omitted ...>`` note attribute rather than ``str(obj)``.
    * Scalars (incl. strings) become attributes.
    * Arrays become datasets; arrays larger than ``max_attr_array_size`` are omitted with a
      shape note so the group never duplicates large buffers or overflows the 64 KB attr limit.
    """
    k = str(key)
    if value is None:
        return
    if _is_non_data_object(value):
        group.attrs[k] = f'<omitted_non_data_object type={type(value).__name__}>'
        return
    if isinstance(value, (str, bytes)):
        group.attrs[k] = value
        return
    try:
        arr = np.asarray(value)
    except Exception:
        group.attrs[k] = f'<omitted_unserializable type={type(value).__name__}>'
        return
    if arr.ndim == 0:
        try:
            group.attrs[k] = arr.item()
        except Exception:
            group.attrs[k] = str(value)[:512]
        return
    if arr.dtype == object:
        group.attrs[k] = f'<omitted_object_array shape={arr.shape}>'
        return
    if arr.size <= max_attr_array_size:
        try:
            group.create_dataset(k, data=arr)
        except Exception:
            # Small non-numeric arrays (e.g. unicode lists) that h5py can't store directly:
            # a compact string attribute is acceptable and readable for short metadata.
            group.attrs[k] = str(value)[:512]
    else:
        group.attrs[k] = f'<omitted_large_array shape={arr.shape}>'


def _store_config_field(group: Any, key: str, value: Any) -> None:
    """Store one **config-file** field into ``group`` for full provenance, additive only.

    Mapping (per approved Task 1 contract):

    * ``dict`` (e.g. ``units`` / ``unit_rules``) → JSON-string attribute written with an
      explicit variable-length string dtype so long dictionaries are never truncated.
    * ``str`` / ``bytes`` / scalars → attributes.
    * ``list`` / ``tuple`` / ``ndarray`` → datasets (string/object arrays are encoded to bytes
      so h5py can store them; if that fails, fall back to a JSON-string attribute).

    Caller must guard against name collisions with existing ``01_Metadata`` entries; this
    function never inspects or rewrites entries other than the one ``key`` it is asked to add.
    """
    k = str(key)
    vlen_str = h5py.special_dtype(vlen=str)
    if isinstance(value, dict):
        try:
            payload = json.dumps(value, default=str, sort_keys=True)
        except Exception:
            payload = str(value)
        group.attrs.create(k, payload, dtype=vlen_str)
        return
    if isinstance(value, (str, bytes)):
        group.attrs[k] = value
        return
    if isinstance(value, (bool, int, float, complex)) or isinstance(value, np.generic):
        group.attrs[k] = value
        return
    if isinstance(value, (list, tuple, np.ndarray)):
        arr = np.asarray(value)
        if arr.dtype.kind in ('U', 'O'):
            try:
                arr = arr.astype('S')
            except Exception:
                try:
                    group.attrs.create(k, json.dumps(list(value), default=str), dtype=vlen_str)
                except Exception:
                    group.attrs[k] = str(value)
                return
        try:
            group.create_dataset(k, data=arr)
        except Exception:
            group.attrs[k] = str(value)
        return
    group.attrs[k] = str(value)


def _iter_public_config_fields(config_module_name: str):
    """Yield ``(name, value)`` for every public, data-like attribute of the config module.

    The module is taken from ``sys.modules`` when already imported (it is, since ``Bender``
    imported it at construction) and re-imported otherwise. Modules, callables, and classes are
    skipped so only actual config values (scalars, strings, lists, dicts) are yielded.
    """
    name = str(config_module_name or '').strip()
    if not name:
        return
    cfg = sys.modules.get(name)
    if cfg is None:
        try:
            cfg = importlib.import_module(name)
        except Exception:
            return
    for attr in sorted(dir(cfg)):
        if attr.startswith('_'):
            continue
        try:
            value = getattr(cfg, attr)
        except Exception:
            continue
        if inspect.ismodule(value) or inspect.isclass(value) or callable(value):
            continue
        yield attr, value


def _read_existing_post_trial_notes(h5_path: str) -> str:
    """Return root ``post-trial notes`` attribute from an existing file, if any."""
    if not h5_path or not os.path.isfile(h5_path):
        return ''
    try:
        with h5py.File(h5_path, 'r') as old:
            raw = old.attrs.get('post-trial notes', '') or ''
            return str(raw).strip()
    except Exception:
        return ''


def unique_filepath(path: str) -> str:
    """
    If ``path`` is already an existing file, return the first unused name
    ``stem_001.ext``, ``stem_002.ext``, … in the same directory so exports never
    overwrite data on disk.
    """
    if not path:
        return path
    p = os.path.normpath(path)
    if not os.path.isfile(p):
        return p
    directory, filename = os.path.split(p)
    stem, ext = os.path.splitext(filename)
    n = 1
    while True:
        cand = os.path.join(directory, f'{stem}_{n:03d}{ext}')
        if not os.path.isfile(cand):
            return cand
        n += 1


def _append_post_trial_notes_to_file_notes(previous: str, addition: str) -> str:
    """
    Build the full ``post-trial notes`` string for the HDF5 root / metadata attrs.

    New non-empty text is appended as a dated block so repeated exports to the same path
    preserve earlier QC / rig comments.
    """
    add = (addition or '').strip()
    prev = (previous or '').strip()
    if not add:
        return prev
    if not prev:
        return add
    ts = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    return f'{prev}\n\n--- QC / post-experiment ({ts}) ---\n{add}'


def export_primary_h5(
    bender: Any,
    *,
    post_trial_notes: Optional[str] = None,
    outputfile: Optional[str] = None,
    append_post_trial_notes: bool = True,
) -> Dict[str, Any]:
    """
    Write schema-v2 H5 from ``bender.trial_records`` (or legacy flat buffers).

    Parameters
    ----------
    bender
        :class:`bender_functions.Bender` instance with ``outputfile`` set.
    post_trial_notes
        New text to store. If ``append_post_trial_notes`` is ``True``, only this
        string is **appended** (``None`` or empty = add nothing new, keep existing
        file note). If ``append_post_trial_notes`` is ``False``, ``None`` falls
        back to ``bender.post_trial_notes`` for the full stored string.
    outputfile
        Full path for the new file; default ``bender.outputfile``.
    append_post_trial_notes
        If ``True`` (default), read any existing ``post-trial notes`` from the
        file that will be written (or, when the chosen path is already taken and
        a ``_001`` / ``_002`` … name is used, from the conflicting path) and
        **append** this call's text as a new timestamped block. If ``False``,
        stored notes are replaced by this call's text only. Empty new text leaves
        prior notes unchanged when appending.

    Returns
    -------
    dict
        ``outputfile``, ``n_trials``, ``test_type``, ``schema_version``, ``message``,
        ``post_trial_notes`` (full string written to the file).
    """
    out_path = outputfile or getattr(bender, 'outputfile', None)
    if not out_path:
        raise ValueError("export_primary_h5 requires bender.outputfile or outputfile=...")

    final_path = unique_filepath(out_path)

    test_type = str(getattr(bender, 'test_type', 'unknown') or 'unknown')
    h5_schema_version = str(getattr(bender, 'h5_schema_version', '2.0'))
    if append_post_trial_notes:
        incoming = '' if post_trial_notes is None else str(post_trial_notes).strip()
        previous = ''
        if os.path.isfile(final_path):
            previous = _read_existing_post_trial_notes(final_path)
        if (
            not previous
            and os.path.normpath(final_path) != os.path.normpath(out_path)
            and os.path.isfile(out_path)
        ):
            previous = _read_existing_post_trial_notes(out_path)
        notes = _append_post_trial_notes_to_file_notes(previous, incoming)
    else:
        if post_trial_notes is not None:
            notes = str(post_trial_notes).strip()
        else:
            notes = str(getattr(bender, 'post_trial_notes', '') or '').strip()

    setattr(bender, 'post_trial_notes', notes)
    setattr(bender, 'outputfile', final_path)

    cal_file = str(getattr(bender, 'inertial_calibration_file', '') or '')
    use_cal = bool(getattr(bender, 'use_inertial_calibration', False))
    cal_available = bool(use_cal and cal_file and os.path.isfile(cal_file))

    trial_records = list(getattr(bender, 'trial_records', []) or [])
    if len(trial_records) == 0:
        trial_records = [{
            'test_type': test_type,
            'trial_index': 0,
            'cycle_index': 0,
            't': np.asarray(getattr(bender, 't', np.array([]))),
            'angle_cmd': np.asarray(getattr(bender, 'angle', np.array([]))),
            'anglevel_cmd': np.asarray(getattr(bender, 'anglevel', np.array([]))),
            'tnorm': np.asarray(getattr(bender, 'tnorm', np.array([]))),
            'S1stimcmd': np.asarray(getattr(bender, 'S1stimcmd', np.array([]))),
            'S2stimcmd': np.asarray(getattr(bender, 'S2stimcmd', np.array([]))),
            'aidata': np.asarray(getattr(bender, 'aidata', np.array([]))),
            'angle_measured': np.asarray(getattr(bender, 'angle_measured', np.array([]))),
            'forcetorque': np.asarray(getattr(bender, 'forcetorque', np.array([]))),
            'forcetorque_raw': np.asarray(getattr(bender, 'forcetorque_raw', np.array([]))),
            'forcetorque_corrected': np.asarray(getattr(bender, 'forcetorque_corrected', np.array([]))),
            'inertial_torque_system_primary': np.asarray(
                getattr(bender, 'inertial_torque_system_primary', np.array([]))
            ),
            'inertial_torque_specimen_primary': np.asarray(
                getattr(bender, 'inertial_torque_specimen_primary', np.array([]))
            ),
            'inertial_torque_total_primary': np.asarray(
                getattr(bender, 'inertial_torque_total_primary', np.array([]))
            ),
            'primary_torque_raw': np.asarray(getattr(bender, 'primary_torque_raw', np.array([]))),
            'primary_torque_corrected': np.asarray(getattr(bender, 'primary_torque_corrected', np.array([]))),
        }]

    # Store the file's own name, not the absolute path: portable across machines and avoids
    # leaking the rig's directory layout into the data file.
    filename_only = os.path.basename(final_path)
    # Experiment start time as ISO-8601. ``bender.timestamp`` is the compact run stamp
    # (``%Y%m%d-%H%M%S``); convert it (falling back to now()) without overwriting it.
    ts_raw = str(getattr(bender, 'timestamp', '') or '').strip()
    try:
        start_dt = datetime.strptime(ts_raw, '%Y%m%d-%H%M%S')
    except (ValueError, TypeError):
        start_dt = datetime.now()
    start_time_iso = start_dt.isoformat()

    with h5py.File(final_path, 'w') as f:
        f.attrs['schema_version'] = h5_schema_version
        f.attrs['test_type'] = test_type
        f.attrs['post-trial notes'] = notes
        f.attrs['filename'] = filename_only
        f.attrs['start_time_iso'] = start_time_iso

        g_meta = f.create_group('01_Metadata')
        g_meta.attrs['test_type'] = test_type
        g_meta.attrs['schema_version'] = h5_schema_version
        g_meta.attrs['post-trial notes'] = notes
        g_meta.attrs['filename'] = filename_only
        g_meta.attrs['start_time_iso'] = start_time_iso

        g_cal_link = g_meta.create_group('calibration_link')
        g_cal_link.attrs['use_inertial_calibration'] = bool(use_cal)
        g_cal_link.attrs['calibration_file'] = cal_file
        g_cal_link.attrs['calibration_available'] = bool(cal_available)

        g_ts = f.create_group('02_TimeSeries')
        manifest_rows = []

        for i, rec in enumerate(trial_records):
            tg = g_ts.create_group(f'trial_{i:04d}')
            tg.attrs['trial_index'] = int(rec.get('trial_index', i))
            tg.attrs['cycle_index'] = int(rec.get('cycle_index', i))
            rec_tt = str(rec.get('test_type', test_type))
            tg.attrs['test_type'] = rec_tt
            manifest = {
                'trial_name': f'trial_{i:04d}',
                'trial_index': int(rec.get('trial_index', i)),
                'cycle_index': int(rec.get('cycle_index', i)),
                'test_type': rec_tt,
            }

            series_keys = [
                't', 'angle_cmd', 'anglevel_cmd', 'tnorm', 'S1stimcmd', 'S2stimcmd',
                'aidata', 'angle_measured', 'forcetorque', 'forcetorque_raw', 'forcetorque_corrected',
                'inertial_torque_system_primary', 'inertial_torque_specimen_primary',
                'inertial_torque_total_primary',
                'primary_torque_raw', 'primary_torque_corrected',
                'cycle_index_by_sample', 'stim_type', 'stim_state', 'stim_side',
            ]

            t_arr = np.asarray(rec.get('t', np.array([]))).reshape(-1)
            s1_arr = np.asarray(rec.get('S1stimcmd', np.array([]))).reshape(-1)
            s2_arr = np.asarray(rec.get('S2stimcmd', np.array([]))).reshape(-1)
            cyc_arr = np.asarray(rec.get('cycle_index_by_sample', np.array([]))).reshape(-1)
            n = int(t_arr.size)
            if n <= 0:
                n = int(max(s1_arr.size, s2_arr.size, cyc_arr.size, 0))
            # Align cycle_index_by_sample to the trial's sample count so the written dataset and
            # the stim-state logic agree with the time base (pad unknown samples with -1, the
            # "not a numbered cycle" tag; truncate any overrun).
            if n > 0 and cyc_arr.size > 0 and cyc_arr.size != n:
                cyc_aligned = np.full(n, -1, dtype=int)
                m = min(n, cyc_arr.size)
                cyc_aligned[:m] = np.asarray(cyc_arr[:m]).astype(int, copy=False)
                cyc_arr = cyc_aligned
                rec['cycle_index_by_sample'] = cyc_aligned
            if n > 0:
                s1 = np.zeros(n, dtype=float)
                s2 = np.zeros(n, dtype=float)
                s1[:min(n, s1_arr.size)] = s1_arr[:min(n, s1_arr.size)]
                s2[:min(n, s2_arr.size)] = s2_arr[:min(n, s2_arr.size)]
                stim_on_mask = (np.abs(s1) > 1e-12) | (np.abs(s2) > 1e-12)
                stim_enabled = bool(rec.get('is_stim', getattr(bender, 'is_stim', False)))
                phase_state = np.full(n, 'passive', dtype='<U8')
                activity_state = np.full(n, 'passive', dtype='<U8')
                if stim_enabled:
                    if cyc_arr.size == n:
                        valid_cyc = np.isfinite(cyc_arr)
                        if np.any(valid_cyc):
                            for cyc in np.unique(cyc_arr[valid_cyc]):
                                m = valid_cyc & (cyc_arr == cyc)
                                if np.any(stim_on_mask[m]):
                                    phase_state[m] = 'off'
                                    phase_state[m & stim_on_mask] = 'on'
                                    activity_state[m] = 'active'
                        else:
                            if np.any(stim_on_mask):
                                phase_state[:] = 'off'
                                phase_state[stim_on_mask] = 'on'
                                activity_state[:] = 'active'
                    else:
                        if np.any(stim_on_mask):
                            phase_state[:] = 'off'
                            phase_state[stim_on_mask] = 'on'
                            activity_state[:] = 'active'
                side_state = np.full(n, 'none', dtype='<U8')
                left_on = np.abs(s1) > 1e-12
                right_on = np.abs(s2) > 1e-12
                side_state[left_on & ~right_on] = 'left'
                side_state[~left_on & right_on] = 'right'
                side_state[left_on & right_on] = 'both'
                rec['stim_type'] = np.asarray(activity_state, dtype='S8')
                rec['stim_state'] = np.asarray(phase_state, dtype='S8')
                rec['stim_side'] = np.asarray(side_state, dtype='S8')
                tg.attrs['stim_enabled'] = bool(stim_enabled)
                tg.attrs['stim_any_on'] = bool(np.any(stim_on_mask))
                manifest['stim_enabled'] = bool(stim_enabled)
                manifest['stim_any_on'] = bool(np.any(stim_on_mask))
            for key in series_keys:
                if key in rec and rec[key] is not None:
                    arr = np.asarray(rec[key])
                    if arr.size > 0:
                        tg.create_dataset(key, data=arr)

            for k, v in rec.items():
                if k in series_keys or k in ('trial_index', 'cycle_index', 'test_type'):
                    continue
                try:
                    arr = np.asarray(v)
                    if arr.ndim == 0:
                        vv = arr.item()
                        tg.attrs[str(k)] = vv
                        manifest[str(k)] = vv
                    elif arr.size == 1:
                        vv = arr.reshape(-1)[0].item()
                        tg.attrs[str(k)] = vv
                        manifest[str(k)] = vv
                except Exception:
                    tg.attrs[str(k)] = str(v)
                    manifest[str(k)] = str(v)
            manifest_rows.append(manifest)

        g_idx = g_meta.create_group('trial_index')
        g_idx.create_dataset(
            'trial_names',
            data=np.array([f'trial_{i:04d}' for i in range(len(trial_records))], dtype='S'),
        )
        g_idx.create_dataset(
            'trial_index',
            data=np.array([int(r.get('trial_index', i)) for i, r in enumerate(manifest_rows)], dtype=np.int64),
        )
        g_idx.create_dataset(
            'cycle_index',
            data=np.array([int(r.get('cycle_index', i)) for i, r in enumerate(manifest_rows)], dtype=np.int64),
        )
        g_idx.create_dataset(
            'test_type',
            data=np.array([str(r.get('test_type', test_type)) for r in manifest_rows], dtype='S'),
        )
        reserved = {'trial_name', 'trial_index', 'cycle_index', 'test_type'}
        all_keys = sorted({k for r in manifest_rows for k in r.keys() if k not in reserved})
        for cond_key in all_keys:
            col = [r.get(cond_key, np.nan) for r in manifest_rows]
            numeric_vals = []
            numeric_ok = True
            has_numeric = False
            for v in col:
                try:
                    fv = float(v)
                    numeric_vals.append(fv)
                    if np.isfinite(fv):
                        has_numeric = True
                except Exception:
                    numeric_ok = False
                    break
            if numeric_ok and has_numeric:
                g_idx.create_dataset(cond_key, data=np.array(numeric_vals, dtype=float))
            else:
                svals = ['' if (v is None) else str(v) for v in col]
                if any(len(s) > 0 for s in svals):
                    g_idx.create_dataset(cond_key, data=np.array(svals, dtype='S'))

        g_meta.attrs['n_trials'] = int(len(trial_records))

        prof = getattr(bender, 'inertial_calibration_profile', None)
        if isinstance(prof, dict):
            g_ic = g_meta.create_group('inertial_calibration_profile')
            for k, v in prof.items():
                g_ic.attrs[str(k)] = v

        h5p = dict(getattr(bender, 'h5_protocol_metadata', {}) or {})
        g_proto = g_meta.create_group('protocol_metadata')
        for k, v in h5p.items():
            # Protocol arrays (e.g. frequency/amplitude lists) stay full datasets regardless of
            # length; opaque objects are skipped rather than stringified.
            _store_h5_metadata_value(g_proto, k, v, max_attr_array_size=2**31)

        g_settings = g_meta.create_group('bender_settings')
        skip_keys = {
            'aidata', 'forcetorque', 'angle', 'anglevel', 'tnorm', 't', 'angledata',
            'S1stimcmd', 'S2stimcmd', 'trial_records',
            # Not data: the MasterLogger object must not be serialized into the file.
            'master_logger',
            # Absolute output path is machine-specific; the portable filename is stored on
            # the root and 01_Metadata as ``filename`` instead.
            'outputfile',
        }
        for k, v in bender.__dict__.items():
            if k.startswith('_') or k in skip_keys or v is None:
                continue
            _store_h5_metadata_value(g_settings, k, v)

        # Full config-file provenance: write every public config attribute directly into
        # 01_Metadata under its EXACT config field name. Additive only — any name that already
        # exists as an attribute, dataset, or subgroup of 01_Metadata is left untouched so the
        # established HDF5 / R-pipeline contract is preserved.
        for cfg_key, cfg_val in _iter_public_config_fields(getattr(bender, 'config_name', '')):
            if cfg_val is None:
                continue
            if cfg_key in g_meta.attrs or cfg_key in g_meta:
                continue
            _store_config_field(g_meta, cfg_key, cfg_val)

    msg = f'EXPORT FINISHED (schema={h5_schema_version}, test_type={test_type}, n_trials={len(trial_records)})'
    return {
        'outputfile': final_path,
        'n_trials': len(trial_records),
        'test_type': test_type,
        'schema_version': h5_schema_version,
        'message': msg,
        'post_trial_notes': notes,
    }


def _concat_trial_records(records):
    """Concatenate per-step trial records into one continuous timeline for the QC figure.

    Each step's ``t`` restarts near 0; we offset each segment by the cumulative duration so the
    saved figure shows every discrete step (matching the live preview), instead of collapsing to a
    single step. Missing per-segment arrays are filled with NaN so trace lengths stay aligned.
    """
    t_parts = []
    keys_1d = ['angle_cmd', 'anglevel_cmd', 'angle_measured', 'S1stimcmd', 'S2stimcmd',
               'inertial_torque_total_primary']
    parts_1d = {k: [] for k in keys_1d}
    ft_raw_parts = []
    ft_corr_parts = []
    t_offset = 0.0
    for r in records:
        tr = np.asarray(r.get('t', np.array([])), dtype=float).reshape(-1)
        if tr.size == 0:
            continue
        n = tr.size
        tr0 = tr - tr[0] + t_offset
        t_parts.append(tr0)
        for k in keys_1d:
            a = np.asarray(r.get(k, np.array([])), dtype=float).reshape(-1)
            parts_1d[k].append(a if a.size == n else np.full(n, np.nan, dtype=float))
        raw = r.get('forcetorque_raw', r.get('forcetorque', None))
        raw = np.asarray(raw, dtype=float) if raw is not None else np.array([])
        ft_raw_parts.append(raw if (raw.ndim == 2 and raw.shape[1] == n and raw.shape[0] >= 6)
                            else np.full((6, n), np.nan, dtype=float))
        corr = r.get('forcetorque_corrected', None)
        corr = np.asarray(corr, dtype=float) if corr is not None else np.array([])
        ft_corr_parts.append(corr if (corr.ndim == 2 and corr.shape[1] == n and corr.shape[0] >= 6)
                             else np.full((6, n), np.nan, dtype=float))
        dt = (tr[-1] - tr[0]) / max(1, n - 1)
        t_offset = float(tr0[-1] + dt)
    combined: dict = {'t': np.concatenate(t_parts) if t_parts else np.array([])}
    for k in keys_1d:
        combined[k] = np.concatenate(parts_1d[k]) if parts_1d[k] else np.array([])
    combined['forcetorque_raw'] = np.concatenate(ft_raw_parts, axis=1) if ft_raw_parts else np.array([])
    combined['forcetorque_corrected'] = (
        np.concatenate(ft_corr_parts, axis=1) if ft_corr_parts else np.array([])
    )
    return combined


def build_universal_qc_figure(bender: Any, qc_trial_index=None):
    """
    Build the multi-row QC Plotly figure.

    With multiple per-step ``trial_records`` and no specific ``qc_trial_index`` (``None`` or the
    string ``'all'``), all steps are concatenated into one timeline so the saved figure shows the
    discrete steps just like the preview. Pass an integer to plot a single step in detail.
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    trial_records = list(getattr(bender, 'trial_records', []) or [])
    if len(trial_records) == 0:
        trial_records = [{
            't': np.asarray(getattr(bender, 't', np.array([]))),
            'angle_cmd': np.asarray(getattr(bender, 'angle', np.array([]))),
            'anglevel_cmd': np.asarray(getattr(bender, 'anglevel', np.array([]))),
            'angle_measured': np.asarray(getattr(bender, 'angle_measured', np.array([]))),
            'S1stimcmd': np.asarray(getattr(bender, 'S1stimcmd', np.array([]))),
            'S2stimcmd': np.asarray(getattr(bender, 'S2stimcmd', np.array([]))),
            'forcetorque': np.asarray(getattr(bender, 'forcetorque', np.array([]))),
            'forcetorque_raw': np.asarray(getattr(bender, 'forcetorque_raw', np.array([]))),
            'forcetorque_corrected': np.asarray(getattr(bender, 'forcetorque_corrected', np.array([]))),
            'inertial_torque_total_primary': np.asarray(
                getattr(bender, 'inertial_torque_total_primary', np.array([]))
            ),
        }]

    combine_all = (qc_trial_index is None) or (
        isinstance(qc_trial_index, str) and qc_trial_index.lower() == 'all'
    )
    if combine_all and len(trial_records) > 1:
        rec = _concat_trial_records(trial_records)
        qc_trial_index = 'all'
    else:
        if not isinstance(qc_trial_index, (int, np.integer)):
            qc_trial_index = len(trial_records) - 1
        qc_trial_index = int(max(0, min(int(qc_trial_index), len(trial_records) - 1)))
        rec = trial_records[qc_trial_index]

    t = np.asarray(rec.get('t', np.array([])))
    angle_cmd = np.asarray(rec.get('angle_cmd', np.array([])))
    anglevel_cmd = np.asarray(rec.get('anglevel_cmd', np.array([])))
    angle_meas = np.asarray(rec.get('angle_measured', np.array([])))
    S1 = np.asarray(rec.get('S1stimcmd', np.array([])))
    S2 = np.asarray(rec.get('S2stimcmd', np.array([])))

    ft_raw = np.asarray(rec.get('forcetorque_raw', rec.get('forcetorque', np.array([]))))
    ft_corr = np.asarray(rec.get('forcetorque_corrected', np.array([])))
    inertial_total = np.asarray(rec.get('inertial_torque_total_primary', np.array([])))

    axis_key = str(
        getattr(bender, 'primary_bending_axis', getattr(bender, 'bending_axis_sensor', 'zTorque'))
    ).lower().strip()
    axis_norm = {
        'x': 'xTorque', 'xtorque': 'xTorque',
        'y': 'yTorque', 'ytorque': 'yTorque',
        'z': 'zTorque', 'ztorque': 'zTorque',
    }.get(axis_key, 'zTorque')
    axis_to_idx = {'xTorque': 3, 'yTorque': 4, 'zTorque': 5}
    primary_idx = axis_to_idx.get(axis_norm, 5)
    all_torque_axes = ['xTorque', 'yTorque', 'zTorque']
    off_axes = [ax for ax in all_torque_axes if ax != axis_norm]

    def _torque_row(arr, idx):
        a = np.asarray(arr)
        if a.ndim == 2 and a.shape[0] >= 6:
            return a[idx, :]
        return np.array([])

    primary_raw = _torque_row(ft_raw, primary_idx)
    primary_corr = _torque_row(ft_corr, primary_idx)
    off1 = _torque_row(ft_raw, axis_to_idx[off_axes[0]]) if len(off_axes) > 0 else np.array([])
    off2 = _torque_row(ft_raw, axis_to_idx[off_axes[1]]) if len(off_axes) > 1 else np.array([])

    fig = make_subplots(
        rows=5, cols=1, shared_xaxes=True, vertical_spacing=0.03,
        specs=[[{'secondary_y': True}], [{}], [{}], [{}], [{}]],
    )

    if t.size > 0 and angle_cmd.size == t.size:
        fig.add_trace(
            go.Scatter(x=t, y=angle_cmd, mode='lines', name='angle_cmd', line=dict(dash='dash', color='black')),
            row=1, col=1, secondary_y=False,
        )
    if t.size > 0 and angle_meas.size == t.size:
        fig.add_trace(
            go.Scatter(x=t, y=angle_meas, mode='lines', name='angle_measured', line=dict(color='royalblue')),
            row=1, col=1, secondary_y=False,
        )
    if t.size > 0 and anglevel_cmd.size == t.size:
        fig.add_trace(
            go.Scatter(x=t, y=anglevel_cmd, mode='lines', name='anglevel_cmd', line=dict(color='orange', width=1)),
            row=1, col=1, secondary_y=True,
        )

    if t.size > 0 and primary_raw.size == t.size:
        fig.add_trace(
            go.Scatter(x=t, y=primary_raw, mode='lines', name=f'{axis_norm} raw', line=dict(color='firebrick')),
            row=2, col=1,
        )
    if t.size > 0 and primary_corr.size == t.size:
        fig.add_trace(
            go.Scatter(x=t, y=primary_corr, mode='lines', name=f'{axis_norm} corrected', line=dict(color='seagreen')),
            row=2, col=1,
        )
    if t.size > 0 and inertial_total.size == t.size and np.any(np.isfinite(inertial_total)):
        fig.add_trace(
            go.Scatter(
                x=t, y=inertial_total, mode='lines',
                name='inertial torque (total, primary)',
                line=dict(color='goldenrod', dash='dot'),
            ),
            row=2, col=1,
        )

    if t.size > 0 and off1.size == t.size:
        fig.add_trace(
            go.Scatter(x=t, y=off1, mode='lines', name=f'{off_axes[0]} raw', line=dict(color='darkorange')),
            row=3, col=1,
        )
    if t.size > 0 and off2.size == t.size:
        fig.add_trace(
            go.Scatter(x=t, y=off2, mode='lines', name=f'{off_axes[1]} raw', line=dict(color='purple')),
            row=4, col=1,
        )

    if t.size > 0 and S1.size == t.size:
        fig.add_trace(go.Scatter(x=t, y=S1, mode='lines', name='S1 stim', line=dict(color='teal')), row=5, col=1)
    if t.size > 0 and S2.size == t.size:
        fig.add_trace(
            go.Scatter(x=t, y=S2, mode='lines', name='S2 stim', line=dict(color='gray', dash='dot')),
            row=5, col=1,
        )

    proc = str(getattr(bender, 'test_type', 'unknown'))
    angle_title = 'Angle (deg)'
    if proc in ('isometric', 'isovelocity'):
        parts = bender.strain_yaxis_title_pct().split(' — ', 1)
        angle_title = 'Angle (deg) — ' + (parts[1] if len(parts) > 1 else parts[0])
    fig.update_yaxes(title_text=angle_title, row=1, col=1, secondary_y=False)
    fig.update_yaxes(title_text='Velocity (deg/s)', row=1, col=1, secondary_y=True, showgrid=False)
    fig.update_yaxes(title_text=f'{axis_norm} (N-m)', row=2, col=1)
    fig.update_yaxes(title_text=f'{off_axes[0]} (N-m)', row=3, col=1)
    fig.update_yaxes(title_text=f'{off_axes[1]} (N-m)', row=4, col=1)
    fig.update_yaxes(title_text='Stim (V)', row=5, col=1)
    fig.update_xaxes(title_text='Time (s)', row=5, col=1)
    qc_cap = f'QC: {proc} | ' + ('all steps' if qc_trial_index == 'all' else f'trial {qc_trial_index}')
    if proc in ('isometric', 'isovelocity'):
        qc_cap = qc_cap + "<br><sup style='font-size:11px'>" + bender.strain_geometry_plot_context() + '</sup>'
    fig.update_layout(height=1400, width=1100, title_text=qc_cap, showlegend=True, hovermode='x unified')
    return fig, qc_trial_index


def save_universal_qc_figure(
    bender: Any,
    qc_trial_index: Optional[int] = None,
    *,
    base_path: Optional[str] = None,
) -> Tuple[str, Optional[str]]:
    """
    Write QC figure to PNG (kaleido) or HTML fallback.

    Returns
    -------
    (png_or_html_path, secondary_path_or_none)
    """
    fig, idx = build_universal_qc_figure(bender, qc_trial_index=qc_trial_index)
    proc = str(getattr(bender, 'test_type', 'unknown'))
    h5p = str(getattr(bender, 'outputfile', 'bender_output.h5'))
    idx_tag = 'all' if idx == 'all' else f'trial{int(idx):03d}'
    if base_path:
        base = base_path
    elif h5p.lower().endswith('.h5'):
        base = h5p[:-3] + f'_{proc}_{idx_tag}_qc'
    else:
        base = h5p + f'_{proc}_{idx_tag}_qc'
    png_path = unique_filepath(f'{base}.png')
    html_path = unique_filepath(f'{base}.html')
    try:
        fig.write_image(png_path, engine='kaleido', scale=2)
        return png_path, None
    except Exception:
        fig.write_html(html_path)
        return html_path, None
