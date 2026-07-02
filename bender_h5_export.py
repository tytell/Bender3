"""
Experiment **data file** (HDF5 / `.h5`) export and universal QC figure — same behavior as ``bender_run.ipynb`` cells.

Used by the Streamlit GUI and callable from scripts/notebooks to avoid duplicating logic.
"""
from __future__ import annotations

import json
import logging
import os
from datetime import datetime
from typing import Any, Dict, Optional, Tuple

import h5py
import numpy as np

from bender_routing_spec import BENDER_ROUTING, EXCLUDED


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


def _store_metadata_canonical(group: Any, key: str, value: Any) -> bool:
    """Store one **canonical metadata** value, fail-loud and lossless.

    Used by the ledger-driven exporter to write a routed Bender attribute (or an unrouted
    one, into ``99_Unrouted``) under its canonical schema key. Mapping:

    * ``None`` / non-data objects (loggers, callables, file handles) → skipped.
    * ``dict`` → JSON-string attribute (variable-length str dtype; never truncated).
    * ``str`` / ``bytes`` / scalars → attributes.
    * ``list`` / ``tuple`` / ``ndarray`` → datasets at full length (string/object arrays are
      encoded to bytes; if that fails, a JSON-string attribute is used as a fallback).

    Returns ``True`` if anything was written, ``False`` if the value was skipped.
    """
    k = str(key)
    if value is None or _is_non_data_object(value):
        return False
    vlen_str = h5py.special_dtype(vlen=str)
    if isinstance(value, dict):
        try:
            payload = json.dumps(value, default=str, sort_keys=True)
        except Exception:
            payload = str(value)
        group.attrs.create(k, payload, dtype=vlen_str)
        return True
    if isinstance(value, (str, bytes)):
        group.attrs[k] = value
        return True
    if isinstance(value, (bool, int, float, complex)) or isinstance(value, np.generic):
        group.attrs[k] = value
        return True
    if isinstance(value, (list, tuple, np.ndarray)):
        arr = np.asarray(value)
        if arr.size == 0:
            return False
        if arr.dtype.kind in ('U', 'O'):
            try:
                arr = arr.astype('S')
            except Exception:
                try:
                    group.attrs.create(k, json.dumps(list(value), default=str), dtype=vlen_str)
                except Exception:
                    group.attrs[k] = str(value)
                return True
        try:
            group.create_dataset(k, data=arr)
        except Exception:
            group.attrs[k] = str(value)
        return True
    group.attrs[k] = str(value)
    return True


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

    # Prefix all sim-mode files so they sort separately and the analysis pipeline can quarantine
    # them by name alone (belt-and-suspenders alongside metadata.attrs['session_simulated']).
    # Idempotent: never double-prefixes if the caller has already named the file sim_*. The GUI
    # auto-namer already composes the sim_ name (so its disk scan and preview match this output);
    # this branch still covers manual-override, legacy, and direct (non-GUI) export paths.
    if bool(getattr(bender, 'session_simulated', False)):
        _d, _base = os.path.split(str(out_path))
        if not _base.startswith('sim_'):
            out_path = os.path.join(_d, 'sim_' + _base) if _d else 'sim_' + _base

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
            'primary_torque_raw': np.asarray(getattr(bender, 'primary_torque_raw', np.array([]))),
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
        # Root holds only schema_version (file-level provenance, D6). All other formerly root-level
        # attrs (test_type, post-trial notes, filename, start_time_iso) are written into metadata:
        # test_type -> protocol_type via ledger, post_trial_notes -> note_posthoc via ledger,
        # filename -> metadata/filename directly, start_time_iso -> metadata/session_date directly.
        f.attrs['schema_version'] = h5_schema_version

        g_meta = f.create_group('metadata')
        # Canonical scalar metadata (protocol_type, session_*, schema_version, note_posthoc,
        # starting_angle, etc.) is written by the ledger-driven loop near the end of this block
        # from BENDER_ROUTING. Portable file name and session date are kept here as direct writes
        # (filename is not a schema field; session_date is computed from bender.timestamp here).
        g_meta.attrs['filename'] = filename_only
        g_meta.attrs['session_date'] = start_time_iso

        # calibration_link subgroup removed (D3, D12): flatten to a single provenance attr.
        # use_inertial_calibration and calibration_available dropped -- no correction runs at export.
        # Always written: schema says "empty string if none" (D12).
        g_meta.attrs['calibration_inertia_file'] = cal_file
        # calibration_inertia_matrix: empty lookup table until apparatus-calibration workflow
        # is implemented (D12, Deferred flag 4). Present-but-empty -- omission not permitted.
        g_meta.create_dataset('calibration_inertia_matrix', data=np.empty((0,), dtype=np.float64))

        # calibration_forcetorque_matrix (D11, M2b): embed the ATI 6x6 VALUES ONLY when a REAL
        # calibration matrix is loaded. loadCalibration() in bender_functions falls back to
        # np.eye(6) when the .cal file is missing (sim / uncalibrated bench). An identity matrix
        # is NOT a calibration; embedding it would masquerade uncalibrated data as calibrated.
        # Refuse it: leave calibration_forcetorque_matrix ABSENT so an absent matrix unambiguously
        # means "not calibrated" (and derived/forcetorque_calibrated is likewise skipped below).
        # 'calibration' is bypassed in the ledger (special_metadata) so it is written only here.
        _cal_mat = getattr(bender, 'calibration', None)
        _cal_is_real = False
        if _cal_mat is not None:
            _cal_arr = np.asarray(_cal_mat, dtype=float)
            if _cal_arr.shape == (6, 6) and not np.array_equal(_cal_arr, np.eye(6)):
                _cal_is_real = True
                g_meta.create_dataset('calibration_forcetorque_matrix', data=_cal_arr)

        g_ts = f.create_group('timeseries')
        manifest_rows = []

        # Canonical per-trial timeseries dataset names (PI decision: full spelled-out rename).
        # timeseries holds ONLY the immutable raw ADC buffer ``aidata`` (D11, M2b). The F/T rows
        # live in aidata rows 0-5 (identity is declared in metadata/daq_ai_channel_map); the
        # calibrated F/T copy is NOT written to timeseries -- it is derived in R (authoritative
        # hub ``derived/``) and, for live RA inspection only, into this file's ``derived/`` group
        # when a REAL calibration matrix is present (see the derived/ block below).
        ts_rename = {
            't': 'time_second',
            'tnorm': 'time_normalized',
            'angle_cmd': 'angle_commanded_degree',
            'angle_measured': 'angle_measured_degree',
            'anglevel_cmd': 'angular_velocity_commanded_degree_per_second',
            'S1stimcmd': 'stim_channel1_command_volt',
            'S2stimcmd': 'stim_channel2_command_volt',
            'aidata': 'aidata',
            'stim_type': 'stim_type',
            'stim_state': 'stim_state',
            'stim_side': 'stim_side',
            'sweep_instantaneous_freq': 'instantaneous_frequency_hertz',
        }
        # Per-sample index family (converted from former per-trial scalar attrs). cycle_index is
        # written for every trial; the rest are written only where the source rec field exists
        # (step protocols carry step/sequence/block; dynamic/frequency_sweep carry only cycles).
        ts_index_int = ('step_index', 'sequence_index', 'block_index')
        ts_index_str = ('block_direction', 'block_stim_sides')
        # trial_index is dropped (redundant with filename + sequence_index, PI decision).
        ts_converted = set(ts_index_int) | set(ts_index_str) | {
            'cycle_index', 'cycle_index_by_sample', 'trial_index',
        }
        # Dropped from the raw file's timeseries (D11, M2b): ``forcetorque`` and ``forcetorque_raw``
        # (both are CALIBRATED F/T copies -- forcetorque_raw is a copy=True of the decoded
        # self.forcetorque, bender_functions L1114 -- so neither belongs in the raw-voltage archive;
        # aidata rows 0-5 are the sole archived F/T source) and ``primary_torque_raw`` (derived;
        # recomputed in R -- belongs in the research-hub ``derived`` group, not the raw file).
        ts_drop = {'forcetorque', 'forcetorque_raw', 'primary_torque_raw'}

        # Branch on daq_collection_type (Steps 2/3).
        # 'continuous' with exactly 1 record (single_finite, dynamic/sweep): write channel
        # datasets flat under timeseries — no per-step subgroup.
        # 'segmented' (isometric/isovelocity/FL/FV): one subgroup per step named step_{step_index:03d}
        # (one-based, 3-digit, step_ prefix); step_NNN suffix == rec['step_index'].
        # The old trial_{i:04d} branch is removed: all protocols are either single_finite
        # (continuous, 1 record) or segmented_finite (per-step run() loop) after Step 4.
        _daq_ct = str(getattr(bender, 'daq_collection_type', '') or '')
        _is_flat = (_daq_ct == 'continuous' and len(trial_records) <= 1)
        _is_segmented = (_daq_ct == 'segmented')

        for i, rec in enumerate(trial_records):
            if _is_flat:
                tg = g_ts
                _tname = 'timeseries'
            elif _is_segmented:
                _si = int(rec.get('step_index', i + 1))
                _tname = f'step_{_si:03d}'
                tg = g_ts.create_group(_tname)
            else:
                raise ValueError(
                    f'export_primary_h5: unrecognized daq_collection_type={_daq_ct!r}. '
                    'Expected "continuous" (single record) or "segmented".'
                )
            rec_tt = str(rec.get('test_type', test_type))
            tg.attrs['test_type'] = rec_tt
            manifest = {
                'trial_name': _tname,
                'test_type': rec_tt,
            }

            # RAW torque only. Inertial correction (system + specimen) is NOT applied or saved by
            # the acquisition GUI: empirical apparatus inertia is fit from an empty calibration run
            # post-hoc in R, and specimen inertia is analytic from geometry. The parameters needed
            # for that post-hoc correction (the use_theoretical_inertial_correction flag, specimen
            # MOI, and the inertial_calibration_profile I_est/bias) are stored in metadata, not
            # baked into a corrected time series here. So forcetorque_corrected,
            # primary_torque_corrected, and the inertial_torque_* traces are intentionally omitted.
            # Source rec keys we know how to canonicalize (written below under ts_rename names).
            series_keys = list(ts_rename.keys())

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

            # Per-sample index arrays (PI decision: convert from per-trial scalar attrs). Within a
            # trial group each step is one record, so the step/sequence/block fields are constant
            # across the trial's samples; they are broadcast to length n so the timeline stays
            # self-describing when trials are later concatenated. cycle_index is the per-sample
            # cycle number (dynamic/frequency_sweep) or -1 for step protocols (no numbered cycles).
            if n > 0:
                if cyc_arr.size == n:
                    cycle_idx_ps = np.asarray(cyc_arr, dtype=np.int64)
                else:
                    cycle_idx_ps = np.full(n, -1, dtype=np.int64)
                tg.create_dataset('cycle_index', data=cycle_idx_ps)

                for key in ts_index_int:
                    if key in rec and rec[key] is not None:
                        try:
                            val = int(np.asarray(rec[key]).reshape(-1)[0])
                        except Exception:
                            continue
                        tg.create_dataset(key, data=np.full(n, val, dtype=np.int64))
                for key in ts_index_str:
                    if key in rec and rec[key] is not None:
                        raw = np.asarray(rec[key]).reshape(-1)[0]
                        sval = raw.decode('utf-8') if isinstance(raw, (bytes, bytearray)) else str(raw)
                        tg.create_dataset(key, data=np.array([sval] * n, dtype='S16'))

            for key in series_keys:
                if key in rec and rec[key] is not None:
                    arr = np.asarray(rec[key])
                    if arr.size > 0:
                        tg.create_dataset(ts_rename[key], data=arr)

            for k, v in rec.items():
                if k in series_keys or k in ts_drop or k in ts_converted or k in ('test_type',):
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

        g_meta.attrs['n_trials'] = int(len(trial_records))

        # protocol_motor_positive_bend_direction / protocol_sensor_positive_bend_direction --
        # flat trial-level metadata attrs (constant within a trial; apparatus + mount choice).
        # Source: positive_motor_direction config constant ('left' or 'right').
        # The encoder (angle sensor) is co-mounted on the gearbox output shaft and shares the
        # same sign convention as the motor, so both fields derive from the same config attr.
        _pos_motor = str(getattr(bender, 'positive_motor_direction', '') or '').strip().lower()
        if _pos_motor not in ('left', 'right'):
            _pos_motor = 'none'
        g_meta.attrs['protocol_motor_positive_bend_direction'] = _pos_motor
        g_meta.attrs['protocol_sensor_positive_bend_direction'] = _pos_motor

        # index_step_* -- flat parallel arrays (one element per step), replacing the trial_index
        # subgroup (M2). Explicit name map: (source key in entry dict, hdf5 dataset name, dtype).
        # dtype 'float'    : missing key or None -> NaN; written only when at least one finite value.
        # dtype 'int'      : missing key or None -> 0; always written.
        # dtype 'int_1based': 0-based source int -> write as 1-based float; absent -> NaN.
        # dtype 'str'      : missing key or None -> ''; written only when at least one non-empty.
        _INDEX_STEP_MAP = [
            # Shared (all protocols)
            ('step_index',                       'index_step_step_number',                         'int'),
            ('wall_clock_start',                 'index_step_wall_clock_start',                    'str'),
            ('duration_second',                  'index_step_duration_second',                     'float'),
            ('rest_before_second',               'index_step_rest_before_second',                  'float'),
            ('operating_point',                  'index_step_operating_point',                     'float'),
            ('operating_point_units',            'index_step_operating_point_units',               'str'),
            ('recruitment',                      'index_step_recruitment',                         'str'),
            # None -> NaN when bilateral (unilateral returns an int, bilateral returns None)
            ('unilateral_posture_lateral_index', 'index_step_unilateral_posture_lateral_index',    'float'),
            ('stim_t0',                          'index_step_stim_t0_second',                      'float'),
            ('stim_t1',                          'index_step_stim_t1_second',                      'float'),
            ('t_pre_baseline_start',             'index_step_t_pre_baseline_start_second',         'float'),
            ('t_pre_baseline_end',               'index_step_t_pre_baseline_end_second',           'float'),
            ('t_active_start',                   'index_step_t_active_start_second',               'float'),
            ('t_active_end',                     'index_step_t_active_end_second',                 'float'),
            ('t_post_baseline_start',            'index_step_t_post_baseline_start_second',        'float'),
            ('t_post_baseline_end',              'index_step_t_post_baseline_end_second',          'float'),
            # Motor-position reference (_record_motor_position_reference).
            # index_step_cumulative_commanded_motor_microstep: Teknic motor-shaft MICROSTEPS
            # (1600/rev on motor shaft; gear ratio 5:1 to specimen frame = 8000/output-shaft rev).
            # NOT output-shaft steps or degrees (D9 rename from cumulative_commanded_steps).
            ('cumulative_commanded_steps',       'index_step_cumulative_commanded_motor_microstep', 'float'),
            ('encoder_cumulative_deg',           'index_step_encoder_cumulative_degree',           'float'),
            # Isometric-only (absent in isovelocity / single_finite -> NaN / '')
            ('target_deg',                       'index_step_target_angle_degree',                 'float'),
            ('ramp_from_deg',                    'index_step_ramp_from_angle_degree',              'float'),
            # Isovelocity-only (absent in isometric / single_finite -> NaN / '')
            ('velocity_deg_s',                   'index_step_velocity_degree_per_second',          'float'),
            ('theta_start_deg',                  'index_step_theta_start_degree',                  'float'),
            ('pre_stim_t0',                      'index_step_pre_stim_t0_second',                  'float'),
            ('pre_stim_t1',                      'index_step_pre_stim_t1_second',                  'float'),
            ('iso_t0',                           'index_step_iso_t0_second',                       'float'),
            ('iso_t1',                           'index_step_iso_t1_second',                       'float'),
            # None when F/T cal absent or no stim-active samples (sim mode -> NaN)
            ('mean_xforce_stim',                 'index_step_mean_xforce_stim_newton',             'float'),
            ('guard_triggered',                  'index_step_guard_triggered',                     'int'),
            ('guard_angle_deg',                  'index_step_guard_angle_degree',                  'float'),
            # Absent from entry dict when mirror mode inactive (not set at all -> NaN)
            ('velocity_seg1_deg_s',              'index_step_velocity_seg1_degree_per_second',     'float'),
            ('velocity_seg2_deg_s',              'index_step_velocity_seg2_degree_per_second',     'float'),
            # Block metadata: absent when use_block_stim=False; block_index is 0-based in entry
            # dict -> written 1-based here per schema convention.
            ('block_index',                      'index_step_block_number',                        'int_1based'),
            ('block_direction',                  'index_step_block_direction',                     'str'),
            ('block_stim_sides',                 'index_step_block_stim_sides',                    'str'),
            ('left_stim_voltage',                'index_step_stim_voltage_left_volt',              'float'),
            ('right_stim_voltage',               'index_step_stim_voltage_right_volt',             'float'),
        ]
        if trial_records:
            for _src, _dst, _dtype in _INDEX_STEP_MAP:
                if _dtype in ('float', 'int_1based'):
                    _col = []
                    for _j, _r in enumerate(trial_records):
                        v = _r.get(_src)
                        if _src == 'step_index' and v is None:
                            v = _j + 1
                        if _dtype == 'int_1based':
                            try:
                                _col.append(float(int(v) + 1) if v is not None else float('nan'))
                            except (TypeError, ValueError):
                                _col.append(float('nan'))
                        else:
                            try:
                                _col.append(float(v) if v is not None else float('nan'))
                            except (TypeError, ValueError):
                                _col.append(float('nan'))
                    if any(np.isfinite(x) for x in _col):
                        g_meta.create_dataset(_dst, data=np.array(_col, dtype=float))
                elif _dtype == 'int':
                    _col = []
                    for _j, _r in enumerate(trial_records):
                        v = _r.get(_src)
                        if _src == 'step_index' and v is None:
                            v = _j + 1
                        try:
                            _col.append(int(v) if v is not None else 0)
                        except (TypeError, ValueError):
                            _col.append(0)
                    g_meta.create_dataset(_dst, data=np.array(_col, dtype=np.int64))
                else:  # 'str'
                    _col = []
                    for _r in trial_records:
                        v = _r.get(_src)
                        if v is None:
                            _col.append('')
                        elif isinstance(v, (bytes, bytearray)):
                            _col.append(v.decode('utf-8', errors='replace'))
                        else:
                            _col.append(str(v))
                    if any(len(s) > 0 for s in _col):
                        g_meta.create_dataset(_dst, data=np.array(_col, dtype='S'))

        # index_cycle_* -- per-cycle design grid, single_finite ONLY (D10, M2c). Parallel arrays,
        # one element per bending cycle (length C); R joins them to samples via the per-sample
        # cycle_index stream (already in timeseries, section 3b-i). The source per-cycle arrays are
        # built by organize_cycles() and live on the Bender (EXCLUDED from the ledger), so they are
        # written here directly. segmented_finite has no per-cycle design grid, so this whole block
        # is gated on _is_flat; frequency_sweep (also single_finite) may not populate these arrays,
        # so each dataset is guarded on a present, non-empty source of consistent length.
        # Does NOT add or touch index_step_* columns (D4 intact).
        if _is_flat:
            _freq_by_cycle = getattr(bender, 'freq_by_cycle', None)
            _cyc = np.asarray(_freq_by_cycle, dtype=float).reshape(-1) if _freq_by_cycle is not None else np.array([])
            _C = int(_cyc.size)
            if _C > 0:
                def _cycle_col(attr, dtype=float):
                    """Per-cycle source array as length-C dtype array, or None if absent/ragged."""
                    v = getattr(bender, attr, None)
                    if v is None:
                        return None
                    a = np.asarray(v).reshape(-1)
                    if a.size != _C:
                        return None
                    return a.astype(dtype)

                g_meta.create_dataset('index_cycle_index', data=np.arange(1, _C + 1, dtype=np.int64))
                g_meta.create_dataset('index_cycle_frequency_hertz', data=_cyc.astype(float))
                _amp_deg = _cycle_col('amp_by_cycle', float)
                if _amp_deg is not None:
                    g_meta.create_dataset('index_cycle_motor_amplitude_degree', data=_amp_deg)
                _active = _cycle_col('is_stim_cycle', bool)
                if _active is not None:
                    g_meta.create_dataset('index_cycle_active', data=_active)
                _duty = _cycle_col('duty_by_cycle', float)
                if _duty is not None:
                    g_meta.create_dataset('index_cycle_activation_duty', data=_duty)
                _phase = _cycle_col('phase_by_cycle', float)
                if _phase is not None:
                    g_meta.create_dataset('index_cycle_activation_phase', data=_phase)

                # operating_point + units: user's native drive amplitude per cycle (mode-approved
                # mapping, PI sign-off 2026-07-02). motor amplitude (deg) is written above regardless
                # so R always has degrees plus the native operating point.
                _mode = str(getattr(bender, 'all_amps_mode', '') or '').strip().lower()
                _strain = _cycle_col('strain_by_cycle', float)
                _op_val = None
                _op_units = None
                if _mode == 'angle':
                    _op_val, _op_units = _amp_deg, 'degree'
                elif _mode == 'strain' and _strain is not None:
                    _op_val, _op_units = _strain, 'strain'
                elif _mode == 'strain_pct' and _strain is not None:
                    _op_val, _op_units = _strain * 100.0, 'percent'
                elif _mode == 'curvature' and _strain is not None:
                    # curvature = strain / lever_arm; lever_arm = (xsec_width/2 - depth)/1000 (m).
                    try:
                        _xw = float(getattr(bender, 'xsec_width', 0.0) or 0.0)
                        _md = float(getattr(bender, 'target_muscle_depth_mm', 0.0) or 0.0)
                        _lever_m = (_xw / 2.0 - _md) / 1000.0
                        if np.isfinite(_lever_m) and _lever_m > 0:
                            _op_val, _op_units = _strain / _lever_m, 'per_meter'
                    except (TypeError, ValueError):
                        _op_val = None
                if _op_val is None:
                    # Unknown/fallback: expose the realized motor amplitude in degrees.
                    _op_val, _op_units = _amp_deg, 'degree'
                if _op_val is not None:
                    g_meta.create_dataset('index_cycle_operating_point', data=np.asarray(_op_val, dtype=float))
                    # units as a [C] parallel array (schema section 4: index_cycle_operating_point_units
                    # is [C] str), matching the index_step_operating_point_units pattern.
                    g_meta.create_dataset(
                        'index_cycle_operating_point_units',
                        data=np.array([_op_units] * _C, dtype='S'),
                    )

                # protocol_ per-cycle scalars (D10). protocol_cycles_per_step is already written by
                # the ledger (cycles_per_step route); protocol_end_cycle_count by the ledger too.
                g_meta.attrs['protocol_cycle_count'] = _C
                _active_arr = _cycle_col('is_stim_cycle', bool)
                if _active_arr is not None and bool(np.any(_active_arr)):
                    # 0-based index of the first activated cycle (schema documents the deliberate
                    # 0-based/1-based asymmetry vs index_cycle_index; do not normalize).
                    _start = int(np.argmax(_active_arr))
                else:
                    _start = -1
                g_meta.attrs['protocol_activation_start_cycle'] = _start

        # step_manifest: per-step timing and operating-point index (schema §4).
        # Scaffold: step_index (1-based), duration_second from t array, rest_before_second
        # from rec field (0 until Steps 3/4 populate wall_clock timing). operating_point
        # fields added for segmented protocols where the value is already in the record.
        daq_ct = str(getattr(bender, 'daq_collection_type', '') or '')
        _ai_rate = float(getattr(bender, 'daq_ai_sample_rate_hz', 1000.0) or 1000.0)
        step_manifest_rows = []
        for _j, _rec in enumerate(trial_records):
            _t_arr = np.asarray(_rec.get('t', []), dtype=float).reshape(-1)
            _dur = float(_t_arr.size) / max(_ai_rate, 1.0) if _t_arr.size > 0 else None
            _row = {
                'step_index': int(_rec.get('step_index', _j + 1)),
                'wall_clock_start': _rec.get('wall_clock_start', None),
                'duration_second': _dur,
                'rest_before_second': float(_rec.get('rest_before_second', 0.0)),
            }
            if daq_ct == 'segmented':
                _op = _rec.get('velocity_deg_s') if 'velocity_deg_s' in _rec else _rec.get('target_deg', None)
                if _op is not None:
                    _row['operating_point'] = float(_op)
                _op_units = _rec.get('operating_point_units', None)
                if _op_units is not None:
                    _row['operating_point_units'] = str(_op_units)
            step_manifest_rows.append(_row)
        _step_manifest_json = json.dumps(step_manifest_rows, default=str)
        g_meta.attrs.create('step_manifest', _step_manifest_json, dtype=h5py.special_dtype(vlen=str))

        # inertial_calibration_profile -> three flat calibration_inertia_* attrs (D12, Point 2).
        # I_est is a MOI in N*m/(deg/s^2); multiply by (180/pi)*1e9 to convert to g*mm^2
        # so it sits in the same units as the geometry-derived specimen MOI.
        # bias_est is already N*m (no conversion). axis_sensor is categorical.
        _prof = getattr(bender, 'inertial_calibration_profile', None)
        if isinstance(_prof, dict):
            _i_est = _prof.get('I_est', None)
            if _i_est is not None:
                g_meta.attrs['calibration_inertia_apparatus_moi_gram_millimeter_squared'] = (
                    float(_i_est) * (180.0 / np.pi) * 1e9
                )
            _bias = _prof.get('bias_est', None)
            if _bias is not None:
                g_meta.attrs['calibration_inertia_bias_newton_meter'] = float(_bias)
            _axis = _prof.get('axis_sensor', None)
            if _axis is not None:
                g_meta.attrs['calibration_inertia_axis_sensor'] = str(_axis)

        # frustum_inputs dict -> four flat measurement_specimen_inertia_frustum_* attrs (v2.7 amended).
        # set_frustum_inertial_model stores values in the dict and does NOT set the four attrs
        # individually, so the ledger pass misses them. This block reads from the dict directly.
        # The ledger still fires for any attr set independently via update_metadata.
        _finp = getattr(bender, 'frustum_inputs', None)
        if isinstance(_finp, dict):
            for _fsrc, _fdst in (
                ('height_mm',         'measurement_specimen_inertia_frustum_height_millimeter'),
                ('width_mm',          'measurement_specimen_inertia_frustum_width_millimeter'),
                ('length_mm',         'measurement_specimen_inertia_frustum_length_millimeter'),
                ('density_g_per_mm3', 'measurement_specimen_inertia_frustum_density_gram_per_cubic_millimeter'),
            ):
                _fv = _finp.get(_fsrc, None)
                if _fv is not None:
                    g_meta.attrs[_fdst] = float(_fv)

        # calibration_inertia_apparatus_aor_to_clamp_millimeter -- many-to-one from two sources.
        # Primary: specimen_profile_clamp_offset_mm (set by live GUI path).
        # Fallback: frustum_inputs['clamp_offset_mm'] (set by frustum path).
        # OMIT entirely if neither source resolves to a positive value; never default-write 0.
        _aor_val = None
        _aor_primary = getattr(bender, 'specimen_profile_clamp_offset_mm', None)
        if _aor_primary is not None and float(_aor_primary) > 0:
            _aor_val = float(_aor_primary)
        if _aor_val is None and isinstance(_finp, dict):
            _aor_fb = _finp.get('clamp_offset_mm', None)
            if _aor_fb is not None and float(_aor_fb) > 0:
                _aor_val = float(_aor_fb)
        if _aor_val is not None:
            g_meta.attrs['calibration_inertia_apparatus_aor_to_clamp_millimeter'] = _aor_val

        # protocol_metadata subgroup REMOVED (PI decision, post-Phase-0 audit). Its attr-mirror keys
        # are written canonically by the ledger pass below; block_sequence is routed as
        # protocol_block_sequence (JSON). Run-computed provenance now routed canonically too
        # (Step 2): protocol_acquisition_mode, protocol_guard_triggered/_angle_degree,
        # inertial_specimen_from_geometry, inertial_system_from_profile,
        # daq_motor_positive_bend_lateral_index. The remaining dict-only summary keys
        # (n_trials/protocol/motion_test_type/frequency_hz/curvature_1_per_m/*movedur*/
        # simulation_model/theoretical_i_total_system/step_order) are dropped as redundant.

        # === Ledger-driven canonical metadata (Phase 1) ===================
        # Replaces BOTH the old ``bender_settings`` catch-all (raw bender.__dict__ dump) AND the
        # verbatim config-provenance writer. Every public Bender attribute is routed via
        # bender_routing_spec.BENDER_ROUTING to its canonical schema key:
        #   * tier == 'metadata'   -> written here into metadata under route['key']
        #   * tier == 'timeseries' -> already written per-trial above (skip here)
        #   * in EXCLUDED / underscore / None -> skipped
        #   * unmapped (neither routed nor excluded) -> preserved under 99_Unrouted + loud warning
        # Canonical-only output, no dual-write. Many-to-one routes write the FIRST present source
        # (subsequent sources for the same key are skipped) to avoid dataset collisions.
        # Pass A -- write canonical metadata. Iterate the LEDGER (not bender.__dict__) so that
        # @property-backed routed attributes (e.g. daq_ai_sample_rate_hz) are captured; instance
        # __dict__ never contains properties. Many-to-one routes write the FIRST present source
        # (ledger insertion order) to avoid dataset collisions.
        # 'calibration' is written ONLY by the real-matrix-gated block above (D11, M2b) -- bypass
        # the ledger so it is never emitted as the np.eye(6) identity fallback.
        special_metadata: set = {'calibration', 'h5_schema_version', 'inertial_calibration_file'}
        written_keys: set = set()
        for name, route in BENDER_ROUTING.items():
            if route.get('tier') != 'metadata' or name in special_metadata:
                continue
            key = route['key']
            if key in written_keys:
                continue
            value = getattr(bender, name, None)
            if value is None:
                continue
            if name == 'input_channels':
                # daq_ai_channel_map: zip channels with their identities -> index:identity map.
                names = list(getattr(bender, 'input_channel_names', []) or [])
                chans = list(value)
                cmap = {str(ix): f'{ch}:{nm}' for ix, (ch, nm) in enumerate(zip(chans, names))}
                if _store_metadata_canonical(g_meta, key, cmap):
                    written_keys.add(key)
            elif name == 'block_sequence':
                # protocol_block_sequence: list-of-dicts block plan -> JSON string (PI decision).
                try:
                    payload = json.dumps(value, default=str)
                except Exception:
                    payload = str(value)
                g_meta.attrs.create(key, payload, dtype=h5py.special_dtype(vlen=str))
                written_keys.add(key)
            elif _store_metadata_canonical(g_meta, key, value):
                written_keys.add(key)

        # note_bench fallback: schema marks this as a required field. The ledger writes it
        # when stim_clamp_notices is non-empty; emit an empty string when absent or empty.
        if 'note_bench' not in written_keys:
            g_meta.attrs['note_bench'] = ''

        # Pass B -- unmapped detection over the instance __dict__: any public attribute that is
        # neither routed nor excluded is preserved (lossless) under metadata/99_Unrouted (D7) with
        # a loud warning.
        unrouted: Dict[str, Any] = {}
        for name, value in list(vars(bender).items()):
            if name.startswith('_') or value is None:
                continue
            if name in BENDER_ROUTING or name in EXCLUDED:
                continue
            unrouted[name] = value
        if unrouted:
            g_un = g_meta.create_group('99_Unrouted')
            for name, value in unrouted.items():
                _store_metadata_canonical(g_un, name, value)
            logging.warning(
                'export_primary_h5: %d unrouted Bender attribute(s) preserved under metadata/99_Unrouted '
                '(add each to BENDER_ROUTING or EXCLUDED in bender_routing_spec.py): %s',
                len(unrouted), sorted(unrouted),
            )

        # Print surviving root-attr set for PI eyeball (D6). After M1 only schema_version
        # should remain at root; all other formerly root attrs live in metadata.
        print('[info] export_primary_h5: surviving root attrs = ' + str(list(f.attrs.keys())))

        # === RA-inspection derived/ group (D11 amendment, Point 4) ==============
        # Written into the raw file for live signal checking during/after a run.
        # R ignores this group entirely and re-derives from raw + embedded cal VALUES.
        # Unit discipline: I_est is N*m/(deg/s^2); specimen MOI is g*mm^2 and must be
        # converted to N*m/(deg/s^2) before adding. Do NOT reuse get_corrected_torque
        # (dead code, unit-inconsistent: multiplies g*mm^2 by deg/s^2).
        # Gated on _cal_is_real (D11, M2b): with no REAL matrix embedded there is no calibrated
        # output to derive -- an identity fallback would just echo raw voltages, so skip entirely.
        _ft_cal = None  # shared between the two derived/ blocks below
        if _cal_is_real:
            # Block A: forcetorque_calibrated -- re-anchored to aidata rows 0-5 via the channel
            # map (D11, M2b). Locate the 6 F/T rows by identity in input_channel_names (the same
            # source zipped into metadata/daq_ai_channel_map) rather than reusing the removed
            # in-memory forcetorque_raw, so derived is computed from the immutable aidata archive
            # x the embedded matrix -- exactly what R re-derives downstream.
            try:
                _combined = _concat_trial_records(trial_records)
                _ai = _combined.get('aidata', None)
                _ft_names = ['xForce', 'yForce', 'zForce', 'xTorque', 'yTorque', 'zTorque']
                _ch_names = list(getattr(bender, 'input_channel_names', []) or [])
                _ft_rows = [_ch_names.index(_nm) for _nm in _ft_names if _nm in _ch_names]
                if (
                    _ai is not None and np.asarray(_ai).size > 0
                    and len(_ft_rows) == 6 and np.asarray(_ai).ndim == 2
                ):
                    _ai = np.asarray(_ai, dtype=float)
                    # apply_calibration_forcetorque orientation: (raw6[6,N].T @ cal[6,6]).T -> [6,N].
                    _ft_cal = np.dot(_ai[_ft_rows, :].T, _cal_arr).T
                    g_der = f.create_group('derived')
                    g_der.attrs['note'] = (
                        'RA inspection only -- not ground truth. '
                        'R re-derives from raw aidata rows 0-5 (daq_ai_channel_map) + '
                        'embedded calibration_forcetorque_matrix VALUES.'
                    )
                    g_der.create_dataset('forcetorque_calibrated', data=_ft_cal)
            except Exception as _der_exc:
                _ft_cal = None
                logging.warning(
                    'export_primary_h5: derived/forcetorque_calibrated skipped: %s', _der_exc
                )

            # Block B: torque_inertia_corrected -- only runs when Block A succeeded.
            # Failure here does not affect forcetorque_calibrated.
            # _i_sum = specimen MOI (converted to N*m/(deg/s^2)) + empirical I_est if present.
            # I_est is the apparatus MOI from the empty-run calibration fit (inertial_calibration_profile);
            # it is NOT the old theoretical CLAMP_DIM baseline (superseded, see Deferred flag 4).
            # The apparatus term is I_est-only until the calibration-matrix workflow (Deferred flag 4)
            # delivers a per-configuration apparatus MOI.
            if _ft_cal is not None and 'derived' in f:
                try:
                    # Gather empirical apparatus MOI: I_est in native N*m/(deg/s^2).
                    _der_prof = getattr(bender, 'inertial_calibration_profile', None)
                    _i_app_native = None
                    if isinstance(_der_prof, dict):
                        _ie = _der_prof.get('I_est', None)
                        if _ie is not None:
                            _i_app_native = float(_ie)
                    # Gather specimen inertia: g*mm^2 -> N*m/(deg/s^2).
                    # Conversion: I[g*mm^2] * 1e-9 [kg*m^2] * (pi/180) [N*m/(deg/s^2)]
                    _i_spec_native = None
                    for _mattr in ('specimen_moi_specimen', 'specimen_moi_profile', 'specimen_moi_frustum'):
                        _sv = getattr(bender, _mattr, None)
                        if _sv is not None and np.isfinite(float(_sv)) and float(_sv) > 0.0:
                            _i_spec_native = float(_sv) * 1e-9 * (np.pi / 180.0)
                            break
                    if _i_app_native is not None or _i_spec_native is not None:
                        _i_sum = (_i_app_native or 0.0) + (_i_spec_native or 0.0)
                        _av_cmd = _combined.get('anglevel_cmd', None)
                        _ai_rate = float(
                            getattr(bender, 'daq_ai_sample_rate_hz', 1000.0) or 1000.0
                        )
                        if _av_cmd is not None and np.asarray(_av_cmd).size == _ft_cal.shape[1]:
                            _alpha = np.gradient(
                                np.asarray(_av_cmd, dtype=float).reshape(-1)
                            ) * _ai_rate  # deg/s^2
                            # Select primary torque row matching configured bending axis.
                            _ax = str(getattr(
                                bender, 'primary_bending_axis',
                                getattr(bender, 'bending_axis_sensor', 'z')
                            )).strip().lower()
                            _ax_row = {
                                'x': 3, 'xtorque': 3,
                                'y': 4, 'ytorque': 4,
                                'z': 5, 'ztorque': 5,
                            }.get(_ax, 5)
                            _tau_corrected = _ft_cal[_ax_row, :] - _i_sum * _alpha
                            f['derived'].create_dataset('torque_inertia_corrected', data=_tau_corrected)
                except Exception as _der_exc:
                    logging.warning(
                        'export_primary_h5: derived/torque_inertia_corrected skipped: %s', _der_exc
                    )

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
    keys_1d = ['angle_cmd', 'anglevel_cmd', 'angle_measured', 'S1stimcmd', 'S2stimcmd']
    parts_1d = {k: [] for k in keys_1d}
    ft_raw_parts = []
    ai_parts = []
    ai_rows = None  # channel count of the first record with a 2D aidata buffer
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
        # aidata: immutable raw ADC buffer [n_channels, n]; concatenated along the sample axis so
        # the derived/ block can re-anchor F/T to aidata rows via daq_ai_channel_map (D11, M2b).
        ai = np.asarray(r.get('aidata', np.array([])), dtype=float)
        if ai.ndim == 2 and ai.shape[1] == n:
            if ai_rows is None:
                ai_rows = ai.shape[0]
            ai_parts.append(ai if ai.shape[0] == ai_rows else None)
        else:
            ai_parts.append(None)
        dt = (tr[-1] - tr[0]) / max(1, n - 1)
        t_offset = float(tr0[-1] + dt)
    combined: dict = {'t': np.concatenate(t_parts) if t_parts else np.array([])}
    for k in keys_1d:
        combined[k] = np.concatenate(parts_1d[k]) if parts_1d[k] else np.array([])
    combined['forcetorque_raw'] = np.concatenate(ft_raw_parts, axis=1) if ft_raw_parts else np.array([])
    # Only expose a concatenated aidata when every recorded segment supplied a consistent buffer;
    # a missing/ragged segment makes a whole-timeline decode meaningless, so fall back to empty.
    if ai_rows is not None and ai_parts and all(p is not None for p in ai_parts):
        combined['aidata'] = np.concatenate(ai_parts, axis=1)
    else:
        combined['aidata'] = np.array([])
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
    off1 = _torque_row(ft_raw, axis_to_idx[off_axes[0]]) if len(off_axes) > 0 else np.array([])
    off2 = _torque_row(ft_raw, axis_to_idx[off_axes[1]]) if len(off_axes) > 1 else np.array([])

    # Sonomicrometry (optional). Calibrated length (mm, from apply_calibration_sono) lives on the
    # bender as sono_left_mm / sono_right_mm, aligned with t for a single-trial run; rec keys are
    # checked first for forward compatibility. A side is plotted only when its length matches t and
    # it has finite samples. When sono is absent (no instrumentation) or length-mismatched (e.g. a
    # multi-step concatenated timeline), no subplot is added — no error, no placeholder.
    def _sono_series(rec_key, attr_name):
        v = rec.get(rec_key, None)
        if v is None:
            v = getattr(bender, attr_name, None)
        if v is None:
            return np.array([])
        return np.asarray(v, dtype=float).reshape(-1)

    sono_left = _sono_series('sono_left', 'sono_left_mm')
    sono_right = _sono_series('sono_right', 'sono_right_mm')
    sono_specs = []  # (legend name, array, line color)
    if t.size > 0 and sono_left.size == t.size and np.any(np.isfinite(sono_left)):
        sono_specs.append(('sono_left', sono_left, 'mediumvioletred'))
    if t.size > 0 and sono_right.size == t.size and np.any(np.isfinite(sono_right)):
        sono_specs.append(('sono_right', sono_right, 'darkcyan'))
    has_sono = len(sono_specs) > 0
    _units_map = getattr(bender, 'units', None) or {}
    sono_unit = str(_units_map.get('sono_left') or _units_map.get('sono_right') or 'mm')

    n_rows = 6 if has_sono else 5
    row_specs = [[{'secondary_y': True}], [{}], [{}], [{}], [{}]]
    if has_sono:
        row_specs.append([{}])
    fig = make_subplots(
        rows=n_rows, cols=1, shared_xaxes=True, vertical_spacing=0.03,
        specs=row_specs,
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

    # RAW primary torque only. The inertial correction is applied post-hoc in R (empirical
    # apparatus inertia from an empty calibration run; analytic specimen inertia from geometry),
    # so the QC figure does not draw a "corrected" trace or an inertial-torque overlay.
    if t.size > 0 and primary_raw.size == t.size:
        fig.add_trace(
            go.Scatter(x=t, y=primary_raw, mode='lines', name=f'{axis_norm} raw', line=dict(color='firebrick')),
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

    sono_row = 6
    if has_sono:
        for _name, _arr, _color in sono_specs:
            fig.add_trace(
                go.Scatter(x=t, y=_arr, mode='lines', name=_name, line=dict(color=_color)),
                row=sono_row, col=1,
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
    bottom_row = sono_row if has_sono else 5
    if has_sono:
        fig.update_yaxes(title_text=f'Sono length ({sono_unit})', row=sono_row, col=1)
    fig.update_xaxes(title_text='Time (s)', row=bottom_row, col=1)
    qc_cap = f'QC: {proc} | ' + ('all steps' if qc_trial_index == 'all' else f'trial {qc_trial_index}')
    if proc in ('isometric', 'isovelocity'):
        qc_cap = qc_cap + "<br><sup style='font-size:11px'>" + bender.strain_geometry_plot_context() + '</sup>'
    fig_height = 1400 + (240 if has_sono else 0)
    fig.update_layout(height=fig_height, width=1100, title_text=qc_cap, showlegend=True, hovermode='x unified')
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
