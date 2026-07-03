"""Phase 1 canonical export schema (lockstep regression).

Runs sim-mode protocols end-to-end and asserts the canonical HDF5 contract:
ledger-driven metadata, canonical timeseries names, per-sample index/stim arrays,
dropped fields (trial_index subgroup / step_order / forcetorque / protocol_metadata),
flat index_step_* parallel arrays (M2), protocol direction attrs,
the block_sequence -> protocol_block_sequence JSON route, and an empty 99_Unrouted.
"""
import json
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))
sys.path.append(str(Path(__file__).resolve().parents[1] / 'templates' / 'configs'))

import h5py  # noqa: E402
from bender_functions import Bender  # noqa: E402
from bender_h5_export import export_primary_h5  # noqa: E402


def _common(b):
    b.session_simulated = True
    b.dclamp = 10.0
    b.xsec_width = 2.0
    b.is_stim = False
    b.stim_pulse_rate = 75.0


def _isometric():
    b = Bender('jimenez_bender_config_A')
    _common(b)
    b.isometric_initial = -5.0
    b.isometric_final = 5.0
    b.isometric_num_steps = 3
    b.isometric_mode = 'angle'
    return b, 'isometric'


def _dynamic():
    b = Bender('jimenez_bender_config_A')
    _common(b)
    b.all_freqs = [1.0, 2.0]
    b.all_amps = [0.05]
    b.all_amps_mode = 'strain'
    b.cycles_per_step = 2
    b.n_end_cycles = 1
    b.randomize = False
    b.stim_cycles_in_step = []
    return b, 'dynamic'


def _run(maker, tmp_path, *, calibration=None):
    b, tt = maker()
    out = str(tmp_path / f'{tt}.h5')
    b.outputfile = out
    b.run_experiment(test_type=tt)
    # Optionally pin the F/T calibration matrix so the identity-refusal branch (M2b) is
    # deterministic regardless of whether the config's .cal file resolved on this machine.
    if calibration is not None:
        b.calibration = np.asarray(calibration, dtype=float)
    res = export_primary_h5(b, outputfile=out)
    return res['outputfile']


# A non-identity 6x6 stands in for a REAL ATI calibration matrix (M2b embeds it; identity is
# refused). Diagonal keeps the calibrated F/T finite and easy to reason about.
_REAL_CAL = np.diag([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]).astype(float)


def test_no_unrouted_and_no_legacy_groups(tmp_path):
    for maker in (_isometric, _dynamic):
        path = _run(maker, tmp_path)
        with h5py.File(path, 'r') as f:
            assert '99_Unrouted' not in f, '99_Unrouted must not be at root (moved under metadata)'
            assert 'protocol_metadata' not in f['metadata'], 'protocol_metadata subgroup removed'
            assert 'metadata' in f, 'root must have metadata group'
            assert 'timeseries' in f, 'root must have timeseries group'
            assert '01_Metadata' not in f, 'old 01_Metadata group must be absent'
            assert '02_TimeSeries' not in f, 'old 02_TimeSeries group must be absent'
            assert 'trial_index' not in f['metadata'], 'trial_index subgroup must be absent (replaced by index_step_* arrays in M2)'


def test_step_protocol_per_sample_index_arrays(tmp_path):
    path = _run(_isometric, tmp_path)
    with h5py.File(path, 'r') as f:
        ts = f['timeseries']
        for tn in ts.keys():
            tg = ts[tn]
            n = tg['time_second'].shape[0]
            for fld in ('cycle_index', 'step_index', 'sequence_index', 'block_index',
                        'block_direction', 'block_stim_sides'):
                assert fld in tg, f'{tn}: missing per-sample {fld}'
                assert tg[fld].shape[0] == n, f'{tn}: {fld} length != N'
            # cycle_index is -1 throughout for step protocols (no numbered cycles)
            assert np.all(tg['cycle_index'][:] == -1)
            # converted fields must NOT remain as scalar trial attrs
            for fld in ('trial_index', 'cycle_index', 'step_index', 'sequence_index',
                        'block_index', 'block_direction', 'block_stim_sides'):
                assert fld not in tg.attrs, f'{tn}: {fld} should not be a scalar attr'


def test_dynamic_has_cycle_index_but_no_step_index(tmp_path):
    path = _run(_dynamic, tmp_path)
    with h5py.File(path, 'r') as f:
        # continuous (single_finite): datasets written flat directly under timeseries, no subgroup
        ts = f['timeseries']
        assert 'trial_0000' not in ts, 'continuous layout must not create a trial_0000 subgroup'
        assert 'time_second' in ts, 'continuous layout must have time_second flat under timeseries'
        n = ts['time_second'].shape[0]
        assert 'cycle_index' in ts and ts['cycle_index'].shape[0] == n
        assert 'step_index' not in ts, 'dynamic has no discrete shuffled steps'


def test_forcetorque_raw_absent_and_aidata_present(tmp_path):
    """M2b (D11): timeseries holds ONLY raw aidata; no calibrated forcetorque copy is written."""
    path = _run(_isometric, tmp_path)
    with h5py.File(path, 'r') as f:
        # Isometric is segmented_finite: per-step subgroups are step_NNN (one-based).
        tg = f['timeseries']['step_001']
        assert 'forcetorque_raw' not in tg, 'forcetorque_raw must not be written to timeseries (D11, M2b)'
        assert 'forcetorque' not in tg, 'calibrated forcetorque must not be in timeseries'
        assert 'aidata' in tg, 'raw aidata must remain the timeseries F/T archive (rows 0-5)'
        assert np.asarray(tg['aidata']).ndim == 2, 'aidata must be the 2D [n_channels, N] buffer'
        for split in ('force_x_newton', 'torque_z_newton_meter'):
            assert split not in tg, 'F/T must not be split into per-channel datasets'


def test_calibration_matrix_and_derived_present_when_real(tmp_path):
    """M2b: a REAL matrix is embedded in metadata and derived/forcetorque_calibrated is written.

    Covers both sampling modes: segmented_finite (isometric) and single_finite (dynamic).
    """
    for maker in (_isometric, _dynamic):
        path = _run(maker, tmp_path, calibration=_REAL_CAL)
        with h5py.File(path, 'r') as f:
            meta = f['metadata']
            assert 'calibration_forcetorque_matrix' in meta, 'real matrix must be embedded in metadata'
            assert tuple(np.asarray(meta['calibration_forcetorque_matrix']).shape) == (6, 6)
            assert 'derived' in f and 'forcetorque_calibrated' in f['derived'], (
                'derived/forcetorque_calibrated must be present when the matrix is real'
            )
            ftc = np.asarray(f['derived']['forcetorque_calibrated'])
            assert ftc.ndim == 2 and ftc.shape[0] == 6, 'calibrated F/T must be [6, N]'


def test_calibration_matrix_and_derived_absent_on_identity(tmp_path):
    """M2b: an np.eye(6) identity fallback is REFUSED -- no matrix, no calibrated derived output."""
    path = _run(_isometric, tmp_path, calibration=np.eye(6))
    with h5py.File(path, 'r') as f:
        assert 'calibration_forcetorque_matrix' not in f['metadata'], (
            'identity matrix must NOT be embedded (absent = "not calibrated", M2b)'
        )
        if 'derived' in f:
            assert 'forcetorque_calibrated' not in f['derived'], (
                'no calibrated output may be written without a real matrix'
            )


def test_index_cycle_arrays_single_finite_only(tmp_path):
    """M2c (D10): single_finite writes the index_cycle_* per-cycle design grid; segmented does not."""
    # single_finite (dynamic): index_cycle_* present as parallel arrays of equal length C.
    path = _run(_dynamic, tmp_path)
    with h5py.File(path, 'r') as f:
        meta = f['metadata']
        assert 'index_cycle_index' in meta, 'single_finite must write index_cycle_index'
        C = int(meta['index_cycle_index'].shape[0])
        assert C > 0
        assert int(meta['index_cycle_index'][0]) == 1, 'index_cycle_index is 1-based'
        for key in (
            'index_cycle_frequency_hertz',
            'index_cycle_motor_amplitude_degree',
            'index_cycle_active',
            'index_cycle_activation_duty',
            'index_cycle_activation_phase',
            'index_cycle_operating_point',
            'index_cycle_operating_point_units',
        ):
            assert key in meta, f'missing per-cycle array: {key}'
            assert meta[key].shape[0] == C, f'{key} length != C'
        # protocol scalars
        assert int(meta.attrs['protocol_cycle_count']) == C
        assert 'protocol_activation_start_cycle' in meta.attrs
        # protocol_end_cycle_count (D10 name; routed from n_end_cycles)
        assert 'protocol_end_cycle_count' in meta.attrs
        assert 'protocol_end_cycles_count' not in meta.attrs, 'old plural key must be gone'
        # D4 intact: the per-cycle block adds no index_step_operating_point-style CYCLE columns
        # into the step family (no index_step_cycle_* cross-contamination).
        assert not any(k.startswith('index_step_cycle') for k in meta.keys()), (
            'index_cycle_* must not create index_step_cycle_* columns (D4 intact)'
        )

    # segmented_finite (isometric): no per-cycle design grid.
    seg = _run(_isometric, tmp_path)
    with h5py.File(seg, 'r') as f:
        meta = f['metadata']
        assert not any(k.startswith('index_cycle_') for k in meta.keys()), (
            'segmented_finite must not write index_cycle_* arrays'
        )
        assert 'protocol_cycle_count' not in meta.attrs, 'segmented has no per-cycle scalar block'


def test_block_sequence_is_json_and_step_order_dropped(tmp_path):
    path = _run(_isometric, tmp_path)
    with h5py.File(path, 'r') as f:
        meta = f['metadata']
        assert 'protocol_block_sequence' in meta.attrs, 'block_sequence must be routed canonically'
        parsed = json.loads(meta.attrs['protocol_block_sequence'])
        assert isinstance(parsed, list) and isinstance(parsed[0], dict)
        # step_order must be gone everywhere (no protocol_metadata, no metadata key)
        assert 'step_order' not in meta and 'step_order' not in meta.attrs


def test_index_step_arrays_replace_trial_index_subgroup(tmp_path):
    """M2: trial_index subgroup removed; flat index_step_* parallel arrays written to metadata."""
    path = _run(_isometric, tmp_path)
    with h5py.File(path, 'r') as f:
        meta = f['metadata']
        assert 'trial_index' not in meta, 'trial_index subgroup must be absent after M2'
        # Shared fields present as datasets (not attrs)
        for key in (
            'index_step_number',
            'index_step_duration_second',
            'index_step_rest_before_second',
            'index_step_wall_clock_start',
            'index_step_recruitment',
            'index_step_stim_t0_second',
            'index_step_stim_t1_second',
        ):
            assert key in meta, f'missing expected index_step dataset: {key}'
        # Isometric-specific fields
        assert 'index_step_target_angle_degree' in meta, 'isometric target angle missing'
        assert 'index_step_ramp_from_angle_degree' in meta, 'isometric ramp-from angle missing'
        # Array length matches n_trials
        n = int(meta.attrs['session_step_count'])
        assert meta['index_step_number'].shape[0] == n, 'index_step_number length mismatch'
        assert meta['index_step_duration_second'].shape[0] == n, 'index_step_duration_second length mismatch'
        # step_number is 1-based
        assert int(meta['index_step_number'][0]) == 1, 'step_number must be 1-based'
        # Protocol direction attrs
        assert 'protocol_motor_positive_bend_direction' in meta.attrs
        assert 'protocol_sensor_positive_bend_direction' in meta.attrs
        v = meta.attrs['protocol_motor_positive_bend_direction']
        assert v in ('left', 'right', 'none'), f'unexpected protocol_motor_positive_bend_direction value: {v}'


def test_root_holds_only_schema_version(tmp_path):
    for maker in (_isometric, _dynamic):
        path = _run(maker, tmp_path)
        with h5py.File(path, 'r') as f:
            root_attrs = list(f.attrs.keys())
            assert root_attrs == ['schema_version'], (
                f'Root must have only schema_version; got {root_attrs}'
            )
            assert 'test_type' not in f.attrs, 'test_type must not be at root (moved to metadata/protocol_type)'
            assert 'filename' not in f.attrs, 'filename must not be at root (in metadata)'
            assert 'start_time_iso' not in f.attrs, 'start_time_iso must not be at root (metadata/session_date)'


def test_metadata_has_session_date_and_no_calibration_link_subgroup(tmp_path):
    path = _run(_isometric, tmp_path)
    with h5py.File(path, 'r') as f:
        meta = f['metadata']
        assert 'session_date' in meta.attrs, 'session_date must be a flat metadata attr'
        assert 'calibration_link' not in meta, 'calibration_link subgroup must be absent (D3/D12: replaced by flat attr)'
