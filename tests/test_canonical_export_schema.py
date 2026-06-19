"""Phase 0 canonical export schema (lockstep regression).

Runs sim-mode protocols end-to-end and asserts the canonical HDF5 contract:
ledger-driven metadata, canonical timeseries names, per-sample index/stim arrays,
dropped fields (trial_index / step_order / forcetorque / protocol_metadata), the
block_sequence -> protocol_block_sequence JSON route, and an empty 99_Unrouted.
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


def _run(maker, tmp_path):
    b, tt = maker()
    out = str(tmp_path / f'{tt}.h5')
    b.outputfile = out
    b.run_experiment(test_type=tt)
    res = export_primary_h5(b, outputfile=out)
    return res['outputfile']


def test_no_unrouted_and_no_legacy_groups(tmp_path):
    for maker in (_isometric, _dynamic):
        path = _run(maker, tmp_path)
        with h5py.File(path, 'r') as f:
            assert '99_Unrouted' not in f, '99_Unrouted must be empty/absent'
            assert 'protocol_metadata' not in f['01_Metadata'], 'protocol_metadata subgroup removed'


def test_step_protocol_per_sample_index_arrays(tmp_path):
    path = _run(_isometric, tmp_path)
    with h5py.File(path, 'r') as f:
        ts = f['02_TimeSeries']
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
        # continuous (single_finite): datasets written flat directly under 02_TimeSeries, no subgroup
        ts = f['02_TimeSeries']
        assert 'trial_0000' not in ts, 'continuous layout must not create a trial_0000 subgroup'
        assert 'time_second' in ts, 'continuous layout must have time_second flat under 02_TimeSeries'
        n = ts['time_second'].shape[0]
        assert 'cycle_index' in ts and ts['cycle_index'].shape[0] == n
        assert 'step_index' not in ts, 'dynamic has no discrete shuffled steps'


def test_forcetorque_raw_kept_not_split(tmp_path):
    path = _run(_isometric, tmp_path)
    with h5py.File(path, 'r') as f:
        # Isometric is segmented_finite: per-step subgroups are step_NNN (one-based), not trial_XXXX.
        tg = f['02_TimeSeries']['step_001']
        assert 'forcetorque_raw' in tg
        assert tg['forcetorque_raw'].shape[0] == 6
        assert 'forcetorque' not in tg, 'duplicate forcetorque dropped'
        for split in ('force_x_newton', 'torque_z_newton_meter'):
            assert split not in tg, 'F/T must not be split into per-channel datasets'


def test_block_sequence_is_json_and_step_order_dropped(tmp_path):
    path = _run(_isometric, tmp_path)
    with h5py.File(path, 'r') as f:
        meta = f['01_Metadata']
        assert 'protocol_block_sequence' in meta.attrs, 'block_sequence must be routed canonically'
        parsed = json.loads(meta.attrs['protocol_block_sequence'])
        assert isinstance(parsed, list) and isinstance(parsed[0], dict)
        # step_order must be gone everywhere (no protocol_metadata, no metadata key)
        assert 'step_order' not in meta and 'step_order' not in meta.attrs


def test_manifest_drops_trial_index_column(tmp_path):
    path = _run(_isometric, tmp_path)
    with h5py.File(path, 'r') as f:
        gi = f['01_Metadata/trial_index']
        assert 'trial_names' in gi
        assert 'trial_index' not in gi, 'trial_index column dropped'
        assert 'cycle_index' not in gi, 'cycle_index column dropped (now per-sample)'
