"""Unit tests for strict GUI JSON persistence."""
import json
import os
import tempfile

import numpy as np
import pytest

from bender_json_persistent import JsonPersistTypeError, to_json_persistent
import bender_protocol_templates as pt


def test_to_json_persistent_ndarray_becomes_list():
    arr = np.array([1.0, 2.5, 3.0])
    out = to_json_persistent(arr)
    assert out == [1.0, 2.5, 3.0]
    json.dumps(out)


def test_to_json_persistent_numpy_scalar():
    assert to_json_persistent(np.float64(2.5)) == 2.5
    assert to_json_persistent(np.int32(7)) == 7


def test_to_json_persistent_rejects_opaque_object():
    class _Handle:
        pass

    with pytest.raises(JsonPersistTypeError, match='Non-JSON-persistable'):
        to_json_persistent(_Handle(), path='gui_bad')


def test_snapshot_bender_procedure_serializes_ndarray_fields():
    class _B:
        all_freqs = np.array([1.0, 2.0, 3.0])
        all_amps = np.array([0.05, 0.1])

    schema = {'isometric_required': [], 'isometric_optional': [], 'isovelocity_required': [], 'isovelocity_optional': [], 'calibration_optional': []}
    snap = pt.snapshot_bender_procedure(_B(), schema, 'dynamic')
    assert snap['all_freqs'] == [1.0, 2.0, 3.0]
    assert snap['all_amps'] == [0.05, 0.1]


def test_save_protocol_template_roundtrip_ndarray_procedure():
    fd, path = tempfile.mkstemp(suffix='.json')
    os.close(fd)
    try:
        proc = {'all_freqs': np.array([5.0, 10.0]), 'all_amps': [0.02, 0.04]}
        pt.save_protocol_template(
            path,
            name='t',
            description='',
            test_type='dynamic',
            procedure=proc,
        )
        loaded = pt.load_protocol_template(path)
        assert loaded['procedure']['all_freqs'] == [5.0, 10.0]
        assert loaded['procedure']['all_amps'] == [0.02, 0.04]
        assert isinstance(loaded['procedure']['all_freqs'], list)
        assert not isinstance(loaded['procedure']['all_freqs'], str)
    finally:
        os.remove(path)
