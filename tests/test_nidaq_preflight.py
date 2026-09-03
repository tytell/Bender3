"""Unit tests for NI-DAQ preflight probe (no hardware required)."""
from __future__ import annotations

import sys
from types import ModuleType, SimpleNamespace
from unittest import mock

import bender_daq_kill as daq_kill


def _clear_nidaqmx_modules():
    return {
        k: sys.modules.pop(k)
        for k in list(sys.modules)
        if k == 'nidaqmx' or k.startswith('nidaqmx.')
    }


def test_probe_package_missing():
    saved = _clear_nidaqmx_modules()
    try:
        # Sentinel None makes ``import nidaqmx...`` fail as if the package is absent.
        with mock.patch.dict(sys.modules, {'nidaqmx': None, 'nidaqmx.system': None}):
            status, message, devices = daq_kill.probe_nidaq_status()
    finally:
        sys.modules.update(saved)
    assert status == 'package_missing'
    assert 'nidaqmx' in message.lower()
    assert devices == []


def test_probe_drivers_missing_on_system_local():
    fake_nidaqmx = ModuleType('nidaqmx')
    fake_system = ModuleType('nidaqmx.system')

    class _Boom:
        @staticmethod
        def local():
            raise OSError('Failed to load nicaiu.dll — NI-DAQmx driver not installed')

    fake_system.System = _Boom
    fake_nidaqmx.system = fake_system
    saved = _clear_nidaqmx_modules()
    try:
        with mock.patch.dict(
            sys.modules,
            {'nidaqmx': fake_nidaqmx, 'nidaqmx.system': fake_system},
        ):
            status, message, devices = daq_kill.probe_nidaq_status()
    finally:
        sys.modules.update(saved)
    assert status == 'drivers_missing'
    assert 'driver' in message.lower()
    assert devices == []


def test_probe_no_devices():
    fake_nidaqmx = ModuleType('nidaqmx')
    fake_system = ModuleType('nidaqmx.system')

    class _Empty:
        @staticmethod
        def local():
            return SimpleNamespace(devices=[])

    fake_system.System = _Empty
    fake_nidaqmx.system = fake_system
    saved = _clear_nidaqmx_modules()
    try:
        with mock.patch.dict(
            sys.modules,
            {'nidaqmx': fake_nidaqmx, 'nidaqmx.system': fake_system},
        ):
            status, message, devices = daq_kill.probe_nidaq_status()
    finally:
        sys.modules.update(saved)
    assert status == 'no_devices'
    assert 'no local devices' in message.lower()
    assert devices == []


def test_probe_ok_lists_devices():
    fake_nidaqmx = ModuleType('nidaqmx')
    fake_system = ModuleType('nidaqmx.system')

    class _Ok:
        @staticmethod
        def local():
            return SimpleNamespace(
                devices=[SimpleNamespace(name='Dev1'), SimpleNamespace(name='Dev2')]
            )

    fake_system.System = _Ok
    fake_nidaqmx.system = fake_system
    saved = _clear_nidaqmx_modules()
    try:
        with mock.patch.dict(
            sys.modules,
            {'nidaqmx': fake_nidaqmx, 'nidaqmx.system': fake_system},
        ):
            status, message, devices = daq_kill.probe_nidaq_status()
    finally:
        sys.modules.update(saved)
    assert status == 'ok'
    assert devices == ['Dev1', 'Dev2']
    assert 'Dev1' in message
