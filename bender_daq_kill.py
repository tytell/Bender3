"""
Emergency hardware stop for NI-DAQ: reset device(s) to cease AI/AO/DO activity.

Used by the Streamlit **kill switch**. Safe to call when nidaqmx is missing (no-op message).

Also provides ``probe_nidaq_status`` for GUI preflight (package vs drivers vs devices).
"""
from __future__ import annotations

from typing import List, Optional, Tuple


def probe_nidaq_status() -> Tuple[str, str, List[str]]:
    """Probe NI-DAQmx availability without resetting hardware.

    Distinguishes:
      - ``package_missing``: Python ``nidaqmx`` import failed
      - ``drivers_missing``: package imports but the NI-DAQmx driver/runtime is not usable
      - ``no_devices``: drivers respond but no local devices are listed
      - ``ok``: drivers respond and at least one device name is present
      - ``error``: unexpected failure while probing

    Returns ``(status, message, device_names)``. ``device_names`` is empty unless status is
    ``ok`` (or rarely ``no_devices`` stays empty by definition).
    """
    try:
        import nidaqmx.system as nisys
    except ImportError:
        return (
            'package_missing',
            'Python package nidaqmx is not installed. Install it in this environment, '
            'or use simulation / non-DAQ routes. NI-DAQmx Runtime drivers are separate '
            'system software from National Instruments.',
            [],
        )
    except Exception as e:
        # Rare: import-time DLL failure can surface as OSError / DaqError, not ImportError.
        msg = f'{type(e).__name__}: {e}'
        low = msg.lower()
        if any(
            token in low
            for token in (
                'daqmx',
                'nidaq',
                'nicaiu',
                'dll',
                'driver',
                'library',
                'shared object',
            )
        ):
            return (
                'drivers_missing',
                'NI-DAQmx drivers were not detected (Python package present, but the '
                f'driver/runtime failed to load: {msg}). Install NI-DAQmx Runtime from '
                'National Instruments, then restart the app. Simulation / non-DAQ routes '
                'still work without drivers.',
                [],
            )
        return 'error', f'Unexpected NI-DAQ probe failure during import: {msg}', []

    try:
        local = nisys.System.local()
        names = [str(d.name) for d in local.devices]
    except Exception as e:
        msg = f'{type(e).__name__}: {e}'
        low = msg.lower()
        if any(
            token in low
            for token in (
                'daqmx',
                'nidaq',
                'nicaiu',
                'dll',
                'driver',
                'library',
                'not installed',
                'not found',
                'shared object',
            )
        ):
            return (
                'drivers_missing',
                'NI-DAQmx drivers were not detected (Python package present, but the '
                f'driver/runtime is not usable: {msg}). Install NI-DAQmx Runtime from '
                'National Instruments, then restart the app. Simulation / non-DAQ routes '
                'still work without drivers.',
                [],
            )
        return 'error', f'Could not query local NI-DAQ system: {msg}', []

    if not names:
        return (
            'no_devices',
            'NI-DAQmx drivers responded, but no local devices were found. Check USB '
            'connection / device power, then confirm the device appears in NI MAX.',
            [],
        )
    return 'ok', f'NI-DAQmx OK — devices: {", ".join(names)}', names


def daq_emergency_stop(
    device_name: Optional[str] = None,
    release_motor_enable_line: Optional[str] = None,
) -> Tuple[bool, str]:
    """
    Reset the named NI-DAQ device, or all local devices if ``device_name`` is empty.

    ``release_motor_enable_line`` (e.g. ``'port0/line2'``): when given with a named device, the
    persistent digital power-up state of the ENTIRE port containing that line is forced to
    TRISTATE *before* the device reset, so an emergency stop always RELEASES (de-energizes) the
    motor and it stays released -- overriding any power-up HIGH configured for normal runs
    (Bender._ensure_motor_enable_power_up_state). Without this, reset_device() would return the
    ENABLE line to a power-up HIGH state and leave the motor energized/holding after an e-stop.
    An e-stop must never leave the motor powered.

    Whole port, not just the one line: NI devices like the USB-6361 reject per-line power-up
    states with error -200652 ("must specify programmable powerup state for entire ports"), so a
    single-line set silently failed and the release depended on the factory-default TRISTATE.
    All lines of the motor port go TRISTATE here -- STEP/DIR floating is harmless with the driver
    de-energized.

    Returns ``(ok, message)``. ``ok`` is True if at least one reset was attempted without
    fatal import errors; peripheral errors are summarized in the message.
    """
    try:
        import nidaqmx.system as nisys
    except ImportError:
        return False, 'nidaqmx is not installed — no hardware reset was performed.'

    try:
        local = nisys.System.local()
        name = (device_name or '').strip()
        if name:
            release_msg = ''
            if release_motor_enable_line:
                # Release the motor: make the device's power-up state de-energize ENABLE so the
                # reset below lands on (and stays at) TRISTATE, not the run-time energized HIGH.
                # The device requires the WHOLE port (-200652), so every line of the motor port
                # is specified.
                try:
                    from nidaqmx.types import DOPowerUpState
                    from nidaqmx.constants import PowerUpStates
                    line = f"{name}/{str(release_motor_enable_line).strip().lstrip('/')}"
                    port_path = line.rsplit('/', 1)[0]            # e.g. 'Dev1/port0'
                    port_lines = [
                        ch.name for ch in local.devices[name].do_lines
                        if ch.name.startswith(port_path + '/')
                    ]
                    if not port_lines:
                        raise ValueError(f'no DO lines found under {port_path!r}')
                    local.set_digital_power_up_states(
                        name,
                        [
                            DOPowerUpState(physical_channel=ch, power_up_state=PowerUpStates.TRISTATE)
                            for ch in port_lines
                        ],
                    )
                    release_msg = (
                        f' Motor port {port_path!r} power-up forced TRISTATE '
                        f'({len(port_lines)} lines) -- motor released.'
                    )
                except Exception as e:
                    release_msg = f' WARNING: could not release motor ENABLE line: {e}'
            dev = local.devices[name]
            dev.reset_device()
            return True, f'NI device {name!r} was reset (tasks stopped, outputs cleared).{release_msg}'

        errs: list[str] = []
        n = 0
        for d in local.devices:
            try:
                d.reset_device()
                n += 1
            except Exception as e:
                errs.append(f'{d.name}: {e}')
        if n == 0 and not errs:
            return False, 'No NI devices found in the local system.'
        msg = f'Reset {n} NI device(s).'
        if errs:
            msg += ' Some issues: ' + '; '.join(errs[:5])
            if len(errs) > 5:
                msg += f' … (+{len(errs) - 5} more)'
        return True, msg
    except Exception as e:
        return False, f'{type(e).__name__}: {e}'
