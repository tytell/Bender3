"""
Emergency hardware stop for NI-DAQ: reset device(s) to cease AI/AO/DO activity.

Used by the Streamlit **kill switch**. Safe to call when nidaqmx is missing (no-op message).
"""
from __future__ import annotations

from typing import Optional, Tuple


def daq_emergency_stop(
    device_name: Optional[str] = None,
    release_motor_enable_line: Optional[str] = None,
) -> Tuple[bool, str]:
    """
    Reset the named NI-DAQ device, or all local devices if ``device_name`` is empty.

    ``release_motor_enable_line`` (e.g. ``'port0/line2'``): when given with a named device, the
    motor ENABLE line's persistent digital power-up state is forced to TRISTATE *before* the device
    reset, so an emergency stop always RELEASES (de-energizes) the motor and it stays released --
    overriding any power-up HIGH configured for normal runs (Bender._ensure_motor_enable_power_up_state).
    Without this, reset_device() would return the line to a power-up HIGH state and leave the motor
    energized/holding after an e-stop. An e-stop must never leave the motor powered.

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
                try:
                    from nidaqmx.types import DOPowerUpState
                    from nidaqmx.constants import PowerUpStates
                    line = f"{name}/{str(release_motor_enable_line).strip().lstrip('/')}"
                    local.set_digital_power_up_states(
                        name,
                        [DOPowerUpState(physical_channel=line, power_up_state=PowerUpStates.TRISTATE)],
                    )
                    release_msg = ' Motor ENABLE released (tristate).'
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
