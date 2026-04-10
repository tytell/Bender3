"""
Emergency hardware stop for NI-DAQ: reset device(s) to cease AI/AO/DO activity.

Used by the Streamlit **kill switch**. Safe to call when nidaqmx is missing (no-op message).
"""
from __future__ import annotations

from typing import Optional, Tuple


def daq_emergency_stop(device_name: Optional[str] = None) -> Tuple[bool, str]:
    """
    Reset the named NI-DAQ device, or all local devices if ``device_name`` is empty.

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
            dev = local.devices[name]
            dev.reset_device()
            return True, f'NI device {name!r} was reset (tasks stopped, outputs cleared).'

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
