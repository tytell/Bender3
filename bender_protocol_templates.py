"""
Named protocol templates: saved **procedure** parameters (test_type + form fields), not biometrics.

Supports every GUI experiment type: **dynamic**, **frequency_sweep**, **frequency_step**,
**curvature_step**, **step_change**, **isometric**, **isovelocity**, and **calibration**.
**calibration** templates may include an optional ``base_protocol`` block (motion parameters
for ``calibration_base_test_type``) so load + **Apply** restores both layers.

Templates are JSON files you can share or version-control. Load fills Streamlit ``fld_*`` session
keys and ``test_type_select``; the user then **Apply**es to the ``Bender`` instance. Biometrics
(geometry, profile, species) are separate—use **Load / Save biometrics** in the Streamlit GUI
(section 2) or enter them manually after loading a protocol template.
"""
from __future__ import annotations

import json
import os
import re
from datetime import datetime, timezone
from typing import Any, Callable, Dict, List, Optional, Tuple

PROTOCOL_TEMPLATE_VERSION = 2

# Motion-series keys per test_type (must match ``_motion_parameter_rows`` in ``bender_streamlit_gui``).
_MOTION_STIM = (
    'is_stim',
    'stim_pulse_rate',
    'S1volts',
    'S2volts',
)
_MOTION_FREQ_AMP_MODE = ('all_freqs', 'all_amps', 'all_amps_mode')
_MOTION_ORGANIZE = (
    'cycles_per_step',
    'n_end_cycles',
    'randomize',
    'random_seed',
    'stim_cycles_in_step',
)

MOTION_PROCEDURE_KEYS: Dict[str, List[str]] = {
    'dynamic': [
        *_MOTION_FREQ_AMP_MODE,
        *_MOTION_ORGANIZE,
        *_MOTION_STIM,
    ],
    'frequency_sweep': [
        'duration',
        *_MOTION_FREQ_AMP_MODE,
        'amplitude_frequency_exponent',
        *_MOTION_STIM,
    ],
    'frequency_step': ['duration', *_MOTION_FREQ_AMP_MODE, *_MOTION_STIM],
    'curvature_step': ['duration', *_MOTION_FREQ_AMP_MODE, *_MOTION_STIM],
    'step_change': [
        'step_change_frequencies',
        'step_change_curves',
        'step_change_cycles_per_step',
        *_MOTION_STIM,
    ],
}

# Keys edited as JSON text areas in the GUI
JSON_DICT_KEYS = frozenset(
    {
        'isometric_stim_params',
        'isometric_stim_overrides',
        'isovelocity_stim_params',
        'isovelocity_stim_overrides',
    }
)

LIST_FLOAT_KEYS = frozenset(
    {
        'all_freqs',
        'all_amps',
        'step_change_frequencies',
        'step_change_curves',
    }
)

LIST_INT_KEYS = frozenset(
    {
        'stim_cycles_in_step',
        'step_change_cycles_per_step',
    }
)


def default_templates_dir(project_root: str) -> str:
    return os.path.normpath(os.path.join(project_root, 'ProtocolTemplates'))


def sanitize_template_filename_stem(name: str) -> str:
    raw = (name or '').strip()
    if not raw:
        return 'protocol'
    s = re.sub(r'[^\w\s\-]', '_', raw, flags=re.UNICODE)
    s = re.sub(r'\s+', '_', s.strip())
    s = s.strip('_') or 'protocol'
    return s[:120]


def list_template_files(folder: str) -> List[str]:
    """Sorted paths to ``*.json`` in ``folder`` (non-recursive)."""
    d = str(folder or '').strip()
    if not d or not os.path.isdir(d):
        return []
    out = []
    for n in sorted(os.listdir(d)):
        if not n.lower().endswith('.json'):
            continue
        fp = os.path.join(d, n)
        if os.path.isfile(fp):
            out.append(fp)
    return out


def template_display_label(path: str) -> str:
    try:
        with open(path, encoding='utf-8') as f:
            data = json.load(f)
        name = str(data.get('name') or '').strip()
        tt = str(data.get('test_type') or '').strip()
        base = os.path.splitext(os.path.basename(path))[0]
        if name:
            return f'{name}  ({tt})' if tt else name
        return f'{base}  ({tt})' if tt else base
    except (OSError, json.JSONDecodeError, TypeError):
        return os.path.basename(path)


def _to_json_safe(v: Any) -> Any:
    if hasattr(v, 'item') and callable(getattr(v, 'item', None)):
        try:
            return _to_json_safe(v.item())
        except Exception:
            pass
    if v is None or isinstance(v, (str, bool)):
        return v
    if isinstance(v, (int, float)):
        if isinstance(v, float) and not (v == v):  # NaN
            return None
        return v
    if isinstance(v, dict):
        return {str(kk): _to_json_safe(vv) for kk, vv in v.items()}
    if isinstance(v, (list, tuple)):
        return [_to_json_safe(x) for x in v]
    return str(v)


def save_protocol_template(
    path: str,
    *,
    name: str,
    description: str,
    test_type: str,
    procedure: Dict[str, Any],
    base_protocol: Optional[Dict[str, Any]] = None,
) -> None:
    payload = {
        'version': PROTOCOL_TEMPLATE_VERSION,
        'name': str(name or '').strip() or os.path.splitext(os.path.basename(path))[0],
        'description': str(description or '').strip(),
        'saved_utc': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
        'test_type': str(test_type),
        'procedure': {str(k): _to_json_safe(v) for k, v in procedure.items()},
    }
    if base_protocol is not None:
        payload['base_protocol'] = _to_json_safe(base_protocol)
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(payload, f, indent=2, allow_nan=False)


def load_protocol_template(path: str) -> Dict[str, Any]:
    with open(path, encoding='utf-8') as f:
        return json.load(f)


def procedure_field_names_for_test_type(schema: Dict[str, Any], tt: str) -> List[str]:
    """Bender attribute names used by the GUI for one experiment ``test_type``."""
    if tt == 'isometric':
        return list(schema['isometric_required']) + list(schema['isometric_optional'])
    if tt == 'isovelocity':
        return list(schema['isovelocity_required']) + list(schema['isovelocity_optional'])
    if tt == 'calibration':
        return ['calibration_base_test_type'] + list(schema['calibration_optional'])
    if tt in MOTION_PROCEDURE_KEYS:
        return list(MOTION_PROCEDURE_KEYS[tt])
    return []


def snapshot_bender_procedure(b: Any, schema: Dict[str, Any], tt: str) -> Dict[str, Any]:
    """Read procedure-relevant attributes from a ``Bender`` for embedding (e.g. calibration base protocol)."""
    out: Dict[str, Any] = {}
    for name in procedure_field_names_for_test_type(schema, tt):
        if hasattr(b, name):
            out[name] = _to_json_safe(getattr(b, name))
    return out


def inject_procedure_value_into_session_state(
    session_state: Any,
    key: str,
    value: Any,
    *,
    widget_key: Callable[[str], str],
) -> None:
    """Set one ``fld_<key>`` (and helpers) from a JSON-loaded value."""
    sk = widget_key(key)

    if key in JSON_DICT_KEYS:
        if value is None or value == {}:
            session_state[sk] = '{}'
        elif isinstance(value, dict):
            session_state[sk] = json.dumps(value, indent=2)
        else:
            session_state[sk] = '{}'
        return

    if key in LIST_FLOAT_KEYS:
        if value is None:
            session_state[sk] = ''
        elif isinstance(value, (list, tuple)):
            session_state[sk] = ', '.join(str(float(x)) for x in value)
        else:
            session_state[sk] = str(value)
        return

    if key in LIST_INT_KEYS:
        if value is None:
            session_state[sk] = ''
        elif isinstance(value, (list, tuple)):
            session_state[sk] = ', '.join(str(int(round(float(x)))) for x in value)
        else:
            session_state[sk] = str(value)
        return

    if key == 'random_seed':
        use_k = f'{sk}_use'
        if value is None:
            session_state[use_k] = False
            session_state[sk] = 0
        else:
            session_state[use_k] = True
            try:
                session_state[sk] = int(value)
            except (TypeError, ValueError):
                session_state[use_k] = False
                session_state[sk] = 0
        return

    if key in ('isometric_random_seed', 'isovelocity_random_seed'):
        if value is None:
            session_state[sk] = ''
        else:
            session_state[sk] = str(int(value)) if isinstance(value, (int, float)) else str(value)
        return

    if isinstance(value, bool):
        session_state[sk] = value
        return

    if isinstance(value, float):
        session_state[sk] = float(value)
        return

    if isinstance(value, int) and not isinstance(value, bool):
        session_state[sk] = int(value)
        return

    if value is None:
        session_state[sk] = ''
        return

    session_state[sk] = str(value)


def apply_template_to_session_state(
    session_state: Any,
    template: Dict[str, Any],
    *,
    widget_key: Callable[[str], str],
    valid_test_types: List[str],
) -> Tuple[bool, str]:
    """
    Write template fields into Streamlit session_state (``fld_*`` + ``test_type_select``).

    Returns (ok, message). On failure, ok is False and message explains why.
    """
    ver = template.get('version', 1)
    try:
        ver = int(ver)
    except (TypeError, ValueError):
        ver = 0
    if ver not in (1, PROTOCOL_TEMPLATE_VERSION):
        return False, f'Unsupported template version {ver!r} (expected 1 or {PROTOCOL_TEMPLATE_VERSION}).'

    tt = str(template.get('test_type') or '').strip()
    if not tt:
        return False, 'Template is missing test_type.'
    if tt not in valid_test_types:
        return False, f'Template test_type {tt!r} is not in this app schema.'

    proc = template.get('procedure')
    if proc is None:
        proc = {}
    if not isinstance(proc, dict):
        return False, 'Template "procedure" must be a JSON object.'

    session_state['test_type_select'] = tt

    for k, v in proc.items():
        inject_procedure_value_into_session_state(
            session_state, str(k), v, widget_key=widget_key
        )

    bp = template.get('base_protocol')
    if isinstance(bp, dict) and tt == 'calibration':
        btt = str(bp.get('test_type') or '').strip()
        bproc = bp.get('procedure')
        if btt and btt in valid_test_types and isinstance(bproc, dict):
            for k, v in bproc.items():
                inject_procedure_value_into_session_state(
                    session_state, str(k), v, widget_key=widget_key
                )
            session_state['gui_calibration_embedded_base'] = {
                'test_type': btt,
                'procedure': dict(bproc),
            }
        else:
            session_state.pop('gui_calibration_embedded_base', None)
    else:
        session_state.pop('gui_calibration_embedded_base', None)

    return True, f'Loaded template {template.get("name")!r} ({tt}).'


def build_procedure_dict_from_updates(updates: Dict[str, Any]) -> Dict[str, Any]:
    """Copy procedure keys for JSON export."""
    return dict(updates)
