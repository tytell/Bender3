"""
Morphometrics templates: JSON snapshots of section **3 · Morphometrics** session keys (geometry, profile, conditions, flags).

Unlike protocol templates, these do not change ``test_type`` or procedure fields. Load fills
Streamlit session state; use **Apply** on geometry and profile blocks to copy onto ``Bender``.
"""
from __future__ import annotations

import json
import os
import re
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

from bender_json_persistent import JsonPersistTypeError, to_json_persistent

MORPHOMETRICS_TEMPLATE_VERSION = 1
_MORPHOMETRICS_TEMPLATE_DIR = os.path.join('templates', 'specimens')
# Older flat layouts read for backward compatibility (newest first).
_LEGACY_MORPHOMETRICS_TEMPLATE_DIRS = ('TemplateMorphometrics', 'MorphometricsTemplates')

# Protocol templates (section 4) use this version and a ``procedure`` block — not morphometrics.
_PROTOCOL_TEMPLATE_VERSION = 2

# Keys saved/loaded (Streamlit widget session keys).
MORPHOMETRICS_SESSION_KEYS = (
    'morpho_segment',
    'morpho_fishmass',
    'morpho_fishlen_TL',
    'morpho_fishlen_SL',
    'morpho_xsec_height',
    'morpho_dvert',
    'morpho_dhoriz',
    'morpho_dclamp',
    'morpho_dbend',
    'morpho_xsec',
    'morpho_muscle_depth',
    'morpho_temp_room',
    'morpho_temp_tank',
    'morpho_prep_condition',
    'morpho_prof_ph',
    'morpho_prof_pw',
    'morpho_prof_dh',
    'morpho_prof_dw',
    'morpho_prof_L',
    'morpho_prof_rho_preset',
    'morpho_prof_rho',
    'morpho_prof_clamp',
    'morpho_prof_samples',
    'morpho_use_theoretical_inertial',
    'gui_genus_species',
    'gui_specimen_id',
)


def _coerce_morphometrics_version(raw: Any) -> Optional[int]:
    """Accept ``1``, ``1.0``, ``"1"`` from JSON; return ``None`` if missing or unusable."""
    if raw is None or raw is False:
        return None
    if isinstance(raw, bool):
        return None
    if isinstance(raw, int):
        return int(raw)
    if isinstance(raw, float) and raw == raw and abs(raw - round(raw)) < 1e-9:
        return int(round(raw))
    if isinstance(raw, str):
        s = raw.strip()
        if s.isdigit() or (s.startswith('-') and s[1:].isdigit()):
            return int(s)
    try:
        return int(raw)
    except (TypeError, ValueError):
        return None


def _looks_like_protocol_template(data: Dict[str, Any]) -> bool:
    """Heuristic: saved **procedure** JSON, not a morphometrics snapshot."""
    if not isinstance(data, dict):
        return False
    if 'procedure' in data and 'test_type' in data:
        return True
    ver = _coerce_morphometrics_version(data.get('version'))
    if ver == _PROTOCOL_TEMPLATE_VERSION:
        return True
    return False


def _extract_morphometrics_session_dict(data: Dict[str, Any]) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
    """
    Return ``(session_dict, None)`` or ``(None, error_message)``.
    Supports normal ``{"session": {...}}``, legacy files without ``version``, and flat top-level ``morpho_*`` keys.
    """
    sess = data.get('session')
    if isinstance(sess, dict):
        return sess, None
    flat = {k: data[k] for k in MORPHOMETRICS_SESSION_KEYS if k in data}
    if flat:
        return flat, None
    return None, 'Template missing a "session" object (and no morphometrics fields like morpho_dclamp at the top level).'


def default_morphometrics_templates_dir(project_root: str) -> str:
    return os.path.normpath(os.path.join(project_root, _MORPHOMETRICS_TEMPLATE_DIR))


def _legacy_morphometrics_dir_for(canonical_dir: str) -> Optional[str]:
    """If ``canonical_dir`` is the templates/specimens dir, return an existing legacy dir (or None)."""
    norm = os.path.normpath(canonical_dir)
    suffix = os.path.normpath(_MORPHOMETRICS_TEMPLATE_DIR)
    if not (norm == suffix or norm.endswith(os.sep + suffix)):
        return None
    project_root = norm[: len(norm) - len(suffix)].rstrip(os.sep) or os.sep
    for legacy_name in _LEGACY_MORPHOMETRICS_TEMPLATE_DIRS:
        legacy = os.path.join(project_root, legacy_name)
        if os.path.isdir(legacy):
            return legacy
    return None


def sanitize_morphometrics_filename_stem(name: str) -> str:
    raw = (name or '').strip()
    if not raw:
        return 'morphometrics'
    s = re.sub(r'[^\w\s\-]', '_', raw, flags=re.UNICODE)
    s = re.sub(r'\s+', '_', s.strip())
    s = s.strip('_') or 'morphometrics'
    return s[:120]


def morphometrics_template_display_label(path: str) -> str:
    try:
        with open(path, encoding='utf-8-sig') as f:
            d = json.load(f)
        n = str(d.get('name') or '').strip()
        return f'{n}  ·  {os.path.basename(path)}' if n else os.path.basename(path)
    except Exception:
        return os.path.basename(path)


def list_morphometrics_template_files(folder: str) -> List[str]:
    d = str(folder or '').strip()
    if not d:
        return []
    if not os.path.isdir(d):
        # Backward-compatible read path for repos still using an old folder layout.
        legacy = _legacy_morphometrics_dir_for(d)
        if legacy is not None:
            d = legacy
        else:
            return []
    out = []
    for n in sorted(os.listdir(d)):
        if not n.lower().endswith('.json'):
            continue
        fp = os.path.join(d, n)
        if os.path.isfile(fp):
            out.append(fp)
    return out


def build_morphometrics_payload(session_state: Any) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for k in MORPHOMETRICS_SESSION_KEYS:
        if k in session_state:
            out[k] = to_json_persistent(session_state[k], path=k)
    return out


def save_morphometrics_template(
    path: str,
    *,
    name: str,
    description: str,
    session_state: Any,
) -> None:
    payload = {
        'version': MORPHOMETRICS_TEMPLATE_VERSION,
        'name': str(name or '').strip(),
        'description': str(description or '').strip(),
        'saved_at': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
        'session': build_morphometrics_payload(session_state),
    }
    try:
        json.dumps(payload, allow_nan=False)
    except (TypeError, ValueError) as e:
        raise JsonPersistTypeError(f'Morphometrics template not JSON-serializable: {e}') from e
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(payload, f, indent=2, allow_nan=False)


def load_morphometrics_template(path: str) -> Dict[str, Any]:
    # utf-8-sig strips a UTF-8 BOM if present (common for Notepad / Excel-exported JSON).
    with open(path, encoding='utf-8-sig') as f:
        return json.load(f)


def apply_morphometrics_template_to_session(
    session_state: Any, data: Dict[str, Any]
) -> Tuple[bool, str]:
    if not isinstance(data, dict):
        return False, 'Morphometrics file root must be a JSON object `{ ... }`, not a list or bare value.'

    if _looks_like_protocol_template(data):
        return (
            False,
            'This JSON is a **protocol** template (experiment type + procedure), not a morphometrics file. '
            'Load it from **Protocol** using **Load template into form** (section 4), not from **Morphometrics**. '
            'Tip: files like `Mandytest.json` with `"version": 2`, `"test_type"`, and `"procedure"` are protocol files. '
            'Morphometrics files include `"version": 1` and a `"session"` object with keys like `morpho_dclamp` and `morpho_xsec`.',
        )

    sess, err = _extract_morphometrics_session_dict(data)
    if err or sess is None:
        return False, err or 'Could not read morphometrics fields from file.'

    raw_ver = _coerce_morphometrics_version(data.get('version'))
    # Legacy exports without version, or version 0: accept if we have a session dict.
    if raw_ver is None or raw_ver == 0:
        effective_ver = MORPHOMETRICS_TEMPLATE_VERSION
    else:
        effective_ver = raw_ver

    if effective_ver != MORPHOMETRICS_TEMPLATE_VERSION:
        return (
            False,
            f'Unsupported morphometrics template version {data.get("version")!r} (expected {MORPHOMETRICS_TEMPLATE_VERSION}). '
            'If this is a protocol file, load it from section 4 instead.',
        )

    n = 0
    for k, v in sess.items():
        if k in MORPHOMETRICS_SESSION_KEYS:
            session_state[k] = v
            n += 1
        elif k == 'morpho_fishcode' and v is not None and str(v).strip():
            # Legacy templates: former fish code → single specimen ID field
            if not str(session_state.get('gui_specimen_id') or '').strip():
                session_state['gui_specimen_id'] = str(v).strip()
                n += 1
    if n == 0:
        return (
            False,
            'No recognized morphometrics keys in the file. Expected at least one of: '
            + ', '.join(sorted(MORPHOMETRICS_SESSION_KEYS)[:6])
            + ', … inside "session".',
        )
    return True, f'Loaded {n} morphometrics field(s) into the form. Use **Apply** in section 2 to update the experiment object.'
