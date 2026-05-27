"""
Biometrics templates: JSON snapshots of section **3 · Biometrics** session keys (geometry, profile, conditions, flags).

Unlike protocol templates, these do not change ``test_type`` or procedure fields. Load fills
Streamlit session state; use **Apply** on geometry and profile blocks to copy onto ``Bender``.
"""
from __future__ import annotations

import json
import os
import re
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

BIOMETRICS_TEMPLATE_VERSION = 1
_BIOMETRICS_TEMPLATE_DIR = 'TemplateBiometrics'
_LEGACY_BIOMETRICS_TEMPLATE_DIR = 'BiometricsTemplates'

# Protocol templates (section 4) use this version and a ``procedure`` block — not biometrics.
_PROTOCOL_TEMPLATE_VERSION = 2

# Keys saved/loaded (Streamlit widget session keys).
BIOMETRICS_SESSION_KEYS = (
    'bio_segment',
    'bio_fishmass',
    'bio_fishlen_TL',
    'bio_fishlen_SL',
    'bio_xsec_height',
    'bio_dvert',
    'bio_dhoriz',
    'bio_dclamp',
    'bio_dbend',
    'bio_xsec',
    'bio_temp_room',
    'bio_temp_tank',
    'bio_prep_condition',
    'bio_prof_ph',
    'bio_prof_pw',
    'bio_prof_dh',
    'bio_prof_dw',
    'bio_prof_L',
    'bio_prof_rho_preset',
    'bio_prof_rho',
    'bio_prof_clamp',
    'bio_prof_samples',
    'bio_use_theoretical_inertial',
    'gui_genus_species',
    'gui_specimen_id',
)


def _coerce_biometrics_version(raw: Any) -> Optional[int]:
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
    """Heuristic: saved **procedure** JSON, not a biometrics snapshot."""
    if not isinstance(data, dict):
        return False
    if 'procedure' in data and 'test_type' in data:
        return True
    ver = _coerce_biometrics_version(data.get('version'))
    if ver == _PROTOCOL_TEMPLATE_VERSION:
        return True
    return False


def _extract_biometrics_session_dict(data: Dict[str, Any]) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
    """
    Return ``(session_dict, None)`` or ``(None, error_message)``.
    Supports normal ``{"session": {...}}``, legacy files without ``version``, and flat top-level ``bio_*`` keys.
    """
    sess = data.get('session')
    if isinstance(sess, dict):
        return sess, None
    flat = {k: data[k] for k in BIOMETRICS_SESSION_KEYS if k in data}
    if flat:
        return flat, None
    return None, 'Template missing a "session" object (and no biometrics fields like bio_dclamp at the top level).'


def default_biometrics_templates_dir(project_root: str) -> str:
    return os.path.normpath(os.path.join(project_root, _BIOMETRICS_TEMPLATE_DIR))


def sanitize_biometrics_filename_stem(name: str) -> str:
    raw = (name or '').strip()
    if not raw:
        return 'biometrics'
    s = re.sub(r'[^\w\s\-]', '_', raw, flags=re.UNICODE)
    s = re.sub(r'\s+', '_', s.strip())
    s = s.strip('_') or 'biometrics'
    return s[:120]


def biometrics_template_display_label(path: str) -> str:
    try:
        with open(path, encoding='utf-8-sig') as f:
            d = json.load(f)
        n = str(d.get('name') or '').strip()
        return f'{n}  ·  {os.path.basename(path)}' if n else os.path.basename(path)
    except Exception:
        return os.path.basename(path)


def list_biometrics_template_files(folder: str) -> List[str]:
    d = str(folder or '').strip()
    if not d:
        return []
    if not os.path.isdir(d):
        # Backward-compatible read path for repos still using the old folder name.
        if os.path.basename(os.path.normpath(d)) == _BIOMETRICS_TEMPLATE_DIR:
            legacy = os.path.join(os.path.dirname(os.path.normpath(d)), _LEGACY_BIOMETRICS_TEMPLATE_DIR)
            if os.path.isdir(legacy):
                d = legacy
            else:
                return []
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


def build_biometrics_payload(session_state: Any) -> Dict[str, Any]:
    out: Dict[str, Any] = {}
    for k in BIOMETRICS_SESSION_KEYS:
        if k in session_state:
            v = session_state[k]
            if isinstance(v, (bool, int, float, str)):
                out[k] = v
            else:
                out[k] = v
    return out


def save_biometrics_template(
    path: str,
    *,
    name: str,
    description: str,
    session_state: Any,
) -> None:
    payload = {
        'version': BIOMETRICS_TEMPLATE_VERSION,
        'name': str(name or '').strip(),
        'description': str(description or '').strip(),
        'saved_at': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
        'session': build_biometrics_payload(session_state),
    }
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(payload, f, indent=2)


def load_biometrics_template(path: str) -> Dict[str, Any]:
    # utf-8-sig strips a UTF-8 BOM if present (common for Notepad / Excel-exported JSON).
    with open(path, encoding='utf-8-sig') as f:
        return json.load(f)


def apply_biometrics_template_to_session(
    session_state: Any, data: Dict[str, Any]
) -> Tuple[bool, str]:
    if not isinstance(data, dict):
        return False, 'Biometrics file root must be a JSON object `{ ... }`, not a list or bare value.'

    if _looks_like_protocol_template(data):
        return (
            False,
            'This JSON is a **protocol** template (experiment type + procedure), not a biometrics file. '
            'Load it from **Protocol** using **Load template into form** (section 4), not from **Biometrics**. '
            'Tip: files like `Mandytest.json` with `"version": 2`, `"test_type"`, and `"procedure"` are protocol files. '
            'Biometrics files include `"version": 1` and a `"session"` object with keys like `bio_dclamp` and `bio_xsec`.',
        )

    sess, err = _extract_biometrics_session_dict(data)
    if err or sess is None:
        return False, err or 'Could not read biometrics fields from file.'

    raw_ver = _coerce_biometrics_version(data.get('version'))
    # Legacy exports without version, or version 0: accept if we have a session dict.
    if raw_ver is None or raw_ver == 0:
        effective_ver = BIOMETRICS_TEMPLATE_VERSION
    else:
        effective_ver = raw_ver

    if effective_ver != BIOMETRICS_TEMPLATE_VERSION:
        return (
            False,
            f'Unsupported biometrics template version {data.get("version")!r} (expected {BIOMETRICS_TEMPLATE_VERSION}). '
            'If this is a protocol file, load it from section 4 instead.',
        )

    n = 0
    for k, v in sess.items():
        if k in BIOMETRICS_SESSION_KEYS:
            session_state[k] = v
            n += 1
        elif k == 'bio_fishcode' and v is not None and str(v).strip():
            # Legacy templates: former fish code → single specimen ID field
            if not str(session_state.get('gui_specimen_id') or '').strip():
                session_state['gui_specimen_id'] = str(v).strip()
                n += 1
    if n == 0:
        return (
            False,
            'No recognized biometrics keys in the file. Expected at least one of: '
            + ', '.join(sorted(BIOMETRICS_SESSION_KEYS)[:6])
            + ', … inside "session".',
        )
    return True, f'Loaded {n} biometrics field(s) into the form. Use **Apply** in section 2 to update the experiment object.'
