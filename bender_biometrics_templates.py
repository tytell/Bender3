"""
Biometrics templates: JSON snapshots of **section 2** session keys (geometry, profile, flags).

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
    'bio_prof_ph',
    'bio_prof_pw',
    'bio_prof_dh',
    'bio_prof_dw',
    'bio_prof_L',
    'bio_prof_rho',
    'bio_prof_clamp',
    'bio_prof_samples',
    'bio_use_theoretical_inertial',
    'gui_genus_species',
    'gui_specimen_id',
)


def default_biometrics_templates_dir(project_root: str) -> str:
    return os.path.normpath(os.path.join(project_root, 'BiometricsTemplates'))


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
        with open(path, encoding='utf-8') as f:
            d = json.load(f)
        n = str(d.get('name') or '').strip()
        return f'{n}  ·  {os.path.basename(path)}' if n else os.path.basename(path)
    except Exception:
        return os.path.basename(path)


def list_biometrics_template_files(folder: str) -> List[str]:
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
    with open(path, encoding='utf-8') as f:
        return json.load(f)


def apply_biometrics_template_to_session(
    session_state: Any, data: Dict[str, Any]
) -> Tuple[bool, str]:
    ver = data.get('version', 0)
    if ver != BIOMETRICS_TEMPLATE_VERSION:
        return False, f'Unsupported biometrics template version {ver!r} (expected {BIOMETRICS_TEMPLATE_VERSION}).'
    sess = data.get('session')
    if not isinstance(sess, dict):
        return False, 'Template missing "session" object.'
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
    return True, f'Loaded {n} biometrics field(s) into the form. Use **Apply** in section 2 to update the experiment object.'
