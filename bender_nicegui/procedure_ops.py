"""Procedure field seeding and gathering (``fld_*`` session keys)."""
from __future__ import annotations

import json
from typing import Any, Dict, List, Optional, Tuple

from .constants import (
    ALL_AMPS_MODE_OPTIONS,
    BILATERAL_MIRROR_LABEL,
    ISOVELOCITY_WIDGET_LABEL,
    LATERAL_MODE_LABEL,
    MOTION_TYPES,
    RECRUITMENT_OPTIONS,
    motion_parameter_rows,
)


def widget_key(name: str) -> str:
    return f'fld_{name}'


def get_session_value(b: Any, name: str, default=None):
    if hasattr(b, name):
        v = getattr(b, name)
        if v is not None:
            return v
    return default


def parse_float_list(s: str):
    s = (s or '').strip()
    if not s:
        return None
    out = []
    for part in s.split(','):
        part = part.strip()
        if part:
            out.append(float(part))
    return out if out else None


def parse_int_list(s: str):
    s = (s or '').strip()
    if not s:
        return None
    out = []
    for part in s.split(','):
        part = part.strip()
        if part:
            out.append(int(round(float(part))))
    return out if out else None


def _seed_typed(b: Any, sess: Dict[str, Any], name: str, kind: str) -> None:
    sk = widget_key(name)
    cur = get_session_value(b, name)
    if kind == 'float':
        if sk not in sess:
            sess[sk] = float(cur) if cur is not None else 0.0
    elif kind == 'int':
        if sk not in sess:
            sess[sk] = int(cur) if cur is not None else 0
    elif kind == 'optional_int':
        use_k = f'{sk}_use'
        if use_k not in sess:
            sess[use_k] = cur is not None
        if sk not in sess:
            sess[sk] = int(cur) if cur is not None else 0
    elif kind == 'optional_float':
        use_k = f'{sk}_use'
        if use_k not in sess:
            sess[use_k] = cur is not None
        if sk not in sess:
            sess[sk] = float(cur) if cur is not None else 0.0
    elif kind == 'bool':
        if sk not in sess:
            sess[sk] = bool(cur) if cur is not None else False
    elif kind == 'str':
        if sk not in sess:
            sess[sk] = str(cur or '')
    elif kind == 'select':
        opts = list(ALL_AMPS_MODE_OPTIONS)
        dv = str(cur or 'strain')
        if dv not in opts:
            dv = 'strain'
        if sk not in sess:
            sess[sk] = dv
    elif kind == 'list_float':
        if sk not in sess:
            if cur is not None:
                try:
                    sess[sk] = ', '.join(str(float(x)) for x in list(cur))
                except (TypeError, ValueError):
                    sess[sk] = ''
            else:
                sess[sk] = ''
    elif kind == 'list_int':
        if sk not in sess:
            if cur is not None:
                try:
                    sess[sk] = ', '.join(str(int(x)) for x in list(cur))
                except (TypeError, ValueError):
                    sess[sk] = ''
            else:
                sess[sk] = ''
    elif kind == 'json_dict':
        if sk not in sess:
            sess[sk] = json.dumps(cur, indent=2) if isinstance(cur, dict) else '{}'


def _read_typed(sess: Dict[str, Any], name: str, kind: str) -> Any:
    sk = widget_key(name)
    if kind == 'float':
        return float(sess[sk])
    if kind == 'int':
        return int(round(float(sess[sk])))
    if kind == 'optional_int':
        use_k = f'{sk}_use'
        if not sess.get(use_k):
            return None
        return int(sess[sk])
    if kind == 'optional_float':
        use_k = f'{sk}_use'
        if not sess.get(use_k):
            return None
        return float(sess[sk])
    if kind == 'bool':
        return bool(sess[sk])
    if kind == 'str':
        return str(sess[sk])
    if kind == 'select':
        return str(sess[sk])
    if kind == 'list_float':
        return parse_float_list(str(sess[sk]))
    if kind == 'list_int':
        return parse_int_list(str(sess[sk]))
    if kind == 'json_dict':
        return json.loads(str(sess[sk]))
    return None


def seed_procedure_fields(b: Any, tt: str, schema: Dict[str, Any], test_types: List[str], sess: Dict[str, Any]) -> None:
    if tt == 'isometric':
        for key in schema['isometric_required']:
            kind = 'float' if 'steps' not in key else 'int'
            _seed_typed(b, sess, key, kind)
        for key in schema['isometric_optional']:
            if key in ('isometric_stim_params', 'isometric_stim_overrides'):
                _seed_typed(b, sess, key, 'json_dict')
            elif key in ('bilateral_mirror_motor',):
                _seed_typed(b, sess, key, 'bool')
            elif key in ('isometric_mirror_target_left', 'isometric_mirror_target_right'):
                _seed_typed(b, sess, key, 'optional_float')
            elif key == 'recruitment':
                skr = widget_key('recruitment')
                cur_r = get_session_value(b, 'recruitment', 'bilateral_simultaneous')
                if skr not in sess:
                    sess[skr] = cur_r if cur_r in RECRUITMENT_OPTIONS else RECRUITMENT_OPTIONS[0]
            elif key == 'lateral_mode':
                skl = widget_key('lateral_mode')
                if skl not in sess:
                    sess[skl] = str(get_session_value(b, key) or '')
            elif key == 'bilateral_sequential_left_frac':
                _seed_typed(b, sess, key, 'float')
            elif key == 'isometric_mode':
                skm = widget_key('isometric_mode')
                modes = list(ALL_AMPS_MODE_OPTIONS)
                cur_m = str(get_session_value(b, key, 'strain'))
                if skm not in sess:
                    sess[skm] = cur_m if cur_m in modes else 'strain'
            elif key == 'isometric_inter_step_interval_s':
                _seed_typed(b, sess, key, 'float')
            elif 'random_seed' in key:
                sks = widget_key(key)
                if sks not in sess:
                    v0 = get_session_value(b, key)
                    sess[sks] = '' if v0 is None else str(v0)
            else:
                kind = 'bool' if 'randomize' in key else 'str'
                _seed_typed(b, sess, key, kind)
    elif tt == 'isovelocity':
        for key in schema['isovelocity_required']:
            kind = 'int' if 'num_steps' in key else 'float'
            _seed_typed(b, sess, key, kind)
        for key in schema['isovelocity_optional']:
            if key in ('isovelocity_stim_params', 'isovelocity_stim_overrides'):
                _seed_typed(b, sess, key, 'json_dict')
            elif key in ('bilateral_mirror_motor',):
                _seed_typed(b, sess, key, 'bool')
            elif key == 'recruitment':
                skr = widget_key('recruitment')
                cur_r = get_session_value(b, 'recruitment', 'bilateral_simultaneous')
                if skr not in sess:
                    sess[skr] = cur_r if cur_r in RECRUITMENT_OPTIONS else RECRUITMENT_OPTIONS[0]
            elif key == 'lateral_mode':
                skl = widget_key('lateral_mode')
                if skl not in sess:
                    sess[skl] = str(get_session_value(b, key) or '')
            elif key == 'bilateral_sequential_left_frac':
                _seed_typed(b, sess, key, 'float')
            elif key == 'isovelocity_starting_strain_mode':
                skm = widget_key('isovelocity_starting_strain_mode')
                modes = list(ALL_AMPS_MODE_OPTIONS)
                cur_m = str(get_session_value(b, key, 'strain'))
                if skm not in sess:
                    sess[skm] = cur_m if cur_m in modes else 'strain'
            elif 'random_seed' in key:
                sks = widget_key(key)
                if sks not in sess:
                    v0 = get_session_value(b, key)
                    sess[sks] = '' if v0 is None else str(v0)
            else:
                kind = 'bool' if 'randomize' in key else 'float'
                _seed_typed(b, sess, key, kind)
    elif tt == 'calibration':
        bases = [x for x in test_types if x != 'calibration']
        cur = get_session_value(b, 'calibration_base_test_type', 'dynamic')
        if cur not in bases:
            cur = bases[0] if bases else 'dynamic'
        sk_cal = widget_key('calibration_base_test_type')
        if sk_cal not in sess:
            sess[sk_cal] = cur
        if sess[sk_cal] not in bases and bases:
            sess[sk_cal] = bases[0]
        for key in schema['calibration_optional']:
            sko = widget_key(key)
            if sko not in sess:
                sess[sko] = str(get_session_value(b, key) or '')
    elif tt in MOTION_TYPES:
        for name, kind, _label in motion_parameter_rows(tt):
            _seed_typed(b, sess, name, kind)


def gather_procedure_updates(
    b: Any, tt: str, schema: Dict[str, Any], test_types: List[str], sess: Dict[str, Any]
) -> Tuple[Dict[str, Any], Optional[str]]:
    updates: Dict[str, Any] = {}
    try:
        if tt == 'isometric':
            for key in schema['isometric_required']:
                kind = 'float' if 'steps' not in key else 'int'
                updates[key] = _read_typed(sess, key, kind)
            if 'isometric_num_steps' in updates and updates['isometric_num_steps'] is not None:
                updates['isometric_num_steps'] = int(updates['isometric_num_steps'])
            for key in schema['isometric_optional']:
                if key in ('isometric_stim_params', 'isometric_stim_overrides'):
                    updates[key] = _read_typed(sess, key, 'json_dict')
                elif key in ('bilateral_mirror_motor',):
                    updates[key] = _read_typed(sess, key, 'bool')
                elif key in ('isometric_mirror_target_left', 'isometric_mirror_target_right'):
                    updates[key] = _read_typed(sess, key, 'optional_float')
                elif key == 'recruitment':
                    updates[key] = str(sess[widget_key('recruitment')])
                elif key == 'lateral_mode':
                    updates[key] = str(sess[widget_key('lateral_mode')])
                elif key == 'bilateral_sequential_left_frac':
                    updates[key] = _read_typed(sess, key, 'float')
                elif key == 'isometric_mode':
                    updates[key] = str(sess[widget_key('isometric_mode')])
                elif key == 'isometric_inter_step_interval_s':
                    updates[key] = _read_typed(sess, key, 'float')
                elif 'random_seed' in key:
                    s = str(sess.get(widget_key(key), '')).strip()
                    updates[key] = None if not s else int(s)
                else:
                    kind = 'bool' if 'randomize' in key else 'str'
                    updates[key] = _read_typed(sess, key, kind)
        elif tt == 'isovelocity':
            for key in schema['isovelocity_required']:
                kind = 'int' if 'num_steps' in key else 'float'
                updates[key] = _read_typed(sess, key, kind)
            if 'isovelocity_num_steps' in updates and updates['isovelocity_num_steps'] is not None:
                updates['isovelocity_num_steps'] = int(updates['isovelocity_num_steps'])
            for key in schema['isovelocity_optional']:
                if key in ('isovelocity_stim_params', 'isovelocity_stim_overrides'):
                    updates[key] = _read_typed(sess, key, 'json_dict')
                elif key in ('bilateral_mirror_motor',):
                    updates[key] = _read_typed(sess, key, 'bool')
                elif key == 'recruitment':
                    updates[key] = str(sess[widget_key('recruitment')])
                elif key == 'lateral_mode':
                    updates[key] = str(sess[widget_key('lateral_mode')])
                elif key == 'bilateral_sequential_left_frac':
                    updates[key] = _read_typed(sess, key, 'float')
                elif key == 'isovelocity_starting_strain_mode':
                    updates[key] = str(sess[widget_key('isovelocity_starting_strain_mode')])
                elif 'random_seed' in key:
                    s = str(sess.get(widget_key(key), '')).strip()
                    updates[key] = None if not s else int(s)
                else:
                    kind = 'bool' if 'randomize' in key else 'float'
                    updates[key] = _read_typed(sess, key, kind)
        elif tt == 'calibration':
            updates['calibration_base_test_type'] = str(sess[widget_key('calibration_base_test_type')])
            for key in schema['calibration_optional']:
                updates[key] = str(sess[widget_key(key)])
        elif tt in MOTION_TYPES:
            for name, kind, _lbl in motion_parameter_rows(tt):
                updates[name] = _read_typed(sess, name, kind)
        else:
            return {}, f'No dedicated field panel for {tt!r} yet.'
    except (json.JSONDecodeError, KeyError, TypeError, ValueError) as e:
        return {}, str(e)
    return updates, None
