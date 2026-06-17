"""
Helpers to discover config modules and emit a small generated ``.py`` that does
``from <base> import *`` plus user overrides (CritterGripper **Build** path).
"""
from __future__ import annotations

import importlib
import os
import re
import sys
from typing import Any, Dict, List, Optional

# Hardware config modules live here (relative to the project root) after the templates reorg.
_CONFIG_TEMPLATE_DIR = os.path.join('templates')

# Skip framework / app modules when scanning the project folder.
_SKIP_NAME_PREFIXES = (
    'bender_',
    'test_',
)
_SKIP_EXACT = frozenset(
    {
        'setup',
        'conftest',
    }
)


def default_configs_dir(project_root: str) -> str:
    """Canonical folder for hardware config modules: ``<project_root>/templates/configs``."""
    return os.path.normpath(os.path.join(str(project_root or ''), _CONFIG_TEMPLATE_DIR))


def discover_config_modules(project_root: str) -> List[str]:
    """Return importable module stems for likely hardware config files.

    Scans ``<project_root>/templates/configs`` (the canonical location after the templates
    reorg) and ensures that folder is on ``sys.path`` so the stems import by bare name.
    Falls back to scanning ``project_root`` itself for older flat layouts.
    """
    configs_dir = default_configs_dir(project_root)
    if os.path.isdir(configs_dir):
        if configs_dir not in sys.path:
            sys.path.insert(0, configs_dir)
        root = configs_dir
    else:
        root = os.path.abspath(str(project_root or ''))
    if not os.path.isdir(root):
        return ['jimenez_bender_config_A']
    out: List[str] = []
    for n in sorted(os.listdir(root)):
        if not n.endswith('.py'):
            continue
        stem = n[:-3]
        low = stem.lower()
        if low in _SKIP_EXACT:
            continue
        if any(low.startswith(p) for p in _SKIP_NAME_PREFIXES):
            continue
        if 'config' in low or stem.startswith('jimenez_bender'):
            out.append(stem)
    return out or ['jimenez_bender_config_A']


def sanitize_config_module_stem(name: str) -> str:
    raw = (name or '').strip()
    if not raw:
        return 'user_bender_config'
    s = re.sub(r'[^\w]', '_', raw, flags=re.UNICODE)
    s = re.sub(r'_+', '_', s).strip('_')
    if not s:
        return 'user_bender_config'
    if s[0].isdigit():
        s = 'cfg_' + s
    return s[:120]


def read_base_defaults(base_module: str) -> Dict[str, Any]:
    """Pull editable fields from an existing config module (must import)."""
    m = importlib.import_module(base_module)
    sono_l = list(getattr(m, 'sono_cal_left', [1.1, 4.5, 11.8, 47.0]))
    sono_r = list(getattr(m, 'sono_cal_right', [1.1, 4.5, 11.8, 47.0]))
    return {
        'apparatus_id': str(getattr(m, 'apparatus_id', '') or ''),
        'forcetorque_calibration_file': getattr(m, 'forcetorque_calibration_file', 'FT56491.cal'),
        'positive_motor_direction': getattr(m, 'positive_motor_direction', 'left'),
        'specimen_lateral_index_on_positive_motor_side': int(
            getattr(m, 'specimen_lateral_index_on_positive_motor_side', -1)
        ),
        'motor_axis': getattr(m, 'motor_axis', 'z'),
        'bending_axis_sensor': getattr(m, 'bending_axis_sensor', 'z'),
        'primary_bending_axis': getattr(m, 'primary_bending_axis', 'zTorque'),
        'bending_axis_specimen': getattr(m, 'bending_axis_specimen', 'dorsoventral'),
        'device_name': getattr(m, 'device_name', 'Dev1'),
        'daq_ai_sample_rate_hz': float(getattr(m, 'daq_ai_sample_rate_hz', 1000.0)),
        'daq_ao_do_sample_rate_hz': float(getattr(m, 'daq_ao_do_sample_rate_hz', 60000.0)),
        'motor_full_steps_per_rev': int(getattr(m, 'motor_full_steps_per_rev', 1600)),
        'motor_gear_ratio': int(getattr(m, 'motor_gear_ratio', 5)),
        'stim_channels': list(getattr(m, 'stim_channels', ['ao0', 'ao1'])),
        'motor_port': getattr(m, 'motor_port', 'port0'),
        'encoder_chan': getattr(m, 'encoder_chan', 'ctr0'),
        'SG_chan': list(getattr(m, 'SG_chan', ['ai0', 'ai1', 'ai2', 'ai3', 'ai4', 'ai5'])),
        'SG_name': list(
            getattr(m, 'SG_name', ['xForce', 'yForce', 'zForce', 'xTorque', 'yTorque', 'zTorque'])
        ),
        'stim_monitor_chan': list(getattr(m, 'stim_monitor_chan', [])),
        'stim_monitor_name': list(getattr(m, 'stim_monitor_name', [])),
        'S1side': getattr(m, 'S1side', 'left'),
        'S2side': getattr(m, 'S2side', 'right'),
        'use_sono': bool(getattr(m, 'use_sono', True)),
        'sono_channel': list(getattr(m, 'sono_channel', ['ai6', 'ai7'])),
        'sono_name': list(getattr(m, 'sono_name', ['sono_left', 'sono_right'])),
        'sono_internal_samplefreq': int(getattr(m, 'sono_internal_samplefreq', 241)),
        'sono_cal_left': sono_l,
        'sono_cal_right': sono_r,
        # Sono acquisition parameters. ``sono_distance`` is a comma-separated string (one value per
        # crystal pair, e.g. '12.5,14.2'), parsed downstream. Defaults keep old configs loadable.
        'sono_transmit_pulse': float(getattr(m, 'sono_transmit_pulse', 0.0)),
        'sono_inhibit_delay': float(getattr(m, 'sono_inhibit_delay', 0.0)),
        'sono_distance': str(getattr(m, 'sono_distance', '') or ''),
        'encoder_pulses_per_rev': int(getattr(m, 'encoder_pulses_per_rev', 10000)),
        'amp_step_vel': int(getattr(m, 'amp_step_vel', 10)),
        'ramp_mode_default': getattr(m, 'ramp_mode_default', 'linear'),
        'waitbefore': float(getattr(m, 'waitbefore', 3.0)),
        'waitafter': float(getattr(m, 'waitafter', 4.0)),
        'rampdur': float(getattr(m, 'rampdur', 0.25)),
        'prepoststim_dur': float(getattr(m, 'prepoststim_dur', 0.06)),
        'prepoststim_sep': float(getattr(m, 'prepoststim_sep', 1.0)),
        'prestim_time': float(getattr(m, 'prestim_time', -2.0)),
        'poststim_time': float(getattr(m, 'poststim_time', 2.0)),
    }


def render_generated_config(base_module: str, assignments: Dict[str, Any]) -> str:
    """Return ``.py`` source: import * from base, assignments, then rebuild input channel lists."""
    lines = [
        '# -*- coding: utf-8 -*-',
        f'"""Auto-generated by CritterGripper. Template: {base_module}. Edit freely after load."""',
        f'from {base_module} import *',
        '',
    ]
    for k in sorted(assignments.keys()):
        v = assignments[k]
        lines.append(f'{k} = {repr(v)}')
    lines.extend(
        [
            '',
            '# Combined AI lists after overrides (stim monitor appended when listed)',
            'input_channels = SG_chan + (sono_channel if use_sono else []) + list(stim_monitor_chan)',
            'input_channel_names = SG_name + (sono_name if use_sono else []) + list(stim_monitor_name)',
            '',
        ]
    )
    return '\n'.join(lines)


def parse_comma_list(s: str) -> List[str]:
    s = (s or '').strip()
    if not s:
        return []
    return [p.strip() for p in s.split(',') if p.strip()]


def parse_n_floats(s: str, n: int) -> List[float]:
    """Parse exactly ``n`` comma-separated floats (e.g. sono calibration quads)."""
    parts = parse_comma_list(s)
    if len(parts) != n:
        raise ValueError(f'need exactly {n} comma-separated numbers, got {len(parts)}')
    return [float(x) for x in parts]


def parse_sono_calibration(s: str) -> List[float]:
    """Parse a sonomicrometer calibration list.

    The list is grouped: the first half are voltage values and the second
    half are the matching mm values, i.e. ``[v_1, ..., v_N, mm_1, ..., mm_N]``
    for N >= 2 calibration points. Any even count of at least 4 numbers is
    accepted (>= 2 points). The 2-point case ``[v_lo, v_hi, mm_lo, mm_hi]``
    is unchanged. Fewer than 2 points is rejected.
    """
    parts = parse_comma_list(s)
    if len(parts) < 4 or len(parts) % 2 != 0:
        raise ValueError(
            'need an even count of at least 4 comma-separated numbers '
            '(>= 2 calibration points: all volts first, then all mm), '
            f'got {len(parts)}'
        )
    return [float(x) for x in parts]


def effective_load_module_name(*, typed: str, selected: str) -> Optional[str]:
    t = (typed or '').strip()
    if t:
        return t
    s = (selected or '').strip()
    return s or None
