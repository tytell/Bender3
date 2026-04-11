"""NiceGUI widgets for procedure fields."""
from __future__ import annotations

from typing import Any, Callable, Dict, List

from nicegui import ui

from .constants import (
    ALL_AMPS_MODE_OPTIONS,
    BILATERAL_MIRROR_LABEL,
    ISOMETRIC_FIELD_HELP,
    ISOVELOCITY_FIELD_HELP,
    ISOVELOCITY_WIDGET_LABEL,
    LATERAL_MODE_LABEL,
    MOTION_FIELD_HELP,
    MOTION_TYPES,
    RECRUITMENT_OPTIONS,
    format_strain_or_amp_mode,
    motion_parameter_rows,
)
from .procedure_ops import seed_procedure_fields, widget_key


def _label_block(text: str, cls: str) -> None:
    ui.markdown(text).classes(cls)


def render_procedure_fields(
    b: Any,
    tt: str,
    schema: Dict[str, Any],
    test_types: List[str],
    sess: Dict[str, Any],
    *,
    caption_class: str,
    label_class: str,
) -> None:
    seed_procedure_fields(b, tt, schema, test_types, sess)

    if tt == 'isometric':
        _label_block(
            '**Isometric** uses **section 3** segment length and width for strain geometry after **Apply**.',
            caption_class,
        )
        ui.markdown('**Required**').classes(label_class)
        for key in schema['isometric_required']:
            kind = 'float' if 'steps' not in key else 'int'
            lbl = key.replace('_', ' ')
            sk = widget_key(key)
            h = ISOMETRIC_FIELD_HELP.get(key)
            if kind == 'float':
                ui.number(lbl, format='%.6g').classes('w-full').bind_value(sess, sk).props(f'title="{h}"' if h else '')
            else:
                ui.number(lbl, format='%.0f').classes('w-full').bind_value(sess, sk).props(f'title="{h}"' if h else '')
        ui.markdown('**Optional**').classes(label_class + ' mt-2')
        for key in schema['isometric_optional']:
            sk = widget_key(key)
            h = ISOMETRIC_FIELD_HELP.get(key)
            if key in ('isometric_stim_params', 'isometric_stim_overrides'):
                ui.textarea(key.replace('_', ' '), value=sess.get(sk, '{}')).classes('w-full').bind_value(sess, sk).props(
                    'rows=4'
                )
            elif key == 'recruitment':
                ui.select(
                    {x: x for x in RECRUITMENT_OPTIONS},
                    label='recruitment',
                    value=sess[sk],
                ).classes('w-full').bind_value(sess, sk)
            elif key == 'lateral_mode':
                ui.input(LATERAL_MODE_LABEL).classes('w-full').bind_value(sess, sk)
            elif key in ('bilateral_mirror_motor',):
                ui.checkbox(BILATERAL_MIRROR_LABEL).bind_value(sess, sk)
            elif key == 'bilateral_sequential_left_frac':
                ui.number(key, format='%.6g').classes('w-full').bind_value(sess, sk)
            elif key == 'isometric_mode':
                opts = {m: format_strain_or_amp_mode(m) for m in ALL_AMPS_MODE_OPTIONS}
                ui.select(opts, label='isometric_mode').classes('w-full').bind_value(sess, sk)
            elif key == 'isometric_inter_step_interval_s':
                ui.number('Time between steps (s)', min=0.0, format='%.6g').classes('w-full').bind_value(sess, sk)
            elif 'random_seed' in key:
                ui.input('Random seed (optional)').classes('w-full').bind_value(sess, sk)
            else:
                if 'randomize' in key:
                    ui.checkbox(key.replace('_', ' ')).bind_value(sess, sk)
                else:
                    ui.input(key.replace('_', ' ')).classes('w-full').bind_value(sess, sk)

    elif tt == 'isovelocity':
        ui.markdown('**Required**').classes(label_class)
        for key in schema['isovelocity_required']:
            kind = 'int' if 'num_steps' in key else 'float'
            lbl = ISOVELOCITY_WIDGET_LABEL.get(key, key.replace('_', ' '))
            sk = widget_key(key)
            if kind == 'float':
                ui.number(lbl, format='%.6g').classes('w-full').bind_value(sess, sk)
            else:
                ui.number(lbl, format='%.0f').classes('w-full').bind_value(sess, sk)
        ui.markdown('**Optional**').classes(label_class + ' mt-2')
        for key in schema['isovelocity_optional']:
            sk = widget_key(key)
            if key in ('isovelocity_stim_params', 'isovelocity_stim_overrides'):
                ui.textarea(key.replace('_', ' '), value=sess.get(sk, '{}')).classes('w-full').bind_value(sess, sk).props(
                    'rows=4'
                )
            elif key == 'recruitment':
                ui.select({x: x for x in RECRUITMENT_OPTIONS}, label='recruitment').classes('w-full').bind_value(sess, sk)
            elif key == 'lateral_mode':
                ui.input(LATERAL_MODE_LABEL).classes('w-full').bind_value(sess, sk)
            elif key in ('bilateral_mirror_motor',):
                ui.checkbox(BILATERAL_MIRROR_LABEL).bind_value(sess, sk)
            elif key == 'bilateral_sequential_left_frac':
                ui.number(key, format='%.6g').classes('w-full').bind_value(sess, sk)
            elif key == 'isovelocity_starting_strain_mode':
                opts = {m: format_strain_or_amp_mode(m) for m in ALL_AMPS_MODE_OPTIONS}
                lbl = ISOVELOCITY_WIDGET_LABEL.get(key, key)
                ui.select(opts, label=lbl).classes('w-full').bind_value(sess, sk)
            elif 'random_seed' in key:
                ui.input('Random seed (optional)').classes('w-full').bind_value(sess, sk)
            else:
                kind = 'bool' if 'randomize' in key else 'float'
                lbl = ISOVELOCITY_WIDGET_LABEL.get(key, key.replace('_', ' '))
                if kind == 'bool':
                    ui.checkbox(lbl).bind_value(sess, sk)
                else:
                    ui.number(lbl, format='%.6g').classes('w-full').bind_value(sess, sk)

    elif tt == 'calibration':
        bases = [x for x in test_types if x != 'calibration']
        ui.markdown('**Required**').classes(label_class)
        sk_cal = widget_key('calibration_base_test_type')
        ui.select({x: x for x in bases}, label='calibration_base_test_type').classes('w-full').bind_value(sess, sk_cal)
        ui.markdown(
            'Calibration runs the **base** motion protocol. Set **experiment type** to that base, **Apply procedure**, '
            'then switch back to **calibration** before **Run**.'
        ).classes(caption_class)
        ui.markdown('**Optional**').classes(label_class + ' mt-2')
        for key in schema['calibration_optional']:
            ui.input(key).classes('w-full').bind_value(sess, widget_key(key))

    elif tt in MOTION_TYPES:
        ui.markdown('**Motion-series parameters**').classes(label_class)
        for name, kind, label in motion_parameter_rows(tt):
            sk = widget_key(name)
            h = MOTION_FIELD_HELP.get(name)
            if kind == 'float':
                ui.number(label, format='%.6g').classes('w-full').bind_value(sess, sk)
            elif kind == 'int':
                ui.number(label, format='%.0f').classes('w-full').bind_value(sess, sk)
            elif kind == 'optional_int':
                use_k = f'{sk}_use'
                with ui.row().classes('items-center gap-2 w-full'):
                    ui.checkbox(f'{label} (use custom)').bind_value(sess, use_k)
                ui.number(label, format='%.0f').classes('w-full').bind_value(sess, sk)
            elif kind == 'bool':
                ui.checkbox(label).bind_value(sess, sk)
            elif kind == 'select':
                opts = {m: format_strain_or_amp_mode(m) for m in ALL_AMPS_MODE_OPTIONS}
                ui.select(opts, label=label).classes('w-full').bind_value(sess, sk)
            elif kind == 'list_float':
                ui.input(label).classes('w-full').bind_value(sess, sk)
            elif kind == 'list_int':
                ui.input(label).classes('w-full').bind_value(sess, sk)
    else:
        ui.label(f'No field panel for {tt!r}').classes('text-amber-800')
