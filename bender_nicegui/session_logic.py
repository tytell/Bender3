"""Session-dict operations mirroring ``bender_streamlit_gui`` (no Streamlit)."""
from __future__ import annotations

import json
import math
import os
from typing import Any, Dict, Optional, TYPE_CHECKING

from .constants import BIO_DENSITY_PRESET_G_PER_MM3

if TYPE_CHECKING:
    from bender_functions import Bender


def compose_output_h5_path(sess: Dict[str, Any]) -> str:
    folder = str(sess.get('gui_data_folder') or '').strip()
    fn = str(sess.get('gui_data_filename') or '').strip()
    if not fn:
        return ''
    if not fn.lower().endswith('.h5'):
        fn = fn + '.h5'
    if folder:
        return os.path.normpath(os.path.join(folder, fn))
    return os.path.normpath(fn)


def section2_destination_incomplete(sess: Dict[str, Any]) -> bool:
    folder = str(sess.get('gui_data_folder') or '').strip()
    fn = str(sess.get('gui_data_filename') or '').strip()
    return (not folder) or (not fn)


def needs_missing_calibration_confirmation(b: Any) -> bool:
    use_cal = bool(getattr(b, 'use_inertial_calibration', False))
    cal_file = str(getattr(b, 'inertial_calibration_file', '') or '').strip()
    if not use_cal:
        return False
    return (not cal_file) or (not os.path.isfile(cal_file))


def read_qc_trial_index(sess: Dict[str, Any], n_trials: int) -> int:
    if n_trials <= 0:
        return 0
    try:
        ix = int(sess.get('gui_qc_trial_index'))
    except (TypeError, ValueError):
        ix = n_trials - 1
    return int(max(0, min(ix, n_trials - 1)))


def qc_figure_base_path(b: Any, selected_path: str, trial_idx: int) -> Optional[str]:
    p = str(selected_path or '').strip()
    if not p.lower().endswith('.h5'):
        return None
    proc = str(getattr(b, 'test_type', 'unknown'))
    stem = p[:-3]
    return stem + f'_{proc}_trial{int(trial_idx):03d}_qc'


def shared_experiment_dir(sess: Dict[str, Any], project_root: str) -> str:
    from bender_protocol_templates import default_templates_dir

    d = str(sess.get('gui_data_folder') or '').strip()
    if d:
        norm = os.path.normpath(d)
        if os.path.isdir(norm):
            return norm
    return os.path.normpath(default_templates_dir(project_root))


def sync_bio_prof_rho_from_density_preset(sess: Dict[str, Any]) -> None:
    label = str(sess.get('bio_prof_rho_preset') or '')
    v = BIO_DENSITY_PRESET_G_PER_MM3.get(label)
    if v is not None:
        sess['bio_prof_rho'] = float(v)


def sync_biometric_flags_from_session(b: Any, sess: Dict[str, Any]) -> None:
    if 'bio_use_theoretical_inertial' in sess:
        b.use_theoretical_inertial_correction = bool(sess['bio_use_theoretical_inertial'])
    if 'bio_dclamp' in sess:
        v = float(sess['bio_dclamp'])
        b.dclamp = v
        b.test_segment_length_mm = v
    if 'bio_dbend' in sess:
        v = float(sess['bio_dbend'])
        b.dbend = v
        b.test_segment_position_mm = v
    if 'bio_xsec' in sess:
        b.xsec_width = float(sess['bio_xsec'])
    if 'bio_temp_room' in sess:
        b.temp_C_room = float(sess['bio_temp_room'])
    if 'bio_temp_tank' in sess:
        b.temp_C_tank = float(sess['bio_temp_tank'])
    if 'bio_xsec_height' in sess:
        b.xsec_height = float(sess['bio_xsec_height'])
    if 'bio_dvert' in sess:
        b.dvert = float(sess['bio_dvert'])
    if 'bio_dhoriz' in sess:
        b.dhoriz = float(sess['bio_dhoriz'])
    if 'gui_specimen_id' in sess:
        b.fishcode = str(sess.get('gui_specimen_id') or '')
    if 'bio_segment' in sess:
        b.segment = str(sess['bio_segment'] or '')
    if 'bio_fishmass' in sess:
        b.fishmass = float(sess['bio_fishmass'])
    if 'bio_fishlen_TL' in sess:
        b.fishlen_TL = float(sess['bio_fishlen_TL'])
    if 'bio_fishlen_SL' in sess:
        b.fishlen_SL = float(sess['bio_fishlen_SL'])


def sync_genus_species_to_bender(b: Any, sess: Dict[str, Any]) -> None:
    meta = dict(getattr(b, 'h5_protocol_metadata', {}) or {})
    meta['genus_species'] = str(sess.get('gui_genus_species') or '').strip()
    sid = str(sess.get('gui_specimen_id') or '').strip()
    meta['specimen_id'] = sid
    setattr(b, 'fishcode', sid)

    def _str_attr(sk: str, bender_attr: str) -> None:
        if sk not in sess:
            return
        s = str(sess[sk] or '').strip()
        meta[bender_attr] = s
        setattr(b, bender_attr, s)

    def _float_attr(sk: str, bender_attr: str) -> None:
        if sk not in sess:
            return
        try:
            v = float(sess[sk])
        except (TypeError, ValueError):
            return
        meta[bender_attr] = v
        setattr(b, bender_attr, v)

    _str_attr('bio_segment', 'segment')
    _float_attr('bio_fishmass', 'fishmass')
    _float_attr('bio_fishlen_TL', 'fishlen_TL')
    _float_attr('bio_fishlen_SL', 'fishlen_SL')
    _float_attr('bio_xsec_height', 'xsec_height')
    if 'bio_prep_condition' in sess:
        meta['prep_condition'] = str(sess.get('bio_prep_condition') or '').strip()
    if 'bio_temp_room' in sess:
        try:
            meta['temp_C_room'] = float(sess['bio_temp_room'])
        except (TypeError, ValueError):
            pass
    if 'bio_temp_tank' in sess:
        try:
            meta['temp_C_tank'] = float(sess['bio_temp_tank'])
        except (TypeError, ValueError):
            pass
    if 'bio_dvert' in sess:
        try:
            meta['dvert'] = float(sess['bio_dvert'])
        except (TypeError, ValueError):
            pass
    if 'bio_dhoriz' in sess:
        try:
            meta['dhoriz'] = float(sess['bio_dhoriz'])
        except (TypeError, ValueError):
            pass

    b.h5_protocol_metadata = meta


def apply_specimen_identity_to_bender(b: Any, sess: Dict[str, Any]) -> None:
    meta = dict(getattr(b, 'h5_protocol_metadata', {}) or {})
    meta['genus_species'] = str(sess.get('gui_genus_species') or '').strip()
    sid = str(sess.get('gui_specimen_id') or '').strip()
    meta['specimen_id'] = sid
    b.fishcode = sid
    seg = str(sess.get('bio_segment') or '').strip()
    meta['segment'] = seg
    b.segment = seg
    b.h5_protocol_metadata = meta


def apply_pair(b: Any, name: str, value: Any) -> None:
    if value is None and name == 'random_seed':
        setattr(b, name, None)
        return
    if value is None:
        return
    if name == 'lateral_mode' and isinstance(value, str) and not value.strip():
        setattr(b, 'lateral_mode', None)
        return
    setattr(b, name, value)


def apply_form_updates(b: Any, updates: Dict[str, Any], tt: str, sess: Dict[str, Any]) -> None:
    if tt == 'calibration':
        emb = sess.get('gui_calibration_embedded_base')
        if isinstance(emb, dict):
            proc = emb.get('procedure')
            if isinstance(proc, dict):
                for k, v in proc.items():
                    apply_pair(b, k, v)
    for k, v in updates.items():
        apply_pair(b, k, v)
    b.test_type = tt


def apply_procedure_form_to_bender(b: Any, updates: Dict[str, Any], tt: str, sess: Dict[str, Any]) -> None:
    sync_biometric_flags_from_session(b, sess)
    apply_form_updates(b, updates, tt, sess)
    sync_genus_species_to_bender(b, sess)
    _pn = str(sess.get('gui_post_notes') or '').strip()
    if _pn:
        b.post_trial_notes = _pn


def apply_intrinsic_biometrics_to_bender(b: Any, sess: Dict[str, Any]) -> None:
    b.fishlen_TL = float(sess['bio_fishlen_TL'])
    b.fishlen_SL = float(sess['bio_fishlen_SL'])
    b.fishmass = float(sess['bio_fishmass'])
    sync_genus_species_to_bender(b, sess)


def apply_experimental_conditions_to_bender(b: Any, sess: Dict[str, Any]) -> None:
    b.temp_C_room = float(sess['bio_temp_room'])
    b.temp_C_tank = float(sess['bio_temp_tank'])
    meta = dict(getattr(b, 'h5_protocol_metadata', {}) or {})
    meta['temp_C_room'] = float(sess['bio_temp_room'])
    meta['temp_C_tank'] = float(sess['bio_temp_tank'])
    meta['prep_condition'] = str(sess.get('bio_prep_condition') or '').strip()
    b.h5_protocol_metadata = meta


def apply_clamp_geometry_to_bender(b: Any, sess: Dict[str, Any]) -> None:
    b.dclamp = float(sess['bio_dclamp'])
    b.test_segment_length_mm = float(sess['bio_dclamp'])
    b.dbend = float(sess['bio_dbend'])
    b.test_segment_position_mm = float(sess['bio_dbend'])
    b.xsec_width = float(sess['bio_xsec'])
    b.xsec_height = float(sess['bio_xsec_height'])
    b.dvert = float(sess['bio_dvert'])
    b.dhoriz = float(sess['bio_dhoriz'])
    meta = dict(getattr(b, 'h5_protocol_metadata', {}) or {})
    meta['dvert'] = float(sess['bio_dvert'])
    meta['dhoriz'] = float(sess['bio_dhoriz'])
    b.h5_protocol_metadata = meta


def apply_mounted_profile_inertial_to_bender(b: Any, sess: Dict[str, Any]) -> None:
    sync_bio_prof_rho_from_density_preset(sess)
    stations = b.make_profile_stations(
        sess['bio_prof_ph'],
        sess['bio_prof_pw'],
        sess['bio_prof_dh'],
        sess['bio_prof_dw'],
    )
    b.set_profiled_specimen_inertial_model(
        stations,
        sess['bio_prof_L'],
        sess['bio_prof_rho'],
        clamp_offset_mm=float(sess['bio_prof_clamp']),
        num_samples=int(sess['bio_prof_samples']),
    )


def apply_all_biometrics_to_bender(b: Any, sess: Dict[str, Any]) -> None:
    sync_biometric_flags_from_session(b, sess)
    apply_intrinsic_biometrics_to_bender(b, sess)
    apply_experimental_conditions_to_bender(b, sess)
    apply_clamp_geometry_to_bender(b, sess)
    apply_mounted_profile_inertial_to_bender(b, sess)


def init_biometrics_session_state(b: Any, sess: Dict[str, Any], *, force: bool = False) -> None:
    def _put(key: str, val: Any) -> None:
        if force or key not in sess:
            sess[key] = val

    dc = getattr(b, 'dclamp', None)
    xw = getattr(b, 'xsec_width', None)
    _meta_b = getattr(b, 'h5_protocol_metadata', {}) or {}
    _put(
        'gui_specimen_id',
        str(_meta_b.get('specimen_id') or getattr(b, 'fishcode', '') or '').strip(),
    )
    _put('bio_segment', str(getattr(b, 'segment', '') or ''))
    _fm = getattr(b, 'fishmass', None)
    _put('bio_fishmass', float(_fm) if _fm is not None and math.isfinite(float(_fm)) else 0.0)
    _ftl = getattr(b, 'fishlen_TL', None)
    _put('bio_fishlen_TL', float(_ftl) if _ftl is not None and math.isfinite(float(_ftl)) else 0.0)
    _fsl = getattr(b, 'fishlen_SL', None)
    _put('bio_fishlen_SL', float(_fsl) if _fsl is not None and math.isfinite(float(_fsl)) else 0.0)
    _xh = getattr(b, 'xsec_height', None)
    _put(
        'bio_xsec_height',
        float(_xh) if _xh is not None and math.isfinite(float(_xh)) else (float(xw) if xw is not None else 8.0),
    )
    _put('bio_dvert', float(getattr(b, 'dvert', 0.0) or 0.0))
    _put('bio_dhoriz', float(getattr(b, 'dhoriz', 0.0) or 0.0))
    _put('bio_dclamp', float(dc) if dc is not None else 10.0)
    _put('bio_xsec', float(xw) if xw is not None else 8.0)
    _put('bio_dbend', float(getattr(b, 'dbend', 0.0) or 0.0))
    _put('bio_temp_room', float(getattr(b, 'temp_C_room', 22.0) or 22.0))
    _put('bio_temp_tank', float(getattr(b, 'temp_C_tank', 22.0) or 22.0))
    _put('bio_prep_condition', str(_meta_b.get('prep_condition', '') or ''))
    _put(
        'bio_use_theoretical_inertial',
        bool(getattr(b, 'use_theoretical_inertial_correction', False)),
    )
    _put('bio_prof_L', float(getattr(b, 'specimen_profile_length_mm', 25.0) or 25.0))
    _put('bio_prof_rho', float(getattr(b, 'specimen_profile_density_g_per_mm3', 1.03e-3) or 1.03e-3))
    _put('bio_prof_ph', float(4.0))
    _put('bio_prof_pw', float(6.0))
    _put('bio_prof_dh', float(3.5))
    _put('bio_prof_dw', float(5.0))
    _put('bio_prof_clamp', float(getattr(b, 'specimen_profile_clamp_offset_mm', 20.0) or 20.0))
    _put('bio_prof_samples', int(getattr(b, 'specimen_profile_num_samples', 120) or 120))
    _put('bio_prof_rho_preset', 'Custom — edit the number below')


def clear_fld_session_keys(sess: Dict[str, Any]) -> None:
    for k in list(sess.keys()):
        if k.startswith('fld_'):
            del sess[k]


def ensure_gui_defaults(sess: Dict[str, Any]) -> None:
    sess.setdefault('gui_genus_species', '')
    sess.setdefault('gui_post_notes', '')
    sess.setdefault('gui_protocol_new_name', '')
    sess.setdefault('gui_protocol_new_desc', '')
    sess.setdefault('gui_protocol_overwrite', False)
    sess.setdefault('gui_biometrics_new_name', '')
    sess.setdefault('gui_biometrics_new_desc', '')
    sess.setdefault('gui_biometrics_overwrite', False)
    sess.setdefault('gui_qc_notes_append', True)
    sess.setdefault('gui_qc_trial_index', 0)
    sess.setdefault('gui_exp_hide', False)
    sess.setdefault('gui_sec4_hide', False)
    sess.setdefault('gui_sec5_hide', False)
    sess.setdefault('gui_show_preview', True)
    sess.setdefault('gui_preview_pts', 6000)
    sess.setdefault('gui_run_biometrics_confirm', False)
    sess.setdefault('gui_confirm_run_without_calibration', False)
    sess.setdefault('gui_confirm_run_without_destination', False)
    sess.setdefault('test_type_select', '')
    sess.setdefault('gui_calibration_embedded_base', None)
    sess.setdefault('gui_review_selected', '')
    sess.setdefault('gui_h5_explore_path', '')
    sess.setdefault('gui_h5_explore_trial', '')
    sess.setdefault('gui_h5_explore_x', 'Time (s)')
    sess.setdefault('gui_h5_explore_y', 'Primary torque corrected (N·m)')
    sess.setdefault('gui_h5_explore_loaded_path', '')
    sess.setdefault('gui_h5_explore_upload_path', '')
    sess.setdefault('gui_h5_explore_catalog', None)
    sess.setdefault('gui_h5_explore_notes', None)
    sess.setdefault('gui_h5_explore_schema', None)
    sess.setdefault('gui_run_daq_ok', False)
    sess.setdefault('gui_protocol_template_pick', '')
    sess.setdefault('gui_biometrics_template_pick', '')
