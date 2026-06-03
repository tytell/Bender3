"""
CritterGripper — Streamlit UI for Bender experiment dispatch.

Run from the project directory:
    streamlit run bender_streamlit_gui.py

Select ``test_type`` first; edit fields, use **Apply** to copy them onto the
``Bender`` instance, optionally **save / load protocol** and **biometrics** files in **section 3**,
**Check required fields** to validate,
**Refresh experiment preview** (in **Procedure fields**) for a table/plot of commanded motion (no DAQ), then **Run experiment**
when hardware is ready.
"""
from __future__ import annotations

import base64
from datetime import datetime
import hashlib
import importlib
import json
import math
import ntpath
import os
import posixpath
import sys
import tempfile
import time
from typing import Any, Optional, Tuple

import h5py
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st

# Project root on path for `bender_functions` and config modules
_ROOT = os.path.dirname(os.path.abspath(__file__))
_LOGO_PATH = os.path.join(_ROOT, 'assets', 'jimenez_biomechanics_logo.png')
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

_DEBUG_LOG_PATH = os.path.join(_ROOT, 'debug-4506cd.log')


def _agent_debug_log(*, hypothesis_id: str, location: str, message: str, data: dict) -> None:
    # #region agent log
    try:
        import json as _json

        _payload = {
            'sessionId': '4506cd',
            'hypothesisId': hypothesis_id,
            'location': location,
            'message': message,
            'data': data,
            'timestamp': int(time.time() * 1000),
        }
        _safe_data = to_json_persistent(data, path='data')
        _payload['data'] = _safe_data
        with open(_DEBUG_LOG_PATH, 'a', encoding='utf-8') as _df:
            _df.write(_json.dumps(_payload) + '\n')
    except Exception:
        pass
    # #endregion


def _nsf_logo_path() -> Optional[str]:
    for name in ('nsf_logo.png', 'nsf_logo.svg'):
        p = os.path.join(_ROOT, 'assets', name)
        if os.path.isfile(p):
            return p
    return None


def _img_data_uri(path: str) -> str:
    ext = os.path.splitext(path)[1].lower()
    with open(path, 'rb') as f:
        raw = f.read()
    b64 = base64.b64encode(raw).decode('ascii')
    mime = 'image/svg+xml' if ext == '.svg' else 'image/png'
    return f'data:{mime};base64,{b64}'


from bender_json_persistent import JsonPersistTypeError, to_json_persistent  # noqa: E402
from bender_daq_kill import daq_emergency_stop  # noqa: E402
from bender_config_builder import (  # noqa: E402
    default_configs_dir,
    discover_config_modules,
    parse_comma_list,
    parse_n_floats,
    read_base_defaults,
    render_generated_config,
    sanitize_config_module_stem,
)
from bender_functions import Bender, _strain_lever_arm_m  # noqa: E402
from bender_gui_preview import build_protocol_preview  # noqa: E402
from bender_simulation import (  # noqa: E402
    cantilever_stiffness_N_per_m_for_geometry,
    force_displacement_comparison_curve,
    force_displacement_series,
    oscillatory_viscoelastic_timeseries,
    static_stiffness_comparison_delta_grid,
)
from bender_biometrics_templates import (  # noqa: E402
    apply_biometrics_template_to_session,
    biometrics_template_display_label,
    default_biometrics_templates_dir,
    list_biometrics_template_files,
    load_biometrics_template,
    save_biometrics_template,
    sanitize_biometrics_filename_stem,
)
from bender_protocol_templates import (  # noqa: E402
    apply_template_to_session_state,
    build_procedure_dict_from_updates,
    default_templates_dir,
    list_template_files,
    load_protocol_template,
    save_protocol_template,
    sanitize_template_filename_stem,
    snapshot_bender_procedure,
    template_display_label,
)
from bender_h5_explore import (  # noqa: E402
    FT_ROW_LABELS,
    align_xy as align_h5_catalog_xy,
    build_series_catalog_generic,
    build_series_catalog_legacy,
    build_series_catalog_v2,
    detect_h5_schema,
    h5_join_internal_path,
    list_h5_attribute_rows,
    list_h5_group_children,
    list_v2_trials,
    read_h5_series_1d,
    write_h5_user_attributes,
)

_BND_LS_ACTION_MARK = '<div class="bnd-ls-action" aria-hidden="true"></div>'

# App chrome presets (CSS markers: `body:has(.bnd-theme-*)` in `_inject_accessibility_theme`).
GUI_UI_THEME_OPTIONS = (
    'Default (Streamlit)',
    'Warm paper',
    'Cool gray',
    'Slate & ivory',
)
_GUI_THEME_TO_MARKER = {
    'Default (Streamlit)': '',
    'Warm paper': 'bnd-theme-warm',
    'Cool gray': 'bnd-theme-cool',
    'Slate & ivory': 'bnd-theme-slateivory',
}


def _migrate_gui_ui_theme_session() -> None:
    """Map removed options and unknown values onto ``GUI_UI_THEME_OPTIONS``."""
    v = str(st.session_state.get('gui_ui_theme', '') or '').strip()
    if v in ('Default', 'High contrast'):
        st.session_state['gui_ui_theme'] = GUI_UI_THEME_OPTIONS[0]
        return
    if v and v not in GUI_UI_THEME_OPTIONS:
        st.session_state['gui_ui_theme'] = GUI_UI_THEME_OPTIONS[0]


def _render_display_preferences_sidebar() -> None:
    """Theme and text size inside **Settings** (persists in session; applies via `_inject_accessibility_theme`)."""
    st.session_state.setdefault('gui_ui_theme', GUI_UI_THEME_OPTIONS[0])
    _migrate_gui_ui_theme_session()
    st.session_state.setdefault('gui_ui_large_text', False)
    st.markdown('**Display**')
    st.selectbox(
        'Theme',
        options=list(GUI_UI_THEME_OPTIONS),
        key='gui_ui_theme',
        help=(
            '**Default** uses Streamlit’s built-in look. Other options tint the app background, sidebar, and bordered panels '
            'for different reading preferences.'
        ),
    )
    st.checkbox(
        'Larger text',
        key='gui_ui_large_text',
        help='Slightly increases base font size in the main panel and sidebar.',
    )


def _render_sidebar_settings_expander(*, leading_divider: bool = True) -> None:
    """Bottom of sidebar: collapsed **Settings** panel (display / accessibility)."""
    if leading_divider:
        st.divider()
    with st.expander('Settings', expanded=False):
        st.caption('Theme and readability for this browser session.')
        _render_display_preferences_sidebar()


def _inject_accessibility_theme() -> None:
    """Skip link, optional theme tint and large-text (`:has` markers must be in DOM)."""
    _migrate_gui_ui_theme_session()
    theme_label = str(st.session_state.get('gui_ui_theme', GUI_UI_THEME_OPTIONS[0]) or '')
    theme_cls = _GUI_THEME_TO_MARKER.get(theme_label, '')
    lt = bool(st.session_state.get('gui_ui_large_text'))
    markers = []
    if theme_cls:
        markers.append(f'<div class="bnd-ui-theme {theme_cls}" aria-hidden="true"></div>')
    if lt:
        markers.append('<div class="bnd-a11y-large-text" aria-hidden="true"></div>')
    marker_html = ''.join(markers)
    st.markdown(
        f'''
<a href="#bnd-main-content" class="bnd-skip-link">Skip to main content</a>
{marker_html}
<style>
.bnd-skip-link {{
  position: absolute;
  left: -9999px;
  top: 0;
  z-index: 100000;
  padding: 0.5rem 1rem;
  background: #000000;
  color: #ffffff !important;
  font-weight: 600;
  text-decoration: underline;
  border-radius: 4px;
}}
.bnd-skip-link:focus {{
  left: 0.75rem;
  top: 0.75rem;
  width: auto;
  height: auto;
  overflow: visible;
  outline: 3px solid #fbbf24;
  outline-offset: 2px;
}}
/* Focus visibility */
[data-testid="stMain"] button:focus-visible,
[data-testid="stMain"] a:focus-visible,
[data-testid="stSidebar"] button:focus-visible,
[data-testid="stSidebar"] a:focus-visible {{
  outline: 3px solid #2563eb !important;
  outline-offset: 2px !important;
}}
/* —— Warm paper —— */
body:has(.bnd-theme-warm) .stApp {{
  background-color: #faf8f5 !important;
}}
body:has(.bnd-theme-warm) [data-testid="stSidebar"] {{
  background-color: #f3eee6 !important;
  border-right: 1px solid #d6cfc4 !important;
}}
body:has(.bnd-theme-warm) [data-testid="stMain"] .block-container {{
  background-color: #fffefb !important;
}}
body:has(.bnd-theme-warm) [data-testid="stMain"] p,
body:has(.bnd-theme-warm) [data-testid="stMain"] li,
body:has(.bnd-theme-warm) [data-testid="stMain"] h1,
body:has(.bnd-theme-warm) [data-testid="stMain"] h2,
body:has(.bnd-theme-warm) [data-testid="stMain"] h3,
body:has(.bnd-theme-warm) [data-testid="stCaption"],
body:has(.bnd-theme-warm) [data-testid="stWidgetLabel"] p {{
  color: #1c1917 !important;
}}
body:has(.bnd-theme-warm) [data-testid="stSidebar"] p,
body:has(.bnd-theme-warm) [data-testid="stSidebar"] span,
body:has(.bnd-theme-warm) [data-testid="stSidebar"] label {{
  color: #292524 !important;
}}
/* —— Cool gray —— */
body:has(.bnd-theme-cool) .stApp {{
  background-color: #f1f5f9 !important;
}}
body:has(.bnd-theme-cool) [data-testid="stSidebar"] {{
  background-color: #e2e8f0 !important;
  border-right: 1px solid #cbd5e1 !important;
}}
body:has(.bnd-theme-cool) [data-testid="stMain"] .block-container {{
  background-color: #f8fafc !important;
}}
body:has(.bnd-theme-cool) [data-testid="stMain"] p,
body:has(.bnd-theme-cool) [data-testid="stMain"] li,
body:has(.bnd-theme-cool) [data-testid="stMain"] h1,
body:has(.bnd-theme-cool) [data-testid="stMain"] h2,
body:has(.bnd-theme-cool) [data-testid="stMain"] h3,
body:has(.bnd-theme-cool) [data-testid="stCaption"],
body:has(.bnd-theme-cool) [data-testid="stWidgetLabel"] p {{
  color: #0f172a !important;
}}
body:has(.bnd-theme-cool) [data-testid="stSidebar"] p,
body:has(.bnd-theme-cool) [data-testid="stSidebar"] span,
body:has(.bnd-theme-cool) [data-testid="stSidebar"] label {{
  color: #1e293b !important;
}}
/* —— Slate sidebar + ivory main —— */
body:has(.bnd-theme-slateivory) .stApp {{
  background-color: #fefce8 !important;
}}
body:has(.bnd-theme-slateivory) [data-testid="stSidebar"] {{
  background-color: #1e293b !important;
  border-right: 1px solid #0f172a !important;
}}
body:has(.bnd-theme-slateivory) [data-testid="stMain"] .block-container {{
  background-color: #fffff7 !important;
}}
body:has(.bnd-theme-slateivory) [data-testid="stMain"] p,
body:has(.bnd-theme-slateivory) [data-testid="stMain"] li,
body:has(.bnd-theme-slateivory) [data-testid="stMain"] h1,
body:has(.bnd-theme-slateivory) [data-testid="stMain"] h2,
body:has(.bnd-theme-slateivory) [data-testid="stMain"] h3,
body:has(.bnd-theme-slateivory) [data-testid="stCaption"],
body:has(.bnd-theme-slateivory) [data-testid="stWidgetLabel"] p {{
  color: #0f172a !important;
}}
body:has(.bnd-theme-slateivory) [data-testid="stSidebar"] p,
body:has(.bnd-theme-slateivory) [data-testid="stSidebar"] span,
body:has(.bnd-theme-slateivory) [data-testid="stSidebar"] label {{
  color: #e2e8f0 !important;
}}
/* —— Larger text —— */
body:has(.bnd-a11y-large-text) [data-testid="stMain"] .block-container {{
  font-size: 1.125rem !important;
}}
body:has(.bnd-a11y-large-text) [data-testid="stSidebar"] {{
  font-size: 1.0625rem !important;
}}
/* —— Workflow: main column white (all themes); logo strip matches —— */
body:has(.bnd-workflow-active) [data-testid="stMain"],
body:has(.bnd-workflow-active) [data-testid="stMain"] .block-container {{
  background-color: #ffffff !important;
}}
body:has(.bnd-workflow-active) [data-testid="stMain"] [data-testid="stImage"],
body:has(.bnd-workflow-active) [data-testid="stMain"] [data-testid="stImage"] > div {{
  background-color: #ffffff !important;
}}
/* Light touch: main panel labels and subheaders read a bit clearer without heavy chrome */
body:has(.bnd-workflow-active) [data-testid="stMain"] [data-testid="stWidgetLabel"] p {{
  color: #334155 !important;
}}
body:has(.bnd-workflow-active):has(.bnd-theme-warm) [data-testid="stMain"] [data-testid="stWidgetLabel"] p {{
  color: #44403c !important;
}}
body:has(.bnd-workflow-active) [data-testid="stMain"] h2,
body:has(.bnd-workflow-active) [data-testid="stMain"] h3 {{
  color: #334155 !important;
}}
body:has(.bnd-workflow-active):has(.bnd-theme-warm) [data-testid="stMain"] h2,
body:has(.bnd-workflow-active):has(.bnd-theme-warm) [data-testid="stMain"] h3 {{
  color: #44403c !important;
}}
</style>
''',
        unsafe_allow_html=True,
    )


def _inject_load_save_button_theme() -> None:
    """Inject CSS for persistent action buttons.

    Visual tiers (main panel, workflow):
    - **Primary:** red fill, white label; hover brightens and adds a soft red glow.
    - **Secondary:** white fill, slate border and text; hover tints toward light red.
    - **Load/Save:** `_load_save_button` uses primary styling (red).
    """
    st.markdown(
        """
<style>
/* Workflow: white shell (overrides Streamlit default gray and optional theme tints) */
body:has(.bnd-workflow-active) .stApp,
body:has(.bnd-workflow-active) [data-testid="stAppViewContainer"],
body:has(.bnd-workflow-active) [data-testid="stHeader"] {
    background-color: #ffffff !important;
}
body:has(.bnd-workflow-active) [data-testid="stSidebar"] {
    background-color: #ffffff !important;
    border-right: 1px solid #e2e8f0 !important;
}
/* Sidebar — primary: red / white; secondary: white / slate (hover lights up) */
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="baseButton-primary"],
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="stBaseButton-primary"] {
    background-color: #dc2626 !important;
    background-image: none !important;
    color: #ffffff !important;
    border: 1px solid #b91c1c !important;
    font-weight: 600 !important;
    transition: background-color 0.15s ease, border-color 0.15s ease, box-shadow 0.15s ease !important;
}
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="baseButton-primary"] p,
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="baseButton-primary"] span,
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="stBaseButton-primary"] p,
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="stBaseButton-primary"] span {
    color: #ffffff !important;
}
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="baseButton-primary"]:hover,
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="stBaseButton-primary"]:hover {
    background-color: #f87171 !important;
    border-color: #fca5a5 !important;
    box-shadow: 0 0 0 3px rgba(248, 113, 113, 0.45) !important;
    background-image: none !important;
    filter: none !important;
}
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="baseButton-secondary"],
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="stBaseButton-secondary"] {
    background-color: #ffffff !important;
    background-image: none !important;
    color: #334155 !important;
    border: 1px solid #cbd5e1 !important;
    font-weight: 600 !important;
    transition: background-color 0.15s ease, border-color 0.15s ease, box-shadow 0.15s ease, color 0.15s ease !important;
}
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="baseButton-secondary"] p,
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="baseButton-secondary"] span,
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="stBaseButton-secondary"] p,
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="stBaseButton-secondary"] span {
    color: #334155 !important;
}
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="baseButton-secondary"]:hover,
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="stBaseButton-secondary"]:hover {
    background-color: #fff1f2 !important;
    border-color: #fb7185 !important;
    color: #991b1b !important;
    box-shadow: 0 2px 10px rgba(220, 38, 38, 0.18) !important;
}
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="baseButton-secondary"]:hover p,
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="baseButton-secondary"]:hover span,
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="stBaseButton-secondary"]:hover p,
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="stBaseButton-secondary"]:hover span {
    color: #991b1b !important;
}
/* Main panel — match Streamlit primary / secondary */
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-primary"],
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="baseButton-primary"],
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-primary"],
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-primary"],
body:has(.bnd-workflow-active) .stApp div.row-widget.stButton[data-testid="stButton"] button[data-testid="baseButton-primary"],
body:has(.bnd-workflow-active) .stApp div.row-widget.stButton[data-testid="stButton"] button[data-testid="stBaseButton-primary"],
body:has(.bnd-workflow-active) div.row-widget.stButton[data-testid="stButton"] button[data-testid="baseButton-primary"],
body:has(.bnd-workflow-active) div.row-widget.stButton[data-testid="stButton"] button[data-testid="stBaseButton-primary"] {
    background-color: #dc2626 !important;
    background-image: none !important;
    color: #ffffff !important;
    border: 1px solid #b91c1c !important;
    font-weight: 600 !important;
    transition: background-color 0.15s ease, border-color 0.15s ease, box-shadow 0.15s ease !important;
}
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-primary"] p,
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-primary"] span,
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-primary"] p,
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-primary"] span,
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="baseButton-primary"] p,
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="baseButton-primary"] span,
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-primary"] p,
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-primary"] span {
    color: #ffffff !important;
}
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-primary"]:hover,
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-primary"]:hover,
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="baseButton-primary"]:hover,
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-primary"]:hover {
    background-color: #f87171 !important;
    border-color: #fca5a5 !important;
    box-shadow: 0 0 0 3px rgba(248, 113, 113, 0.45) !important;
    color: #ffffff !important;
    background-image: none !important;
    filter: none !important;
}
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-secondary"],
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="baseButton-secondary"],
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-secondary"],
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-secondary"],
body:has(.bnd-workflow-active) .stApp div.row-widget.stButton[data-testid="stButton"] button[data-testid="baseButton-secondary"],
body:has(.bnd-workflow-active) .stApp div.row-widget.stButton[data-testid="stButton"] button[data-testid="stBaseButton-secondary"],
body:has(.bnd-workflow-active) div.row-widget.stButton[data-testid="stButton"] button[data-testid="baseButton-secondary"],
body:has(.bnd-workflow-active) div.row-widget.stButton[data-testid="stButton"] button[data-testid="stBaseButton-secondary"] {
    background-color: #ffffff !important;
    background-image: none !important;
    color: #334155 !important;
    border: 1px solid #cbd5e1 !important;
    font-weight: 600 !important;
    transition: background-color 0.15s ease, border-color 0.15s ease, box-shadow 0.15s ease, color 0.15s ease !important;
}
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-secondary"] p,
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-secondary"] span,
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-secondary"] p,
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-secondary"] span,
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="baseButton-secondary"] p,
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="baseButton-secondary"] span,
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-secondary"] p,
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-secondary"] span {
    color: #334155 !important;
}
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-secondary"]:hover,
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-secondary"]:hover,
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="baseButton-secondary"]:hover,
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-secondary"]:hover {
    background-color: #fff1f2 !important;
    border-color: #fb7185 !important;
    color: #991b1b !important;
    box-shadow: 0 2px 10px rgba(220, 38, 38, 0.18) !important;
}
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-secondary"]:hover p,
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-secondary"]:hover span,
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-secondary"]:hover p,
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-secondary"]:hover span,
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="baseButton-secondary"]:hover p,
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="baseButton-secondary"]:hover span,
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-secondary"]:hover p,
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-secondary"]:hover span {
    color: #991b1b !important;
}
/* role=button fallbacks (match primary) */
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] [role="button"] {
    background-color: #dc2626 !important;
    background-image: none !important;
    color: #ffffff !important;
    border: 1px solid #b91c1c !important;
    font-weight: 600 !important;
    transition: background-color 0.15s ease, border-color 0.15s ease, box-shadow 0.15s ease !important;
}
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] [role="button"] p,
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] [role="button"] span {
    color: #ffffff !important;
}
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] [role="button"]:hover {
    background-color: #f87171 !important;
    border-color: #fca5a5 !important;
    box-shadow: 0 0 0 3px rgba(248, 113, 113, 0.45) !important;
    background-image: none !important;
    filter: none !important;
}
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="baseButton-primary"],
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="stBaseButton-primary"] {
    background-color: #dc2626 !important;
    background-image: none !important;
    color: #ffffff !important;
    border: 1px solid #b91c1c !important;
    font-weight: 600 !important;
    transition: background-color 0.15s ease, border-color 0.15s ease, box-shadow 0.15s ease !important;
}
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="baseButton-primary"] p,
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="baseButton-primary"] span,
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="stBaseButton-primary"] p,
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="stBaseButton-primary"] span {
    color: #ffffff !important;
}
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="baseButton-primary"]:hover,
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="stBaseButton-primary"]:hover {
    background-color: #f87171 !important;
    border-color: #fca5a5 !important;
    box-shadow: 0 0 0 3px rgba(248, 113, 113, 0.45) !important;
    background-image: none !important;
    filter: none !important;
}
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="baseButton-secondary"],
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="stBaseButton-secondary"] {
    background-color: #ffffff !important;
    background-image: none !important;
    color: #334155 !important;
    border: 1px solid #cbd5e1 !important;
    font-weight: 600 !important;
    transition: background-color 0.15s ease, border-color 0.15s ease, box-shadow 0.15s ease, color 0.15s ease !important;
}
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="baseButton-secondary"] p,
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="baseButton-secondary"] span,
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="stBaseButton-secondary"] p,
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="stBaseButton-secondary"] span {
    color: #334155 !important;
}
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="baseButton-secondary"]:hover,
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="stBaseButton-secondary"]:hover {
    background-color: #fff1f2 !important;
    border-color: #fb7185 !important;
    color: #991b1b !important;
    box-shadow: 0 2px 10px rgba(220, 38, 38, 0.18) !important;
}
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="baseButton-secondary"]:hover p,
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="baseButton-secondary"]:hover span,
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="stBaseButton-secondary"]:hover p,
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="stBaseButton-secondary"]:hover span {
    color: #991b1b !important;
}
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="baseButton-primary"]:disabled,
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="stBaseButton-primary"]:disabled {
    background-color: #fca5a5 !important;
    border-color: #f87171 !important;
    color: #fef2f2 !important;
    box-shadow: none !important;
    opacity: 0.9 !important;
}
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="baseButton-secondary"]:disabled,
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="stBaseButton-secondary"]:disabled {
    background-color: #f1f5f9 !important;
    border-color: #cbd5e1 !important;
    color: #94a3b8 !important;
    box-shadow: none !important;
    opacity: 0.9 !important;
}
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="baseButton-primary"]:disabled,
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="stBaseButton-primary"]:disabled {
    background-color: #fca5a5 !important;
    border-color: #f87171 !important;
    opacity: 0.9 !important;
}
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="baseButton-secondary"]:disabled,
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="stBaseButton-secondary"]:disabled {
    background-color: #f1f5f9 !important;
    border-color: #cbd5e1 !important;
    opacity: 0.9 !important;
}
/* Primary labels: Streamlit nests div/label under <button>; theme can color those. */
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="baseButton-primary"] *,
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="stBaseButton-primary"] *,
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="baseButton-primary"] *,
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="stBaseButton-primary"] *,
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-primary"] *,
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-primary"] *,
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="baseButton-primary"] *,
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-primary"] *,
body:has(.bnd-workflow-active) [data-testid="stMain"] button[kind="primary"] *,
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[kind="primary"] *,
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[kind="primary"] * {
    color: #ffffff !important;
    fill: #ffffff !important;
    stroke: #ffffff !important;
}
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="baseButton-primary"]:hover *,
body:has(.bnd-workflow-active) [data-testid="stMain"] button[data-testid="stBaseButton-primary"]:hover *,
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-primary"]:hover *,
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-primary"]:hover *,
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="baseButton-primary"]:hover *,
body:has(.bnd-workflow-active) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-primary"]:hover *,
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="baseButton-primary"]:hover *,
body:has(.bnd-workflow-active) [data-testid="stSidebar"] button[data-testid="stBaseButton-primary"]:hover * {
    color: #ffffff !important;
    fill: #ffffff !important;
}
/* Hover when pointer is on wrapper / column gutter: match inner button from parent :hover */
body:has(.bnd-workflow-active) [data-testid="stMain"] [data-testid="stButton"]:hover button[data-testid="baseButton-primary"],
body:has(.bnd-workflow-active) [data-testid="stMain"] [data-testid="stButton"]:hover button[data-testid="stBaseButton-primary"],
body:has(.bnd-workflow-active) [data-testid="stMain"] [data-testid="stButton"]:hover button[kind="primary"],
body:has(.bnd-workflow-active) .stApp div.row-widget.stButton[data-testid="stButton"]:hover button[data-testid="baseButton-primary"],
body:has(.bnd-workflow-active) .stApp div.row-widget.stButton[data-testid="stButton"]:hover button[data-testid="stBaseButton-primary"],
body:has(.bnd-workflow-active) div.row-widget.stButton[data-testid="stButton"]:hover button[data-testid="baseButton-primary"],
body:has(.bnd-workflow-active) div.row-widget.stButton[data-testid="stButton"]:hover button[data-testid="stBaseButton-primary"],
body:has(.bnd-workflow-active) [data-testid="stSidebar"] [data-testid="stButton"]:hover button[data-testid="baseButton-primary"],
body:has(.bnd-workflow-active) [data-testid="stSidebar"] [data-testid="stButton"]:hover button[data-testid="stBaseButton-primary"] {
    background-color: #f87171 !important;
    border-color: #fca5a5 !important;
    box-shadow: 0 0 0 3px rgba(248, 113, 113, 0.45) !important;
    background-image: none !important;
    color: #ffffff !important;
    filter: none !important;
}
body:has(.bnd-workflow-active) [data-testid="stMain"] [data-testid="stButton"]:hover button[data-testid="baseButton-primary"] *,
body:has(.bnd-workflow-active) [data-testid="stMain"] [data-testid="stButton"]:hover button[data-testid="stBaseButton-primary"] *,
body:has(.bnd-workflow-active) [data-testid="stMain"] [data-testid="stButton"]:hover button[kind="primary"] *,
body:has(.bnd-workflow-active) .stApp div.row-widget.stButton[data-testid="stButton"]:hover button[data-testid="baseButton-primary"] *,
body:has(.bnd-workflow-active) .stApp div.row-widget.stButton[data-testid="stButton"]:hover button[data-testid="stBaseButton-primary"] *,
body:has(.bnd-workflow-active) div.row-widget.stButton[data-testid="stButton"]:hover button[data-testid="baseButton-primary"] *,
body:has(.bnd-workflow-active) div.row-widget.stButton[data-testid="stButton"]:hover button[data-testid="stBaseButton-primary"] *,
body:has(.bnd-workflow-active) [data-testid="stSidebar"] [data-testid="stButton"]:hover button[data-testid="baseButton-primary"] *,
body:has(.bnd-workflow-active) [data-testid="stSidebar"] [data-testid="stButton"]:hover button[data-testid="stBaseButton-primary"] * {
    color: #ffffff !important;
    fill: #ffffff !important;
}
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"]:hover [role="button"] {
    background-color: #f87171 !important;
    border-color: #fca5a5 !important;
    box-shadow: 0 0 0 3px rgba(248, 113, 113, 0.45) !important;
    background-image: none !important;
    filter: none !important;
}
body:has(.bnd-workflow-active) [data-testid="stMain"] div[data-testid="stButton"]:hover [role="button"] * {
    color: #ffffff !important;
    fill: #ffffff !important;
}
/* Load/Save full-width actions (marker div .bnd-ls-action); default is primary = red */
[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action) button[data-testid="baseButton-primary"],
[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action) button[data-testid="stBaseButton-primary"] {
    background-color: #dc2626 !important;
    background-image: none !important;
    color: #ffffff !important;
    border: 1px solid #b91c1c !important;
    transition: background-color 0.15s ease, border-color 0.15s ease, box-shadow 0.15s ease !important;
}
[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action) button[data-testid="baseButton-primary"]:hover,
[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action) button[data-testid="stBaseButton-primary"]:hover,
[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action):hover button[data-testid="baseButton-primary"],
[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action):hover button[data-testid="stBaseButton-primary"] {
    background-color: #f87171 !important;
    border-color: #fca5a5 !important;
    color: #ffffff !important;
    box-shadow: 0 0 0 3px rgba(248, 113, 113, 0.45) !important;
    background-image: none !important;
    filter: none !important;
}
[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action) button[data-testid="baseButton-primary"] *,
[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action) button[data-testid="stBaseButton-primary"] *,
[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action):hover button[data-testid="baseButton-primary"] *,
[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action):hover button[data-testid="stBaseButton-primary"] * {
    color: #ffffff !important;
    fill: #ffffff !important;
}
[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action) button[data-testid="baseButton-secondary"],
[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action) button[data-testid="stBaseButton-secondary"] {
    background-color: #ffffff !important;
    background-image: none !important;
    color: #334155 !important;
    border: 1px solid #cbd5e1 !important;
    transition: background-color 0.15s ease, border-color 0.15s ease, box-shadow 0.15s ease, color 0.15s ease !important;
}
[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action) button[data-testid="baseButton-secondary"]:hover,
[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action) button[data-testid="stBaseButton-secondary"]:hover {
    background-color: #fff1f2 !important;
    border-color: #fb7185 !important;
    color: #991b1b !important;
    box-shadow: 0 2px 10px rgba(220, 38, 38, 0.18) !important;
}
[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action) button[data-testid="baseButton-primary"]:focus-visible,
[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action) button[data-testid="stBaseButton-primary"]:focus-visible {
    box-shadow: 0 0 0 2px #ffffff, 0 0 0 4px #dc2626 !important;
}
[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action) button[data-testid="baseButton-secondary"]:focus-visible,
[data-testid="stMain"] div[data-testid="stVerticalBlock"]:has(.bnd-ls-action) button[data-testid="stBaseButton-secondary"]:focus-visible {
    box-shadow: 0 0 0 2px #ffffff, 0 0 0 4px #94a3b8 !important;
}
/* Bordered layout panels (Streamlit container border=True) — white on white page */
[data-testid="stMain"] [data-testid="stVerticalBlockBorderWrapper"] {
    border-width: 2px !important;
    border-color: #94a3b8 !important;
    border-radius: 10px !important;
    padding: 0.65rem 0.85rem 0.85rem 0.85rem !important;
    background: #ffffff !important;
    margin-bottom: 0.75rem !important;
}
/* Alert blocks: clearer frame */
[data-testid="stMain"] div[data-testid="stAlert"] {
    border: 1px solid #cbd5e1 !important;
    border-radius: 8px !important;
}
/* Theme-tinted panel borders (background stays white in workflow via rule below) */
body:has(.bnd-theme-warm) [data-testid="stMain"] [data-testid="stVerticalBlockBorderWrapper"] {
    background: #fffefb !important;
    border-color: #d6cfc4 !important;
}
body:has(.bnd-theme-warm) [data-testid="stMain"] div[data-testid="stAlert"] {
    border-color: #c4b8a8 !important;
}
body:has(.bnd-theme-cool) [data-testid="stMain"] [data-testid="stVerticalBlockBorderWrapper"] {
    background: #f8fafc !important;
    border-color: #94a3b8 !important;
}
body:has(.bnd-theme-cool) [data-testid="stMain"] div[data-testid="stAlert"] {
    border-color: #94a3b8 !important;
}
body:has(.bnd-theme-slateivory) [data-testid="stMain"] [data-testid="stVerticalBlockBorderWrapper"] {
    background: #fffff7 !important;
    border-color: #64748b !important;
}
body:has(.bnd-theme-slateivory) [data-testid="stMain"] div[data-testid="stAlert"] {
    border-color: #64748b !important;
}
/* Workflow: always white bordered panels (wins over theme panel tint) */
body:has(.bnd-workflow-active) [data-testid="stMain"] [data-testid="stVerticalBlockBorderWrapper"] {
    background: #ffffff !important;
    border-color: #cbd5e1 !important;
}
/* Simulation & Comparison — tighter controls beside plots (workbench layout) */
body:has(.bnd-sim-compare-active) [data-testid="stMain"] [data-testid="stVerticalBlockBorderWrapper"] {
    padding: 0.45rem 0.6rem 0.5rem 0.6rem !important;
    margin-bottom: 0.4rem !important;
}
body:has(.bnd-sim-compare-active) [data-testid="stMain"] [data-testid="stSlider"] {
    margin-top: 0 !important;
    margin-bottom: 0.15rem !important;
}
body:has(.bnd-sim-compare-active) [data-testid="stMain"] [data-testid="stSelectbox"] {
    margin-bottom: 0.25rem !important;
}
body:has(.bnd-sim-compare-active) [data-testid="stMain"] [data-testid="stCheckbox"] {
    margin-bottom: 0.35rem !important;
}
/* Simulation & Comparison — oscillatory mode banner (navy accent, distinct from live hardware CTAs) */
body:has(.bnd-workflow-active) [data-testid="stMain"] .bnd-sim-osc-banner {
    border-left: 4px solid #1e3a5f !important;
    background: linear-gradient(90deg, #f1f5f9 0%, #ffffff 100%) !important;
    padding: 0.65rem 1rem !important;
    border-radius: 10px !important;
    margin: 0 0 0.85rem 0 !important;
    color: #0f172a !important;
    font-weight: 600 !important;
    font-size: 0.98rem !important;
    border: 1px solid #cbd5e1 !important;
    border-left-width: 4px !important;
}
</style>
""",
        unsafe_allow_html=True,
    )


def _inject_stepwise_compact_layout_css() -> None:
    """Tighten main padding and step rail when ``.bnd-stepwise-active`` is present (see chrome)."""
    st.markdown(
        """
<style>
body:has(.bnd-stepwise-active) .main .block-container {
    padding-top: 0.35rem !important;
    padding-bottom: 0.65rem !important;
}
body:has(.bnd-stepwise-active) [data-testid="stMain"] [data-testid="stVerticalBlockBorderWrapper"] {
    padding: 0.4rem 0.55rem 0.5rem 0.55rem !important;
    margin-bottom: 0.35rem !important;
}
body:has(.bnd-stepwise-active) [data-testid="stMain"] [data-testid="stProgress"] {
    margin-top: 0 !important;
    margin-bottom: 0.3rem !important;
}
body:has(.bnd-stepwise-active) [data-testid="stMain"] [data-testid="stImage"] {
    margin-top: 0 !important;
    margin-bottom: 0 !important;
}
body:has(.bnd-stepwise-active) [data-testid="stMain"] h1,
body:has(.bnd-stepwise-active) [data-testid="stMain"] h2,
body:has(.bnd-stepwise-active) [data-testid="stMain"] h3 {
    margin-top: 0.4rem !important;
    margin-bottom: 0.25rem !important;
}
body:has(.bnd-stepwise-active) [data-testid="stMain"] [data-testid="stHeader"] {
    margin-top: 0.25rem !important;
}
@media (max-width: 1100px) {
  body:has(.bnd-stepwise-active) [data-testid="stMain"] .bnd-stepwise-tab-marker ~ div [data-testid="column"] {
    min-width: 48% !important;
    flex: 1 1 48% !important;
  }
  body:has(.bnd-stepwise-active) [data-testid="stMain"] .bnd-stepwise-tab-marker ~ div [data-testid="stButton"] button {
    min-height: 2.55rem !important;
    white-space: normal !important;
  }
}
</style>
""",
        unsafe_allow_html=True,
    )


def _load_save_button(
    label: str, *, key: str, help: Optional[str] = None, button_type: str = 'primary', **kwargs
) -> bool:
    """Full-width Load/Save-style action; styled black/white via `_inject_load_save_button_theme`."""
    kwargs.pop('use_container_width', None)
    with st.container():
        st.markdown(_BND_LS_ACTION_MARK, unsafe_allow_html=True)
        return bool(
            st.button(label, key=key, help=help, type=button_type, use_container_width=True, **kwargs)
        )


def _hardware_configuration_mode_toggle() -> str:
    """Two-action setup: load existing config or switch to build new."""
    st.session_state.setdefault('gui_config_setup_mode', 'Load existing')
    mode = str(st.session_state.get('gui_config_setup_mode', 'Load existing'))
    if mode not in ('Load existing', 'Build new'):
        mode = 'Load existing'
        st.session_state['gui_config_setup_mode'] = mode
    c1, c2 = st.columns(2)
    with c1:
        if st.button(
            'Load existing',
            key='gui_hw_cfg_mode_load_existing',
            use_container_width=True,
            type='primary' if mode == 'Load existing' else 'secondary',
        ):
            st.session_state['gui_config_setup_mode'] = 'Load existing'
            st.rerun()
    with c2:
        if st.button(
            'Build new',
            key='gui_hw_cfg_mode_build',
            use_container_width=True,
            type='primary' if mode == 'Build new' else 'secondary',
        ):
            st.session_state['gui_config_setup_mode'] = 'Build new'
            st.rerun()
    return str(st.session_state.get('gui_config_setup_mode', 'Load existing'))


def _shared_experiment_dir() -> str:
    """Protocol & biometrics JSON files use the same folder as **Data folder** (section 2) when it exists."""
    d = str(st.session_state.get('gui_data_folder') or '').strip()
    if d:
        norm = os.path.normpath(d)
        if os.path.isdir(norm):
            return norm
    return os.path.normpath(default_templates_dir(_ROOT))


def _session_snapshots_dir() -> str:
    return os.path.normpath(os.path.join(_ROOT, 'SessionSnapshots'))


_AUTOSAVE_SCHEMA_VERSION = 1
_AUTOSAVE_PREFIXES = ('gui_', 'bio_', 'fld_')
_AUTOSAVE_EXCLUDE_KEYS = {
    'gui_autosave_last_saved_at',
    'gui_autosave_last_sig',
    'gui_autosave_bootstrapped',
    'gui_autosave_last_error',
    'gui_state_origin_map',
    'gui_default_state_baseline',
    'gui_recovered_state_baseline',
    'gui_recovery_banner_message',
    'gui_recovery_banner_level',
    'gui_recovery_summary',
    'gui_session_source',
    'gui_setup_confirmed',
    'gui_measurements_confirmed',
    'gui_protocol_confirmed',
    # Action button keys cannot be assigned via session_state.
    'gui_nav_home_stepwise',
    'gui_nav_home_main',
    'gui_kill_daq',
    # Derived caches — recompute after restore; never stringify ndarrays into JSON.
    'gui_last_preview',
    'gui_last_preview_tt',
    'gui_h5_explore_catalog',
    'gui_h5_explore_schema',
    'gui_h5_explore_notes',
    'gui_h5_explore_loaded_path',
    'gui_h5_explore_upload_path',
}


def _autosave_latest_path() -> str:
    return os.path.join(_session_snapshots_dir(), 'autosave_latest.json')


def _autosave_history_path(stamp: str) -> str:
    return os.path.join(_session_snapshots_dir(), f'autosave_{stamp}.json')


def _collect_persistable_state() -> dict[str, object]:
    keys = [
        k
        for k in st.session_state.keys()
        if k.startswith(_AUTOSAVE_PREFIXES) and k not in _AUTOSAVE_EXCLUDE_KEYS and _is_restore_safe_key(k)
    ]
    out: dict[str, object] = {}
    for k in sorted(keys):
        try:
            out[k] = to_json_persistent(st.session_state.get(k), path=k)
        except JsonPersistTypeError as e:
            raise JsonPersistTypeError(f'Autosave cannot serialize session key {k!r}: {e}') from e
    return out


def _is_restore_safe_key(key: str) -> bool:
    if key in _AUTOSAVE_EXCLUDE_KEYS:
        return False
    if key in {'gui_ui_large_text', 'gui_ui_theme'}:
        return False
    if key in {'gui_setup_confirmed', 'gui_measurements_confirmed', 'gui_protocol_confirmed', 'gui_session_source'}:
        return False
    # Buttons / transient actions should never be restored into session state.
    if key.startswith('gui_'):
        # Ephemeral Streamlit widget keys for data-folder dropdown (must follow ``gui_data_folder``).
        if key.startswith('gui_data_folder_dd'):
            return False
        if '_apply_' in key and not key.startswith('gui_apply_tracking_'):
            return False
        action_suffixes = (
            '_btn',
            '_load',
            '_save',
            '_browse',
            '_open',
            '_run',
            '_next',
            '_back',
            '_home',
            '_bottom',
            '_kill_daq',
        )
        if (
            key.startswith('gui_btn_')
            or key.startswith('gui_save_')
            or key.startswith('gui_go_')
            or key.startswith('gui_biometrics_btn_')
            or key.startswith('bio_btn_')
            or key.startswith('gui_nav_')
            or key.startswith('gui_sw_')
            or key.startswith('gui_hw_cfg_mode_')
            or key.startswith('gui_recovery_')
            or key.startswith('gui_kill_')
            or 'upload' in key
            or key.endswith(action_suffixes)
        ):
            return False
    return True


def _build_state_origin_map(state: dict[str, object]) -> dict[str, str]:
    default_base = dict(st.session_state.get('gui_default_state_baseline') or {})
    recovered_base = dict(st.session_state.get('gui_recovered_state_baseline') or {})
    origin: dict[str, str] = {}
    for k, v in state.items():
        if k in recovered_base and to_json_persistent(recovered_base.get(k), path=k) == v:
            origin[k] = 'recovered'
        elif k in default_base and to_json_persistent(default_base.get(k), path=k) == v and k not in recovered_base:
            origin[k] = 'default'
        else:
            origin[k] = 'user'
    return origin


def _build_snapshot_payload(*, source: str) -> dict[str, object]:
    state = _collect_persistable_state()
    provenance = _build_state_origin_map(state)
    return {
        'schema_version': _AUTOSAVE_SCHEMA_VERSION,
        'source': source,
        'saved_at': datetime.now().isoformat(timespec='seconds'),
        'app_route': str(st.session_state.get('gui_app_route') or ''),
        'state': state,
        'provenance': provenance,
    }


def _state_signature(state: dict[str, object]) -> str:
    raw = json.dumps(state, sort_keys=True, ensure_ascii=True, separators=(',', ':'))
    return hashlib.sha256(raw.encode('utf-8')).hexdigest()


def _atomic_write_json(path: str, payload: dict[str, object]) -> None:
    parent = os.path.dirname(path)
    os.makedirs(parent, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix='autosave_', suffix='.tmp', dir=parent)
    try:
        with os.fdopen(fd, 'w', encoding='utf-8') as f:
            json.dump(payload, f, indent=2)
        os.replace(tmp, path)
    except Exception:
        try:
            os.remove(tmp)
        except OSError:
            pass
        raise


def _read_json_file(path: str) -> tuple[Optional[dict], Optional[str]]:
    if not os.path.isfile(path):
        return None, None
    try:
        with open(path, 'r', encoding='utf-8') as f:
            raw = json.load(f)
    except Exception as e:
        return None, f'{type(e).__name__}: {e}'
    if not isinstance(raw, dict):
        return None, 'JSON root must be an object.'
    return raw, None


def _write_snapshot_payload(*, source: str, update_latest: bool) -> tuple[bool, str]:
    try:
        payload = _build_snapshot_payload(source=source)
    except JsonPersistTypeError as e:
        msg = f'{type(e).__name__}: {e}'
        st.session_state['gui_autosave_last_error'] = msg
        return False, msg
    stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    hist_path = _autosave_history_path(stamp) if source == 'autosave' else os.path.join(
        _session_snapshots_dir(), f'session_snapshot_{stamp}.json'
    )
    try:
        _atomic_write_json(hist_path, payload)
        if update_latest:
            _atomic_write_json(_autosave_latest_path(), payload)
        st.session_state['gui_autosave_last_saved_at'] = payload['saved_at']
        st.session_state['gui_state_origin_map'] = dict(payload.get('provenance') or {})
        st.session_state['gui_autosave_last_sig'] = _state_signature(dict(payload.get('state') or {}))
        st.session_state.pop('gui_autosave_last_error', None)
        return True, hist_path
    except Exception as e:
        msg = f'{type(e).__name__}: {e}'
        st.session_state['gui_autosave_last_error'] = msg
        return False, msg


def _save_progress_snapshot() -> tuple[bool, str]:
    """Save current in-app workflow state to JSON (separate from templates)."""
    return _write_snapshot_payload(source='manual_snapshot', update_latest=False)


def _list_manual_snapshot_files(*, max_entries: int = 50) -> list[str]:
    """Most-recent manual session snapshots from ``SessionSnapshots``."""
    root = _session_snapshots_dir()
    if not os.path.isdir(root):
        return []
    out: list[tuple[float, str]] = []
    try:
        for name in os.listdir(root):
            if not (name.startswith('session_snapshot_') and name.endswith('.json')):
                continue
            p = os.path.normpath(os.path.join(root, name))
            try:
                mt = float(os.path.getmtime(p))
            except OSError:
                mt = 0.0
            out.append((mt, p))
    except OSError:
        return []
    out.sort(key=lambda row: row[0], reverse=True)
    return [p for _, p in out[: max(1, int(max_entries))]]


def _validate_snapshot_payload(payload: dict, *, context: str) -> Optional[str]:
    schema = int(payload.get('schema_version', 0) or 0)
    if schema != _AUTOSAVE_SCHEMA_VERSION:
        return f'Unsupported {context} schema version: {schema}.'
    state = payload.get('state')
    if not isinstance(state, dict):
        return f'{context.capitalize()} payload missing valid `state` object.'
    return None


def _restore_snapshot_payload(payload: dict, *, source_label: str) -> None:
    state = dict(payload.get('state') or {})
    for k, v in state.items():
        if _is_restore_safe_key(k):
            try:
                st.session_state[k] = v
            except Exception as e:
                # Streamlit forbids assigning certain widget-managed keys (buttons/uploaders).
                # Ignore those stale keys from older payloads.
                if 'cannot be set using `st.session_state`' not in str(e):
                    raise
    st.session_state['gui_recovered_state_baseline'] = {
        k: to_json_persistent(v, path=k) for k, v in state.items()
    }
    st.session_state['gui_recovery_summary'] = {
        'saved_at': str(payload.get('saved_at') or ''),
        'app_route': str(payload.get('app_route') or ''),
        'source': str(payload.get('source') or source_label),
    }
    st.session_state['gui_recovery_banner_message'] = (
        f"Recovered from {source_label} ({st.session_state['gui_recovery_summary']['saved_at'] or 'unknown time'})."
    )
    st.session_state['gui_recovery_banner_level'] = 'success'
    st.session_state['gui_session_source'] = 'restored'
    _repair_data_path_fields_from_session()


def _load_manual_snapshot(path: str) -> tuple[bool, str]:
    snap_path = os.path.normpath(str(path or '').strip())
    if not snap_path:
        return False, 'Choose a snapshot file first.'
    payload, err = _read_json_file(snap_path)
    if err:
        return False, err
    if payload is None:
        return False, f'Snapshot file not found: `{snap_path}`'
    verr = _validate_snapshot_payload(payload, context='snapshot')
    if verr:
        return False, verr
    _restore_snapshot_payload(payload, source_label='snapshot')
    cfg_sel = _normalize_config_module_name(str(st.session_state.get('gui_load_cfg_select') or ''))
    if cfg_sel:
        mods = discover_config_modules(_ROOT)
        if cfg_sel not in mods:
            st.session_state['gui_recovery_banner_message'] = (
                f"Recovered snapshot, but config module `{cfg_sel}` is missing. Choose a config in Setup."
            )
            st.session_state['gui_recovery_banner_level'] = 'warning'
    return True, snap_path


def _load_latest_autosave() -> tuple[Optional[dict], Optional[str]]:
    payload, err = _read_json_file(_autosave_latest_path())
    if err:
        return None, err
    if payload is None:
        return None, None
    verr = _validate_snapshot_payload(payload, context='autosave')
    if verr:
        return None, verr
    return payload, None


def _restore_autosave_payload(payload: dict) -> None:
    _restore_snapshot_payload(payload, source_label='autosave')


def _bootstrap_autosave_recovery() -> None:
    if st.session_state.get('gui_autosave_bootstrapped'):
        return
    st.session_state['gui_autosave_bootstrapped'] = True
    payload, err = _load_latest_autosave()
    if err:
        st.session_state['gui_recovery_banner_message'] = f'Autosave ignored: {err}'
        st.session_state['gui_recovery_banner_level'] = 'warning'
        return
    if payload is None:
        st.session_state.setdefault('gui_session_source', 'fresh')
        return
    _restore_autosave_payload(payload)
    cfg_sel = _normalize_config_module_name(str(st.session_state.get('gui_load_cfg_select') or ''))
    if cfg_sel:
        mods = discover_config_modules(_ROOT)
        if cfg_sel not in mods:
            st.session_state['gui_recovery_banner_message'] = (
                f"Recovered session, but config module `{cfg_sel}` is missing. Choose a config in Setup."
            )
            st.session_state['gui_recovery_banner_level'] = 'warning'


def _clear_latest_autosave_marker() -> None:
    p = _autosave_latest_path()
    if os.path.isfile(p):
        try:
            os.remove(p)
        except OSError:
            pass


def _autosave_tick(*, force: bool = False) -> None:
    if _nav_route() == 'landing' and not force:
        return
    # Do not call _repair_data_path_fields_from_session() here: it mutates widget keys
    # (e.g. gui_data_filename) after those widgets may already exist this run.
    state = _collect_persistable_state()
    sig = _state_signature(state)
    if not force and sig == str(st.session_state.get('gui_autosave_last_sig') or ''):
        return
    now = float(time.time())
    last = float(st.session_state.get('gui_autosave_last_write_ts') or 0.0)
    if not force and (now - last) < 1.0:
        return
    ok, _ = _write_snapshot_payload(source='autosave', update_latest=True)
    if ok:
        st.session_state['gui_autosave_last_write_ts'] = now


def _render_recovery_status_ui() -> None:
    lvl = str(st.session_state.get('gui_recovery_banner_level') or '')
    msg = str(st.session_state.get('gui_recovery_banner_message') or '').strip()
    if not msg:
        return
    with st.container(border=True):
        if lvl == 'warning':
            st.warning(msg)
        else:
            st.success(msg)
        # Recovery is preserved here as a status banner; the single start-fresh / save-progress
        # control lives in the setup "Start fresh / save progress" panel (no duplicate button).
        c1, c2 = st.columns(2, gap='small')
        with c1:
            with st.expander('Recovery details', expanded=False):
                info = dict(st.session_state.get('gui_recovery_summary') or {})
                st.caption(f"Saved at: `{info.get('saved_at', '(unknown)')}`")
                st.caption(f"Route: `{info.get('app_route', '(unknown)')}`")
        with c2:
            if st.button('Dismiss', key='gui_recovery_dismiss', use_container_width=True):
                st.session_state['gui_recovery_banner_message'] = ''


def _update_state_origin_summary() -> None:
    state = _collect_persistable_state()
    origin = _build_state_origin_map(state)
    st.session_state['gui_state_origin_map'] = origin
    counts = {'default': 0, 'recovered': 0, 'user': 0}
    for v in origin.values():
        if v in counts:
            counts[v] += 1
    st.session_state['gui_state_origin_counts'] = counts


def _reset_workflow_session_to_home(*, clear_autosave: bool = False, target_route: str = 'landing') -> None:
    """Clear workflow state and route to target module."""
    keep = {'gui_ui_theme', 'gui_ui_large_text', 'gui_data_folder', 'gui_data_filename', 'gui_autosave_bootstrapped'}
    keep_prefixes = ('bio_', 'gui_biometrics_', 'gui_genus_', 'gui_specimen_')
    for k in list(st.session_state.keys()):
        if k not in keep and not k.startswith(keep_prefixes):
            st.session_state.pop(k, None)
    st.session_state['gui_app_route'] = str(target_route or 'landing')
    st.session_state['gui_session_source'] = 'fresh'
    if clear_autosave:
        _clear_latest_autosave_marker()


def _protocol_template_option_label(p: Optional[str]):
    if p is None:
        return '— Choose a template —'
    return f'{template_display_label(p)}  ·  {os.path.basename(p)}'


def _biometrics_template_option_label(p: Optional[str]):
    if p is None:
        return '— Choose a biometrics file —'
    return biometrics_template_display_label(p)
from bender_h5_export import build_universal_qc_figure, export_primary_h5, save_universal_qc_figure  # noqa: E402
from bender_h5_plot_helpers import (  # noqa: E402
    align_xy,
    h5_custom_plot_summary,
    list_h5_plot_variables,
    read_h5_series,
)

MOTION_TYPES = frozenset(
    {'dynamic', 'frequency_sweep', 'frequency_step', 'curvature_step', 'step_change'}
)

# Extra fields for motion-series protocols (not fully enumerated in get_dispatch_schema).
MOTION_GUI_FIELDS = [
    ('duration', 'float', 'Duration (s)'),
    ('all_freqs', 'list_float', 'Frequencies (Hz, comma-separated)'),
    (
        'all_amps',
        'list_float',
        'Amplitudes (comma-separated)',
    ),
    ('all_amps_mode', 'select', 'Amplitude mode'),
    ('cycles_per_step', 'int', 'Cycles per step'),
    ('n_end_cycles', 'int', 'End cycles'),
    ('randomize', 'bool', 'Randomize order'),
    ('random_seed', 'optional_int', 'Random seed (optional)'),
    ('stim_cycles_in_step', 'list_int', 'Stim cycles per step (e.g., 2, 3)'),
    ('is_stim', 'bool', 'Enable stimulation'),
    ('stim_pulse_rate', 'float', 'Stim pulse rate (Hz)'),
    ('S1volts', 'float', 'S1 voltage (V)'),
    ('S2volts', 'float', 'S2 voltage (V)'),
    (
        'all_stimduties',
        'list_float',
        'Stim duty (fraction of cycle, comma-separated)',
    ),
    (
        'all_stimphases',
        'list_float',
        'Stim phase (fraction of cycle, comma-separated)',
    ),
]

STEP_CHANGE_EXTRA = [
    ('step_change_frequencies', 'list_float', 'Step-change frequencies'),
    ('step_change_curves', 'list_float', 'Step-change curves/amplitudes'),
    ('step_change_cycles_per_step', 'list_int', 'Step-change cycles per step'),
]

# Required by :meth:`Bender.make_cycles_frequency_sweep` but not used by other motion modes.
FREQUENCY_SWEEP_ONLY_FIELDS = [
    (
        'amplitude_frequency_exponent',
        'float',
        'Amplitude-frequency exponent (alpha)',
    ),
]

_MOTION_ROW_BY_NAME = {row[0]: row for row in MOTION_GUI_FIELDS}


def _motion_parameter_rows(test_type: str):
    """
    Rows (name, kind, label) shown for each motion ``test_type``, aligned with
    :meth:`Bender.run_experiment` branches (dynamic vs sweep/step vs step_change).
    """
    # Shared blocks (see ``organize_cycles`` + stim wiring in ``run_experiment``).
    freq_amp_mode = [
        _MOTION_ROW_BY_NAME['all_freqs'],
        _MOTION_ROW_BY_NAME['all_amps'],
        _MOTION_ROW_BY_NAME['all_amps_mode'],
    ]
    organize_block = [
        _MOTION_ROW_BY_NAME['cycles_per_step'],
        _MOTION_ROW_BY_NAME['n_end_cycles'],
        _MOTION_ROW_BY_NAME['randomize'],
        _MOTION_ROW_BY_NAME['random_seed'],
        _MOTION_ROW_BY_NAME['stim_cycles_in_step'],
    ]
    stim_block = [
        _MOTION_ROW_BY_NAME['is_stim'],
        _MOTION_ROW_BY_NAME['stim_pulse_rate'],
        _MOTION_ROW_BY_NAME['S1volts'],
        _MOTION_ROW_BY_NAME['S2volts'],
        _MOTION_ROW_BY_NAME['all_stimduties'],
        _MOTION_ROW_BY_NAME['all_stimphases'],
    ]
    if test_type == 'dynamic':
        # Duration comes from ``sum(period_by_cycle)`` after ``organize_cycles``, not ``self.duration``.
        return [*freq_amp_mode, *organize_block, *stim_block]
    if test_type == 'frequency_sweep':
        return [
            _MOTION_ROW_BY_NAME['duration'],
            *freq_amp_mode,
            *FREQUENCY_SWEEP_ONLY_FIELDS,
            *stim_block,
        ]
    if test_type in ('frequency_step', 'curvature_step'):
        return [_MOTION_ROW_BY_NAME['duration'], *freq_amp_mode, *stim_block]
    if test_type == 'step_change':
        return [*list(STEP_CHANGE_EXTRA), *stim_block]
    return list(MOTION_GUI_FIELDS)


ALL_AMPS_MODE_OPTIONS = ('strain', 'strain_pct', 'curvature', 'angle')
VELOCITY_MODE_OPTIONS = ('strain_rate', 'strain_pct_rate', 'curvature_rate', 'angle_vel')

DATA_FOLDER_HELP = (
    'Choose the folder where experiment files live. Enter the folder path only — do not put the file name here. '
    'Runs, exports, protocol templates, and biometrics files can all use this folder.'
)
DATA_FILE_NAME_HELP = (
    'This is the name of your saved measurements file (HDF5). Enter only the file name, not the full path. '
    'The app joins it with the data folder field above to build where data is saved. You may type .h5 or leave it off. '
    'If that exact file already exists, the app uses a new name like my_run_001.h5 so nothing is overwritten.'
)

# Specimen density presets: ρ (g/cm³) × 1e-3 → g/mm³ (same mass per volume, different unit).
BIO_DENSITY_PRESET_LABELS = (
    'Custom — edit the number below',
    'Water-like (~1.00 g/cm³)',
    'Skeletal muscle / soft tissue (~1.06 g/cm³)',
    'Cortical bone (~1.9 g/cm³)',
)
BIO_DENSITY_PRESET_G_PER_MM3 = {
    'Water-like (~1.00 g/cm³)': 1.0e-3,
    'Skeletal muscle / soft tissue (~1.06 g/cm³)': 1.06e-3,
    'Cortical bone (~1.9 g/cm³)': 1.9e-3,
}

BIO_DBEND_FIELD_HELP = (
    'Distance **along the body** (mm) from your length reference (same as TL/SL, often snout or a fixed landmark) to the '
    '**midpoint between the two clamps** — i.e. where the bending test is centered. Use **0** only if that reference is '
    'already at the segment center. Saved as `dbend` / `test_segment_position_mm`.'
)

BIO_PROF_CLAMP_FIELD_HELP = (
    'Used only for **rotating hardware** mass/MOI in the profiled inertial model: offset from the **bend / rotation axis** '
    'to the clamps (mm). The code adds half of its built-in clamp depth to this value when estimating clamp contribution. '
    'Saved as `specimen_profile_clamp_offset_mm`.'
)


def _resolved_bio_prof_rho_g_per_mm3() -> float:
    """Density for profile/inertial apply: preset table value, or the custom ``bio_prof_rho`` number input."""
    label = str(st.session_state.get('bio_prof_rho_preset') or '')
    v = BIO_DENSITY_PRESET_G_PER_MM3.get(label)
    if v is not None:
        return float(v)
    return float(st.session_state['bio_prof_rho'])


def _queue_bio_prof_rho_widget_sync_from_preset() -> None:
    """
    Preset-driven density cannot assign ``st.session_state['bio_prof_rho']`` during the same run as the
    ``number_input(..., key='bio_prof_rho')`` widget (Streamlit blocks writes to widget-bound keys after mount).
    Queue the value and flush at the start of the next run (see ``_flush_pending_bio_prof_rho_sync``).
    """
    label = str(st.session_state.get('bio_prof_rho_preset') or '')
    v = BIO_DENSITY_PRESET_G_PER_MM3.get(label)
    if v is not None:
        st.session_state['_pending_sync_bio_prof_rho'] = float(v)


def _flush_pending_bio_prof_rho_sync() -> None:
    """Copy queued preset density into the widget key before ``bio_prof_rho`` number_input is created."""
    if '_pending_sync_bio_prof_rho' in st.session_state:
        st.session_state['bio_prof_rho'] = float(st.session_state.pop('_pending_sync_bio_prof_rho'))


ISOMETRIC_STIM_JSON_HELP = (
    'Optional per-step timing and stimulation. **Leave `{}`** unless your protocol needs custom ramps, holds, or stim.\n\n'
    '**Advanced (JSON text):** use `{ "key": value }` with commas. Common keys: **ramp_duration_s**, **hold_duration_s**, '
    '**stim_onset_s** (signed, relative to hold start), **stim_duration_s**, **inter_step_interval_s** (0 = back-to-back), '
    '**is_stim**, **stim_pulse_rate** (Hz), **device_name** (null = NI config). '
    'Example: `{"ramp_duration_s": 2, "hold_duration_s": 5, "stim_onset_s": 0.5, "stim_duration_s": 4, "stim_pulse_rate": 75}`'
)
ISOVELOCITY_STIM_JSON_HELP = (
    'Optional segment timing and stimulation. **Leave `{}`** unless you need overrides.\n\n'
    '**Advanced (JSON text):** **stim_onset_s** (signed, relative to constant-velocity start), **stim_duration_s**, '
    '**is_stim**, **stim_pulse_rate**, **device_name**; '
    'optionally **iso_duration_s** / **pre_hold_s** to override the main fields.'
)

RANDOM_SEED_HELP = (
    'Only used when **randomize order** is checked. Enter an integer for a **reproducible** shuffle (same order '
    'every run); leave empty for a **different** random order each time.'
)

RECRUITMENT_FIELD_HELP = (
    'How stimulation and/or motor commands are routed across left vs right: simultaneous, sequential halves of a '
    'step, or unilateral (one side). If **Perform test on both sides (bilateral)** is checked, simultaneous or '
    'unilateral recruitment is **upgraded to bilateral sequential** at run time so the motor can bend toward left '
    'then right with stim matched to each hold.'
)

ISOVELOCITY_WIDGET_LABEL = {
    'isovelocity_min_vel': 'Minimum velocity',
    'isovelocity_max_vel': 'Maximum velocity',
    'isovelocity_starting_strain': 'Starting posture',
    'isovelocity_num_steps': 'Number of velocity steps',
    'isovelocity_starting_strain_mode': 'Unit for starting posture',
    'isovelocity_velocity_mode': 'Unit for min/max velocity',
    'isovelocity_randomize': 'Randomize order of velocity steps',
    'isovelocity_random_seed': 'Random seed (optional)',
    'isovelocity_iso_duration_s': 'Constant-velocity bend duration (s)',
    'isovelocity_pre_hold_s': 'Pre-hold at starting angle (s)',
}

ISOVELOCITY_FIELD_HELP = {
    'isovelocity_min_vel': 'Lower end of the velocity sweep in the units selected below (e.g. deg/s, dε/dt, or dκ/dt).',
    'isovelocity_max_vel': 'Upper end of the velocity sweep (same units as minimum).',
    'isovelocity_velocity_mode': (
        'How **minimum** / **maximum velocity** are interpreted before conversion to motor deg/s '
        '(strain_rate, curvature_rate, or angle_vel).'
    ),
    'isovelocity_starting_strain': 'Initial posture before each velocity step; interpreted with **starting strain mode** (e.g. decimal ε vs motor angle).',
    'isovelocity_num_steps': 'How many commanded angular velocities between min and max (inclusive). Each is a separate trial from the same start posture.',
    'isovelocity_starting_strain_mode': 'Whether **starting posture** is decimal ε, percent strain, curvature κ, or motor angle (deg).',
    'isovelocity_randomize': 'Shuffle the order of velocity steps. Use **random seed** if you want the same order every time.',
    'isovelocity_random_seed': RANDOM_SEED_HELP,
    'isovelocity_iso_duration_s': 'After the quiet **pre-hold**, how long (seconds) the motor holds **constant angular velocity** before the step ends. This is not the pre-hold time—that is **pre-hold**.',
    'isovelocity_pre_hold_s': 'Time (s) at the starting posture before each trial’s constant-velocity segment begins.',
    'recruitment': RECRUITMENT_FIELD_HELP,
    'lateral_mode': (
        'Expert only. Leave **blank** unless you need a custom stim-router label. Normal left/right/bilateral behavior '
        'is set with **recruitment** above; this overrides that name inside merged stim settings.'
    ),
    'bilateral_mirror_motor': (
        'Two constant-velocity bends in one trial: first toward **left**, then toward **right**, with stimulation on '
        'each side in turn. If recruitment is **bilateral simultaneous** or **left/right**, it is upgraded to '
        '**bilateral sequential** for this run.'
    ),
    'bilateral_sequential_left_frac': 'Share of each step spent on the **left** side before the right (0–1).',
}

ISOMETRIC_FIELD_HELP = {
    'isometric_initial': 'First target in the sweep; units follow **isometric mode** (strain, curvature, or angle).',
    'isometric_final': 'Last target in the sweep (same units as initial).',
    'isometric_num_steps': 'Number of steps between initial and final (endpoints included).',
    'isometric_mode': 'How **initial** / **final** are interpreted: decimal ε, percent strain, κ, or motor angle.',
    'isometric_randomize': 'Shuffle step order. Use **random seed** for a fixed order across runs.',
    'isometric_random_seed': RANDOM_SEED_HELP,
    'isometric_inter_step_interval_s': 'Idle time after each step’s acquisition before the next ramp (0 = back-to-back).',
    'recruitment': RECRUITMENT_FIELD_HELP,
    'lateral_mode': (
        'Expert only. Leave **blank** unless you need a custom stim-router label. **Recruitment** sets routine '
        'left/right/bilateral behavior; this overrides the router name in merged stim params.'
    ),
    'bilateral_mirror_motor': (
        'Left posture then right, with stim on each side in turn. Recruitment is coerced to **bilateral sequential** '
        'when needed so both bends run.'
    ),
    'isometric_mirror_target_left': (
        'Optional. When **both** left and right values are set, each bilateral step uses this magnitude for the **first** '
        'hold (bend toward **left** specimen side) instead of the sweep step value. Same units as **isometric mode** '
        '(e.g. decimal ε 0.1 = 10% strain, or use **strain_pct** with 10).'
    ),
    'isometric_mirror_target_right': (
        'Optional. Magnitude for the **second** hold (bend toward **right** specimen side). Must be set together with '
        'the left target; same units as **isometric mode**.'
    ),
    'bilateral_sequential_left_frac': 'Share of each step on the **left** side before the right (0–1).',
}

MOTION_FIELD_HELP = {
    'random_seed': RANDOM_SEED_HELP,
    'duration': 'Total motion timeline length (s) for protocols that use a single continuous run.',
    'all_freqs': 'Comma-separated list of frequencies (Hz) for each step or segment.',
    'all_amps': 'Comma-separated amplitudes; meaning depends on **amplitude interpretation** (strain vs angle, etc.).',
    'all_amps_mode': 'How each entry in **amplitudes** is converted to curvature / motor angle.',
    'cycles_per_step': 'Oscillation cycles at each frequency step before advancing.',
    'n_end_cycles': 'Extra cycles at the last frequency at the end of the run.',
    'randomize': 'Shuffle the order of frequency/amplitude steps when applicable.',
    'stim_cycles_in_step': 'Which cycle indices (1-based) within each step receive stimulation.',
    'is_stim': 'Enable patterned stimulation during the motion.',
    'stim_pulse_rate': 'Carrier pulse rate for stimulation (Hz).',
    'S1volts': 'Stimulus channel 1 amplitude (V).',
    'S2volts': 'Stimulus channel 2 amplitude (V).',
    'all_stimduties': (
        'Comma-separated **duty** values (each 0–1 as a **fraction of one motion cycle**). '
        'Used with frequencies and amplitudes to build stimulation timing via ``organize_cycles``; '
        'leave blank to use defaults (e.g. 0.3) in preview.'
    ),
    'all_stimphases': (
        'Comma-separated **phase** values (each 0–1 as a **fraction of one motion cycle**, same convention as duty). '
        'Leave blank to use defaults (e.g. 0.5) in preview.'
    ),
    'amplitude_frequency_exponent': 'For frequency sweep: exponent α so amplitude ∝ f^α relative to the sweep start.',
    'step_change_frequencies': 'Comma-separated frequencies for each step-change segment.',
    'step_change_curves': 'Comma-separated curvature / amplitude targets per segment.',
    'step_change_cycles_per_step': 'Comma-separated integers: cycles per step-change segment.',
}


def _format_strain_or_amp_mode(opt: str) -> str:
    """User-facing labels: make decimal vs percent strain explicit."""
    o = str(opt)
    if o == 'strain':
        return 'strain — decimal ε (0.05 = 5%)'
    if o == 'strain_pct':
        return 'strain_pct — percent scale (5 = 5%)'
    if o == 'curvature':
        return 'curvature — κ (1/m)'
    if o == 'angle':
        return 'angle — motor (deg)'
    return o


def _format_velocity_mode(opt: str) -> str:
    o = str(opt)
    if o == 'strain_rate':
        return 'strain_rate — dε/dt (1/s, decimal ε)'
    if o == 'strain_pct_rate':
        return 'strain_pct_rate — d(% strain)/dt (%/s)'
    if o == 'curvature_rate':
        return 'curvature_rate — dκ/dt (1/m/s)'
    if o == 'angle_vel':
        return 'angle_vel — motor angular rate (deg/s)'
    return o


RECRUITMENT_OPTIONS = (
    'bilateral_simultaneous',
    'bilateral_sequential',
    'left',
    'right',
)

BILATERAL_MIRROR_LABEL = 'Perform test on both sides (bilateral)'
LATERAL_MODE_LABEL = 'Stim routing override (optional; experts only)'

BLOCK_DIRECTION_OPTIONS = ('left', 'right')
BLOCK_STIM_SIDES_OPTIONS = ('left', 'right', 'both')
BLOCK_DIRECTION_LABELS = {'left': 'Bend LEFT', 'right': 'Bend RIGHT'}
BLOCK_STIM_SIDES_LABELS = {
    'left': 'Stim LEFT',
    'right': 'Stim RIGHT',
    'both': 'Stim BOTH',
}

_BLOCK_SEQUENCE_PROCEDURE_KEYS = frozenset({
    'block_sequence',
    'left_stim_voltage',
    'right_stim_voltage',
    'block_reset_ramp_duration_s',
})


def _seed_block_sequence_widget_state(b: Bender) -> None:
    """Initialize block-sequence widget keys before first render."""
    seq = getattr(b, 'block_sequence', None)
    if not isinstance(seq, list) or not seq:
        seq = [{'direction': 'left', 'stim_sides': 'left'}]
    if 'gui_block_seq_count' not in st.session_state:
        n = len(seq)
        st.session_state['gui_block_seq_count'] = max(1, min(12, int(n)))
    count = int(st.session_state.get('gui_block_seq_count', 1))
    default_blocks = seq
    for i in range(max(count, len(default_blocks))):
        block = default_blocks[i] if i < len(default_blocks) else {'direction': 'left', 'stim_sides': 'left'}
        d_sk = _widget_key(f'block_{i}_direction')
        s_sk = _widget_key(f'block_{i}_stim_sides')
        if d_sk not in st.session_state:
            d = str(block.get('direction', 'left')).lower()
            st.session_state[d_sk] = d if d in BLOCK_DIRECTION_OPTIONS else 'left'
        if s_sk not in st.session_state:
            s = str(block.get('stim_sides', 'left')).lower()
            st.session_state[s_sk] = s if s in BLOCK_STIM_SIDES_OPTIONS else 'left'
    if _widget_key('left_stim_voltage') not in st.session_state:
        st.session_state[_widget_key('left_stim_voltage')] = float(
            getattr(b, 'left_stim_voltage', 5.0) or 5.0
        )
    if _widget_key('right_stim_voltage') not in st.session_state:
        st.session_state[_widget_key('right_stim_voltage')] = float(
            getattr(b, 'right_stim_voltage', 5.0) or 5.0
        )
    if _widget_key('block_reset_ramp_duration_s') not in st.session_state:
        st.session_state[_widget_key('block_reset_ramp_duration_s')] = float(
            getattr(b, 'block_reset_ramp_duration_s', 2.0) or 2.0
        )


def _render_block_sequence_fields(b: Bender) -> Optional[dict]:
    """
    Render block-sequence UI; return updates dict for Apply, or None if validation fails.
    Block sequence is always active (default: one block, bend LEFT / stim LEFT).
    """
    _seed_block_sequence_widget_state(b)

    count = int(
        st.number_input(
            'Number of blocks',
            min_value=1,
            max_value=12,
            step=1,
            key='gui_block_seq_count',
            help='Each block resets to neutral, then runs the full step protocol.',
        )
    )
    blocks = []
    for i in range(count):
        c1, c2 = st.columns(2)
        with c1:
            direction = st.selectbox(
                f'Block {i + 1} — bend direction',
                list(BLOCK_DIRECTION_OPTIONS),
                key=_widget_key(f'block_{i}_direction'),
                format_func=lambda x, _l=BLOCK_DIRECTION_LABELS: _l.get(x, x),
            )
        with c2:
            stim_sides = st.selectbox(
                f'Block {i + 1} — stim sides',
                list(BLOCK_STIM_SIDES_OPTIONS),
                key=_widget_key(f'block_{i}_stim_sides'),
                format_func=lambda x, _l=BLOCK_STIM_SIDES_LABELS: _l.get(x, x),
            )
        blocks.append({'direction': direction, 'stim_sides': stim_sides})

    left_v = float(
        st.number_input(
            'Left stim voltage (V)',
            key=_widget_key('left_stim_voltage'),
            format='%.6g',
            min_value=0.0,
            help='Voltage on the LEFT stim channel when a block includes LEFT or BOTH.',
        )
    )
    right_v = float(
        st.number_input(
            'Right stim voltage (V)',
            key=_widget_key('right_stim_voltage'),
            format='%.6g',
            min_value=0.0,
            help='Voltage on the RIGHT stim channel when a block includes RIGHT or BOTH.',
        )
    )
    reset_ramp = float(
        st.number_input(
            'Neutral reset ramp duration (s)',
            key=_widget_key('block_reset_ramp_duration_s'),
            format='%.6g',
            min_value=0.0,
            help='Motor ramp time from the previous angle back to 0° before each block.',
        )
    )

    sides_used = {b['stim_sides'] for b in blocks}
    if ('left' in sides_used or 'both' in sides_used) and not (np.isfinite(left_v) and left_v > 0):
        st.error('Left stim voltage must be finite and > 0 when any block uses LEFT or BOTH stim.')
        return None
    if ('right' in sides_used or 'both' in sides_used) and not (np.isfinite(right_v) and right_v > 0):
        st.error('Right stim voltage must be finite and > 0 when any block uses RIGHT or BOTH stim.')
        return None
    if not (np.isfinite(reset_ramp) and reset_ramp >= 0):
        st.error('Neutral reset ramp duration must be finite and ≥ 0 s.')
        return None

    try:
        b._validate_block_sequence_voltages(blocks, left_v, right_v)
        b._normalize_block_sequence(blocks)
    except ValueError as exc:
        st.error(str(exc))
        return None

    return {
        'block_sequence': blocks,
        'left_stim_voltage': left_v,
        'right_stim_voltage': right_v,
        'block_reset_ramp_duration_s': reset_ramp,
    }


def _parse_float_list(s: str):
    s = (s or '').strip()
    if not s:
        return None
    out = []
    for part in s.split(','):
        part = part.strip()
        if part:
            out.append(float(part))
    return out if out else None


def _parse_int_list(s: str):
    s = (s or '').strip()
    if not s:
        return None
    out = []
    for part in s.split(','):
        part = part.strip()
        if part:
            out.append(int(round(float(part))))
    return out if out else None


def _widget_key(name: str) -> str:
    return f'fld_{name}'


def _stim_params_dict_from_bender(b: Bender, attr_name: str) -> dict:
    sp = getattr(b, attr_name, None)
    return dict(sp) if isinstance(sp, dict) else {}


def _stim_onset_duration_seed_from_sp(sp: dict, *, default_onset: float, default_duration: float) -> tuple[float, float]:
    """Seed widget values from stim params, migrating legacy settle/pre-iso keys."""
    if sp.get('stim_onset_s') is not None:
        onset = float(sp['stim_onset_s'])
    elif float(sp.get('pre_iso_stim_duration_s', 0) or 0) > 0:
        onset = -float(sp['pre_iso_stim_duration_s'])
    else:
        onset = float(sp.get('settle_before_stim_s', default_onset) or default_onset)
    if sp.get('stim_duration_s') is not None:
        duration = float(sp['stim_duration_s'])
    else:
        duration = float(default_duration)
    return onset, duration


def _seed_isovelocity_stim_widget_state(b: Bender) -> None:
    """Initialize isovelocity stim widget keys before first render (Apply-button pattern)."""
    sp = _stim_params_dict_from_bender(b, 'isovelocity_stim_params')
    spr_raw = sp.get('stim_pulse_rate', None)
    spr_default = float(
        spr_raw if spr_raw is not None else getattr(b, 'stim_pulse_rate', 75.0) or 75.0
    )
    iso_seg = float(getattr(b, 'isovelocity_iso_duration_s', 0.2) or 0.2)
    onset, duration = _stim_onset_duration_seed_from_sp(sp, default_onset=0.0, default_duration=iso_seg)
    defaults = {
        'isovelocity_stim_enable': bool(sp.get('is_stim', False)),
        'isovelocity_stim_pulse_rate': spr_default,
        'isovelocity_stim_onset_s': onset,
        'isovelocity_stim_duration_s': duration,
    }
    for name, val in defaults.items():
        sk = _widget_key(name)
        if sk not in st.session_state:
            st.session_state[sk] = val


def _render_isovelocity_stim_fields(b: Bender) -> Optional[dict]:
    """
    Render isovelocity stimulation widgets and return assembled ``isovelocity_stim_params``.

    Returns None when validation fails (caller must skip setattr on Apply).
    """
    _seed_isovelocity_stim_widget_state(b)
    enable_sk = _widget_key('isovelocity_stim_enable')
    enable = bool(
        st.checkbox(
            'Enable stimulation',
            key=enable_sk,
            help=(
                'Gates only whether stimulation fires during the run. The stim settings below stay '
                'visible and editable either way, so values are preserved while toggling and you can '
                'pre-condition (run with stim off, then enable) or tune.'
            ),
        )
    )
    onset_sk = _widget_key('isovelocity_stim_onset_s')
    duration_sk = _widget_key('isovelocity_stim_duration_s')
    pr_sk = _widget_key('isovelocity_stim_pulse_rate')

    pulse_rate = float(
        st.number_input(
            'Stim pulse rate (Hz)',
            key=pr_sk,
            format='%.6g',
            min_value=0.0,
            help=MOTION_FIELD_HELP.get('stim_pulse_rate'),
        )
    )
    onset = float(
        st.number_input(
            'Stim onset (s, relative to active-segment start)',
            key=onset_sk,
            format='%.6g',
            help=(
                'Signed seconds from constant-velocity segment start. Negative = pre-activation during '
                'pre-hold; 0 = at segment start; positive = after.'
            ),
        )
    )
    duration = float(
        st.number_input(
            'Stim duration (s)',
            key=duration_sk,
            format='%.6g',
            min_value=0.0,
            help='How long stimulation runs from onset.',
        )
    )

    # Validation guards apply only when stim will fire; disabled stim never blocks Apply so
    # values are preserved for pre-conditioning and tuning.
    if enable:
        if not (np.isfinite(pulse_rate) and pulse_rate > 0):
            _st_error_actions(
                'Stim pulse rate invalid.',
                ['Enter a value > 0 Hz', 'Or uncheck Enable stimulation'],
            )
            return None
        if not (np.isfinite(duration) and duration > 0):
            _st_error_actions('Stim duration invalid.', ['Enter a value > 0 s'])
            return None
        if not np.isfinite(onset):
            _st_error_actions('Stim onset invalid.', ['Enter a finite value in seconds'])
            return None

    params: dict = {
        'is_stim': enable,
        'stim_onset_s': onset if np.isfinite(onset) else 0.0,
        'stim_duration_s': duration,
        'stim_pulse_rate': pulse_rate,
    }
    return params


def _seed_isometric_stim_widget_state(b: Bender) -> None:
    """Initialize isometric stim widget keys before first render."""
    sp = _stim_params_dict_from_bender(b, 'isometric_stim_params')
    spr_raw = sp.get('stim_pulse_rate', None)
    spr_default = float(
        spr_raw if spr_raw is not None else getattr(b, 'stim_pulse_rate', 75.0) or 75.0
    )
    onset, duration = _stim_onset_duration_seed_from_sp(sp, default_onset=0.5, default_duration=4.5)
    defaults = {
        'isometric_stim_enable': bool(sp.get('is_stim', False)),
        'isometric_stim_pulse_rate': spr_default,
        'isometric_stim_onset_s': onset,
        'isometric_stim_duration_s': duration,
        'isometric_hold_duration_s': float(sp.get('hold_duration_s', 5.0)),
    }
    for name, val in defaults.items():
        sk = _widget_key(name)
        if sk not in st.session_state:
            st.session_state[sk] = val


def _render_isometric_stim_fields(b: Bender) -> Optional[dict]:
    """Render isometric stimulation widgets; return assembled ``isometric_stim_params`` or None."""
    _seed_isometric_stim_widget_state(b)
    enable_sk = _widget_key('isometric_stim_enable')
    enable = bool(
        st.checkbox(
            'Enable stimulation',
            key=enable_sk,
            help=(
                'Gates only whether stimulation fires during the run. The stim settings below stay '
                'visible and editable either way, so values are preserved while toggling and you can '
                'pre-condition (run with stim off, then enable) or tune.'
            ),
        )
    )
    onset_sk = _widget_key('isometric_stim_onset_s')
    duration_sk = _widget_key('isometric_stim_duration_s')
    pr_sk = _widget_key('isometric_stim_pulse_rate')

    hold_sk = _widget_key('isometric_hold_duration_s')
    hold_duration = float(
        st.number_input(
            'Hold duration (s)',
            key=hold_sk,
            format='%.6g',
            min_value=0.0,
            help=(
                'How long (seconds) the motor holds at each target angle (the active segment). '
                'Stim onset and duration are measured relative to the start of this hold.'
            ),
        )
    )

    pulse_rate = float(
        st.number_input(
            'Stim pulse rate (Hz)',
            key=pr_sk,
            format='%.6g',
            min_value=0.0,
            help=MOTION_FIELD_HELP.get('stim_pulse_rate'),
        )
    )
    onset = float(
        st.number_input(
            'Stim onset (s, relative to active-segment start)',
            key=onset_sk,
            format='%.6g',
            help=(
                'Signed seconds from hold start (active segment). Negative = stim begins during the '
                'pre-hold ramp; it cannot start before the ramp begins (limited by ramp duration). '
                '0 = at hold start; positive = later in the hold.'
            ),
        )
    )
    duration = float(
        st.number_input(
            'Stim duration (s)',
            key=duration_sk,
            format='%.6g',
            min_value=0.0,
            help='How long stimulation runs from onset.',
        )
    )

    # Validation guards apply only when stim will fire; disabled stim never blocks Apply so
    # values are preserved for pre-conditioning and tuning.
    if enable:
        if not (np.isfinite(pulse_rate) and pulse_rate > 0):
            _st_error_actions(
                'Stim pulse rate invalid.',
                ['Enter a value > 0 Hz', 'Or uncheck Enable stimulation'],
            )
            return None
        if not (np.isfinite(duration) and duration > 0):
            _st_error_actions('Stim duration invalid.', ['Enter a value > 0 s'])
            return None
        if not np.isfinite(onset):
            _st_error_actions('Stim onset invalid.', ['Enter a finite value in seconds'])
            return None
    if not (np.isfinite(hold_duration) and hold_duration > 0):
        _st_error_actions('Hold duration invalid.', ['Enter a value > 0 s'])
        return None

    params: dict = {
        'is_stim': enable,
        'stim_onset_s': onset if np.isfinite(onset) else 0.0,
        'stim_duration_s': duration,
        'hold_duration_s': hold_duration,
        'stim_pulse_rate': pulse_rate,
    }
    return params


def _validate_procedure_stim_timing(b: Bender, updates: dict, tt: str) -> Optional[str]:
    """Return an error message when stim timing would bleed outside segment bounds."""
    try:
        if tt == 'isometric':
            sp = updates.get('isometric_stim_params')
            if not isinstance(sp, dict) or not sp.get('is_stim'):
                return None
            num_steps = updates.get('isometric_num_steps', getattr(b, 'isometric_num_steps', 1))
            hold_s = float(sp.get('hold_duration_s', 5.0))
            ramp_s = float(sp.get('ramp_duration_s', 2.0))
            b._validate_stim_timing_for_steps(
                sp,
                test_type='isometric',
                num_steps=int(num_steps),
                pre_hold_at_start_s=ramp_s,
                segment_duration_s=hold_s,
            )
        elif tt == 'isovelocity':
            sp = updates.get('isovelocity_stim_params')
            if not isinstance(sp, dict) or not sp.get('is_stim'):
                return None
            num_steps = updates.get('isovelocity_num_steps', getattr(b, 'isovelocity_num_steps', 1))
            pre_hold_s = float(
                updates.get('isovelocity_pre_hold_s', getattr(b, 'isovelocity_pre_hold_s', 0.3)) or 0.3
            )
            iso_s = float(
                updates.get('isovelocity_iso_duration_s', getattr(b, 'isovelocity_iso_duration_s', 0.2)) or 0.2
            )
            b._validate_stim_timing_for_steps(
                sp,
                test_type='isovelocity',
                num_steps=int(num_steps),
                pre_hold_at_start_s=pre_hold_s,
                segment_duration_s=iso_s,
            )
    except ValueError as exc:
        return str(exc)
    return None


def _get_session_value(b: Bender, name: str, default=None):
    if hasattr(b, name):
        v = getattr(b, name)
        if v is not None:
            return v
    return default


def _clip_action_words(text: str, max_words: int = 5) -> str:
    words = str(text).strip().split()
    if len(words) <= max_words:
        return str(text).strip()
    return ' '.join(words[:max_words])


def _st_error_actions(headline: str, actions: list[str], *, max_words: int = 5) -> None:
    """Show ``st.error`` with a headline and bulleted next-step phrases (≤ ``max_words`` each)."""
    parts = [f'**{headline.strip()}**']
    for a in actions:
        s = str(a).strip()
        if s:
            parts.append(f'- {_clip_action_words(s, max_words)}')
    st.error('\n\n'.join(parts))


def _st_error_detail(headline: str, actions: list[str], detail: str) -> None:
    """Like `_st_error_actions` plus an expander with full text (e.g. server or Python message)."""
    _st_error_actions(headline, actions)
    with st.expander('Details'):
        st.code(str(detail).strip() or '(empty)', language=None)


def _missing_fields_to_actions(missing: list[str]) -> list[str]:
    """Turn validator missing names into short imperative lines (capped for UI size)."""
    names = [str(n).strip() for n in missing if str(n).strip()]
    lines = [_clip_action_words(f'Set {n}', 5) for n in names[:12]]
    if len(names) > 12:
        lines.append('Review other required fields')
    return lines


def _render_field(b: Bender, name: str, kind: str, label: str, *, help_text: Optional[str] = None):
    """
    Render one widget; return value to setattr, or None to skip.

    Streamlit keeps keyed widget values in ``st.session_state``. Passing ``value=`` / ``index=``
    every rerun is ignored once the key exists, which can leave the UI stale (e.g. empty text)
    while ``Bender`` still holds good data. We only seed session state when the key is missing.
    """
    sk = _widget_key(name)
    cur = _get_session_value(b, name)
    h = help_text if help_text else None

    if kind == 'float':
        if sk not in st.session_state:
            st.session_state[sk] = float(cur) if cur is not None else 0.0
        return float(st.number_input(label, key=sk, format='%.6g', help=h))

    if kind == 'int':
        if sk not in st.session_state:
            seed = int(cur) if cur is not None else 0
            if name == 'n_end_cycles' and seed < 0:
                seed = 0
            if name == 'cycles_per_step' and seed < 1:
                seed = 1
            st.session_state[sk] = seed
        elif name == 'n_end_cycles' and int(st.session_state.get(sk, 0) or 0) < 0:
            st.session_state[sk] = 0
        elif name == 'cycles_per_step' and int(st.session_state.get(sk, 1) or 1) < 1:
            st.session_state[sk] = 1
        min_v = 0 if name == 'n_end_cycles' else (1 if name == 'cycles_per_step' else None)
        return int(st.number_input(label, key=sk, step=1, min_value=min_v, help=h))

    if kind == 'optional_int':
        use_key = f'{sk}_use'
        if use_key not in st.session_state:
            st.session_state[use_key] = cur is not None
        use = st.checkbox(f'{label} (use custom)', key=use_key, help=h)
        if not use:
            return None
        if sk not in st.session_state:
            st.session_state[sk] = int(cur) if cur is not None else 0
        return int(st.number_input(label, key=sk, step=1))

    if kind == 'optional_float':
        use_key = f'{sk}_use'
        if use_key not in st.session_state:
            st.session_state[use_key] = cur is not None
        use = st.checkbox(f'{label} (enable)', key=use_key, help=h)
        if not use:
            return None
        if sk not in st.session_state:
            st.session_state[sk] = float(cur) if cur is not None else 0.0
        return float(st.number_input(label, key=sk, format='%.6g'))

    if kind == 'bool':
        if sk not in st.session_state:
            st.session_state[sk] = bool(cur) if cur is not None else False
        return bool(st.checkbox(label, key=sk, help=h))

    if kind == 'str':
        if sk not in st.session_state:
            st.session_state[sk] = str(cur or '')
        return str(st.text_input(label, key=sk, help=h))

    if kind == 'select':
        opts = list(ALL_AMPS_MODE_OPTIONS)
        dv = str(cur or 'strain')
        if dv not in opts:
            dv = 'strain'
        if sk not in st.session_state:
            st.session_state[sk] = dv
        return str(st.selectbox(label, opts, key=sk, format_func=_format_strain_or_amp_mode, help=h))

    if kind == 'list_float':
        if sk not in st.session_state:
            if cur is not None:
                try:
                    st.session_state[sk] = ', '.join(str(float(x)) for x in list(cur))
                except (TypeError, ValueError):
                    st.session_state[sk] = ''
            else:
                st.session_state[sk] = ''
        s = str(st.text_input(label, key=sk, help=h))
        return _parse_float_list(s)

    if kind == 'list_int':
        if sk not in st.session_state:
            if cur is not None:
                try:
                    st.session_state[sk] = ', '.join(str(int(x)) for x in list(cur))
                except (TypeError, ValueError):
                    st.session_state[sk] = ''
            else:
                st.session_state[sk] = ''
        s = str(st.text_input(label, key=sk, help=h))
        return _parse_int_list(s)

    if kind == 'json_dict':
        if sk not in st.session_state:
            st.session_state[sk] = json.dumps(cur, indent=2) if isinstance(cur, dict) else '{}'
        jh = (
            help_text
            if help_text is not None
            else (
                f'Advanced optional settings for `{name}`. **Leave `{{}}`** unless directed. '
                'Uses **JSON** text (`{{ "key": value }}`); see protocol docs or your PI.'
            )
        )
        s = str(st.text_area(label, height=120, key=sk, help=jh))
        try:
            return json.loads(s)
        except json.JSONDecodeError:
            _st_error_actions(
                f'Could not read `{name}`',
                ['Fix commas and brackets', 'Check quotes and braces', 'See ? help or protocol docs'],
            )
            return None

    return None


def _clear_fld_session_keys():
    """Drop procedure widget keys so the next render re-seeds from ``Bender`` (e.g. new config)."""
    for k in list(st.session_state.keys()):
        if k.startswith('fld_'):
            del st.session_state[k]


def _friendly_error_actions(err: Exception, *, action: str) -> tuple[str, list[str]]:
    """Headline and bulleted next steps (≤5 words each) for common failures."""
    msg = str(err or '').strip()
    low = msg.lower()
    if 'export_primary_h5 requires bender.outputfile or outputfile' in low:
        return (
            'HDF5 path not set.',
            ['Open section 2', 'Set data folder', 'Set file name', 'Apply data path'],
        )
    if 'please enter a note first' in low:
        return ('No note text.', ['Type your note', 'Click append again'])
    if 'file not found' in low and 'h5' in low:
        return ('File missing.', ['Pick file from list', 'Check folder path'])
    if 'selected file is not .h5' in low:
        return ('Not an HDF5 file.', ['Select a `.h5` file from the list'])
    if 'kaleido' in low and action == 'save_qc':
        return (
            'QC PNG export failed.',
            ['Run pip install kaleido', 'Or use HTML output'],
        )
    if action == 'save_h5':
        return (
            'HDF5 save failed.',
            ['Open Technical details', 'Fix error shown there'],
        )
    if action == 'save_qc':
        return (
            'QC figure save failed.',
            ['Open Technical details', 'Fix error shown there'],
        )
    if action == 'run_experiment':
        if '-200077' in msg or 'sampclk_rate' in low or 'daqmx_sampclk_rate' in low:
            return (
                'DAQ rejected sample rate.',
                ['Check DAQ rates in config', 'Check motion timeline dt', 'Open Technical details'],
            )
        if '-200560' in msg or 'wait until done' in low:
            return (
                'DAQ acquisition timed out.',
                ['Check USB connection', 'Close other DAQ software', 'Retry run experiment'],
            )
        if '-50103' in msg or 'resource reserved' in low or 'specified resource is reserved' in low:
            return (
                'DAQ still reserved elsewhere.',
                ['Close NI MAX panels', 'Close extra Streamlit tabs', 'Wait then retry'],
            )
        return (
            'Experiment did not run.',
            ['Open Technical details', 'Fix error shown there'],
        )
    return (
        'Something went wrong.',
        ['Open Technical details', 'Fix error shown there'],
    )


def _show_friendly_error(err: Exception, *, action: str):
    h, bullets = _friendly_error_actions(err, action=action)
    _st_error_actions(h, bullets)
    with st.expander('Technical details'):
        st.exception(err)


def _preview_error_actions(err_text: str) -> tuple[str, list[str]]:
    e = str(err_text or '').strip()
    low = e.lower()
    if 'set all_freqs' in low:
        return (
            'Preview needs frequencies.',
            ['Enter Frequencies (Hz) in motion controls', 'Click Refresh experiment preview'],
        )
    if 'set all_amps' in low:
        return (
            'Preview needs amplitudes.',
            ['Enter Amplitudes in motion controls', 'Click Refresh experiment preview'],
        )
    if 'needs duration' in low:
        return (
            'Preview needs duration.',
            ['Set Duration (s) for this protocol', 'Click Refresh experiment preview'],
        )
    if 'isometric_initial' in low or 'isometric_final' in low or 'isometric_num_steps' in low:
        return (
            'Preview needs isometric start/end/steps values.',
            ['Set isometric start, end, and number of steps', 'Click Refresh experiment preview'],
        )
    if 'step_change' in low and 'needs' in low:
        return (
            'Preview needs step-change parameters.',
            ['Set step-change frequencies, amplitudes, and cycles/step', 'Click Refresh experiment preview'],
        )
    return ('Preview failed.', ['Review preview error detail', 'Update required fields and refresh preview'])


def _needs_missing_calibration_confirmation(b: Bender) -> bool:
    """True when inertial-calibration use is requested but file path is absent or missing."""
    use_cal = bool(getattr(b, 'use_inertial_calibration', False))
    cal_file = str(getattr(b, 'inertial_calibration_file', '') or '').strip()
    if not use_cal:
        return False
    return (not cal_file) or (not os.path.isfile(cal_file))


def _append_note_to_h5_file(h5_path: str, note_text: str):
    """Append a dated note block to an existing .h5 file's post-trial notes attrs."""
    note = str(note_text or '').strip()
    if not note:
        raise ValueError('Please enter a note first.')
    if not os.path.isfile(h5_path):
        raise ValueError(f'File not found: {h5_path}')
    ts = pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')
    with h5py.File(h5_path, 'r+') as f:
        prev = str(f.attrs.get('post-trial notes', '') or '').strip()
        merged = note if not prev else f"{prev}\n\n--- Post-Experiment QC note ({ts}) ---\n{note}"
        f.attrs['post-trial notes'] = merged
        if '01_Metadata' in f:
            f['01_Metadata'].attrs['post-trial notes'] = merged


def _candidate_review_files(data_file_path: str):
    """Files in the data directory likely useful for review (h5/images/html)."""
    p = str(data_file_path or '').strip()
    if not p:
        return []
    d = os.path.dirname(p) or '.'
    if not os.path.isdir(d):
        return []
    exts = {'.h5', '.png', '.jpg', '.jpeg', '.webp', '.html'}
    out = []
    for n in sorted(os.listdir(d)):
        fp = os.path.join(d, n)
        if os.path.isfile(fp) and os.path.splitext(n.lower())[1] in exts:
            out.append(fp)
    return out


def _ensure_review_file_selection(files: list) -> None:
    if not files:
        return
    sel = st.session_state.get('gui_review_selected')
    if sel not in files:
        st.session_state['gui_review_selected'] = files[0]


def _read_qc_trial_index(b: Bender) -> int:
    """Which in-memory `trial_records` entry the QC figure uses (after Run / export)."""
    tr = list(getattr(b, 'trial_records', []) or [])
    n = len(tr)
    if n <= 1:
        return 0
    try:
        ix = int(st.session_state.get('gui_qc_trial_index'))
    except (TypeError, ValueError):
        ix = n - 1
    return int(max(0, min(ix, n - 1)))


def _qc_figure_base_path(b: Bender, selected_path: str, trial_idx: int) -> Optional[str]:
    """Output stem for `save_universal_qc_figure` when tying the QC file name to a chosen `.h5`."""
    p = str(selected_path or '').strip()
    if not p.lower().endswith('.h5'):
        return None
    proc = str(getattr(b, 'test_type', 'unknown'))
    stem = p[:-3]
    return stem + f'_{proc}_trial{int(trial_idx):03d}_qc'


def _compose_output_h5_path() -> str:
    """Full `.h5` path from **section 2** **Data folder** + **Data file name**."""
    folder = str(st.session_state.get('gui_data_folder') or '').strip()
    fn = str(st.session_state.get('gui_data_filename') or '').strip()
    # When a folder is selected, treat filename input as basename-only.
    if folder and fn:
        fn = os.path.basename(os.path.normpath(fn))
    if not fn:
        return ''
    if not fn.lower().endswith('.h5'):
        fn = fn + '.h5'
    if folder:
        return os.path.normpath(os.path.join(folder, fn))
    return os.path.normpath(fn)


def _normalize_data_filename_field() -> None:
    """Normalize filename field to basename when users paste/type a full path."""
    raw = str(st.session_state.get('gui_data_filename') or '').strip()
    if not raw:
        return
    base = os.path.basename(os.path.normpath(raw))
    if base and base != raw:
        st.session_state['gui_data_filename'] = base


def _split_path_cross_platform(p: str) -> Tuple[str, str]:
    """Return ``(dirname, basename)`` for a path that may use Windows or POSIX separators.

    Acquisition runs on Windows, but persisted signatures may be re-read on macOS where
    ``os.path`` would not split ``C:\\...\\file.h5`` (backslashes aren't separators on POSIX).
    Detect Windows-style paths (drive letter or backslash) and split with :mod:`ntpath`;
    otherwise use :mod:`posixpath`. Returns ``('', '')`` for empty input.
    """
    s = str(p or '').strip()
    if not s:
        return '', ''
    mod = _pathmod_cross_platform(s)
    nv = mod.normpath(s)
    return mod.dirname(nv), mod.basename(nv)


def _pathmod_cross_platform(p: str):
    """Return :mod:`ntpath` for Windows-style paths (drive letter or backslash), else :mod:`posixpath`."""
    s = str(p or '')
    is_windows = ('\\' in s) or (len(s) >= 2 and s[1] == ':' and s[0].isalpha())
    return ntpath if is_windows else posixpath


def _normpath_cross_platform(p: str) -> str:
    s = str(p or '').strip()
    if not s:
        return ''
    return _pathmod_cross_platform(s).normpath(s)


def _repair_data_path_fields_from_session() -> None:
    """Best-effort repair for data-path fields from persisted session artifacts."""
    folder = str(st.session_state.get('gui_data_folder') or '').strip()
    fn = str(st.session_state.get('gui_data_filename') or '').strip()
    if folder and fn:
        _, fn_base = _split_path_cross_platform(fn)
        st.session_state['gui_data_filename'] = fn_base or fn
        return

    sig = st.session_state.get('gui_data_path_applied_sig')
    if not isinstance(sig, (list, tuple)) or len(sig) != 2:
        return
    vals = [str(v or '').strip() for v in sig]
    vals = [v for v in vals if v]
    if not vals:
        return

    for v in vals:
        v_dir, v_base = _split_path_cross_platform(v)
        if v_base.lower().endswith('.h5') and (v_dir or ('\\' in v or '/' in v)):
            st.session_state['gui_data_folder'] = v_dir or folder
            st.session_state['gui_data_filename'] = v_base
            return

    if len(vals) == 2:
        a, b = vals
        a_dir, a_base = _split_path_cross_platform(a)
        b_dir, b_base = _split_path_cross_platform(b)
        if a_dir and b_base:
            st.session_state['gui_data_folder'] = _normpath_cross_platform(a)
            st.session_state['gui_data_filename'] = b_base
            return
        if b_dir and a_base:
            st.session_state['gui_data_folder'] = _normpath_cross_platform(b)
            st.session_state['gui_data_filename'] = a_base


def _section2_destination_incomplete() -> bool:
    """True when no usable output `.h5` destination is currently available."""
    folder = str(st.session_state.get('gui_data_folder') or '').strip()
    fn = str(st.session_state.get('gui_data_filename') or '').strip()
    if fn:
        # If filename includes a path, or folder+filename compose cleanly, we have a destination.
        if folder or os.path.dirname(fn):
            return False
        # Filename-only still resolves to a concrete relative `.h5` path.
        return False
    # Fall back to any already-bound experiment output path.
    b = st.session_state.get('bender')
    if b is not None:
        outp = str(getattr(b, 'outputfile', '') or '').strip()
        if outp:
            return False
    return True


def _output_path_anchor_for_review(b: Optional[Bender] = None) -> str:
    """Path used to locate the data directory (composed section 2 path, else ``b.outputfile``)."""
    p = _compose_output_h5_path().strip()
    if p:
        return p
    inst = b if b is not None else st.session_state.get('bender')
    if inst is not None:
        return str(getattr(inst, 'outputfile', '') or '').strip()
    return ''


def _normalize_config_module_name(raw: str) -> str:
    """Strip whitespace and optional ``.py`` suffix for ``importlib.import_module``."""
    s = str(raw or '').strip()
    if s.lower().endswith('.py'):
        s = s[:-3].strip()
    s = s.replace('\\', '.').replace('/', '.')
    s = '.'.join(part for part in s.split('.') if part)
    return s


def _project_config_py_path_for_stem(stem: str) -> str:
    """Absolute path to ``<stem>.py`` under the app folder when that layout applies.

    Checks the project root first (legacy flat layout, also what the resolver tests use),
    then ``templates/configs/`` where config modules live after the templates reorg.
    """
    mod = _normalize_config_module_name(stem)
    if not mod:
        return ''
    rel = mod.replace('.', os.sep) + '.py'
    root_path = os.path.normpath(os.path.join(_ROOT, rel))
    if os.path.isfile(root_path):
        return root_path
    configs_path = os.path.normpath(os.path.join(default_configs_dir(_ROOT), rel))
    if os.path.isfile(configs_path):
        return configs_path
    return root_path


def _resolve_hardware_config_import_target(raw: str) -> tuple[str, str]:
    """Return ``(module_stem, absolute_.py_path)`` for ``Bender(stem)`` import."""
    s = str(raw or '').strip()
    if not s:
        raise ValueError('Enter a config module name or `.py` path.')
    candidate = os.path.normpath(s)
    if s.lower().endswith('.py') or os.path.isfile(candidate):
        if not candidate.lower().endswith('.py'):
            raise ValueError('Hardware config path must be a `.py` file.')
        abs_path = os.path.normpath(os.path.abspath(candidate))
        if not os.path.isfile(abs_path):
            raise FileNotFoundError(f'Config file not found: `{abs_path}`')
        stem = os.path.splitext(os.path.basename(abs_path))[0]
        if not stem:
            raise ValueError('Invalid config file name.')
        return stem, abs_path
    stem = _normalize_config_module_name(s)
    if not stem:
        raise ValueError('Enter a config module name or `.py` path.')
    project_py = _project_config_py_path_for_stem(stem)
    if os.path.isfile(project_py):
        return stem, os.path.normpath(project_py)
    return stem, ''


def _ensure_config_dir_on_syspath_for_file(abs_py_path: str) -> None:
    if not abs_py_path:
        return
    cfg_dir = os.path.dirname(os.path.normpath(abs_py_path))
    if cfg_dir and cfg_dir not in sys.path:
        sys.path.insert(0, cfg_dir)


def _raw_mod_for_hardware_config_load(*, module_stem: str) -> str:
    """Prefer ``gui_load_cfg_file_path`` when it points at the intended module stem."""
    stem = _normalize_config_module_name(module_stem)
    cfg_path = str(st.session_state.get('gui_load_cfg_file_path') or '').strip()
    if cfg_path:
        norm = os.path.normpath(cfg_path)
        if os.path.isfile(norm):
            try:
                path_stem, _ = _resolve_hardware_config_import_target(norm)
                if path_stem == stem:
                    return norm
            except (OSError, ValueError):
                pass
    return stem


def _ensure_hw_config_session_defaults() -> None:
    """Keep ``gui_load_cfg_select`` valid when **section 1** is off-screen or session state is stale."""
    mods = discover_config_modules(_ROOT)
    cur_sel = str(st.session_state.get('gui_load_cfg_select') or '').strip()
    if cur_sel:
        return
    if not mods:
        return
    pick = str(st.session_state.get('cfg_mod') or '')
    if pick not in mods:
        pick = mods[0]
    st.session_state['gui_load_cfg_select'] = pick


def _selected_config_matches_bender(b: Bender, eff_raw: str) -> bool:
    """True if ``b`` matches the current selection (module stem and path when known)."""
    eff = _normalize_config_module_name(eff_raw)
    if not eff:
        return False
    loaded = _normalize_config_module_name(getattr(b, 'config_name', '') or '')
    if not loaded or loaded != eff:
        return False
    loaded_path = str(st.session_state.get('gui_loaded_cfg_abs_path') or '').strip()
    sel_path = str(st.session_state.get('gui_load_cfg_file_path') or '').strip()
    if sel_path:
        sel_norm = os.path.normpath(sel_path)
        if not os.path.isfile(sel_norm):
            return True
        if loaded_path:
            return _paths_equal_norm(sel_norm, loaded_path)
        project_py = _project_config_py_path_for_stem(eff)
        if project_py and os.path.isfile(project_py):
            return _paths_equal_norm(sel_norm, project_py)
        return True
    if loaded_path:
        project_py = _project_config_py_path_for_stem(eff)
        if project_py and os.path.isfile(project_py):
            return _paths_equal_norm(loaded_path, project_py)
    return True


def _apply_loaded_config_module(raw_mod: str) -> Optional[str]:
    """Instantiate ``Bender`` from a config module and refresh session. Returns an error message or ``None``."""
    try:
        _cm, _cfg_abs = _resolve_hardware_config_import_target(raw_mod)
    except FileNotFoundError as e:
        st.session_state.pop('bender', None)
        st.session_state.pop('gui_loaded_cfg_abs_path', None)
        return str(e)
    except ValueError as e:
        st.session_state.pop('bender', None)
        st.session_state.pop('gui_loaded_cfg_abs_path', None)
        return str(e)
    if _ROOT not in sys.path:
        sys.path.insert(0, _ROOT)
    _ensure_config_dir_on_syspath_for_file(_cfg_abs)
    try:
        st.session_state['bender'] = Bender(_cm)
        st.session_state['cfg_mod'] = _cm
        if _cfg_abs:
            st.session_state['gui_loaded_cfg_abs_path'] = _cfg_abs
        else:
            st.session_state.pop('gui_loaded_cfg_abs_path', None)
        b0 = st.session_state['bender']
        outp0 = str(getattr(b0, 'outputfile', '') or '').strip()
        if outp0:
            n0 = os.path.normpath(outp0)
            st.session_state['gui_pending_data_folder'] = os.path.dirname(n0) or ''
            st.session_state['gui_pending_data_filename'] = os.path.basename(n0)
        _meta0 = getattr(b0, 'h5_protocol_metadata', {}) or {}
        st.session_state['gui_pending_genus_species'] = str(_meta0.get('genus_species', '') or '')
        st.session_state['gui_pending_specimen_id'] = str(_meta0.get('specimen_id', '') or '')
        st.session_state['gui_pending_post_notes'] = str(getattr(b0, 'post_trial_notes', '') or '')
        _init_biometrics_session_state(b0, force=True)
        _clear_fld_session_keys()
        st.session_state.pop('gui_tpl_bio_done', None)
        return None
    except ImportError as e:
        st.session_state.pop('bender', None)
        st.session_state.pop('gui_loaded_cfg_abs_path', None)
        _path_hint = f'- Config file: `{_cfg_abs}`\n' if _cfg_abs else f'- App folder: `{_ROOT}`\n'
        return (
            f'Could not import config module `{_cm}`.\n\n'
            f'- Use a **module name** (e.g. `jimenez_bender_config_A`) or a full path to a `.py` file.\n'
            f'{_path_hint}'
            f'- Current working directory: `{os.getcwd()}`\n\n'
            f'Detail: {e}'
        )
    except Exception as e:
        st.session_state.pop('bender', None)
        st.session_state.pop('gui_loaded_cfg_abs_path', None)
        return f'{type(e).__name__}: {e}'


def _sec1_apply_composed_path_to_bender() -> Optional[str]:
    """Copy **Data folder** + **Data file name** onto ``bender.outputfile``. Returns an error message or ``None``."""
    b1 = st.session_state.get('bender')
    if b1 is None:
        return 'Load hardware configuration first.'
    outp = _compose_output_h5_path().strip()
    if not outp:
        return 'Set **Data file name** first.'
    b1.outputfile = outp
    _mark_data_path_applied()
    return None


def _paths_equal_norm(a: str, b: str) -> bool:
    if not str(a).strip() or not str(b).strip():
        return False
    try:
        return os.path.normcase(os.path.normpath(a)) == os.path.normcase(os.path.normpath(b))
    except Exception:
        return str(a).strip() == str(b).strip()


_BIO_APPLY_SESSION_KEYS = (
    'gui_genus_species',
    'gui_specimen_id',
    'bio_segment',
    'bio_fishmass',
    'bio_fishlen_TL',
    'bio_fishlen_SL',
    'bio_xsec_height',
    'bio_dvert',
    'bio_dhoriz',
    'bio_dclamp',
    'bio_xsec',
    'bio_muscle_depth',
    'bio_dbend',
    'bio_temp_room',
    'bio_temp_tank',
    'bio_prep_condition',
    'bio_use_theoretical_inertial',
    'bio_prof_L',
    'bio_prof_rho_preset',
    'bio_prof_rho',
    'bio_prof_ph',
    'bio_prof_pw',
    'bio_prof_dh',
    'bio_prof_dw',
    'bio_prof_clamp',
    'bio_prof_samples',
)


def _bio_fingerprint() -> tuple:
    return tuple(st.session_state.get(k) for k in _BIO_APPLY_SESSION_KEYS)


def _data_path_fingerprint() -> tuple:
    return (st.session_state.get('gui_data_folder'), st.session_state.get('gui_data_filename'))


def _procedure_fingerprint() -> tuple:
    tt = st.session_state.get('test_type_select')
    pairs = [(k, st.session_state.get(k)) for k in sorted(st.session_state.keys()) if k.startswith('fld_')]
    return (tt, tuple(pairs))


def _ensure_apply_tracking_bender(b: Bender) -> None:
    """Reset apply/dirty baselines when the in-memory ``Bender`` instance is replaced."""
    bid = id(b)
    if st.session_state.get('gui_apply_tracking_bender_id') != bid:
        st.session_state['gui_apply_tracking_bender_id'] = bid
        for _k in (
            'gui_bio_applied_sig',
            'gui_data_path_applied_sig',
            'gui_proc_applied_sig',
            'gui_bio_apply_invalidated',
            'gui_proc_apply_invalidated',
        ):
            st.session_state.pop(_k, None)


def _touch_bio_apply_baseline_if_clean() -> None:
    _inv = bool(st.session_state.get('gui_bio_apply_invalidated'))
    _has_sig = 'gui_bio_applied_sig' in st.session_state
    if _inv:
        # #region agent log
        _agent_debug_log(
            hypothesis_id='B',
            location='bender_streamlit_gui.py:_touch_bio_apply_baseline_if_clean',
            message='skip_touch_invalidated',
            data={'invalidated': True, 'has_sig': _has_sig},
        )
        # #endregion
        return
    if not _has_sig:
        st.session_state['gui_bio_applied_sig'] = _bio_fingerprint()
        # #region agent log
        _agent_debug_log(
            hypothesis_id='B',
            location='bender_streamlit_gui.py:_touch_bio_apply_baseline_if_clean',
            message='baseline_sig_set_without_apply',
            data={
                'fishmass': st.session_state.get('bio_fishmass'),
                'dclamp': st.session_state.get('bio_dclamp'),
                'xsec': st.session_state.get('bio_xsec'),
            },
        )
        # #endregion


def _touch_data_path_baseline_if_clean() -> None:
    if 'gui_data_path_applied_sig' not in st.session_state:
        st.session_state['gui_data_path_applied_sig'] = _data_path_fingerprint()


def _touch_proc_apply_baseline_if_clean() -> None:
    if st.session_state.get('gui_proc_apply_invalidated'):
        return
    if 'gui_proc_applied_sig' not in st.session_state:
        st.session_state['gui_proc_applied_sig'] = _procedure_fingerprint()


def _mark_bio_applied() -> None:
    st.session_state['gui_bio_apply_invalidated'] = False
    st.session_state['gui_bio_applied_sig'] = _bio_fingerprint()


def _mark_data_path_applied() -> None:
    st.session_state['gui_data_path_applied_sig'] = _data_path_fingerprint()


def _mark_procedure_applied() -> None:
    st.session_state['gui_proc_apply_invalidated'] = False
    st.session_state['gui_proc_applied_sig'] = _procedure_fingerprint()


def _bio_apply_dirty() -> bool:
    if st.session_state.get('bender') is None:
        return False
    if st.session_state.get('gui_bio_apply_invalidated'):
        return True
    if 'gui_bio_applied_sig' not in st.session_state:
        return False
    return _bio_fingerprint() != st.session_state['gui_bio_applied_sig']


def _bio_apply_dirty_reason() -> str:
    if st.session_state.get('bender') is None:
        return 'no_bender'
    if st.session_state.get('gui_bio_apply_invalidated'):
        return 'invalidated'
    if 'gui_bio_applied_sig' not in st.session_state:
        return 'no_applied_sig'
    if _bio_fingerprint() != st.session_state['gui_bio_applied_sig']:
        return 'fingerprint_mismatch'
    return 'clean'


def _data_path_apply_dirty() -> bool:
    if 'gui_data_path_applied_sig' not in st.session_state:
        return False
    return _data_path_fingerprint() != st.session_state['gui_data_path_applied_sig']


def _procedure_apply_dirty() -> bool:
    if st.session_state.get('gui_proc_apply_invalidated'):
        return True
    if 'gui_proc_applied_sig' not in st.session_state:
        return False
    return _procedure_fingerprint() != st.session_state['gui_proc_applied_sig']


def _procedure_ready_for_run() -> tuple[bool, str]:
    """Procedure committed via Apply / preview — not silently on Run."""
    if _procedure_apply_dirty():
        return (
            False,
            'Procedure fields changed — click **Apply procedure** or **Refresh experiment preview** first.',
        )
    if not bool(st.session_state.get('gui_protocol_confirmed')):
        return (
            False,
            'Click **Apply procedure** or **Refresh experiment preview** before **Run experiment**.',
        )
    return True, ''


def _soft_apply_reminder() -> None:
    # Intentionally silent: setup readiness is shown in the sidebar checklist.
    return


def _session_float(key: str) -> Optional[float]:
    try:
        v = st.session_state.get(key)
        if v is None:
            return None
        x = float(v)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def _measurements_confirmed_for_checklist() -> bool:
    """True when biometrics/measurements look intentionally confirmed, not just default-populated."""
    if not bool(st.session_state.get('gui_measurements_confirmed')):
        return False
    if 'gui_bio_applied_sig' not in st.session_state:
        return False
    if _bio_apply_dirty():
        return False
    if not str(st.session_state.get('gui_specimen_id') or '').strip():
        return False
    if not str(st.session_state.get('gui_genus_species') or '').strip():
        return False
    m = _session_float('bio_fishmass')
    if m is None or m <= 0:
        return False
    dc = _session_float('bio_dclamp')
    if dc is None or dc <= 0:
        return False
    xw = _session_float('bio_xsec')
    if xw is None or xw <= 0:
        return False
    return True


def _sanitize_stale_run_state() -> None:
    """Clear stuck run flags when no acquisition can be active (e.g. after refresh)."""
    if st.session_state.get('bender') is None:
        st.session_state['gui_run_in_progress'] = False
        st.session_state['gui_run_pending_confirm'] = False
        st.session_state['gui_run_soft_warnings'] = []


def _run_button_state() -> tuple[bool, str]:
    """Whether **Run experiment** should be disabled and short help text."""
    if st.session_state.get('bender') is None:
        return True, 'Load hardware configuration in Setup (section 1) first.'
    if bool(st.session_state.get('gui_run_in_progress', False)):
        return True, 'Run in progress — wait for acquisition to finish.'
    return False, 'Starts NI-DAQ acquisition (simulation when enabled on other routes).'


def _measurements_fields_ok() -> bool:
    if not str(st.session_state.get('gui_specimen_id') or '').strip():
        return False
    if not str(st.session_state.get('gui_genus_species') or '').strip():
        return False
    m = _session_float('bio_fishmass')
    if m is None or m <= 0:
        return False
    dc = _session_float('bio_dclamp')
    if dc is None or dc <= 0:
        return False
    xw = _session_float('bio_xsec')
    if xw is None or xw <= 0:
        return False
    return True


def _setup_ready(b: Optional[Bender]) -> bool:
    if b is None:
        return False
    if _section2_destination_incomplete():
        return False
    outp = str(getattr(b, 'outputfile', '') or '').strip()
    composed = _compose_output_h5_path().strip()
    if not outp and not composed:
        return False
    if composed and outp and not _paths_equal_norm(outp, composed) and _data_path_apply_dirty():
        return False
    _cfg_sel = _normalize_config_module_name(str(st.session_state.get('gui_load_cfg_select') or ''))
    if not _cfg_sel:
        return False
    _cfg_mode = str(st.session_state.get('gui_config_setup_mode', 'Load existing'))
    if _cfg_mode == 'Load existing' and not _selected_config_matches_bender(b, _cfg_sel):
        return False
    return True


def _protocol_ready(b: Optional[Bender], tt: str) -> bool:
    if b is None:
        return False
    if _procedure_apply_dirty():
        return False
    try:
        rep = b.validate_dispatch_setup(test_type=tt)
        return bool(rep.get('ok', False))
    except Exception:
        return False


def _workflow_ready_state(b: Optional[Bender], tt: str) -> dict[str, Any]:
    """Unified readiness for sidebar checklist and Run UX (scratch-oriented)."""
    run_disabled, run_reason = _run_button_state()
    setup_ok = _setup_ready(b)
    measurements_ok = _measurements_fields_ok() and not _bio_apply_dirty()
    # #region agent log
    _agent_debug_log(
        hypothesis_id='A',
        location='bender_streamlit_gui.py:_workflow_ready_state',
        message='checklist_measurements',
        data={
            'route': _nav_route(),
            'stepwise_step': _stepwise_step() if _nav_route() == 'stepwise' else None,
            'fields_ok': _measurements_fields_ok(),
            'bio_dirty': _bio_apply_dirty(),
            'dirty_reason': _bio_apply_dirty_reason(),
            'measurements_ok': measurements_ok,
            'fishmass': st.session_state.get('bio_fishmass'),
            'dclamp': st.session_state.get('bio_dclamp'),
            'bender_dclamp': getattr(b, 'dclamp', None) if b is not None else None,
        },
    )
    # #endregion
    if measurements_ok and b is not None:
        st.session_state.setdefault('gui_measurements_confirmed', True)
    protocol_ok = _protocol_ready(b, tt)
    if protocol_ok:
        st.session_state.setdefault('gui_protocol_confirmed', True)
    if setup_ok:
        st.session_state.setdefault('gui_setup_confirmed', True)
    return {
        'setup_ok': setup_ok,
        'measurements_ok': measurements_ok,
        'protocol_ok': protocol_ok,
        'run_disabled': run_disabled,
        'run_disabled_reason': run_reason,
    }


def _setup_confirmed_for_checklist(b: Optional[Bender]) -> bool:
    tt = str(st.session_state.get('test_type_select') or getattr(b, 'test_type', '') or 'dynamic')
    return _workflow_ready_state(b, tt)['setup_ok']


def _measurements_confirmed_for_checklist() -> bool:
    b = st.session_state.get('bender')
    tt = str(st.session_state.get('test_type_select') or getattr(b, 'test_type', '') or 'dynamic')
    return _workflow_ready_state(b, tt)['measurements_ok']


_CHK_SEC_DATA = '2 · Data path'
_CHK_SEC_BIO = '3 · Measurements'
_CHK_SEC_EXP = '4–6 · Experiment'


def _fld_raw_str(name: str) -> str:
    v = st.session_state.get(_widget_key(name), '')
    if v is None:
        return ''
    return str(v).strip()


def _collect_experiment_form_status_messages(tt: str) -> list[str]:
    """Warnings from **Procedure fields** widgets (`fld_*`), including blank frequencies before **Apply**."""
    msgs: list[str] = []
    if tt in MOTION_TYPES and tt != 'step_change':
        raw_f = _fld_raw_str('all_freqs')
        if not raw_f:
            msgs.append('Frequencies field is blank.')
        else:
            parsed = _parse_float_list(raw_f)
            if not parsed:
                msgs.append('Frequencies not parseable.')
            else:
                if any(not math.isfinite(x) for x in parsed):
                    msgs.append('Frequencies must be finite.')
                if any(x <= 0 for x in parsed):
                    msgs.append('Frequencies must be > 0 Hz.')
        raw_a = _fld_raw_str('all_amps')
        if not raw_a:
            msgs.append('Amplitudes field is blank.')
        else:
            ap = _parse_float_list(raw_a)
            if not ap:
                msgs.append('Amplitudes not parseable.')
        if tt in ('frequency_sweep', 'frequency_step', 'curvature_step'):
            skd = _widget_key('duration')
            dv = None
            if skd in st.session_state:
                try:
                    dv = float(st.session_state[skd])
                except (TypeError, ValueError):
                    dv = None
            if dv is None or not math.isfinite(dv) or dv <= 0:
                msgs.append('Duration (s) must be > 0.')
        if tt == 'dynamic':
            sk = _widget_key('cycles_per_step')
            if sk in st.session_state:
                try:
                    cps = int(st.session_state[sk])
                    if cps <= 0:
                        msgs.append('Cycles per step must be ≥ 1.')
                except (TypeError, ValueError):
                    msgs.append('Cycles per step must be an integer ≥ 1.')
    elif tt == 'step_change':
        for fname, label in (
            ('step_change_frequencies', 'Step-change frequencies'),
            ('step_change_curves', 'Step-change amplitudes'),
            ('step_change_cycles_per_step', 'Cycles per step'),
        ):
            raw = _fld_raw_str(fname)
            if not raw:
                msgs.append(f'{label} field is blank.')
                continue
            if fname == 'step_change_cycles_per_step':
                pi = _parse_int_list(raw)
                if not pi:
                    msgs.append(f'{label} not parseable.')
                elif any(x < 1 for x in pi):
                    msgs.append(f'{label}: values must be ≥ 1.')
            else:
                pf = _parse_float_list(raw)
                if not pf:
                    msgs.append(f'{label} not parseable.')
                elif fname == 'step_change_frequencies' and any(not math.isfinite(x) or x <= 0 for x in pf):
                    msgs.append('Step-change frequencies must be finite and > 0 Hz.')
    if tt == 'isometric':
        sk = _widget_key('isometric_num_steps')
        if sk in st.session_state:
            try:
                if int(st.session_state[sk]) < 1:
                    msgs.append('Isometric steps must be ≥ 1.')
            except (TypeError, ValueError):
                pass
    if tt == 'isovelocity':
        sk = _widget_key('isovelocity_num_steps')
        if sk in st.session_state:
            try:
                if int(st.session_state[sk]) < 1:
                    msgs.append('Isovelocity steps must be ≥ 1.')
            except (TypeError, ValueError):
                pass
    return msgs


def _collect_check_tuples(b: Bender) -> list[tuple[str, str]]:
    """Return ``(section_label, message)`` for sidebar **Status check**: path, biometrics, experiment form + Bender validation."""
    out: list[tuple[str, str]] = []
    tt = str(st.session_state.get('test_type_select') or getattr(b, 'test_type', '') or 'dynamic')

    df = str(st.session_state.get('gui_data_folder') or '').strip()
    fn = str(st.session_state.get('gui_data_filename') or '').strip()
    if not df:
        out.append((_CHK_SEC_DATA, 'Data folder is empty.'))
    else:
        try:
            if not os.path.isdir(os.path.normpath(df)):
                out.append((_CHK_SEC_DATA, 'Data folder is not a valid directory.'))
        except Exception:
            out.append((_CHK_SEC_DATA, 'Data folder path is invalid.'))
    if not fn:
        out.append((_CHK_SEC_DATA, 'Data file name is empty.'))

    composed = _compose_output_h5_path().strip()
    applied = str(getattr(b, 'outputfile', '') or '').strip()
    if composed and applied and not _paths_equal_norm(applied, composed):
        out.append((_CHK_SEC_DATA, 'Form path ≠ experiment object (section 2).'))
    elif composed and not applied:
        out.append((_CHK_SEC_DATA, 'Data path not set on experiment object (section 2).'))

    if not str(st.session_state.get('gui_specimen_id') or '').strip():
        out.append((_CHK_SEC_BIO, 'Specimen ID is blank.'))
    if not str(st.session_state.get('gui_genus_species') or '').strip():
        out.append((_CHK_SEC_BIO, 'Genus-species is blank.'))

    m = _session_float('bio_fishmass')
    if m is not None:
        if m <= 0:
            out.append((_CHK_SEC_BIO, 'Mass is zero or negative.'))
        elif m < 1.0:
            out.append((_CHK_SEC_BIO, f'Mass {m:g} g (check units).'))

    _BIO_INTRINSIC_MM = (
        ('Whole-body TL', 'bio_fishlen_TL'),
        ('Whole-body SL', 'bio_fishlen_SL'),
    )
    _BIO_CLAMP_GEOMETRY_MM = (
        ('Clamp spacing (dclamp)', 'bio_dclamp'),
        ('Cross-section width', 'bio_xsec'),
        ('Muscle depth from surface', 'bio_muscle_depth'),
        ('Cross-section height', 'bio_xsec_height'),
    )
    _BIO_PROFILE_OUTLINE_MM = (('Profile outline length', 'bio_prof_L'),)
    for label, key in _BIO_INTRINSIC_MM + _BIO_CLAMP_GEOMETRY_MM + _BIO_PROFILE_OUTLINE_MM:
        v = _session_float(key)
        if v is None:
            continue
        if v <= 0:
            out.append((_CHK_SEC_BIO, f'{label}: invalid ({v:g} mm).'))
        elif v < 1.0 and key != 'bio_muscle_depth':
            out.append((_CHK_SEC_BIO, f'{label}: {v:g} mm (check units).'))

    xw_chk = _session_float('bio_xsec')
    md_chk = _session_float('bio_muscle_depth')
    if md_chk is not None and xw_chk is not None and xw_chk > 0:
        if md_chk < 0:
            out.append((_CHK_SEC_BIO, f'Muscle depth: invalid ({md_chk:g} mm).'))
        elif md_chk >= xw_chk / 2.0:
            out.append(
                (
                    _CHK_SEC_BIO,
                    f'Muscle depth ({md_chk:g} mm) must be < half width ({xw_chk / 2.0:g} mm).',
                )
            )

    v_dbend = _session_float('bio_dbend')
    if v_dbend is not None:
        if v_dbend < 0:
            out.append((_CHK_SEC_BIO, f'Segment center distance: invalid ({v_dbend:g} mm).'))
        elif 0 < v_dbend < 1.0:
            out.append((_CHK_SEC_BIO, f'Segment center distance: {v_dbend:g} mm (check units).'))

    for label, key in (
        ('Vertical offset (dvert)', 'bio_dvert'),
        ('Horizontal offset (dhoriz)', 'bio_dhoriz'),
    ):
        v = _session_float(key)
        if v is None:
            continue
        if v < 0:
            out.append((_CHK_SEC_BIO, f'{label}: invalid ({v:g} mm).'))
        elif 0 < v < 1.0:
            out.append((_CHK_SEC_BIO, f'{label}: {v:g} mm (check units).'))

    for label, key in (
        ('Profile proximal H', 'bio_prof_ph'),
        ('Profile proximal W', 'bio_prof_pw'),
        ('Profile distal H', 'bio_prof_dh'),
        ('Profile distal W', 'bio_prof_dw'),
    ):
        v = _session_float(key)
        if v is None:
            continue
        if v <= 0:
            out.append((_CHK_SEC_BIO, f'{label}: invalid ({v:g} mm).'))
        elif v < 1.0:
            out.append((_CHK_SEC_BIO, f'{label}: {v:g} mm (check units).'))

    v_pclamp = _session_float('bio_prof_clamp')
    if v_pclamp is not None:
        if v_pclamp < 0:
            out.append((_CHK_SEC_BIO, f'Axis-clamp distance (profile): invalid ({v_pclamp:g} mm).'))
        elif 0 < v_pclamp < 1.0:
            out.append((_CHK_SEC_BIO, f'Axis-clamp distance (profile): {v_pclamp:g} mm (check units).'))

    for msg in _collect_experiment_form_status_messages(tt):
        out.append((_CHK_SEC_EXP, msg))

    try:
        rep = b.validate_dispatch_setup(test_type=tt)
        for m in rep.get('missing') or []:
            out.append((_CHK_SEC_EXP, f'{tt}: {m}'))
    except Exception:
        pass

    # Dedupe (section, message) while keeping order
    seen: set[tuple[str, str]] = set()
    uniq: list[tuple[str, str]] = []
    for t in out:
        if t not in seen:
            seen.add(t)
            uniq.append(t)
    return uniq


def _setup_confirmed_for_checklist(b: Optional[Bender]) -> bool:
    if b is None:
        return False
    if not bool(st.session_state.get('gui_setup_confirmed')):
        return False
    if _section2_destination_incomplete() or _data_path_apply_dirty():
        return False
    _cfg_sel = _normalize_config_module_name(str(st.session_state.get('gui_load_cfg_select') or ''))
    if not _cfg_sel:
        return False
    _cfg_mode = str(st.session_state.get('gui_config_setup_mode', 'Load existing'))
    if _cfg_mode == 'Load existing' and not _selected_config_matches_bender(b, _cfg_sel):
        return False
    return True


def _protocol_confirmed_for_checklist(b: Optional[Bender], checks_by_sec: dict[str, list[str]]) -> bool:
    """Protocol step complete when Bender validates; form-widget warnings do not block."""
    del checks_by_sec  # retained for call-site compatibility
    tt = str(st.session_state.get('test_type_select') or getattr(b, 'test_type', '') or 'dynamic')
    return _workflow_ready_state(b, tt)['protocol_ok']


def _refresh_confirmation_flags() -> None:
    if st.session_state.get('bender') is None:
        st.session_state['gui_setup_confirmed'] = False
        st.session_state['gui_measurements_confirmed'] = False
        st.session_state['gui_protocol_confirmed'] = False
        return
    if _data_path_apply_dirty() or _section2_destination_incomplete():
        st.session_state['gui_setup_confirmed'] = False
    if _bio_apply_dirty():
        st.session_state['gui_measurements_confirmed'] = False
    if _procedure_apply_dirty():
        st.session_state['gui_protocol_confirmed'] = False


def _mark_review_data_used() -> None:
    st.session_state['gui_review_data_used'] = True


def _arm_acquired_trial_review(h5_path, qc_path) -> None:
    """Flag a just-acquired trial for the inline keep/delete review.

    The raw ``.h5`` is already persisted on disk by the run; this only stores the paths so
    the review plot + Keep/Delete controls render. The save is never gated on this review.
    """
    st.session_state['gui_review_pending'] = True
    st.session_state['gui_review_h5_path'] = str(h5_path or '')
    st.session_state['gui_review_qc_path'] = str(qc_path or '')


def _render_acquired_trial_review(b: Bender) -> None:
    """Inline review of the last acquired trial: raw (and corrected, when present) torque, then Keep/Delete."""
    if not bool(st.session_state.get('gui_review_pending', False)):
        return
    st.divider()
    st.subheader('Review acquired trial')
    rev_h5 = str(st.session_state.get('gui_review_h5_path') or '')
    if rev_h5:
        st.caption(f'Saved data file (already on disk): `{rev_h5}`')
    try:
        rev_fig, _ = build_universal_qc_figure(b)
        st.plotly_chart(rev_fig, use_container_width=True)
        st.caption('Raw torque shown; inertia-corrected torque appears when available.')
    except Exception as rev_exc:
        st.warning(f'Could not build review plot: {type(rev_exc).__name__}: {rev_exc}')
    keep_col, del_col = st.columns(2)
    with keep_col:
        if st.button('Keep data', key='gui_review_keep', type='primary', use_container_width=True):
            st.session_state['gui_review_pending'] = False
            st.toast('Kept acquired data.')
            st.rerun()
    with del_col:
        if st.button('Delete data', key='gui_review_delete', type='secondary', use_container_width=True):
            deleted = []
            for p in (rev_h5, str(st.session_state.get('gui_review_qc_path') or '')):
                try:
                    if p and os.path.isfile(p):
                        os.remove(p)
                        deleted.append(p)
                except OSError as del_exc:
                    st.warning(f'Could not delete `{p}`: {del_exc}')
            st.session_state['gui_review_pending'] = False
            if deleted:
                st.warning('Deleted: ' + ', '.join(f'`{os.path.basename(p)}`' for p in deleted))
            else:
                st.info('Nothing to delete (file already removed).')
            st.rerun()


def _build_checklist_fix_lines(
    *,
    b: Optional[Bender],
    setup_ok: bool,
    measurements_ok: bool,
    protocol_ok: bool,
    checks_by_sec: Optional[dict[str, list[str]]] = None,
    tt: str = 'dynamic',
) -> list[str]:
    lines: list[str] = []
    if not setup_ok:
        lines.append('Setup: complete Step 1 and click Apply setup.')
    if not measurements_ok:
        lines.append('Measurements: enter required values and click Apply specimen and Apply clamp geometry & inertial correction.')
    if not protocol_ok:
        lines.append('Protocol/Run: fill required procedure fields and click Apply procedure or Refresh experiment preview.')
    if b is None:
        lines.append('Hardware: load a hardware configuration in Setup.')
    if checks_by_sec:
        for sec in (_CHK_SEC_DATA, _CHK_SEC_BIO, _CHK_SEC_EXP):
            for msg in (checks_by_sec.get(sec) or [])[:2]:
                lines.append(f'{sec}: {msg}')
    elif b is not None and not protocol_ok:
        try:
            for msg in (b.validate_dispatch_setup(test_type=tt).get('missing') or [])[:3]:
                lines.append(f'Protocol: missing {msg}')
        except Exception:
            pass
    return lines


def _check_sections_for_sidebar() -> Optional[frozenset[str]]:
    """Which section labels to show in **Status check**. ``None`` = all sections (full workflow)."""
    if _nav_route() != 'stepwise':
        return None
    s = _stepwise_step()
    if s <= 0:
        return frozenset()
    if s == 1:
        return frozenset({_CHK_SEC_DATA})
    if s == 2:
        return frozenset({_CHK_SEC_BIO})
    if s == 3:
        return frozenset({_CHK_SEC_EXP})
    # Step 5 (index 4): full session review before plots
    return None


def _render_sidebar_stepwise_progress() -> None:
    """Stepwise HW + path summary (sidebar)."""
    b = st.session_state.get('bender')
    st.markdown('**Progress**')
    if b is None:
        st.caption('Step 1: load hardware.')
        return
    st.caption(f"Module `{getattr(b, 'config_name', '?')}`")
    composed = _compose_output_h5_path().strip()
    applied = str(getattr(b, 'outputfile', '') or '').strip()
    if not applied:
        if composed:
            st.success(f"Selected HDF5 → `{os.path.basename(composed)}` (not applied yet — click **Apply setup**)")
        else:
            st.caption('Step 2: apply data path.')
    elif not composed:
        st.caption('Step 2: set a file name.')
    elif _paths_equal_norm(applied, composed):
        st.caption(f"HDF5 → `{os.path.basename(applied)}`")
    else:
        st.success(
            f"Selected HDF5 → `{os.path.basename(composed)}` (selected, not applied yet — click **Apply setup**)"
        )
    st.divider()


def _render_sidebar_input_checks() -> None:
    """Status check: suspicious or missing inputs when hardware is loaded; grouped by section, filtered in stepwise."""
    b = st.session_state.get('bender')
    if b is None:
        return
    items = _collect_check_tuples(b)
    allowed = _check_sections_for_sidebar()
    if allowed is not None:
        items = [(sec, msg) for sec, msg in items if sec in allowed]
    if not items:
        return
    st.markdown('**Status check**')
    if allowed is None:
        st.caption('Path, biometrics, and experiment settings worth a second look (full workflow).')
    else:
        st.caption(f'For **{_STEPWISE_TAB_LABELS[_stepwise_step()]}** only.')
    groups: dict[str, list[str]] = {}
    for sec, msg in items[:24]:
        groups.setdefault(sec, []).append(msg)
    parts: list[str] = []
    for sec, msgs in groups.items():
        parts.append(f'**{sec}**\n' + '\n'.join(f'- {m}' for m in msgs))
    st.warning('\n\n'.join(parts))
    if len(items) > 24:
        st.caption(f'… and {len(items) - 24} more.')
    if _CHK_SEC_EXP in groups:
        st.caption(
            'Experiment rows combine **Procedure fields** (form) and values on the **Bender** object after **Apply**.'
        )
    st.divider()


def _render_simulation_sidebar() -> None:
    """Simulation mode: numpy stand-in for DAQ; material affects cantilever stiffness in :meth:`Bender.run`."""
    st.session_state.setdefault('gui_simulation_mode', False)
    st.session_state.setdefault('gui_simulation_material', 'polyurethane')
    st.markdown('### Simulation mode')
    st.caption(
        'Run the app **without NI hardware**: a solid 25.4 mm OD tube cantilever model drives synthetic '
        'force/torque channels (polyurethane vs silicone Young\'s modulus).'
    )
    with st.container(border=True):
        st.markdown('<div class="bnd-simulation-panel" aria-hidden="true"></div>', unsafe_allow_html=True)
        st.checkbox(
            'Enable simulation mode (no DAQ)',
            key='gui_simulation_mode',
            help='Bypasses NI-DAQmx in **Run experiment**; use **Refresh experiment preview** in **Procedure fields** for force–displacement plots.',
        )
        _sim_on = bool(st.session_state.get('gui_simulation_mode', False))
        st.selectbox(
            'Simulated tube material',
            options=['polyurethane', 'silicone'],
            format_func=lambda k: {
                'polyurethane': 'Polyurethane (E ≈ 35 MPa — stiffer)',
                'silicone': 'Silicone (E ≈ 3 MPa — softer)',
            }[k],
            key='gui_simulation_material',
            disabled=not _sim_on,
            help='Sets E in k = 3EI/L³ (viscous term η·dδ/dt adds mild rate dependence).',
        )
    b = st.session_state.get('bender')
    if b is not None:
        b.simulation_mode = bool(st.session_state.get('gui_simulation_mode', False))
        b.simulation_material = str(st.session_state.get('gui_simulation_material', 'polyurethane'))
    st.divider()


def _seed_cfg_build_from_source_config(source_config: str) -> None:
    """Fill ``gui_cfg_bld_*`` widget defaults from an existing config module."""
    try:
        if _ROOT not in sys.path:
            sys.path.insert(0, _ROOT)
        base_stem = _normalize_config_module_name(source_config)
        build_path = str(st.session_state.get('gui_cfg_build_base_path') or '').strip()
        if build_path and os.path.isfile(os.path.normpath(build_path)):
            try:
                path_stem, abs_path = _resolve_hardware_config_import_target(build_path)
                if path_stem == base_stem:
                    _ensure_config_dir_on_syspath_for_file(abs_path)
            except (OSError, ValueError):
                pass
        d = read_base_defaults(base_stem)
    except Exception:
        return
    st.session_state['gui_cfg_bld_forcetorque_calibration_file'] = str(d['forcetorque_calibration_file'])
    st.session_state['gui_cfg_bld_positive_motor_direction'] = str(d['positive_motor_direction'])
    st.session_state['gui_cfg_bld_specimen_lateral_index'] = int(d['specimen_lateral_index_on_positive_motor_side'])
    st.session_state['gui_cfg_bld_primary_bending_axis'] = str(d['primary_bending_axis'])
    st.session_state['gui_cfg_bld_device_name'] = str(d['device_name'])
    st.session_state['gui_cfg_bld_daq_ai_sr'] = float(d['daq_ai_sample_rate_hz'])
    st.session_state['gui_cfg_bld_daq_ao_sr'] = float(d['daq_ao_do_sample_rate_hz'])
    st.session_state['gui_cfg_bld_motor_steps'] = int(d['motor_full_steps_per_rev'])
    st.session_state['gui_cfg_bld_motor_gear'] = int(d['motor_gear_ratio'])
    st.session_state['gui_cfg_bld_stim_channels'] = ', '.join(str(x) for x in d['stim_channels'])
    st.session_state['gui_cfg_bld_motor_port'] = str(d['motor_port'])
    st.session_state['gui_cfg_bld_encoder_chan'] = str(d['encoder_chan'])
    st.session_state['gui_cfg_bld_SG_chan'] = ', '.join(str(x) for x in d['SG_chan'])
    st.session_state['gui_cfg_bld_SG_name'] = ', '.join(str(x) for x in d['SG_name'])
    st.session_state['gui_cfg_bld_S1side'] = str(d['S1side'])
    st.session_state['gui_cfg_bld_S2side'] = str(d['S2side'])
    st.session_state['gui_cfg_bld_use_sono'] = bool(d['use_sono'])
    st.session_state['gui_cfg_bld_sono_channel'] = ', '.join(str(x) for x in d['sono_channel'])
    st.session_state['gui_cfg_bld_sono_name'] = ', '.join(str(x) for x in d['sono_name'])
    st.session_state['gui_cfg_bld_encoder_ppr'] = int(d['encoder_pulses_per_rev'])
    st.session_state['gui_cfg_bld_motor_axis'] = str(d['motor_axis'])
    st.session_state['gui_cfg_bld_bending_axis_sensor'] = str(d['bending_axis_sensor'])
    st.session_state['gui_cfg_bld_bending_axis_specimen'] = str(d['bending_axis_specimen'])
    st.session_state['gui_cfg_bld_stim_monitor_chan'] = ', '.join(str(x) for x in d['stim_monitor_chan'])
    st.session_state['gui_cfg_bld_stim_monitor_name'] = ', '.join(str(x) for x in d['stim_monitor_name'])
    st.session_state['gui_cfg_bld_sono_internal_samplefreq'] = int(d['sono_internal_samplefreq'])
    st.session_state['gui_cfg_bld_sono_cal_left'] = ', '.join(str(x) for x in d['sono_cal_left'])
    st.session_state['gui_cfg_bld_sono_cal_right'] = ', '.join(str(x) for x in d['sono_cal_right'])
    st.session_state['gui_cfg_bld_amp_step_vel'] = int(d['amp_step_vel'])
    rm = str(d['ramp_mode_default'])
    st.session_state['gui_cfg_bld_ramp_mode_default'] = rm if rm in ('linear', 'exponential') else 'linear'
    st.session_state['gui_cfg_bld_waitbefore'] = float(d['waitbefore'])
    st.session_state['gui_cfg_bld_waitafter'] = float(d['waitafter'])
    st.session_state['gui_cfg_bld_rampdur'] = float(d['rampdur'])
    st.session_state['gui_cfg_bld_prepoststim_dur'] = float(d['prepoststim_dur'])
    st.session_state['gui_cfg_bld_prepoststim_sep'] = float(d['prepoststim_sep'])
    st.session_state['gui_cfg_bld_prestim_time'] = float(d['prestim_time'])
    st.session_state['gui_cfg_bld_poststim_time'] = float(d['poststim_time'])


def _maybe_seed_cfg_build_fields() -> None:
    source_config = str(st.session_state.get('gui_cfg_build_base') or 'jimenez_bender_config_A')
    if st.session_state.get('gui_cfg_build_seeded_for') == source_config:
        return
    _seed_cfg_build_from_source_config(source_config)
    st.session_state['gui_cfg_build_seeded_for'] = source_config


def _flush_pending_cfg_build_base() -> None:
    """Apply deferred base-module selection before build-form widgets mount."""
    if 'gui_pending_cfg_build_base' not in st.session_state:
        return
    base = str(st.session_state.pop('gui_pending_cfg_build_base') or '').strip()
    if not base:
        return
    if base != str(st.session_state.get('gui_cfg_build_base') or ''):
        st.session_state['gui_cfg_build_base'] = base
        st.session_state.pop('gui_cfg_build_seeded_for', None)


def _flush_pending_load_config_session():
    """Apply deferred session updates before widgets bind to these keys (Streamlit restriction)."""
    if 'gui_pending_genus_species' in st.session_state:
        st.session_state['gui_genus_species'] = st.session_state.pop('gui_pending_genus_species')
    if 'gui_pending_specimen_id' in st.session_state:
        st.session_state['gui_specimen_id'] = st.session_state.pop('gui_pending_specimen_id')
    if 'gui_pending_data_folder' in st.session_state:
        _pending_folder = str(st.session_state.pop('gui_pending_data_folder') or '').strip()
        _cur_folder = str(st.session_state.get('gui_data_folder') or '').strip()
        # Do not clobber a user-picked folder with an empty pending value.
        if _pending_folder or not _cur_folder:
            st.session_state['gui_data_folder'] = _pending_folder
    if 'gui_pending_data_filename' in st.session_state:
        _pending_file = str(st.session_state.pop('gui_pending_data_filename') or '').strip()
        _cur_file = str(st.session_state.get('gui_data_filename') or '').strip()
        if _pending_file or not _cur_file:
            st.session_state['gui_data_filename'] = _pending_file
    if 'gui_pending_post_notes' in st.session_state:
        st.session_state['gui_post_notes'] = st.session_state.pop('gui_pending_post_notes')
    if 'gui_pending_load_cfg_file_path' in st.session_state:
        st.session_state['gui_load_cfg_file_path'] = st.session_state.pop('gui_pending_load_cfg_file_path')


def _ensure_gui_data_path_session_keys():
    """Seed folder/filename from legacy single-field ``gui_outputfile`` if needed."""
    if 'gui_data_folder' in st.session_state and 'gui_data_filename' in st.session_state:
        return
    leg = str(st.session_state.get('gui_outputfile', '') or '').strip()
    if leg:
        norm = os.path.normpath(leg)
        folder = os.path.dirname(norm) or ''
        filename = os.path.basename(norm)
        if folder:
            st.session_state['gui_data_folder'] = folder
        if filename:
            st.session_state['gui_data_filename'] = filename
    # else: leave existing session_state values untouched


def _sync_genus_species_to_bender(b: Bender) -> None:
    """Store identity, notebook-style specimen metadata, and protocol IDs for HDF5 export."""
    meta = dict(getattr(b, 'h5_protocol_metadata', {}) or {})
    meta['genus_species'] = str(st.session_state.get('gui_genus_species') or '').strip()
    sid = str(st.session_state.get('gui_specimen_id') or '').strip()
    meta['specimen_id'] = sid
    setattr(b, 'fishcode', sid)

    def _str_attr(sess_key: str, bender_attr: str) -> None:
        if sess_key not in st.session_state:
            return
        s = str(st.session_state[sess_key] or '').strip()
        meta[bender_attr] = s
        setattr(b, bender_attr, s)

    def _float_attr(sess_key: str, bender_attr: str) -> None:
        if sess_key not in st.session_state:
            return
        try:
            v = float(st.session_state[sess_key])
        except (TypeError, ValueError):
            return
        meta[bender_attr] = v
        setattr(b, bender_attr, v)

    _str_attr('bio_segment', 'segment')
    _float_attr('bio_fishmass', 'fishmass')
    _float_attr('bio_fishlen_TL', 'fishlen_TL')
    _float_attr('bio_fishlen_SL', 'fishlen_SL')
    _float_attr('bio_xsec_height', 'xsec_height')
    if 'bio_prep_condition' in st.session_state:
        meta['prep_condition'] = str(st.session_state.get('bio_prep_condition') or '').strip()
    if 'bio_temp_room' in st.session_state:
        try:
            meta['temp_C_room'] = float(st.session_state['bio_temp_room'])
        except (TypeError, ValueError):
            pass
    if 'bio_temp_tank' in st.session_state:
        try:
            meta['temp_C_tank'] = float(st.session_state['bio_temp_tank'])
        except (TypeError, ValueError):
            pass
    if 'bio_dvert' in st.session_state:
        try:
            meta['dvert'] = float(st.session_state['bio_dvert'])
        except (TypeError, ValueError):
            pass
    if 'bio_dhoriz' in st.session_state:
        try:
            meta['dhoriz'] = float(st.session_state['bio_dhoriz'])
        except (TypeError, ValueError):
            pass

    b.h5_protocol_metadata = meta


def _apply_specimen_identity_to_bender(b: Bender) -> None:
    """Copy genus/species, specimen ID, ``fishcode``, and ``segment`` from section 3 onto ``b`` and ``h5_protocol_metadata``."""
    meta = dict(getattr(b, 'h5_protocol_metadata', {}) or {})
    meta['genus_species'] = str(st.session_state.get('gui_genus_species') or '').strip()
    sid = str(st.session_state.get('gui_specimen_id') or '').strip()
    meta['specimen_id'] = sid
    b.fishcode = sid
    seg = str(st.session_state.get('bio_segment') or '').strip()
    meta['segment'] = seg
    b.segment = seg
    b.h5_protocol_metadata = meta
    _mark_bio_applied()


def _apply_pair(b: Bender, name: str, value):
    if value is None and name in ('random_seed', 'block_sequence'):
        setattr(b, name, None)
        return
    if value is None:
        return
    if name == 'lateral_mode' and isinstance(value, str) and not value.strip():
        setattr(b, 'lateral_mode', None)
        return
    if name == 'n_end_cycles':
        try:
            value = max(0, int(value))
        except (TypeError, ValueError):
            value = 0
    if name == 'cycles_per_step':
        try:
            value = max(1, int(value))
        except (TypeError, ValueError):
            value = 1
    setattr(b, name, value)


def _apply_form_updates(b: Bender, updates: dict, tt: str):
    if tt == 'calibration':
        emb = st.session_state.get('gui_calibration_embedded_base')
        if isinstance(emb, dict):
            proc = emb.get('procedure')
            if isinstance(proc, dict):
                for k, v in proc.items():
                    _apply_pair(b, k, v)
    for k, v in updates.items():
        _apply_pair(b, k, v)
    b.test_type = tt


def _apply_procedure_form_to_bender(b: Bender, updates: dict, tt: str) -> None:
    """Sync biometrics flags, copy procedure fields onto ``b``, and mirror any QC note text from session."""
    _sync_biometric_flags_from_session(b)
    _apply_form_updates(b, updates, tt)
    _sync_genus_species_to_bender(b)
    _pn = str(st.session_state.get('gui_post_notes') or '').strip()
    if _pn:
        b.post_trial_notes = _pn
    _mark_procedure_applied()
    st.session_state['gui_protocol_confirmed'] = True


def _consume_pending_protocol_template(valid_test_types: list) -> None:
    """Apply a template queued by **Load template** (runs before procedure widgets render)."""
    path = st.session_state.pop('gui_pending_protocol_template_path', None)
    if not path:
        return
    try:
        data = load_protocol_template(path)
        ok, msg = apply_template_to_session_state(
            st.session_state,
            data,
            widget_key=_widget_key,
            valid_test_types=list(valid_test_types),
        )
        st.session_state['gui_protocol_load_feedback'] = (ok, msg)
        if ok:
            st.session_state['gui_proc_apply_invalidated'] = True
    except OSError as e:
        st.session_state['gui_protocol_load_feedback'] = (False, f'Could not read template: {e}')
    except json.JSONDecodeError as e:
        st.session_state['gui_protocol_load_feedback'] = (False, f'Could not read template file: {e}')
    except Exception as e:
        st.session_state['gui_protocol_load_feedback'] = (False, f'{type(e).__name__}: {e}')


def _consume_pending_biometrics_template() -> None:
    path = st.session_state.pop('gui_pending_biometrics_path', None)
    if not path:
        return
    try:
        data = load_biometrics_template(path)
        ok, msg = apply_biometrics_template_to_session(st.session_state, data)
        st.session_state['gui_biometrics_load_feedback'] = (ok, msg)
        if ok:
            st.session_state['gui_bio_apply_invalidated'] = True
            st.session_state.pop('gui_tpl_bio_done', None)
            # #region agent log
            _agent_debug_log(
                hypothesis_id='A',
                location='bender_streamlit_gui.py:_consume_pending_biometrics_template',
                message='template_loaded',
                data={
                    'path': os.path.basename(str(path)),
                    'fishmass': st.session_state.get('bio_fishmass'),
                    'dclamp': st.session_state.get('bio_dclamp'),
                    'xsec': st.session_state.get('bio_xsec'),
                    'invalidated': True,
                },
            )
            # #endregion
            st.rerun()
    except OSError as e:
        st.session_state['gui_biometrics_load_feedback'] = (False, f'Could not read file: {e}')
    except json.JSONDecodeError as e:
        st.session_state['gui_biometrics_load_feedback'] = (False, f'Could not read biometrics file: {e}')
    except Exception as e:
        st.session_state['gui_biometrics_load_feedback'] = (False, f'{type(e).__name__}: {e}')


def _apply_intrinsic_biometrics_to_bender(b: Bender) -> None:
    """Whole-body TL/SL and mass → ``b``; identity/metadata via ``_sync_genus_species_to_bender``."""
    b.fishlen_TL = float(st.session_state['bio_fishlen_TL'])
    b.fishlen_SL = float(st.session_state['bio_fishlen_SL'])
    b.fishmass = float(st.session_state['bio_fishmass'])
    _sync_genus_species_to_bender(b)
    _mark_bio_applied()


def _apply_experimental_conditions_to_bender(b: Bender) -> None:
    """Room/tank temperatures and prep condition → ``b`` and ``h5_protocol_metadata``."""
    b.temp_C_room = float(st.session_state['bio_temp_room'])
    b.temp_C_tank = float(st.session_state['bio_temp_tank'])
    meta = dict(getattr(b, 'h5_protocol_metadata', {}) or {})
    meta['temp_C_room'] = float(st.session_state['bio_temp_room'])
    meta['temp_C_tank'] = float(st.session_state['bio_temp_tank'])
    meta['prep_condition'] = str(st.session_state.get('bio_prep_condition') or '').strip()
    b.h5_protocol_metadata = meta
    _mark_bio_applied()


def _apply_clamp_geometry_to_bender(b: Bender) -> bool:
    """Clamp spacing, bend position, cross-section, and vertical/horizontal offsets → ``b``."""
    xw = float(st.session_state['bio_xsec'])
    md = float(st.session_state.get('bio_muscle_depth', 0.0) or 0.0)
    try:
        _strain_lever_arm_m(xw, md)
    except ValueError as exc:
        st.error(str(exc))
        return False
    b.dclamp = float(st.session_state['bio_dclamp'])
    b.test_segment_length_mm = float(st.session_state['bio_dclamp'])
    b.dbend = float(st.session_state['bio_dbend'])
    b.test_segment_position_mm = float(st.session_state['bio_dbend'])
    b.xsec_width = xw
    b.target_muscle_depth_mm = md
    b.xsec_height = float(st.session_state['bio_xsec_height'])
    b.dvert = float(st.session_state['bio_dvert'])
    b.dhoriz = float(st.session_state['bio_dhoriz'])
    meta = dict(getattr(b, 'h5_protocol_metadata', {}) or {})
    meta['dvert'] = float(st.session_state['bio_dvert'])
    meta['dhoriz'] = float(st.session_state['bio_dhoriz'])
    meta['target_muscle_depth_mm'] = md
    b.h5_protocol_metadata = meta
    _mark_bio_applied()
    return True


def _apply_mounted_profile_inertial_to_bender(b: Bender) -> None:
    """Tapered outline + specimen density + length for profiled / inertial (and frustum-style) corrections → ``b``."""
    rho = _resolved_bio_prof_rho_g_per_mm3()
    _queue_bio_prof_rho_widget_sync_from_preset()
    stations = b.make_profile_stations(
        st.session_state['bio_prof_ph'],
        st.session_state['bio_prof_pw'],
        st.session_state['bio_prof_dh'],
        st.session_state['bio_prof_dw'],
    )
    b.set_profiled_specimen_inertial_model(
        stations,
        st.session_state['bio_prof_L'],
        rho,
        clamp_offset_mm=float(st.session_state['bio_prof_clamp']),
        num_samples=int(st.session_state['bio_prof_samples']),
    )
    _mark_bio_applied()


def _apply_all_biometrics_to_bender(b: Bender) -> None:
    """Copy intrinsic, experimental conditions, clamp, profile (density + outline), and inertial flags from section 3 onto ``b``."""
    _sync_biometric_flags_from_session(b)
    _apply_intrinsic_biometrics_to_bender(b)
    _apply_experimental_conditions_to_bender(b)
    if not _apply_clamp_geometry_to_bender(b):
        return
    _apply_mounted_profile_inertial_to_bender(b)
    st.session_state['gui_measurements_confirmed'] = True


def _sync_biometric_flags_from_session(b: Bender):
    """Copy biometrics panel session state onto ``b`` (flags + specimen geometry aliases)."""
    if 'bio_use_theoretical_inertial' in st.session_state:
        b.use_theoretical_inertial_correction = bool(st.session_state['bio_use_theoretical_inertial'])
    if 'bio_dclamp' in st.session_state:
        v = float(st.session_state['bio_dclamp'])
        b.dclamp = v
        b.test_segment_length_mm = v
    if 'bio_dbend' in st.session_state:
        v = float(st.session_state['bio_dbend'])
        b.dbend = v
        b.test_segment_position_mm = v
    if 'bio_xsec' in st.session_state:
        b.xsec_width = float(st.session_state['bio_xsec'])
    if 'bio_muscle_depth' in st.session_state:
        b.target_muscle_depth_mm = float(st.session_state['bio_muscle_depth'] or 0.0)
    if 'bio_temp_room' in st.session_state:
        b.temp_C_room = float(st.session_state['bio_temp_room'])
    if 'bio_temp_tank' in st.session_state:
        b.temp_C_tank = float(st.session_state['bio_temp_tank'])
    if 'bio_xsec_height' in st.session_state:
        b.xsec_height = float(st.session_state['bio_xsec_height'])
    if 'bio_dvert' in st.session_state:
        b.dvert = float(st.session_state['bio_dvert'])
    if 'bio_dhoriz' in st.session_state:
        b.dhoriz = float(st.session_state['bio_dhoriz'])
    if 'gui_specimen_id' in st.session_state:
        b.fishcode = str(st.session_state.get('gui_specimen_id') or '')
    if 'bio_segment' in st.session_state:
        b.segment = str(st.session_state['bio_segment'] or '')
    if 'bio_fishmass' in st.session_state:
        b.fishmass = float(st.session_state['bio_fishmass'])
    if 'bio_fishlen_TL' in st.session_state:
        b.fishlen_TL = float(st.session_state['bio_fishlen_TL'])
    if 'bio_fishlen_SL' in st.session_state:
        b.fishlen_SL = float(st.session_state['bio_fishlen_SL'])


def _bio_prof_outline_mm_from_bender(b: Bender) -> tuple[float, float, float, float]:
    """Proximal H/W and distal H/W (mm) from ``specimen_profile_stations`` when present; else GUI placeholders."""
    st_list = getattr(b, 'specimen_profile_stations', None)
    if not isinstance(st_list, list) or len(st_list) < 2:
        return (4.0, 6.0, 3.5, 5.0)
    valid: list = []
    for s in st_list:
        if not isinstance(s, dict):
            continue
        try:
            float(s.get('position', 0.0))
            float(s.get('height_mm', 0.0))
            float(s.get('width_mm', 0.0))
        except (TypeError, ValueError):
            continue
        valid.append(s)
    if len(valid) < 2:
        return (4.0, 6.0, 3.5, 5.0)
    by_pos = sorted(valid, key=lambda s: float(s['position']))
    p0, p1 = by_pos[0], by_pos[-1]
    return (
        float(p0['height_mm']),
        float(p0['width_mm']),
        float(p1['height_mm']),
        float(p1['width_mm']),
    )


def _bio_widget_defaults_from_bender(b: Bender) -> dict[str, Any]:
    """Map ``Bender`` + protocol metadata to Streamlit biometrics widget keys (same source as **Apply all**)."""
    dc = getattr(b, 'dclamp', None)
    xw = getattr(b, 'xsec_width', None)
    meta = getattr(b, 'h5_protocol_metadata', {}) or {}
    ph, pw, dh, dw = _bio_prof_outline_mm_from_bender(b)
    _fm = getattr(b, 'fishmass', None)
    _ftl = getattr(b, 'fishlen_TL', None)
    _fsl = getattr(b, 'fishlen_SL', None)
    _xh = getattr(b, 'xsec_height', None)
    return {
        'gui_genus_species': str(meta.get('genus_species', '') or '').strip(),
        'gui_specimen_id': str(meta.get('specimen_id') or getattr(b, 'fishcode', '') or '').strip(),
        'bio_segment': str(getattr(b, 'segment', '') or ''),
        'bio_fishmass': float(_fm) if _fm is not None and math.isfinite(float(_fm)) else 0.0,
        'bio_fishlen_TL': float(_ftl) if _ftl is not None and math.isfinite(float(_ftl)) else 0.0,
        'bio_fishlen_SL': float(_fsl) if _fsl is not None and math.isfinite(float(_fsl)) else 0.0,
        'bio_xsec_height': float(_xh)
        if _xh is not None and math.isfinite(float(_xh))
        else (float(xw) if xw is not None else 8.0),
        'bio_dvert': float(getattr(b, 'dvert', 0.0) or 0.0),
        'bio_dhoriz': float(getattr(b, 'dhoriz', 0.0) or 0.0),
        'bio_dclamp': float(dc) if dc is not None else 10.0,
        'bio_xsec': float(xw) if xw is not None else 8.0,
        'bio_dbend': float(getattr(b, 'dbend', 0.0) or 0.0),
        'bio_temp_room': float(getattr(b, 'temp_C_room', 22.0) or 22.0),
        'bio_temp_tank': float(getattr(b, 'temp_C_tank', 22.0) or 22.0),
        'bio_prep_condition': str(meta.get('prep_condition', '') or ''),
        'bio_use_theoretical_inertial': bool(getattr(b, 'use_theoretical_inertial_correction', False)),
        'bio_prof_L': float(getattr(b, 'specimen_profile_length_mm', 25.0) or 25.0),
        'bio_prof_rho': float(getattr(b, 'specimen_profile_density_g_per_mm3', 1.03e-3) or 1.03e-3),
        'bio_prof_ph': ph,
        'bio_prof_pw': pw,
        'bio_prof_dh': dh,
        'bio_prof_dw': dw,
        'bio_prof_clamp': float(getattr(b, 'specimen_profile_clamp_offset_mm', 20.0) or 20.0),
        'bio_prof_samples': int(getattr(b, 'specimen_profile_num_samples', 120) or 120),
        'bio_prof_rho_preset': BIO_DENSITY_PRESET_LABELS[0],
    }


def _init_biometrics_session_state(b: Bender, *, force: bool = False):
    """Seed Streamlit widget keys from ``b`` (``force`` overwrites after config reload)."""
    defaults = _bio_widget_defaults_from_bender(b)
    _seeded: list[str] = []
    for key, val in defaults.items():
        if force or key not in st.session_state:
            st.session_state[key] = val
            _seeded.append(key)
    if _seeded:
        # #region agent log
        _agent_debug_log(
            hypothesis_id='C',
            location='bender_streamlit_gui.py:_init_biometrics_session_state',
            message='seeded_from_bender',
            data={
                'force': force,
                'seeded_keys': _seeded[:8],
                'seeded_count': len(_seeded),
                'route': _nav_route(),
                'step': _stepwise_step() if _nav_route() == 'stepwise' else None,
                'dclamp_after': st.session_state.get('bio_dclamp'),
            },
        )
        # #endregion


def _rehydrate_missing_biometrics_from_bender(b: Bender) -> None:
    """Restore dropped or stale-empty ``bio_*`` / identity keys from ``B`` when measurements were applied.

    Streamlit may omit widget keys between runs; keys can also stay as ``''``/0 while ``Bender`` still holds
    values from **Apply all biometrics** (``setdefault`` keeps empty keys from being overwritten by ``_init``).
    """
    if not bool(st.session_state.get('gui_measurements_confirmed')):
        return
    defaults = _bio_widget_defaults_from_bender(b)
    for key, val in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = val
            continue
        if key in ('gui_specimen_id', 'gui_genus_species'):
            if not str(st.session_state.get(key) or '').strip():
                if str(val or '').strip():
                    st.session_state[key] = val
            continue
        if key == 'bio_fishmass':
            m = _session_float('bio_fishmass')
            if m is None or m <= 0:
                fm = getattr(b, 'fishmass', None)
                if fm is not None:
                    try:
                        fmf = float(fm)
                    except (TypeError, ValueError):
                        continue
                    if math.isfinite(fmf) and fmf > 0:
                        st.session_state[key] = fmf
            continue
        if key == 'bio_dclamp':
            dc = _session_float('bio_dclamp')
            if dc is None or dc <= 0:
                bd = getattr(b, 'dclamp', None)
                if bd is None:
                    bd = getattr(b, 'test_segment_length_mm', None)
                if bd is not None:
                    try:
                        bdf = float(bd)
                    except (TypeError, ValueError):
                        continue
                    if math.isfinite(bdf) and bdf > 0:
                        st.session_state[key] = bdf
            continue
        if key == 'bio_xsec':
            xw = _session_float('bio_xsec')
            if xw is None or xw <= 0:
                bx = getattr(b, 'xsec_width', None)
                if bx is not None:
                    try:
                        bxf = float(bx)
                    except (TypeError, ValueError):
                        continue
                    if math.isfinite(bxf) and bxf > 0:
                        st.session_state[key] = bxf
            continue


def _nav_route() -> str:
    return str(st.session_state.get('gui_app_route') or 'landing')


def _clear_pending_bio_nav_warning() -> None:
    for _k in (
        'gui_pending_bio_nav_route',
        'gui_pending_bio_nav_stepwise_step',
        'gui_pending_bio_nav_clear_stepwise',
        'gui_pending_bio_nav_origin',
    ):
        st.session_state.pop(_k, None)


def _apply_route_switch(*, target_route: str, stepwise_step: Optional[int], clear_stepwise: bool) -> None:
    st.session_state['gui_app_route'] = str(target_route or 'landing')
    if clear_stepwise:
        st.session_state.pop('gui_stepwise_step', None)
    elif stepwise_step is not None:
        st.session_state['gui_stepwise_step'] = int(stepwise_step)


def _attempt_route_switch_with_bio_warning(
    *,
    target_route: str,
    origin: str,
    stepwise_step: Optional[int] = None,
    clear_stepwise: bool = False,
) -> bool:
    if _bio_apply_dirty():
        st.session_state['gui_pending_bio_nav_route'] = str(target_route or 'landing')
        st.session_state['gui_pending_bio_nav_stepwise_step'] = (
            None if stepwise_step is None else int(stepwise_step)
        )
        st.session_state['gui_pending_bio_nav_clear_stepwise'] = bool(clear_stepwise)
        st.session_state['gui_pending_bio_nav_origin'] = str(origin or '')
        return False
    _apply_route_switch(target_route=target_route, stepwise_step=stepwise_step, clear_stepwise=clear_stepwise)
    _clear_pending_bio_nav_warning()
    return True


def _render_pending_bio_nav_warning(origin: str) -> None:
    if str(st.session_state.get('gui_pending_bio_nav_origin') or '') != str(origin or ''):
        return
    _target = str(st.session_state.get('gui_pending_bio_nav_route') or '').strip()
    if not _target:
        return
    _step = st.session_state.get('gui_pending_bio_nav_stepwise_step', None)
    _clear_step = bool(st.session_state.get('gui_pending_bio_nav_clear_stepwise', False))
    st.warning(
        'You have unsaved biometrics form edits. Click Apply in Measurements before switching workflows, '
        'or continue and keep edits only in this session state.'
    )
    st.caption('- Apply specimen')
    st.caption('- Apply clamp geometry & inertial correction (or Apply all biometrics)')
    _w1, _w2 = st.columns(2, gap='small')
    with _w1:
        if st.button('Switch anyway', key=f'gui_bio_nav_switch_anyway_{origin}', type='primary', use_container_width=True):
            _apply_route_switch(target_route=_target, stepwise_step=_step, clear_stepwise=_clear_step)
            _clear_pending_bio_nav_warning()
            st.rerun()
    with _w2:
        if st.button('Stay here', key=f'gui_bio_nav_stay_{origin}', use_container_width=True):
            _clear_pending_bio_nav_warning()
            st.rerun()


def _stepwise_step() -> int:
    v = st.session_state.get('gui_stepwise_step', 0)
    try:
        return int(v)
    except (TypeError, ValueError):
        return 0


def _stepwise_on_data_file_path_step() -> bool:
    # Single-page workflow no longer uses stepwise tab position to gate setup behavior.
    return False


def _template_hide_config_build_new() -> bool:
    """In template mode, user said they have a saved config — skip the *Build new* path in section 1."""
    return _nav_route() == 'templates' and bool(st.session_state.get('gui_tpl_chk_config', True))


def _tpl_only_procedure() -> bool:
    """True when user has config + biometrics files but still needs to set the procedure in the app (no protocol JSON yet)."""
    if _nav_route() != 'templates':
        return False
    if not bool(st.session_state.get('gui_tpl_chk_config', True)):
        return False
    if not bool(st.session_state.get('gui_tpl_chk_biometrics', True)):
        return False
    # Third box: "I already have a protocol template" — when checked, use full sections (load protocol from file).
    return not bool(st.session_state.get('gui_tpl_have_protocol_template', False))


def _show_hw_config_section() -> bool:
    """Section 1 · hardware configuration (no data-folder block)."""
    if _nav_route() == 'templates' and _tpl_only_procedure():
        return False
    return True


def _show_data_path_section() -> bool:
    """Section 2 · data folder & file name (+ **Load hardware configuration and data path** for load-existing)."""
    if _nav_route() == 'templates' and _tpl_only_procedure():
        return False
    return True


def _show_full_sec2() -> bool:
    """Section 3 · biometrics."""
    if _nav_route() == 'templates' and _tpl_only_procedure():
        return False
    return True


def _show_sec3_through_6() -> bool:
    return True


def _show_sec7_and_8() -> bool:
    return True


_LANDING_STIM_HZ = 75.0


def _landing_stimulus_75hz_gated(t: np.ndarray, windows: tuple[tuple[float, float], ...]) -> np.ndarray:
    """50%-duty square wave at ``_LANDING_STIM_HZ``, gated on during ``windows`` (schematic bursts)."""
    carrier = (np.sin(2 * np.pi * _LANDING_STIM_HZ * t) >= 0).astype(np.float64)
    env = np.zeros_like(t)
    for t0, t1 in windows:
        env[(t >= t0) & (t <= t1)] = 1.0
    return env * carrier


def _landing_plotly_demo_layout(fig: go.Figure, *, x_title: str, y_title: str, chart_title: str) -> None:
    fig.update_layout(
        margin=dict(l=48, r=12, t=36, b=40),
        height=240,
        title=dict(text=chart_title, font=dict(size=13, color='#1e293b')),
        xaxis_title=x_title,
        yaxis_title=y_title,
        paper_bgcolor='#ffffff',
        plot_bgcolor='#f8fafc',
        font=dict(size=11, color='#1e293b'),
        xaxis=dict(showgrid=True, gridcolor='#e2e8f0', zeroline=False),
        yaxis=dict(showgrid=True, gridcolor='#e2e8f0', zeroline=True, zerolinecolor='#cbd5e1'),
        showlegend=False,
    )


def _landing_demo_figure(test_type: str) -> go.Figure:
    """Schematic Plotly figures for landing-page education (not real trial data)."""
    tt = str(test_type)
    if tt == 'dynamic':
        t = np.linspace(0, 2.5, 500)
        y = 0.85 * np.sin(2 * np.pi * 2.0 * t)
        fig = go.Figure(go.Scatter(x=t, y=y, mode='lines', line=dict(color='#0ea5e9', width=2)))
        _landing_plotly_demo_layout(
            fig, x_title='Time (s)', y_title='Bending / strain (a.u.)', chart_title='Steady cyclic command'
        )
        return fig
    if tt == 'frequency_sweep':
        t = np.linspace(0, 3, 800)
        dt = float(t[1] - t[0])
        inst_f = 0.5 + 1.5 * (t / t[-1]) ** 1.1
        phase = np.cumsum(2 * np.pi * inst_f * dt)
        y = 0.55 * np.sin(phase)
        fig = go.Figure(go.Scatter(x=t, y=y, mode='lines', line=dict(color='#6366f1', width=2)))
        _landing_plotly_demo_layout(
            fig, x_title='Time (s)', y_title='Bending / strain (a.u.)', chart_title='Frequency increases over time'
        )
        return fig
    if tt == 'frequency_step':
        t = np.linspace(0, 2.4, 600)
        seg = (t * 3 / 2.4).astype(int)
        seg = np.clip(seg, 0, 2)
        freqs = np.array([0.8, 1.6, 2.4])[seg]
        dt = float(t[1] - t[0])
        phase = np.cumsum(2 * np.pi * freqs * dt)
        y = 0.6 * np.sin(phase)
        fig = go.Figure(go.Scatter(x=t, y=y, mode='lines', line=dict(color='#8b5cf6', width=2)))
        _landing_plotly_demo_layout(
            fig, x_title='Time (s)', y_title='Bending / strain (a.u.)', chart_title='Plateau segments at fixed frequencies'
        )
        return fig
    if tt == 'curvature_step':
        t = np.linspace(0, 2.2, 500)
        k = (t * 2.2 / 2.2).astype(int)
        k = np.clip(k, 0, 3)
        levels = np.array([-0.5, 0.2, 0.75, 1.1])[k]
        fig = go.Figure(go.Scatter(x=t, y=levels, mode='lines', line=dict(color='#0d9488', width=2.5)))
        _landing_plotly_demo_layout(
            fig, x_title='Time (s)', y_title='Target curvature (a.u.)', chart_title='Hold at each curvature level'
        )
        return fig
    if tt == 'step_change':
        t = np.linspace(0, 2.0, 400)
        amp = np.where(t < 1.0, 0.35, 0.9)
        y = amp * np.sin(2 * np.pi * 1.5 * t)
        fig = go.Figure(go.Scatter(x=t, y=y, mode='lines', line=dict(color='#ea580c', width=2)))
        _landing_plotly_demo_layout(
            fig, x_title='Time (s)', y_title='Bending / strain (a.u.)', chart_title='Abrupt amplitude / condition change'
        )
        return fig
    if tt == 'isometric':
        t = np.linspace(0, 2.6, int(round(2.6 * 1200)) + 1)
        stim = _landing_stimulus_75hz_gated(t, ((0.35, 0.52), (1.05, 1.22), (1.85, 2.02)))
        strain = np.ones_like(t) * 0.55
        torque = 0.12 + 0.62 * (1 - np.exp(-2.0 * t)) + 0.07 * stim
        fig = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.07)
        fig.add_trace(
            go.Scatter(x=t, y=stim, mode='lines', name='Stimulus (75 Hz)', line=dict(color='#a855f7', width=1.2)),
            row=1,
            col=1,
        )
        fig.add_trace(
            go.Scatter(x=t, y=strain, mode='lines', name='Strain (held)', line=dict(color='#0ea5e9', width=2)),
            row=2,
            col=1,
        )
        fig.add_trace(
            go.Scatter(x=t, y=torque, mode='lines', name='Torque', line=dict(color='#be123c', width=2)),
            row=3,
            col=1,
        )
        fig.update_yaxes(title_text='Stim (75 Hz, a.u.)', row=1, col=1, showgrid=True, gridcolor='#e2e8f0')
        fig.update_yaxes(title_text='Strain (a.u.)', row=2, col=1, showgrid=True, gridcolor='#e2e8f0')
        fig.update_yaxes(title_text='Torque (a.u.)', row=3, col=1, showgrid=True, gridcolor='#e2e8f0')
        fig.update_xaxes(title_text='Time (s)', row=3, col=1, showgrid=True, gridcolor='#e2e8f0')
        fig.update_layout(
            margin=dict(l=52, r=18, t=40, b=36),
            height=320,
            title=dict(
                text='Isometric: 75 Hz stimulus bursts, held strain, evolving torque',
                font=dict(size=13, color='#1e293b'),
            ),
            paper_bgcolor='#ffffff',
            plot_bgcolor='#f8fafc',
            font=dict(size=11, color='#1e293b'),
            showlegend=True,
            legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='left', x=0),
        )
        return fig
    if tt == 'isovelocity':
        t = np.linspace(0, 2.2, int(round(2.2 * 1200)) + 1)
        m1, m2, m3 = 0.42, 0.28, -0.35
        m1a = t < 0.75
        m2a = (t >= 0.75) & (t < 1.45)
        m3a = t >= 1.45
        strain = np.zeros_like(t)
        strain[m1a] = m1 * t[m1a]
        strain[m2a] = m1 * 0.75 + m2 * (t[m2a] - 0.75)
        strain[m3a] = m1 * 0.75 + m2 * 0.7 + m3 * (t[m3a] - 1.45)
        vel = np.gradient(strain, t, edge_order=2)
        stim = _landing_stimulus_75hz_gated(t, ((0.2, 0.38), (0.9, 1.05), (1.55, 1.72)))
        fig = make_subplots(rows=3, cols=1, shared_xaxes=True, vertical_spacing=0.07)
        fig.add_trace(
            go.Scatter(x=t, y=stim, mode='lines', name='Stimulus (75 Hz)', line=dict(color='#a855f7', width=1.2)),
            row=1,
            col=1,
        )
        fig.add_trace(
            go.Scatter(x=t, y=vel, mode='lines', name='Strain rate', line=dict(color='#d97706', width=2)),
            row=2,
            col=1,
        )
        fig.add_trace(
            go.Scatter(x=t, y=strain, mode='lines', name='Strain', line=dict(color='#15803d', width=2)),
            row=3,
            col=1,
        )
        fig.update_yaxes(title_text='Stim (75 Hz, a.u.)', row=1, col=1, showgrid=True, gridcolor='#e2e8f0')
        fig.update_yaxes(title_text='dε/dt (a.u.)', row=2, col=1, showgrid=True, gridcolor='#e2e8f0')
        fig.update_yaxes(title_text='Strain (a.u.)', row=3, col=1, showgrid=True, gridcolor='#e2e8f0')
        fig.update_xaxes(title_text='Time (s)', row=3, col=1, showgrid=True, gridcolor='#e2e8f0')
        fig.update_layout(
            margin=dict(l=52, r=18, t=40, b=36),
            height=320,
            title=dict(
                text='Isovelocity: 75 Hz stimulus bursts, strain rate, integrated strain',
                font=dict(size=13, color='#1e293b'),
            ),
            paper_bgcolor='#ffffff',
            plot_bgcolor='#f8fafc',
            font=dict(size=11, color='#1e293b'),
            showlegend=True,
            legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='left', x=0),
        )
        return fig
    # calibration — commanded motion vs raw vs corrected torque (schematic)
    t = np.linspace(0, 2.2, 550)
    theta_cmd = 0.32 * np.sin(2 * np.pi * 1.1 * t)
    theta_meas = 0.31 * np.sin(2 * np.pi * 1.1 * t - 0.07) + 0.018 * np.sin(2 * np.pi * 16.0 * t)
    tau_raw = 0.55 * np.sin(2 * np.pi * 1.1 * t) + 0.22 * np.sin(2 * np.pi * 16.0 * t + 0.4)
    tau_model = 0.20 * np.sin(2 * np.pi * 16.0 * t + 0.4)
    tau_corr = tau_raw - tau_model
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.11)
    fig.add_trace(
        go.Scatter(x=t, y=theta_cmd, mode='lines', name='θ command', line=dict(color='#2563eb', width=2)),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(x=t, y=theta_meas, mode='lines', name='θ measured', line=dict(color='#f97316', width=2, dash='dot')),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Scatter(x=t, y=tau_raw, mode='lines', name='Torque raw', line=dict(color='#94a3b8', width=2)),
        row=2,
        col=1,
    )
    fig.add_trace(
        go.Scatter(x=t, y=tau_corr, mode='lines', name='Torque corrected', line=dict(color='#0d9488', width=2)),
        row=2,
        col=1,
    )
    fig.update_yaxes(title_text='Angle (a.u.)', row=1, col=1, showgrid=True, gridcolor='#e2e8f0')
    fig.update_yaxes(title_text='Torque (a.u.)', row=2, col=1, showgrid=True, gridcolor='#e2e8f0')
    fig.update_xaxes(title_text='Time (s)', row=2, col=1, showgrid=True, gridcolor='#e2e8f0')
    fig.update_layout(
        margin=dict(l=52, r=18, t=44, b=40),
        height=300,
        title=dict(
            text='Calibration: repeat a known motion; subtract modeled inertial ripple from torque',
            font=dict(size=12, color='#1e293b'),
        ),
        paper_bgcolor='#ffffff',
        plot_bgcolor='#f8fafc',
        font=dict(size=11, color='#1e293b'),
        showlegend=True,
        legend=dict(orientation='h', yanchor='bottom', y=1.05, xanchor='left', x=0),
    )
    return fig


_LANDING_EXPERIMENT_BLURBS: dict[str, tuple[str, str]] = {
    'dynamic': (
        'Dynamic (cyclic)',
        'Repeat the same bend cycle at fixed frequency and amplitude to measure cyclic mechanics, energy dissipation, '
        'and responses comparable across specimens — useful for stiffness–damping characterization and fatigue-style loading.',
    ),
    'frequency_sweep': (
        'Frequency sweep',
        'Vary excitation frequency over time to map how mechanical response changes with frequency — helpful for '
        'resonance-like behavior and frequency-dependent stiffness in soft tissues.',
    ),
    'frequency_step': (
        'Frequency step',
        'Hold several discrete frequencies in separate segments so each band is sampled clearly — separates steady '
        'response at each frequency without the ambiguity of a continuous sweep.',
    ),
    'curvature_step': (
        'Curvature step',
        'Step through held curvature levels to read torque (or force) at specific bending states — supports '
        'quasi-static–like curves and stress-relaxation–style interpretation between plateaus.',
    ),
    'step_change': (
        'Step change',
        'Suddenly change amplitude, frequency, or related parameters mid-protocol to probe transient adjustment, '
        'history dependence, and nonlinear responses to a new loading condition.',
    ),
    'isometric': (
        'Isometric hold',
        'Geometry is held constant while you may deliver **electrical or other stimuli** in pulses; torque evolves with '
        'activation and passive stress relaxation — analogous to isometric contraction, useful for separating active and '
        'passive mechanics without ongoing bending.',
    ),
    'isovelocity': (
        'Isovelocity ramp',
        'Bending follows a **commanded strain rate** (piecewise ramps here); optional **stimulus trains** can be aligned '
        'with loading segments so velocity-dependent effects are not confounded with steady cyclic loops.',
    ),
    'calibration': (
        'Calibration',
        'Run a **repeatable reference motion** (often the same `test_type` you use later). **Raw torque** mixes tissue load '
        'with high-frequency inertial “ripple” from acceleration; after fitting a dynamic model, **corrected torque** '
        'tracks the biology-linked component more cleanly for subsequent trials.',
    ),
}


def _render_landing_equation_of_motion_box() -> None:
    th = '\N{GREEK SMALL LETTER THETA}'
    thdd = th + '\N{COMBINING DIAERESIS}'
    thd = th + '\N{COMBINING DOT ABOVE}'
    st.markdown(
        f'<div class="bnd-landing-eom-box">'
        f'<p class="bnd-eom-heading"><strong>Fundamental equation of motion (schematic)</strong></p>'
        f'<p class="bnd-eom-eq"><em>J</em>{thdd} + <em>B</em>{thd} + <em>K</em>{th} ≈ τ<sub>motor</sub> − '
        f'τ<sub>passive</sub>({th})</p>'
        f'<p class="bnd-eom-note">Inertial (<em>J{thdd}</em>), viscous (<em>B{thd}</em>), and elastic (<em>K{th}</em>) '
        f'terms balance the actuator against passive tissue. Calibration and correction in the workflow help separate '
        f'biology-linked torque from acceleration-dependent artifacts in recorded data.</p></div>',
        unsafe_allow_html=True,
    )


def _render_landing_learn_section() -> None:
    _render_landing_equation_of_motion_box()
    st.markdown(
        '<p class="bnd-landing-learn-intro">'
        'Each run uses one <strong>experiment type</strong> (<code>test_type</code>); configure hardware and biometrics, '
        'preview commanded motion, then run and save. Curves below are <strong>schematics</strong> — not real specimen '
        'data — to show the idea behind each type.</p>',
        unsafe_allow_html=True,
    )
    tabs = st.tabs(
        [
            'Oscillatory',
            'Steps & sweeps',
            'Quasi-static',
            'Calibration',
        ]
    )
    with tabs[0]:
        c1, c2 = st.columns(2)
        with c1:
            title, blurb = _LANDING_EXPERIMENT_BLURBS['dynamic']
            st.markdown(f'**{title}**')
            st.caption(blurb)
            st.plotly_chart(_landing_demo_figure('dynamic'), use_container_width=True, config={'displayModeBar': False})
        with c2:
            title, blurb = _LANDING_EXPERIMENT_BLURBS['frequency_sweep']
            st.markdown(f'**{title}**')
            st.caption(blurb)
            st.plotly_chart(
                _landing_demo_figure('frequency_sweep'), use_container_width=True, config={'displayModeBar': False}
            )
    with tabs[1]:
        c1, c2 = st.columns(2)
        with c1:
            title, blurb = _LANDING_EXPERIMENT_BLURBS['frequency_step']
            st.markdown(f'**{title}**')
            st.caption(blurb)
            st.plotly_chart(
                _landing_demo_figure('frequency_step'), use_container_width=True, config={'displayModeBar': False}
            )
        with c2:
            title, blurb = _LANDING_EXPERIMENT_BLURBS['curvature_step']
            st.markdown(f'**{title}**')
            st.caption(blurb)
            st.plotly_chart(
                _landing_demo_figure('curvature_step'), use_container_width=True, config={'displayModeBar': False}
            )
        c3, c4 = st.columns(2)
        with c3:
            title, blurb = _LANDING_EXPERIMENT_BLURBS['step_change']
            st.markdown(f'**{title}**')
            st.caption(blurb)
            st.plotly_chart(
                _landing_demo_figure('step_change'), use_container_width=True, config={'displayModeBar': False}
            )
        with c4:
            st.markdown('**How to pick a type**')
            st.caption(
                'In the full workflow, choose **Experiment type (test_type)** first; the form shows only the motion '
                'fields that protocol needs. Preview plots approximate commanded motion before you run on hardware.'
            )
    with tabs[2]:
        c1, c2 = st.columns(2)
        with c1:
            title, blurb = _LANDING_EXPERIMENT_BLURBS['isometric']
            st.markdown(f'**{title}**')
            st.caption(blurb)
            st.plotly_chart(_landing_demo_figure('isometric'), use_container_width=True, config={'displayModeBar': False})
        with c2:
            title, blurb = _LANDING_EXPERIMENT_BLURBS['isovelocity']
            st.markdown(f'**{title}**')
            st.caption(blurb)
            st.plotly_chart(
                _landing_demo_figure('isovelocity'), use_container_width=True, config={'displayModeBar': False}
            )
    with tabs[3]:
        title, blurb = _LANDING_EXPERIMENT_BLURBS['calibration']
        st.markdown(f'**{title}**')
        st.caption(blurb)
        st.plotly_chart(_landing_demo_figure('calibration'), use_container_width=True, config={'displayModeBar': False})


def _render_landing_page() -> None:
    st.markdown(
        '<div id="bnd-main-content" tabindex="-1"></div>',
        unsafe_allow_html=True,
    )
    st.markdown(
        '<div class="bnd-landing-page" aria-hidden="true"></div>',
        unsafe_allow_html=True,
    )
    # Hero/logo layout must load before the <img> or the PNG paints at intrinsic size and then jumps.
    st.markdown(
        """
<style>
:is(#root, body:has(.bnd-landing-page)) .bnd-landing-hero-wrap {
    position: relative;
    width: 100%;
    box-sizing: border-box;
    padding-right: min(38%, 15rem);
    min-height: 11rem;
    margin-bottom: 0.35rem;
    padding-bottom: 0.35rem;
}
:is(#root, body:has(.bnd-landing-page)) .bnd-landing-hero-text { display: block; }
:is(#root, body:has(.bnd-landing-page)) .bnd-landing-hero-logo {
    position: absolute;
    top: 0;
    right: 0;
    bottom: 0;
    width: min(34%, 13.5rem);
    max-width: 13.5rem;
    display: flex;
    align-items: center;
    justify-content: center;
    pointer-events: none;
    background: #ffffff;
    border-radius: 12px;
}
:is(#root, body:has(.bnd-landing-page)) .bnd-landing-hero-logo-img {
    display: block;
    max-width: 100%;
    max-height: 100%;
    width: auto;
    height: auto;
    object-fit: contain;
    pointer-events: auto;
}
</style>
""",
        unsafe_allow_html=True,
    )
    if os.path.isfile(_LOGO_PATH):
        _logo_uri = _img_data_uri(_LOGO_PATH)
        st.markdown(
            f'<div class="bnd-landing-hero-wrap">'
            f'<div class="bnd-landing-hero-text">'
            f'<h1>The CritterGripper App</h1>'
            f'<p class="bnd-landing-tagline">Non-destructive bending biomechanics</p>'
            f'</div>'
            f'<div class="bnd-landing-hero-logo">'
            f'<img class="bnd-landing-hero-logo-img" src="{_logo_uri}" alt="Laboratory logo" '
            f'style="max-width:min(100%,13.5rem);max-height:100%;width:auto;height:auto;object-fit:contain;display:block;" '
            f'decoding="async" fetchpriority="high"/>'
            f'</div>'
            f'</div>',
            unsafe_allow_html=True,
        )
    else:
        st.title('The CritterGripper App')
        st.markdown(
            '<p class="bnd-landing-tagline">Non-destructive bending biomechanics</p>',
            unsafe_allow_html=True,
        )

    st.markdown(
        '<div class="bnd-landing-hero-rule" aria-hidden="true"></div>',
        unsafe_allow_html=True,
    )

    st.markdown('<div class="bnd-landing-between-sections"></div>', unsafe_allow_html=True)

    with st.container(border=True):
        st.markdown(
            '<h3 class="bnd-landing-section-title">Run experiment</h3>',
            unsafe_allow_html=True,
        )
        st.markdown(
            '<p class="bnd-landing-section-sub">The single-page experiment engine: hardware, data path, specimen, clamp '
            'geometry &amp; inertial correction, protocol, preview, run, save, and notes — all on one scrolling page.</p>',
            unsafe_allow_html=True,
        )
        st.markdown('<div class="bnd-land-workflow-btn-row-marker" aria-hidden="true"></div>', unsafe_allow_html=True)
        if st.button('Run experiment', key='land_scratch', use_container_width=True, type='primary'):
            _autosave_tick(force=True)
            if _attempt_route_switch_with_bio_warning(
                target_route='scratch', origin='landing', clear_stepwise=True
            ):
                st.rerun()
        _render_pending_bio_nav_warning('landing')

        st.markdown(
            '<p class="bnd-landing-sim-cta-sub">Offline exploration (no NI hardware)</p>',
            unsafe_allow_html=True,
        )
        st.markdown('<div class="bnd-land-sim-btn-marker" aria-hidden="true"></div>', unsafe_allow_html=True)
        if st.button(
            'Simulate Bending Mechanics',
            key='land_sim_compare',
            use_container_width=True,
            type='secondary',
            help='Compare two cantilever specimens (material + geometry) using numpy only — no NI-DAQ.',
        ):
            _autosave_tick(force=True)
            st.session_state['gui_app_route'] = 'sim_compare'
            st.session_state.pop('gui_stepwise_step', None)
            st.rerun()

    st.markdown('<div class="bnd-landing-between-sections"></div>', unsafe_allow_html=True)

    with st.container(border=True):
        st.markdown(
            '<h3 class="bnd-landing-section-title">Visualize Data</h3>',
            unsafe_allow_html=True,
        )
        st.markdown(
            '<p class="bnd-landing-section-sub">Open exported <code>.h5</code> files: <strong>Standard</strong> mode '
            'lists CritterGripper series automatically; <strong>Custom</strong> mode uses dropdowns to navigate '
            'the file tree for non-standard layouts. You can edit HDF5 attributes when metadata is wrong.</p>',
            unsafe_allow_html=True,
        )
        if st.button('Open Visualize Data', key='land_h5_explorer', use_container_width=True, type='primary'):
            _autosave_tick(force=True)
            st.session_state['gui_app_route'] = 'h5_explorer'
            st.rerun()

    st.markdown('<div class="bnd-landing-between-sections"></div>', unsafe_allow_html=True)

    with st.container(border=True):
        st.markdown(
            '<h3 class="bnd-landing-section-title">What this App Does</h3>',
            unsafe_allow_html=True,
        )
        st.markdown(
            '<p class="bnd-landing-section-sub bnd-landing-what-section-sub">Physics overview, experiment types, '
            'schematic plots, and why you might run each one.</p>',
            unsafe_allow_html=True,
        )
        _render_landing_learn_section()

    st.markdown('<div class="bnd-landing-between-sections"></div>', unsafe_allow_html=True)
    with st.expander('Settings', expanded=False):
        st.caption('The sidebar is hidden on this page. Adjust theme and text size here; they apply after you open a workflow.')
        _render_display_preferences_sidebar()

    # Landing-only overrides (this function runs only on the home route). Selectors use
    # :is(#root, body:has(.bnd-landing-page)) because some Streamlit builds have no #root id.
    # st.html keeps <style> in the live document (markdown path can drop or scope styles).
    _landing_style = """
<style>
:is(#root, body:has(.bnd-landing-page)) .stApp,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stAppViewContainer"],
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"],
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] > div,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] .block-container {
    background: #ffffff !important;
    background-color: #ffffff !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] .block-container {
    padding-top: 0.35rem !important;
    padding-bottom: 2rem !important;
    max-width: 52rem !important;
    margin-left: auto !important;
    margin-right: auto !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] [data-testid="stVerticalBlock"] > div:first-child {
    margin-top: 0 !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] h1 {
    font-weight: 700 !important;
    letter-spacing: -0.03em !important;
    color: #334155 !important;
    line-height: 1.12 !important;
    margin-top: 0 !important;
    margin-bottom: 0.35rem !important;
    font-size: clamp(2.45rem, 5.2vw, 3.45rem) !important;
}
:is(#root, body:has(.bnd-landing-page)) .bnd-landing-tagline {
    font-size: 1.32rem !important;
    line-height: 1.48 !important;
    color: #64748b !important;
    margin: 0.1rem 0 0.2rem 0 !important;
    font-weight: 450;
}
:is(#root, body:has(.bnd-landing-page)) .bnd-landing-hero-rule {
    height: 0;
    margin: 0.85rem 0 0.55rem 0;
    clear: both;
    padding: 0;
    border: none;
    border-top: 4px solid #c2410c;
    background: none;
}
:is(#root, body:has(.bnd-landing-page)) .bnd-landing-eom-box {
    background: linear-gradient(145deg, #e2f4f0 0%, #dff4f8 50%, #e4e8f5 100%);
    border-radius: 14px;
    padding: 1rem 1.2rem 1.05rem;
    margin: 1.95rem 0 1.05rem 0;
    border-left: 4px solid #0d9488;
}
:is(#root, body:has(.bnd-landing-page)) .bnd-eom-heading { margin: 0 0 0.45rem 0; font-size: 0.98rem; color: #475569; }
:is(#root, body:has(.bnd-landing-page)) .bnd-eom-eq {
    margin: 0 0 0.55rem 0;
    font-family: Georgia, 'Times New Roman', serif;
    font-size: 1.22rem;
    line-height: 1.45;
    color: #334155;
}
:is(#root, body:has(.bnd-landing-page)) .bnd-eom-note { margin: 0; font-size: 0.92rem; line-height: 1.5; color: #64748b; }
:is(#root, body:has(.bnd-landing-page)) .bnd-landing-section-sub.bnd-landing-what-section-sub { margin-bottom: 1.15rem !important; }
:is(#root, body:has(.bnd-landing-page)) .bnd-landing-section-title {
    text-align: center !important;
    width: 100%;
    margin: 0.1rem 0 0.35rem 0 !important;
    font-size: 1.45rem !important;
    font-weight: 600 !important;
    color: #475569 !important;
    letter-spacing: -0.02em;
}
:is(#root, body:has(.bnd-landing-page)) .bnd-landing-section-sub {
    text-align: center;
    color: #64748b;
    font-size: 1.08rem;
    line-height: 1.5;
    margin: 0 0 0.9rem 0;
}
:is(#root, body:has(.bnd-landing-page)) .bnd-landing-learn-intro {
    font-size: 1.02rem;
    line-height: 1.55;
    color: #64748b !important;
    margin: 0 0 0.65rem 0;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] [data-testid="stCaption"] { color: #64748b !important; }
:is(#root, body:has(.bnd-landing-page)) .bnd-landing-between-sections { height: 1.65rem; }
@media (max-width: 980px) {
  :is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] .bnd-land-workflow-cards-marker ~ div [data-testid="column"],
  :is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] .bnd-land-workflow-btn-row-marker ~ div [data-testid="column"] {
    min-width: 100% !important;
    flex: 1 1 100% !important;
  }
  :is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] .bnd-land-workflow-btn-row-marker ~ div [data-testid="stButton"] button {
    min-height: 2.85rem !important;
  }
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stVerticalBlockBorderWrapper"] {
    padding: 1.05rem 1.05rem 1.2rem !important;
    margin-top: 0 !important;
    margin-bottom: 0 !important;
    background: #ffffff !important;
    border: 1px solid #cbd5e1 !important;
    border-radius: 14px !important;
    box-shadow: 0 1px 3px rgba(15, 23, 42, 0.06) !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stVerticalBlockBorderWrapper"] [data-testid="stMarkdownContainer"] p,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stVerticalBlockBorderWrapper"] [data-testid="stMarkdownContainer"] li {
    color: #475569 !important;
}
/* st.button — primary red / secondary white; hover brightens (see workflow theme). */
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-primary"],
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="baseButton-primary"],
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-primary"],
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-primary"],
:is(#root, body:has(.bnd-landing-page)) .stApp div.row-widget.stButton[data-testid="stButton"] button[data-testid="baseButton-primary"],
:is(#root, body:has(.bnd-landing-page)) .stApp div.row-widget.stButton[data-testid="stButton"] button[data-testid="stBaseButton-primary"],
:is(#root, body:has(.bnd-landing-page)) div.row-widget.stButton[data-testid="stButton"] button[data-testid="baseButton-primary"],
:is(#root, body:has(.bnd-landing-page)) div.row-widget.stButton[data-testid="stButton"] button[data-testid="stBaseButton-primary"],
:is(#root, body:has(.bnd-landing-page)) div[data-testid="stButton"] button[data-testid="baseButton-primary"],
:is(#root, body:has(.bnd-landing-page)) div[data-testid="stButton"] button[data-testid="stBaseButton-primary"] {
    min-height: 3rem !important;
    background-color: #dc2626 !important;
    background-image: none !important;
    color: #ffffff !important;
    border: 1px solid #b91c1c !important;
    font-weight: 600 !important;
    transition: background-color 0.15s ease, border-color 0.15s ease, box-shadow 0.15s ease !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-primary"] p,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-primary"] span,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-primary"] p,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-primary"] span,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="baseButton-primary"] p,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="baseButton-primary"] span,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-primary"] p,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-primary"] span {
    color: #ffffff !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-primary"]:hover,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-primary"]:hover,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="baseButton-primary"]:hover,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-primary"]:hover {
    background-color: #f87171 !important;
    border-color: #fca5a5 !important;
    box-shadow: 0 0 0 3px rgba(248, 113, 113, 0.45) !important;
    color: #ffffff !important;
    background-image: none !important;
    filter: none !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-secondary"],
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="baseButton-secondary"],
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-secondary"],
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-secondary"],
:is(#root, body:has(.bnd-landing-page)) .stApp div.row-widget.stButton[data-testid="stButton"] button[data-testid="baseButton-secondary"],
:is(#root, body:has(.bnd-landing-page)) .stApp div.row-widget.stButton[data-testid="stButton"] button[data-testid="stBaseButton-secondary"],
:is(#root, body:has(.bnd-landing-page)) div.row-widget.stButton[data-testid="stButton"] button[data-testid="baseButton-secondary"],
:is(#root, body:has(.bnd-landing-page)) div.row-widget.stButton[data-testid="stButton"] button[data-testid="stBaseButton-secondary"],
:is(#root, body:has(.bnd-landing-page)) div[data-testid="stButton"] button[data-testid="baseButton-secondary"],
:is(#root, body:has(.bnd-landing-page)) div[data-testid="stButton"] button[data-testid="stBaseButton-secondary"] {
    min-height: 3rem !important;
    background-color: #ffffff !important;
    background-image: none !important;
    color: #334155 !important;
    border: 1px solid #cbd5e1 !important;
    font-weight: 600 !important;
    transition: background-color 0.15s ease, border-color 0.15s ease, box-shadow 0.15s ease, color 0.15s ease !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-secondary"] p,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-secondary"] span,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-secondary"] p,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-secondary"] span,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="baseButton-secondary"] p,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="baseButton-secondary"] span,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-secondary"] p,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-secondary"] span {
    color: #334155 !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-secondary"]:hover,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-secondary"]:hover,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="baseButton-secondary"]:hover,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-secondary"]:hover {
    background-color: #fff1f2 !important;
    border-color: #fb7185 !important;
    color: #991b1b !important;
    box-shadow: 0 2px 10px rgba(220, 38, 38, 0.18) !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-secondary"]:hover p,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-secondary"]:hover span,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-secondary"]:hover p,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-secondary"]:hover span,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="baseButton-secondary"]:hover p,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="baseButton-secondary"]:hover span,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-secondary"]:hover p,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-secondary"]:hover span {
    color: #991b1b !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] [role="button"] {
    min-height: 3rem !important;
    background-color: #dc2626 !important;
    background-image: none !important;
    color: #ffffff !important;
    border: 1px solid #b91c1c !important;
    font-weight: 600 !important;
    transition: background-color 0.15s ease, border-color 0.15s ease, box-shadow 0.15s ease !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] [role="button"] p,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] [role="button"] span {
    color: #ffffff !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] [role="button"]:hover {
    background-color: #f87171 !important;
    border-color: #fca5a5 !important;
    box-shadow: 0 0 0 3px rgba(248, 113, 113, 0.45) !important;
    background-image: none !important;
    filter: none !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-primary"]:focus-visible,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-primary"]:focus-visible,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="baseButton-primary"]:focus-visible,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-primary"]:focus-visible {
    box-shadow: 0 0 0 2px #ffffff, 0 0 0 4px #dc2626 !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-secondary"]:focus-visible,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-secondary"]:focus-visible,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="baseButton-secondary"]:focus-visible,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-secondary"]:focus-visible {
    box-shadow: 0 0 0 2px #ffffff, 0 0 0 4px #94a3b8 !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="baseButton-primary"],
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="stBaseButton-primary"] {
    min-height: 3rem !important;
    background-color: #dc2626 !important;
    background-image: none !important;
    color: #ffffff !important;
    border: 1px solid #b91c1c !important;
    font-weight: 600 !important;
    transition: background-color 0.15s ease, border-color 0.15s ease, box-shadow 0.15s ease !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="baseButton-primary"]:hover,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="stBaseButton-primary"]:hover {
    background-color: #f87171 !important;
    border-color: #fca5a5 !important;
    box-shadow: 0 0 0 3px rgba(248, 113, 113, 0.45) !important;
    background-image: none !important;
    filter: none !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="baseButton-secondary"],
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="stBaseButton-secondary"] {
    min-height: 3rem !important;
    background-color: #ffffff !important;
    background-image: none !important;
    color: #334155 !important;
    border: 1px solid #cbd5e1 !important;
    font-weight: 600 !important;
    transition: background-color 0.15s ease, border-color 0.15s ease, box-shadow 0.15s ease, color 0.15s ease !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="baseButton-secondary"]:hover,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="stBaseButton-secondary"]:hover {
    background-color: #fff1f2 !important;
    border-color: #fb7185 !important;
    color: #991b1b !important;
    box-shadow: 0 2px 10px rgba(220, 38, 38, 0.18) !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="baseButton-primary"] p,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="baseButton-primary"] span,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="stBaseButton-primary"] p,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="stBaseButton-primary"] span {
    color: #ffffff !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="baseButton-secondary"] p,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="baseButton-secondary"] span,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="stBaseButton-secondary"] p,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="stBaseButton-secondary"] span {
    color: #334155 !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="baseButton-secondary"]:hover p,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="baseButton-secondary"]:hover span,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="stBaseButton-secondary"]:hover p,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="stBaseButton-secondary"]:hover span {
    color: #991b1b !important;
}
/* Landing primary: white nested labels; parent-hover + kind=primary (Streamlit variance) */
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-primary"] *,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-primary"] *,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="baseButton-primary"] *,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-primary"] *,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="baseButton-primary"] *,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="stBaseButton-primary"] *,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[kind="primary"] *,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[kind="primary"] *,
:is(#root, body:has(.bnd-landing-page)) div[data-testid="stButton"] button[data-testid="baseButton-primary"] *,
:is(#root, body:has(.bnd-landing-page)) div[data-testid="stButton"] button[data-testid="stBaseButton-primary"] *,
:is(#root, body:has(.bnd-landing-page)) div.row-widget.stButton[data-testid="stButton"] button[data-testid="baseButton-primary"] *,
:is(#root, body:has(.bnd-landing-page)) div.row-widget.stButton[data-testid="stButton"] button[data-testid="stBaseButton-primary"] * {
    color: #ffffff !important;
    fill: #ffffff !important;
    stroke: #ffffff !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="baseButton-primary"]:hover *,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] button[data-testid="stBaseButton-primary"]:hover *,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="baseButton-primary"]:hover *,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div.stButton button[data-testid="stBaseButton-primary"]:hover *,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="baseButton-primary"]:hover *,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] button[data-testid="stBaseButton-primary"]:hover * {
    color: #ffffff !important;
    fill: #ffffff !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] [data-testid="stButton"]:hover button[data-testid="baseButton-primary"],
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] [data-testid="stButton"]:hover button[data-testid="stBaseButton-primary"],
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] [data-testid="stButton"]:hover button[kind="primary"],
:is(#root, body:has(.bnd-landing-page)) .stApp div.row-widget.stButton[data-testid="stButton"]:hover button[data-testid="baseButton-primary"],
:is(#root, body:has(.bnd-landing-page)) .stApp div.row-widget.stButton[data-testid="stButton"]:hover button[data-testid="stBaseButton-primary"],
:is(#root, body:has(.bnd-landing-page)) div.row-widget.stButton[data-testid="stButton"]:hover button[data-testid="baseButton-primary"],
:is(#root, body:has(.bnd-landing-page)) div.row-widget.stButton[data-testid="stButton"]:hover button[data-testid="stBaseButton-primary"],
:is(#root, body:has(.bnd-landing-page)) div[data-testid="stButton"]:hover button[data-testid="baseButton-primary"],
:is(#root, body:has(.bnd-landing-page)) div[data-testid="stButton"]:hover button[data-testid="stBaseButton-primary"] {
    background-color: #f87171 !important;
    border-color: #fca5a5 !important;
    box-shadow: 0 0 0 3px rgba(248, 113, 113, 0.45) !important;
    background-image: none !important;
    color: #ffffff !important;
    filter: none !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] [data-testid="stButton"]:hover button[data-testid="baseButton-primary"] *,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] [data-testid="stButton"]:hover button[data-testid="stBaseButton-primary"] *,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] [data-testid="stButton"]:hover button[kind="primary"] *,
:is(#root, body:has(.bnd-landing-page)) .stApp div.row-widget.stButton[data-testid="stButton"]:hover button[data-testid="baseButton-primary"] *,
:is(#root, body:has(.bnd-landing-page)) .stApp div.row-widget.stButton[data-testid="stButton"]:hover button[data-testid="stBaseButton-primary"] *,
:is(#root, body:has(.bnd-landing-page)) div.row-widget.stButton[data-testid="stButton"]:hover button[data-testid="baseButton-primary"] *,
:is(#root, body:has(.bnd-landing-page)) div.row-widget.stButton[data-testid="stButton"]:hover button[data-testid="stBaseButton-primary"] *,
:is(#root, body:has(.bnd-landing-page)) div[data-testid="stButton"]:hover button[data-testid="baseButton-primary"] *,
:is(#root, body:has(.bnd-landing-page)) div[data-testid="stButton"]:hover button[data-testid="stBaseButton-primary"] * {
    color: #ffffff !important;
    fill: #ffffff !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"]:hover [role="button"] {
    background-color: #f87171 !important;
    border-color: #fca5a5 !important;
    box-shadow: 0 0 0 3px rgba(248, 113, 113, 0.45) !important;
    background-image: none !important;
    filter: none !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"]:hover [role="button"] * {
    color: #ffffff !important;
    fill: #ffffff !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] div[data-testid="stButton"] [role="button"] * {
    color: #ffffff !important;
    fill: #ffffff !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stSidebar"],
:is(#root, body:has(.bnd-landing-page)) [data-testid="stSidebarCollapsedControl"],
:is(#root, body:has(.bnd-landing-page)) [data-testid="collapsedControl"] {
    display: none !important;
}
/* Fourth workflow CTA — deep slate / navy (distinct from live red primaries) */
:is(#root, body:has(.bnd-landing-page)) .bnd-landing-sim-cta-sub {
    text-align: center;
    font-size: 0.95rem;
    color: #475569;
    margin: 1.65rem 0 0.5rem 0;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] .bnd-land-sim-btn-marker ~ div [data-testid="stButton"] button[data-testid="baseButton-secondary"],
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] .bnd-land-sim-btn-marker ~ div [data-testid="stButton"] button[data-testid="stBaseButton-secondary"],
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] .bnd-land-sim-btn-marker ~ div [data-testid="stButton"] [role="button"] {
    min-height: 3rem !important;
    background: linear-gradient(180deg, #1e3a5f 0%, #0f172a 100%) !important;
    background-image: none !important;
    color: #f8fafc !important;
    border: 1px solid #334155 !important;
    font-weight: 600 !important;
    transition: background-color 0.15s ease, border-color 0.15s ease, box-shadow 0.15s ease !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] .bnd-land-sim-btn-marker ~ div [data-testid="stButton"] button[data-testid="baseButton-secondary"] *,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] .bnd-land-sim-btn-marker ~ div [data-testid="stButton"] button[data-testid="stBaseButton-secondary"] *,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] .bnd-land-sim-btn-marker ~ div [data-testid="stButton"] [role="button"] * {
    color: #f8fafc !important;
    fill: #f8fafc !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] .bnd-land-sim-btn-marker ~ div [data-testid="stButton"] button[data-testid="baseButton-secondary"]:hover,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] .bnd-land-sim-btn-marker ~ div [data-testid="stButton"] button[data-testid="stBaseButton-secondary"]:hover,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] .bnd-land-sim-btn-marker ~ div [data-testid="stButton"]:hover [role="button"] {
    background: linear-gradient(180deg, #2563eb 0%, #1e3a5f 55%, #0f172a 100%) !important;
    border-color: #38bdf8 !important;
    box-shadow: 0 0 0 3px rgba(56, 189, 248, 0.35) !important;
    filter: none !important;
}
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] .bnd-land-sim-btn-marker ~ div [data-testid="stButton"] button[data-testid="baseButton-secondary"]:focus-visible,
:is(#root, body:has(.bnd-landing-page)) [data-testid="stMain"] .bnd-land-sim-btn-marker ~ div [data-testid="stButton"] button[data-testid="stBaseButton-secondary"]:focus-visible {
    box-shadow: 0 0 0 2px #f8fafc, 0 0 0 4px #38bdf8 !important;
}
</style>
"""
    st.markdown(_landing_style, unsafe_allow_html=True)
    if hasattr(st, 'html'):
        st.html(_landing_style)


def _cb_h5_explorer_trial_changed() -> None:
    st.session_state.pop('gui_h5_explore_catalog', None)
    st.session_state.pop('gui_h5_explore_notes', None)


_SIM_CMP_COLOR_PAIRS = (
    ('Navy', '#0f172a'),
    ('Blue', '#2563eb'),
    ('Teal', '#0d9488'),
    ('Orange', '#ea580c'),
    ('Green', '#16a34a'),
    ('Purple', '#7c3aed'),
    ('Coral', '#f43f5e'),
    ('Slate', '#475569'),
)
_SIM_CMP_COLOR_LABELS = [p[0] for p in _SIM_CMP_COLOR_PAIRS]
_SIM_CMP_COLOR_HEX = {p[0]: p[1] for p in _SIM_CMP_COLOR_PAIRS}

_SIM_CMP_OSC_VAR_OPTS = (
    ('t', 'Time (s)'),
    ('delta_mm', 'Tip displacement δ (mm)'),
    ('F', 'Force F (N)'),
    ('theta_deg', 'Bending angle θ (deg)'),
    ('M', 'Root moment M (N·m)'),
    ('kappa', 'Curvature κ (1/m)'),
    ('delta_dot_mm_s', 'Tip velocity δ̇ (mm/s)'),
)
_SIM_CMP_OSC_VAR_LABELS = {k: v for k, v in _SIM_CMP_OSC_VAR_OPTS}

_SIM_CMP_QS_VAR_OPTS = (
    ('index', 'Sample index'),
    ('delta_mm', 'Tip displacement δ (mm)'),
    ('F', 'Force F (N)'),
    ('theta_deg', 'Bending angle θ (deg)'),
    ('M', 'Root moment M (N·m)'),
    ('kappa', 'Curvature κ (1/m)'),
)
_SIM_CMP_QS_VAR_LABELS = {k: v for k, v in _SIM_CMP_QS_VAR_OPTS}


def _sim_cmp_resolve_osc_series(s: dict, var_id: str) -> np.ndarray:
    if var_id == 'delta_dot_mm_s':
        return np.asarray(s['delta_dot_m_s'], dtype=float) * 1000.0
    return np.asarray(s[var_id], dtype=float)


def _sim_cmp_qs_series_bundle(
    *,
    length_mm: float,
    width_mm: float,
    material_key: str,
    d_mm: np.ndarray,
    F: np.ndarray,
) -> dict:
    d_mm = np.asarray(d_mm, dtype=float)
    F = np.asarray(F, dtype=float)
    L_m = max(float(length_mm), 1.0) / 1000.0
    k, E, I = cantilever_stiffness_N_per_m_for_geometry(
        length_mm=length_mm,
        width_mm=width_mm,
        material_key=material_key,
    )
    d_m = d_mm / 1000.0
    theta_rad = d_m / max(L_m, 1e-12)
    theta_deg = np.rad2deg(theta_rad)
    M = F * L_m
    kappa = M / (E * I + 1e-30)
    return {
        'index': np.arange(len(d_mm), dtype=float),
        'delta_mm': d_mm,
        'F': F,
        'theta_deg': theta_deg,
        'M': M,
        'kappa': kappa,
    }


def _render_simulation_comparison_page() -> None:
    """
    Community / teaching view: two independent cantilever simulations (Specimen 1 & 2),
    same axes, numpy-only (no Bender / NI-DAQ).
    """
    st.session_state.setdefault('gui_sim_cmp_a_mat', 'polyurethane')
    st.session_state.setdefault('gui_sim_cmp_b_mat', 'silicone')
    st.session_state.setdefault('gui_sim_cmp_L_a', 30)
    st.session_state.setdefault('gui_sim_cmp_L_b', 45)
    st.session_state.setdefault('gui_sim_cmp_w_a', 25.4)
    st.session_state.setdefault('gui_sim_cmp_w_b', 25.4)
    st.session_state.setdefault('gui_sim_cmp_delta_max', 10.0)
    st.session_state.setdefault('gui_sim_cmp_noise', 0.025)
    st.session_state.setdefault('gui_sim_cmp_oscillatory', False)
    st.session_state.setdefault('gui_sim_cmp_freq_hz', 1.0)
    st.session_state.setdefault('gui_sim_cmp_amp_mm', 2.0)
    st.session_state.setdefault('gui_sim_cmp_cycles', 3.5)
    st.session_state.setdefault('gui_sim_cmp_ppc', 200)
    st.session_state.setdefault('gui_sim_cmp_color_a', 'Navy')
    st.session_state.setdefault('gui_sim_cmp_color_b', 'Blue')
    st.session_state.setdefault('gui_sim_cmp_osc_x', 'delta_mm')
    st.session_state.setdefault('gui_sim_cmp_osc_y', 'F')
    st.session_state.setdefault('gui_sim_cmp_qs_x', 'delta_mm')
    st.session_state.setdefault('gui_sim_cmp_qs_y', 'F')

    st.subheader('Simulation & Comparison')
    st.caption(
        '**Workbench layout:** parameters on the **left**, plots on the **right** — they update together on each slider '
        'change. Widen the browser so both columns stay side-by-side.'
    )

    _mat_labels = {
        'polyurethane': 'Polyurethane (E ≈ 35 MPa)',
        'silicone': 'Silicone (E ≈ 3 MPa)',
    }

    st.checkbox(
        'Oscillatory viscoelastic mode',
        key='gui_sim_cmp_oscillatory',
        help='Sinusoidal δ(t); Kelvin–Voigt F = kδ + cδ̇ with material-dependent damping (PU = wide loop, silicone = thin).',
    )
    _osc = bool(st.session_state.get('gui_sim_cmp_oscillatory', False))
    if _osc:
        st.markdown(
            '<div class="bnd-sim-osc-banner">Oscillatory mode — try **Frequency** and **Amplitude** on the left; use **X / Y** '
            'below the table to compare variables. Higher f tends to open the loop (more viscous work per cycle).</div>',
            unsafe_allow_html=True,
        )
    else:
        st.caption('Quasi-static ramp: **F ≈ k·δ** with noise. Stiffness **k ∝ D⁴/L³** (solid circular section).')

    try:
        c_lo, c_hi = st.columns([0.36, 0.64], gap='small')
    except TypeError:
        c_lo, c_hi = st.columns([0.36, 0.64])

    with c_lo:
        st.markdown('<div class="bnd-sim-compare-workbench" aria-hidden="true"></div>', unsafe_allow_html=True)
        with st.container(border=True):
            st.markdown('**Specimen 1**')
            st.selectbox(
                'Material',
                options=['polyurethane', 'silicone'],
                key='gui_sim_cmp_a_mat',
                format_func=lambda k: _mat_labels[k],
            )
            st.slider('L₁ length (mm)', min_value=5.0, max_value=150.0, step=0.5, key='gui_sim_cmp_L_a')
            st.slider('D₁ tube OD (mm)', min_value=2.0, max_value=50.0, step=0.1, key='gui_sim_cmp_w_a')
            st.selectbox(
                'Line color (Specimen 1)',
                options=_SIM_CMP_COLOR_LABELS,
                key='gui_sim_cmp_color_a',
                format_func=lambda x: x,
            )
        with st.container(border=True):
            st.markdown('**Specimen 2**')
            st.selectbox(
                'Material ',
                options=['polyurethane', 'silicone'],
                key='gui_sim_cmp_b_mat',
                format_func=lambda k: _mat_labels[k],
            )
            st.slider('L₂ length (mm)', min_value=5.0, max_value=150.0, step=0.5, key='gui_sim_cmp_L_b')
            st.slider('D₂ tube OD (mm)', min_value=2.0, max_value=50.0, step=0.1, key='gui_sim_cmp_w_b')
            st.selectbox(
                'Line color (Specimen 2)',
                options=_SIM_CMP_COLOR_LABELS,
                key='gui_sim_cmp_color_b',
                format_func=lambda x: x,
            )
        with st.container(border=True):
            if _osc:
                st.markdown('**Oscillatory drive**')
                st.slider('Frequency f (Hz)', min_value=0.1, max_value=8.0, step=0.05, key='gui_sim_cmp_freq_hz')
                st.slider('Amplitude A (mm, peak)', min_value=0.2, max_value=12.0, step=0.1, key='gui_sim_cmp_amp_mm')
                st.slider('Cycles shown', min_value=2.0, max_value=10.0, step=0.5, key='gui_sim_cmp_cycles')
                st.slider('Samples / cycle', min_value=80, max_value=400, step=20, key='gui_sim_cmp_ppc')
                st.slider('Force noise σ (N)', min_value=0.0, max_value=0.2, step=0.005, key='gui_sim_cmp_noise')
            else:
                st.markdown('**Quasi-static sweep**')
                st.slider('Max tip δ (mm)', min_value=0.5, max_value=25.0, step=0.5, key='gui_sim_cmp_delta_max')
                st.slider('Force noise σ (N)', min_value=0.0, max_value=0.2, step=0.005, key='gui_sim_cmp_noise')

    La = float(st.session_state['gui_sim_cmp_L_a'])
    Lb = float(st.session_state['gui_sim_cmp_L_b'])
    wa = float(st.session_state['gui_sim_cmp_w_a'])
    wb = float(st.session_state['gui_sim_cmp_w_b'])
    ma = str(st.session_state['gui_sim_cmp_a_mat'])
    mb = str(st.session_state['gui_sim_cmp_b_mat'])
    sigma = float(st.session_state['gui_sim_cmp_noise'])

    k_a, _, _ = cantilever_stiffness_N_per_m_for_geometry(length_mm=La, width_mm=wa, material_key=ma)
    k_b, _, _ = cantilever_stiffness_N_per_m_for_geometry(length_mm=Lb, width_mm=wb, material_key=mb)
    _kr = (k_a / k_b) if k_b > 1e-30 else float('nan')

    with c_hi:
        if _osc:
            f_hz = float(st.session_state['gui_sim_cmp_freq_hz'])
            amp_mm = float(st.session_state['gui_sim_cmp_amp_mm'])
            n_cyc = float(st.session_state['gui_sim_cmp_cycles'])
            ppc = int(st.session_state['gui_sim_cmp_ppc'])

            rng_a = np.random.default_rng(201)
            rng_b = np.random.default_rng(202)
            sa = oscillatory_viscoelastic_timeseries(
                length_mm=La,
                width_mm=wa,
                material_key=ma,
                freq_hz=f_hz,
                amplitude_mm=amp_mm,
                n_cycles_shown=n_cyc,
                points_per_cycle=ppc,
                noise_std_N=max(sigma, 0.0),
                rng=rng_a,
            )
            sb = oscillatory_viscoelastic_timeseries(
                length_mm=Lb,
                width_mm=wb,
                material_key=mb,
                freq_hz=f_hz,
                amplitude_mm=amp_mm,
                n_cycles_shown=n_cyc,
                points_per_cycle=ppc,
                noise_std_N=max(sigma, 0.0),
                rng=rng_b,
            )

            _e1 = sa['energy_dissipated_last_cycle_J']
            _e2 = sb['energy_dissipated_last_cycle_J']
            st.markdown('**Results** (updates with sliders)')
            st.dataframe(
                pd.DataFrame(
                    [
                        {
                            'Spec': '1',
                            'k (N/m)': f'{k_a:,.0f}',
                            'c (N·s/m)': f'{sa["c"]:,.0f}',
                            'φ (°)': f'{sa["phase_lag_deg"]:.1f}',
                            'J/cycle': f'{_e1:.4f}' if np.isfinite(_e1) else '—',
                            'tan φ': f'{sa["tan_delta"]:.3f}' if np.isfinite(sa['tan_delta']) else '—',
                        },
                        {
                            'Spec': '2',
                            'k (N/m)': f'{k_b:,.0f}',
                            'c (N·s/m)': f'{sb["c"]:,.0f}',
                            'φ (°)': f'{sb["phase_lag_deg"]:.1f}',
                            'J/cycle': f'{_e2:.4f}' if np.isfinite(_e2) else '—',
                            'tan φ': f'{sb["tan_delta"]:.3f}' if np.isfinite(sb['tan_delta']) else '—',
                        },
                        {'Spec': 'k1/k2', 'k (N/m)': f'{_kr:.3f}', 'c (N·s/m)': '—', 'φ (°)': '—', 'J/cycle': '—', 'tan φ': '—'},
                    ]
                ),
                use_container_width=True,
                hide_index=True,
            )

            _osc_ids = [p[0] for p in _SIM_CMP_OSC_VAR_OPTS]
            _px, _py = st.columns(2, gap='small')
            with _px:
                st.selectbox(
                    'X axis',
                    options=_osc_ids,
                    key='gui_sim_cmp_osc_x',
                    format_func=lambda i: _SIM_CMP_OSC_VAR_LABELS[i],
                )
            with _py:
                st.selectbox(
                    'Y axis',
                    options=_osc_ids,
                    key='gui_sim_cmp_osc_y',
                    format_func=lambda i: _SIM_CMP_OSC_VAR_LABELS[i],
                )

            _xid = str(st.session_state.get('gui_sim_cmp_osc_x') or 'delta_mm')
            _yid = str(st.session_state.get('gui_sim_cmp_osc_y') or 'F')
            if _xid not in _SIM_CMP_OSC_VAR_LABELS:
                _xid = 'delta_mm'
            if _yid not in _SIM_CMP_OSC_VAR_LABELS:
                _yid = 'F'

            _c1 = _SIM_CMP_COLOR_HEX.get(str(st.session_state.get('gui_sim_cmp_color_a', 'Navy')), '#0f172a')
            _c2 = _SIM_CMP_COLOR_HEX.get(str(st.session_state.get('gui_sim_cmp_color_b', 'Blue')), '#2563eb')
            _xa = _sim_cmp_resolve_osc_series(sa, _xid)
            _ya = _sim_cmp_resolve_osc_series(sa, _yid)
            _xb = _sim_cmp_resolve_osc_series(sb, _xid)
            _yb = _sim_cmp_resolve_osc_series(sb, _yid)

            fig = go.Figure()
            fig.add_trace(
                go.Scatter(
                    x=_xa,
                    y=_ya,
                    mode='lines',
                    name='Specimen 1',
                    line=dict(color=_c1, width=2.2),
                )
            )
            fig.add_trace(
                go.Scatter(
                    x=_xb,
                    y=_yb,
                    mode='lines',
                    name='Specimen 2',
                    line=dict(color=_c2, width=2.2),
                )
            )
            fig.update_layout(
                title=dict(
                    text=f'Oscillatory · f = {f_hz:.2f} Hz, A = {amp_mm:.2f} mm',
                    font=dict(size=14, color='#0f172a'),
                ),
                xaxis_title=_SIM_CMP_OSC_VAR_LABELS[_xid],
                yaxis_title=_SIM_CMP_OSC_VAR_LABELS[_yid],
                height=480,
                margin=dict(l=52, r=28, t=52, b=44),
                paper_bgcolor='#ffffff',
                plot_bgcolor='#f8fafc',
                hovermode='closest',
                legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='right', x=1),
            )
            st.plotly_chart(fig, use_container_width=True, config={'scrollZoom': True, 'displayModeBar': True})

            with st.expander('Equations & model detail (oscillatory)', expanded=False):
                st.markdown(
                    '- **δ(t) = A sin(2πft)**. **F = kδ + cδ̇** (Kelvin–Voigt); **c** is much larger for polyurethane than '
                    'silicone (wider hysteresis).\n'
                    '- **Phase lag:** **F** lags **δ** by **φ = arctan(cω/k)**.\n'
                    '- **Dissipation:** **J/cycle** ≈ ∫ F δ̇ dt over the last full period.\n'
                    '- **M = F·L**, **θ ≈ δ/L**, **κ = M/(EI)**.'
                )

        else:
            dmax = float(st.session_state['gui_sim_cmp_delta_max'])
            st.markdown('**Results** (updates with sliders)')
            st.caption(f'k1 = {k_a:,.0f} N/m · k2 = {k_b:,.0f} N/m · k1/k2 = {_kr:.3f}')

            rng_a = np.random.default_rng(201)
            rng_b = np.random.default_rng(202)
            d1, F1 = force_displacement_comparison_curve(
                length_mm=La,
                width_mm=wa,
                material_key=ma,
                delta_max_mm=dmax,
                n_points=180,
                noise_std_N=max(sigma, 0.0),
                rng=rng_a,
            )
            d2, F2 = force_displacement_comparison_curve(
                length_mm=Lb,
                width_mm=wb,
                material_key=mb,
                delta_max_mm=dmax,
                n_points=180,
                noise_std_N=max(sigma, 0.0),
                rng=rng_b,
            )

            _ba = _sim_cmp_qs_series_bundle(
                length_mm=La, width_mm=wa, material_key=ma, d_mm=d1, F=F1
            )
            _bb = _sim_cmp_qs_series_bundle(
                length_mm=Lb, width_mm=wb, material_key=mb, d_mm=d2, F=F2
            )

            _qs_ids = [p[0] for p in _SIM_CMP_QS_VAR_OPTS]
            _qx, _qy = st.columns(2, gap='small')
            with _qx:
                st.selectbox(
                    'X axis',
                    options=_qs_ids,
                    key='gui_sim_cmp_qs_x',
                    format_func=lambda i: _SIM_CMP_QS_VAR_LABELS[i],
                )
            with _qy:
                st.selectbox(
                    'Y axis',
                    options=_qs_ids,
                    key='gui_sim_cmp_qs_y',
                    format_func=lambda i: _SIM_CMP_QS_VAR_LABELS[i],
                )

            _qxid = str(st.session_state.get('gui_sim_cmp_qs_x') or 'delta_mm')
            _qyid = str(st.session_state.get('gui_sim_cmp_qs_y') or 'F')
            if _qxid not in _SIM_CMP_QS_VAR_LABELS:
                _qxid = 'delta_mm'
            if _qyid not in _SIM_CMP_QS_VAR_LABELS:
                _qyid = 'F'

            _qc1 = _SIM_CMP_COLOR_HEX.get(str(st.session_state.get('gui_sim_cmp_color_a', 'Navy')), '#0f172a')
            _qc2 = _SIM_CMP_COLOR_HEX.get(str(st.session_state.get('gui_sim_cmp_color_b', 'Blue')), '#2563eb')
            _qxa = np.asarray(_ba[_qxid], dtype=float)
            _qya = np.asarray(_ba[_qyid], dtype=float)
            _qxb = np.asarray(_bb[_qxid], dtype=float)
            _qyb = np.asarray(_bb[_qyid], dtype=float)

            fig = go.Figure()
            fig.add_trace(
                go.Scatter(
                    x=_qxa,
                    y=_qya,
                    mode='lines',
                    name='Specimen 1',
                    line=dict(color=_qc1, width=2.2),
                )
            )
            fig.add_trace(
                go.Scatter(
                    x=_qxb,
                    y=_qyb,
                    mode='lines',
                    name='Specimen 2',
                    line=dict(color=_qc2, width=2.2),
                )
            )
            fig.update_layout(
                title='Quasi-static comparison (choose X / Y above)',
                xaxis_title=_SIM_CMP_QS_VAR_LABELS[_qxid],
                yaxis_title=_SIM_CMP_QS_VAR_LABELS[_qyid],
                height=440,
                margin=dict(l=52, r=28, t=48, b=44),
                hovermode='x unified',
                legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='right', x=1),
                paper_bgcolor='#ffffff',
                plot_bgcolor='#f8fafc',
            )
            st.plotly_chart(fig, use_container_width=True, config={'scrollZoom': True, 'displayModeBar': True})

            with st.expander('Model detail (quasi-static)', expanded=False):
                st.markdown(
                    '- **Geometry:** solid circle, diameter D (mm) → I = πD⁴/64.\n'
                    '- **Material:** E from polyurethane (~35 MPa) or silicone (~3 MPa).\n'
                    '- **Stiffness:** k = 3EI/L³ with L in meters.\n'
                    '- **Curve:** F = k·δ (δ in m) + Normal(0, σ) on F in newtons.'
                )


def _path_ch_axis_label(path: str, ch) -> str:
    if not path:
        return ''
    if ch is None or ch == '':
        return path
    try:
        i = int(ch)
        if 0 <= i < len(FT_ROW_LABELS):
            return f'{path} — {FT_ROW_LABELS[i]}'
    except (TypeError, ValueError):
        pass
    return path


def _h5_attr_row_token(name: str) -> str:
    return hashlib.sha256(name.encode('utf-8', errors='replace')).hexdigest()[:16]


def _h5_custom_axis_folder_suffix(axis: str, cwd: str) -> str:
    """Per-axis + per-folder suffix so selectbox keys stay valid when folder or axis context changes."""
    return hashlib.sha256(f'{axis}|{cwd or ""}'.encode('utf-8', errors='replace')).hexdigest()[:16]


def _flush_pending_h5_custom_axis_nav(axis: str) -> None:
    """Apply folder navigation queued by Enter / Up / Root before axis widgets mount."""
    pk = f'gui_pending_h5_custom_nav_{axis}'
    if pk not in st.session_state:
        return
    op = st.session_state.pop(pk)
    cwd_key = f'gui_h5_custom_cwd_{axis}'
    cwd = str(st.session_state.get(cwd_key) or '')
    if op == 'root':
        st.session_state[cwd_key] = ''
    elif op == 'up':
        parts = cwd.split('/')
        st.session_state[cwd_key] = '/'.join(parts[:-1]) if parts else ''
    elif isinstance(op, str) and op.startswith('enter:'):
        pick = op[6:]
        if pick and pick != '—':
            st.session_state[cwd_key] = h5_join_internal_path(cwd, pick)


def _render_h5_custom_axis_browser(loaded: str, axis: str, heading: str) -> None:
    """One column: independent cwd, subfolder nav, dataset pick, optional 6×N channel."""
    _flush_pending_h5_custom_axis_nav(axis)
    cwd_key = f'gui_h5_custom_cwd_{axis}'
    st.session_state.setdefault(cwd_key, '')
    cwd = str(st.session_state.get(cwd_key) or '')

    kids = list_h5_group_children(loaded, cwd)
    if not kids and cwd:
        st.session_state[cwd_key] = ''
        cwd = ''
        kids = list_h5_group_children(loaded, cwd)

    groups = sorted([e['name'] for e in kids if e['kind'] == 'group'], key=str.lower)
    ds_list = [e for e in kids if e['kind'] == 'dataset']
    _ks = _h5_custom_axis_folder_suffix(axis, cwd)

    st.markdown(f'**{heading}**')
    st.caption('Folder: `' + (cwd or '/') + '`')

    nav1, nav2, nav3 = st.columns(3)
    with nav1:
        pick = st.selectbox(
            'Subfolder',
            options=['—'] + groups,
            key=f'gui_h5_custom_{axis}_sub_{_ks}',
        )
        if st.button('Enter', key=f'gui_h5_custom_{axis}_enter_{_ks}', disabled=pick == '—'):
            st.session_state[f'gui_pending_h5_custom_nav_{axis}'] = f'enter:{pick}'
            st.rerun()
    with nav2:
        if st.button('Up', key=f'gui_h5_custom_{axis}_up_{_ks}', disabled=not cwd):
            st.session_state[f'gui_pending_h5_custom_nav_{axis}'] = 'up'
            st.rerun()
    with nav3:
        if st.button('Root', key=f'gui_h5_custom_{axis}_root_{_ks}'):
            st.session_state[f'gui_pending_h5_custom_nav_{axis}'] = 'root'
            st.rerun()

    plottable = [e for e in ds_list if e['plot'] in ('1d', 'six')]
    opt = ['—'] + [e['name'] for e in plottable]

    st.selectbox(
        'Dataset',
        options=opt,
        key=f'gui_h5_custom_{axis}_ds_{_ks}',
    )
    ds_sel = str(st.session_state.get(f'gui_h5_custom_{axis}_ds_{_ks}', '—'))
    plot_mode = next((e['plot'] for e in plottable if e['name'] == ds_sel), None) if ds_sel != '—' else None
    if plot_mode == 'six':
        ch_key = f'gui_h5_custom_{axis}_ch_{_ks}'
        st.session_state.setdefault(ch_key, 5)
        st.selectbox(
            '6×N channel',
            options=list(range(6)),
            format_func=lambda i: FT_ROW_LABELS[i],
            key=ch_key,
        )


def _h5_custom_swap_xy_session() -> None:
    """Exchange X/Y folder, dataset, and channel widget state (custom hierarchy mode)."""
    cx = str(st.session_state.get('gui_h5_custom_cwd_x') or '')
    cy = str(st.session_state.get('gui_h5_custom_cwd_y') or '')
    ksx = _h5_custom_axis_folder_suffix('x', cx)
    ksy = _h5_custom_axis_folder_suffix('y', cy)
    ds_x = str(st.session_state.get(f'gui_h5_custom_x_ds_{ksx}', '—'))
    ds_y = str(st.session_state.get(f'gui_h5_custom_y_ds_{ksy}', '—'))
    ch_x = st.session_state.get(f'gui_h5_custom_x_ch_{ksx}', 5)
    ch_y = st.session_state.get(f'gui_h5_custom_y_ch_{ksy}', 5)
    try:
        ch_x_i = int(ch_x)
    except (TypeError, ValueError):
        ch_x_i = 5
    try:
        ch_y_i = int(ch_y)
    except (TypeError, ValueError):
        ch_y_i = 5
    ch_x_i = max(0, min(5, ch_x_i))
    ch_y_i = max(0, min(5, ch_y_i))

    st.session_state['gui_h5_custom_cwd_x'] = cy
    st.session_state['gui_h5_custom_cwd_y'] = cx

    ksx_new = _h5_custom_axis_folder_suffix('x', cy)
    ksy_new = _h5_custom_axis_folder_suffix('y', cx)
    st.session_state[f'gui_h5_custom_x_ds_{ksx_new}'] = ds_y
    st.session_state[f'gui_h5_custom_y_ds_{ksy_new}'] = ds_x
    st.session_state[f'gui_h5_custom_x_ch_{ksx_new}'] = ch_y_i
    st.session_state[f'gui_h5_custom_y_ch_{ksy_new}'] = ch_x_i


def _h5_custom_resolve_axis_path(loaded: str, axis: str) -> Tuple[str, Optional[int], Optional[str]]:
    """Return ``(internal_dataset_path, channel_or_none, plot_mode)`` for axis ``'x'`` or ``'y'``."""
    cwd_key = f'gui_h5_custom_cwd_{axis}'
    cwd = str(st.session_state.get(cwd_key) or '')
    kids = list_h5_group_children(loaded, cwd)
    if not kids and cwd:
        cwd = ''
        kids = list_h5_group_children(loaded, cwd)
    ds_list = [e for e in kids if e['kind'] == 'dataset']
    plottable = [e for e in ds_list if e['plot'] in ('1d', 'six')]
    opt_names = [e['name'] for e in plottable]
    _ks = _h5_custom_axis_folder_suffix(axis, cwd)
    ds_sel = str(st.session_state.get(f'gui_h5_custom_{axis}_ds_{_ks}', '—'))
    if ds_sel != '—' and ds_sel not in opt_names:
        ds_sel = '—'
    plot_mode = next((e['plot'] for e in plottable if e['name'] == ds_sel), None) if ds_sel != '—' else None
    ch: Optional[int] = None
    if plot_mode == 'six':
        ch = int(st.session_state.get(f'gui_h5_custom_{axis}_ch_{_ks}', 5))
    path = h5_join_internal_path(cwd, ds_sel) if ds_sel != '—' else ''
    return path, ch, plot_mode


def _render_h5_custom_hierarchy_dropdowns(loaded: str) -> None:
    """Non-standard HDF5: **X** and **Y** each have their own folder navigation and dataset pick."""
    if st.session_state.pop('gui_pending_h5_custom_swap_xy', False):
        _h5_custom_swap_xy_session()
    st.caption(
        '**X** and **Y** use **separate folders**. In each column: open subfolders with **Enter**, pick a **Dataset**, '
        'and set **6×N channel** when needed (Fx…Tz).'
    )
    cx, cy = st.columns(2)
    with cx:
        _render_h5_custom_axis_browser(loaded, 'x', 'X axis (horizontal)')
    with cy:
        _render_h5_custom_axis_browser(loaded, 'y', 'Y axis (vertical)')

    _sw1, _sw2, _sw3 = st.columns([2, 1, 2])
    with _sw2:
        if st.button(
            'Swap X ↔ Y',
            key='gui_h5_custom_swap_xy',
            help='Exchange X and Y (each column’s folder, dataset, and 6×N channel).',
        ):
            st.session_state['gui_pending_h5_custom_swap_xy'] = True
            st.rerun()

    st.divider()
    x_path, xch, _ = _h5_custom_resolve_axis_path(loaded, 'x')
    y_path, ych, _ = _h5_custom_resolve_axis_path(loaded, 'y')
    if not x_path or not y_path:
        st.info('Choose a **Dataset** for both **X** and **Y** (each in its own folder column above).')
        return

    x_title = _path_ch_axis_label(x_path, xch)
    y_title = _path_ch_axis_label(y_path, ych)
    try:
        xa = read_h5_series_1d(loaded, x_path, xch)
        ya = read_h5_series_1d(loaded, y_path, ych)
    except Exception as e:
        st.error(str(e))
        return

    x_data, y_data = align_xy(xa, ya)
    n = int(min(x_data.size, y_data.size))
    if n <= 0:
        st.warning('Selected series are empty or trimming produced zero samples.')
        return

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=x_data,
            y=y_data,
            mode='lines',
            name=y_title,
            line=dict(width=1.2),
        )
    )
    fig.update_layout(
        title=f'{y_title} vs {x_title} (n={n})',
        xaxis_title=x_title,
        yaxis_title=y_title,
        margin=dict(l=48, r=24, t=48, b=48),
        hovermode='x unified',
        height=520,
    )
    st.plotly_chart(fig, use_container_width=True)


def _render_h5_attribute_editor(loaded: str) -> None:
    """Edit scalar HDF5 attributes on file root, a group, or a dataset (writable file)."""
    with st.expander('Edit HDF5 attributes', expanded=False):
        st.caption(
            'Fix wrong metadata (e.g. **dclamp**, sample rate). Requires **write permission** on the file path. '
            'Back up important data first; arrays and opaque types are read-only here.'
        )
        if 'gui_h5_attr_path_typed' not in st.session_state:
            st.session_state['gui_h5_attr_path_typed'] = ''
        if 'gui_pending_h5_attr_edit_target' in st.session_state:
            st.session_state['gui_h5_attr_edit_target'] = str(
                st.session_state.pop('gui_pending_h5_attr_edit_target') or ''
            ).strip()
        st.text_input(
            'Path inside the file (blank = file root)',
            key='gui_h5_attr_path_typed',
            placeholder='e.g. Calibrated or 02_TimeSeries/trial_1',
            help='Use forward slashes. Leave empty to edit attributes on the file itself.',
        )
        if st.button('Load attributes', key='gui_h5_attr_load_btn', type='secondary'):
            st.session_state['gui_pending_h5_attr_edit_target'] = str(
                st.session_state.get('gui_h5_attr_path_typed') or ''
            ).strip()
            st.rerun()

        if 'gui_h5_attr_edit_target' not in st.session_state:
            st.info('Enter a path (or leave blank for file root) and click **Load attributes**.')
            return

        ap = str(st.session_state.get('gui_h5_attr_edit_target') or '').strip()
        try:
            rows = list_h5_attribute_rows(loaded, ap)
        except Exception as e:
            st.error(f'Could not read attributes: {e}')
            return

        st.markdown(f'Editing **`{ap or "/"}`** — {len(rows)} attribute(s)')
        for row in rows:
            name = row['name']
            tok = _h5_attr_row_token(f'{ap}::{name}')
            r1, r2 = st.columns([5, 1])
            with r1:
                if row['editable']:
                    sk = f'gui_h5_attr_v_{tok}'
                    if sk not in st.session_state:
                        st.session_state[sk] = row['value_text']
                    st.text_input(name, key=sk)
                else:
                    st.text_input(
                        f'{name} (read-only: {row["kind"]})',
                        value=row['value_text'],
                        disabled=True,
                        key=f'gui_h5_attr_ro_{tok}',
                    )
            with r2:
                st.checkbox('Delete', key=f'gui_h5_attr_d_{tok}', value=False)

        st.divider()
        st.markdown('**Add attribute**')
        a1, a2, a3 = st.columns(3)
        with a1:
            if 'gui_h5_attr_new_name' not in st.session_state:
                st.session_state['gui_h5_attr_new_name'] = ''
            st.text_input('New name', key='gui_h5_attr_new_name')
        with a2:
            if 'gui_h5_attr_new_val' not in st.session_state:
                st.session_state['gui_h5_attr_new_val'] = ''
            st.text_input('New value', key='gui_h5_attr_new_val')
        with a3:
            if 'gui_h5_attr_new_kind' not in st.session_state:
                st.session_state['gui_h5_attr_new_kind'] = None
            st.selectbox(
                'Type',
                options=['str', 'float', 'int', 'bool'],
                key='gui_h5_attr_new_kind',
            )

        if not st.button('Apply changes to file', key='gui_h5_attr_apply', type='primary'):
            return

        deletes: list = []
        final_updates: dict = {}
        for row in rows:
            name = row['name']
            tok = _h5_attr_row_token(f'{ap}::{name}')
            if st.session_state.get(f'gui_h5_attr_d_{tok}', False):
                deletes.append(name)
            elif row['editable']:
                final_updates[name] = (
                    row['kind'],
                    str(st.session_state.get(f'gui_h5_attr_v_{tok}', row['value_text'])),
                )

        nn = str(st.session_state.get('gui_h5_attr_new_name') or '').strip()
        adds: Optional[list] = None
        if nn:
            adds = [
                (
                    nn,
                    str(st.session_state.get('gui_h5_attr_new_kind') or 'str'),
                    str(st.session_state.get('gui_h5_attr_new_val') or ''),
                )
            ]

        try:
            notes = write_h5_user_attributes(
                loaded,
                ap,
                updates=final_updates,
                delete_names=deletes,
                additions=adds,
            )
        except OSError as e:
            st.error(f'Could not write (read-only or locked file?): {e}')
            return
        except Exception as e:
            st.error(str(e))
            return

        st.success('Attributes written. Refreshing series from file…')
        for line in notes:
            st.caption(line)
        try:
            fresh = list_h5_attribute_rows(loaded, ap)
            for row in fresh:
                if row['editable']:
                    tt = _h5_attr_row_token(f'{ap}::{row["name"]}')
                    st.session_state[f'gui_h5_attr_v_{tt}'] = row['value_text']
        except OSError:
            pass
        st.session_state.pop('gui_h5_explore_schema', None)
        st.session_state.pop('gui_h5_explore_catalog', None)
        st.session_state.pop('gui_h5_explore_notes', None)
        st.rerun()


def _render_h5_explorer() -> None:
    """Visualize Data: HDF5 path or upload; standard auto-series or custom hierarchy; optional attr edit."""
    st.subheader('Visualize Data')
    st.caption(
        'Paths are read on this machine (Streamlit server). Use **Standard** for CritterGripper exports, '
        'or **Custom** to drill into any layout. Edit attributes at the bottom if metadata is wrong.'
    )

    st.session_state.setdefault('gui_h5_explore_path', '')
    st.session_state.setdefault('gui_h5_explore_trial', '')
    st.session_state.setdefault('gui_h5_explore_x', 'Time (s)')
    st.session_state.setdefault('gui_h5_explore_y', 'Primary torque corrected (N·m)')

    st.session_state.setdefault('gui_h5_explore_ui_mode', 'standard')
    st.radio(
        'File Type',
        options=['standard', 'custom'],
        format_func=lambda x: (
            'Standard — auto series (v2 / legacy / generic)'
            if x == 'standard'
            else 'Custom — navigate folders + X/Y dropdowns'
        ),
        horizontal=True,
        key='gui_h5_explore_ui_mode',
    )

    with st.container(border=True):
        path_in = st.text_input(
            'Path to `.h5` file',
            key='gui_h5_explore_path',
            placeholder=r'Example: C:\Data\experiment_001.h5',
            help='File must exist on the computer running Streamlit.',
        )
        up = st.file_uploader('Or upload an `.h5` file', type=['h5'], key='gui_h5_explore_upload')
        if up is not None:
            try:
                suf = os.path.splitext(up.name)[1] or '.h5'
                fd, tmp_path = tempfile.mkstemp(suffix=suf, prefix='crittergripper_h5_')
                os.close(fd)
                with open(tmp_path, 'wb') as wf:
                    wf.write(up.getbuffer())
                st.session_state['gui_h5_explore_upload_path'] = tmp_path
                st.caption('Temporary upload saved on the server — click **Load file** below.')
                with st.expander('Show temporary file path', expanded=False):
                    st.code(str(tmp_path))
            except Exception as e:
                st.error(f'Could not save upload: {e}')

        if st.button('Load file', key='gui_h5_explore_load', type='primary'):
            pin = str(st.session_state.get('gui_h5_explore_path') or '').strip()
            upl = str(st.session_state.get('gui_h5_explore_upload_path') or '').strip()
            chosen = ''
            if pin and os.path.isfile(pin):
                chosen = os.path.normpath(pin)
                st.session_state.pop('gui_h5_explore_upload_path', None)
            elif upl and os.path.isfile(upl):
                chosen = upl
            if chosen:
                st.session_state['gui_h5_explore_loaded_path'] = chosen
            else:
                st.session_state.pop('gui_h5_explore_loaded_path', None)
            st.session_state.pop('gui_h5_explore_catalog', None)
            st.session_state.pop('gui_h5_explore_notes', None)
            st.session_state.pop('gui_h5_explore_schema', None)
            for k in (
                'gui_h5_custom_cwd_x',
                'gui_h5_custom_cwd_y',
                'gui_h5_attr_edit_target',
            ):
                st.session_state.pop(k, None)
            st.rerun()

    loaded = str(st.session_state.get('gui_h5_explore_loaded_path') or '').strip()
    if not loaded or not os.path.isfile(loaded):
        if loaded:
            st.warning('File not found. Check the path or upload again, then **Load file**.')
        st.info('Load an `.h5` file to list series and plot.')
        return

    mode = str(st.session_state.get('gui_h5_explore_ui_mode') or 'standard')
    prev_mode = st.session_state.get('_h5_explore_mode_prev')
    if prev_mode is not None and prev_mode != mode:
        st.session_state.pop('gui_h5_explore_schema', None)
        st.session_state.pop('gui_h5_explore_catalog', None)
        st.session_state.pop('gui_h5_explore_notes', None)
    st.session_state['_h5_explore_mode_prev'] = mode

    if st.session_state.get('gui_h5_explore_schema') is None or st.session_state.get('gui_h5_explore_catalog') is None:
        if mode == 'custom':
            try:
                with h5py.File(loaded, 'r'):
                    pass
            except OSError:
                st.error('Could not open this file as HDF5 (missing, corrupted, or invalid format).')
                return
            st.session_state['gui_h5_explore_schema'] = 'custom'
            st.session_state['gui_h5_explore_catalog'] = {}
            st.session_state['gui_h5_explore_notes'] = []
        else:
            schema = detect_h5_schema(loaded)
            st.session_state['gui_h5_explore_schema'] = schema
            notes: list = []
            if schema == 'v2':
                trials = list_v2_trials(loaded)
                if not trials:
                    st.error('Schema v2 file has no `02_TimeSeries/trial_*` groups.')
                    return
                tid = str(st.session_state.get('gui_h5_explore_trial') or trials[0])
                if tid not in trials:
                    tid = trials[0]
                st.session_state['gui_h5_explore_trial'] = tid
                cat, n2 = build_series_catalog_v2(loaded, tid)
                notes.extend(n2)
            elif schema == 'legacy':
                st.session_state['gui_h5_explore_trial'] = ''
                cat, n2 = build_series_catalog_legacy(loaded)
                notes.extend(n2)
            elif schema == 'generic':
                st.session_state['gui_h5_explore_trial'] = ''
                cat, n2 = build_series_catalog_generic(loaded)
                notes.extend(n2)
            elif schema == 'browse':
                st.session_state['gui_h5_explore_trial'] = ''
                cat, n2 = {}, [
                    'No structured series list detected — use **Custom hierarchy** mode to plot.'
                ]
                notes.extend(n2)
            else:
                st.error('Could not open this file as HDF5 (missing, corrupted, or invalid format).')
                return
            st.session_state['gui_h5_explore_catalog'] = cat
            st.session_state['gui_h5_explore_notes'] = notes

    schema = str(st.session_state.get('gui_h5_explore_schema') or '')
    cat = dict(st.session_state.get('gui_h5_explore_catalog') or {})
    notes = list(st.session_state.get('gui_h5_explore_notes') or [])

    if mode == 'custom':
        st.success(f'Loaded `{os.path.basename(loaded)}` — **custom hierarchy** mode.')
        with st.expander('File notes', expanded=False):
            st.caption('Use the **X** and **Y** columns to open different folders and pick one dataset each.')
        _render_h5_custom_hierarchy_dropdowns(loaded)
        _render_h5_attribute_editor(loaded)
        return

    if not cat:
        st.warning(
            'No auto-detected series for this file. Switch to **Custom — navigate folders + X/Y dropdowns** above.'
        )
        _render_h5_attribute_editor(loaded)
        return

    st.success(f'Loaded `{os.path.basename(loaded)}` — **{schema}** schema, **{len(cat)}** series.')

    with st.expander('File notes', expanded=False):
        if notes:
            for line in notes:
                st.caption(line)
        else:
            st.caption('No loader warnings.')

    if schema == 'v2':
        trials = list_v2_trials(loaded)
        if trials:
            if st.session_state.get('gui_h5_explore_trial') not in trials:
                st.session_state['gui_h5_explore_trial'] = trials[0]
            st.selectbox(
                'Trial group (`02_TimeSeries/…`)',
                options=trials,
                key='gui_h5_explore_trial',
                on_change=_cb_h5_explorer_trial_changed,
                help='Each group is one exported segment. Changing trial reloads available series.',
            )

    keys = sorted(cat.keys())
    if st.session_state.get('gui_h5_explore_x') not in keys:
        st.session_state['gui_h5_explore_x'] = keys[0]
    if st.session_state.get('gui_h5_explore_y') not in keys:
        for pref in (
            'Primary torque corrected (N·m)',
            'Tz (N·m)',
            'z Torque (N·m)',
            'x Torque (N·m)',
        ):
            if pref in keys:
                st.session_state['gui_h5_explore_y'] = pref
                break
        else:
            st.session_state['gui_h5_explore_y'] = keys[min(1, len(keys) - 1)]

    st.caption('Quick presets — click to set X and Y (axis dropdowns below).')
    pr1, pr2, pr3, pr4 = st.columns(4)
    with pr1:
        if st.button('Time → torque (primary)', key='gui_h5_preset_tt'):
            st.session_state['gui_h5_explore_x'] = 'Time (s)' if 'Time (s)' in keys else keys[0]
            if 'Primary torque corrected (N·m)' in keys:
                st.session_state['gui_h5_explore_y'] = 'Primary torque corrected (N·m)'
            elif 'Tz (N·m)' in keys:
                st.session_state['gui_h5_explore_y'] = 'Tz (N·m)'
            elif 'z Torque (N·m)' in keys:
                st.session_state['gui_h5_explore_y'] = 'z Torque (N·m)'
            st.rerun()
    with pr2:
        if st.button('Curvature → torque', key='gui_h5_preset_ct'):
            if 'Curvature κ (1/m) from angle' in keys:
                st.session_state['gui_h5_explore_x'] = 'Curvature κ (1/m) from angle'
            elif 'Angle measured (deg)' in keys:
                st.session_state['gui_h5_explore_x'] = 'Angle measured (deg)'
            elif 'Encoder / angle (deg)' in keys:
                st.session_state['gui_h5_explore_x'] = 'Encoder / angle (deg)'
            if 'Primary torque corrected (N·m)' in keys:
                st.session_state['gui_h5_explore_y'] = 'Primary torque corrected (N·m)'
            elif 'Tz (N·m)' in keys:
                st.session_state['gui_h5_explore_y'] = 'Tz (N·m)'
            elif 'z Torque (N·m)' in keys:
                st.session_state['gui_h5_explore_y'] = 'z Torque (N·m)'
            st.rerun()
    with pr3:
        if st.button('Time → Tz', key='gui_h5_preset_tz'):
            st.session_state['gui_h5_explore_x'] = 'Time (s)' if 'Time (s)' in keys else keys[0]
            if 'Tz (N·m)' in keys:
                st.session_state['gui_h5_explore_y'] = 'Tz (N·m)'
            elif 'z Torque (N·m)' in keys:
                st.session_state['gui_h5_explore_y'] = 'z Torque (N·m)'
            st.rerun()
    with pr4:
        if st.button(
            'Swap X ↔ Y',
            key='gui_h5_preset_swap_xy',
            help='Exchange the X and Y axis series.',
        ):
            sx = st.session_state.get('gui_h5_explore_x')
            sy = st.session_state.get('gui_h5_explore_y')
            st.session_state['gui_h5_explore_x'] = sy
            st.session_state['gui_h5_explore_y'] = sx
            st.rerun()

    c1, c2 = st.columns(2)
    with c1:
        x_key = st.selectbox('X axis', options=keys, key='gui_h5_explore_x')
    with c2:
        y_key = st.selectbox('Y axis', options=keys, key='gui_h5_explore_y')

    try:
        x_data, y_data, n = align_h5_catalog_xy(cat, x_key, y_key)
    except KeyError as e:
        st.error(str(e))
        return

    if n <= 0:
        st.warning('Selected series are empty or length mismatch produced zero samples.')
        return

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=x_data,
            y=y_data,
            mode='lines',
            name=y_key,
            line=dict(width=1.2),
        )
    )
    fig.update_layout(
        title=f'{y_key} vs {x_key} (n={n})',
        xaxis_title=x_key,
        yaxis_title=y_key,
        margin=dict(l=48, r=24, t=48, b=48),
        hovermode='x unified',
        height=520,
    )
    st.plotly_chart(fig, use_container_width=True)

    with st.expander('Series lengths', expanded=False):
        rows = [{'Series': k, 'Length': int(v.size)} for k, v in sorted(cat.items(), key=lambda kv: kv[0])]
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    _render_h5_attribute_editor(loaded)


def _trigger_emergency_stop() -> tuple[bool, str]:
    """Run NI-DAQ emergency stop and return `(ok, message)`."""
    dev = None
    if st.session_state.get('bender') is not None:
        dev = getattr(st.session_state['bender'], 'device_name', None)
    return daq_emergency_stop(dev)


def _render_app_chrome() -> None:
    st.markdown(
        '<div id="bnd-main-content" tabindex="-1"></div>'
        '<div class="bnd-workflow-active" aria-hidden="true"></div>',
        unsafe_allow_html=True,
    )
    if _nav_route() == 'sim_compare':
        st.markdown('<div class="bnd-sim-compare-active" aria-hidden="true"></div>', unsafe_allow_html=True)
    if _nav_route() == 'stepwise' and bool(st.session_state.get('gui_use_legacy_stepwise_chrome', False)):
        st.markdown('<div class="bnd-stepwise-active" aria-hidden="true"></div>', unsafe_allow_html=True)
        _inject_stepwise_compact_layout_css()
        # Wide-enough columns so "Home" stays one line; align row vertically.
        try:
            c_home, c_spacer, c_logo = st.columns([2.2, 7.6, 2.2], vertical_alignment='center')
        except TypeError:
            c_home, c_spacer, c_logo = st.columns([2.2, 7.6, 2.2])
        with c_home:
            if st.button('Home', key='gui_nav_home_stepwise', type='secondary'):
                _autosave_tick(force=True)
                st.session_state['gui_app_route'] = 'landing'
                st.rerun()
        with c_spacer:
            pass
        with c_logo:
            lg0, lg1 = st.columns([1, 1])
            with lg1:
                if os.path.isfile(_LOGO_PATH):
                    try:
                        st.image(_LOGO_PATH, width=40)
                    except Exception:
                        pass
        _inject_load_save_button_theme()
        _render_recovery_status_ui()
        return

    h0, h1 = st.columns([1.1, 8.9], vertical_alignment='center')
    with h0:
        if st.button('← Home', key='gui_nav_home_main'):
            _autosave_tick(force=True)
            st.session_state['gui_app_route'] = 'landing'
            st.rerun()
    with h1:
        st.markdown('### The CritterGripper App')
        _mode = _nav_route()
        if _mode == 'scratch':
            st.caption('Full workflow — all sections below are visible.')
        elif _mode == 'templates':
            st.caption('Template mode — checklist at top controls which sections appear.')
        elif _mode == 'sim_compare':
            st.caption('Simulation & Comparison — numpy cantilever mechanics only; no NI-DAQ.')
        elif _mode == 'h5_explorer':
            st.caption('Visualize Data — plot HDF5 series; standard or custom layout; optional attribute edits.')
    lg_l, lg_r = st.columns([1, 1])
    with lg_l:
        if os.path.isfile(_LOGO_PATH):
            try:
                st.image(_LOGO_PATH, width=120)
            except Exception:
                pass
    with lg_r:
        nsf = _nsf_logo_path()
        if nsf:
            st.markdown(
                f'<img src="{_img_data_uri(nsf)}" style="max-height:56px;width:auto;" alt="NSF"/>',
                unsafe_allow_html=True,
            )
    st.divider()
    _inject_load_save_button_theme()
    _render_recovery_status_ui()


def _render_sidebar() -> None:
    """Sidebar: unified workflow checklist + emergency + settings."""
    # #region agent log
    _agent_debug_log(
        hypothesis_id='F',
        location='bender_streamlit_gui.py:_render_sidebar_entry',
        message='pre_sidebar_state',
        data={
            'sess_fishmass': st.session_state.get('bio_fishmass'),
            'sess_dclamp': st.session_state.get('bio_dclamp'),
            'bender_fishmass': getattr(st.session_state.get('bender'), 'fishmass', None) if st.session_state.get('bender') is not None else None,
            'bender_dclamp': getattr(st.session_state.get('bender'), 'dclamp', None) if st.session_state.get('bender') is not None else None,
            'meas_conf': bool(st.session_state.get('gui_measurements_confirmed')),
            'invalidated': bool(st.session_state.get('gui_bio_apply_invalidated')),
            'route': _nav_route(),
            'step': _stepwise_step() if _nav_route() == 'stepwise' else None,
        },
    )
    # #endregion
    with st.sidebar:
        if _nav_route() == 'sim_compare':
            st.markdown('### Simulation & Comparison')
            st.caption('This path does not use NI-DAQ or the live experiment stack. Use **Settings** for display options.')
            _render_sidebar_settings_expander(leading_divider=False)
            return
        b = st.session_state.get('bender')
        checks = _collect_check_tuples(b) if b is not None else []
        checks_by_sec: dict[str, list[str]] = {}
        for sec, msg in checks:
            checks_by_sec.setdefault(sec, []).append(msg)
        st.markdown('**Workflow checklist**')
        _session_src = str(st.session_state.get('gui_session_source') or 'fresh')
        st.caption(f"Session: {'Restored' if _session_src == 'restored' else 'Fresh start'}")
        tt_chk = str(st.session_state.get('test_type_select') or getattr(b, 'test_type', '') or 'dynamic')
        _ready = _workflow_ready_state(b, tt_chk)
        _setup_ok = _ready['setup_ok']
        _meas_ok = _ready['measurements_ok']
        _proto_ok = _ready['protocol_ok']
        _review_used = bool(st.session_state.get('gui_review_data_used'))
        st.caption(f"1. Setup {'✅' if _setup_ok else '⚠️'}")
        st.caption(f"2. Measurements {'✅' if _meas_ok else '⚠️'}")
        st.caption(f"3. Protocol/Run {'✅' if _proto_ok else '⚠️'}")
        st.caption(f"4. Review Data {'✅' if _review_used else '⚪'}")
        st.session_state.setdefault('gui_fix_panel_seen', False)
        if st.button('Click here to fix warning', key='gui_fix_warning_open', use_container_width=True, type='secondary'):
            st.session_state['gui_fix_panel_open'] = True
        _show_fix_panel = bool(st.session_state.get('gui_fix_panel_open'))
        if _show_fix_panel:
            fix_lines = _build_checklist_fix_lines(
                b=b,
                setup_ok=_setup_ok,
                measurements_ok=_meas_ok,
                protocol_ok=_proto_ok,
                checks_by_sec=checks_by_sec,
                tt=tt_chk,
            )
            with st.container(border=True):
                if fix_lines:
                    for row in fix_lines:
                        st.caption(f'- {row}')
                else:
                    st.caption('All required checklist steps are complete.')
                if st.button('Hide fixes', key='gui_fix_warning_hide', use_container_width=True):
                    st.session_state['gui_fix_panel_open'] = False
                    st.session_state['gui_fix_panel_seen'] = True
                    st.rerun()
            st.session_state['gui_fix_panel_seen'] = True
        with st.expander('Advanced diagnostics', expanded=False):
            if checks_by_sec:
                for sec, rows in checks_by_sec.items():
                    st.caption(f'{sec}')
                    for row in rows:
                        st.caption(f'- {row}')
            else:
                st.caption('No diagnostics.')
        st.divider()
        # No simulation-mode toggle on NI-DAQ experiment workflows (real hardware path only).
        _route = _nav_route()
        if _route in ('stepwise', 'scratch', 'templates'):
            st.session_state['gui_simulation_mode'] = False
            _b = st.session_state.get('bender')
            if _b is not None:
                _b.simulation_mode = False
        else:
            _render_simulation_sidebar()
        with st.container(border=True):
            st.markdown('### Emergency stop')
            if st.button(
                'Stop DAQ & reset NI device',
                key='gui_kill_daq',
                type='primary',
                use_container_width=True,
            ):
                ok, msg = _trigger_emergency_stop()
                if ok:
                    st.success(msg)
                else:
                    st.warning(msg)
        _render_sidebar_settings_expander(leading_divider=True)


_STEPWISE_TAB_LABELS = (
    '1 · Setup (hardware + data path)',
    '2 · Measurements',
    '3 · Protocol & run',
    '4 · Visualize & notes',
    '5 · Review',
)

# Short labels for tab buttons (keep narrow columns readable)
_STEPWISE_TAB_SHORT = (
    'Setup',
    'Measurements',
    'Protocol/Run',
    'Visualize',
    'Review',
)


def _render_stepwise_rail() -> None:
    st.session_state.setdefault('gui_stepwise_step', 0)
    cur = int(st.session_state.get('gui_stepwise_step', 0) or 0)
    cur = max(0, min(4, cur))
    if cur != st.session_state.get('gui_stepwise_step'):
        st.session_state['gui_stepwise_step'] = cur

    if err := st.session_state.pop('gui_sw_default_hw_err', None):
        _st_error_detail(
            'Default hardware load failed.',
            ['Use step 1 to pick module', 'Read Details'],
            err,
        )

    with st.container(border=True):
        st.progress(min(1.0, (cur + 1) / 5.0), text=f'Step {cur + 1} / 5 · {_STEPWISE_TAB_SHORT[cur]}')

        st.markdown('<div class="bnd-stepwise-tab-marker" aria-hidden="true"></div>', unsafe_allow_html=True)
        tab_cols = st.columns(5, gap='small')
        for i in range(5):
            with tab_cols[i]:
                _lbl = f'{i + 1}. {_STEPWISE_TAB_SHORT[i]}'
                if st.button(
                    _lbl,
                    key=f'gui_sw_tab_{i}',
                    use_container_width=True,
                    type='primary' if i == cur else 'secondary',
                    disabled=i == cur,
                    help=_STEPWISE_TAB_LABELS[i],
                ):
                    st.session_state['gui_stepwise_step'] = i
                    st.rerun()

        n1, n2, _ = st.columns([1, 1, 4])
        with n1:
            if st.button(
                '← Previous',
                key='gui_sw_back',
                disabled=cur <= 0,
                use_container_width=True,
                help='Move to the previous step',
            ):
                st.session_state['gui_stepwise_step'] = cur - 1
                st.rerun()
        with n2:
            if st.button(
                'Next →',
                key='gui_sw_next',
                disabled=cur >= 4,
                use_container_width=True,
                help='Move to the next step',
            ):
                st.session_state['gui_stepwise_step'] = cur + 1
                st.rerun()

    # Prerequisite reminders are centralized in the sidebar checklist + fix panel.


def _cb_tpl_config_module_changed() -> None:
    st.session_state['gui_tpl_reload_config'] = True


def _pick_folder_with_dialog(initial_dir: str) -> tuple[Optional[str], Optional[str]]:
    """Return (folder_path, error_text) from native folder picker."""
    try:
        import tkinter as tk
        from tkinter import filedialog

        root = tk.Tk()
        root.withdraw()
        root.attributes('-topmost', True)
        picked = filedialog.askdirectory(initialdir=initial_dir or _ROOT, title='Select folder')
        root.destroy()
        if not picked:
            return None, None
        return os.path.normpath(str(picked)), None
    except Exception as e:
        return None, f'{type(e).__name__}: {e}'


_DATA_FOLDER_DD_SENTINEL = '— Select folder —'


def _immediate_subdirs(root: str, *, max_entries: int = 50) -> list[str]:
    """List immediate child directories under ``root`` (stable, browser-safe)."""
    root = os.path.normpath(str(root or '').strip())
    if not root or not os.path.isdir(root):
        return []
    out: list[str] = []
    try:
        for name in sorted(os.listdir(root), key=str.casefold):
            p = os.path.join(root, name)
            if os.path.isdir(p):
                out.append(os.path.normpath(p))
            if len(out) >= max_entries:
                break
    except OSError:
        return []
    return out


def _data_folder_dropdown_choice_list(current_folder: str) -> list[str]:
    """
    Candidate data folders for a Streamlit dropdown (no tkinter).

    Includes the project root, common subfolders, and one level of child directories.
    If ``current_folder`` is a valid directory, it is included even when not under _ROOT.
    """
    cur = os.path.normpath(str(current_folder or '').strip()) if str(current_folder or '').strip() else ''
    seeds = [
        _ROOT,
        os.path.join(_ROOT, 'TestData'),
        os.path.join(_ROOT, 'SessionSnapshots'),
        os.path.join(_ROOT, 'templates', 'protocols'),
        os.path.join(_ROOT, 'templates', 'specimens'),
        os.path.join(_ROOT, 'templates', 'configs'),
    ]
    seen: set[str] = set()
    ordered: list[str] = []
    for seed in seeds:
        norm = os.path.normpath(seed)
        if not norm or norm in seen or not os.path.isdir(norm):
            continue
        seen.add(norm)
        ordered.append(norm)
        for sub in _immediate_subdirs(norm, max_entries=40):
            if sub not in seen:
                seen.add(sub)
                ordered.append(sub)
    if cur and os.path.isdir(cur) and cur not in seen:
        ordered.insert(0, cur)
    return ordered


def _render_data_folder_dropdown(*, key_suffix: str) -> None:
    """
    Browser-native folder entry: text input accepts pasted folder paths.
    """
    _ = key_suffix
    if 'gui_data_folder' not in st.session_state:
        st.session_state['gui_data_folder'] = ''

    def _on_folder_text_change() -> None:
        raw = str(st.session_state.get('gui_data_folder') or '').strip()
        if not raw:
            return
        picked = os.path.normpath(raw)
        _store_selected_data_folder(picked)

    folder_path = st.text_input(
        'Data folder',
        key='gui_data_folder',
        on_change=_on_folder_text_change,
        help=(
            'Paste a folder path. The full HDF5 path is **folder + file name** '
            'in the next column. Native **Browse…** may not work on remote desktops or hosted Streamlit.'
        ),
    )
    folder_path = str(folder_path or '').strip()
    if folder_path:
        if os.path.isdir(folder_path):
            st.success('✅ Folder found')
        else:
            st.error('❌ Folder not found — check path')


def _pick_file_with_dialog(initial_dir: str, *, title: str, filetypes: list[tuple[str, str]]) -> tuple[Optional[str], Optional[str]]:
    """Return (file_path, error_text) from native open-file picker."""
    try:
        import tkinter as tk
        from tkinter import filedialog

        root = tk.Tk()
        root.withdraw()
        root.attributes('-topmost', True)
        picked = filedialog.askopenfilename(
            initialdir=initial_dir or _ROOT,
            title=title,
            filetypes=filetypes,
        )
        root.destroy()
        if not picked:
            return None, None
        return os.path.normpath(str(picked)), None
    except Exception as e:
        return None, f'{type(e).__name__}: {e}'


def _store_selected_data_folder(picked_dir: str) -> None:
    """Persist selected folder and keep a usable filename for composed output paths."""
    folder = os.path.normpath(str(picked_dir or '').strip())
    if not folder:
        return
    # User-picked folder is authoritative; discard deferred config-load path seeds.
    st.session_state.pop('gui_pending_data_folder', None)
    st.session_state.pop('gui_pending_data_filename', None)
    st.session_state['gui_data_folder'] = folder
    cur_fn = str(st.session_state.get('gui_data_filename') or '').strip()
    if cur_fn:
        return
    # If filename is blank, preserve the current experiment basename so saves move to selected folder.
    b = st.session_state.get('bender')
    fallback = str(getattr(b, 'outputfile', '') or '').strip() if b is not None else ''
    if fallback:
        base = os.path.basename(os.path.normpath(fallback))
        if base:
            st.session_state['gui_data_filename'] = base


def _render_config_module_navigator(*, key_prefix: str, label: str = 'Hardware configuration module') -> None:
    """Hardware config picker: paste a `.py` path, confirm, then click Load."""
    _sel = _normalize_config_module_name(str(st.session_state.get('gui_load_cfg_select') or ''))
    if _sel:
        _sel_path = os.path.join(_ROOT, _sel.replace('.', os.sep) + '.py')
    else:
        _sel_path = ''
    if 'gui_load_cfg_file_path' not in st.session_state:
        st.session_state['gui_load_cfg_file_path'] = _sel_path if _sel and os.path.isfile(_sel_path) else ''

    _path_col, _load_col = st.columns([5, 1])
    with _path_col:
        cfg_path = str(
            st.text_input(
                label,
                key='gui_load_cfg_file_path',
                placeholder='Paste full path to a hardware config .py file',
                help='Paste any `.py` hardware config path, then click **Load**.',
            )
            or ''
        ).strip()
    with _load_col:
        load_clicked = st.button(
            'Load',
            key=f'gui_btn_load_hw_cfg_{key_prefix}',
            use_container_width=True,
            type='primary',
            help='Import the pasted hardware configuration module into the app.',
        )

    _eff_mod = None
    if cfg_path:
        norm_cfg_path = os.path.normpath(cfg_path)
        if os.path.isfile(norm_cfg_path):
            st.success('✅ File found')
            try:
                _eff_mod, _ = _resolve_hardware_config_import_target(norm_cfg_path)
            except ValueError:
                pass
        else:
            st.error('❌ File not found — check path')

    if load_clicked:
        if not cfg_path:
            st.warning('Paste a hardware config `.py` path first.')
        elif not os.path.isfile(os.path.normpath(cfg_path)):
            st.error('❌ File not found — check path')
        else:
            norm_cfg_path = os.path.normpath(cfg_path)
            try:
                eff, _ = _resolve_hardware_config_import_target(norm_cfg_path)
            except ValueError as e:
                st.error(str(e))
                eff = ''
            if eff:
                err = _apply_loaded_config_module(norm_cfg_path)
                if err:
                    _st_error_detail(
                        'Hardware config load failed.',
                        ['Check module path', 'Read Details'],
                        err,
                    )
                else:
                    st.session_state['gui_load_cfg_select'] = eff
                    st.session_state['gui_pending_load_cfg_file_path'] = norm_cfg_path
                    st.success(f'Loaded `{eff}`')
                    st.rerun()

    _display_mod = _eff_mod or _normalize_config_module_name(str(st.session_state.get('gui_load_cfg_select') or ''))
    st.caption(f'{label}: `{_display_mod or "(none selected)"}`')


def _render_template_procedure_strip() -> None:
    """Compact loaders: config + biometrics from disk, then experiment sections (procedure edited in the app)."""
    _ensure_hw_config_session_defaults()
    st.subheader('Load saved files')
    st.caption('Order: hardware `.py` → data path → biometrics file. Reload module after changing the selection.')

    if st.session_state.pop('gui_tpl_reload_config', False):
        eff = _normalize_config_module_name(str(st.session_state.get('gui_load_cfg_select') or ''))
        if not eff:
            st.session_state['gui_tpl_config_load_err'] = 'Choose a config module with Browse.'
        else:
            err = _apply_loaded_config_module(_raw_mod_for_hardware_config_load(module_stem=eff))
            if err:
                st.session_state['gui_tpl_config_load_err'] = err
            else:
                st.rerun()
    if err_msg := st.session_state.pop('gui_tpl_config_load_err', None):
        _st_error_detail(
            'Hardware config load failed.',
            ['Check module name', 'Fix errors in Details'],
            err_msg,
        )

    st.subheader('1 · Hardware configuration')
    with st.container(border=True):
        _cfg_mods = discover_config_modules(_ROOT)
        _default_pick = str(st.session_state.get('cfg_mod') or '')
        if _default_pick not in _cfg_mods:
            _default_pick = _cfg_mods[0]
        st.session_state.setdefault('gui_load_cfg_select', _default_pick)
        _render_config_module_navigator(key_prefix='tpl', label='Hardware configuration module')
        if _load_save_button(
            'Load hardware configuration',
            key='gui_btn_load_config_tpl',
            help='Reload the selected hardware `.py` module from disk.',
        ):
            st.session_state['gui_tpl_reload_config'] = True
            st.rerun()
        if st.session_state.get('bender') is not None:
            st.success(f"Hardware configuration active: `{getattr(st.session_state['bender'], 'config_name', '?')}`")

    st.divider()
    st.subheader('2 · Data file path')
    with st.container(border=True):
        _dfc, _fnc = st.columns(2)
        with _dfc:
            _render_data_folder_dropdown(key_suffix='tpl')
            _preview_out = _compose_output_h5_path().strip()
            _preview_folder = str(st.session_state.get('gui_data_folder') or '').strip()
            if _preview_out:
                st.caption(f'Save path: `{_preview_out}`')
            elif _preview_folder:
                st.caption(f'Save path: `{_preview_folder}`')
        with _fnc:
            st.text_input(
                ' ',
                key='gui_data_filename',
                placeholder='datafilename.h5',
                help=DATA_FILE_NAME_HELP,
                label_visibility='collapsed',
            )

    st.divider()
    _bio_tpl_dir = _shared_experiment_dir()
    _bio_tpl_list = list_biometrics_template_files(_bio_tpl_dir)
    _bio_opts: list = [None] + _bio_tpl_list
    if 'gui_biometrics_template_select' not in st.session_state:
        st.session_state['gui_biometrics_template_select'] = None
    st.selectbox(
        'Biometrics file',
        _bio_opts,
        format_func=_biometrics_template_option_label,
        key='gui_biometrics_template_select',
        help='Saved specimen/setup files in your data folder (`.json` on disk).',
    )
    if _load_save_button('Load biometrics into form', key='gui_biometrics_btn_load_tpl'):
        _bp = st.session_state.get('gui_biometrics_template_select')
        if not _bp:
            st.session_state['gui_biometrics_load_feedback'] = (False, 'Choose a biometrics file first.')
        else:
            st.session_state['gui_pending_biometrics_path'] = _bp
        st.rerun()
    if bf := st.session_state.pop('gui_biometrics_load_feedback', None):
        ok_bf, txt_bf = bf
        if ok_bf:
            st.success(txt_bf)
        elif str(txt_bf).strip() == 'Choose a biometrics file first.':
            st.warning(txt_bf)
        else:
            _st_error_detail(
                'Biometrics load failed.',
                ['Check file format', 'Read Details'],
                txt_bf,
            )
    if st.session_state.get('bender') is not None:
        if _load_save_button(
            'Apply loaded biometrics to experiment',
            key='gui_tpl_apply_bio_to_bender',
            help='Same as **section 3 · Apply all biometrics**: identity, intrinsic, experimental conditions, clamp, profile, inertial flag.',
        ):
            _apply_all_biometrics_to_bender(st.session_state['bender'])
            st.session_state['gui_tpl_bio_done'] = True
            st.success('Biometrics applied. Procedure sections are below.')
            st.rerun()


def _template_procedure_gate() -> None:
    if _nav_route() != 'templates' or not _tpl_only_procedure():
        return
    _render_template_procedure_strip()
    if not st.session_state.get('bender'):
        st.info('Load **hardware configuration** above to continue.')
        st.stop()
    if not st.session_state.get('gui_tpl_bio_done'):
        st.info('Load a **biometrics** file and click **Apply loaded biometrics to experiment**.')
        st.stop()


def main():
    st.set_page_config(
        page_title='CritterGripper',
        layout='wide',
        initial_sidebar_state='expanded',
        menu_items={
            'Get Help': None,
            'Report a bug': None,
            'About': (
                'CritterGripper is the browser-based console for Bender bending experiments. It walks operators through '
                'hardware configuration, where HDF5 files will be written, specimen biometrics, and protocol parameters; '
                'offers a no-DAQ preview of commanded motion; validates required fields; and runs acquisition and export '
                'when you deliberately enable it. A separate H5 data explorer plots saved trials without loading the full '
                'workflow.\n\n'
                'Reliable biomechanics depends on consistent metadata, a clear picture of what the rig will do before it '
                'moves, and separation between “configure / preview” and “run.” This app supports that workflow, surfaces '
                'common setup gaps (paths, biometrics, calibration), and exposes an NI-DAQ hardware reset for situations '
                'where outputs should halt immediately — alongside, not instead of, your lab’s normal safety practices.'
            ),
        },
    )
    # #region agent log
    _agent_debug_log(
        hypothesis_id='F',
        location='bender_streamlit_gui.py:main_start',
        message='entry_state',
        data={
            'sess_fishmass': st.session_state.get('bio_fishmass'),
            'sess_dclamp': st.session_state.get('bio_dclamp'),
            'sess_xsec': st.session_state.get('bio_xsec'),
            'fishmass_in_keys': 'bio_fishmass' in st.session_state,
            'dclamp_in_keys': 'bio_dclamp' in st.session_state,
            'bender_fishmass': getattr(st.session_state.get('bender'), 'fishmass', None) if st.session_state.get('bender') is not None else None,
            'bender_dclamp': getattr(st.session_state.get('bender'), 'dclamp', None) if st.session_state.get('bender') is not None else None,
            'applied_sig_fm': st.session_state.get('gui_bio_applied_sig')[3] if isinstance(st.session_state.get('gui_bio_applied_sig'), tuple) and len(st.session_state.get('gui_bio_applied_sig')) > 3 else None,
            'applied_sig_dc': st.session_state.get('gui_bio_applied_sig')[9] if isinstance(st.session_state.get('gui_bio_applied_sig'), tuple) and len(st.session_state.get('gui_bio_applied_sig')) > 9 else None,
            'invalidated': bool(st.session_state.get('gui_bio_apply_invalidated')),
            'route': str(st.session_state.get('gui_app_route') or 'landing'),
            'step': int(st.session_state.get('gui_stepwise_step', 0) or 0),
        },
    )
    # #endregion
    st.session_state.setdefault('gui_ui_theme', GUI_UI_THEME_OPTIONS[0])
    _migrate_gui_ui_theme_session()
    st.session_state.setdefault('gui_ui_large_text', False)
    st.session_state.setdefault('gui_app_route', 'landing')
    st.session_state.setdefault('gui_session_source', 'fresh')
    st.session_state.setdefault('gui_setup_confirmed', False)
    st.session_state.setdefault('gui_measurements_confirmed', False)
    st.session_state.setdefault('gui_protocol_confirmed', False)
    _bootstrap_autosave_recovery()
    _inject_accessibility_theme()
    _ensure_gui_data_path_session_keys()
    if 'gui_post_notes' not in st.session_state:
        st.session_state['gui_post_notes'] = ''
    if 'gui_genus_species' not in st.session_state:
        st.session_state['gui_genus_species'] = ''
    if 'gui_specimen_id' not in st.session_state:
        st.session_state['gui_specimen_id'] = ''
    if 'bio_prep_condition' not in st.session_state:
        st.session_state['bio_prep_condition'] = ''
    st.session_state.setdefault('bio_prof_rho_preset', 'Custom — edit the number below')
    if 'gui_default_state_baseline' not in st.session_state and 'gui_recovered_state_baseline' not in st.session_state:
        st.session_state['gui_default_state_baseline'] = _collect_persistable_state()
    _update_state_origin_summary()
    if _nav_route() == 'landing':
        _render_landing_page()
        _autosave_tick()
        return

    _flush_pending_load_config_session()
    _consume_pending_biometrics_template()
    _refresh_confirmation_flags()
    _sanitize_stale_run_state()
    # Repair data-path fields from persisted signatures before any widget binds to
    # gui_data_folder / gui_data_filename (Streamlit forbids mutating those keys later).
    _repair_data_path_fields_from_session()

    _render_app_chrome()
    _render_sidebar()

    if _nav_route() == 'h5_explorer':
        _render_h5_explorer()
        _autosave_tick()
        return

    if _nav_route() == 'sim_compare':
        _render_simulation_comparison_page()
        _autosave_tick()
        return

    if _nav_route() == 'templates':
        st.session_state.pop('gui_tpl_need_procedure', None)
        st.markdown('**Template mode** — check what you already have on disk.')
        if 'gui_tpl_chk_config' not in st.session_state:
            st.session_state['gui_tpl_chk_config'] = False
        st.checkbox('I have a saved hardware **config**', value=True, key='gui_tpl_chk_config')
        if 'gui_tpl_chk_biometrics' not in st.session_state:
            st.session_state['gui_tpl_chk_biometrics'] = False
        st.checkbox(
            'I have a saved **biometrics** file',
            value=True,
            key='gui_tpl_chk_biometrics',
            help='Usually a `.json` file in your data folder with lengths, clamp spacing, etc.',
        )
        if 'gui_tpl_have_protocol_template' not in st.session_state:
            st.session_state['gui_tpl_have_protocol_template'] = False
        st.checkbox(
            'I already have a **protocol** template',
            value=False,
            key='gui_tpl_have_protocol_template',
            help='A saved procedure file (`.json`) listing experiment type and parameters.',
        )
        if not _tpl_only_procedure():
            st.info(
                'All sections below match **Build from scratch**, or use **Load template into form** in **section 4** when you '
                'already have a protocol file. Uncheck the third box above if you only need to configure the procedure after '
                'loading config and biometrics.'
            )

    st.caption('Sections follow the numbered order on the page.')

    _cfg_mods = discover_config_modules(_ROOT)
    _template_procedure_gate()

    if _show_hw_config_section() or _show_data_path_section():
        st.session_state.setdefault('gui_setup_actions_show', False)
        if st.button(
            'Start fresh / save progress',
            key='gui_setup_actions_btn',
            type='secondary',
            use_container_width=True,
        ):
            st.session_state['gui_setup_actions_show'] = not bool(st.session_state.get('gui_setup_actions_show', False))
        if st.session_state.get('gui_setup_actions_show'):
            st.caption('Choose what to do with your current form state before starting over.')
            _a_save, _a_home, _a_tpl = st.columns(3, gap='small')
            with _a_save:
                if st.button('Save progress snapshot', key='gui_save_progress_snapshot', use_container_width=True):
                    ok, msg = _save_progress_snapshot()
                    if ok:
                        st.success(f'Snapshot saved: `{msg}`')
                    else:
                        _st_error_detail(
                            'Could not save progress snapshot.',
                            ['Check write permissions', 'Verify workspace is writable'],
                            msg,
                        )
                _snap_files = _list_manual_snapshot_files(max_entries=40)
                st.selectbox(
                    'Recent snapshots',
                    options=[''] + _snap_files,
                    key='gui_snapshot_file_pick',
                    format_func=lambda p: '— Select saved snapshot —' if not p else os.path.basename(str(p)),
                    help='Manual snapshots saved in SessionSnapshots.',
                )
                st.text_input(
                    'Snapshot file path',
                    key='gui_manual_snapshot_path',
                    placeholder='Paste full path to session_snapshot_*.json',
                    help='Use a recent snapshot above or paste a full path.',
                )
                if _load_save_button('Load snapshot', key='gui_btn_load_snapshot', button_type='secondary'):
                    _picked = str(st.session_state.get('gui_snapshot_file_pick') or '').strip()
                    _typed = str(st.session_state.get('gui_manual_snapshot_path') or '').strip()
                    _path = _typed or _picked
                    ok, msg = _load_manual_snapshot(_path)
                    if ok:
                        st.success(f'Snapshot loaded: `{msg}`')
                        st.info('Re-apply setup, biometrics, and procedure before running hardware.')
                        st.rerun()
                    elif msg == 'Choose a snapshot file first.':
                        st.warning(msg)
                    else:
                        _st_error_detail(
                            'Could not load snapshot.',
                            ['Check file path', 'Verify JSON schema and state object'],
                            msg,
                        )
            with _a_home:
                if st.button(
                    'Start fresh (discard & go to Home)',
                    key='gui_go_home_discard',
                    type='primary',
                    use_container_width=True,
                    help='Discard current unsaved workflow state and return to landing page.',
                ):
                    _autosave_tick(force=True)
                    _reset_workflow_session_to_home(clear_autosave=True)
                    st.rerun()
            with _a_tpl:
                if st.button(
                    'Go to Templates',
                    key='gui_go_templates_from_setup',
                    use_container_width=True,
                    help='Use template workflow instead of continuing with current unsaved setup.',
                ):
                    _autosave_tick(force=True)
                    if _attempt_route_switch_with_bio_warning(target_route='templates', origin='setup'):
                        st.rerun()
            _render_pending_bio_nav_warning('setup')

    _show_hw = _show_hw_config_section()
    _show_data = _show_data_path_section()
    _setup_left = _setup_right = None
    if _show_hw and _show_data:
        _setup_left, _setup_right = st.columns(2, gap='large')

    def _apply_setup_action(*, sw_dp: bool) -> None:
        _cfg_mode = str(st.session_state.get('gui_config_setup_mode', 'Load existing'))
        if _cfg_mode == 'Load existing':
            _ensure_hw_config_session_defaults()
            eff = _normalize_config_module_name(str(st.session_state.get('gui_load_cfg_select') or ''))
            if not eff:
                _mods_here = discover_config_modules(_ROOT)
                if not _mods_here:
                    _st_error_actions(
                        'No hardware config modules found.',
                        ['Add a `.py` config module in the project folder'],
                    )
                else:
                    _st_error_actions(
                        'No config module resolved.',
                        ['Use Browse in Step 1 to select module'],
                    )
                return
            b0 = st.session_state.get('bender')
            need_hw_reload = b0 is None or not _selected_config_matches_bender(b0, eff)
            err = (
                _apply_loaded_config_module(_raw_mod_for_hardware_config_load(module_stem=eff))
                if need_hw_reload
                else None
            )
            if err:
                _st_error_detail(
                    'Hardware config load failed.',
                    ['Check module name', 'Fix errors in Details'],
                    err,
                )
                return
            perr = _sec1_apply_composed_path_to_bender()
            if perr:
                if need_hw_reload:
                    _st_error_detail(
                        'Hardware OK; path failed.',
                        ['Fix folder name', 'Fix file name'],
                        perr,
                    )
                else:
                    _st_error_detail(
                        'Data path not applied.',
                        ['Fix folder name', 'Fix file name'],
                        perr,
                    )
                return
            if sw_dp and not need_hw_reload:
                st.toast('Data file path set.')
            elif sw_dp and need_hw_reload:
                st.toast('Hardware configuration loaded and data file path set.')
            else:
                st.toast('Hardware configuration and data file path applied.')
            if need_hw_reload:
                st.success(f'Loaded `{_normalize_config_module_name(eff)}`')
            else:
                st.success('Data file path set on the experiment object.')
            st.session_state['gui_setup_confirmed'] = True
            st.rerun()
            return

        if st.session_state.get('bender') is None:
            _st_error_actions('Setup apply blocked.', ['Load or build hardware config first'])
            return
        perr = _sec1_apply_composed_path_to_bender()
        if perr:
            _st_error_detail(
                'Data path not applied.',
                ['Fix folder name', 'Fix file name'],
                perr,
            )
        else:
            st.toast('Data file path set.' if sw_dp else 'Data file path applied.')
            st.session_state['gui_setup_confirmed'] = True
            st.rerun()

    if _show_hw:
        _hw_host = _setup_left if _setup_left is not None else st
        _hw_host.subheader('1 · Hardware configuration')
        with _hw_host.container(border=True):
            _default_pick = str(st.session_state.get('cfg_mod') or '')
            if _default_pick not in _cfg_mods:
                _default_pick = _cfg_mods[0]
            st.session_state.setdefault('gui_load_cfg_select', _default_pick)
            st.session_state.setdefault('gui_cfg_build_base', _cfg_mods[0])
            st.session_state.setdefault('gui_cfg_build_out', '')
            st.session_state.setdefault('gui_cfg_build_overwrite', False)

            if _template_hide_config_build_new():
                st.caption('Template mode: load-only here. Uncheck "I have a saved config" above to use **Build new**.')
                st.session_state['gui_config_setup_mode'] = 'Load existing'
                mode = 'Load existing'
            else:
                mode = _hardware_configuration_mode_toggle()

            if mode == 'Load existing':
                _render_config_module_navigator(key_prefix='main', label='Hardware configuration module')
            else:
                _flush_pending_cfg_build_base()
                _maybe_seed_cfg_build_fields()
                c_top_l, c_top_r = st.columns(2, gap='large')
                with c_top_l:
                    if 'gui_cfg_build_base_path' not in st.session_state:
                        _base_mod = str(st.session_state.get('gui_cfg_build_base') or '')
                        _base_path = os.path.join(_ROOT, _base_mod.replace('.', os.sep) + '.py') if _base_mod else ''
                        st.session_state['gui_cfg_build_base_path'] = _base_path if _base_path and os.path.isfile(_base_path) else ''
                    _base_cfg_path = str(
                        st.text_input(
                            'Start from template (base module for import *)',
                            key='gui_cfg_build_base_path',
                            placeholder='Paste full path to a base config .py file',
                            help='Other settings from the template stay unless you override them below.',
                        )
                        or ''
                    ).strip()
                    _resolved_base = None
                    if _base_cfg_path:
                        _base_cfg_norm = os.path.normpath(_base_cfg_path)
                        if os.path.isfile(_base_cfg_norm):
                            st.success('✅ File found')
                            try:
                                _resolved_base, _ = _resolve_hardware_config_import_target(_base_cfg_norm)
                            except ValueError:
                                pass
                        else:
                            st.error('❌ File not found — check path')
                    if _resolved_base and _resolved_base != str(
                        st.session_state.get('gui_cfg_build_base') or ''
                    ).strip():
                        st.session_state['gui_pending_cfg_build_base'] = _resolved_base
                        st.rerun()
                with c_top_r:
                    st.text_input(
                        'Save new config as (module name, no `.py`)',
                        key='gui_cfg_build_out',
                        placeholder='e.g. lab_setup_2026',
                        help='Writes a new `.py` file in this folder and loads it.',
                    )
                c_cfg_l, c_cfg_r = st.columns(2, gap='large')
                with c_cfg_l:
                    with st.expander('Calibration, direction & axis labels', expanded=False):
                        st.text_input('Force/torque calibration file', key='gui_cfg_bld_forcetorque_calibration_file')
                        st.text_input('Positive motor direction (`left` / `right`)', key='gui_cfg_bld_positive_motor_direction')
                        st.number_input(
                            'Specimen lateral index on positive motor side',
                            key='gui_cfg_bld_specimen_lateral_index',
                            step=1,
                            format='%d',
                        )
                        st.text_input('Motor axis (config `motor_axis`)', key='gui_cfg_bld_motor_axis')
                        st.text_input('Bending axis — sensor (`bending_axis_sensor`)', key='gui_cfg_bld_bending_axis_sensor')
                        st.text_input('Primary bending axis (e.g. zTorque)', key='gui_cfg_bld_primary_bending_axis')
                        st.text_input(
                            'Bending axis — specimen (`bending_axis_specimen`: dorsoventral / lateral / anteroposterior)',
                            key='gui_cfg_bld_bending_axis_specimen',
                        )
                        st.text_input('S1 side label', key='gui_cfg_bld_S1side')
                        st.text_input('S2 side label', key='gui_cfg_bld_S2side')
                    with st.expander('Motor, encoder & stim channels', expanded=False):
                        st.number_input('Motor full steps per rev', key='gui_cfg_bld_motor_steps', min_value=1, step=1, format='%d')
                        st.number_input('Motor gear ratio', key='gui_cfg_bld_motor_gear', min_value=1, step=1, format='%d')
                        st.number_input('Encoder pulses per rev', key='gui_cfg_bld_encoder_ppr', min_value=1, step=1, format='%d')
                        st.text_input('Motor port', key='gui_cfg_bld_motor_port')
                        st.text_input('Encoder counter channel', key='gui_cfg_bld_encoder_chan')
                        st.text_input('Stim AO channels (comma-separated)', key='gui_cfg_bld_stim_channels')
                    with st.expander('Sonomicrometry', expanded=False):
                        st.checkbox('Use sonomicrometry', key='gui_cfg_bld_use_sono')
                        st.text_input('Sono AI channels (comma-separated)', key='gui_cfg_bld_sono_channel')
                        st.text_input('Sono names (comma-separated)', key='gui_cfg_bld_sono_name')
                        st.number_input(
                            'Sono internal sample frequency (`sono_internal_samplefreq`)',
                            key='gui_cfg_bld_sono_internal_samplefreq',
                            min_value=1,
                            step=1,
                            format='%d',
                        )
                        st.text_input(
                            'Sono cal left [V_lo, V_hi, mm_lo, mm_hi] comma-separated',
                            key='gui_cfg_bld_sono_cal_left',
                            help='Four numbers, same order as the config file.',
                        )
                        st.text_input(
                            'Sono cal right [V_lo, V_hi, mm_lo, mm_hi] comma-separated',
                            key='gui_cfg_bld_sono_cal_right',
                        )
                with c_cfg_r:
                    with st.expander('DAQ rates & device', expanded=True):
                        st.text_input('NI-DAQ device name', key='gui_cfg_bld_device_name')
                        st.number_input('AI + encoder sample rate (Hz)', key='gui_cfg_bld_daq_ai_sr', min_value=1.0, format='%.f')
                        st.number_input('AO + DO sample rate (Hz)', key='gui_cfg_bld_daq_ao_sr', min_value=1.0, format='%.f')
                    with st.expander('Strain / ATI channels', expanded=False):
                        st.text_input('SG AI channels (comma-separated)', key='gui_cfg_bld_SG_chan')
                        st.text_input('SG channel names (comma-separated)', key='gui_cfg_bld_SG_name')
                    with st.expander('Stim monitor AI (optional)', expanded=False):
                        st.text_input(
                            'Stim monitor AI channels (comma-separated; empty = none)',
                            key='gui_cfg_bld_stim_monitor_chan',
                        )
                        st.text_input('Stim monitor channel names (comma-separated)', key='gui_cfg_bld_stim_monitor_name')
                    with st.expander('Default timing, ramp & motion', expanded=False):
                        st.number_input('waitbefore (s)', key='gui_cfg_bld_waitbefore', min_value=0.0, format='%.6g')
                        st.number_input('waitafter (s)', key='gui_cfg_bld_waitafter', min_value=0.0, format='%.6g')
                        st.number_input('rampdur (s)', key='gui_cfg_bld_rampdur', min_value=0.0, format='%.6g')
                        st.number_input('prepoststim_dur (s)', key='gui_cfg_bld_prepoststim_dur', min_value=0.0, format='%.6g')
                        st.number_input('prepoststim_sep (s)', key='gui_cfg_bld_prepoststim_sep', format='%.6g')
                        st.number_input('prestim_time', key='gui_cfg_bld_prestim_time', format='%.6g')
                        st.number_input('poststim_time', key='gui_cfg_bld_poststim_time', format='%.6g')
                        st.number_input('amp_step_vel', key='gui_cfg_bld_amp_step_vel', min_value=1, step=1, format='%d')
                        st.selectbox(
                            'ramp_mode_default',
                            options=['linear', 'exponential'],
                            key='gui_cfg_bld_ramp_mode_default',
                        )
                        st.caption('`units` / `unit_rules` dicts stay from the template unless you edit the generated `.py` file.')
                st.checkbox(
                    'Overwrite if a `.py` file with that name already exists',
                    key='gui_cfg_build_overwrite',
                )
                if _load_save_button('Write config file and load', key='gui_btn_write_load_config'):
                    base = str(st.session_state.get('gui_cfg_build_base') or '').strip()
                    out_raw = str(st.session_state.get('gui_cfg_build_out') or '').strip()
                    out_stem = sanitize_config_module_stem(out_raw)
                    if not base:
                        _st_error_actions('Build new blocked.', ['Pick template module above'])
                    elif not out_raw:
                        _st_error_actions('Build new blocked.', ['Enter new module file name'])
                    else:
                        stim_ch = parse_comma_list(str(st.session_state.get('gui_cfg_bld_stim_channels') or ''))
                        sg_ch = parse_comma_list(str(st.session_state.get('gui_cfg_bld_SG_chan') or ''))
                        sg_nm = parse_comma_list(str(st.session_state.get('gui_cfg_bld_SG_name') or ''))
                        sm_ch = parse_comma_list(str(st.session_state.get('gui_cfg_bld_stim_monitor_chan') or ''))
                        sm_nm = parse_comma_list(str(st.session_state.get('gui_cfg_bld_stim_monitor_name') or ''))
                        sono_cal_ok = True
                        sono_lf: list[float] = []
                        sono_rf: list[float] = []
                        try:
                            sono_lf = parse_n_floats(str(st.session_state.get('gui_cfg_bld_sono_cal_left') or ''), 4)
                            sono_rf = parse_n_floats(str(st.session_state.get('gui_cfg_bld_sono_cal_right') or ''), 4)
                        except ValueError as e:
                            sono_cal_ok = False
                            _st_error_detail(
                                'Sono calibration invalid.',
                                ['Enter four comma-separated numbers', 'Match V and mm pairs'],
                                str(e),
                            )
                        if not stim_ch:
                            _st_error_actions('Stim channels missing.', ['Add comma-separated port lines'])
                        elif not sg_ch or not sg_nm or len(sg_ch) != len(sg_nm):
                            _st_error_actions('SG lists mismatch.', ['Match SG channel name counts'])
                        elif (bool(sm_ch) ^ bool(sm_nm)) or (sm_ch and sm_nm and len(sm_ch) != len(sm_nm)):
                            _st_error_actions('Stim monitor mismatch.', ['Fill both lists or clear both'])
                        elif not sono_cal_ok:
                            pass
                        else:
                            sono_ch = parse_comma_list(str(st.session_state.get('gui_cfg_bld_sono_channel') or ''))
                            sono_nm = parse_comma_list(str(st.session_state.get('gui_cfg_bld_sono_name') or ''))
                            use_sono = bool(st.session_state.get('gui_cfg_bld_use_sono'))
                            if use_sono and (not sono_ch or not sono_nm or len(sono_ch) != len(sono_nm)):
                                _st_error_actions('Sono lists mismatch.', ['Match sono name count', 'Or disable sonomicrometry'])
                            else:
                                _configs_dir = default_configs_dir(_ROOT)
                                os.makedirs(_configs_dir, exist_ok=True)
                                if _configs_dir not in sys.path:
                                    sys.path.insert(0, _configs_dir)
                                path = os.path.join(_configs_dir, out_stem + '.py')
                                if os.path.isfile(path) and not st.session_state.get('gui_cfg_build_overwrite'):
                                    _st_error_actions('Config file exists.', ['Enable overwrite checkbox', 'Or pick new name'])
                                else:
                                    rm = str(st.session_state.get('gui_cfg_bld_ramp_mode_default') or 'linear')
                                    if rm not in ('linear', 'exponential'):
                                        rm = 'linear'
                                    assignments = {
                                        'forcetorque_calibration_file': str(
                                            st.session_state.get('gui_cfg_bld_forcetorque_calibration_file') or 'FT56491.cal'
                                        ),
                                        'positive_motor_direction': str(
                                            st.session_state.get('gui_cfg_bld_positive_motor_direction') or 'left'
                                        ),
                                        'specimen_lateral_index_on_positive_motor_side': int(
                                            st.session_state.get('gui_cfg_bld_specimen_lateral_index') or -1
                                        ),
                                        'motor_axis': str(st.session_state.get('gui_cfg_bld_motor_axis') or 'z'),
                                        'bending_axis_sensor': str(
                                            st.session_state.get('gui_cfg_bld_bending_axis_sensor') or 'z'
                                        ),
                                        'primary_bending_axis': str(
                                            st.session_state.get('gui_cfg_bld_primary_bending_axis') or 'zTorque'
                                        ),
                                        'bending_axis_specimen': str(
                                            st.session_state.get('gui_cfg_bld_bending_axis_specimen') or 'dorsoventral'
                                        ),
                                        'device_name': str(st.session_state.get('gui_cfg_bld_device_name') or 'Dev1'),
                                        'daq_ai_sample_rate_hz': float(st.session_state.get('gui_cfg_bld_daq_ai_sr') or 1000.0),
                                        'daq_ao_do_sample_rate_hz': float(
                                            st.session_state.get('gui_cfg_bld_daq_ao_sr') or 60000.0
                                        ),
                                        'motor_full_steps_per_rev': int(st.session_state.get('gui_cfg_bld_motor_steps') or 1600),
                                        'motor_gear_ratio': int(st.session_state.get('gui_cfg_bld_motor_gear') or 5),
                                        'encoder_pulses_per_rev': int(st.session_state.get('gui_cfg_bld_encoder_ppr') or 10000),
                                        'stim_channels': stim_ch,
                                        'motor_port': str(st.session_state.get('gui_cfg_bld_motor_port') or 'port0'),
                                        'encoder_chan': str(st.session_state.get('gui_cfg_bld_encoder_chan') or 'ctr0'),
                                        'SG_chan': sg_ch,
                                        'SG_name': sg_nm,
                                        'stim_monitor_chan': sm_ch,
                                        'stim_monitor_name': sm_nm,
                                        'S1side': str(st.session_state.get('gui_cfg_bld_S1side') or 'left'),
                                        'S2side': str(st.session_state.get('gui_cfg_bld_S2side') or 'right'),
                                        'use_sono': use_sono,
                                        'sono_channel': sono_ch if use_sono else [],
                                        'sono_name': sono_nm if use_sono else [],
                                        'sono_internal_samplefreq': int(
                                            st.session_state.get('gui_cfg_bld_sono_internal_samplefreq') or 241
                                        ),
                                        'sono_cal_left': sono_lf,
                                        'sono_cal_right': sono_rf,
                                        'amp_step_vel': int(st.session_state.get('gui_cfg_bld_amp_step_vel') or 10),
                                        'ramp_mode_default': rm,
                                        'waitbefore': float(st.session_state.get('gui_cfg_bld_waitbefore') or 3.0),
                                        'waitafter': float(st.session_state.get('gui_cfg_bld_waitafter') or 4.0),
                                        'rampdur': float(st.session_state.get('gui_cfg_bld_rampdur') or 0.25),
                                        'prepoststim_dur': float(st.session_state.get('gui_cfg_bld_prepoststim_dur') or 0.06),
                                        'prepoststim_sep': float(st.session_state.get('gui_cfg_bld_prepoststim_sep') or 1.0),
                                        'prestim_time': float(st.session_state.get('gui_cfg_bld_prestim_time') or -2.0),
                                        'poststim_time': float(st.session_state.get('gui_cfg_bld_poststim_time') or 2.0),
                                    }
                                    src = render_generated_config(base, assignments)
                                    with open(path, 'w', encoding='utf-8') as f:
                                        f.write(src)
                                    importlib.invalidate_caches()
                                    err = _apply_loaded_config_module(out_stem)
                                    if err:
                                        _st_error_detail(
                                            'Generated config did not load.',
                                            ['Fix assignments above', 'Read Details message'],
                                            err,
                                        )
                                    else:
                                        st.session_state['gui_load_cfg_select'] = out_stem
                                        st.success(f'Wrote and loaded `{out_stem}`')
                                        st.rerun()

    if _show_data:
        _data_host = _setup_right if _setup_right is not None else st
        _ensure_hw_config_session_defaults()
        _data_host.subheader('2 · Data file path')
        _sw_dp = _stepwise_on_data_file_path_step()
        with _data_host.container(border=True):
            df_col, fn_col = st.columns(2)
            with df_col:
                _render_data_folder_dropdown(key_suffix='main')
                _preview_out = _compose_output_h5_path().strip()
                _preview_folder = str(st.session_state.get('gui_data_folder') or '').strip()
                if _preview_out:
                    st.caption(f'Save path: `{_preview_out}`')
                elif _preview_folder:
                    st.caption(f'Save path: `{_preview_folder}`')
            with fn_col:
                st.text_input(
                    ' ',
                    key='gui_data_filename',
                    placeholder='datafilename.h5',
                    help=DATA_FILE_NAME_HELP,
                    label_visibility='collapsed',
                )
            full_out = _compose_output_h5_path()
            if full_out:
                # Show "selected" state immediately; committing to the in-memory experiment still requires "Apply setup".
                b0 = st.session_state.get('bender')
                applied = str(getattr(b0, 'outputfile', '') or '').strip() if b0 is not None else ''
                if applied and _paths_equal_norm(applied, full_out):
                    st.success(f'**Save path:** `{full_out}`')
                else:
                    st.success(f'**Save path:** `{full_out}` (selected, not applied yet — click **Apply setup**)')
            if _data_path_apply_dirty():
                _soft_apply_reminder()
            _touch_data_path_baseline_if_clean()

        _apply_hlp = (
            'Applies selected hardware config and current data path to the experiment object. Does not start DAQ.'
            if _sw_dp
            else 'Applies hardware config and data path from setup sections. Does not start DAQ.'
        )
        if _load_save_button('Apply setup', key='gui_setup_apply_bottom', help=_apply_hlp):
            _apply_setup_action(sw_dp=bool(_sw_dp))

    if 'bender' not in st.session_state:
        st.stop()

    b: Bender = st.session_state['bender']
    _ensure_apply_tracking_bender(b)
    _init_biometrics_session_state(b, force=False)
    _rehydrate_missing_biometrics_from_bender(b)
    _sync_biometric_flags_from_session(b)
    # #region agent log
    _agent_debug_log(
        hypothesis_id='D',
        location='bender_streamlit_gui.py:main_post_bio_sync',
        message='session_vs_bender',
        data={
            'route': _nav_route(),
            'step': _stepwise_step() if _nav_route() == 'stepwise' else None,
            'show_sec2': _show_full_sec2(),
            'sess_dclamp': st.session_state.get('bio_dclamp'),
            'bender_dclamp': getattr(b, 'dclamp', None),
            'sess_fishmass': st.session_state.get('bio_fishmass'),
            'bender_fishmass': getattr(b, 'fishmass', None),
        },
    )
    # #endregion
    _ensure_review_file_selection(
        _candidate_review_files(_output_path_anchor_for_review(b)) if _output_path_anchor_for_review(b) else []
    )
    schema = b.get_dispatch_schema()
    test_types = list(schema['test_types'])
    _consume_pending_protocol_template(test_types)

    if _show_full_sec2():
        st.subheader('3 · Specimen')
        if _bio_apply_dirty():
            _soft_apply_reminder()

        with st.expander('About the measurements', expanded=False):
            st.markdown(
                '- Use this section to set specimen identity, morphometrics, experimental conditions, clamp geometry, and profile/inertial settings.\n'
                '- **Apply specimen** commits identity, morphometrics, and session conditions; **Apply clamp geometry & inertial correction** commits clamp geometry, profile, and the inertial flag.\n'
                '- Biometrics templates are separate from protocol templates.'
            )

        if bf := st.session_state.pop('gui_biometrics_load_feedback', None):
            ok_bf, txt_bf = bf
            if ok_bf:
                st.success(txt_bf)
            elif str(txt_bf).strip() == 'Choose a biometrics file first.':
                st.warning(txt_bf)
            else:
                _st_error_detail(
                    'Biometrics load failed.',
                    ['Check file format', 'Read Details'],
                    txt_bf,
                )
        if bfs := st.session_state.pop('gui_biometrics_save_feedback', None):
            ok_s, txt_s = bfs
            if ok_s:
                st.success(txt_s)
            else:
                _st_error_detail(
                    'Biometrics save failed.',
                    ['Check name and folder', 'Read Details'],
                    txt_s,
                )
        with st.expander('Load / save biometrics templates', expanded=False):
            st.caption('Optional: reuse biometrics snapshots (`.json`) in your selected data folder.')
            _df_check = str(st.session_state.get('gui_data_folder') or '').strip()
            _bio_tpl_dir = _shared_experiment_dir()
            if not (_df_check and os.path.isdir(os.path.normpath(_df_check))):
                st.caption(
                    f'**Data folder** is not set or not found on disk—listing saved template files from `{_bio_tpl_dir}` until '
                    '**section 2** points to a valid folder.'
                )
            _bio_tpl_list = list_biometrics_template_files(_bio_tpl_dir)
            _bio_opts: list = [None] + _bio_tpl_list
            if 'gui_biometrics_template_select' not in st.session_state:
                st.session_state['gui_biometrics_template_select'] = None
            _bio_pick = st.selectbox(
                'Biometrics file to load',
                _bio_opts,
                format_func=_biometrics_template_option_label,
                key='gui_biometrics_template_select',
            )
            if _load_save_button(
                'Load biometrics into form',
                key='gui_biometrics_btn_load',
                help='Fills biometrics widgets from the file. Then click **Apply specimen** and **Apply clamp geometry & inertial correction** (or **Apply all biometrics**).',
            ):
                if not _bio_pick:
                    st.session_state['gui_biometrics_load_feedback'] = (False, 'Choose a biometrics file first.')
                else:
                    st.session_state['gui_pending_biometrics_path'] = _bio_pick
                st.rerun()

            if 'gui_biometrics_new_name' not in st.session_state:
                st.session_state['gui_biometrics_new_name'] = ''
            st.text_input('Save biometrics as (name)', key='gui_biometrics_new_name', placeholder='e.g. Zebrafish adult default')
            st.text_area(
                'Description (optional)',
                key='gui_biometrics_new_desc',
                height=68,
                placeholder='Optional note saved inside the file.',
                help='Stored in the file metadata when you save.',
            )
            if 'gui_biometrics_overwrite' not in st.session_state:
                st.session_state['gui_biometrics_overwrite'] = False
            st.checkbox('Overwrite if same file name exists', key='gui_biometrics_overwrite')
            if _load_save_button('Save biometrics', key='gui_biometrics_btn_save'):
                _bn = str(st.session_state.get('gui_biometrics_new_name') or '').strip()
                _bd = str(st.session_state.get('gui_biometrics_new_desc') or '').strip()
                _bst = sanitize_biometrics_filename_stem(_bn or 'biometrics')
                _bout = os.path.normpath(os.path.join(_shared_experiment_dir(), f'{_bst}.json'))
                try:
                    if os.path.isfile(_bout) and not bool(st.session_state.get('gui_biometrics_overwrite')):
                        st.session_state['gui_biometrics_save_feedback'] = (
                            False,
                            f'File already exists: `{_bout}`. Enable **Overwrite** or change the name.',
                        )
                    else:
                        os.makedirs(os.path.dirname(_bout) or '.', exist_ok=True)
                        save_biometrics_template(
                            _bout,
                            name=_bn or _bst,
                            description=_bd,
                            session_state=st.session_state,
                        )
                        st.session_state['gui_biometrics_save_feedback'] = (True, f'Saved `{_bout}`')
                except Exception as e:
                    st.session_state['gui_biometrics_save_feedback'] = (False, f'{type(e).__name__}: {e}')
                st.rerun()

        st.divider()
        sub_specimen = sub_clamp_inertial = False
        with st.container():
            # Section 2 (Specimen): identity + morphometrics + session temperature + prep condition.
            # One Apply commits only these fields (per 4-section model, ux_spec §2.1/§3).
            with st.form('bio_form_specimen', clear_on_submit=False):
                id1, id2 = st.columns(2)
                with id1:
                    st.text_input(
                        'Genus-species',
                        key='gui_genus_species',
                        placeholder='e.g. Danio rerio',
                        help='Stored in the exported `.h5` under protocol metadata (`genus_species`) when you run or export.',
                    )
                with id2:
                    st.text_input(
                        'Specimen ID',
                        key='gui_specimen_id',
                        placeholder='e.g. fish-042 or prep code',
                        help='Primary specimen label; also written to `fishcode` on the experiment object for notebook compatibility.',
                    )
                if 'bio_segment' not in st.session_state:
                    st.session_state['bio_segment'] = ''
                st.text_input('Segment / preparation label (`segment`)', key='bio_segment', placeholder='e.g. whole body, hemi')

                st.divider()
                st.markdown('**Morphometrics**')
                if 'bio_fishlen_TL' not in st.session_state:
                    st.session_state['bio_fishlen_TL'] = 0.0
                st.number_input(
                    'Total Length (`fishlen_TL`, mm)',
                    min_value=0.0,
                    format='%.6g',
                    key='bio_fishlen_TL',
                )
                if 'bio_fishlen_SL' not in st.session_state:
                    st.session_state['bio_fishlen_SL'] = 0.0
                st.number_input(
                    'Standard Length (`fishlen_SL`, mm)',
                    min_value=0.0,
                    format='%.6g',
                    key='bio_fishlen_SL',
                )
                if 'bio_fishmass' not in st.session_state:
                    st.session_state['bio_fishmass'] = 0.0
                st.number_input('Mass `fishmass` (g)', min_value=0.0, format='%.6g', key='bio_fishmass')

                st.divider()
                st.markdown('**Session conditions**')
                if 'bio_temp_room' not in st.session_state:
                    st.session_state['bio_temp_room'] = 0.0
                st.number_input(
                    'Room temperature (`temp_C_room`, °C)',
                    min_value=-5.0,
                    max_value=60.0,
                    format='%.3f',
                    key='bio_temp_room',
                )
                if 'bio_temp_tank' not in st.session_state:
                    st.session_state['bio_temp_tank'] = 0.0
                st.number_input(
                    'Tank / bath temperature (`temp_C_tank`, °C)',
                    min_value=-5.0,
                    max_value=60.0,
                    format='%.3f',
                    key='bio_temp_tank',
                )
                st.text_input(
                    'Prep condition',
                    key='bio_prep_condition',
                    placeholder='e.g. anesthetized, recovered 24 h, fasted',
                    help='Free text (e.g. handling, anesthesia, recovery). Saved as `prep_condition` in protocol metadata on export.',
                )
                sub_specimen = st.form_submit_button(
                    'Apply specimen',
                    type='primary',
                    use_container_width=True,
                    help='Commits identity, morphometrics, session temperatures, and prep condition onto the experiment object.',
                )

            st.divider()
            st.subheader('4 · Clamp geometry & inertial correction')
            # Section 3: clamp geometry + mounted profile + inertial flag. One merged Apply
            # commits only these fields (replaces the former separate clamp / profile Applies).
            with st.form('bio_form_clamp_inertial', clear_on_submit=False):
                st.markdown('**Clamp geometry**')
                if 'bio_dclamp' not in st.session_state:
                    st.session_state['bio_dclamp'] = 0.0
                st.number_input(
                    'Test segment length = clamp spacing (`dclamp` / `test_segment_length_mm`, mm)',
                    min_value=0.001,
                    format='%.6g',
                    key='bio_dclamp',
                )
                if 'bio_dbend' not in st.session_state:
                    st.session_state['bio_dbend'] = 0.0
                st.number_input(
                    'Along-body distance to center of clamped test segment (mm)',
                    min_value=0.0,
                    format='%.6g',
                    key='bio_dbend',
                    help=BIO_DBEND_FIELD_HELP,
                )
                cw1, cw2 = st.columns(2)
                with cw1:
                    if 'bio_xsec' not in st.session_state:
                        st.session_state['bio_xsec'] = 0.0
                    st.number_input('Width `xsec_width` (mm)', min_value=0.001, format='%.6g', key='bio_xsec')
                    if 'bio_muscle_depth' not in st.session_state:
                        st.session_state['bio_muscle_depth'] = 0.0
                    st.number_input(
                        'Muscle depth `target_muscle_depth_mm` (mm)',
                        min_value=0.0,
                        format='%.6g',
                        key='bio_muscle_depth',
                        help=(
                            'Distance from the outer surface to the muscle/fiber layer used for strain↔curvature '
                            '(effective lever = xsec_width/2 − muscle_depth). 0 = surface strain (legacy).'
                        ),
                    )
                with cw2:
                    if 'bio_xsec_height' not in st.session_state:
                        st.session_state['bio_xsec_height'] = 0.0
                    st.number_input('Height `xsec_height` (mm)', min_value=0.001, format='%.6g', key='bio_xsec_height')
                co1, co2 = st.columns(2)
                with co1:
                    if 'bio_dvert' not in st.session_state:
                        st.session_state['bio_dvert'] = 0.0
                    st.number_input('Vertical offset `dvert` (mm)', min_value=0.0, format='%.6g', key='bio_dvert')
                with co2:
                    if 'bio_dhoriz' not in st.session_state:
                        st.session_state['bio_dhoriz'] = 0.0
                    st.number_input('Horizontal offset `dhoriz` (mm)', min_value=0.0, format='%.6g', key='bio_dhoriz')

                st.divider()
                st.markdown('**Mounted body profile (inertial model)**')
                st.selectbox(
                    'Typical density (sets g/mm³ on Apply)',
                    BIO_DENSITY_PRESET_LABELS,
                    key='bio_prof_rho_preset',
                    help=(
                        'Quick picks from literature-scale values: ~1.00 g/cm³ water-like, ~1.06 g/cm³ muscle/soft tissue, '
                        '~1.9 g/cm³ cortical bone. Values are copied into **Specimen density** when you click **Apply clamp '
                        'geometry & inertial correction** (or choose **Custom** and edit the number).'
                    ),
                )
                _flush_pending_bio_prof_rho_sync()
                st.number_input(
                    'Specimen density (g / mm³)',
                    min_value=1e-9,
                    format='%.6g',
                    key='bio_prof_rho',
                    help=(
                        'Mass density for the inertial model (`specimen_profile_density_g_per_mm3`). '
                        '1 g/cm³ = 1×10⁻³ g/mm³. Adjust after a preset or type your own.'
                    ),
                )
                if 'bio_prof_L' not in st.session_state:
                    st.session_state['bio_prof_L'] = 0.0
                st.number_input(
                    'Specimen outline length for profile model (mm)',
                    min_value=0.001,
                    format='%.6g',
                    key='bio_prof_L',
                    help='Length along the specimen used with the tapered outline below (`specimen_profile_length_mm`).',
                )
                p1, p2 = st.columns(2)
                with p1:
                    if 'bio_prof_ph' not in st.session_state:
                        st.session_state['bio_prof_ph'] = 0.0
                    st.number_input('Proximal height (mm)', min_value=0.001, format='%.6g', key='bio_prof_ph')
                    if 'bio_prof_pw' not in st.session_state:
                        st.session_state['bio_prof_pw'] = 0.0
                    st.number_input('Proximal width (mm)', min_value=0.001, format='%.6g', key='bio_prof_pw')
                with p2:
                    if 'bio_prof_dh' not in st.session_state:
                        st.session_state['bio_prof_dh'] = 0.0
                    st.number_input('Distal height (mm)', min_value=0.001, format='%.6g', key='bio_prof_dh')
                    if 'bio_prof_dw' not in st.session_state:
                        st.session_state['bio_prof_dw'] = 0.0
                    st.number_input('Distal width (mm)', min_value=0.001, format='%.6g', key='bio_prof_dw')
                p3, p4 = st.columns(2)
                with p3:
                    if 'bio_prof_clamp' not in st.session_state:
                        st.session_state['bio_prof_clamp'] = 0.0
                    st.number_input(
                        'Distance from rotation axis to clamps (mm)',
                        min_value=0.0,
                        format='%.6g',
                        key='bio_prof_clamp',
                        help=BIO_PROF_CLAMP_FIELD_HELP,
                    )
                with p4:
                    if 'bio_prof_samples' not in st.session_state:
                        st.session_state['bio_prof_samples'] = 0.0
                    st.number_input('Profile integration samples', min_value=20, max_value=400, step=10, key='bio_prof_samples')

                st.divider()
                if 'bio_use_theoretical_inertial' not in st.session_state:
                    st.session_state['bio_use_theoretical_inertial'] = False
                st.checkbox(
                    'Check here to perform inertial correction',
                    key='bio_use_theoretical_inertial',
                    help=(
                        'Subtracts model **system** inertia (from calibration, if loaded) and **specimen** inertia from the '
                        'profile above when correcting measured torque.'
                    ),
                )
                sub_clamp_inertial = st.form_submit_button(
                    'Apply clamp geometry & inertial correction',
                    type='primary',
                    use_container_width=True,
                    help='Commits clamp spacing/offsets, cross-section, mounted profile, density, and the inertial-correction flag.',
                )
        if sub_specimen:
            _apply_specimen_identity_to_bender(b)
            _apply_intrinsic_biometrics_to_bender(b)
            _apply_experimental_conditions_to_bender(b)
            st.toast('Specimen applied.')
        if sub_clamp_inertial:
            _sync_biometric_flags_from_session(b)
            if _apply_clamp_geometry_to_bender(b):
                _apply_mounted_profile_inertial_to_bender(b)
                st.session_state['gui_measurements_confirmed'] = True
                st.toast('Clamp geometry & inertial correction applied.')

        _touch_bio_apply_baseline_if_clean()

    if _show_sec3_through_6():

        st.divider()
        st.subheader('7 · Protocol / Run')

        with st.expander('Load protocol template (optional)', expanded=False):
            st.caption(
                'Templates set experiment type + procedure fields only (not biometrics). '
                'Use **Apply procedure** after loading.'
            )
            _tpl_folder_top = _shared_experiment_dir()
            _tpl_files_top = list_template_files(_tpl_folder_top)
            _tpl_options_top: list = [None] + _tpl_files_top
            if 'gui_protocol_template_select' not in st.session_state:
                st.session_state['gui_protocol_template_select'] = None
            _tpl_pick_top = st.selectbox(
                'Template to load',
                _tpl_options_top,
                format_func=_protocol_template_option_label,
                key='gui_protocol_template_select',
                help='Procedure files saved from this app (`.json` in the data folder).',
            )
            if _load_save_button(
                'Load template into form',
                key='gui_protocol_btn_load',
                help='Sets **Experiment type** and procedure widgets from the file. Then click **Apply** to copy onto the Bender object.',
            ):
                if not _tpl_pick_top:
                    st.session_state['gui_protocol_load_feedback'] = (False, 'Choose a template file first.')
                else:
                    st.session_state['gui_pending_protocol_template_path'] = _tpl_pick_top
                st.rerun()
        if fb := st.session_state.pop('gui_protocol_load_feedback', None):
            ok_fb, txt_fb = fb
            if ok_fb:
                st.success(txt_fb)
            else:
                _st_error_detail(
                    'Protocol load failed.',
                    ['Check template file', 'Read Details'],
                    txt_fb,
                )
        if 'test_type_select' not in st.session_state:
            st.session_state['test_type_select'] = None
        tt = st.selectbox('Experiment type (test_type)', test_types, key='test_type_select')
        b.test_type = tt

        st.session_state.setdefault('gui_exp_hide', False)
        st.caption(
            'Set procedure fields below, then **Apply procedure** or **Refresh experiment preview** (both buttons are at the '
            'bottom of **Procedure fields**).'
        )

        updates = {}
        sub_proc_apply = False
        sub_proc_save = False
        sub_proc_preview = False
        pv_pts = 6000

        with st.expander('Procedure fields', expanded=not bool(st.session_state.get('gui_exp_hide'))):
            if sf := st.session_state.pop('gui_protocol_save_feedback', None):
                ok_sf, txt_sf = sf
                if ok_sf:
                    st.success(txt_sf)
                else:
                    _st_error_detail(
                        'Protocol save failed.',
                        ['Check name and folder', 'Read Details'],
                        txt_sf,
                    )
            with st.form('gui_procedure_form', clear_on_submit=False):
                if _procedure_apply_dirty():
                    _soft_apply_reminder()
                if tt == 'isometric':
                    st.caption(
                        '**Isometric** turns strain or curvature targets into motor angles using **test segment length** '
                        'and **cross-section width** from **section 3** (same as clamp spacing `dclamp`). '
                        'Those values are copied when you use **Apply** in **section 3** (clamp / intrinsic / experimental / **Apply all**) or when you **Run**.'
                    )
                    st.markdown('**Required**')
                    for key in schema['isometric_required']:
                        if key in _BLOCK_SEQUENCE_PROCEDURE_KEYS:
                            continue
                        label = key.replace('_', ' ')
                        updates[key] = _render_field(
                            b,
                            key,
                            'float' if 'steps' not in key else 'int',
                            label,
                            help_text=ISOMETRIC_FIELD_HELP.get(key),
                        )
                    if 'isometric_num_steps' in updates and updates['isometric_num_steps'] is not None:
                        updates['isometric_num_steps'] = int(updates['isometric_num_steps'])
                    st.markdown('**Optional**')
                    for key in schema['isometric_optional']:
                        if key == 'isometric_stim_params':
                            st.markdown('**Stimulation**')
                            assembled = _render_isometric_stim_fields(b)
                            updates[key] = assembled
                        elif key == 'isometric_mode':
                            modes = list(ALL_AMPS_MODE_OPTIONS)
                            skm = _widget_key('isometric_mode')
                            cur_m = str(_get_session_value(b, key, 'strain'))
                            if skm not in st.session_state:
                                st.session_state[skm] = cur_m if cur_m in modes else 'strain'
                            updates[key] = st.selectbox(
                                'Isometric mode (units for initial/final)',
                                modes,
                                key=skm,
                                format_func=_format_strain_or_amp_mode,
                                help=ISOMETRIC_FIELD_HELP.get(key),
                            )
                        elif key == 'isometric_inter_step_interval_s':
                            skg = _widget_key('isometric_inter_step_interval_s')
                            if skg not in st.session_state:
                                v0 = _get_session_value(b, key, 0.0)
                                st.session_state[skg] = float(v0) if v0 is not None else 0.0
                            updates[key] = float(
                                st.number_input(
                                    'Time between steps (s)',
                                    min_value=0.0,
                                    format='%.6g',
                                    key=skg,
                                    help=(
                                        'Seconds to wait after each step finishes (after acquisition) before the next '
                                        'ramp/stimulus period. Use **0** for back-to-back steps.'
                                    ),
                                )
                            )
                        elif 'random_seed' in key:
                            sks = _widget_key(key)
                            if sks not in st.session_state:
                                v0 = _get_session_value(b, key)
                                st.session_state[sks] = '' if v0 is None else str(v0)
                            s = st.text_input('Random seed (optional)', key=sks, help=RANDOM_SEED_HELP)
                            if not str(s).strip():
                                updates[key] = None
                            else:
                                try:
                                    updates[key] = int(s)
                                except ValueError:
                                    _st_error_actions(
                                        'Random seed invalid.',
                                        ['Use whole number only', 'Or leave field blank'],
                                    )
                                    updates[key] = None
                        else:
                            kind = 'bool' if 'randomize' in key else 'str'
                            lbl = key.replace('_', ' ')
                            updates[key] = _render_field(
                                b, key, kind, lbl, help_text=ISOMETRIC_FIELD_HELP.get(key)
                            )
                    with st.container(border=True):
                        st.markdown('**Block sequence**')
                        _block_up = _render_block_sequence_fields(b)
                        if _block_up is not None:
                            updates.update(_block_up)

                elif tt == 'isovelocity':
                    st.markdown('**Required**')
                    for key in schema['isovelocity_required']:
                        if key in _BLOCK_SEQUENCE_PROCEDURE_KEYS:
                            continue
                        if key == 'isovelocity_starting_strain_mode':
                            modes = list(ALL_AMPS_MODE_OPTIONS)
                            skm = _widget_key('isovelocity_starting_strain_mode')
                            cur_m = str(_get_session_value(b, key, 'strain'))
                            if skm not in st.session_state:
                                st.session_state[skm] = cur_m if cur_m in modes else 'strain'
                            updates[key] = st.selectbox(
                                ISOVELOCITY_WIDGET_LABEL.get(key, 'Unit for starting posture'),
                                modes,
                                key=skm,
                                format_func=_format_strain_or_amp_mode,
                                help=ISOVELOCITY_FIELD_HELP.get(key),
                            )
                        elif key == 'isovelocity_velocity_mode':
                            vmodes = list(VELOCITY_MODE_OPTIONS)
                            skv = _widget_key('isovelocity_velocity_mode')
                            cur_v = str(_get_session_value(b, key, 'angle_vel'))
                            if skv not in st.session_state:
                                st.session_state[skv] = cur_v if cur_v in vmodes else 'angle_vel'
                            updates[key] = st.selectbox(
                                ISOVELOCITY_WIDGET_LABEL.get(key, 'Unit for min/max velocity'),
                                vmodes,
                                key=skv,
                                format_func=_format_velocity_mode,
                                help=ISOVELOCITY_FIELD_HELP.get(key),
                            )
                        else:
                            kind = 'int' if 'num_steps' in key else 'float'
                            lbl = ISOVELOCITY_WIDGET_LABEL.get(key, key.replace('_', ' '))
                            updates[key] = _render_field(
                                b, key, kind, lbl, help_text=ISOVELOCITY_FIELD_HELP.get(key)
                            )
                    if 'isovelocity_num_steps' in updates and updates['isovelocity_num_steps'] is not None:
                        updates['isovelocity_num_steps'] = int(updates['isovelocity_num_steps'])
                    st.markdown('**Optional**')
                    for key in schema['isovelocity_optional']:
                        if key == 'isovelocity_stim_params':
                            st.markdown('**Stimulation**')
                            assembled = _render_isovelocity_stim_fields(b)
                            updates[key] = assembled
                        elif key == 'isovelocity_velocity_mode':
                            pass  # required-only; rendered in Required section
                        elif 'random_seed' in key:
                            sks = _widget_key(key)
                            if sks not in st.session_state:
                                v0 = _get_session_value(b, key)
                                st.session_state[sks] = '' if v0 is None else str(v0)
                            s = st.text_input('Random seed (optional)', key=sks, help=RANDOM_SEED_HELP)
                            if not str(s).strip():
                                updates[key] = None
                            else:
                                try:
                                    updates[key] = int(s)
                                except ValueError:
                                    _st_error_actions(
                                        'Random seed invalid.',
                                        ['Use whole number only', 'Or leave field blank'],
                                    )
                                    updates[key] = None
                        else:
                            kind = 'bool' if 'randomize' in key else 'float'
                            lbl = ISOVELOCITY_WIDGET_LABEL.get(key, key.replace('_', ' '))
                            updates[key] = _render_field(
                                b, key, kind, lbl, help_text=ISOVELOCITY_FIELD_HELP.get(key)
                            )
                    with st.container(border=True):
                        st.markdown('**Block sequence**')
                        _block_up = _render_block_sequence_fields(b)
                        if _block_up is not None:
                            updates.update(_block_up)

                elif tt == 'calibration':
                    st.markdown('**Required**')
                    bases = [x for x in test_types if x != 'calibration']
                    cur = _get_session_value(b, 'calibration_base_test_type', 'dynamic')
                    if cur not in bases:
                        cur = 'dynamic'
                    sk_cal = _widget_key('calibration_base_test_type')
                    if sk_cal not in st.session_state:
                        st.session_state[sk_cal] = cur
                    if st.session_state[sk_cal] not in bases:
                        st.session_state[sk_cal] = bases[0]
                    updates['calibration_base_test_type'] = st.selectbox(
                        'Calibration base test type', bases, key=sk_cal
                    )
                    st.info(
                        'Calibration runs the **base** motion protocol. Set **test_type** to that base '
                        '(e.g. dynamic), click **Apply** in **section 4** or **6** (and **Refresh experiment preview** in **Procedure fields** if you use preview), '
                        'then switch back to **calibration** before running.'
                    )
                    st.markdown('**Optional**')
                    for key in schema['calibration_optional']:
                        sko = _widget_key(key)
                        if sko not in st.session_state:
                            st.session_state[sko] = str(_get_session_value(b, key) or '')
                        updates[key] = st.text_input(key.replace('_', ' '), key=sko)

                elif tt in MOTION_TYPES:
                    st.markdown('**Motion-series parameters** (procedure-specific)')
                    if tt == 'dynamic':
                        st.info('Dynamic timing uses cycles (not Duration).')
                    elif tt == 'step_change':
                        st.info('Step-change timing is derived from step-change cycle parameters.')
                    fields = _motion_parameter_rows(tt)
                    _stim_field_names = {
                        'stim_cycles_in_step',
                        'is_stim',
                        'stim_pulse_rate',
                        'S1volts',
                        'S2volts',
                        'all_stimduties',
                        'all_stimphases',
                    }
                    _motion_fields = [row for row in fields if row[0] not in _stim_field_names]
                    _stim_fields = [row for row in fields if row[0] in _stim_field_names]
                    mcol, scol = st.columns([1.35, 1.0], gap='large')
                    with mcol:
                        st.markdown('**Motion controls**')
                        for name, kind, label in _motion_fields:
                            updates[name] = _render_field(
                                b, name, kind, label, help_text=MOTION_FIELD_HELP.get(name)
                            )
                    with scol:
                        st.markdown('**Stimulation**')
                        for name, kind, label in _stim_fields:
                            updates[name] = _render_field(
                                b, name, kind, label, help_text=MOTION_FIELD_HELP.get(name)
                            )

                else:
                    st.warning(f'No dedicated field panel for {tt!r} yet; use notebook or extend this script.')

                # Persistent anti-bleed warning: re-evaluated every render so it stays visible
                # while the stim-timing condition holds, instead of flashing for one rerun
                # after Apply. Returns None when stim is disabled or timing is valid.
                _stim_timing_warn = _validate_procedure_stim_timing(b, updates, tt)
                if _stim_timing_warn:
                    st.warning(_stim_timing_warn)

                st.divider()
                if 'gui_protocol_show_save_template' not in st.session_state:
                    st.session_state['gui_protocol_show_save_template'] = False
                _show_tpl_save = st.checkbox(
                    'Show "Save procedure as template"',
                    key='gui_protocol_show_save_template',
                    value=False,
                )
                if _show_tpl_save:
                    st.markdown('**Save procedure as template**')
                    if 'gui_protocol_new_name' not in st.session_state:
                        st.session_state['gui_protocol_new_name'] = ''
                    st.text_input(
                        'Template name',
                        key='gui_protocol_new_name',
                        placeholder='e.g. Protocol A (any test_type)',
                    )
                    st.text_area(
                        'Description (optional)',
                        key='gui_protocol_new_desc',
                        height=70,
                        placeholder='e.g. Isometric 5 steps; or dynamic 1/3/5 Hz x strains; or calibration + base',
                    )
                    if 'gui_protocol_overwrite' not in st.session_state:
                        st.session_state['gui_protocol_overwrite'] = False
                    st.checkbox('Overwrite if a file with the same name already exists', key='gui_protocol_overwrite')
                    _pc1, _pc2 = st.columns(2)
                    with _pc1:
                        sub_proc_apply = st.form_submit_button(
                            'Apply procedure',
                            use_container_width=True,
                            help='Copy procedure fields onto the experiment object (not **Run experiment**).',
                        )
                    with _pc2:
                        sub_proc_save = st.form_submit_button('Save template', use_container_width=True)
                else:
                    sub_proc_apply = st.form_submit_button(
                        'Apply procedure',
                        use_container_width=True,
                        help='Copy procedure fields onto the experiment object (not **Run experiment**).',
                    )
                sub_proc_preview = st.form_submit_button(
                    'Refresh experiment preview',
                    use_container_width=True,
                    help=(
                        'Submit this form: copy procedure fields onto the experiment object and rebuild the preview plot. '
                        'Use after editing fields here so preview matches your inputs (same data path as **Apply procedure**).'
                    ),
                )

            if sub_proc_apply:
                _form_invalid = (
                    (tt == 'isometric' and updates.get('isometric_stim_params') is None)
                    or (tt == 'isovelocity' and updates.get('isovelocity_stim_params') is None)
                    or (tt in ('isometric', 'isovelocity') and 'block_sequence' not in updates)
                )
                if _form_invalid:
                    pass  # widget validation already surfaced st.error
                else:
                    _stim_err = _validate_procedure_stim_timing(b, updates, tt)
                    if _stim_err:
                        st.error(_stim_err)
                    else:
                        _apply_procedure_form_to_bender(b, updates, tt)
                        st.toast('Settings applied.')
            if sub_proc_save:
                _name = str(st.session_state.get('gui_protocol_new_name') or '').strip()
                _desc = str(st.session_state.get('gui_protocol_new_desc') or '').strip()
                _stem = sanitize_template_filename_stem(_name or 'protocol')
                _out = os.path.normpath(os.path.join(_shared_experiment_dir(), f'{_stem}.json'))
                try:
                    if os.path.isfile(_out) and not bool(st.session_state.get('gui_protocol_overwrite')):
                        st.session_state['gui_protocol_save_feedback'] = (
                            False,
                            f'File already exists: `{_out}`. Enable **Overwrite** or pick a different name.',
                        )
                    else:
                        os.makedirs(os.path.dirname(_out) or '.', exist_ok=True)
                        _proc = build_procedure_dict_from_updates(updates)
                        _base = None
                        if tt == 'calibration':
                            _btt = _proc.get('calibration_base_test_type')
                            if _btt:
                                _snap = snapshot_bender_procedure(b, schema, str(_btt))
                                if _snap:
                                    _base = {'test_type': str(_btt), 'procedure': _snap}
                        save_protocol_template(
                            _out,
                            name=_name or _stem,
                            description=_desc,
                            test_type=tt,
                            procedure=_proc,
                            base_protocol=_base,
                        )
                        st.session_state['gui_protocol_save_feedback'] = (True, f'Saved `{_out}`')
                except Exception as e:
                    st.session_state['gui_protocol_save_feedback'] = (False, f'{type(e).__name__}: {e}')
                st.rerun()

            if sub_proc_preview:
                _form_invalid = (
                    (tt == 'isometric' and updates.get('isometric_stim_params') is None)
                    or (tt == 'isovelocity' and updates.get('isovelocity_stim_params') is None)
                    or (tt in ('isometric', 'isovelocity') and 'block_sequence' not in updates)
                )
                if _form_invalid:
                    pass
                else:
                    _stim_err = _validate_procedure_stim_timing(b, updates, tt)
                    if _stim_err:
                        st.error(_stim_err)
                    else:
                        _sync_biometric_flags_from_session(b)
                        _sync_genus_species_to_bender(b)
                        _apply_form_updates(b, updates, tt)
                        _mark_procedure_applied()
                        st.session_state['gui_last_preview'] = build_protocol_preview(
                            b, requested_test_type=tt, max_plot_points=int(pv_pts)
                        )
                        st.session_state['gui_last_preview_tt'] = tt
                        if st.session_state['gui_last_preview'].get('ok'):
                            st.session_state['gui_protocol_confirmed'] = True
                        st.toast('Preview updated.')
                        st.rerun()

            _touch_proc_apply_baseline_if_clean()

        st.checkbox(
            'Hide section (values stay; unhide to edit)',
            key='gui_exp_hide',
            help='Collapse **Procedure fields** after you finish editing, saving a template, or using **Apply procedure** / **Refresh experiment preview**.',
        )

        st.divider()
        st.subheader('8 · Experiment preview')
        if _procedure_apply_dirty() or _bio_apply_dirty():
            _soft_apply_reminder()

        def _render_current_settings_table() -> None:
            _sync_biometric_flags_from_session(b)
            _sync_genus_species_to_bender(b)
            _apply_form_updates(b, updates, tt)
            _mark_procedure_applied()
            settings_rows = [
                {'group': 'experiment', 'name': 'test_type', 'value': tt},
                {
                    'group': 'export',
                    'name': 'data_file_target_h5',
                    'value': _compose_output_h5_path() or getattr(b, 'outputfile', None),
                },
                {
                    'group': 'specimen',
                    'name': 'genus_species',
                    'value': (getattr(b, 'h5_protocol_metadata', {}) or {}).get('genus_species', ''),
                },
                {
                    'group': 'specimen',
                    'name': 'specimen_id',
                    'value': (getattr(b, 'h5_protocol_metadata', {}) or {}).get('specimen_id', ''),
                },
                {'group': 'biometric', 'name': 'test_segment_length_mm', 'value': getattr(b, 'dclamp', None)},
                {'group': 'biometric', 'name': 'test_segment_position_mm', 'value': getattr(b, 'dbend', None)},
                {'group': 'biometric', 'name': 'xsec_width', 'value': getattr(b, 'xsec_width', None)},
                {'group': 'biometric', 'name': 'target_muscle_depth_mm', 'value': getattr(b, 'target_muscle_depth_mm', None)},
                {'group': 'biometric', 'name': 'dvert', 'value': getattr(b, 'dvert', None)},
                {'group': 'biometric', 'name': 'dhoriz', 'value': getattr(b, 'dhoriz', None)},
                {
                    'group': 'conditions',
                    'name': 'temp_C_room',
                    'value': getattr(b, 'temp_C_room', None),
                },
                {
                    'group': 'conditions',
                    'name': 'temp_C_tank',
                    'value': getattr(b, 'temp_C_tank', None),
                },
                {
                    'group': 'conditions',
                    'name': 'prep_condition',
                    'value': (getattr(b, 'h5_protocol_metadata', {}) or {}).get('prep_condition', ''),
                },
            ]
            for k, v in sorted(updates.items(), key=lambda kv: kv[0]):
                settings_rows.append({'group': 'parameter', 'name': k, 'value': str(v)})
            _settings_df = pd.DataFrame(settings_rows)
            # Keep the value column a single dtype (string): mixing str/float/None triggers a
            # pyarrow "mixed-type column" warning when Streamlit serializes the DataFrame.
            if 'value' in _settings_df.columns:
                _settings_df['value'] = _settings_df['value'].apply(
                    lambda x: '' if x is None else str(x)
                )
            st.dataframe(_settings_df, use_container_width=True, hide_index=True)

        if st.session_state.get('gui_last_preview') is not None:
            prev = st.session_state['gui_last_preview']
            if st.session_state.get('gui_last_preview_tt') != tt:
                st.warning(
                    'Test type changed since the last preview — open **Procedure fields** and click **Refresh experiment preview**.'
                )
            if prev.get('error'):
                _ph, _pb = _preview_error_actions(str(prev.get('error') or ''))
                _st_error_actions(
                    _ph,
                    _pb,
                )
                with st.expander('Preview error detail'):
                    st.code(str(prev['error']))
            elif prev.get('ok'):
                    if tt == 'calibration':
                        st.info(
                            f"Calibration uses base protocol **{prev.get('motion_test_type')}** for motion "
                            '(same as **Run experiment**).'
                        )
                    if prev.get('table'):
                        st.markdown('**Summary table**')
                        st.dataframe(pd.DataFrame(prev['table']), use_container_width=True, hide_index=True)
                    tp = prev.get('t_plot')
                    ap = prev.get('angle_plot')
                    vp = prev.get('anglevel_plot')
                    if tp is not None and ap is not None and len(tp) > 0:
                        if tt == 'isometric' and prev.get('preview_isometric'):
                            st.markdown('**Command preview** (isometric: ramp–hold per step, same timing as run)')
                        elif tt == 'isovelocity' and prev.get('preview_isovelocity'):
                            st.markdown(
                                '**Command preview** (isovelocity: pre-hold at start angle + constant ω segment per step, '
                                'same timing as run)'
                            )
                        else:
                            st.markdown('**Command preview** (motor angle and angular velocity)')
                        fig = go.Figure()
                        fig.add_trace(go.Scatter(x=tp, y=ap, mode='lines', name='Commanded angle (deg)'))
                        if vp is not None and len(vp) > 0:
                            fig.add_trace(
                                go.Scatter(
                                    x=tp,
                                    y=vp,
                                    mode='lines',
                                    name='Commanded anglevel (deg/s)',
                                    yaxis='y2',
                                )
                            )
                        if vp is not None and len(vp) > 0:
                            fig.update_layout(
                                height=420,
                                margin=dict(l=48, r=48, t=40, b=40),
                                xaxis_title='Time (s)',
                                yaxis=dict(title='Angle (deg)'),
                                yaxis2=dict(title='Anglevel (deg/s)', overlaying='y', side='right'),
                                legend=dict(yanchor='top', y=0.99, xanchor='left', x=0.01),
                            )
                        else:
                            fig.update_layout(
                                height=420,
                                margin=dict(l=48, r=48, t=40, b=40),
                                xaxis_title='Time (s)',
                                yaxis=dict(title='Angle (deg)'),
                                legend=dict(yanchor='top', y=0.99, xanchor='left', x=0.01),
                            )
                        st.plotly_chart(fig, use_container_width=True)
                        sp_plot = prev.get('strain_plot')
                        kp_plot = prev.get('curvature_plot')
                        if sp_plot is not None or kp_plot is not None:
                            st.markdown('**Native units** (derived from motor angle)')
                            fig_nu = go.Figure()
                            if kp_plot is not None and len(kp_plot) > 0:
                                fig_nu.add_trace(
                                    go.Scatter(x=tp, y=kp_plot, mode='lines', name='Curvature κ (1/m)')
                                )
                            if sp_plot is not None and len(sp_plot) > 0:
                                fig_nu.add_trace(
                                    go.Scatter(
                                        x=tp,
                                        y=sp_plot,
                                        mode='lines',
                                        name='Surface strain ε',
                                        yaxis='y2' if kp_plot is not None else 'y',
                                    )
                                )
                            _nu_layout = dict(
                                height=360,
                                margin=dict(l=48, r=48, t=40, b=40),
                                xaxis_title='Time (s)',
                                legend=dict(yanchor='top', y=0.99, xanchor='left', x=0.01),
                            )
                            if sp_plot is not None and kp_plot is not None:
                                _nu_layout['yaxis'] = dict(title='κ (1/m)')
                                _nu_layout['yaxis2'] = dict(title='ε', overlaying='y', side='right')
                            elif kp_plot is not None:
                                _nu_layout['yaxis'] = dict(title='κ (1/m)')
                            else:
                                _nu_layout['yaxis'] = dict(title='ε')
                            fig_nu.update_layout(**_nu_layout)
                            st.plotly_chart(fig_nu, use_container_width=True)
                        stim_plot = prev.get('stim_plot')
                        if stim_plot is not None and len(stim_plot) > 0:
                            st.markdown('**Stimulation preview** (commanded AO voltage)')
                            fig_st = go.Figure()
                            fig_st.add_trace(
                                go.Scatter(x=tp, y=stim_plot, mode='lines', name='Total stim (S1+S2) (V)')
                            )
                            fig_st.update_layout(
                                height=320,
                                margin=dict(l=48, r=48, t=40, b=40),
                                xaxis_title='Time (s)',
                                yaxis_title='Voltage (V)',
                                legend=dict(yanchor='top', y=0.99, xanchor='left', x=0.01),
                            )
                            st.plotly_chart(fig_st, use_container_width=True)
                        elif prev.get('table') and any(
                            row.get('metric') == 'stimulation enabled' and row.get('value') is False
                            for row in (prev.get('table') or [])
                            if isinstance(row, dict)
                        ):
                            st.caption('Stimulation disabled for this protocol.')
                        if bool(st.session_state.get('gui_simulation_mode', False)) and prev.get('t') is not None:
                            t_prev = np.asarray(prev['t'], dtype=float).reshape(-1)
                            ang_prev = np.asarray(prev['angle'], dtype=float).reshape(-1)
                            av_prev = prev.get('anglevel')
                            if av_prev is None:
                                av_prev = np.zeros_like(ang_prev)
                            else:
                                av_prev = np.asarray(av_prev, dtype=float).reshape(-1)
                            n_prev = min(t_prev.size, ang_prev.size, av_prev.size)
                            if n_prev >= 2:
                                Lpv = getattr(b, 'dclamp', None) or getattr(b, 'test_segment_length_mm', None)
                                mat = str(st.session_state.get('gui_simulation_material', 'polyurethane'))
                                d_mm, F_N = force_displacement_series(
                                    ang_prev[:n_prev],
                                    av_prev[:n_prev],
                                    length_mm=Lpv,
                                    material_key=mat,
                                    max_points=int(pv_pts),
                                )
                                fig_fd = go.Figure()
                                fig_fd.add_trace(
                                    go.Scatter(
                                        x=d_mm,
                                        y=F_N,
                                        mode='lines',
                                        name='Simulated F(δ) along protocol',
                                    )
                                )
                                fig_fd.update_layout(
                                    title='Simulation preview: bending force vs tip displacement (25.4 mm OD cantilever)',
                                    xaxis_title='Tip displacement δ (mm)',
                                    yaxis_title='Bending reaction F (N)',
                                    height=420,
                                    margin=dict(l=48, r=48, t=48, b=40),
                                    legend=dict(yanchor='top', y=0.99, xanchor='left', x=0.01),
                                )
                                st.markdown(
                                    '**Simulation · force–displacement** (material from sidebar; span from **section 3 · clamp length**). '
                                    'Cyclic protocols trace a path (hysteresis-like from η·dδ/dt).'
                                )
                                st.plotly_chart(fig_fd, use_container_width=True)
                                dm, Fpu, Fsi = static_stiffness_comparison_delta_grid(Lpv)
                                fig_cmp = go.Figure()
                                fig_cmp.add_trace(
                                    go.Scatter(x=dm, y=Fpu, mode='lines', name='Polyurethane (E ≈ 35 MPa)')
                                )
                                fig_cmp.add_trace(
                                    go.Scatter(x=dm, y=Fsi, mode='lines', name='Silicone (E ≈ 3 MPa)')
                                )
                                fig_cmp.update_layout(
                                    title='Quasi-static comparison (same δ; polyurethane is much steeper)',
                                    xaxis_title='Tip displacement δ (mm)',
                                    yaxis_title='Elastic force F = kδ (N)',
                                    height=380,
                                    margin=dict(l=48, r=48, t=48, b=40),
                                    legend=dict(yanchor='top', y=0.99, xanchor='left', x=0.01),
                                )
                                st.plotly_chart(fig_cmp, use_container_width=True)
                    elif prev.get('table') and tt in ('isometric', 'isovelocity'):
                        st.caption(
                            'Step protocols: table lists setpoints; refresh preview after fixing errors to see the plot.'
                        )
            else:
                st.warning(
                    'Preview incomplete — open **Procedure fields** and click **Refresh experiment preview** (bottom of that form).'
                )

        st.divider()
        st.markdown('### Run controls')
        if bool(st.session_state.get('gui_simulation_mode', False)):
            st.info('Simulation mode active: run uses numpy only (no NI-DAQ).')
        if _procedure_apply_dirty() or _bio_apply_dirty():
            _soft_apply_reminder()

        def _execute_run() -> None:
            if bool(st.session_state.get('gui_run_in_progress', False)):
                st.warning('A run is already in progress.')
                return
            _proc_ok, _proc_msg = _procedure_ready_for_run()
            if not _proc_ok:
                st.error(_proc_msg)
                return
            st.session_state['gui_run_in_progress'] = True
            _rehydrate_missing_biometrics_from_bender(b)
            _sync_genus_species_to_bender(b)
            outp = _compose_output_h5_path().strip()
            if outp:
                b.outputfile = outp
                _mark_data_path_applied()
            notes_in = str(st.session_state.get('gui_post_notes') or '').strip()
            try:
                b.simulation_mode = bool(st.session_state.get('gui_simulation_mode', False))
                b.simulation_material = str(st.session_state.get('gui_simulation_material', 'polyurethane'))
                _acq_label = (
                    'Simulating acquisition (no DAQ)…'
                    if b.simulation_mode
                    else 'Acquiring (DAQ)…'
                )
                _status_factory = getattr(st, 'status', None)
                if callable(_status_factory):
                    with _status_factory('Run in progress…', expanded=True) as run_status:
                        run_status.write('Acquisition…')
                        with st.spinner(_acq_label):
                            b.run_experiment(test_type=tt)
                        st.success('Acquisition finished.')
                        run_status.write('HDF5…')
                        with st.spinner('Writing data file (.h5)…'):
                            rep = export_primary_h5(
                                b,
                                post_trial_notes=notes_in if notes_in else None,
                                outputfile=outp or None,
                                append_post_trial_notes=bool(st.session_state.get('gui_qc_notes_append', True)),
                            )
                        qix = _read_qc_trial_index(b)
                        sel_h5 = str(st.session_state.get('gui_review_selected') or '').strip()
                        qc_base = _qc_figure_base_path(b, sel_h5, qix)
                        run_status.write('QC plot…')
                        with st.spinner('Saving QC plot…'):
                            qc_path, _ = save_universal_qc_figure(b, qc_trial_index=qix, base_path=qc_base)
                        _st_done = getattr(run_status, 'update', None)
                        if callable(_st_done):
                            _st_done(label='Run finished', state='complete')
                        st.success('Data has been saved! Check data folder to confirm before proceeding.')
                        st.info(f"Data file: `{rep['outputfile']}`  |  QC plot: `{qc_path}`")
                        _arm_acquired_trial_review(rep.get('outputfile'), qc_path)
                        if bool(st.session_state.get('gui_qc_notes_append', True)):
                            st.session_state['gui_post_notes'] = ''
                        else:
                            st.session_state['gui_post_notes'] = str(rep.get('post_trial_notes') or '')
                else:
                    with st.spinner(_acq_label):
                        b.run_experiment(test_type=tt)
                    st.success('Acquisition finished.')
                    with st.spinner('Writing data file (.h5)…'):
                        rep = export_primary_h5(
                            b,
                            post_trial_notes=notes_in if notes_in else None,
                            outputfile=outp or None,
                            append_post_trial_notes=bool(st.session_state.get('gui_qc_notes_append', True)),
                        )
                    qix = _read_qc_trial_index(b)
                    sel_h5 = str(st.session_state.get('gui_review_selected') or '').strip()
                    qc_base = _qc_figure_base_path(b, sel_h5, qix)
                    with st.spinner('Saving QC plot…'):
                        qc_path, _ = save_universal_qc_figure(b, qc_trial_index=qix, base_path=qc_base)
                    st.success('Data has been saved! Check data folder to confirm before proceeding.')
                    st.info(f"Data file: `{rep['outputfile']}`  |  QC plot: `{qc_path}`")
                    _arm_acquired_trial_review(rep.get('outputfile'), qc_path)
                    if bool(st.session_state.get('gui_qc_notes_append', True)):
                        st.session_state['gui_post_notes'] = ''
                    else:
                        st.session_state['gui_post_notes'] = str(rep.get('post_trial_notes') or '')
            except Exception as e:
                _show_friendly_error(e, action='run_experiment')
            finally:
                st.session_state['gui_run_in_progress'] = False

        if st.button('View experiment settings', use_container_width=True):
            _render_current_settings_table()

        st.session_state.setdefault('gui_run_soft_warnings', [])
        st.session_state.setdefault('gui_run_in_progress', False)
        _pending_run_confirm = bool(st.session_state.get('gui_run_pending_confirm', False))
        _run_disabled, _run_help = _run_button_state()
        _ready_run = _workflow_ready_state(b, tt)
        if bool(st.session_state.get('gui_run_in_progress', False)):
            st.warning(_run_help)
            if st.button('Reset run state', key='gui_run_reset_inprogress', use_container_width=True, type='secondary'):
                st.session_state['gui_run_in_progress'] = False
                st.session_state['gui_run_pending_confirm'] = False
                st.session_state['gui_run_soft_warnings'] = []
                st.rerun()
        elif _run_disabled:
            st.warning(_run_help)
            if st.button('Reset run state', key='gui_run_reset_stuck', use_container_width=True, type='secondary'):
                st.session_state['gui_run_in_progress'] = False
                st.session_state['gui_run_pending_confirm'] = False
                st.session_state['gui_run_soft_warnings'] = []
                st.rerun()
        elif _ready_run['protocol_ok'] and _ready_run['setup_ok'] and _ready_run['measurements_ok']:
            st.caption('Checklist complete — review warnings below if any, then run or click Proceed.')
        if st.button(
            'Run experiment',
            type='primary',
            use_container_width=True,
            disabled=_run_disabled,
            help=_run_help,
        ):
            _proc_ok, _proc_msg = _procedure_ready_for_run()
            if not _proc_ok:
                st.error(_proc_msg)
            else:
                _rehydrate_missing_biometrics_from_bender(b)
            run_warnings: list[str] = []
            bio_soft_missing: list[str] = []
            if not str(st.session_state.get('gui_specimen_id') or '').strip():
                bio_soft_missing.append('specimen ID')
            if not str(st.session_state.get('gui_genus_species') or '').strip():
                bio_soft_missing.append('genus-species')
            if _session_float('bio_dclamp') is None or _session_float('bio_dclamp') <= 0:
                bio_soft_missing.append('clamp spacing')
            if _session_float('bio_xsec') is None or _session_float('bio_xsec') <= 0:
                bio_soft_missing.append('cross-section width')
            if _session_float('bio_fishmass') is None or _session_float('bio_fishmass') <= 0:
                bio_soft_missing.append('mass')
            if bio_soft_missing:
                _missing_txt = ', '.join(bio_soft_missing[:4])
                if len(bio_soft_missing) > 4:
                    _missing_txt += ', ...'
                run_warnings.append(f'Biometrics look incomplete ({_missing_txt}).')
            if _needs_missing_calibration_confirmation(b) and not bool(st.session_state.get('gui_simulation_mode', False)):
                run_warnings.append('No calibration file detected.')
            if _section2_destination_incomplete():
                run_warnings.append('No designated file destination in section 2.')
            if _proc_ok:
                rep = b.validate_dispatch_setup(test_type=tt)
                if not rep.get('ok', False):
                    run_warnings.append('Required protocol fields are missing.')
                if run_warnings:
                    st.session_state['gui_run_soft_warnings'] = run_warnings
                    st.session_state['gui_run_pending_confirm'] = True
                else:
                    st.session_state['gui_run_soft_warnings'] = []
                    st.session_state['gui_run_pending_confirm'] = False
                    _execute_run()
            _pending_run_confirm = bool(st.session_state.get('gui_run_pending_confirm', False))

        if _pending_run_confirm:
            _warns = list(st.session_state.get('gui_run_soft_warnings') or [])
            if _warns:
                st.warning('\n'.join([f'- {w}' for w in _warns]))
            c_go, c_stop = st.columns(2)
            with c_go:
                if st.button(
                    'Proceed',
                    type='primary',
                    use_container_width=True,
                    key='gui_run_proceed',
                    disabled=bool(st.session_state.get('gui_run_in_progress', False)),
                ):
                    st.session_state['gui_run_pending_confirm'] = False
                    st.session_state['gui_run_soft_warnings'] = []
                    _execute_run()
            with c_stop:
                if st.button('Abort', use_container_width=True, key='gui_run_abort'):
                    st.session_state['gui_run_pending_confirm'] = False
                    st.session_state['gui_run_soft_warnings'] = []
                    st.info('Run aborted.')

        # Persist-then-review: the raw .h5 is already written by the run above. Show the
        # acquired trial now and let the operator keep or delete it (save is not gated on this).
        _render_acquired_trial_review(b)

        st.divider()
        st.session_state.setdefault('gui_sec6_hide', False)
        st.subheader('7 · Save data here')
        st.caption('Writes from current in-memory **trial_records** (after **Run** or a prior save). No new DAQ.')
        if st.session_state.get('gui_sec6_hide'):
            st.caption('Save controls hidden. Uncheck **Hide section** below.')
        else:

            def _export_h5_from_session():
                _sync_genus_species_to_bender(b)
                outp = _compose_output_h5_path().strip()
                if not outp and not getattr(b, 'outputfile', None):
                    _st_error_actions(
                        'No HDF5 path set.',
                        ['Open section 2', 'Set data folder', 'Set file name'],
                    )
                    return None
                notes_in = str(st.session_state.get('gui_post_notes') or '').strip()
                if outp:
                    b.outputfile = outp
                return export_primary_h5(
                    b,
                    post_trial_notes=notes_in if notes_in else None,
                    outputfile=outp or None,
                    append_post_trial_notes=bool(st.session_state.get('gui_qc_notes_append', True)),
                )

            def _save_qc_plot_only():
                qix = _read_qc_trial_index(b)
                sel_h5 = str(st.session_state.get('gui_review_selected') or '').strip()
                qc_base = _qc_figure_base_path(b, sel_h5, qix)
                return save_universal_qc_figure(b, qc_trial_index=qix, base_path=qc_base)

            e1, e2 = st.columns(2)
            with e1:
                if _load_save_button('Only save Data File (.h5)', key='gui_save_h5_only', button_type='secondary'):
                    try:
                        rep = _export_h5_from_session()
                        if rep is not None:
                            st.success(f"{rep['message']}  →  `{rep['outputfile']}`")
                            if bool(st.session_state.get('gui_qc_notes_append', True)):
                                st.session_state['gui_post_notes'] = ''
                            else:
                                st.session_state['gui_post_notes'] = str(rep.get('post_trial_notes') or '')
                    except Exception as e:
                        _show_friendly_error(e, action='save_h5')
            with e2:
                if _load_save_button('Only Save QC Plot', key='gui_save_qc_only', button_type='secondary'):
                    try:
                        qc_path, _ = _save_qc_plot_only()
                        st.success(f'QC plot saved: `{qc_path}`')
                    except Exception as e:
                        _show_friendly_error(e, action='save_qc')

            if _load_save_button('Save Data File (.h5) and QC Plot', key='gui_save_h5_and_qc'):
                try:
                    with st.spinner('Writing data file (.h5)…'):
                        rep = _export_h5_from_session()
                    if rep is None:
                        pass
                    else:
                        qix = _read_qc_trial_index(b)
                        sel_h5 = str(st.session_state.get('gui_review_selected') or '').strip()
                        qc_base = _qc_figure_base_path(b, sel_h5, qix)
                        try:
                            with st.spinner('Saving QC plot…'):
                                qc_path, _ = save_universal_qc_figure(b, qc_trial_index=qix, base_path=qc_base)
                        except Exception as e:
                            _show_friendly_error(e, action='save_qc')
                            st.warning(f"The data file was saved: `{rep['outputfile']}`")
                        else:
                            st.success('Data file and QC plot saved.')
                            st.info(f"Data file: `{rep['outputfile']}`  |  QC plot: `{qc_path}`")
                        if bool(st.session_state.get('gui_qc_notes_append', True)):
                            st.session_state['gui_post_notes'] = ''
                        else:
                            st.session_state['gui_post_notes'] = str(rep.get('post_trial_notes') or '')
                except Exception as e:
                    _show_friendly_error(e, action='save_h5')

        st.checkbox(
            'Hide section (values stay; unhide to edit)',
            key='gui_sec6_hide',
            help='Collapse save buttons when finished.',
        )

    if _show_sec7_and_8():

        st.divider()
        st.session_state.setdefault('gui_sec7_hide', False)
        st.subheader('9 · Review data')
        if st.session_state.get('gui_sec7_hide'):
            st.caption('Visualization panel hidden. Uncheck **Hide section** below.')
        if not st.session_state.get('gui_sec7_hide'):
            data_path = _output_path_anchor_for_review(b)
            review_files = _candidate_review_files(data_path) if data_path else []
            if not review_files:
                st.info('No matching files in the data folder yet. Set **Data folder** and **Data file name** in Step 2, then run or export.')
            else:
                selected_file = st.session_state.get('gui_review_selected')
                if selected_file not in review_files:
                    st.session_state['gui_review_selected'] = review_files[0]
                    selected_file = review_files[0]
                st.caption(f'**Selected file:** `{os.path.basename(selected_file)}`')
                ext = os.path.splitext(str(selected_file).lower())[1]
                if ext in ('.png', '.jpg', '.jpeg', '.webp'):
                    st.image(selected_file, caption=os.path.basename(selected_file))
                else:
                    st.caption(f"Selected file: `{selected_file}`")
    
                st.markdown('**Custom plots from data file**')
                if ext != '.h5':
                    st.info('Select a **`.h5`** file above to plot saved time series (torque, angle, stim, etc.).')
                else:
                    summ = h5_custom_plot_summary(selected_file)
                    if not summ['ok']:
                        st.warning(summ.get('error') or 'Could not read this HDF5 file.')
                    else:
                        m1, m2, m3, m4 = st.columns(4)
                        with m1:
                            st.metric('test_type (file)', summ['test_type'] or '—')
                        with m2:
                            st.metric('Trials', str(summ['n_trials']))
                        with m3:
                            st.metric('schema', summ['schema_version'] or '—')
                        with m4:
                            ax = summ.get('primary_bending_axis') or ''
                            st.metric('primary axis (saved)', ax or '—')
    
                        trial_names = summ['trial_names']
                        if not trial_names:
                            st.info('No `trial_*` groups found under `02_TimeSeries`.')
                        else:
                            def_trial = trial_names[0]
                            if st.session_state.get('gui_h5_trial_select') not in trial_names:
                                st.session_state['gui_h5_trial_select'] = def_trial
                            trial_sel = st.selectbox(
                                'Trial',
                                trial_names,
                                key='gui_h5_trial_select',
                                help='Each trial is one acquisition segment stored in the file.',
                            )
                            plot_vars = list_h5_plot_variables(
                                selected_file,
                                trial_sel,
                                channel_names=summ.get('channel_names'),
                            )
                            if not plot_vars:
                                st.info('No plottable numeric series in this trial (need length ≥ 2).')
                            else:
                                ids = [v['id'] for v in plot_vars]
                                id_to_label = {v['id']: v['label'] for v in plot_vars}
                                tt0 = plot_vars[0].get('trial_test_type') or ''
                                if tt0:
                                    st.caption(f'**test_type (this trial):** `{tt0}`')
    
                                if 'gui_h5_n_panels' not in st.session_state:
                                    st.session_state['gui_h5_n_panels'] = 0.0
                                n_panel = st.number_input(
                                    'Number of figure panels',
                                    min_value=1,
                                    max_value=4,
                                    value=1,
                                    step=1,
                                    key='gui_h5_n_panels',
                                    help='Each panel is a separate Plotly figure. Use multiple Y traces in one panel for several lines on the same axes.',
                                )
                                for p in range(int(n_panel)):
                                    st.markdown(f'**Panel {p + 1}**')
                                    cxa, cya = st.columns(2)
                                    with cxa:
                                        st.selectbox(
                                            'X (time, angle, stim, …)',
                                            ids,
                                            key=f'gui_h5_x_{p}',
                                            format_func=lambda i: id_to_label.get(i, i),
                                        )
                                    with cya:
                                        st.multiselect(
                                            'Y — one or more traces (same plot)',
                                            ids,
                                            key=f'gui_h5_y_{p}',
                                            format_func=lambda i: id_to_label.get(i, i),
                                            help='Pick several series to overlay on the same axes (e.g. raw vs corrected torque).',
                                        )
    
                                if st.button('Generate plots from file', key='gui_h5_gen_plots'):
                                    _any_panel_ok = False
                                    for p in range(int(n_panel)):
                                        x_id = st.session_state.get(f'gui_h5_x_{p}')
                                        y_ids = st.session_state.get(f'gui_h5_y_{p}') or []
                                        if not x_id or not y_ids:
                                            _st_error_actions(
                                                f'Panel {p + 1} incomplete.',
                                                ['Choose X axis', 'Choose one or more Y'],
                                            )
                                            continue
                                        try:
                                            x_raw = read_h5_series(selected_file, trial_sel, x_id)
                                            fig_p = go.Figure()
                                            x_label = id_to_label.get(x_id, x_id)
                                            for y_id in y_ids:
                                                y_raw = read_h5_series(selected_file, trial_sel, y_id)
                                                x, y = align_xy(x_raw, y_raw)
                                                y_label = id_to_label.get(y_id, y_id)
                                                fig_p.add_trace(
                                                    go.Scatter(x=x, y=y, mode='lines', name=y_label),
                                                )
                                            fig_p.update_layout(
                                                title=f'Panel {p + 1}: {os.path.basename(selected_file)} · {trial_sel}',
                                                xaxis_title=x_label,
                                                yaxis_title='Y (multiple traces)',
                                                height=420,
                                                legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='right', x=1),
                                            )
                                            st.plotly_chart(fig_p, use_container_width=True)
                                            _any_panel_ok = True
                                        except Exception as e:
                                            _st_error_actions(
                                                f'Panel {p + 1} plot failed.',
                                                ['Check trial selection', 'Open Technical details'],
                                            )
                                            with st.expander('Technical details'):
                                                st.exception(e)
                                    if _any_panel_ok:
                                        _mark_review_data_used()
    
                st.markdown('**Add post-experiment note to selected file**')
                note_file = st.text_area(
                    'New note text (appended for .h5 files)',
                    key='gui_selected_file_note',
                    height=90,
                )
                if _load_save_button('Append note to selected file', key='gui_append_note_to_h5'):
                    try:
                        if ext != '.h5':
                            raise ValueError('Selected file is not .h5. Choose a data file to append notes.')
                        _append_note_to_h5_file(selected_file, note_file)
                        _mark_review_data_used()
                        st.success('Note appended to selected data file.')
                    except Exception as e:
                        _show_friendly_error(e, action='save_h5')

        st.checkbox(
            'Hide section (values stay; unhide to edit)',
            key='gui_sec7_hide',
            help='Collapse plots and file preview when finished.',
        )

        st.divider()
        st.session_state.setdefault('gui_sec8_hide', False)
        st.subheader('9 · Review data notes')
        st.caption(
            'Optional notes (specimen, setup, data quality). Pick a file below; notes append to the chosen `.h5`. '
            'QC plot export may use **kaleido** for PNG (`pip install kaleido`); otherwise HTML.'
        )
        if st.session_state.get('gui_sec8_hide'):
            st.caption('Note controls hidden. Uncheck **Hide section** below.')
        if not st.session_state.get('gui_sec8_hide'):
            data_path_qc = _output_path_anchor_for_review(b)
            review_files_qc = _candidate_review_files(data_path_qc) if data_path_qc else []
            if not review_files_qc:
                st.info('Set **Data folder** and **Data file name** in Step 2 (or save once) so files appear here.')
            else:
                _ensure_review_file_selection(review_files_qc)
                st.selectbox(
                    'Data folder file',
                    review_files_qc,
                    key='gui_review_selected',
                    format_func=lambda fp: os.path.basename(fp),
                    help='Files in your data folder. An `.h5` selection sets the base name for the next QC plot export.',
                )

            tr_qc = list(getattr(b, 'trial_records', []) or [])
            n_tr = len(tr_qc)
            if n_tr > 1:
                opts_ix = list(range(n_tr))
                _cur_ix = st.session_state.get('gui_qc_trial_index')
                if _cur_ix not in opts_ix:
                    st.session_state['gui_qc_trial_index'] = opts_ix[-1]

                def _qc_trial_label(i: int) -> str:
                    r = tr_qc[i]
                    stp = r.get('step_index', r.get('trial_index', r.get('cycle_index')))
                    ttp = str(r.get('test_type', '') or '')
                    parts = [f'Trial {i}']
                    if ttp:
                        parts.append(ttp)
                    if stp is not None:
                        parts.append(f'step {stp}')
                    return ' — '.join(parts) if len(parts) > 1 else parts[0]

                st.selectbox(
                    'In-memory trial for QC plot (this session)',
                    opts_ix,
                    format_func=_qc_trial_label,
                    key='gui_qc_trial_index',
                    help=(
                        'QC plot uses **trial_records** in memory after Run or save. '
                        'Pick which segment when the protocol produced several (e.g. multiple isometric steps).'
                    ),
                )
            else:
                st.caption(
                    '**QC trial:** only one in-memory segment — no trial choice needed. '
                    '(After multi-step runs, a trial list appears here.)'
                )

            st.markdown('**Note text for next save (stored on `.h5`)**')
            st.text_area(
                'Comment for export (optional)',
                key='gui_post_notes',
                height=100,
                help=(
                    'On save/export, this text is **appended** to any existing note with a timestamp. Leave empty to skip. '
                    'Uncheck **Append…** below to replace the stored note entirely.'
                ),
            )
            if 'gui_qc_notes_append' not in st.session_state:
                st.session_state['gui_qc_notes_append'] = False
            st.checkbox(
                'Append new text to existing note in this data file path (when the file already exists)',
                value=True,
                key='gui_qc_notes_append',
            )

        st.checkbox(
            'Hide section (values stay; unhide to edit)',
            key='gui_sec8_hide',
            help='Collapse file picker and note fields when finished.',
        )

    _autosave_tick()


if __name__ == '__main__':
    main()
