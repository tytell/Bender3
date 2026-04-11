"""
CritterGripper — NiceGUI application entry.

Run from the project directory:
    python main.py

Full workflow lives in ``bender_nicegui.workflow`` (sections 1–6, H5 explorer, landing).
"""
from __future__ import annotations

import os
import sys

_ROOT = os.path.dirname(os.path.abspath(__file__))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from bender_nicegui.workflow import run_app  # noqa: E402

if __name__ in {'__main__', '__mp_main__'}:
    run_app(port=8765)
