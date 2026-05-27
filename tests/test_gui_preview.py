"""Preview builder tests (stim + native units, no DAQ)."""
from pathlib import Path
import sys
from unittest.mock import MagicMock

import numpy as np
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

import bender_gui_preview as pv


class _PreviewBender:
    specimen_side_index_left = 0
    specimen_side_index_right = 1
    stim_pulse_rate = 75.0
    daq_ai_sample_rate_hz = 1000.0
    recruitment = 'bilateral_simultaneous'
    bilateral_sequential_left_frac = 0.5
    bilateral_mirror_motor = False

    def __init__(self):
        self.dclamp = 10.0
        self.xsec_width = 2.0
        self.isometric_initial = 0.02
        self.isometric_final = 0.06
        self.isometric_num_steps = 3
        self.isometric_mode = 'strain'
        self.isovelocity_min_vel = 1.0
        self.isovelocity_max_vel = 5.0
        self.isovelocity_starting_strain = 0.03
        self.isovelocity_starting_strain_mode = 'strain'
        self.isovelocity_velocity_mode = 'angle_vel'
        self.isovelocity_num_steps = 2
        self.isovelocity_iso_duration_s = 0.2
        self.isovelocity_pre_hold_s = 0.1
        self.isometric_stim_params = {'is_stim': True, 'stim_voltage': 5.0}
        self.isovelocity_stim_params = {'is_stim': True, 'stim_voltage': 5.0}

    def _effective_dclamp_mm(self):
        return float(self.dclamp)

    def _normalize_recruitment(self, r):
        return str(r)

    def _recruitment_with_bilateral_mirror_motor(self, r, _bm):
        return r

    def recruitment_unilateral_lateral_index(self, _rec):
        return None

    def motor_command_sign_for_bend_toward_index(self, _idx):
        return 1.0

    def _stim_params_with_lateral(self, sp):
        return dict(sp)

    def _timeline_ramp_hold(self, a0, a1, ramp, hold, hz):
        n = max(4, int((ramp + hold) * hz))
        t = np.linspace(0, ramp + hold, n)
        ang = np.linspace(a0, a1, n)
        w = np.gradient(ang, t)
        return t, ang, w

    def _timeline_mirror_two_holds(self, *args, **kwargs):
        t, ang, w = self._timeline_ramp_hold(0, 1, 0.1, 0.2, 1000)
        return t, ang, w, (0.1, 0.2), (0.2, 0.3)

    def _pulse_carrier_volts(self, t, mask, spr, volt):
        m = np.asarray(mask, dtype=bool)
        out = np.zeros_like(np.asarray(t, dtype=float))
        if np.any(m):
            out[m] = float(volt)
        return out

    def _route_recruitment_stim(self, pulse, _rec, sequential_left_frac=0.5):
        p = np.asarray(pulse, dtype=float)
        return p, np.zeros_like(p)

    def _deposit_stim_on_side(self, pulse, side, s1, s2):
        pass

    def _isovelocity_one_block(self, theta0, vel, **kw):
        n = 50
        t = np.linspace(0, 0.5, n)
        ang = theta0 + vel * t
        w = np.full(n, vel)
        active = (t >= 0.2) & (t < 0.45)
        pulse = self._pulse_carrier_volts(t, active, 75.0, kw.get('stim_voltage', 5.0))
        s1, s2 = self._route_recruitment_stim(pulse, kw.get('recruitment', 'bilateral_simultaneous'))
        return {'t': t, 'angle': ang, 'anglevel': w, 's1': s1, 's2': s2}


def test_isometric_preview_includes_stim_and_native_units():
    b = _PreviewBender()
    out = pv.build_protocol_preview(b, requested_test_type='isometric', max_plot_points=200)
    assert out.get('ok') is True
    assert out.get('stim_s1') is not None
    assert out.get('stim_plot') is not None
    assert out.get('strain') is not None or out.get('curvature') is not None


def test_isovelocity_preview_includes_stim():
    b = _PreviewBender()
    out = pv.build_protocol_preview(b, requested_test_type='isovelocity', max_plot_points=200)
    assert out.get('ok') is True
    assert out.get('stim_total') is not None or out.get('stim_plot') is not None


def test_isometric_preview_errors_without_clamp():
    b = _PreviewBender()
    b.dclamp = None
    b._effective_dclamp_mm = lambda: None
    out = pv.build_protocol_preview(b, requested_test_type='isometric')
    assert out.get('ok') is False
    assert out.get('error')
