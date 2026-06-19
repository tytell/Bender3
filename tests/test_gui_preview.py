"""Preview builder tests (stim + native units, no DAQ)."""
from pathlib import Path
import sys
from unittest.mock import MagicMock

import numpy as np
import pytest

sys.path.append(str(Path(__file__).resolve().parents[1]))

import bender_gui_preview as pv
from bender_functions import Bender as _RealBender


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
        self.block_sequence = [{'direction': 'left', 'stim_sides': 'left'}]
        self.left_stim_voltage = 5.0
        self.right_stim_voltage = 5.0
        self.isometric_stim_params = {'is_stim': True, 'settle_before_stim_s': 0.5}
        self.isovelocity_stim_params = {'is_stim': True, 'settle_before_stim_s': 0.02}

    def _effective_dclamp_mm(self):
        return float(self.dclamp)

    def _normalize_recruitment(self, r):
        return str(r)

    def _normalize_block_sequence(self, block_sequence):
        if not block_sequence:
            return None
        return block_sequence

    def _lateral_index_for_block_direction(self, direction):
        return 0 if str(direction).lower() == 'left' else 1

    def _route_stim_sides_volts(self, t, active, spr, stim_sides, left_v, right_v):
        p = self._pulse_carrier_volts(t, active, spr, left_v)
        s1 = np.zeros_like(np.asarray(t, dtype=float))
        s2 = np.zeros_like(s1)
        sides = str(stim_sides).lower()
        if sides in ('left', 'both'):
            s1 = np.where(active, float(left_v), 0.0)
        if sides in ('right', 'both'):
            s2 = np.where(active, float(right_v), 0.0)
        return s1, s2

    def _recruitment_with_bilateral_mirror_motor(self, r, _bm):
        return r

    def recruitment_unilateral_lateral_index(self, _rec):
        return None

    def motor_command_sign_for_bend_toward_index(self, _idx):
        return 1.0

    def _stim_params_with_lateral(self, sp):
        return dict(sp)

    def _resolve_stim_onset_duration_s(self, sp, *, segment_duration_s):
        return _RealBender._resolve_stim_onset_duration_s(self, sp, segment_duration_s=segment_duration_s)

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

    def _neutral_reset_ramp_duration_s(self, from_deg):
        max_speed = 15.0
        ramp = abs(float(from_deg)) / max_speed
        return max(ramp, 2.0 / float(self.daq_ai_sample_rate_hz))

    def _build_isometric_one_step(self, target_deg, *, prev_deg=0.0, ramp_duration_s,
                                  hold_duration_s, pre_baseline_s=0.0, post_baseline_s=0.0,
                                  is_stim=False, spr=75.0, stim_voltage=5.0, daq_hz=1000.0,
                                  post_buffer_s=1.0, **kw):
        # Stub mirroring the run-path single-step shape: ramp 0->target, baselines+hold at target,
        # speed-capped return to home, trailing post_buffer hold at 0. Start/end at home.
        total_hold = float(pre_baseline_s) + float(hold_duration_s) + float(post_baseline_s)
        t1, a1, w1 = self._timeline_ramp_hold(
            float(prev_deg), float(target_deg), float(ramp_duration_s), total_hold, daq_hz
        )
        ret_s = self._neutral_reset_ramp_duration_s(float(target_deg))
        t2, a2, w2 = self._timeline_ramp_hold(float(target_deg), 0.0, ret_s, float(post_buffer_s), daq_hz)
        t = np.concatenate([t1, float(t1[-1]) + np.asarray(t2, dtype=float)[1:]])
        ang = np.concatenate([a1, np.asarray(a2, dtype=float)[1:]])
        w = np.concatenate([w1, np.asarray(w2, dtype=float)[1:]])
        stim_start = float(ramp_duration_s) + float(pre_baseline_s)
        active = (t >= stim_start) & (t < stim_start + float(hold_duration_s))
        if is_stim and np.any(active):
            pulse = self._pulse_carrier_volts(t, active, spr, stim_voltage)
            s1, s2 = self._route_recruitment_stim(pulse, kw.get('recruitment', 'bilateral_simultaneous'))
        else:
            s1 = np.zeros_like(t)
            s2 = np.zeros_like(t)
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


def test_isometric_preview_baselines_extend_timeline():
    # Issue 12: pre/post baselines should lengthen each isometric step's ramp+hold segment.
    b = _PreviewBender()
    targets = np.array([5.0, 10.0])
    sp_none = pv._isometric_stim_params_from_b(b)
    sp_none.update({'ramp_duration_s': 2.0, 'hold_duration_s': 5.0,
                    'pre_baseline_s': 0.0, 'post_baseline_s': 0.0, 'inter_step_interval_s': 0.0})
    t0, *_ = pv._preview_concat_isometric_timeline(b, targets, sp_none, mode='strain')

    sp_base = dict(sp_none)
    sp_base.update({'pre_baseline_s': 1.5, 'post_baseline_s': 2.0})
    t1, *_ = pv._preview_concat_isometric_timeline(b, targets, sp_base, mode='strain')

    # Two steps each gain (1.5 + 2.0) s of silent baseline hold.
    assert float(t1[-1]) == pytest.approx(float(t0[-1]) + 2 * (1.5 + 2.0), rel=0.05)


def test_isometric_preview_errors_without_clamp():
    b = _PreviewBender()
    b.dclamp = None
    b._effective_dclamp_mm = lambda: None
    out = pv.build_protocol_preview(b, requested_test_type='isometric')
    assert out.get('ok') is False
    assert out.get('error')


def _dynamic_stim_bender():
    b = _RealBender('jimenez_bender_config_A')
    b.dclamp = 10.0
    b.xsec_width = 2.0
    b.all_freqs = [1.0]
    b.all_amps = [0.05]
    b.all_amps_mode = 'strain'
    b.cycles_per_step = 3
    b.n_end_cycles = 0
    b.randomize = False
    b.stim_cycles_in_step = []
    b.stim_pulse_rate = 75.0
    b.is_stim = True
    # Deterministic conditioning-pulse timing for the assertions below.
    b.prestim_time = -2.0
    b.poststim_time = 2.0
    b.prepoststim_dur = 0.06
    b.prepoststim_sep = 1.0
    return b


def _pulse_window_starts(t, sig):
    t = np.asarray(t, dtype=float)
    on = np.asarray(sig, dtype=float) > 0
    idx = np.where(on)[0]
    starts = []
    for i in idx:
        if i == 0 or not on[i - 1]:
            starts.append(float(t[i]))
    return starts


def test_dynamic_preview_post_stim_relative_to_active_window():
    """Post-conditioning pulses must land after the active motion ends, not at a fixed
    absolute time inside it (regression: make_stimuli movedur not passed in preview)."""
    b = _dynamic_stim_bender()
    out = pv.build_protocol_preview(b, requested_test_type='dynamic', max_plot_points=100000)
    assert out.get('ok') is True
    t = np.asarray(out['t'], dtype=float)
    s1 = np.asarray(out['stim_s1'], dtype=float)
    s2 = np.asarray(out['stim_s2'], dtype=float)
    movedur = float(np.sum(np.asarray(b.period_by_cycle, dtype=float)))

    s1_starts = _pulse_window_starts(t, s1)
    s2_starts = _pulse_window_starts(t, s2)
    assert s1_starts and s2_starts

    # Each conditioning burst is a train of carrier pulses spanning prepoststim_dur; the burst
    # ONSET is the first pulse of the relevant group. Pre-stim fires before motion (explicit
    # negative onset); post-stim fires after the active motion window ends.
    s1_pre = [s for s in s1_starts if s < 0.0]
    s1_post = [s for s in s1_starts if s >= movedur - 1e-6]
    s2_post = [s for s in s2_starts if s >= movedur - 1e-6]
    assert s1_pre and s1_post and s2_post

    assert min(s1_pre) == pytest.approx(b.prestim_time, abs=2e-2)
    assert min(s1_post) == pytest.approx(movedur + b.poststim_time, abs=2e-2)
    assert min(s2_post) == pytest.approx(movedur + b.poststim_time + b.prepoststim_sep, abs=2e-2)

    # With no in-cycle stim configured, nothing should fire strictly inside the motion window.
    in_motion = [s for s in (s1_starts + s2_starts) if 1e-6 < s < movedur - 1e-6]
    assert in_motion == []
