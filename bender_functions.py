import numpy as np
import importlib
from scipy import interpolate
from datetime import datetime
from copy import copy
import time
import re
import os
import sys
import xml.etree.ElementTree as ElementTree
import json
print(f"DEBUG: Loading bender_functions.py from: {os.path.abspath(__file__)}")
import logging

# Hardware config modules live in ``templates/configs/`` and are imported by bare module
# name (e.g. ``importlib.import_module('jimenez_bender_config_A')``). Put that folder on
# ``sys.path`` at import time so both the app and the test suite resolve them after the
# templates reorg, without changing how callers reference configs by stem.
_CONFIGS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'templates', 'configs')
if os.path.isdir(_CONFIGS_DIR) and _CONFIGS_DIR not in sys.path:
    sys.path.insert(0, _CONFIGS_DIR)

try:
    import nidaqmx.constants as daq
    from nidaqmx import Task
    from nidaqmx.stream_writers import AnalogMultiChannelWriter, DigitalSingleChannelWriter
    from nidaqmx.stream_readers import AnalogMultiChannelReader, CounterReader
    from nidaqmx.errors import DaqError
    from nidaqmx.constants import TerminalConfiguration
except ImportError:
    logging.warning('No DAQmx available')
    Task = None  # type: ignore
    daq = None  # type: ignore

import xml.etree.ElementTree as ElementTree

import h5py


# --- Isovelocity safety limit ------------------------------------------------
# Hard limit on commanded/encoder motor angle during an isovelocity ramp. If the ramp would carry
# the angle to or past this magnitude, the ramp is stopped and the motor is returned to 0° (2.1).
ISOVELOCITY_ANGLE_LIMIT_DEG = 45.0

# --- Specimen lateral frame (signed indices) ---------------------------------
# Live values come from config ``specimen_lateral_index_on_positive_motor_side`` (non-zero)
# paired with ``positive_motor_direction``; assigned on :class:`Bender` construction.
# Default -1 / +1 below match the default config and are only for offline reference.
SPECIMEN_SIDE_INDEX_LEFT = -1
SPECIMEN_SIDE_INDEX_RIGHT = 1


def _normalize_start_end(x):
    """Accept float or sequence [start, end]; return (start, end)."""
    if x is None:
        return None
    if np.isscalar(x):
        xf = float(x)
        return xf, xf
    arr = np.asarray(x, dtype=float).reshape(-1)
    if arr.size == 0:
        raise ValueError("nominal frequency/curvature must be non-empty")
    if arr.size == 1:
        v = float(arr[0])
        return v, v
    return float(arr[0]), float(arr[-1])


def _ramp_progress(u, ramp_mode, velocity_exponent):
    """Map normalized time u in [0,1] to progress in [0,1].

    NOTE: ``velocity_exponent`` is NOT ``amplitude_frequency_exponent``. This helper only
    warps *when* along a generic ramp (ease-in/out of normalized time). It does not set
    how bend amplitude scales with frequency during a log-frequency sweep; see
    :meth:`Bender.make_cycles_frequency_sweep` for that.
    """
    u = np.clip(np.asarray(u, dtype=float), 0.0, 1.0)
    if ramp_mode is None or ramp_mode == 'linear':
        return u
    if ramp_mode == 'exponential':
        exp_x = float(velocity_exponent)
        if abs(exp_x) < 1e-12:
            return u
        den = np.exp(np.clip(exp_x, -50.0, 50.0)) - 1.0
        if abs(den) < 1e-12:
            return u
        exu = np.clip(exp_x * u, -50.0, 50.0)
        return (np.exp(exu) - 1.0) / den
    raise ValueError(f"Unknown ramp_mode: {ramp_mode!r}")


def _scalar_or_pair(f0, f1):
    """Single float if endpoints match, else [start, end] (Hz or 1/m)."""
    a, b = float(f0), float(f1)
    return a if abs(b - a) < 1e-12 else [a, b]


def _strain_lever_arm_m(xsec_width_mm, target_muscle_depth_mm=0.0) -> float:
    """
    Effective lever arm (m) for engineering strain fraction: ε = κ * lever_arm.

    ``lever_arm = (xsec_width_mm/2 - target_muscle_depth_mm) / 1000`` with ``target_muscle_depth_mm=0``
    giving surface strain (legacy behavior).
    """
    xw = float(xsec_width_mm)
    md = float(target_muscle_depth_mm or 0.0)
    if not np.isfinite(xw) or xw <= 0:
        raise ValueError(f"xsec_width_mm must be a finite value > 0 mm; got {xsec_width_mm!r}.")
    if not np.isfinite(md) or md < 0:
        raise ValueError(f"target_muscle_depth_mm must be a finite value >= 0 mm; got {target_muscle_depth_mm!r}.")
    half_w_mm = xw / 2.0
    if md >= half_w_mm:
        raise ValueError(
            f"target_muscle_depth_mm ({md} mm) must be < xsec_width_mm/2 ({half_w_mm} mm) "
            f"so effective lever arm stays positive (xsec_width_mm={xw})."
        )
    return (half_w_mm - md) / 1000.0


def curvature_to_strain_fraction(kappa, *, xsec_width_mm, target_muscle_depth_mm=0.0):
    """Inverse of strain modes in :func:`convert_to_curvature`: ε = κ * lever_arm."""
    lever_m = _strain_lever_arm_m(xsec_width_mm, target_muscle_depth_mm)
    return np.asarray(kappa, dtype=float) * lever_m


def convert_to_curvature(value, mode, *, dclamp_mm=None, xsec_width_mm=None, target_muscle_depth_mm=None):
    """
    Convert specimen inputs to curvature κ (1/m) or curvature rate dκ/dt (1/m/s).

    Position modes
    --------------
    'curvature'
        ``value`` is κ (1/m); returned unchanged.
    'strain'
        ``value`` is **engineering strain as a fraction** (0.05 = 5 %). Same relation as
        ``organize_cycles``: strain = κ * lever_arm_m with
        lever_arm_m = (xsec_width_mm/2 - target_muscle_depth_mm)/1000.

        **Interpretation (bending):** this scalar is κ times distance from the neutral axis to
        the fiber layer (``target_muscle_depth_mm=0`` → **surface** half-width). Linear bending gives **equal magnitude** tensile strain on
        one surface and compressive (shortening) strain on the other; the sign of κ (and thus
        ε) follows the commanded bend direction. It is **not** automatically “ipsilateral
        muscle shortening”; map sign to left/right fibers using anatomy + mounting, or flip
        displayed traces with :attr:`Bender.strain_shortening_positive_display_sign`.
    'strain_pct'
        ``value`` is strain in **percent** (5 = 5 %).
    'angle'
        ``value`` is commanded motor / specimen bend angle in **degrees**; uses
        κ = rad(angle) * 1000 / dclamp_mm (inverse of ``rad2deg(κ * dclamp/1000)``).

    Rate modes (curvature velocity)
    -------------------------------
    These are **time derivatives** (e.g. during a linear ramp), not peak rates in cyclic bending.
    For sinusoidal motion at frequency ``f`` Hz, peak |dε/dt| is ``2π f ε_max`` (strain fraction);
    see :func:`peak_strain_rate_sinusoidal` — same relation as ``organize_cycles`` strain rates.

    'curvature_rate' — dκ/dt in (1/m)/s.
    'strain_rate' — d(strain_fraction)/dt in 1/s.
    'strain_pct_rate' — d(percent strain)/dt in %/s.
    'angle_vel' — angular rate in deg/s.
    """
    mode = str(mode).lower().strip()
    v = np.asarray(value, dtype=float)

    if mode == 'curvature':
        return v

    md = 0.0 if target_muscle_depth_mm is None else float(target_muscle_depth_mm)

    if mode in ('strain', 'strain_rate'):
        if xsec_width_mm is None:
            raise ValueError(f"convert_to_curvature(mode={mode!r}) requires xsec_width_mm (mm).")
        lever_m = _strain_lever_arm_m(xsec_width_mm, md)
        return v / lever_m

    if mode in ('strain_pct', 'strain_pct_rate'):
        if xsec_width_mm is None:
            raise ValueError(f"convert_to_curvature(mode={mode!r}) requires xsec_width_mm (mm).")
        lever_m = _strain_lever_arm_m(xsec_width_mm, md)
        frac = v / 100.0
        return frac / lever_m

    if mode in ('angle', 'angle_vel'):
        if dclamp_mm is None:
            raise ValueError(f"convert_to_curvature(mode={mode!r}) requires dclamp_mm (mm).")
        dc = float(dclamp_mm)
        return np.deg2rad(v) * 1000.0 / dc

    if mode == 'curvature_rate':
        return v

    raise ValueError(
        f"Unknown convert_to_curvature mode {mode!r}. "
        "Use curvature, strain, strain_pct, angle, or *_rate / angle_vel / curvature_rate."
    )


def peak_strain_rate_sinusoidal(strain_fraction, frequency_hz):
    """
    Peak |dε/dt| for ε(t) = ε_max sin(2π f t) with ε_max = ``strain_fraction``.

    Matches ``organize_cycles``: ``strainrate = 2 * π * strain * freq``. Use this when you care
    about cyclic strain rate at a chosen **frequency**; it is not the same as a constant
    ``strain_rate`` passed to :func:`convert_to_curvature` (which is dε/dt for a ramp, in 1/s).
    """
    f = np.asarray(frequency_hz, dtype=float)
    e = np.asarray(strain_fraction, dtype=float)
    return 2.0 * np.pi * f * e


class MasterLogger:
    """Accumulates protocol metadata for HDF5 export (merge into saver attrs/datasets)."""

    def __init__(self):
        self.entries = {}

    def record(self, **kwargs):
        for k, v in kwargs.items():
            if isinstance(v, (np.integer, np.floating)):
                v = float(v) if getattr(v, 'shape', ()) == () else np.asarray(v).tolist()
            elif isinstance(v, np.ndarray):
                v = v.tolist()
            elif isinstance(v, (list, tuple)) and len(v) == 2 and all(np.isscalar(s) for s in v):
                v = [float(v[0]), float(v[1])]
            self.entries[k] = v

    def clear(self):
        self.entries = {}

    def as_dict(self):
        return dict(self.entries)


class Bender:
    def __init__(self, config_module_name):
        # 1. Dynamically import the python config file
        # This turns 'bender_config_A' into an object called 'cfg'
        try:
            cfg = importlib.import_module(config_module_name)
        except ImportError:
            raise ImportError(f"Could not find {config_module_name}.py in the current folder.")

        self.config_name = config_module_name


        # 2. Assign Hardware settings from cfg
        self.device_name = cfg.device_name
        self.motor_port = cfg.motor_port 
        self.encoder_chan = cfg.encoder_chan
        self.stim_channels = cfg.stim_channels  
        self.S1stim_chan, self.S2stim_chan = self.stim_channels # Map names for run_experiment and run
        # Canonical AI acquisition rate from the hardware config. The instance
        # daq_ai_sample_rate_hz (a validated property) is checked against this so a
        # stale/garbage value (e.g. 0.5 Hz) can never under-sample dynamic sine
        # commands into a ramp and make the motor "walk".
        self._config_daq_ai_sample_rate_hz = float(cfg.daq_ai_sample_rate_hz)
        self.daq_ai_sample_rate_hz = cfg.daq_ai_sample_rate_hz
        self.daq_ao_do_sample_rate_hz = cfg.daq_ao_do_sample_rate_hz
        self.motor_gear_ratio = cfg.motor_gear_ratio
        self.motor_full_steps_per_rev = cfg.motor_full_steps_per_rev
        self.encoder_pulses_per_rev = cfg.encoder_pulses_per_rev
          
          
        # Grab the sensor lists
        self._cfg_input_channels = list(getattr(cfg, 'input_channels', []))
        self._cfg_input_channel_names = list(getattr(cfg, 'input_channel_names', []))
        self.sg_channels = list(getattr(cfg, 'SG_chan', []))
        self.sg_names = list(getattr(cfg, 'SG_name', []))
        self.sono_channels = list(getattr(cfg, 'sono_channel', []))
        self.sono_names = list(getattr(cfg, 'sono_name', []))
        self.input_channels = list(self._cfg_input_channels)
        self.input_channel_names = list(self._cfg_input_channel_names)
        self.use_sono = bool(getattr(cfg, 'use_sono', True))
        if not self.use_sono:
            # Defensive filter in case config still includes sono channels.
            pairs = [
                (ch, nm) for ch, nm in zip(self.input_channels, self.input_channel_names)
                if not str(nm).lower().startswith('sono_')
            ]
            self.input_channels = [p[0] for p in pairs]
            self.input_channel_names = [p[1] for p in pairs]

        # 3. Assign Calibration/Directionality from cfg
        self.forcetorque_calibration_file = cfg.forcetorque_calibration_file
        self.positive_motor_direction = cfg.positive_motor_direction
        anchor = getattr(cfg, 'specimen_lateral_index_on_positive_motor_side', None)
        if anchor is None:
            anchor = -1
        self.specimen_lateral_index_on_positive_motor_side = int(anchor)
        if self.specimen_lateral_index_on_positive_motor_side == 0:
            raise ValueError(
                "config specimen_lateral_index_on_positive_motor_side must be non-zero: "
                "signed lateral index for the anatomical side named in positive_motor_direction "
                "(opposite side is the negated index)."
            )
        pm0 = str(self.positive_motor_direction).lower()
        a = self.specimen_lateral_index_on_positive_motor_side
        if pm0 == 'left':
            self.specimen_side_index_left = a
            self.specimen_side_index_right = -a
        elif pm0 == 'right':
            self.specimen_side_index_right = a
            self.specimen_side_index_left = -a
        else:
            logging.warning(
                "positive_motor_direction is not 'left' or 'right'; "
                "defaulting specimen_side_index_left=-1, specimen_side_index_right=1"
            )
            self.specimen_side_index_left = -1
            self.specimen_side_index_right = 1
        self.S1side = cfg.S1side
        self.S2side = cfg.S2side
        self.motor_axis = cfg.motor_axis
        self.bending_axis_sensor = cfg.bending_axis_sensor
        self.primary_bending_axis = getattr(cfg, 'primary_bending_axis', self.bending_axis_sensor)

        # Sonomicrometer specific calibration
        self.sono_cal_left = getattr(cfg, 'sono_cal_left', None)
        self.sono_cal_right = getattr(cfg, 'sono_cal_right', None)
        self.sono_internal_rate = getattr(cfg, 'sono_internal_samplefreq', 241)
        
        # 4. Assign Experiment Defaults from cfg
        self.waitbefore = cfg.waitbefore
        self.waitafter = cfg.waitafter
        self.rampdur = cfg.rampdur
        self.prepoststim_dur = cfg.prepoststim_dur
        self.prepoststim_sep = cfg.prepoststim_sep
        self.prestim_time = cfg.prestim_time
        self.poststim_time = cfg.poststim_time
        
         # 4. AUTO-CONFIGURE HARDWARE
        self.set_stim_channels(*self.stim_channels) 
        self.set_motor_channel(self.motor_port)
        self.set_encoder_channel(self.encoder_chan, encoder_pulses_per_rev=self.encoder_pulses_per_rev)

        # 5. NEW: AUTO-LOAD CALIBRATION
        # This uses the method you already wrote!
        try:
            self.loadCalibration(self.forcetorque_calibration_file)
        except Exception as e:
            print(f"⚠️ WARNING: Calibration failed to load: {e}")
        if not hasattr(self, 'calibration') or self.calibration is None:
            self.calibration = np.eye(6, dtype=float)

        # Simulation mode (set by GUI): :meth:`run` uses numpy instead of NI-DAQmx.
        self.simulation_mode = False
        self.simulation_material = 'polyurethane'

        # Last commanded motor angle (deg). Updated whenever a motion timeline is recorded; used by
        # command_start_position_zero() so the start-of-protocol re-zero ramps from the real last
        # commanded position rather than blindly assuming the motor is already at 0 (1.5).
        self._last_commanded_angle_deg = 0.0

        # 5. Placeholders (to prevent NoneType/0-channel errors)
        self.t = np.array([0.0, 1.0/self.daq_ai_sample_rate_hz])
        self.S1stimcmd = np.zeros(len(self.t))
        self.S2stimcmd = np.zeros(len(self.t))
        self.i_total_system = 0.0
        self.total_mass = 0.0
        
        self.angle = np.array([0.0, 0.0])
        self.anglevel = np.array([0.0, 0.0])
        self.tnorm = np.array([0.0, 0.0])
        self.amp_step_vel = cfg.amp_step_vel
        # NOTE: velocity_exponent vs amplitude_frequency_exponent — different physics:
        #   velocity_exponent feeds _ramp_progress (exponential ramp_mode): warps normalized
        #   time u∈[0,1]; does not appear in the frequency-sweep sine law.
        #   amplitude_frequency_exponent is only for make_cycles_frequency_sweep: amplitude
        #   scales as (f/f_start)^α while f sweeps (α=0 constant amp, α=-1 constant peak |ω|).
        self.velocity_exponent = 1.0
        self.ramp_mode_default = getattr(cfg, 'ramp_mode_default', 'linear')
        self.master_logger = MasterLogger()
        self.h5_protocol_metadata = {}
        self.h5_schema_version = "2.0"
        self.post_trial_notes = ""
        # GUI-friendly aliases for geometry naming
        self.test_segment_length_mm = None  # alias for dclamp
        self.test_segment_position_mm = None  # alias for dbend (metadata only)
        self.target_muscle_depth_mm = 0.0  # depth from outer surface to fiber layer (0 = surface strain)
        # Optional inertial-calibration linkage; safe defaults keep legacy runs working.
        self.inertial_calibration_file = None
        self.use_inertial_calibration = False
        self.inertial_calibration_profile = None
        self.inertial_torque_system_primary = np.array([])
        self.inertial_torque_specimen_primary = np.array([])
        self.inertial_torque_total_primary = np.array([])
        self.primary_torque_raw = np.array([])
        self.primary_torque_corrected = np.array([])
        self.trial_records = []
        # Sweep / step motion length (seconds); set via update_metadata(duration=...) or bender.duration = ...
        self.duration = None
        self.amplitude_frequency_exponent = None
        self.test_type = None  # set in notebook via bender.test_type + update_metadata
        # Recruitment / lateral stimulation for isometric & isovelocity (see _normalize_recruitment).
        self.recruitment = 'bilateral_simultaneous'
        self.lateral_mode = None  # optional alias; merged into stim_params as recruitment
        self.bilateral_mirror_motor = False
        self.bilateral_sequential_left_frac = 0.5
        # Optional native targets (same units as isometric_mode) for LEFT vs RIGHT hold when bilateral mirror is on.
        self.isometric_mirror_target_left = None
        self.isometric_mirror_target_right = None
        # Multiply geometric strain (κ·w/2) in plots/post-hoc by this after sono validation if
        # you want “positive = shortening on the recruited side” (typically ±1).
        self.strain_shortening_positive_display_sign = 1.0
        # Canonical rest (s) the motor holds position after each step finishes before the next
        # step begins, for the stepped protocols (FL = isometric, FV = isovelocity). 0 = back-to-back.
        self.rest_between_steps_s = 2.0
        # Canonical toggle: shuffle the step sequence (curvature levels for FL, velocity levels
        # for FV) before running. OFF = run in the generated linspace order. When a block
        # sequence is used, each block is shuffled independently.
        self.randomize_step_order = False
        # Canonical toggle: after each step's rest, drive the motor back to angle = 0° (resting
        # length) before the next step. OFF = motor stays at its current position. FL + FV only.
        self.reset_between_steps = False
        # Width (ms) of each stimulation pulse (the high time of every carrier pulse), applied to
        # all stimulated protocols. Combined with stim_pulse_rate this sets the pulse high fraction.
        self.pulse_width_ms = 2.0
        # Legacy isometric-only alias for rest_between_steps_s. None = not set; mirrored onto
        # rest_between_steps_s in _normalize_dispatch_aliases so old saved templates still apply.
        self.isometric_inter_step_interval_s = None
        # Block sequence (default: one block, bend LEFT / stim LEFT).
        self.block_sequence = [{'direction': 'left', 'stim_sides': 'left'}]
        self.left_stim_voltage = 5.0
        self.right_stim_voltage = 5.0
        self.block_reset_ramp_duration_s = 2.0

        # Standard 2D shapes for NI-DAQmx (Channels, Samples)
        self.stimcmdhi = np.zeros((2, 2))
        self.dig = np.zeros((1, 2), dtype='uint32')

        print(f"Bender initialized using: {config_module_name}.py")

    def _build_trial_record(self, test_type, trial_index, cycle_index=0, extra=None):
        """Create a standardized per-trial record for H5 export."""
        tt = str(test_type)
        rec = {
            'test_type': tt,
            'trial_index': int(trial_index),
            'cycle_index': int(cycle_index),
            't': np.array(getattr(self, 't', np.array([])), copy=True),
            'angle_cmd': np.array(getattr(self, 'angle', np.array([])), copy=True),
            'anglevel_cmd': np.array(getattr(self, 'anglevel', np.array([])), copy=True),
            'tnorm': np.array(getattr(self, 'tnorm', np.array([])), copy=True),
            'S1stimcmd': np.array(getattr(self, 'S1stimcmd', np.array([])), copy=True),
            'S2stimcmd': np.array(getattr(self, 'S2stimcmd', np.array([])), copy=True),
            'aidata': np.array(getattr(self, 'aidata', np.array([])), copy=True),
            'angle_measured': np.array(getattr(self, 'angle_measured', np.array([])), copy=True),
            'forcetorque': np.array(getattr(self, 'forcetorque', np.array([])), copy=True),
            'forcetorque_raw': np.array(getattr(self, 'forcetorque_raw', np.array([])), copy=True),
            'forcetorque_corrected': np.array(getattr(self, 'forcetorque_corrected', np.array([])), copy=True),
            'inertial_torque_system_primary': np.array(getattr(self, 'inertial_torque_system_primary', np.array([])), copy=True),
            'inertial_torque_specimen_primary': np.array(getattr(self, 'inertial_torque_specimen_primary', np.array([])), copy=True),
            'inertial_torque_total_primary': np.array(getattr(self, 'inertial_torque_total_primary', np.array([])), copy=True),
            'primary_torque_raw': np.array(getattr(self, 'primary_torque_raw', np.array([])), copy=True),
            'primary_torque_corrected': np.array(getattr(self, 'primary_torque_corrected', np.array([])), copy=True),
        }
        if extra:
            rec.update(extra)
        return rec

    def loadCalibration(self, calibrationFile):
        if not os.path.exists(calibrationFile):
            raise IOError("Calibration file %s not found", calibrationFile)

        try:
            tree = ElementTree.parse(calibrationFile)
            cal = tree.getroot().find('Calibration')
            if cal is None:
                raise IOError('Not a calibration XML file')

            mat = []
            for ax in cal.findall('UserAxis'):
                txt = ax.get('values')
                row = [float(v) for v in txt.split()]
                mat.append(row)

        except IOError:
          raise IOError(f"Calibration file {calibrationFile} not found")

        self.calibration = np.array(mat).T

    def apply_calibration_forcetorque(self, rawdata):
        self.forcetorque = np.dot(rawdata[:6, :].T, self.calibration)
        self.forcetorque = self.forcetorque.T

        return self.forcetorque
    
    def apply_calibration_sono(self, raw_volts, cal_list):
        """Unpacks [v_low, v_high, mm_low, mm_high] and returns mm."""
        if raw_volts is None or cal_list is None:
            return None
            
        v_low, v_high, mm_low, mm_high = cal_list
        
        slope = (mm_high - mm_low) / (v_high - v_low)
        intercept = mm_low - (slope * v_low)
        return (raw_volts * slope) + intercept

    def _primary_torque_index(self):
        """Return forcetorque row index for configured primary bending torque axis."""
        axis = str(getattr(self, 'primary_bending_axis', getattr(self, 'bending_axis_sensor', 'zTorque'))).strip()
        axis_l = axis.lower()
        norm = {
            'x': 'xTorque', 'xtorque': 'xTorque',
            'y': 'yTorque', 'ytorque': 'yTorque',
            'z': 'zTorque', 'ztorque': 'zTorque',
        }.get(axis_l, axis)
        lut = {'xTorque': 3, 'yTorque': 4, 'zTorque': 5}
        if norm not in lut:
            raise ValueError(f"Unsupported primary_bending_axis={axis!r}; expected xTorque/yTorque/zTorque (or x/y/z).")
        return lut[norm]

    def _estimate_inertial_profile(self, torque_primary, anglevel_deg_s):
        """Fit torque = I*alpha + bias from a calibration run."""
        tq = np.asarray(torque_primary, dtype=float).reshape(-1)
        av = np.asarray(anglevel_deg_s, dtype=float).reshape(-1)
        if tq.size < 3 or av.size < 3:
            return None
        n = min(tq.size, av.size)
        tq = tq[:n]
        av = av[:n]
        alpha = np.gradient(av) * float(self.daq_ai_sample_rate_hz)
        m = np.isfinite(tq) & np.isfinite(alpha)
        if np.count_nonzero(m) < 3:
            return None
        X = np.column_stack([alpha[m], np.ones(np.count_nonzero(m))])
        y = tq[m]
        coef, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        return {
            'axis_sensor': str(getattr(self, 'bending_axis_sensor', 'zTorque')),
            'I_est': float(coef[0]),
            'bias_est': float(coef[1]),
        }

    def _specimen_moi_for_inertial_torque(self):
        """Return specimen-only MOI estimate (g*mm^2) from active geometry model, else 0."""
        for key in ('specimen_moi_specimen', 'specimen_moi_profile', 'specimen_moi_frustum'):
            v = getattr(self, key, None)
            if v is not None and np.isfinite(v) and float(v) > 0:
                return float(v)
        return 0.0

    def _compute_primary_torque_components(self, raw_primary_torque, anglevel_deg_s):
        """
        Compute raw/system/specimen/total/corrected primary torque vectors.

        corrected = raw_primary - (system_inertial + specimen_inertial)
        """
        raw = np.asarray(raw_primary_torque, dtype=float).reshape(-1)
        av = np.asarray(anglevel_deg_s, dtype=float).reshape(-1)
        n = min(raw.size, av.size)
        if n < 2:
            z = np.zeros_like(raw)
            return {'raw': raw, 'system': z.copy(), 'specimen': z.copy(), 'total': z.copy(), 'corrected': raw.copy()}

        raw = raw[:n]
        av = av[:n]
        alpha = np.gradient(av) * float(self.daq_ai_sample_rate_hz)

        system = np.zeros_like(raw)
        use_cal = bool(getattr(self, 'use_inertial_calibration', False))
        prof = getattr(self, 'inertial_calibration_profile', None)
        if use_cal and isinstance(prof, dict):
            system = float(prof.get('I_est', 0.0)) * alpha + float(prof.get('bias_est', 0.0))

        specimen = np.zeros_like(raw)
        use_specimen = bool(getattr(self, 'use_theoretical_inertial_correction', False))
        if use_specimen:
            i_spec = self._specimen_moi_for_inertial_torque()
            if i_spec > 0:
                specimen = i_spec * alpha

        total = system + specimen
        corrected = raw - total
        return {'raw': raw, 'system': system, 'specimen': specimen, 'total': total, 'corrected': corrected}

    def save_inertial_calibration_file(self, calibration_h5_path):
        """Save inertial calibration profile as a standalone H5 file."""
        prof = getattr(self, 'inertial_calibration_profile', None)
        if not isinstance(prof, dict):
            raise ValueError("No inertial_calibration_profile available to save.")
        with h5py.File(calibration_h5_path, 'w') as f:
            f.attrs['schema'] = 'bender_inertial_calibration_v1'
            f.attrs['axis_sensor'] = str(prof.get('axis_sensor', getattr(self, 'bending_axis_sensor', 'zTorque')))
            f.attrs['I_est'] = float(prof.get('I_est', 0.0))
            f.attrs['bias_est'] = float(prof.get('bias_est', 0.0))
            if hasattr(self, 'timestamp'):
                f.attrs['timestamp'] = str(self.timestamp)
        self.inertial_calibration_file = str(calibration_h5_path)

    def load_inertial_calibration_file(self, calibration_h5_path):
        """Load inertial calibration profile from H5 (non-destructive)."""
        with h5py.File(calibration_h5_path, 'r') as f:
            prof = {
                'axis_sensor': str(f.attrs.get('axis_sensor', getattr(self, 'bending_axis_sensor', 'zTorque'))),
                'I_est': float(f.attrs.get('I_est', 0.0)),
                'bias_est': float(f.attrs.get('bias_est', 0.0)),
            }
        self.inertial_calibration_profile = prof
        self.inertial_calibration_file = str(calibration_h5_path)
        return prof

    def _normalize_dispatch_aliases(self):
        """
        Copy legacy short attribute names onto canonical ``isometric_*`` / ``isovelocity_*`` names.

        After this, dispatcher and validation use one consistent naming scheme; legacy names
        remain on the instance for backward compatibility but are mirrored to canonical keys.
        """
        # --- isometric ---
        if getattr(self, 'isometric_initial', None) is None and getattr(self, 'initial', None) is not None:
            self.isometric_initial = self.initial
        if getattr(self, 'isometric_final', None) is None and getattr(self, 'final', None) is not None:
            self.isometric_final = self.final
        if getattr(self, 'isometric_num_steps', None) is None and getattr(self, 'num_steps', None) is not None:
            self.isometric_num_steps = self.num_steps
        if getattr(self, 'isometric_mode', None) is None and getattr(self, 'mode', None) is not None:
            self.isometric_mode = self.mode
        if getattr(self, 'isometric_randomize', None) is None and hasattr(self, 'randomize'):
            self.isometric_randomize = bool(self.randomize)
        if getattr(self, 'isometric_random_seed', None) is None and getattr(self, 'random_seed', None) is not None:
            self.isometric_random_seed = self.random_seed
        # Mirror the legacy isometric-only inter-step pause onto the canonical rest_between_steps_s
        # so older saved templates that carry isometric_inter_step_interval_s keep their pause.
        if getattr(self, 'isometric_inter_step_interval_s', None) is not None:
            self.rest_between_steps_s = float(self.isometric_inter_step_interval_s)
        # Dormant: the GUI no longer exposes a stim-routing-overrides field (blocks own routing).
        # Kept only so older saved templates carrying isometric_stim_overrides still load + merge.
        sp_iso = dict(getattr(self, 'isometric_stim_params', {}) or {})
        ov_iso = dict(getattr(self, 'isometric_stim_overrides', {}) or {})
        if ov_iso:
            sp_iso.update(ov_iso)
            self.isometric_stim_params = sp_iso

        # --- isovelocity ---
        if getattr(self, 'isovelocity_min_vel', None) is None and getattr(self, 'min_vel', None) is not None:
            self.isovelocity_min_vel = self.min_vel
        if getattr(self, 'isovelocity_max_vel', None) is None and getattr(self, 'max_vel', None) is not None:
            self.isovelocity_max_vel = self.max_vel
        if getattr(self, 'isovelocity_starting_strain', None) is None and getattr(self, 'starting_strain', None) is not None:
            self.isovelocity_starting_strain = self.starting_strain
        if getattr(self, 'isovelocity_num_steps', None) is None and getattr(self, 'num_steps', None) is not None:
            self.isovelocity_num_steps = self.num_steps
        if getattr(self, 'isovelocity_randomize', None) is None and hasattr(self, 'randomize'):
            self.isovelocity_randomize = bool(self.randomize)
        if getattr(self, 'isovelocity_random_seed', None) is None and getattr(self, 'random_seed', None) is not None:
            self.isovelocity_random_seed = self.random_seed
        # Dormant: see isometric note above (template back-compat only; no GUI field).
        sp_iv = dict(getattr(self, 'isovelocity_stim_params', {}) or {})
        ov_iv = dict(getattr(self, 'isovelocity_stim_overrides', {}) or {})
        if ov_iv:
            sp_iv.update(ov_iv)
            self.isovelocity_stim_params = sp_iv

        # Mirror legacy per-protocol randomize flags onto the canonical randomize_step_order so
        # older saved templates that carry isometric_randomize / isovelocity_randomize still shuffle.
        if getattr(self, 'isometric_randomize', None) or getattr(self, 'isovelocity_randomize', None):
            self.randomize_step_order = True

    def _sync_dclamp_test_segment_aliases(self):
        """
        Keep ``dclamp`` and ``test_segment_length_mm`` in sync when only one is set.

        GUI and configs may populate only the newer name; isometric/isovelocity need a
        finite clamp spacing in mm.
        """
        dc = getattr(self, 'dclamp', None)
        ts = getattr(self, 'test_segment_length_mm', None)
        if dc is not None and ts is None:
            self.test_segment_length_mm = float(dc)
        elif ts is not None and dc is None:
            self.dclamp = float(ts)

    def _effective_dclamp_mm(self):
        """Clamp spacing (mm) from ``dclamp`` or alias ``test_segment_length_mm``."""
        self._sync_dclamp_test_segment_aliases()
        dc = getattr(self, 'dclamp', None)
        if dc is not None:
            return float(dc)
        ts = getattr(self, 'test_segment_length_mm', None)
        if ts is not None:
            return float(ts)
        return None

    def _clamp_spacing_mm_valid(self) -> bool:
        """True if ``dclamp`` / ``test_segment_length_mm`` resolves to a finite value > 0 mm."""
        dc = self._effective_dclamp_mm()
        if dc is None:
            return False
        try:
            dcf = float(dc)
        except (TypeError, ValueError):
            return False
        return bool(np.isfinite(dcf) and dcf > 0)

    def _xsec_width_mm_valid(self) -> bool:
        """True if ``xsec_width`` is a finite value > 0 mm."""
        xw = getattr(self, 'xsec_width', None)
        if xw is None:
            return False
        try:
            xwf = float(xw)
        except (TypeError, ValueError):
            return False
        return bool(np.isfinite(xwf) and xwf > 0)

    def _normalize_primary_bending_axis(self):
        """Normalize primary axis to xTorque/yTorque/zTorque."""
        axis = str(getattr(self, 'primary_bending_axis', getattr(self, 'bending_axis_sensor', 'zTorque'))).strip()
        axis_l = axis.lower()
        norm = {
            'x': 'xTorque', 'xtorque': 'xTorque',
            'y': 'yTorque', 'ytorque': 'yTorque',
            'z': 'zTorque', 'ztorque': 'zTorque',
        }.get(axis_l, axis)
        self.primary_bending_axis = norm
        return norm

    def _apply_use_sono(self):
        """Apply `use_sono` flag to channel lists from config defaults."""
        use_sono = bool(getattr(self, 'use_sono', True))
        if self.sg_channels and self.sg_names:
            channels = list(self.sg_channels)
            names = list(self.sg_names)
            if use_sono:
                channels += list(self.sono_channels)
                names += list(self.sono_names)
        else:
            pairs = list(zip(self._cfg_input_channels, self._cfg_input_channel_names))
            if not use_sono:
                pairs = [(ch, nm) for ch, nm in pairs if not str(nm).lower().startswith('sono_')]
            channels = [p[0] for p in pairs]
            names = [p[1] for p in pairs]
        self.input_channels = channels
        self.input_channel_names = names

    def _stim_params_with_lateral(self, base_dict):
        """Merge top-level lateral/recruitment attributes into a stim_params dict for dispatch."""
        sp = dict(base_dict or {})
        if 'recruitment' not in sp:
            lat = getattr(self, 'lateral_mode', None)
            if lat is not None and str(lat).strip() != '':
                sp['recruitment'] = lat
            elif getattr(self, 'recruitment', None) is not None:
                sp['recruitment'] = self.recruitment
        for k in ('bilateral_mirror_motor', 'bilateral_sequential_left_frac'):
            if k not in sp and getattr(self, k, None) is not None:
                sp[k] = getattr(self, k)
        bs = sp.get('block_sequence', getattr(self, 'block_sequence', None))
        if bs is not None:
            sp['block_sequence'] = bs
        for k in ('left_stim_voltage', 'right_stim_voltage', 'block_reset_ramp_duration_s'):
            if k not in sp and getattr(self, k, None) is not None:
                sp[k] = getattr(self, k)
        return sp

    def run_experiment(self, test_type=None):
        self._sync_dclamp_test_segment_aliases()
        if test_type is None:
            test_type = getattr(self, 'test_type', None) or 'dynamic'
        self.test_type = test_type
        requested_test_type = str(test_type)
        motion_test_type = requested_test_type
        if requested_test_type == 'calibration':
            motion_test_type = str(getattr(self, 'calibration_base_test_type', 'dynamic'))
            if motion_test_type in ('isometric', 'isovelocity', 'calibration'):
                raise ValueError(
                    "calibration_base_test_type must be one of dynamic/frequency_sweep/frequency_step/"
                    "curvature_step/step_change."
                )

        # --- THE CLEANING CREW ---
        self.aidata = None
        self.forcetorque = None
        self.sono_left_mm = None
        self.sono_right_mm = None
        self.angledata = None
        self.master_logger.clear()
        self.trial_records = []
        self.acquisition_start = None
        self.stim_clamp_notices = []

        # Now, if the DAQ fails, you won't accidentally save old data!
        # Set input and output names/channels
        self.set_input_channels(input_channels=self.input_channels, input_channel_names=self.input_channel_names)
        self.set_stim_channels(*self.stim_channels)

        self._announce_pre_run_max_rotation(requested_test_type)

        # Starting position is zero: explicitly drive the motor to angle = 0° (and confirm via the
        # encoder) before any protocol runs. Never assume the motor is already at neutral (1.5).
        self.command_start_position_zero(device_name=self.device_name)

        if requested_test_type in ('isometric', 'isovelocity'):
            self._normalize_dispatch_aliases()

        # Unified dispatcher for dedicated step-protocols.
        if requested_test_type == 'isometric':
            initial = getattr(self, 'isometric_initial', None)
            final = getattr(self, 'isometric_final', None)
            num_steps = getattr(self, 'isometric_num_steps', None)
            if initial is None or final is None or num_steps is None:
                raise AttributeError(
                    "isometric dispatch requires isometric_initial, isometric_final, and isometric_num_steps "
                    "on the Bender instance (legacy names initial/final/num_steps are copied automatically)."
                )
            mode = getattr(self, 'isometric_mode', None) or 'strain'
            randomize = bool(getattr(self, 'randomize_step_order', False)) or bool(getattr(self, 'isometric_randomize', False))
            random_seed = getattr(self, 'isometric_random_seed', None)
            stim_params = self._stim_params_with_lateral(getattr(self, 'isometric_stim_params', {}) or {})
            out = self.isometric(
                initial,
                final,
                num_steps,
                mode=mode,
                randomize=randomize,
                random_seed=random_seed,
                stim_params=stim_params if len(stim_params) > 0 else None,
            )
            self.timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
            self.h5_protocol_metadata.update({'test_type': 'isometric'})
            return out

        if requested_test_type == 'isovelocity':
            min_vel = getattr(self, 'isovelocity_min_vel', None)
            max_vel = getattr(self, 'isovelocity_max_vel', None)
            starting_strain = getattr(self, 'isovelocity_starting_strain', None)
            num_steps = getattr(self, 'isovelocity_num_steps', None)
            if min_vel is None or max_vel is None or starting_strain is None or num_steps is None:
                raise AttributeError(
                    "isovelocity dispatch requires isovelocity_min_vel, isovelocity_max_vel, "
                    "isovelocity_starting_strain, and isovelocity_num_steps on the Bender instance "
                    "(legacy names min_vel/max_vel/starting_strain/num_steps are copied automatically)."
                )
            starting_strain_mode = getattr(self, 'isovelocity_starting_strain_mode', None) or 'strain'
            randomize = bool(getattr(self, 'randomize_step_order', False)) or bool(getattr(self, 'isovelocity_randomize', False))
            random_seed = getattr(self, 'isovelocity_random_seed', None)
            stim_params = self._stim_params_with_lateral(getattr(self, 'isovelocity_stim_params', {}) or {})
            out = self.isovelocity(
                min_vel,
                max_vel,
                starting_strain,
                num_steps,
                starting_strain_mode=starting_strain_mode,
                randomize=randomize,
                random_seed=random_seed,
                iso_duration_s=float(getattr(self, 'isovelocity_iso_duration_s', 0.2)),
                pre_hold_s=float(getattr(self, 'isovelocity_pre_hold_s', 0.3)),
                stim_params=stim_params if len(stim_params) > 0 else None,
            )
            self.timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
            self.h5_protocol_metadata.update({'test_type': 'isovelocity'})
            return out

        # self.calibration is assumed to be loaded in the notebook using bender.loadCalibration(...)

        # General parameters (duration/sample_rate passed from notebook)
        sample_rate = self.daq_ai_sample_rate_hz

        # Define angle/anglevel/tnorm variables for scope consistency
        angle, anglevel, tnorm = None, None, None
        duration = None #Initialize, but will be set below

        # THIS LEVEL IS ABOUT CREATING MOTOR ANGLES
        if motion_test_type == 'dynamic':
            self._organize_cycles_for_dynamic_run()
            duration = np.sum(self.period_by_cycle)

            angle, anglevel, tnorm, t = self.make_cycles_dynamic(
                self.period_by_cycle,
                self.freq_by_cycle,
                self.amp_by_cycle,
            )

        elif motion_test_type == 'frequency_sweep':
            if self.duration is None:
                raise AttributeError("frequency_sweep requires 'self.duration' to be set in the notebook first.")
            self._ensure_all_curves_for_run()
            self._ensure_cycle_stim_arrays_for_sweep()

            duration = self.duration + self.waitbefore + self.waitafter

            angle, anglevel, tnorm, sweep_freq, t = self.make_cycles_frequency_sweep(
                self.all_freqs,
                self.all_curves,
                self.amplitude_frequency_exponent,
                self.duration,
                self.waitbefore,
                nominal_frequency=getattr(self, 'sweep_nominal_frequency', None),
                nominal_curvature=getattr(self, 'sweep_nominal_curvature', None),
            )
            self.sweep_instantaneous_freq = sweep_freq

        elif motion_test_type == 'frequency_step':
            if self.duration is None:
                raise AttributeError("frequency_step requires 'self.duration' to be set in the notebook first.")

            duration = self.duration + self.waitbefore + self.waitafter

            angle, anglevel, tnorm, _, t = self.make_cycles_frequency_step(
                self.all_freqs,
                self.all_curves,
                self.duration,
                self.waitbefore,
                nominal_frequency=getattr(self, 'step_nominal_frequency', None),
                nominal_curvature=getattr(self, 'step_nominal_curvature', None),
            )

        elif motion_test_type == 'curvature_step':
            if self.duration is None:
                raise AttributeError("curvature_step requires 'self.duration' to be set in the notebook first.")

            duration = self.duration + self.waitbefore + self.waitafter

            angle, anglevel, tnorm, _, t = self.make_cycles_curvature_step(
                self.all_freqs,
                self.all_curves,
                self.duration,
                self.waitbefore,
                nominal_frequency=getattr(self, 'step_nominal_frequency', None),
                nominal_curvature=getattr(self, 'step_nominal_curvature', None),
            )

        elif motion_test_type == 'step_change':
            freqs = getattr(self, 'step_change_frequencies', None)
            curves = getattr(self, 'step_change_curves', None)
            cps = getattr(self, 'step_change_cycles_per_step', None)
            if freqs is None or curves is None or cps is None:
                raise AttributeError(
                    "step_change requires step_change_frequencies, step_change_curves, "
                    "and step_change_cycles_per_step on the Bender instance."
                )
            angle, anglevel, tnorm, t, movedur = self.make_cycles_step_change(
                freqs,
                curves,
                cps,
                dclamp=getattr(self, 'dclamp', None),
                amp_step_vel=getattr(self, 'step_change_amp_step_vel', None),
            )
            duration = movedur

        else:
            raise ValueError(
                f"Unknown test type: {test_type!r} (use 'dynamic', 'frequency_sweep', "
                f"'frequency_step', 'curvature_step', 'step_change', 'isometric', "
                "'isovelocity', or 'calibration')"
            )

        # Access parameters for electrical stimuli stored in 'self' (set by organize_cycles)
        stimulation_params = {
            'is_stim': self.is_stim, # Must be set in notebook
            'stim_pulse_rate': self.stim_pulse_rate, # Must be set in notebook
            'prestim_time': self.prestim_time, # Must be set in notebook
            'poststim_time': self.poststim_time, # Must be set in notebook
            'prepoststim_dur': self.prepoststim_dur, # Must be set in notebook
            'prepoststim_sep': self.prepoststim_sep, # Must be set in notebook
            'movedur': duration, 
            'stimburstdur': self.stimburstdur,        
            'duty_by_cycle': self.duty_by_cycle,     
            'freq_by_cycle': self.freq_by_cycle,    
            'phase_by_cycle': self.phase_by_cycle, 
            }

        # 1. Make electrical stimuli
        S1stimcmd, S2stimcmd = self.make_stimuli(
            **stimulation_params, 
            t_basis=t,         # Pass the new long time array
            tnorm_basis=tnorm  # Pass the new normalized time array
        )
            
        # 2. Record stimulation (CRITICAL: Ensure these update self.S1stimcmd/S2stimcmd)
        self.record_stim_signal(S1stimcmd, S2stimcmd)

        # 3. Record (and SET) the generated signals
        # This MUST update self.t to match the new duration
        self.record_motor_signal(t, angle, anglevel, tnorm)

        # --- ADD THESE THREE LINES HERE FOR SAFETY ---
        # Force the class attributes to match the newly generated data
        self.t = t
        self.angle = angle
        self.anglevel = anglevel
        self.S1stimcmd = S1stimcmd
        self.S2stimcmd = S2stimcmd

        self.master_logger.record(test_type=requested_test_type, motion_test_type=motion_test_type)
        self.h5_protocol_metadata = self.master_logger.as_dict()
        # Pulse width applies only when stim fires; log null otherwise (1.4).
        self.h5_protocol_metadata['pulse_width_ms'] = (
            float(getattr(self, 'pulse_width_ms', 2.0)) if bool(self.is_stim) else None
        )

        # Create motor stepper pulses based on the generated angle/anglevel signals (MOTION ONLY)
        self.make_motor_stepper_pulses(
                        daq_ao_do_sample_rate_hz=self.daq_ao_do_sample_rate_hz,
                        motor_gear_ratio=self.motor_gear_ratio,
                        motor_full_steps_per_rev=self.motor_full_steps_per_rev
                    )
        self._validate_pre_run_signals()
                            

        # Print file save location
        filename = self.increment_file_name(f'experiment_data_{requested_test_type}_000.h5')
        print(f"Data will be saved to: {filename}")

        # Run the experiment using 'self'
        self.aidata = self.run(device_name=self.device_name)
        if getattr(self, 'simulation_mode', False):
            self.h5_protocol_metadata['simulation_mode'] = True
            self.h5_protocol_metadata['simulation_material'] = str(getattr(self, 'simulation_material', ''))
            self.h5_protocol_metadata['simulation_model'] = 'cantilever_solid_tube_OD_25.4mm'

        # --- 1. Process Force/Torque (Always assumed to exist) ---
        # Find exactly where the 6 SG channels are, regardless of order
        sg_names = ['xForce', 'yForce', 'zForce', 'xTorque', 'yTorque', 'zTorque']
        forcetorque_indices = [self.input_channel_names.index(n) for n in sg_names if n in self.input_channel_names]

        if len(forcetorque_indices) == 6:
            self.forcetorque = self.apply_calibration_forcetorque(self.aidata[forcetorque_indices, :])
            self.forcetorque_raw = np.array(self.forcetorque, copy=True)
            self.forcetorque_corrected = None
            ns = int(self.forcetorque.shape[1]) if np.ndim(self.forcetorque) == 2 else 0
            self.inertial_torque_system_primary = np.zeros(ns, dtype=float)
            self.inertial_torque_specimen_primary = np.zeros(ns, dtype=float)
            self.inertial_torque_total_primary = np.zeros(ns, dtype=float)
            self.primary_torque_raw = np.zeros(ns, dtype=float)
            self.primary_torque_corrected = np.zeros(ns, dtype=float)
            use_cal = bool(getattr(self, 'use_inertial_calibration', False))
            cal_file = getattr(self, 'inertial_calibration_file', None)
            if use_cal and getattr(self, 'inertial_calibration_profile', None) is None and cal_file and os.path.exists(cal_file):
                try:
                    self.load_inertial_calibration_file(cal_file)
                except Exception as e:
                    print(f"⚠️ Inertial calibration file could not be loaded: {e}")
            prof = getattr(self, 'inertial_calibration_profile', None)
            if requested_test_type != 'calibration':
                idx_t = self._primary_torque_index()
                comp = self._compute_primary_torque_components(self.forcetorque[idx_t, :], self.anglevel)
                self.inertial_torque_system_primary = np.array(comp['system'], copy=True)
                self.inertial_torque_specimen_primary = np.array(comp['specimen'], copy=True)
                self.inertial_torque_total_primary = np.array(comp['total'], copy=True)
                self.primary_torque_raw = np.array(comp['raw'], copy=True)
                self.primary_torque_corrected = np.array(comp['corrected'], copy=True)
                corr = np.array(self.forcetorque, copy=True)
                corr[idx_t, :len(self.primary_torque_corrected)] = self.primary_torque_corrected
                self.forcetorque_corrected = corr
                self.h5_protocol_metadata.update({
                    'system_inertial_from_profile': bool(use_cal and isinstance(prof, dict)),
                    'specimen_inertial_from_geometry': bool(self._specimen_moi_for_inertial_torque() > 0),
                    'theoretical_i_total_system_g_mm2': float(getattr(self, 'i_total_system', 0.0)),
                })
            if requested_test_type == 'calibration':
                idx_t = self._primary_torque_index()
                prof = self._estimate_inertial_profile(self.forcetorque[idx_t, :], self.anglevel)
                if prof is not None:
                    self.inertial_calibration_profile = prof
                    self.h5_protocol_metadata.update({
                        'inertial_axis_sensor': prof['axis_sensor'],
                        'inertial_I_est': prof['I_est'],
                        'inertial_bias_est': prof['bias_est'],
                    })
        else:
            print("⚠️ Warning: Could not find all 6 SG channels for Force/Torque calibration.")

          # --- Process Sonometer (Checks if they exist in config first) ---
        self.sono_left_mm = None
        self.sono_right_mm = None

        # Look for 'sono_left' in master name list
        if 'sono_left' in self.input_channel_names:
            raw_l = self.get_data_by_name('sono_left')
            self.sono_left_mm = self.apply_calibration_sono(raw_l, self.sono_cal_left)
            print(f"📏 Sonometer (Left) Calibrated: {np.mean(self.sono_left_mm):.2f} mm")

        # Look for 'sono_right' in master name list
        if 'sono_right' in self.input_channel_names:
            raw_r = self.get_data_by_name('sono_right')
            self.sono_right_mm = self.apply_calibration_sono(raw_r, self.sono_cal_right)
            print(f"📏 Sonometer (Right) Calibrated: {np.mean(self.sono_right_mm):.2f} mm")

        # --- 3. Process Stim Monitor (ONLY if it exists) ---
        if 'stim_monitor' in self.input_channel_names:
            self.stim_monitor = self.get_data_by_name('stim_monitor')

        self.timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        try:
            self.make_cycle_tags()
            cyc = np.array(getattr(self, 'cycle_index_history', np.array([])), copy=True)
        except Exception as exc:
            if not getattr(self, 'simulation_mode', False):
                import warnings
                warnings.warn(f'make_cycle_tags failed: {exc}', UserWarning, stacklevel=2)
            cyc = np.array([], dtype=int)
        entry = self._build_trial_record(
            test_type=requested_test_type,
            trial_index=0,
            cycle_index=0,
            extra={'cycle_index_by_sample': cyc},
        )
        self.trial_records = [entry]
        self.h5_protocol_metadata.update({
            'test_type': str(requested_test_type),
            'n_trials': 1,
            'motion_test_type': str(motion_test_type),
        })

        # Return to resting length: command angle = 0° and confirm before the trial completes (1.7).
        self.return_to_resting_length(device_name=self.device_name)

   
    def get_data_by_name(self, name):
        """Returns the data row for a specific channel name."""
        try:
            idx = self.input_channel_names.index(name)
            return self.aidata[idx, :]
        except (ValueError, AttributeError):
            return None

    def convert_to_curvature(self, value, mode):
        """Instance wrapper for :func:`convert_to_curvature` using ``self.dclamp`` / ``self.xsec_width``."""
        dc = getattr(self, 'dclamp', None)
        if dc is None:
            dc = getattr(self, 'test_segment_length_mm', None)
        if dc is not None:
            dc = float(dc)
        return convert_to_curvature(
            value,
            mode,
            dclamp_mm=dc,
            xsec_width_mm=getattr(self, 'xsec_width', None),
            target_muscle_depth_mm=getattr(self, 'target_muscle_depth_mm', 0.0),
        )

    def get_all_amps(self, target_value, mode='curvature', rate=None, rate_mode=None,
                     dclamp_mm=None, xsec_width_mm=None, target_muscle_depth_mm=None):
        """
        Map target specimen inputs to curvature (1/m) and motor amplitudes (deg), matching
        ``organize_cycles`` geometry.

        Parameters
        ----------
        target_value : float or sequence
            Interpreted per ``mode`` (see :func:`convert_to_curvature`).
        mode : str
            ``curvature``, ``strain`` (fraction), ``strain_pct``, or ``angle`` (deg).
        rate : float or sequence, optional
            Must be accompanied by ``rate_mode`` (no default): this is a **kinematic** rate
            (ramp derivative), not peak cyclic rate. For cyclic peak strain rate at frequency
            ``f``, compute ``peak_strain_rate_sinusoidal(strain, f)`` and pass with
            ``rate_mode='strain_rate'``.
        rate_mode : str, required if ``rate`` is given
            One of ``curvature_rate``, ``strain_rate``, ``strain_pct_rate``, ``angle_vel``.
        dclamp_mm, xsec_width_mm, target_muscle_depth_mm
            Overrides; defaults from ``self``.

        Returns
        -------
        dict
            ``curvature_1_per_m``, ``amp_deg`` (arrays); if ``rate`` given, also
            ``curvature_rate_1_per_m_s`` and ``angle_vel_deg_s``.
        """
        dc = float(dclamp_mm) if dclamp_mm is not None else getattr(self, 'dclamp', None)
        if dc is None:
            dc = getattr(self, 'test_segment_length_mm', None)
        xw = float(xsec_width_mm) if xsec_width_mm is not None else getattr(self, 'xsec_width', None)
        md = float(target_muscle_depth_mm) if target_muscle_depth_mm is not None else float(
            getattr(self, 'target_muscle_depth_mm', 0.0) or 0.0
        )
        if dc is None:
            raise ValueError("get_all_amps requires dclamp (argument or self.dclamp from metadata).")
        tv = np.atleast_1d(np.asarray(target_value, dtype=float)).reshape(-1)
        kappa = convert_to_curvature(tv, mode, dclamp_mm=dc, xsec_width_mm=xw, target_muscle_depth_mm=md)
        amps = np.rad2deg(kappa * (dc / 1000.0))
        out = {
            'curvature_1_per_m': kappa,
            'amp_deg': amps,
        }
        if rate is not None:
            if rate_mode is None:
                raise ValueError(
                    "get_all_amps(..., rate=...) requires rate_mode= "
                    "('curvature_rate', 'strain_rate', 'strain_pct_rate', or 'angle_vel'). "
                    "Position mode does not imply a unique rate: cyclic peak strain rate depends "
                    "on frequency (use peak_strain_rate_sinusoidal(strain, f) with rate_mode='strain_rate')."
                )
            rm = rate_mode
            rv = np.atleast_1d(np.asarray(rate, dtype=float)).reshape(-1)
            if rv.size == 1 and kappa.size > 1:
                rv = np.full_like(kappa, rv.flat[0])
            kdot = convert_to_curvature(rv, rm, dclamp_mm=dc, xsec_width_mm=xw, target_muscle_depth_mm=md)
            out['curvature_rate_1_per_m_s'] = kdot
            out['angle_vel_deg_s'] = np.rad2deg(kdot * (dc / 1000.0))
        return out

    def _organize_cycles_for_dynamic_run(self):
        """
        Build per-cycle motion/stim arrays from current instance fields.

        Called at dynamic run start and by the GUI preview so Apply+Run and preview
        share the same cycle organization (no stale ``period_by_cycle``).
        """
        dc = getattr(self, 'dclamp', None)
        if dc is None:
            raise ValueError(
                'Dynamic run needs test_segment_length_mm (internally `dclamp`) on the Bender — '
                'usually set in the morphometrics section.'
            )
        xw = getattr(self, 'xsec_width', None)
        if xw is None:
            raise ValueError(
                'Dynamic run needs xsec_width (mm) on the Bender (organize_cycles uses it for strain metadata).'
            )

        af_raw = getattr(self, 'all_freqs', None)
        try:
            af_arr = np.asarray(af_raw, dtype=float).reshape(-1)
        except (TypeError, ValueError):
            af_arr = np.array([], dtype=float)
        if af_arr.size == 0 or not np.all(np.isfinite(af_arr)) or np.any(af_arr <= 0):
            raise ValueError('Set all_freqs (Hz list) with finite values > 0 Hz for dynamic run.')
        af = [float(v) for v in af_arr.tolist()]

        aa = getattr(self, 'all_amps', None)
        mode = getattr(self, 'all_amps_mode', None) or getattr(self, 'curve_input_mode', None)
        ac = getattr(self, 'all_curves', None)
        if aa is not None and mode is not None:
            conv = self.get_all_amps(aa, mode=mode)
            all_curves = np.asarray(conv['curvature_1_per_m'], dtype=float).reshape(-1)
        elif ac is not None:
            all_curves = np.asarray(ac, dtype=float).reshape(-1)
            if all_curves.size == 0:
                raise ValueError('Set all_amps (or all_curves) for dynamic run.')
            if aa is not None and mode is None:
                logging.warning(
                    'Dynamic run: using existing all_curves; set all_amps_mode so all_amps '
                    'is converted with the same interpretation as the GUI.'
                )
        elif aa is not None:
            raise ValueError(
                'Dynamic run: set all_amps_mode (strain, strain_pct, angle, etc.) when using '
                'all_amps — otherwise amplitudes may be misread and the motor will move too far.'
            )
        else:
            raise ValueError('Set all_amps (or all_curves) for dynamic run.')

        randomize = bool(getattr(self, 'randomize', False))
        rs = getattr(self, 'random_seed', None)
        if randomize and rs is not None:
            np.random.seed(int(rs))

        cps = int(getattr(self, 'cycles_per_step', 0) or 0)
        nec = int(getattr(self, 'n_end_cycles', 0) or 0)
        if cps <= 0:
            raise ValueError('cycles_per_step must be a positive integer for dynamic run.')

        stim_ix = getattr(self, 'stim_cycles_in_step', None)
        if stim_ix is None:
            stim_ix = []
        stim_ix = np.asarray(stim_ix, dtype=int).reshape(-1).tolist()

        d = getattr(self, 'all_stimduties', None)
        p = getattr(self, 'all_stimphases', None)
        if d is None or (isinstance(d, (list, tuple, np.ndarray)) and len(d) == 0):
            d = [0.3]
        if p is None or (isinstance(p, (list, tuple, np.ndarray)) and len(p) == 0):
            p = [0.5]
        duties = np.asarray(d, dtype=float).reshape(-1).tolist()
        phases = np.asarray(p, dtype=float).reshape(-1).tolist()
        spr = float(getattr(self, 'stim_pulse_rate', 0.0) or 0.0)

        md = float(getattr(self, 'target_muscle_depth_mm', 0.0) or 0.0)
        self.organize_cycles(
            list(all_curves),
            af,
            randomize,
            cps,
            nec,
            float(dc),
            float(xw),
            stim_ix,
            duties,
            phases,
            spr,
            target_muscle_depth_mm=md,
        )

    def _ensure_all_curves_for_run(self):
        """Build ``self.all_curves`` (κ, 1/m) from ``all_amps`` + mode if it is absent.

        Mirrors the curve auto-conversion done elsewhere (``get_all_amps``) so that the
        ``frequency_sweep`` branch of :meth:`run_experiment` does not read a missing
        ``all_curves`` attribute — the analogue of ``period_by_cycle`` being built by
        :meth:`_organize_cycles_for_dynamic_run` before the dynamic branch reads it.
        """
        ac = getattr(self, 'all_curves', None)
        if ac is not None:
            try:
                if np.asarray(ac, dtype=float).reshape(-1).size > 0:
                    return
            except (TypeError, ValueError):
                pass
        aa = getattr(self, 'all_amps', None)
        vals = aa if aa is not None else getattr(self, 'curve_input_values', None)
        mode = getattr(self, 'all_amps_mode', None) or getattr(self, 'curve_input_mode', None)
        if vals is not None and mode is not None:
            conv = self.get_all_amps(vals, mode=mode)
            self.all_curves = np.asarray(conv['curvature_1_per_m'], dtype=float).tolist()
            return
        raise AttributeError(
            "frequency_sweep requires all_curves (or all_amps + all_amps_mode) to be set "
            "before run_experiment."
        )

    def _ensure_cycle_stim_arrays_for_sweep(self):
        """Guarantee the per-cycle stim arrays exist before the shared ``make_stimuli`` call.

        The ``frequency_sweep`` motion path never runs :meth:`_organize_cycles_for_dynamic_run`,
        so the per-cycle arrays it populates (``freq_by_cycle``, ``phase_by_cycle``,
        ``duty_by_cycle``, ``stimburstdur``, ``period_by_cycle``) may be missing and the shared
        stimulation-params dict in :meth:`run_experiment` would raise ``AttributeError``.

        A sweep is a continuous chirp with no discrete cycles, and its ``tnorm`` is in radians
        rather than cycle counts, so the cycle-indexed stim loop in :meth:`make_stimuli` is not
        meaningful here. Set these arrays to empty (clearing any stale values left over from a
        prior dynamic run) so that loop is a no-op while the absolute-time pre/post stim bursts
        still fire.
        """
        for name in ('freq_by_cycle', 'phase_by_cycle', 'duty_by_cycle', 'stimburstdur', 'period_by_cycle'):
            setattr(self, name, np.array([], dtype=float))

    def organize_cycles(self, all_curves, all_freqs, randomize, cycles_per_step, n_end_cycles, dclamp, xsec_width, stim_cycles_in_step, all_stimduties, all_stimphases, stim_pulse_rate, target_muscle_depth_mm=0.0):
        """Build per-cycle arrays. ``all_curves`` is κ (1/m); use :meth:`get_all_amps` to build it from strain/angle."""
        start = time.time()
        self.dclamp = float(dclamp)
        if not np.isfinite(cycles_per_step) or cycles_per_step <= 0:
            raise ValueError(
                f"organize_cycles: cycles_per_step must be a positive finite number; got {cycles_per_step!r}."
            )
        cycles_per_step = int(cycles_per_step)
        if not np.isfinite(n_end_cycles) or n_end_cycles < 0:
            raise ValueError(
                f"organize_cycles: n_end_cycles must be a non-negative finite integer; got {n_end_cycles!r}."
            )
        n_end_cycles = int(n_end_cycles)
        # 1. Build combinations
        combos = []
        for c1 in all_curves:
            for f1 in all_freqs:
                for d1 in all_stimduties:
                    for p1 in all_stimphases:
                        combos.append((c1, f1, d1, p1))
        
        all_curves_arr, all_freqs_arr, all_stimduties_arr, all_stimphases_arr = map(np.array, zip(*combos))

        # --- 2. CALCULATE SECONDARY VARIABLES BEFORE RANDOMIZING ---
        # This ensures they are "locked" to the specific curve/freq combo
        lever_m = _strain_lever_arm_m(xsec_width, target_muscle_depth_mm)
        all_strains_arr = lever_m * all_curves_arr
        all_strainrates_arr = 2 * np.pi * all_strains_arr * all_freqs_arr
        all_degs_arr = np.rad2deg(all_curves_arr * (dclamp/1000))

        # --- 3. RANDOMIZE EVERYTHING TOGETHER ---
        if randomize:
            order = np.arange(len(all_freqs_arr))
            np.random.shuffle(order)
            
            # Shuffle primary variables
            all_curves_arr = all_curves_arr[order]
            all_freqs_arr = all_freqs_arr[order]
            all_stimduties_arr = all_stimduties_arr[order]
            all_stimphases_arr = all_stimphases_arr[order]
            
            # SHUFFLE SECONDARY VARIABLES TOO!
            all_strains_arr = all_strains_arr[order]
            all_strainrates_arr = all_strainrates_arr[order]
            all_degs_arr = all_degs_arr[order]

        # --- 4. CREATE BY-CYCLE ARRAYS ---
        freq_by_cycle = np.repeat(all_freqs_arr, cycles_per_step)
        amp_by_cycle  = np.repeat(all_degs_arr, cycles_per_step)
        duty_by_cycle = np.repeat(all_stimduties_arr, cycles_per_step)
        phase_by_cycle = np.repeat(all_stimphases_arr, cycles_per_step)
        
        # New "By-Cycle" mappings for R
        strain_by_cycle = np.repeat(all_strains_arr, cycles_per_step)
        strainrate_by_cycle = np.repeat(all_strainrates_arr, cycles_per_step)

        # 5. Add end cycles padding
        freq_by_cycle = np.concatenate((freq_by_cycle, [all_freqs_arr[-1]] * n_end_cycles))
        amp_by_cycle  = np.concatenate((amp_by_cycle, [all_degs_arr[-1]] * n_end_cycles))
        duty_by_cycle = np.concatenate((duty_by_cycle, [all_stimduties_arr[-1]] * n_end_cycles))
        phase_by_cycle = np.concatenate((phase_by_cycle, [all_stimphases_arr[-1]] * n_end_cycles))
        
        strain_by_cycle = np.concatenate((strain_by_cycle, [all_strains_arr[-1]] * n_end_cycles))
        strainrate_by_cycle = np.concatenate((strainrate_by_cycle, [all_strainrates_arr[-1]] * n_end_cycles))

        period_by_cycle = 1.0 / freq_by_cycle

        # Only validate stim-cycle indices when stimulation is enabled. With stim off, the
        # preview should render from passive cycles alone; stim_cycles_in_step is irrelevant
        # (is_stim_cycle / stimburstdur below handle out-of-range indices harmlessly).
        if getattr(self, 'is_stim', False) and np.any(np.array(stim_cycles_in_step) >= cycles_per_step):
            raise IndexError("stim_cycles_in_step have to be less than cycles_in_step")

        c = np.arange(0, cycles_per_step)
        is_stim_cycle = np.isin(c, stim_cycles_in_step)
        is_stim_cycle = np.tile(is_stim_cycle, len(all_freqs_arr))
        is_stim_cycle = np.concatenate((is_stim_cycle, [False] * n_end_cycles))

        # When stim is disabled / stim_pulse_rate <= 0, the quantization below would be
        # 0/0 -> NaN (RuntimeWarning). Non-stim cycles are zeroed anyway, so just skip it.
        if np.isfinite(stim_pulse_rate) and stim_pulse_rate > 0:
            stimburstdur = duty_by_cycle / freq_by_cycle
            stimburstdur = np.floor(stimburstdur * stim_pulse_rate * 2) / (stim_pulse_rate * 2)
        else:
            stimburstdur = np.zeros_like(freq_by_cycle, dtype=float)
        stimburstdur[is_stim_cycle == False] = 0

        # --- 6. STORE RESULTS (For BOTH Motor Control & H5 Metadata) ---
        self.period_by_cycle = 1.0 / freq_by_cycle
        self.freq_by_cycle   = freq_by_cycle
        self.amp_by_cycle    = amp_by_cycle

        # Use the original names your motor controller expects
        self.all_degs        = all_degs_arr
        self.all_strains     = all_strains_arr
        self.all_strainrates = all_strainrates_arr
        self.duty_by_cycle   = duty_by_cycle
        self.phase_by_cycle  = phase_by_cycle
        # Keep the timeline mapping for R
        self.strain_by_cycle     = strain_by_cycle
        self.strainrate_by_cycle = strainrate_by_cycle

        # Organized outputs for MetaData (matches the shuffled experiment order)
        self.organized_freqs       = all_freqs_arr
        self.organized_curves      = all_curves_arr
        self.organized_strains     = all_strains_arr
        self.organized_strainrates = all_strainrates_arr
        self.organized_stimduties  = all_stimduties_arr
        self.organized_stimphases = all_stimphases_arr
        self.is_stim_cycle = is_stim_cycle
        self.stimburstdur = stimburstdur    
        # Create an array [0, 0, ..., 1, 1, ..., N, N]
        # len(all_freqs_arr) is the number of randomized trials
        step_indices = np.arange(len(all_freqs_arr))
        step_by_cycle = np.repeat(step_indices, cycles_per_step)

        # Add padding for the "end cycles" (label them as a special step, e.g., -1 or max+1)
        padding_steps = np.full(n_end_cycles, -1) 
        self.step_by_cycle = np.concatenate((step_by_cycle, padding_steps))
                

    # Floor below which an AI/command-timeline sample rate is treated as corruption.
    # Real acquisition for this rig runs far above this; dynamic sine commands need
    # many samples per cycle, so anything this low is a stale/garbage value.
    DAQ_AI_RATE_FLOOR_HZ = 50.0

    @property
    def daq_ai_sample_rate_hz(self):
        return self._daq_ai_sample_rate_hz

    @daq_ai_sample_rate_hz.setter
    def daq_ai_sample_rate_hz(self, value):
        """Validate the AI/command sample rate on every assignment.

        record_motor_signal() writes this from 1/dt of the last timeline, so a single
        bad value (e.g. 0.5 Hz, which equals a typical drive frequency) would otherwise
        re-save itself every run, under-sampling dynamic sine commands into a ramp and
        making the motor walk. Reject implausibly low/invalid values and snap back to
        the hardware config rate.
        """
        cfg_hz = float(getattr(self, '_config_daq_ai_sample_rate_hz', 0.0) or 0.0)
        try:
            v = float(value)
        except (TypeError, ValueError):
            v = float('nan')
        if (not np.isfinite(v) or v < self.DAQ_AI_RATE_FLOOR_HZ) and cfg_hz >= self.DAQ_AI_RATE_FLOOR_HZ:
            if getattr(self, '_daq_ai_sample_rate_hz', None) != cfg_hz:
                logging.warning(
                    'daq_ai_sample_rate_hz=%r Hz is invalid/implausibly low; resetting to '
                    'config rate %r Hz so dynamic commands stay correctly sampled.',
                    v, cfg_hz,
                )
            self._daq_ai_sample_rate_hz = cfg_hz
        else:
            self._daq_ai_sample_rate_hz = v

    def record_motor_signal(self, t, angle, anglevel, tnorm=None):
        t_arr = np.asarray(t, dtype=float)
        if t_arr.ndim != 1:
            raise ValueError(
                f"record_motor_signal: t must be a 1-D sequence; got shape {t_arr.shape}."
            )
        if t_arr.size < 2:
            raise ValueError(
                "record_motor_signal: t must have at least 2 samples to infer dt. "
                "This usually means make_cycles_* produced an empty or single-sample timeline "
                "(check duration, waitbefore/waitafter, daq_ai_sample_rate_hz, and that all "
                "cycle frequencies are finite and > 0)."
            )
        if not np.all(np.isfinite(t_arr)):
            raise ValueError(
                "record_motor_signal: t contains NaN or Inf. "
                "Check for zero/invalid Hz in all_freqs (curvature_step / dynamic) or non-finite motion duration."
            )
        dt = float(t_arr[1] - t_arr[0])
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError(
                f"record_motor_signal: invalid spacing dt = t[1]-t[0] = {dt!r} (need finite dt > 0)."
            )
        # Must match uniform timeline spacing. Do not use round(1/dt): for dt >= 2 s, 1/dt <= 0.5 and
        # round() can return 0 (e.g. round(0.5)==0), which breaks DAQmx (SampClk_Rate 0.0).
        hz = float(1.0 / dt)
        if not np.isfinite(hz) or hz <= 0:
            raise ValueError(
                f"record_motor_signal: inferred AI sample rate 1/dt = {hz!r} Hz is invalid (dt = {dt!r} s)."
            )
        self.daq_ai_sample_rate_hz = hz

        self.t = t_arr
        self.angle = np.asarray(angle, dtype=float)
        self.anglevel = np.asarray(anglevel, dtype=float)
        if self.angle.size:
            self._last_commanded_angle_deg = float(self.angle[-1])

        if tnorm is None:
            self.tnorm = copy(t_arr)
        else:
            self.tnorm = np.asarray(tnorm, dtype=float)
    
    def record_stim_signal(self, S1stimcmd, S2stimcmd):
        self.S1stimcmd = S1stimcmd
        self.S2stimcmd = S2stimcmd

    def _validate_pre_run_signals(self):
        """Fail fast on invalid motion/stimulus buffers before any NI tasks start."""
        t = np.asarray(getattr(self, 't', np.array([])), dtype=float).reshape(-1)
        tout = np.asarray(getattr(self, 'tout', np.array([])), dtype=float).reshape(-1)
        ang = np.asarray(getattr(self, 'angle', np.array([])), dtype=float).reshape(-1)
        av = np.asarray(getattr(self, 'anglevel', np.array([])), dtype=float).reshape(-1)
        stim_hi = np.asarray(getattr(self, 'stimcmdhi', np.array([])), dtype=float)
        if t.size < 2:
            raise ValueError('Run preflight failed: motion timeline must have at least 2 samples.')
        if tout.size < 2:
            raise ValueError('Run preflight failed: AO/DO timeline must have at least 2 samples.')
        if ang.size != t.size or av.size != t.size:
            raise ValueError(
                f'Run preflight failed: angle/velocity lengths must match t ({t.size}); '
                f'got angle={ang.size}, anglevel={av.size}.'
            )
        if not np.all(np.isfinite(t)) or not np.all(np.isfinite(tout)):
            raise ValueError('Run preflight failed: timeline contains NaN/Inf.')
        if not np.all(np.isfinite(ang)) or not np.all(np.isfinite(av)):
            raise ValueError('Run preflight failed: motion command contains NaN/Inf.')
        if stim_hi.ndim != 2 or stim_hi.shape[0] != 2 or stim_hi.shape[1] != tout.size:
            raise ValueError(
                'Run preflight failed: stimulus matrix must be shape (2, len(tout)); '
                f'got {stim_hi.shape}.'
            )
        if not np.all(np.isfinite(stim_hi)):
            raise ValueError('Run preflight failed: stimulus contains NaN/Inf.')
        dig = np.asarray(getattr(self, 'dig', np.array([]))).reshape(-1)
        if dig.size != tout.size:
            raise ValueError(
                f'Run preflight failed: digital buffer length must match tout ({tout.size}); got {dig.size}.'
            )
        ch = list(getattr(self, 'input_channels', []) or [])
        nm = list(getattr(self, 'input_channel_names', []) or [])
        if len(ch) != len(nm):
            raise ValueError(
                f'Run preflight failed: input_channels ({len(ch)}) and input_channel_names ({len(nm)}) '
                'must be the same length (check sono_channel vs sono_name in the hardware config).'
            )

    def _protocol_log(self, protocol, f0=None, f1=None, c0=None, c1=None, **extra):
        """Write compact protocol fields: frequency_hz, curvature_1_per_m (float or [lo, hi])."""
        payload = {'protocol': protocol}
        if f0 is not None and f1 is not None:
            payload['frequency_hz'] = _scalar_or_pair(f0, f1)
        if c0 is not None and c1 is not None:
            payload['curvature_1_per_m'] = _scalar_or_pair(c0, c1)
        payload.update(extra)
        self.master_logger.record(**payload)

    def increment_file_name(self, filename):
        m = re.search('(d+).h5', filename)
        if m is None:
            basename, ext = os.path.splitext(filename)
            num = 1
        else:
            basename = filename[:m.start(1)]
            num = int(m.group(1))
            ext = filename[m.end(1):]

        done = False
        while not done:
            filename = '{}{:03d}{}'.format(basename, num, ext)
            done = not os.path.exists(filename)
            num += 1

        self.filename = filename
        return filename

   
    def make_motor_stepper_pulses(self, daq_ao_do_sample_rate_hz=1000,
                                motor_gear_ratio=5,
                                motor_full_steps_per_rev=6400.0):

        self.daq_ao_do_sample_rate_hz = daq_ao_do_sample_rate_hz
        self.motor_gear_ratio = motor_gear_ratio

        tout = np.arange(self.t[0], self.t[-1], 1.0/daq_ao_do_sample_rate_hz)

        poshi = interpolate.interp1d(self.t, self.angle, kind='linear', assume_sorted=True, bounds_error=False,
                                    fill_value=0.0)(tout)
        velhi = interpolate.interp1d(self.t, self.anglevel, kind='linear', assume_sorted=True, bounds_error=False,
                                    fill_value=0.0)(tout)

        poshi *= motor_gear_ratio
        velhi *= motor_gear_ratio


        if daq_ao_do_sample_rate_hz == 0 or motor_full_steps_per_rev == 0:
            raise ValueError('Problems with parameters')

        stepsize = 360.0 / motor_full_steps_per_rev
        maxspeed = stepsize * daq_ao_do_sample_rate_hz / 2

        if np.any(np.abs(self.anglevel) > maxspeed):
            raise ValueError('Motion is too fast!')

        # np.round gives symmetric step counts for equal-magnitude positive and
        # negative displacements. np.floor would yield one extra step on every
        # negative ramp because it rounds toward −∞ (e.g. floor(−111.11)=−112
        # but floor(+111.11)=+111), causing a 1-step-per-trial systematic bias.
        stepnum = np.round(poshi / stepsize)
        dstep = np.diff(stepnum)
        motorstep = np.concatenate((np.array([0], dtype='uint8'), (dstep != 0).astype('uint8')))
        motordirection = (velhi <= 0).astype('uint8')

        # High = enable on driver; last sample low so DO idles off after FINITE playback.
        motorenable = np.ones_like(motordirection, dtype='uint8')
        if motorenable.size:
            motorenable[-1] = 0

        # Ensure the columns match your wires:
        dig = np.packbits(np.column_stack((
            np.zeros((len(motorstep), 5), dtype=np.uint8), 
            motorenable,    # Goes to P0.2 (BLUE)
            motorstep,      # Goes to P0.1 (BLACK)
            motordirection  # Goes to P0.0 (WHITE)
        )))
        # np.packbits always returns a uint8, so we need to convert to a uint32
        dig = dig.astype('uint32')

        self.tout = tout
        self.dig = dig

        if self.S1stimcmd is None:
            self.stimcmdhi = np.zeros((2, len(tout)))
        else:
            S1stimcmdhi = interpolate.interp1d(self.t, self.S1stimcmd, kind='linear', assume_sorted=True, bounds_error=False,
                                        fill_value=0.0)(tout)
            S2stimcmdhi = interpolate.interp1d(self.t, self.S2stimcmd, kind='linear', assume_sorted=True, bounds_error=False,
                                        fill_value=0.0)(tout)
            self.stimcmdhi = np.vstack((S1stimcmdhi, S2stimcmdhi))

        return tout, dig, motorstep, motordirection

    def set_input_channels(self, input_channels, input_channel_names):
        self.input_channels = input_channels
        self.input_channel_names = input_channel_names

    def set_stim_channels(self, S1stim_chan, S2stim_chan):
        self.S1stim_chan = S1stim_chan
        self.S2stim_chan = S2stim_chan

    def set_motor_channel(self, motor_port):
        self.motor_port = motor_port

    def set_encoder_channel(self, encoder_chan,
                            encoder_pulses_per_rev=10000):
        self.encoder_chan = encoder_chan
        self.encoder_pulses_per_rev = encoder_pulses_per_rev

    def calculate_moi_clamp(self, H, W, D, rho, offset_x, offset_z):
        # Calculates MOI for a single rectangular prism about a global axis 
        # parallel to H (Y-axis), using the Parallel Axis Theorem (I = I_cm + M*d^2)
        mass = rho * H * W * D
        I_cm_local = (1/12) * mass * (D**2 + W**2) 
        # Distance 'd' from component CM to the GLOBAL axis (d^2 = offset_x^2 + offset_z^2)
        distance_sq = offset_x**2 + offset_z**2
        I_total = I_cm_local + mass * distance_sq
        return I_total, mass
    
    def calculate_moi_specimen(self, rho_eff, obj_depth_length, 
                                front_h_semi, back_h_semi, 
                                front_w_semi, back_w_semi, 
                                num_samples=50, axis_offset_x=0.0, axis_offset_z=0.0):
        """
        Calculates MOI for a tapered ellipsoid specimen using vectorized NumPy math.
        
        Args:
            rho_eff (float): Density (g/mm^3)
            obj_depth_length (float): Total length of specimen along Z-axis (mm)
            front_h_semi, back_h_semi (float): Half-heights at front/back (mm)
            front_w_semi, back_w_semi (float): Half-widths at front/back (mm)
            num_samples (int): Resolution of the grid
            axis_offset_x, axis_offset_z (float): Distance from rotation axis (mm)
        """
        # 1. Create the 2D grid (X and Z planes)
        # We use a large enough X range to cover the widest part of the fish
        max_w = max(front_w_semi, back_w_semi)
        x = np.linspace(-max_w, max_w, num_samples)
        z = np.linspace(0, obj_depth_length, num_samples)
        
        # Create 2D coordinate matrices
        X, Z = np.meshgrid(x, z)
        
        # 2. Calculate the taper at every point on the Z-axis
        # f is the fraction of length from front (0) to back (1)
        f = Z / obj_depth_length
        Rx = front_w_semi * (1 - f) + back_w_semi * f # Local semi-width
        Ry = front_h_semi * (1 - f) + back_h_semi * f # Local semi-height

        # 3. Calculate local height of the oval at every (X, Z) coordinate
        # Using the ellipse equation: (x/Rx)^2 + (y/Ry)^2 = 1 -> solve for y
        # We use np.clip to prevent square roots of negative numbers for points outside the oval
    # Use np.maximum to ensure we never pass a negative number to the square root
        h_sq = 1 - (X**2 / Rx**2)
        H = 2 * Ry * np.sqrt(np.maximum(h_sq, 0)) 

        # 4. Calculate Mass and MOI for every "voxel"
        dx = x[1] - x[0]
        dz = z[1] - z[0]
        
        mass_matrix = rho_eff * H * dx * dz
        # r^2 = (x - offset_x)^2 + (z - offset_z)^2
        r_sq_matrix = (X - axis_offset_x)**2 + (Z - axis_offset_z)**2
        
        # 5. Sum the results
        total_mass = np.sum(mass_matrix)
        total_moi = np.sum(mass_matrix * r_sq_matrix)
        
        return total_moi, total_mass

    def _hardware_inertia_baseline(self, clamp_offset_mm):
        """Return baseline rotating hardware inertia/mass terms in (g*mm^2, g)."""
        RHO_PLA, RHO_STEEL, RHO_NEO = 0.001116, 0.008, 0.0075
        CLAMP_DIM = np.array([50.0, 100.0, 20.0])  # W, H, D (mm)
        M_BOLT = 12.0  # g total hardware per clamp
        m_shaft = (np.pi * (9.525 / 2) ** 2 * 304.8) * RHO_STEEL
        i_shaft = 0.5 * m_shaft * (9.525 / 2) ** 2
        r_cm = float(clamp_offset_mm) + (CLAMP_DIM[2] / 2)
        m_unit = (np.prod(CLAMP_DIM) * RHO_PLA) + (np.pi * 5 ** 2 * 3 * RHO_NEO) + M_BOLT
        i_rotating_clamps = 2 * (m_unit * r_cm ** 2)
        return i_shaft, m_shaft, i_rotating_clamps, (m_unit * 2), r_cm

    def set_frustum_inertial_model(
        self,
        height_mm,
        width_mm,
        length_mm,
        density_g_per_mm3,
        *,
        tip_scale=0.0,
        clamp_offset_mm=20.0,
        num_samples=100,
    ):
        """
        Build a specimen inertia model as an elliptical frustum and fold it into total MOI.

        The frustum is modeled with linearly varying semi-axes:
        front = tip_scale * back, back = (width_mm/2, height_mm/2).
        """
        h = float(height_mm)
        w = float(width_mm)
        L = float(length_mm)
        rho = float(density_g_per_mm3)
        ts = float(tip_scale)
        if h <= 0 or w <= 0 or L <= 0 or rho <= 0:
            raise ValueError("Frustum inputs height/width/length/density must be > 0.")
        if ts < 0 or ts > 1:
            raise ValueError("tip_scale must be in [0, 1].")

        back_h = h / 2.0
        back_w = w / 2.0
        front_h = ts * back_h
        front_w = ts * back_w

        i_spec, m_spec = self.calculate_moi_specimen(
            rho_eff=rho,
            obj_depth_length=L,
            front_h_semi=front_h,
            back_h_semi=back_h,
            front_w_semi=front_w,
            back_w_semi=back_w,
            num_samples=int(num_samples),
            axis_offset_x=0.0,
            axis_offset_z=0.0,
        )
        i_shaft, m_shaft, i_clamps, m_clamps, r_cm = self._hardware_inertia_baseline(clamp_offset_mm)
        self.i_total_system = float(i_shaft + i_clamps + i_spec)
        self.total_mass = float(m_shaft + m_clamps + m_spec)
        self.specimen_moi_frustum = float(i_spec)
        self.specimen_mass_frustum = float(m_spec)
        self.specimen_inertial_model = "elliptical_frustum"
        self.frustum_inputs = {
            'height_mm': h,
            'width_mm': w,
            'length_mm': L,
            'density_g_per_mm3': rho,
            'tip_scale': ts,
            'clamp_offset_mm': float(clamp_offset_mm),
        }
        print(
            f"Frustum model set: specimen mass={m_spec:.3f} g, "
            f"specimen MOI={i_spec:.3f} g*mm^2, total MOI={self.i_total_system:.3f} g*mm^2"
        )
        return self.i_total_system

    def _parse_profile_stations(self, stations):
        """
        Parse flexible station list for specimen cross-sections.

        Expected each item to include at least:
            {'height_mm': <float>, 'width_mm': <float>}
        Optional:
            {'position': <0..1>} or {'position_frac': <0..1>} and {'label': <str>}
        If positions are omitted, they are assigned evenly from 0..1.
        """
        if not isinstance(stations, (list, tuple)) or len(stations) < 2:
            raise ValueError("specimen_profile_stations must be a list with >= 2 station entries.")
        parsed = []
        for i, s in enumerate(stations):
            if not isinstance(s, dict):
                raise ValueError(f"Station #{i} must be a dict; got {type(s).__name__}.")
            h = float(s.get('height_mm'))
            w = float(s.get('width_mm'))
            if h <= 0 or w <= 0:
                raise ValueError(f"Station #{i}: height_mm and width_mm must be > 0.")
            p = s.get('position_frac', s.get('position', None))
            label = str(s.get('label', f'station_{i}'))
            parsed.append({'height_mm': h, 'width_mm': w, 'position': p, 'label': label})
        # Fill missing positions with even spacing
        provided = [x['position'] is not None for x in parsed]
        if not any(provided):
            n = len(parsed)
            for i, x in enumerate(parsed):
                x['position'] = 0.0 if n == 1 else float(i) / float(n - 1)
        else:
            for i, x in enumerate(parsed):
                if x['position'] is None:
                    raise ValueError("Either provide position for all stations or for none.")
                x['position'] = float(x['position'])
                if x['position'] < 0 or x['position'] > 1:
                    raise ValueError(f"Station #{i}: position must be in [0,1].")
        # Sort by position
        parsed = sorted(parsed, key=lambda d: d['position'])
        if parsed[0]['position'] != 0.0 or parsed[-1]['position'] != 1.0:
            # Normalize span so first/last map to 0/1 for convenience.
            p0 = parsed[0]['position']
            p1 = parsed[-1]['position']
            if p1 <= p0:
                raise ValueError("Station positions must span increasing values.")
            for x in parsed:
                x['position'] = (x['position'] - p0) / (p1 - p0)
        return parsed

    def make_profile_stations(
        self,
        proximal_height_mm,
        proximal_width_mm,
        distal_height_mm,
        distal_width_mm,
        *,
        middle_height_mm=None,
        middle_width_mm=None,
        middle_position=0.5,
    ):
        """
        Tiny helper to build proximal/middle/distal station dictionaries.

        Middle station is optional; if omitted, returns two-station profile.
        """
        stations = [
            {
                'label': 'proximal',
                'position': 0.0,
                'height_mm': float(proximal_height_mm),
                'width_mm': float(proximal_width_mm),
            },
            {
                'label': 'distal',
                'position': 1.0,
                'height_mm': float(distal_height_mm),
                'width_mm': float(distal_width_mm),
            },
        ]
        if middle_height_mm is not None and middle_width_mm is not None:
            stations.insert(
                1,
                {
                    'label': 'middle',
                    'position': float(middle_position),
                    'height_mm': float(middle_height_mm),
                    'width_mm': float(middle_width_mm),
                },
            )
        return stations

    def set_profiled_specimen_inertial_model(
        self,
        stations,
        length_mm,
        density_g_per_mm3,
        *,
        clamp_offset_mm=20.0,
        num_samples=120,
    ):
        """
        Build specimen inertia model from an arbitrary profile of width/height stations.

        This generalizes proximal/middle/distal inputs to any number of stations.
        """
        L = float(length_mm)
        rho = float(density_g_per_mm3)
        if L <= 0 or rho <= 0:
            raise ValueError("length_mm and density_g_per_mm3 must be > 0.")
        st = self._parse_profile_stations(stations)
        pos = np.array([s['position'] for s in st], dtype=float)
        rx = np.array([s['width_mm'] / 2.0 for s in st], dtype=float)
        ry = np.array([s['height_mm'] / 2.0 for s in st], dtype=float)

        z = np.linspace(0.0, L, int(num_samples))
        f = z / L
        rx_z = np.interp(f, pos, rx)
        ry_z = np.interp(f, pos, ry)

        # 2D X-Z integration, same style as calculate_moi_specimen.
        max_w = float(np.max(rx_z))
        x = np.linspace(-max_w, max_w, int(num_samples))
        X, Z = np.meshgrid(x, z)
        RX = np.interp(Z / L, pos, rx)
        RY = np.interp(Z / L, pos, ry)
        h_sq = 1.0 - (X ** 2 / np.maximum(RX ** 2, 1e-12))
        H = 2.0 * RY * np.sqrt(np.maximum(h_sq, 0.0))
        dx = x[1] - x[0]
        dz = z[1] - z[0]
        mass_matrix = rho * H * dx * dz
        r_sq = (X ** 2) + (Z ** 2)
        m_spec = float(np.sum(mass_matrix))
        i_spec = float(np.sum(mass_matrix * r_sq))

        i_shaft, m_shaft, i_clamps, m_clamps, _ = self._hardware_inertia_baseline(clamp_offset_mm)
        self.i_total_system = float(i_shaft + i_clamps + i_spec)
        self.total_mass = float(m_shaft + m_clamps + m_spec)
        self.specimen_moi_profile = i_spec
        self.specimen_mass_profile = m_spec
        self.specimen_inertial_model = "profiled_stations"
        self.specimen_profile_stations = st
        self.specimen_profile_length_mm = L
        self.specimen_profile_density_g_per_mm3 = rho
        self.specimen_profile_num_samples = int(num_samples)
        self.specimen_profile_clamp_offset_mm = float(clamp_offset_mm)
        print(
            f"Profiled model set: specimen mass={m_spec:.3f} g, specimen MOI={i_spec:.3f} g*mm^2, "
            f"total MOI={self.i_total_system:.3f} g*mm^2"
        )
        return self.i_total_system

    @staticmethod
    def _parse_geometry_dimension_list(values, name, *, require_positive):
        """Coerce ``values`` (list/tuple/ndarray of numbers) to a list of finite floats.

        ``require_positive`` rejects NaN/Inf and any value <= 0 (cross-section
        dimensions). Position lists allow any finite value (AoR-relative, may be
        negative) but still reject NaN/Inf.
        """
        try:
            seq = list(values)
        except TypeError:
            raise ValueError(f"{name} must be a list of numbers.")
        if len(seq) == 0:
            raise ValueError(f"{name} is empty.")
        out = []
        for i, v in enumerate(seq):
            try:
                fv = float(v)
            except (TypeError, ValueError):
                raise ValueError(f"{name}[{i}] is not a number: {v!r}.")
            if not np.isfinite(fv):
                raise ValueError(f"{name}[{i}] is not finite ({fv}).")
            if require_positive and fv <= 0:
                raise ValueError(f"{name}[{i}] must be > 0 (got {fv:g} mm).")
            out.append(fv)
        return out

    def set_specimen_geometry_inertial_model(
        self,
        heights_mm,
        depths_mm,
        positions_mm,
        density_g_per_mm3,
        *,
        clamp_offset_mm=20.0,
    ):
        """User-defined specimen geometry -> volume, mass, and specimen MOI about the
        CENTER transverse axis (the rotation axis = midpoint between clamps).

        Inputs are three equal-length, variable-length per-station lists:
            heights_mm   : x, cross-section height (mm), full dimension
            depths_mm    : y, cross-section depth/width (mm), full dimension
            positions_mm : station position relative to the AoR (mm); AoR center = 0

        Cross-section convention mirrors :meth:`calculate_moi_specimen`: an ELLIPSE
        with semi-axes (height/2, depth/2). Area per station = (pi/4) * height * depth.
        Between adjacent stations height and depth vary linearly with position, so the
        area integral is evaluated analytically per segment (no integration-count knob).

        Specimen MOI about the center transverse axis, mirroring the existing convention
        (rotation axis parallel to the height direction, ``r^2 = x^2 + z^2``):

            I_spec = rho * integral[ (depth/2)**2 / 4 + s**2 ] * area(s) ds

        where the first term is each elliptical slice's own MOI about the transverse axis
        (uses the depth/width semi-axis, as in ``calculate_moi_specimen``) and ``s`` =
        position relative to the AoR is the parallel-axis distance.

        I_spec is the SPECIMEN term ONLY. It composes downstream with the EMPIRICAL
        apparatus inertia (empty run): ``M_corrected = M_raw - (I_apparatus + I_spec)*alpha``.
        This method subtracts nothing; it stores I_spec and the parsed geometry.
        I_spec is NOT rod-scale validated.
        """
        x = self._parse_geometry_dimension_list(heights_mm, 'heights_mm (x)', require_positive=True)
        y = self._parse_geometry_dimension_list(depths_mm, 'depths_mm (y)', require_positive=True)
        s = self._parse_geometry_dimension_list(positions_mm, 'positions_mm', require_positive=False)
        n = len(x)
        if not (len(y) == n and len(s) == n):
            raise ValueError(
                "Specimen geometry lists must be equal length: got "
                f"heights={len(x)}, depths={len(y)}, positions={len(s)}."
            )
        if n < 2:
            raise ValueError("Specimen geometry needs >= 2 stations to integrate a volume.")
        rho = float(density_g_per_mm3)
        if not np.isfinite(rho) or rho <= 0:
            raise ValueError("density_g_per_mm3 must be a finite value > 0.")

        # Integrate front->back: sort stations by AoR-relative position.
        order = np.argsort(np.asarray(s, dtype=float))
        xs = np.asarray(x, dtype=float)[order]
        ys = np.asarray(y, dtype=float)[order]
        ss = np.asarray(s, dtype=float)[order]

        def _integ_0_1(poly):
            # Exact definite integral over t in [0, 1] of a polynomial.
            return float(poly.integ()(1.0))

        quarter_pi = np.pi / 4.0
        volume = 0.0
        i_spec = 0.0
        for k in range(n - 1):
            ds = float(ss[k + 1] - ss[k])
            if ds == 0.0:
                continue  # coincident stations contribute no length
            # Linear interpolation in t in [0, 1] between station k and k+1.
            xp = np.polynomial.Polynomial([xs[k], xs[k + 1] - xs[k]])
            yp = np.polynomial.Polynomial([ys[k], ys[k + 1] - ys[k]])
            sp = np.polynomial.Polynomial([ss[k], ss[k + 1] - ss[k]])
            area_p = quarter_pi * xp * yp                       # ellipse area(t)
            volume += ds * _integ_0_1(area_p)
            # slice-own MOI per mass = (depth/2)**2 / 4 = y**2 / 16 ; parallel-axis = s**2
            ispec_p = ((yp * yp) * (1.0 / 16.0) + sp * sp) * area_p
            i_spec += rho * ds * _integ_0_1(ispec_p)

        mass = rho * volume

        i_shaft, m_shaft, i_clamps, m_clamps, _ = self._hardware_inertia_baseline(clamp_offset_mm)
        self.i_total_system = float(i_shaft + i_clamps + i_spec)
        self.total_mass = float(m_shaft + m_clamps + mass)
        # Specimen-only term consumed by _specimen_moi_for_inertial_torque (when enabled).
        self.specimen_moi_specimen = float(i_spec)
        self.specimen_mass_specimen = float(mass)
        self.specimen_volume_mm3 = float(volume)
        self.specimen_inertial_model = "user_geometry_center_axis"
        self.specimen_geometry_heights_mm = [float(v) for v in xs]
        self.specimen_geometry_depths_mm = [float(v) for v in ys]
        self.specimen_geometry_positions_mm = [float(v) for v in ss]
        self.specimen_geometry_density_g_per_mm3 = rho
        self.specimen_profile_clamp_offset_mm = float(clamp_offset_mm)
        print(
            f"User-geometry model set: specimen volume={volume:.3f} mm^3, mass={mass:.3f} g, "
            f"specimen MOI (center axis)={i_spec:.3f} g*mm^2 [NOT scale-validated]"
        )
        return self.i_total_system

    def _align_vector_to_t(self, arr, n: int, *, fill: float = 0.0):
        """Reshape to length ``n`` (trim or pad with edge/fill) to match ``self.t``."""
        a = np.asarray(arr, dtype=float).reshape(-1)
        if a.size >= n:
            return np.array(a[:n], dtype=float, copy=True)
        out = np.full(n, float(fill), dtype=float)
        if a.size > 0:
            out[: a.size] = a
            out[a.size :] = float(a[-1])
        return out

    def _simulate_daq_acquisition(self):
        """
        Populate ``aidata`` / encoder traces without NI-DAQmx using cantilever tube physics
        (see ``bender_simulation``). Raw voltages are chosen so existing calibration reproduces
        the simulated wrench.
        """
        if not getattr(self, 'acquisition_start', None):
            self.acquisition_start = datetime.now().replace(microsecond=0).strftime('%Y-%m-%dT%H:%M:%S')
        from bender_simulation import (
            forcetorque_six_from_bending,
            forcetorque_to_raw_voltages,
            simulated_bending_force,
            specimen_effective_length_m,
        )

        n = int(np.asarray(self.t, dtype=float).size)
        if n < 1:
            raise ValueError('Simulation requires a non-empty timeline self.t.')

        ang = self._align_vector_to_t(self.angle, n, fill=0.0)
        av = self._align_vector_to_t(self.anglevel, n, fill=0.0)
        L_mm = getattr(self, 'dclamp', None)
        if L_mm is None:
            L_mm = getattr(self, 'test_segment_length_mm', None)

        rng = np.random.default_rng(
            int(getattr(self, 'simulation_rng_seed', 0) or 0) % (2**32)
        )
        F, _, _ = simulated_bending_force(
            ang,
            av,
            length_mm=L_mm,
            material_key=getattr(self, 'simulation_material', 'polyurethane'),
            rng=rng,
        )
        L_m = specimen_effective_length_m(L_mm)
        ft = forcetorque_six_from_bending(F, L_m)
        raw6 = forcetorque_to_raw_voltages(ft, np.asarray(self.calibration, dtype=float))

        n_ch = len(self.input_channels)
        self.aidata = np.zeros((n_ch, n), dtype=np.float64)
        sg_order = ['xForce', 'yForce', 'zForce', 'xTorque', 'yTorque', 'zTorque']
        for row, name in enumerate(sg_order):
            if name in self.input_channel_names:
                chi = self.input_channel_names.index(name)
                self.aidata[chi, :] = raw6[row, :]

        for i, nm in enumerate(self.input_channel_names):
            low = str(nm).lower()
            if 'sono' in low:
                self.aidata[i, :] = 0.02 * rng.standard_normal(n)
            if low == 'stim_monitor':
                s1 = np.asarray(getattr(self, 'S1stimcmd', np.zeros(n)), dtype=float).reshape(-1)
                s2 = np.asarray(getattr(self, 'S2stimcmd', np.zeros(n)), dtype=float).reshape(-1)
                s1 = self._align_vector_to_t(s1, n, fill=0.0)
                s2 = self._align_vector_to_t(s2, n, fill=0.0)
                self.aidata[i, :] = 0.5 * (s1 + s2) * 0.02

        enc_noise = 0.02 * rng.standard_normal(n)
        self.angledata = ang + enc_noise
        self.angle_measured = np.array(self.angledata, copy=True)
        self.endTime = datetime.now()
        return self.aidata

    def run(self, device_name):
        if getattr(self, 'simulation_mode', False):
            return self._simulate_daq_acquisition()

        if Task is None:
            raise RuntimeError(
                'NI-DAQmx is not available. Enable **Simulation mode** in the Streamlit sidebar to run without hardware, '
                'or install NI-DAQmx and the nidaqmx Python package.'
            )

        ai_hz = float(self.daq_ai_sample_rate_hz)
        ao_hz = float(self.daq_ao_do_sample_rate_hz)
        if not np.isfinite(ai_hz) or ai_hz <= 0:
            raise ValueError(
                f"DAQ AI sample rate daq_ai_sample_rate_hz must be finite and > 0; got {self.daq_ai_sample_rate_hz!r}. "
                "This often follows record_motor_signal from the motion timeline: check uniform dt between t[0] and t[1]."
            )
        if not np.isfinite(ao_hz) or ao_hz <= 0:
            raise ValueError(
                f"DAQ AO/DO sample rate daq_ao_do_sample_rate_hz must be finite and > 0; got {self.daq_ao_do_sample_rate_hz!r}."
            )

        input_channels = ['/'.join((device_name, c1)) for c1 in self.input_channels]
        S1stim_chan = '/'.join((device_name, self.S1stim_chan))
        S2stim_chan = '/'.join((device_name, self.S2stim_chan))
        motor_port = '/'.join((device_name, self.motor_port))
        encoder_chan = '/'.join((device_name, self.encoder_chan))

        with Task() as analog_in, Task() as analog_out, \
                Task() as digital_out, Task() as angle_in:
            def _stop_run_tasks():
                for tsk in (analog_in, angle_in, analog_out, digital_out):
                    try:
                        tsk.stop()
                    except Exception:
                        pass

            # set up the input channels
            for c1, name1 in zip(input_channels, self.input_channel_names):
                low_nm = name1.lower()
                # Sonometrics + stim monitor AI lines are single-ended on this rig (not DIFF).
                if 'sono' in low_nm or low_nm == 'stim_monitor':
                    t_config = TerminalConfiguration.RSE
                else:
                    t_config = TerminalConfiguration.DIFF
                
                # Pass the terminal_config to the channel setup
                analog_in.ai_channels.add_ai_voltage_chan(c1, name1, terminal_config=t_config,
                                                              min_val=-10.0, max_val=10.0)     # Change from 5.0)

            # set up the input sample frequency
            # just records as many samples as are in the output
            analog_in.timing.cfg_samp_clk_timing(self.daq_ai_sample_rate_hz,
                                                sample_mode=daq.AcquisitionType.FINITE,
                                                samps_per_chan=len(self.t))

            # set up the encoder channel
            angle_in.ci_channels.add_ci_ang_encoder_chan(encoder_chan, 'encoder',
                                    pulses_per_rev=self.encoder_pulses_per_rev)
            angle_in.timing.cfg_samp_clk_timing(self.daq_ai_sample_rate_hz,
                                                source="ai/SampleClock",
                                                sample_mode=daq.AcquisitionType.FINITE,
                                                samps_per_chan=len(self.t))

            # set up the analog output channels
            analog_out.ao_channels.add_ao_voltage_chan(S1stim_chan, 'S1stim')
            analog_out.ao_channels.add_ao_voltage_chan(S2stim_chan, 'S2stim')

            # it will run much faster than the input channels, because the digital output is linked
            # to it, and it needs to run fast so that the pulses 
            # are output fast enough for smooth motion
            analog_out.timing.cfg_samp_clk_timing(self.daq_ao_do_sample_rate_hz,
                                                sample_mode=daq.AcquisitionType.FINITE,
                                                samps_per_chan=len(self.tout))    

            # set it to start when the analog input starts
            analog_out.triggers.start_trigger.cfg_dig_edge_start_trig("ai/StartTrigger",
                                                    trigger_edge=daq.Edge.RISING)

            # set up the digital output channel
            digital_out.do_channels.add_do_chan(motor_port, 'motor')
            # use the analog output clock for digital output timing
            digital_out.timing.cfg_samp_clk_timing(self.daq_ao_do_sample_rate_hz, 
                                                source = "ao/SampleClock",
                                                sample_mode=daq.AcquisitionType.FINITE,
                                                samps_per_chan=len(self.tout))
            digital_out.triggers.start_trigger.cfg_dig_edge_start_trig("ai/StartTrigger",
                                                    trigger_edge=daq.Edge.RISING)

            # set up to read the input
            reader = AnalogMultiChannelReader(analog_in.in_stream)
            self.aidata = np.zeros((len(self.input_channels), len(self.t)), dtype=np.float64)
            
            angle_reader = CounterReader(angle_in.in_stream)
            self.angledata = np.zeros((len(self.t),), dtype=np.float64)

            # write the output
            analog_writer = AnalogMultiChannelWriter(analog_out.out_stream, 
                                                    auto_start=False)
            analog_writer.write_many_sample(self.stimcmdhi)

            digital_writer = DigitalSingleChannelWriter(digital_out.out_stream,
                                                        auto_start=False)
            nwritten = digital_writer.write_many_sample_port_uint32(self.dig)
            expected_out = int(len(self.tout))
            if int(nwritten) != expected_out:
                raise RuntimeError(
                    f'DO write incomplete: wrote {nwritten} of {expected_out} samples. '
                    'Try lowering daq_ao_do_sample_rate_hz or shortening the protocol.'
                )

            # start everthing
            # make sure to start the output first, because it'll wait until the input starts
            try:
                if not getattr(self, 'acquisition_start', None):
                    self.acquisition_start = datetime.now().replace(microsecond=0).strftime('%Y-%m-%dT%H:%M:%S')
                digital_out.start()
                analog_out.start()
                angle_in.start()
                analog_in.start()

                # FINITE AI acquires exactly len(t) samples at ai_hz → nominal duration n/ai_hz.
                # Using only t[-1]+10 is fragile (t[0]≠0, clock vs timeline mismatch, long isometric holds)
                # and triggers NI -200560 ("Wait Until Done ... timeout") when margin is too small.
                n_samples = int(len(self.t))
                nominal_duration_s = n_samples / ai_hz
                t_np = np.asarray(self.t, dtype=float)
                if n_samples > 1 and np.all(np.isfinite(t_np[[0, -1]])):
                    timeline_span_s = float(t_np[-1] - t_np[0])
                    if (not np.isfinite(timeline_span_s)) or timeline_span_s <= 0:
                        timeline_span_s = nominal_duration_s
                else:
                    timeline_span_s = nominal_duration_s
                base_duration_s = max(nominal_duration_s, timeline_span_s)
                margin_s = max(30.0, 0.25 * base_duration_s)
                wait_timeout_s = float(base_duration_s + margin_s)
                wait_timeout_s = max(wait_timeout_s, 45.0)

                analog_in.wait_until_done(wait_timeout_s)
                self.endTime = datetime.now()

                reader.read_many_sample(self.aidata)
                angle_reader.read_many_sample_double(self.angledata)
                self.angle_measured = self.angledata
            except Exception:
                time.sleep(0.05)
                try:
                    from bender_daq_kill import daq_emergency_stop
                    daq_emergency_stop(device_name)
                except Exception:
                    pass
                raise
            finally:
                _stop_run_tasks()
                # Reset/free the NI device after EVERY run (success or failure).
                # On Windows/NI a second back-to-back FINITE acquisition wedges in
                # wait_until_done() unless the device is released first. Fully guarded
                # so teardown can never raise and mask the real error.
                try:
                    from bender_daq_kill import daq_emergency_stop
                    daq_emergency_stop(device_name)
                except Exception:
                    pass

        return(self.aidata)

    def _normalize_recruitment(self, name):
        """
        Canonical recruitment modes for isometric / isovelocity.

        - left_unilateral / left: only the stim channel mapped to specimen left (config S1side/S2side).
        - right_unilateral / right: only right.
        - bilateral_simultaneous: left + right channels carry the same pulse train when active.
        - bilateral_sequential: in time order, left channel then right (same motor posture unless mirror).
        """
        if name is None:
            return 'bilateral_simultaneous'
        s = str(name).strip().lower().replace('-', '_').replace(' ', '_')
        if not s:
            return 'bilateral_simultaneous'
        aliases = {
            'bilateral_simultaneous': 'bilateral_simultaneous',
            'bilateral_both': 'bilateral_simultaneous',
            'both': 'bilateral_simultaneous',
            'simultaneous': 'bilateral_simultaneous',
            'bilateral_sequential': 'bilateral_sequential',
            'bilateral_left_then_right': 'bilateral_sequential',
            'sequential': 'bilateral_sequential',
            'left_then_right': 'bilateral_sequential',
            'left_unilateral': 'left_unilateral',
            'left': 'left_unilateral',
            'unilateral_left': 'left_unilateral',
            'unilateral': 'left_unilateral',
            'right_unilateral': 'right_unilateral',
            'right': 'right_unilateral',
            'unilateral_right': 'right_unilateral',
        }
        if s in aliases:
            return aliases[s]
        if s in ('bilateral_sequential', 'bilateral_simultaneous', 'left_unilateral', 'right_unilateral'):
            return s
        raise ValueError(
            f"Unknown recruitment/lateral_mode {name!r}; use left, right, bilateral_sequential, "
            "bilateral_simultaneous."
        )

    def _recruitment_with_bilateral_mirror_motor(self, rec, mirror_bm):
        """
        When ``bilateral_mirror_motor`` is True, the motor timeline uses left-then-right postures
        (see :meth:`_timeline_mirror_two_holds`), which requires **bilateral_sequential** stim routing.
        If the user left recruitment at **bilateral_simultaneous** or a unilateral mode, upgrade to
        **bilateral_sequential** so the checkbox matches the label “both sides”.
        """
        rec = self._normalize_recruitment(rec)
        if not mirror_bm:
            return rec
        if rec == 'bilateral_simultaneous':
            return 'bilateral_sequential'
        if rec in ('left_unilateral', 'right_unilateral'):
            return 'bilateral_sequential'
        return rec

    def _normalize_block_direction(self, direction):
        """Canonical block bend direction: ``left`` or ``right`` (specimen side)."""
        if direction is None:
            raise ValueError("block direction is required; use 'left' or 'right'.")
        s = str(direction).strip().lower()
        if s in ('left', 'right'):
            return s
        raise ValueError(f"block direction must be 'left' or 'right', not {direction!r}.")

    def _normalize_stim_sides(self, stim_sides):
        """Canonical per-block stim routing: ``left``, ``right``, ``both``, or ``off``.

        ``off`` means the block runs the motion with no stimulation: ``_route_stim_sides_volts``
        returns all-zero S1/S2 for it regardless of the global stim-enable flag.
        """
        if stim_sides is None:
            return 'left'
        s = str(stim_sides).strip().lower()
        aliases = {
            'left': 'left',
            'right': 'right',
            'both': 'both',
            'bilateral': 'both',
            'bilateral_simultaneous': 'both',
            'off': 'off',
            'none': 'off',
            'no_stim': 'off',
        }
        if s in aliases:
            return aliases[s]
        raise ValueError(f"stim_sides must be 'left', 'right', 'both', or 'off', not {stim_sides!r}.")

    def _normalize_block_sequence(self, block_sequence):
        """
        Validate and normalize an ordered list of block dicts with ``direction`` and ``stim_sides``.
        Returns ``None`` when ``block_sequence`` is absent/empty (legacy single-protocol mode).
        """
        if block_sequence is None:
            return None
        if isinstance(block_sequence, (list, tuple)) and len(block_sequence) == 0:
            return None
        if not isinstance(block_sequence, (list, tuple)):
            raise ValueError("block_sequence must be a list of {direction, stim_sides} dicts.")
        out = []
        for i, raw in enumerate(block_sequence):
            if not isinstance(raw, dict):
                raise ValueError(f"block_sequence[{i}] must be a dict with direction and stim_sides.")
            out.append({
                'direction': self._normalize_block_direction(raw.get('direction')),
                'stim_sides': self._normalize_stim_sides(raw.get('stim_sides', 'left')),
            })
        return out

    def _lateral_index_for_block_direction(self, direction):
        """Specimen lateral index for a block bend direction."""
        d = self._normalize_block_direction(direction)
        if d == 'left':
            return int(self.specimen_side_index_left)
        return int(self.specimen_side_index_right)

    def _stim_sides_to_recruitment(self, stim_sides):
        """Map block ``stim_sides`` to recruitment label for metadata only."""
        s = self._normalize_stim_sides(stim_sides)
        if s == 'left':
            return 'left_unilateral'
        if s == 'right':
            return 'right_unilateral'
        return 'bilateral_simultaneous'

    def _validate_block_sequence_voltages(self, block_sequence, left_voltage, right_voltage):
        """Require finite voltages > 0 for each stim side used in any block."""
        sides_used = {b['stim_sides'] for b in block_sequence}
        if 'left' in sides_used or 'both' in sides_used:
            lv = float(left_voltage)
            if not np.isfinite(lv) or lv <= 0:
                raise ValueError(
                    f"left_stim_voltage must be finite and > 0 when LEFT stim is used; got {left_voltage!r}."
                )
        if 'right' in sides_used or 'both' in sides_used:
            rv = float(right_voltage)
            if not np.isfinite(rv) or rv <= 0:
                raise ValueError(
                    f"right_stim_voltage must be finite and > 0 when RIGHT stim is used; got {right_voltage!r}."
                )

    def _resolve_stim_onset_duration_s(self, sp, *, segment_duration_s):
        """
        Return ``(stim_onset_s, stim_duration_s)`` relative to active-segment start.

        Migrates legacy ``settle_before_stim_s`` / ``pre_iso_stim_duration_s`` when
        ``stim_onset_s`` is absent. Legacy ``stim_duration_s=None`` means through segment end.
        """
        sp = dict(sp or {})
        if sp.get('stim_onset_s') is not None:
            onset = float(sp['stim_onset_s'])
        elif float(sp.get('pre_iso_stim_duration_s', 0) or 0) > 0:
            onset = -float(sp['pre_iso_stim_duration_s'])
        else:
            onset = float(sp.get('settle_before_stim_s', 0) or 0)
        seg = float(segment_duration_s)
        if sp.get('stim_duration_s') is not None:
            dur = float(sp['stim_duration_s'])
        else:
            dur = seg - max(0.0, onset)
        return onset, dur

    def _validate_stim_timing_bounds(
        self,
        *,
        step_index,
        stim_onset_s,
        stim_duration_s,
        pre_hold_at_start_s,
        segment_duration_s,
        protocol_label='step',
    ):
        """
        Ensure stim window stays within pre-hold and active-segment bounds.

        ``stim_onset_s`` and ``stim_duration_s`` are relative to active-segment start.
        """
        onset = float(stim_onset_s)
        dur = float(stim_duration_s)
        pre_hold = float(pre_hold_at_start_s)
        seg = float(segment_duration_s)
        if not np.isfinite(onset):
            raise ValueError(f"{protocol_label} step {step_index}: stim_onset_s must be finite.")
        if not np.isfinite(dur) or dur <= 0:
            raise ValueError(f"{protocol_label} step {step_index}: stim_duration_s must be finite and > 0.")
        if onset < -pre_hold:
            overrun = float(-(onset + pre_hold))
            raise ValueError(
                f"{protocol_label} step {step_index}: stim onset {onset} s starts {overrun:.6g} s "
                f"before the allowed pre-hold window (earliest onset: {-pre_hold} s)."
            )
        if onset + dur > seg + 1e-9:
            overrun = float(onset + dur - seg)
            raise ValueError(
                f"{protocol_label} step {step_index}: stim ends {overrun:.6g} s past the active segment "
                f"({seg} s); reduce onset and/or duration."
            )

    def _clamp_stim_window_to_segment(
        self,
        onset,
        dur,
        *,
        pre_hold_at_start_s,
        segment_duration_s,
    ):
        """Clamp ``(onset, dur)`` so the stim window can never overrun the active segment.

        ``onset`` and ``dur`` are relative to active-segment start. The window is clamped into
        ``[-pre_hold_at_start_s, segment_duration_s]`` so stim never starts before the pre-hold
        window nor ends past the active-segment end, regardless of user input. Returns
        ``(onset, dur, notices)`` where ``notices`` is a list of human-readable strings describing
        any clamping that occurred (empty when the window already fit).
        """
        onset = float(onset)
        dur = float(dur)
        pre_hold = float(pre_hold_at_start_s)
        seg = float(segment_duration_s)
        notices = []
        if onset < -pre_hold:
            notices.append(
                f"Stim onset clamped from {onset:.4g} s to {-pre_hold:.4g} s "
                f"(cannot begin before the pre-hold window)."
            )
            onset = -pre_hold
        if onset > seg:
            notices.append(
                f"Stim onset clamped from {onset:.4g} s to {seg:.4g} s (active-segment end)."
            )
            onset = seg
        max_dur = seg - onset
        if dur > max_dur + 1e-9:
            notices.append(
                f"Stim duration clamped to fit active segment ({seg:.4g} s)."
            )
            dur = max(0.0, max_dur)
        return onset, dur, notices

    def _validate_stim_timing_for_steps(
        self,
        sp,
        *,
        test_type,
        num_steps,
        pre_hold_at_start_s,
        segment_duration_s,
    ):
        """Auto-clamp stim timing into the pre-hold + active-segment bounds when stim is enabled.

        Stim windows that would overrun the active segment (or begin before the pre-hold window)
        are silently clamped so the run is never blocked. Truly invalid inputs (non-finite onset,
        non-positive / non-finite duration) are still rejected loudly. When clamping occurs the
        clamped ``stim_onset_s`` / ``stim_duration_s`` are written back into ``sp`` so the generated
        stim arrays honor the clamp, and the human-readable notices are returned (and appended to
        ``self.stim_clamp_notices`` for post-run reporting).
        """
        if not bool(sp.get('is_stim', False)):
            return []
        onset, dur = self._resolve_stim_onset_duration_s(sp, segment_duration_s=segment_duration_s)
        if not np.isfinite(onset):
            raise ValueError(f"{test_type}: stim_onset_s must be finite.")
        if not np.isfinite(dur) or dur <= 0:
            raise ValueError(f"{test_type}: stim_duration_s must be finite and > 0.")
        onset_c, dur_c, notices = self._clamp_stim_window_to_segment(
            onset,
            dur,
            pre_hold_at_start_s=pre_hold_at_start_s,
            segment_duration_s=segment_duration_s,
        )
        if notices:
            sp['stim_onset_s'] = onset_c
            sp['stim_duration_s'] = dur_c
            existing = list(getattr(self, 'stim_clamp_notices', []) or [])
            for msg in notices:
                full = f"{test_type}: {msg}"
                if full not in existing:
                    existing.append(full)
            self.stim_clamp_notices = existing
        return [f"{test_type}: {m}" for m in notices]

    def lateral_index_from_side_name(self, side):
        """Map ``'left'`` / ``'right'`` to specimen lateral index (from config)."""
        s = str(side).strip().lower()
        if s == 'left':
            return int(self.specimen_side_index_left)
        if s == 'right':
            return int(self.specimen_side_index_right)
        raise ValueError(f"side must be 'left' or 'right', not {side!r}")

    def side_name_from_lateral_index(self, lateral_index):
        """Map a specimen lateral index to ``'left'`` or ``'right'``."""
        li = int(lateral_index)
        if li == int(self.specimen_side_index_left):
            return 'left'
        if li == int(self.specimen_side_index_right):
            return 'right'
        raise ValueError(
            f"lateral_index {li!r} is not specimen_side_index_left/right "
            f"({self.specimen_side_index_left}, {self.specimen_side_index_right})"
        )

    def recruitment_unilateral_lateral_index(self, recruitment_normalized):
        """
        If ``recruitment_normalized`` is unilateral, return recruited lateral index; else ``None``.
        Expects output of :meth:`Bender._normalize_recruitment`.
        """
        if recruitment_normalized == 'left_unilateral':
            return int(self.specimen_side_index_left)
        if recruitment_normalized == 'right_unilateral':
            return int(self.specimen_side_index_right)
        return None

    def motor_positive_bend_lateral_index(self):
        """
        Specimen lateral index toward which **positive** commanded motor angle bends
        (the side named in ``positive_motor_direction``, index from config anchor).
        """
        pos = str(getattr(self, 'positive_motor_direction', 'left')).lower()
        if pos == 'left':
            return int(self.specimen_side_index_left)
        if pos == 'right':
            return int(self.specimen_side_index_right)
        return int(self.specimen_side_index_left)

    def motor_command_sign_for_bend_toward_index(self, lateral_index):
        """
        Return ``+1`` or ``-1`` to multiply a **positive** bend magnitude (angle, velocity)
        so that the resulting motion bends **toward** the specimen side named by ``lateral_index``
        (must equal ``specimen_side_index_left`` or ``specimen_side_index_right``).

        Relates encoder/motor sign convention (``positive_motor_direction``) to the specimen
        lateral frame without changing hardware geometry.
        """
        li = int(lateral_index)
        if li == 0:
            return 1.0
        pos = str(getattr(self, 'positive_motor_direction', 'left')).lower()
        if pos not in ('left', 'right'):
            return 1.0
        toward_left_side = li == int(self.specimen_side_index_left)
        motor_positive_is_left = pos == 'left'
        return 1.0 if toward_left_side == motor_positive_is_left else -1.0

    def _motor_sign_for_specimen_side(self, side):
        """Backward-compatible wrapper for :meth:`motor_command_sign_for_bend_toward_index`."""
        return self.motor_command_sign_for_bend_toward_index(self.lateral_index_from_side_name(side))

    def specimen_lateral_frame_summary(self):
        """Compact dict for metadata / debugging: indices + config wiring."""
        return {
            'specimen_lateral_index_on_positive_motor_side': int(
                getattr(self, 'specimen_lateral_index_on_positive_motor_side', -1)
            ),
            'specimen_side_index_left': int(self.specimen_side_index_left),
            'specimen_side_index_right': int(self.specimen_side_index_right),
            'positive_motor_bend_toward_lateral_index': int(self.motor_positive_bend_lateral_index()),
            'positive_motor_direction': str(getattr(self, 'positive_motor_direction', '')),
            'S1side': str(getattr(self, 'S1side', '')),
            'S2side': str(getattr(self, 'S2side', '')),
        }

    def unilateral_posture_lateral_index(self, recruitment=None):
        """
        If ``recruitment`` is left/right unilateral (after normalization), return that side's
        lateral index; else ``None``. Pass ``recruitment=None`` to use ``self.recruitment``.
        """
        r = self._normalize_recruitment(
            recruitment if recruitment is not None else getattr(self, 'recruitment', None)
        )
        return self.recruitment_unilateral_lateral_index(r)

    def strain_display_sign(self):
        """±1 scale for geometric strain in plots; see ``strain_shortening_positive_display_sign``."""
        s = float(getattr(self, 'strain_shortening_positive_display_sign', 1.0))
        return -1.0 if s < 0 else 1.0

    def strain_yaxis_title_pct(self, recruitment=None):
        """
        Y-axis label for ε = κ·(w/2) previews, tagged with recruitment / posture side when known.
        """
        rec = self._normalize_recruitment(
            recruitment if recruitment is not None else getattr(self, 'recruitment', None)
        )
        base = 'ε_geom = κ·w/2 (%)'
        tag = {
            'left_unilateral': 'posture toward LEFT',
            'right_unilateral': 'posture toward RIGHT',
            'bilateral_sequential': 'bilateral sequential (L/R stim phases)',
            'bilateral_simultaneous': 'bilateral simultaneous',
        }.get(rec, '')
        if tag:
            return f'{base} — {tag}'
        return base

    def strain_geometry_plot_context(self, recruitment=None):
        """
        Short text for figure titles/annotations: recruitment, motor sign, and fiber strain caveat.
        """
        rec = self._normalize_recruitment(
            recruitment if recruitment is not None else getattr(self, 'recruitment', None)
        )
        pos = str(getattr(self, 'positive_motor_direction', 'left')).lower()
        side_label = {
            'left_unilateral': 'LEFT unilateral (command bends toward LEFT)',
            'right_unilateral': 'RIGHT unilateral (command bends toward RIGHT)',
            'bilateral_sequential': 'bilateral sequential (split stim L then R; mirror motor optional)',
            'bilateral_simultaneous': 'bilateral simultaneous (both stim channels)',
        }.get(rec, str(rec))
        flip = self.strain_display_sign()
        flip_note = (
            ' Display sign flipped (strain_shortening_positive_display_sign<0).'
            if flip < 0
            else ''
        )
        mp_idx = int(self.motor_positive_bend_lateral_index())
        uidx = self.unilateral_posture_lateral_index(rec)
        u_note = f' Unilateral posture lateral index={uidx}.' if uidx is not None else ''
        return (
            f'{side_label}. Lateral axis: LEFT={self.specimen_side_index_left}, '
            f'RIGHT={self.specimen_side_index_right}; '
            f'+motor bends toward index {mp_idx} ({pos.upper()}).{u_note} '
            f'ε_geom = κ·(w/2−muscle_depth) from θ/L; opposite surfaces ± same |ε| (shorten vs lengthen).{flip_note}'
        )

    def _pulse_high_fraction(self, stim_pulse_rate_hz, pulse_width_ms=None):
        """Fraction (0..1) of each carrier period the stim pulse is high, from ``pulse_width_ms``.

        Falls back to 50% duty when no pulse width is configured (legacy behaviour). The fraction
        is clamped to (0, 1]: the high time can never exceed one carrier period.
        """
        pw_ms = pulse_width_ms if pulse_width_ms is not None else getattr(self, 'pulse_width_ms', None)
        pr = float(stim_pulse_rate_hz)
        if pw_ms is None or not np.isfinite(pr) or pr <= 0:
            return 0.5
        try:
            pw_s = float(pw_ms) / 1000.0
        except (TypeError, ValueError):
            return 0.5
        if not np.isfinite(pw_s) or pw_s <= 0:
            return 0.5
        return float(min(max(pw_s * pr, 0.0), 1.0))

    def _pulse_carrier_volts(self, t, active_mask, stim_pulse_rate_hz, stim_voltage):
        """Square carrier at stim_pulse_rate_hz with pulse_width_ms high time, gated by active_mask."""
        pr = float(stim_pulse_rate_hz)
        v = float(stim_voltage)
        frac = self._pulse_high_fraction(pr)
        pulse = (np.mod(t * pr, 1.0) < frac).astype(np.float64) * v
        m = active_mask & np.isfinite(t)
        return np.where(m, pulse, 0.0)

    def _deposit_stim_on_side(self, pulse_vec, side_name, s1, s2):
        """Add pulse_vec to S1 and/or S2 depending on config S1side / S2side."""
        sn = str(side_name).lower()
        if str(self.S1side).lower() == sn:
            s1[:] += np.asarray(pulse_vec, dtype=float).reshape(-1)
        if str(self.S2side).lower() == sn:
            s2[:] += np.asarray(pulse_vec, dtype=float).reshape(-1)

    def _route_stim_sides_volts(
        self,
        t,
        active_mask,
        stim_pulse_rate_hz,
        stim_sides,
        left_voltage,
        right_voltage,
    ):
        """
        Route gated pulse carriers onto S1/S2 using per-block ``stim_sides`` and per-side voltages.
        """
        n = int(np.asarray(t).size)
        s1 = np.zeros(n, dtype=float)
        s2 = np.zeros(n, dtype=float)
        sides = self._normalize_stim_sides(stim_sides)
        # OFF: this block runs the motion with no stimulation on either channel,
        # regardless of the global stim-enable flag.
        if sides == 'off':
            return s1, s2
        active = np.asarray(active_mask, dtype=bool).reshape(-1)
        if not np.any(active):
            return s1, s2
        if sides in ('left', 'both'):
            p_l = self._pulse_carrier_volts(t, active, stim_pulse_rate_hz, left_voltage)
            self._deposit_stim_on_side(p_l, 'left', s1, s2)
        if sides in ('right', 'both'):
            p_r = self._pulse_carrier_volts(t, active, stim_pulse_rate_hz, right_voltage)
            self._deposit_stim_on_side(p_r, 'right', s1, s2)
        return s1, s2

    def _route_recruitment_stim(self, pulse, recruitment, sequential_left_frac=0.5):
        """
        Route a single per-sample pulse train onto S1/S2 using recruitment mode.
        ``pulse`` is length n; returns (s1, s2).
        """
        n = int(np.asarray(pulse).size)
        s1 = np.zeros(n, dtype=float)
        s2 = np.zeros(n, dtype=float)
        rec = self._normalize_recruitment(recruitment)
        p = np.asarray(pulse, dtype=float).reshape(-1)
        active = np.abs(p) > 1e-18
        if not np.any(active):
            return s1, s2
        if rec == 'left_unilateral':
            self._deposit_stim_on_side(p, 'left', s1, s2)
            return s1, s2
        if rec == 'right_unilateral':
            self._deposit_stim_on_side(p, 'right', s1, s2)
            return s1, s2
        if rec == 'bilateral_simultaneous':
            self._deposit_stim_on_side(p, 'left', s1, s2)
            self._deposit_stim_on_side(p, 'right', s1, s2)
            return s1, s2
        if rec == 'bilateral_sequential':
            idx = np.where(active)[0]
            frac = float(sequential_left_frac)
            frac = min(0.999, max(0.001, frac))
            split = int(round(frac * idx.size))
            split = max(1, min(idx.size - 1, split))
            left_m = np.zeros(n, dtype=bool)
            right_m = np.zeros(n, dtype=bool)
            left_m[idx[:split]] = True
            right_m[idx[split:]] = True
            pl = np.where(left_m, p, 0.0)
            pr = np.where(right_m, p, 0.0)
            self._deposit_stim_on_side(pl, 'left', s1, s2)
            self._deposit_stim_on_side(pr, 'right', s1, s2)
            return s1, s2
        return s1, s2

    def _isometric_pulse_stim(self, t, active_mask, stim_pulse_rate_hz, stim_voltage):
        """Same carrier shape as make_stimuli (square at stim_pulse_rate_hz, pulse_width_ms high)."""
        pr = float(stim_pulse_rate_hz)
        v = float(stim_voltage)
        frac = self._pulse_high_fraction(pr)
        pulse = (np.mod(t * pr, 1.0) < frac).astype(np.float64) * v
        m = active_mask & np.isfinite(t)
        s = np.where(m, pulse, 0.0)
        return s, s.copy()

    def _timeline_mirror_two_holds(
        self, prev_deg, mag_deg, ramp_s, hold_s_total, daq_hz, *, mirror_abs_deg_left=None, mirror_abs_deg_right=None
    ):
        """
        Ramp prev->T1, hold T1, ramp T1->T2, hold T2 with signed T1/T2 toward left/right specimen sides.
        Default: |T1|=|T2|=|mag_deg|. Optional ``mirror_abs_deg_*`` set per-side hold magnitudes (unsigned deg).
        Each hold duration is hold_s_total/2. Used for bilateral_sequential + bilateral_mirror_motor.
        """
        h = float(hold_s_total) * 0.5
        mL = abs(float(mirror_abs_deg_left)) if mirror_abs_deg_left is not None else abs(float(mag_deg))
        mR = abs(float(mirror_abs_deg_right)) if mirror_abs_deg_right is not None else abs(float(mag_deg))
        T1 = self.motor_command_sign_for_bend_toward_index(self.specimen_side_index_left) * mL
        T2 = self.motor_command_sign_for_bend_toward_index(self.specimen_side_index_right) * mR
        t1, a1, w1 = self._timeline_ramp_hold(float(prev_deg), T1, float(ramp_s), h, daq_hz)
        t2, a2, w2 = self._timeline_ramp_hold(T1, T2, float(ramp_s), h, daq_hz)
        off = float(t1[-1])
        t2s = t2[1:] + off
        a2s = a2[1:]
        w2s = w2[1:]
        t = np.concatenate([t1, t2s])
        angle = np.concatenate([a1, a2s])
        anglevel = np.concatenate([w1, w2s])
        t_hold1_0 = float(ramp_s)
        t_hold1_1 = float(ramp_s) + h
        t_hold2_0 = off + float(ramp_s)
        t_hold2_1 = off + float(ramp_s) + h
        return t, angle, anglevel, (t_hold1_0, t_hold1_1), (t_hold2_0, t_hold2_1)

    def _mirror_hold_deg_at_step(self, arr, step_index, num_steps):
        """Scalar or length ``num_steps`` (or 1) array of unsigned motor magnitudes (deg) for mirror holds."""
        if arr is None:
            return None
        a = np.atleast_1d(np.asarray(arr, dtype=float).reshape(-1))
        if a.size == 1:
            return float(a[0])
        if int(a.size) != int(num_steps):
            raise ValueError(
                f"mirror hold degree targets: expected length 1 or {num_steps}, got {int(a.size)}."
            )
        return float(a[int(step_index)])

    def _timeline_ramp_hold(self, angle_start_deg, angle_end_deg, ramp_s, hold_s, daq_hz):
        """Piecewise-linear ramp then hold; returns (t, angle, anglevel) at AI rate."""
        dt = 1.0 / float(daq_hz)
        ramp_s = float(ramp_s)
        hold_s = float(hold_s)
        a0, a1 = float(angle_start_deg), float(angle_end_deg)
        if ramp_s <= 0.0 or abs(a1 - a0) < 1e-12:
            t_r = np.array([0.0])
            ang_r = np.array([a1])
        else:
            n_r = max(2, int(round(ramp_s / dt)) + 1)
            t_r = np.linspace(0.0, ramp_s, n_r)
            u = np.linspace(0.0, 1.0, n_r)
            ang_r = a0 + (a1 - a0) * u
        n_h = max(2, int(round(hold_s / dt)) + 1)
        if hold_s <= 0.0:
            t_h = np.array([])
            ang_h = np.array([])
        else:
            t_h = ramp_s + np.linspace(0.0, hold_s, n_h)[1:]
            ang_h = np.full(t_h.size, a1, dtype=float)
        t = np.concatenate([t_r, t_h])
        angle = np.concatenate([ang_r, ang_h])
        if t.size < 2:
            raise ValueError("_timeline_ramp_hold: timeline too short; increase ramp_s, hold_s, or daq rate.")
        anglevel = np.gradient(angle, t, edge_order=1)
        return t, angle, anglevel

    def _run_neutral_reset_segment(self, from_deg, ramp_duration_s, device_name):
        """
        Ramp commanded motor angle to 0° (straight/center), acquire one DAQ segment, no stimulation.
        Returns ``0.0`` (commanded neutral angle).
        """
        from_deg = float(from_deg)
        ramp_s = float(ramp_duration_s)
        if not np.isfinite(ramp_s) or ramp_s < 0:
            raise ValueError(f"block_reset_ramp_duration_s must be finite and >= 0; got {ramp_duration_s!r}.")
        dev = device_name if device_name is not None else getattr(self, 'device_name', None)
        if dev is None:
            raise ValueError("_run_neutral_reset_segment requires device_name or self.device_name.")
        # Already at neutral: no reset needed regardless of ramp (ramping 0°->0° is a no-op that
        # would collapse _timeline_ramp_hold to a single sample). Skip it.
        if abs(from_deg) < 1e-12:
            return 0.0
        # A real reset (motor not at 0°) cannot happen in zero time: require a positive ramp so the
        # motor is slewed back to neutral under control. Refuse to silently skip a needed reset.
        if ramp_s <= 0:
            raise ValueError(
                f"block_reset_ramp_duration_s must be > 0 s to return the motor to neutral from "
                f"{from_deg:.6g}° before the next block; got {ramp_duration_s!r}."
            )
        daq_hz = float(self.daq_ai_sample_rate_hz)
        # Floor the reset ramp so a tiny-but-positive duration still spans at least two
        # AI samples; otherwise _timeline_ramp_hold can collapse and trip "timeline too short".
        min_ramp_s = 2.0 / daq_hz
        if ramp_s < min_ramp_s:
            ramp_s = min_ramp_s
        t, angle, anglevel = self._timeline_ramp_hold(from_deg, 0.0, ramp_s, 0.0, daq_hz)
        s1 = np.zeros_like(t)
        s2 = np.zeros_like(t)
        self.record_motor_signal(t, angle, anglevel, tnorm=np.zeros_like(t))
        self.record_stim_signal(s1, s2)
        self.make_motor_stepper_pulses(
            daq_ao_do_sample_rate_hz=self.daq_ao_do_sample_rate_hz,
            motor_gear_ratio=self.motor_gear_ratio,
            motor_full_steps_per_rev=self.motor_full_steps_per_rev,
        )
        self.aidata = self.run(device_name=dev)
        return 0.0

    def _drive_to_zero_and_confirm(self, device_name=None, ramp_duration_s=None):
        """Ramp the commanded motor angle to 0° from the last commanded angle and confirm the move.

        Shared by :meth:`command_start_position_zero` (1.5) and :meth:`return_to_resting_length`
        (1.7). Ramps from ``_last_commanded_angle_deg`` to 0° via the existing move-to-position
        routine; the encoder read during that acquisition confirms position before returning.
        """
        dev = device_name if device_name is not None else getattr(self, 'device_name', None)
        if dev is None:
            # No DAQ device configured yet (e.g. early init): nothing to command safely.
            return 0.0
        ramp = float(
            ramp_duration_s if ramp_duration_s is not None
            else getattr(self, 'block_reset_ramp_duration_s', 2.0) or 2.0
        )
        from_deg = float(getattr(self, '_last_commanded_angle_deg', 0.0) or 0.0)
        self._run_neutral_reset_segment(from_deg, ramp, dev)
        self._last_commanded_angle_deg = 0.0
        return 0.0

    def command_start_position_zero(self, device_name=None, ramp_duration_s=None):
        """Explicitly drive the commanded motor angle to 0° (resting/start) and wait for the move.

        Call at apparatus initialization and at the start of every protocol. The motor's position
        at software start is never assumed to be 0: we ramp from the last commanded angle tracked on
        the instance to 0° using the existing move-to-position routine, and the encoder read during
        that acquisition confirms the position before the protocol proceeds (1.5).
        """
        return self._drive_to_zero_and_confirm(device_name, ramp_duration_s)

    def return_to_resting_length(self, device_name=None, ramp_duration_s=None):
        """Drive the motor back to angle = 0° (resting length) at the end of a protocol and confirm.

        Called after the final step (and any post-trial rest) of every protocol so the preparation
        is returned to neutral before the trial is marked complete (1.7).
        """
        return self._drive_to_zero_and_confirm(device_name, ramp_duration_s)

    def _tag_block_trial_metadata(self, entry, *, block_index, block_direction, block_stim_sides,
                                  left_stim_voltage, right_stim_voltage):
        entry['block_index'] = int(block_index)
        entry['block_direction'] = str(block_direction)
        entry['block_stim_sides'] = str(block_stim_sides)
        entry['left_stim_voltage'] = float(left_stim_voltage)
        entry['right_stim_voltage'] = float(right_stim_voltage)

    def _run_force_length_steps(
        self,
        targets_deg,
        *,
        ramp_duration_s=2.0,
        hold_duration_s=5.0,
        pre_baseline_s=1.0,
        post_baseline_s=1.0,
        stim_onset_s=None,
        settle_before_stim_s=0.5,
        stim_duration_s=None,
        inter_step_interval_s=0.0,
        reset_between_steps=False,
        reset_ramp_duration_s=None,
        is_stim=False,
        stim_pulse_rate=None,
        stim_voltage=5.0,
        device_name=None,
        recruitment=None,
        bilateral_mirror_motor=False,
        bilateral_sequential_left_frac=0.5,
        mirror_hold_deg_left=None,
        mirror_hold_deg_right=None,
        ramp_from_deg=None,
        stim_sides=None,
        left_stim_voltage=None,
        right_stim_voltage=None,
        block_direction=None,
        block_index=None,
    ):
        """Execute ramp–hold–stim–acquire for each target angle in ``targets_deg`` (degrees)."""
        targets_deg = np.atleast_1d(np.asarray(targets_deg, dtype=float)).reshape(-1)
        num_steps = int(targets_deg.size)
        if num_steps < 1:
            raise ValueError("targets_deg must contain at least one sample.")
        if stim_sides is not None and bool(bilateral_mirror_motor):
            raise ValueError("bilateral_mirror_motor cannot be used with per-block stim_sides.")
        dev = device_name if device_name is not None else getattr(self, 'device_name', None)
        if dev is None:
            raise ValueError("_run_force_length_steps requires device_name or self.device_name.")
        spr = float(stim_pulse_rate) if stim_pulse_rate is not None else float(
            getattr(self, 'stim_pulse_rate', 75.0)
        )
        daq_hz = float(self.daq_ai_sample_rate_hz)
        results = []
        use_block_stim = stim_sides is not None
        if use_block_stim:
            lv = float(left_stim_voltage if left_stim_voltage is not None else stim_voltage)
            rv = float(right_stim_voltage if right_stim_voltage is not None else stim_voltage)
            sides_norm = self._normalize_stim_sides(stim_sides)
            rec = self._stim_sides_to_recruitment(sides_norm)
        else:
            lv = float(stim_voltage)
            rv = float(stim_voltage)
            sides_norm = None
            rec = self._normalize_recruitment(
                recruitment if recruitment is not None else getattr(self, 'recruitment', 'bilateral_simultaneous')
            )
        bm = bool(bilateral_mirror_motor) and not use_block_stim
        rec = self._recruitment_with_bilateral_mirror_motor(rec, bm)
        mirror = bm and rec == 'bilateral_sequential'
        seq_frac = float(
            bilateral_sequential_left_frac
            if bilateral_sequential_left_frac is not None
            else getattr(self, 'bilateral_sequential_left_frac', 0.5)
        )
        self.set_input_channels(input_channels=self.input_channels, input_channel_names=self.input_channel_names)
        self.set_stim_channels(*self.stim_channels)

        gap_s = float(inter_step_interval_s)
        if not np.isfinite(gap_s) or gap_s < 0:
            raise ValueError(f"inter_step_interval_s must be finite and >= 0; got {inter_step_interval_s!r}.")
        reset_steps = bool(reset_between_steps)
        reset_ramp = float(
            reset_ramp_duration_s if reset_ramp_duration_s is not None
            else getattr(self, 'block_reset_ramp_duration_s', 2.0)
        )

        for i in range(num_steps):
            if i > 0 and gap_s > 0:
                time.sleep(gap_s)
            if i > 0 and reset_steps:
                # Reset to resting length: drive the motor back to 0° from the previous hold
                # angle and wait for the move to complete before the next step.
                self._run_neutral_reset_segment(float(targets_deg[i - 1]), reset_ramp, dev)
            target = float(targets_deg[i])
            if i > 0:
                prev = 0.0 if reset_steps else float(targets_deg[i - 1])
            elif ramp_from_deg is not None:
                prev = float(ramp_from_deg)
            else:
                prev = float(targets_deg[0])
            if mirror:
                kw = {}
                ml = self._mirror_hold_deg_at_step(mirror_hold_deg_left, i, num_steps)
                mr = self._mirror_hold_deg_at_step(mirror_hold_deg_right, i, num_steps)
                if ml is not None:
                    kw['mirror_abs_deg_left'] = ml
                if mr is not None:
                    kw['mirror_abs_deg_right'] = mr
                t, angle, anglevel, h1, h2 = self._timeline_mirror_two_holds(
                    prev, target, float(ramp_duration_s), float(hold_duration_s), daq_hz, **kw
                )
                active_l = (t >= h1[0] + float(settle_before_stim_s)) & (t < h1[1])
                active_r = (t >= h2[0] + float(settle_before_stim_s)) & (t < h2[1])
                if stim_duration_s is not None:
                    active_l &= t < (h1[0] + float(settle_before_stim_s) + float(stim_duration_s))
                    active_r &= t < (h2[0] + float(settle_before_stim_s) + float(stim_duration_s))
                if is_stim and (np.any(active_l) or np.any(active_r)):
                    p_l = self._pulse_carrier_volts(t, active_l, spr, stim_voltage)
                    p_r = self._pulse_carrier_volts(t, active_r, spr, stim_voltage)
                    s1 = np.zeros_like(t)
                    s2 = np.zeros_like(t)
                    self._deposit_stim_on_side(p_l, 'left', s1, s2)
                    self._deposit_stim_on_side(p_r, 'right', s1, s2)
                else:
                    s1 = np.zeros_like(t)
                    s2 = np.zeros_like(t)
                t_stim0 = float(h1[0] + float(settle_before_stim_s))
                t_stim1 = float(h2[1])
                # Mirror two-hold mode does not add silent baselines; report the active stim span
                # and zero-length pre/post windows so the metadata keys are always present.
                t_active_start = t_stim0
                t_active_end = t_stim1
                t_pre_baseline_start = t_pre_baseline_end = float(ramp_duration_s)
                t_post_baseline_start = t_post_baseline_end = float(t[-1])
            else:
                pre_b = max(0.0, float(pre_baseline_s))
                post_b = max(0.0, float(post_baseline_s))
                # Bracket the active hold with silent baseline holds at the commanded angle: the
                # motor holds for ``pre_b`` (stim off), then the active window (``hold_duration_s``),
                # then ``post_b`` (stim off). All three are a single continuous hold at ``target``.
                total_hold = pre_b + float(hold_duration_s) + post_b
                t, angle, anglevel = self._timeline_ramp_hold(
                    prev, target, float(ramp_duration_s), total_hold, daq_hz
                )
                onset, dur = self._resolve_stim_onset_duration_s(
                    {
                        'stim_onset_s': stim_onset_s,
                        'stim_duration_s': stim_duration_s,
                        'settle_before_stim_s': settle_before_stim_s,
                    },
                    segment_duration_s=float(hold_duration_s),
                )
                t_active_start = float(ramp_duration_s) + pre_b
                t_active_end = t_active_start + float(hold_duration_s)
                t_stim0 = t_active_start + onset
                t_stim1 = t_stim0 + dur
                t_stim1 = min(t_stim1, float(t[-1]) + 1e-9)
                active = (t >= t_stim0) & (t < t_stim1)
                t_pre_baseline_start = float(ramp_duration_s)
                t_pre_baseline_end = t_active_start
                t_post_baseline_start = t_active_end
                t_post_baseline_end = t_active_end + post_b
                if is_stim and np.any(active):
                    if use_block_stim:
                        s1, s2 = self._route_stim_sides_volts(
                            t, active, spr, sides_norm, lv, rv
                        )
                    else:
                        pulse = self._pulse_carrier_volts(t, active, spr, stim_voltage)
                        s1, s2 = self._route_recruitment_stim(pulse, rec, sequential_left_frac=seq_frac)
                else:
                    s1 = np.zeros_like(t)
                    s2 = np.zeros_like(t)

            self.record_motor_signal(t, angle, anglevel, tnorm=np.zeros_like(t))
            self.record_stim_signal(s1, s2)
            self.make_motor_stepper_pulses(
                daq_ao_do_sample_rate_hz=self.daq_ao_do_sample_rate_hz,
                motor_gear_ratio=self.motor_gear_ratio,
                motor_full_steps_per_rev=self.motor_full_steps_per_rev,
            )
            self.aidata = self.run(device_name=dev)

            entry = {
                'test_type': 'isometric',
                'step_index': i,
                'trial_index': i,
                'cycle_index': i,
                'recruitment': rec,
                'unilateral_posture_lateral_index': self.recruitment_unilateral_lateral_index(rec),
                'motor_positive_bend_toward_lateral_index': int(self.motor_positive_bend_lateral_index()),
                'bilateral_mirror_motor': bool(mirror),
                'target_deg': target,
                'ramp_from_deg': prev,
                't': np.array(t, copy=True),
                'angle_cmd': np.array(angle, copy=True),
                'anglevel_cmd': np.array(anglevel, copy=True),
                'tnorm': np.zeros_like(t),
                'S1stimcmd': np.array(s1, copy=True),
                'S2stimcmd': np.array(s2, copy=True),
                'aidata': np.array(self.aidata, copy=True),
                'angle_measured': np.array(self.angle_measured, copy=True),
                'stim_t0': t_stim0,
                'stim_t1': t_stim1,
                't_pre_baseline_start': float(t_pre_baseline_start),
                't_pre_baseline_end': float(t_pre_baseline_end),
                't_active_start': float(t_active_start),
                't_active_end': float(t_active_end),
                't_post_baseline_start': float(t_post_baseline_start),
                't_post_baseline_end': float(t_post_baseline_end),
                'forcetorque': None,
                'mean_xforce_stim': None,
            }
            if use_block_stim:
                self._tag_block_trial_metadata(
                    entry,
                    block_index=block_index if block_index is not None else 0,
                    block_direction=block_direction or sides_norm,
                    block_stim_sides=sides_norm,
                    left_stim_voltage=lv,
                    right_stim_voltage=rv,
                )
            sg_names = ['xForce', 'yForce', 'zForce', 'xTorque', 'yTorque', 'zTorque']
            try:
                idx = [self.input_channel_names.index(n) for n in sg_names if n in self.input_channel_names]
                if len(idx) == 6:
                    ft = self.apply_calibration_forcetorque(self.aidata[idx, :])
                    entry['forcetorque'] = ft
                    entry['forcetorque_raw'] = np.array(ft, copy=True)
                    idx_t = self._primary_torque_index()
                    comp = self._compute_primary_torque_components(ft[idx_t, :], anglevel)
                    ft_corr = np.array(ft, copy=True)
                    ft_corr[idx_t, :len(comp['corrected'])] = comp['corrected']
                    entry['forcetorque_corrected'] = ft_corr
                    entry['inertial_torque_system_primary'] = np.array(comp['system'], copy=True)
                    entry['inertial_torque_specimen_primary'] = np.array(comp['specimen'], copy=True)
                    entry['inertial_torque_total_primary'] = np.array(comp['total'], copy=True)
                    entry['primary_torque_raw'] = np.array(comp['raw'], copy=True)
                    entry['primary_torque_corrected'] = np.array(comp['corrected'], copy=True)
                    m = active & np.isfinite(t)
                    if np.any(m):
                        entry['mean_xforce_stim'] = float(np.mean(ft[0, m]))
            except (ValueError, AttributeError):
                pass
            results.append(entry)

        self.force_length_results = results
        self.test_type = 'isometric'
        self.trial_records = list(results)
        self.h5_protocol_metadata.update({
            'test_type': 'isometric',
            'n_trials': int(len(results)),
            'rest_between_steps_s': float(gap_s),
            'reset_between_steps': bool(reset_steps),
            'pulse_width_ms': float(getattr(self, 'pulse_width_ms', 2.0)) if bool(is_stim) else None,
            'isometric_inter_step_interval_s': float(gap_s),
        })
        return results

    def run_force_length_series(
        self,
        initial_length,
        max_dist,
        num_steps,
        *,
        randomize=False,
        random_seed=None,
        ramp_duration_s=2.0,
        hold_duration_s=5.0,
        settle_before_stim_s=0.5,
        stim_duration_s=None,
        is_stim=False,
        stim_pulse_rate=None,
        stim_voltage=5.0,
        device_name=None,
        inter_step_interval_s=0.0,
    ):
        """
        Force–length / isometric-style protocol in **commanded motor angle (degrees)**.

        For each target between ``initial_length`` and ``max_dist``, the motor ramps, holds,
        stimulates (same carrier as :meth:`make_stimuli`), and acquires via :meth:`run`.
        For specimen units (strain, curvature), use :meth:`test_force_length` instead.
        """
        num_steps = int(num_steps)
        if num_steps < 1:
            raise ValueError("num_steps must be >= 1.")
        positions = np.linspace(float(initial_length), float(max_dist), num_steps)
        if bool(randomize) and positions.size > 1:
            rng = np.random.default_rng(None if random_seed is None else int(random_seed))
            rng.shuffle(positions)
        return self._run_force_length_steps(
            positions,
            ramp_duration_s=ramp_duration_s,
            hold_duration_s=hold_duration_s,
            settle_before_stim_s=settle_before_stim_s,
            stim_duration_s=stim_duration_s,
            inter_step_interval_s=float(inter_step_interval_s),
            is_stim=is_stim,
            stim_pulse_rate=stim_pulse_rate,
            stim_voltage=stim_voltage,
            device_name=device_name,
        )

    def _shuffled_step_order(self, n, randomize, random_seed=None, block_index=0):
        """Return a permutation of ``range(n)`` for step-order randomization (1.2).

        When ``randomize`` is False (or n < 2) the identity order is returned. With no
        ``random_seed`` the shuffle is non-deterministic (``seed=None``) and independent on each
        call (so block sequences get an independent order per block); a ``random_seed`` makes it
        reproducible while still differing per ``block_index``.
        """
        order = np.arange(int(n), dtype=int)
        if bool(randomize) and int(n) > 1:
            if random_seed is None:
                np.random.shuffle(order)
            else:
                rng = np.random.default_rng(int(random_seed) + int(block_index))
                rng.shuffle(order)
        return order

    def isometric(
        self,
        initial,
        final,
        num_steps,
        mode='strain',
        *,
        randomize=False,
        random_seed=None,
        stim_params=None,
        **kwargs,
    ):
        """
        Force–length steps where ``initial`` / ``final`` are in ``mode`` units (see :func:`convert_to_curvature`).

        Each step converts the target to curvature then to motor degrees (``dclamp``), then
        ramps, stimulates, and records force like :meth:`run_force_length_series`.

        Parameters
        ----------
        initial, final : float
            Endpoints for ``np.linspace(initial, final, num_steps)`` in ``mode`` units.
        num_steps : int
        mode : str
            ``strain`` (fraction), ``strain_pct``, ``angle`` (deg), or ``curvature`` (1/m).
        randomize : bool
            If True, randomize the sequence order after generating the linspace targets.
        random_seed : int or None
            Optional seed for deterministic randomization when ``randomize=True``.
        stim_params : dict or None
            Optional keys: ``ramp_duration_s``, ``hold_duration_s``, ``settle_before_stim_s``,
            ``stim_duration_s``, ``inter_step_interval_s`` (seconds of idle time after each step
            before the next begins; 0 = no pause), ``is_stim``, ``stim_pulse_rate``,
            ``stim_voltage``, ``device_name``, ``recruitment``, ``bilateral_mirror_motor``,
            ``isometric_mirror_target_left`` / ``isometric_mirror_target_right`` (same units as ``mode``;
            when both set with bilateral mirror, per-hold magnitudes instead of the step sweep scalar).
            If ``inter_step_interval_s`` is omitted, :attr:`isometric_inter_step_interval_s` on
            the instance is used (default 0).
        **kwargs
            Merged into ``stim_params`` (stim_params takes precedence).

        Returns
        -------
        list of dict
            Same structure as :meth:`run_force_length_series`, plus ``target_value_native``,
            ``curvature_1_per_m`` per step.
        """
        sp = {
            'ramp_duration_s': 2.0,
            'hold_duration_s': 5.0,
            'pre_baseline_s': 1.0,
            'post_baseline_s': 1.0,
            'stim_onset_s': None,
            'settle_before_stim_s': 0.5,
            'stim_duration_s': None,
            'inter_step_interval_s': None,
            'reset_between_steps': None,
            'is_stim': False,
            'stim_pulse_rate': None,
            'stim_voltage': 5.0,
            'device_name': None,
        }
        if stim_params:
            sp.update(stim_params)
        sp.update(kwargs)
        for _bk in ('pre_baseline_s', 'post_baseline_s'):
            if sp.get(_bk, None) is None and getattr(self, _bk, None) is not None:
                sp[_bk] = getattr(self, _bk)
        if sp.get('inter_step_interval_s', None) is None:
            sp['inter_step_interval_s'] = float(getattr(self, 'rest_between_steps_s', 2.0) or 0.0)
        if sp.get('reset_between_steps', None) is None:
            sp['reset_between_steps'] = bool(getattr(self, 'reset_between_steps', False))
        # Pulse width comes from the stim panel; mirror onto self so the pulse carrier uses it (1.4).
        if sp.get('pulse_width_ms', None) is not None:
            self.pulse_width_ms = float(sp['pulse_width_ms'])
        for k in ('isometric_mirror_target_left', 'isometric_mirror_target_right'):
            if k not in sp and getattr(self, k, None) is not None:
                sp[k] = getattr(self, k)
        # Negative stim onset for isometric is measured back into the pre-hold ramp: the ramp is
        # the only time before the hold (active segment) starts, so it bounds how negative onset
        # can go (stim may begin partway through the ramp, never before t=0).
        self._validate_stim_timing_for_steps(
            sp,
            test_type='isometric',
            num_steps=int(num_steps),
            pre_hold_at_start_s=float(sp.get('ramp_duration_s', 2.0)),
            segment_duration_s=float(sp.get('hold_duration_s', 5.0)),
        )
        # Base linspace in generated (un-shuffled) order. The non-block path below shuffles once
        # and the block path shuffles each block independently (1.2 step-order randomization).
        vals = np.linspace(float(initial), float(final), int(num_steps))
        seq_idx = np.arange(vals.size, dtype=int)
        dc = self._effective_dclamp_mm()
        xw = getattr(self, 'xsec_width', None)
        md = float(getattr(self, 'target_muscle_depth_mm', 0.0) or 0.0)
        if dc is None:
            raise ValueError(
                "isometric requires clamp spacing (mm): set `dclamp` or `test_segment_length_mm` on the Bender instance."
            )
        kappa = convert_to_curvature(vals, mode, dclamp_mm=dc, xsec_width_mm=xw, target_muscle_depth_mm=md)
        targets_deg_raw = np.rad2deg(kappa * (float(dc) / 1000.0))
        block_seq = self._normalize_block_sequence(
            sp.get('block_sequence', getattr(self, 'block_sequence', None))
        )
        rec = self._normalize_recruitment(
            sp.get('recruitment', sp.get('lateral_mode', getattr(self, 'recruitment', 'bilateral_simultaneous')))
        )
        mirror_bm = bool(sp.get('bilateral_mirror_motor', getattr(self, 'bilateral_mirror_motor', False)))
        if block_seq:
            if mirror_bm:
                raise ValueError(
                    "bilateral_mirror_motor cannot be used with block_sequence; set direction per block instead."
                )
            left_v = float(sp.get('left_stim_voltage', getattr(self, 'left_stim_voltage', 5.0)))
            right_v = float(sp.get('right_stim_voltage', getattr(self, 'right_stim_voltage', 5.0)))
            reset_ramp = float(
                sp.get(
                    'block_reset_ramp_duration_s',
                    getattr(self, 'block_reset_ramp_duration_s', sp.get('ramp_duration_s', 2.0)),
                )
            )
            if not np.isfinite(reset_ramp) or reset_ramp < 0:
                raise ValueError(
                    f"block_reset_ramp_duration_s must be finite and >= 0; got {reset_ramp!r}."
                )
            self._validate_block_sequence_voltages(block_seq, left_v, right_v)
            dev = sp.get('device_name', None) or getattr(self, 'device_name', None)
            all_out = []
            last_deg = 0.0
            global_step = 0
            step_order_blocks = []
            base_targets = np.abs(targets_deg_raw)
            for bi, block in enumerate(block_seq):
                last_deg = self._run_neutral_reset_segment(last_deg, reset_ramp, dev)
                dir_idx = self._lateral_index_for_block_direction(block['direction'])
                ms = self.motor_command_sign_for_bend_toward_index(dir_idx)
                order = self._shuffled_step_order(base_targets.size, randomize, random_seed, block_index=bi)
                step_order_blocks.append([int(x) for x in order])
                targets_block = base_targets[order] * ms
                block_out = self._run_force_length_steps(
                    targets_block,
                    ramp_duration_s=float(sp.get('ramp_duration_s', 2.0)),
                    hold_duration_s=float(sp.get('hold_duration_s', 5.0)),
                    pre_baseline_s=float(sp.get('pre_baseline_s', 1.0) or 0.0),
                    post_baseline_s=float(sp.get('post_baseline_s', 1.0) or 0.0),
                    stim_onset_s=sp.get('stim_onset_s'),
                    settle_before_stim_s=float(sp.get('settle_before_stim_s', 0.5)),
                    stim_duration_s=sp.get('stim_duration_s'),
                    inter_step_interval_s=float(sp.get('inter_step_interval_s', 0.0) or 0.0),
                    reset_between_steps=bool(sp.get('reset_between_steps', False)),
                    is_stim=bool(sp.get('is_stim', False)),
                    stim_pulse_rate=sp.get('stim_pulse_rate', None),
                    stim_voltage=float(sp.get('stim_voltage', 5.0)),
                    device_name=dev,
                    ramp_from_deg=0.0,
                    stim_sides=block['stim_sides'],
                    left_stim_voltage=left_v,
                    right_stim_voltage=right_v,
                    block_direction=block['direction'],
                    block_index=bi,
                    bilateral_mirror_motor=False,
                )
                for j, e in enumerate(block_out):
                    orig = int(order[j])
                    e['target_value_native'] = float(vals[orig])
                    e['curvature_1_per_m'] = float(kappa[orig])
                    e['sequence_index'] = orig
                    e['step_index'] = global_step
                    e['trial_index'] = global_step
                    e['cycle_index'] = global_step
                    global_step += 1
                last_deg = float(block_out[-1]['angle_cmd'][-1])
                all_out.extend(block_out)
            # Return the motor to neutral once after the final block: block resets only
            # fire BEFORE each block, so without this the motor is left at the last hold
            # angle after the run ends.
            last_deg = self._run_neutral_reset_segment(last_deg, reset_ramp, dev)
            uidx = None
            seq_frac = float(
                sp.get('bilateral_sequential_left_frac', getattr(self, 'bilateral_sequential_left_frac', 0.5))
            )
            ml_n = mr_n = mhl = mhr = None
            meta = {
                'recruitment': rec,
                'bilateral_mirror_motor': False,
                'bilateral_sequential_left_frac': seq_frac,
                'block_sequence': block_seq,
                'randomize_step_order': bool(randomize),
                'step_order': step_order_blocks,
                'left_stim_voltage': left_v,
                'right_stim_voltage': right_v,
                'block_reset_ramp_duration_s': reset_ramp,
                'specimen_side_index_left': int(self.specimen_side_index_left),
                'specimen_side_index_right': int(self.specimen_side_index_right),
                'specimen_lateral_index_on_positive_motor_side': int(
                    getattr(self, 'specimen_lateral_index_on_positive_motor_side', -1)
                ),
                'motor_positive_bend_toward_lateral_index': int(self.motor_positive_bend_lateral_index()),
                'unilateral_posture_lateral_index': uidx,
                'n_trials': int(len(all_out)),
            }
            self.h5_protocol_metadata.update(meta)
            self.force_length_results = all_out
            self.test_type = 'isometric'
            self.trial_records = list(all_out)
            return all_out
        # Single (non-block) path: shuffle the step order once if requested (1.2).
        order = self._shuffled_step_order(vals.size, randomize, random_seed)
        seq_idx = order
        vals = vals[order]
        kappa = kappa[order]
        targets_deg_raw = targets_deg_raw[order]
        targets_deg = targets_deg_raw.copy()
        rec = self._recruitment_with_bilateral_mirror_motor(rec, mirror_bm)
        uidx = self.recruitment_unilateral_lateral_index(rec)
        if uidx is not None:
            targets_deg = targets_deg * self.motor_command_sign_for_bend_toward_index(uidx)
        seq_frac = float(sp.get('bilateral_sequential_left_frac', getattr(self, 'bilateral_sequential_left_frac', 0.5)))
        ml_n = sp.get('isometric_mirror_target_left')
        mr_n = sp.get('isometric_mirror_target_right')
        mhl = None
        mhr = None
        if mirror_bm and ml_n is not None and mr_n is not None:
            try:
                fL, fR = float(ml_n), float(mr_n)
                if np.isfinite(fL) and np.isfinite(fR):
                    kL = convert_to_curvature(
                        np.asarray([fL]), mode, dclamp_mm=float(dc), xsec_width_mm=xw, target_muscle_depth_mm=md
                    )
                    kR = convert_to_curvature(
                        np.asarray([fR]), mode, dclamp_mm=float(dc), xsec_width_mm=xw, target_muscle_depth_mm=md
                    )
                    mhl = float(
                        np.abs(np.rad2deg(float(np.asarray(kL, dtype=float).reshape(-1)[0]) * (float(dc) / 1000.0)))
                    )
                    mhr = float(
                        np.abs(np.rad2deg(float(np.asarray(kR, dtype=float).reshape(-1)[0]) * (float(dc) / 1000.0)))
                    )
            except (TypeError, ValueError):
                mhl = None
                mhr = None
        out = self._run_force_length_steps(
            targets_deg,
            ramp_duration_s=float(sp.get('ramp_duration_s', 2.0)),
            hold_duration_s=float(sp.get('hold_duration_s', 5.0)),
            pre_baseline_s=float(sp.get('pre_baseline_s', 1.0) or 0.0),
            post_baseline_s=float(sp.get('post_baseline_s', 1.0) or 0.0),
            stim_onset_s=sp.get('stim_onset_s'),
            settle_before_stim_s=float(sp.get('settle_before_stim_s', 0.5)),
            stim_duration_s=sp.get('stim_duration_s'),
            inter_step_interval_s=float(sp.get('inter_step_interval_s', 0.0) or 0.0),
            reset_between_steps=bool(sp.get('reset_between_steps', False)),
            is_stim=bool(sp.get('is_stim', False)),
            stim_pulse_rate=sp.get('stim_pulse_rate', None),
            stim_voltage=float(sp.get('stim_voltage', 5.0)),
            device_name=sp.get('device_name', None),
            recruitment=rec,
            bilateral_mirror_motor=mirror_bm,
            bilateral_sequential_left_frac=seq_frac,
            mirror_hold_deg_left=mhl,
            mirror_hold_deg_right=mhr,
        )
        # Return the motor to neutral once after the final trial. The simple isometric
        # path has no block boundaries, so without this the motor is left at the last
        # hold angle when the run ends.
        if out:
            reset_ramp = float(
                sp.get(
                    'block_reset_ramp_duration_s',
                    getattr(self, 'block_reset_ramp_duration_s', sp.get('ramp_duration_s', 2.0)),
                )
            )
            if not np.isfinite(reset_ramp) or reset_ramp < 0:
                raise ValueError(
                    f"block_reset_ramp_duration_s must be finite and >= 0; got {reset_ramp!r}."
                )
            dev = sp.get('device_name', None) or getattr(self, 'device_name', None)
            self._run_neutral_reset_segment(float(out[-1]['angle_cmd'][-1]), reset_ramp, dev)
        for i, e in enumerate(out):
            e['target_value_native'] = float(vals[i])
            e['curvature_1_per_m'] = float(kappa[i])
            e['sequence_index'] = int(seq_idx[i])
        meta = {
            'recruitment': rec,
            'bilateral_mirror_motor': mirror_bm,
            'bilateral_sequential_left_frac': seq_frac,
            'randomize_step_order': bool(randomize),
            'step_order': [int(x) for x in seq_idx],
            'specimen_side_index_left': int(self.specimen_side_index_left),
            'specimen_side_index_right': int(self.specimen_side_index_right),
            'specimen_lateral_index_on_positive_motor_side': int(
                getattr(self, 'specimen_lateral_index_on_positive_motor_side', -1)
            ),
            'motor_positive_bend_toward_lateral_index': int(self.motor_positive_bend_lateral_index()),
            'unilateral_posture_lateral_index': uidx,
        }
        if ml_n is not None:
            meta['isometric_mirror_target_left'] = ml_n
        if mr_n is not None:
            meta['isometric_mirror_target_right'] = mr_n
        if mhl is not None and mhr is not None:
            meta['isometric_mirror_hold_deg_left'] = mhl
            meta['isometric_mirror_hold_deg_right'] = mhr
        self.h5_protocol_metadata.update(meta)
        return out

    def test_force_length(self, initial, final, num_steps, mode='strain', stim_params=None, **kwargs):
        """Backward-compatible alias for :meth:`isometric`."""
        return self.isometric(
            initial,
            final,
            num_steps,
            mode=mode,
            stim_params=stim_params,
            **kwargs,
        )

    def _timeline_prehold_isovelocity(self, theta_start_deg, vel_deg_per_s, pre_hold_s, iso_duration_s, daq_hz):
        """Hold at fixed angle, then bend at constant ``vel_deg_per_s`` for ``iso_duration_s``."""
        dt = 1.0 / float(daq_hz)
        th0 = float(theta_start_deg)
        v = float(vel_deg_per_s)
        ph = float(pre_hold_s)
        iso = float(iso_duration_s)
        if not np.isfinite(iso) or iso <= 0:
            raise ValueError("iso_duration_s must be finite and > 0.")
        if ph < 0:
            raise ValueError("pre_hold_s must be non-negative.")
        # Use a sub-threshold velocity during the pre-hold so the direction bit is
        # pre-set to match the upcoming ramp. An exact zero would always encode as
        # REVERSE (velhi <= 0), causing 300 ms of wrong coil energisation before
        # every rightward ramp. 1e-10 deg/s is numerically zero for position but
        # gives the sign needed by make_motor_stepper_pulses.
        prehold_vel = np.copysign(1e-10, v) if v != 0.0 else 0.0
        if ph > 0:
            n1 = max(2, int(round(ph * daq_hz)) + 1)
            t1 = np.linspace(0.0, ph, n1)
            a1 = np.full(n1, th0)
            w1 = np.full(n1, prehold_vel)
        else:
            t1 = np.array([0.0])
            a1 = np.array([th0])
            w1 = np.array([prehold_vel])
        n2 = max(2, int(round(iso * daq_hz)) + 1)
        t_edge = ph
        t2 = np.linspace(t_edge, t_edge + iso, n2)[1:]
        if t2.size < 1:
            t2 = np.array([t_edge + iso])
        a2 = th0 + v * (t2 - t_edge)
        w2 = np.full(t2.size, v)
        # Velocity-angle guard (2.1): stop the constant-velocity ramp the moment the commanded
        # angle reaches the hard limit. Open-loop stepper ⇒ commanded angle is the encoder angle.
        # The caller's post-baseline ramp then returns the motor to 0° (controlled deceleration).
        guard_triggered = False
        guard_angle_deg = float('nan')
        limit = float(ISOVELOCITY_ANGLE_LIMIT_DEG)
        over = np.nonzero(np.abs(a2) >= limit)[0]
        if over.size:
            k = int(over[0])
            guard_triggered = True
            guard_angle_deg = float(a2[k])
            t2 = t2[:k + 1]
            a2 = a2[:k + 1]
            w2 = w2[:k + 1]
            import warnings
            warnings.warn(
                f"Isovelocity angle guard: |angle|={abs(guard_angle_deg):.3g}° reached the "
                f"{limit:.3g}° limit; ramp stopped and motor will return to 0°.",
                UserWarning,
                stacklevel=2,
            )
        t = np.concatenate([t1, t2])
        angle = np.concatenate([a1, a2])
        anglevel = np.concatenate([w1, w2])
        if t.size < 2:
            raise ValueError("_timeline_prehold_isovelocity: timeline too short; increase durations or daq rate.")
        return t, angle, anglevel, float(t_edge), guard_triggered, guard_angle_deg

    def _isovelocity_one_block(
        self,
        theta_start_deg,
        vel_deg_per_s,
        *,
        pre_hold_s,
        iso_duration_s,
        stim_onset_s=None,
        stim_duration_s=None,
        is_stim,
        spr,
        stim_voltage,
        daq_hz,
        recruitment,
        sequential_left_frac,
        mirror_stim_side,
        stim_sides=None,
        left_stim_voltage=None,
        right_stim_voltage=None,
        settle_before_stim_s=None,
        pre_iso_stim_duration_s=None,
        post_baseline_s=1.0,
    ):
        """
        Build one isovelocity timeline + stim commands.

        ``mirror_stim_side``: None -> use ``recruitment`` routing; ``'left'`` / ``'right'`` -> that side only
        (for bilateral sequential + bilateral_mirror_motor two-segment runs).

        ``stim_onset_s`` is signed seconds relative to constant-velocity segment start (``t_iso0``).
        Legacy ``settle_before_stim_s`` / ``pre_iso_stim_duration_s`` kwargs are migrated when onset omitted.
        """
        stim_onset_s, stim_duration_s = self._resolve_stim_onset_duration_s(
            {
                'stim_onset_s': stim_onset_s,
                'stim_duration_s': stim_duration_s,
                'settle_before_stim_s': settle_before_stim_s,
                'pre_iso_stim_duration_s': pre_iso_stim_duration_s,
            },
            segment_duration_s=iso_duration_s,
        )
        t, angle, anglevel, t_iso0, guard_triggered, guard_angle_deg = self._timeline_prehold_isovelocity(
            float(theta_start_deg), float(vel_deg_per_s), float(pre_hold_s), float(iso_duration_s), daq_hz
        )
        t_seg_end = t_iso0 + float(iso_duration_s)
        # Post-stimulus baseline: ramp the motor back to neutral (0 deg) with stim off; recording
        # continues through the return ramp so the trial is one flat, continuous time series.
        post_b = max(0.0, float(post_baseline_s))
        if post_b > 0:
            final_ang = float(angle[-1])
            n_post = max(2, int(round(post_b * float(daq_hz))) + 1)
            t_post = float(t[-1]) + np.linspace(0.0, post_b, n_post)[1:]
            if t_post.size < 1:
                t_post = np.array([float(t[-1]) + post_b])
            frac = (t_post - float(t[-1])) / post_b
            ang_post = final_ang + (0.0 - final_ang) * frac
            w_post = np.full(t_post.size, (0.0 - final_ang) / post_b)
            t = np.concatenate([t, t_post])
            angle = np.concatenate([angle, ang_post])
            anglevel = np.concatenate([anglevel, w_post])
            # End-of-trial buffer: after the controlled return to neutral, hold flat at the final
            # angle (0 deg) with zero velocity for ``pre_hold_s`` so the trial does not end abruptly
            # the instant motion stops. Mirrors the dynamic protocol's post-motion hold (waitafter)
            # and matches the isovelocity start hold (pre_hold), so each trial is symmetric. Only the
            # trial-ending segment returns to neutral (post_baseline_s > 0); intermediate
            # bilateral-mirror sub-segments (post_baseline_s == 0) get no spurious hold.
            end_hold_s = max(0.0, float(pre_hold_s))
            if end_hold_s > 0:
                hold_ang = float(angle[-1])
                n_hold = max(2, int(round(end_hold_s * float(daq_hz))) + 1)
                t_hold = float(t[-1]) + np.linspace(0.0, end_hold_s, n_hold)[1:]
                if t_hold.size < 1:
                    t_hold = np.array([float(t[-1]) + end_hold_s])
                a_hold = np.full(t_hold.size, hold_ang)
                w_hold = np.zeros(t_hold.size)
                t = np.concatenate([t, t_hold])
                angle = np.concatenate([angle, a_hold])
                anglevel = np.concatenate([anglevel, w_hold])
        t_post0 = t_seg_end
        t_post1 = t_seg_end + post_b
        t_stim0 = t_iso0 + float(stim_onset_s)
        t_stim1 = t_stim0 + float(stim_duration_s)
        # Stim is confined to the active (constant-velocity) segment; it never bleeds into the
        # post-baseline return ramp.
        t_stim1 = min(t_stim1, t_seg_end + 1e-9)
        active = (t >= t_stim0) & (t < t_stim1) if is_stim else np.zeros_like(t, dtype=bool)
        if is_stim and np.any(active):
            if stim_sides is not None:
                lv = float(left_stim_voltage if left_stim_voltage is not None else stim_voltage)
                rv = float(right_stim_voltage if right_stim_voltage is not None else stim_voltage)
                s1, s2 = self._route_stim_sides_volts(t, active, spr, stim_sides, lv, rv)
            elif mirror_stim_side == 'left':
                s1 = np.zeros_like(t)
                s2 = np.zeros_like(t)
                pulse = self._pulse_carrier_volts(t, active, spr, stim_voltage)
                self._deposit_stim_on_side(pulse, 'left', s1, s2)
            elif mirror_stim_side == 'right':
                s1 = np.zeros_like(t)
                s2 = np.zeros_like(t)
                pulse = self._pulse_carrier_volts(t, active, spr, stim_voltage)
                self._deposit_stim_on_side(pulse, 'right', s1, s2)
            else:
                pulse = self._pulse_carrier_volts(t, active, spr, stim_voltage)
                s1, s2 = self._route_recruitment_stim(pulse, recruitment, sequential_left_frac=sequential_left_frac)
        else:
            s1 = np.zeros_like(t)
            s2 = np.zeros_like(t)
        return {
            't': t,
            'angle': angle,
            'anglevel': anglevel,
            't_iso0': t_iso0,
            't_pre0': 0.0,
            't_pre1': t_iso0,
            't_stim0': t_stim0,
            't_stim1': t_stim1,
            't_post0': t_post0,
            't_post1': t_post1,
            'stim_onset_s': float(stim_onset_s),
            'stim_duration_s': float(stim_duration_s),
            'iso_t1': t_seg_end,
            's1': s1,
            's2': s2,
            'active': active,
            'guard_triggered': bool(guard_triggered),
            'guard_angle_deg': float(guard_angle_deg),
        }

    def _run_isovelocity_steps(
        self,
        theta_start_deg,
        velocities_deg_per_s,
        *,
        pre_hold_s=0.3,
        iso_duration_s=0.2,
        post_baseline_s=1.0,
        stim_onset_s=None,
        settle_before_stim_s=0.02,
        pre_iso_stim_duration_s=0.0,
        stim_duration_s=None,
        inter_step_interval_s=0.0,
        reset_between_steps=False,
        reset_ramp_duration_s=None,
        is_stim=False,
        stim_pulse_rate=None,
        stim_voltage=5.0,
        device_name=None,
        recruitment=None,
        bilateral_mirror_motor=False,
        bilateral_sequential_left_frac=0.5,
        block_direction=None,
        stim_sides=None,
        left_stim_voltage=None,
        right_stim_voltage=None,
        block_index=None,
    ):
        """Acquire one DAQ segment per commanded constant angular velocity (deg/s)."""
        vels = np.atleast_1d(np.asarray(velocities_deg_per_s, dtype=float)).reshape(-1)
        n = int(vels.size)
        if n < 1:
            raise ValueError("velocities_deg_per_s must contain at least one value.")
        if stim_sides is not None and bool(bilateral_mirror_motor):
            raise ValueError("bilateral_mirror_motor cannot be used with per-block stim_sides.")
        use_block_stim = stim_sides is not None
        use_block_direction = block_direction is not None
        dev = device_name if device_name is not None else getattr(self, 'device_name', None)
        if dev is None:
            raise ValueError("_run_isovelocity_steps requires device_name or self.device_name.")
        spr = float(stim_pulse_rate) if stim_pulse_rate is not None else float(
            getattr(self, 'stim_pulse_rate', 75.0)
        )
        daq_hz = float(self.daq_ai_sample_rate_hz)
        results = []
        if use_block_stim:
            lv = float(left_stim_voltage if left_stim_voltage is not None else stim_voltage)
            rv = float(right_stim_voltage if right_stim_voltage is not None else stim_voltage)
            sides_norm = self._normalize_stim_sides(stim_sides)
            rec = self._stim_sides_to_recruitment(sides_norm)
        else:
            lv = float(stim_voltage)
            rv = float(stim_voltage)
            sides_norm = None
            rec = self._normalize_recruitment(
                recruitment if recruitment is not None else getattr(self, 'recruitment', 'bilateral_simultaneous')
            )
        bm = bool(bilateral_mirror_motor) and not use_block_stim and not use_block_direction
        rec = self._recruitment_with_bilateral_mirror_motor(rec, bm)
        mirror = bm and rec == 'bilateral_sequential'
        if use_block_direction:
            dir_idx = self._lateral_index_for_block_direction(block_direction)
            block_dir_sign = self.motor_command_sign_for_bend_toward_index(dir_idx)
        else:
            block_dir_sign = None
        seq_frac = float(
            bilateral_sequential_left_frac
            if bilateral_sequential_left_frac is not None
            else getattr(self, 'bilateral_sequential_left_frac', 0.5)
        )
        self.set_input_channels(input_channels=self.input_channels, input_channel_names=self.input_channel_names)
        self.set_stim_channels(*self.stim_channels)
        th0_fixed = float(theta_start_deg)
        post_b = max(0.0, float(post_baseline_s))
        iso_kw = dict(
            pre_hold_s=pre_hold_s,
            iso_duration_s=iso_duration_s,
            stim_onset_s=stim_onset_s,
            stim_duration_s=stim_duration_s,
            settle_before_stim_s=settle_before_stim_s,
            pre_iso_stim_duration_s=pre_iso_stim_duration_s,
            is_stim=is_stim,
            spr=spr,
            stim_voltage=stim_voltage,
            daq_hz=daq_hz,
            recruitment=rec,
            sequential_left_frac=seq_frac,
        )
        if use_block_stim:
            iso_kw.update(
                stim_sides=sides_norm,
                left_stim_voltage=lv,
                right_stim_voltage=rv,
            )

        gap_s = float(inter_step_interval_s)
        if not np.isfinite(gap_s) or gap_s < 0:
            raise ValueError(f"inter_step_interval_s must be finite and >= 0; got {inter_step_interval_s!r}.")
        reset_steps = bool(reset_between_steps)
        reset_ramp = float(
            reset_ramp_duration_s if reset_ramp_duration_s is not None
            else getattr(self, 'block_reset_ramp_duration_s', 2.0)
        )
        prev_end_deg = float(theta_start_deg)

        for i in range(n):
            if i > 0 and gap_s > 0:
                # Hold the motor at its current position and rest before the next velocity step.
                time.sleep(gap_s)
            if i > 0 and reset_steps:
                # Reset to resting length: drive the motor back to 0° from the previous step's
                # end angle and wait for the move to complete before the next step.
                self._run_neutral_reset_segment(prev_end_deg, reset_ramp, dev)
            v_mag = abs(float(vels[i]))
            th0 = th0_fixed
            if mirror:
                v1 = v_mag * self.motor_command_sign_for_bend_toward_index(self.specimen_side_index_left)
                v2 = v_mag * self.motor_command_sign_for_bend_toward_index(self.specimen_side_index_right)
                d1 = self._isovelocity_one_block(
                    th0, v1,
                    mirror_stim_side='left',
                    post_baseline_s=0.0,
                    **iso_kw,
                )
                th_mid = float(d1['angle'][-1])
                d2 = self._isovelocity_one_block(
                    th_mid, v2,
                    mirror_stim_side='right',
                    post_baseline_s=post_b,
                    **iso_kw,
                )
                off = float(d1['t'][-1])
                t = np.concatenate([d1['t'], d2['t'][1:] + off])
                angle = np.concatenate([d1['angle'], d2['angle'][1:]])
                anglevel = np.concatenate([d1['anglevel'], d2['anglevel'][1:]])
                s1 = np.concatenate([d1['s1'], d2['s1'][1:]])
                s2 = np.concatenate([d1['s2'], d2['s2'][1:]])
                t_iso0 = d1['t_iso0']
                t_pre0 = d1['t_pre0']
                t_pre1 = d1['t_pre1']
                t_stim0 = d1['t_stim0']
                t_stim1 = d2['t_stim1'] + off
                active = np.concatenate([d1['active'], d2['active'][1:]])
                v_report = float(vels[i])
                iso_t1_end = d2['iso_t1'] + off
                t_post0 = d2['t_post0'] + off
                t_post1 = d2['t_post1'] + off
                step_guard_triggered = bool(d1['guard_triggered'] or d2['guard_triggered'])
                step_guard_angle = (
                    d1['guard_angle_deg'] if d1['guard_triggered'] else d2['guard_angle_deg']
                )
            else:
                if block_dir_sign is not None:
                    v_sign = v_mag * block_dir_sign
                else:
                    v_sign = float(vels[i])
                    uidx = self.recruitment_unilateral_lateral_index(rec)
                    if uidx is not None:
                        v_sign = v_mag * self.motor_command_sign_for_bend_toward_index(uidx)
                d0 = self._isovelocity_one_block(
                    th0, v_sign,
                    mirror_stim_side=None,
                    **iso_kw,
                )
                t = d0['t']
                angle = d0['angle']
                anglevel = d0['anglevel']
                t_iso0 = d0['t_iso0']
                t_pre0 = d0['t_pre0']
                t_pre1 = d0['t_pre1']
                t_stim0 = d0['t_stim0']
                t_stim1 = d0['t_stim1']
                s1, s2 = d0['s1'], d0['s2']
                active = d0['active']
                v_report = float(v_sign)
                iso_t1_end = d0['iso_t1']
                t_post0 = d0['t_post0']
                t_post1 = d0['t_post1']
                step_guard_triggered = bool(d0['guard_triggered'])
                step_guard_angle = d0['guard_angle_deg']

            prev_end_deg = float(angle[-1])
            self.record_motor_signal(t, angle, anglevel, tnorm=np.zeros_like(t))
            self.record_stim_signal(s1, s2)
            self.make_motor_stepper_pulses(
                daq_ao_do_sample_rate_hz=self.daq_ao_do_sample_rate_hz,
                motor_gear_ratio=self.motor_gear_ratio,
                motor_full_steps_per_rev=self.motor_full_steps_per_rev,
            )
            self.aidata = self.run(device_name=dev)

            entry = {
                'test_type': 'isovelocity',
                'step_index': i,
                'trial_index': i,
                'cycle_index': i,
                'recruitment': rec,
                'unilateral_posture_lateral_index': self.recruitment_unilateral_lateral_index(rec),
                'motor_positive_bend_toward_lateral_index': int(self.motor_positive_bend_lateral_index()),
                'bilateral_mirror_motor': bool(mirror),
                'velocity_deg_s': v_report,
                'theta_start_deg': th0_fixed,
                't': np.array(t, copy=True),
                'angle_cmd': np.array(angle, copy=True),
                'anglevel_cmd': np.array(anglevel, copy=True),
                'tnorm': np.zeros_like(t),
                'S1stimcmd': np.array(s1, copy=True),
                'S2stimcmd': np.array(s2, copy=True),
                'aidata': np.array(self.aidata, copy=True),
                'angle_measured': np.array(self.angle_measured, copy=True),
                'pre_stim_t0': t_pre0,
                'pre_stim_t1': t_pre1,
                'stim_t0': t_stim0,
                'stim_t1': t_stim1,
                'iso_t0': t_iso0,
                'iso_t1': float(iso_t1_end),
                't_pre_baseline_start': float(t_pre0),
                't_pre_baseline_end': float(t_pre1),
                't_active_start': float(t_iso0),
                't_active_end': float(iso_t1_end),
                't_post_baseline_start': float(t_post0),
                't_post_baseline_end': float(t_post1),
                'forcetorque': None,
                'mean_xforce_stim': None,
                'guard_triggered': bool(step_guard_triggered),
                'guard_angle_deg': float(step_guard_angle),
            }
            if use_block_stim:
                self._tag_block_trial_metadata(
                    entry,
                    block_index=block_index if block_index is not None else 0,
                    block_direction=block_direction or sides_norm,
                    block_stim_sides=sides_norm,
                    left_stim_voltage=lv,
                    right_stim_voltage=rv,
                )
            if mirror:
                entry['velocity_seg1_deg_s'] = float(
                    v_mag * self.motor_command_sign_for_bend_toward_index(self.specimen_side_index_left)
                )
                entry['velocity_seg2_deg_s'] = float(
                    v_mag * self.motor_command_sign_for_bend_toward_index(self.specimen_side_index_right)
                )
            sg_names = ['xForce', 'yForce', 'zForce', 'xTorque', 'yTorque', 'zTorque']
            try:
                idx = [self.input_channel_names.index(n) for n in sg_names if n in self.input_channel_names]
                if len(idx) == 6:
                    ft = self.apply_calibration_forcetorque(self.aidata[idx, :])
                    entry['forcetorque'] = ft
                    entry['forcetorque_raw'] = np.array(ft, copy=True)
                    idx_t = self._primary_torque_index()
                    comp = self._compute_primary_torque_components(ft[idx_t, :], anglevel)
                    ft_corr = np.array(ft, copy=True)
                    ft_corr[idx_t, :len(comp['corrected'])] = comp['corrected']
                    entry['forcetorque_corrected'] = ft_corr
                    entry['inertial_torque_system_primary'] = np.array(comp['system'], copy=True)
                    entry['inertial_torque_specimen_primary'] = np.array(comp['specimen'], copy=True)
                    entry['inertial_torque_total_primary'] = np.array(comp['total'], copy=True)
                    entry['primary_torque_raw'] = np.array(comp['raw'], copy=True)
                    entry['primary_torque_corrected'] = np.array(comp['corrected'], copy=True)
                    m = active & np.isfinite(t)
                    if np.any(m):
                        entry['mean_xforce_stim'] = float(np.mean(ft[0, m]))
            except (ValueError, AttributeError):
                pass
            results.append(entry)

        self.isovelocity_results = results
        self.test_type = 'isovelocity'
        self.trial_records = list(results)
        self.h5_protocol_metadata.update({
            'test_type': 'isovelocity',
            'n_trials': int(len(results)),
            'rest_between_steps_s': float(gap_s),
            'reset_between_steps': bool(reset_steps),
            'pulse_width_ms': float(getattr(self, 'pulse_width_ms', 2.0)) if bool(is_stim) else None,
        })
        return results

    def isovelocity(
        self,
        min_vel,
        max_vel,
        starting_strain,
        num_steps,
        *,
        starting_strain_mode='strain',
        randomize=False,
        random_seed=None,
        iso_duration_s=0.2,
        pre_hold_s=0.3,
        stim_params=None,
        **kwargs,
    ):
        """
        Isovelocity sweep: hold at a strain-defined angle, then command constant angular velocity
        (deg/s) for a short interval, with optional stimulation during the iso segment.

        Velocities are ``np.linspace(min_vel, max_vel, num_steps)`` in the units given by
        ``isovelocity_velocity_mode`` (default **angle_vel** = deg/s), then converted to motor deg/s.
        Each trial starts again from the same starting angle.

        Parameters
        ----------
        min_vel, max_vel : float
            Velocity range in ``isovelocity_velocity_mode`` units (converted to motor deg/s).
        starting_strain : float
            Interpreted with ``starting_strain_mode`` (``strain`` fraction, ``strain_pct``, etc.).
        num_steps : int
        starting_strain_mode : str
            Passed to :func:`convert_to_curvature` for the initial posture.
        isovelocity_velocity_mode : str
            Rate mode for min/max velocity (``angle_vel``, ``strain_rate``, etc.).
        randomize : bool
            If True, randomize the order of the generated velocity steps.
        random_seed : int or None
            Optional seed for deterministic randomization when ``randomize=True``.
        iso_duration_s : float
            Duration of the constant-velocity bend (seconds).
        pre_hold_s : float
            Quiet hold at the starting angle before the iso segment (seconds).
        stim_params : dict or None
            Optional: ``settle_before_stim_s``, ``stim_duration_s``, ``is_stim``,
            ``pre_iso_stim_duration_s``, ``stim_pulse_rate``, ``stim_voltage``,
            ``device_name``, ``iso_duration_s``, ``pre_hold_s``.
        **kwargs
            Merged into ``stim_params`` (``stim_params`` wins on duplicate keys).

        Returns
        -------
        list of dict
            Per-step results (see :meth:`_run_isovelocity_steps`); each includes ``velocity_deg_s``,
            ``starting_strain``, ``curvature_1_per_m`` at start.
        """
        sp = {
            'iso_duration_s': iso_duration_s,
            'pre_hold_s': pre_hold_s,
            'post_baseline_s': 1.0,
            'stim_onset_s': None,
            'settle_before_stim_s': 0.02,
            'pre_iso_stim_duration_s': 0.0,
            'stim_duration_s': None,
            'inter_step_interval_s': None,
            'reset_between_steps': None,
            'is_stim': False,
            'stim_pulse_rate': None,
            'stim_voltage': 5.0,
            'device_name': None,
        }
        if stim_params:
            sp.update(stim_params)
        sp.update(kwargs)
        if sp.get('post_baseline_s', None) is None and getattr(self, 'post_baseline_s', None) is not None:
            sp['post_baseline_s'] = getattr(self, 'post_baseline_s')
        if sp.get('inter_step_interval_s', None) is None:
            sp['inter_step_interval_s'] = float(getattr(self, 'rest_between_steps_s', 2.0) or 0.0)
        if sp.get('reset_between_steps', None) is None:
            sp['reset_between_steps'] = bool(getattr(self, 'reset_between_steps', False))
        # Pulse width comes from the stim panel; mirror onto self so the pulse carrier uses it (1.4).
        if sp.get('pulse_width_ms', None) is not None:
            self.pulse_width_ms = float(sp['pulse_width_ms'])

        self._validate_stim_timing_for_steps(
            sp,
            test_type='isovelocity',
            num_steps=int(num_steps),
            pre_hold_at_start_s=float(sp.get('pre_hold_s', pre_hold_s)),
            segment_duration_s=float(sp.get('iso_duration_s', iso_duration_s)),
        )

        dc = self._effective_dclamp_mm()
        xw = getattr(self, 'xsec_width', None)
        md = float(getattr(self, 'target_muscle_depth_mm', 0.0) or 0.0)
        if dc is None:
            raise ValueError(
                "isovelocity requires clamp spacing (mm): set `dclamp` or `test_segment_length_mm` on the Bender instance."
            )
        k0 = convert_to_curvature(
            float(starting_strain),
            starting_strain_mode,
            dclamp_mm=float(dc),
            xsec_width_mm=xw,
            target_muscle_depth_mm=md,
        )
        k0 = float(np.asarray(k0).reshape(-1)[0])
        theta0_raw = float(np.rad2deg(k0 * (float(dc) / 1000.0)))
        block_seq = self._normalize_block_sequence(
            sp.get('block_sequence', getattr(self, 'block_sequence', None))
        )
        rec_iso = self._normalize_recruitment(
            sp.get('recruitment', sp.get('lateral_mode', getattr(self, 'recruitment', 'bilateral_simultaneous')))
        )
        mirror_bm = bool(sp.get('bilateral_mirror_motor', getattr(self, 'bilateral_mirror_motor', False)))

        velocity_mode = str(
            getattr(self, 'isovelocity_velocity_mode', None) or kwargs.get('isovelocity_velocity_mode') or 'angle_vel'
        ).lower()
        vels_native = np.linspace(float(min_vel), float(max_vel), int(num_steps))
        kdot = convert_to_curvature(
            vels_native,
            velocity_mode,
            dclamp_mm=float(dc),
            xsec_width_mm=xw,
            target_muscle_depth_mm=md,
        )
        # Base velocity steps in generated (un-shuffled) order. The non-block path below shuffles
        # once and the block path shuffles each block independently (1.2 step-order randomization).
        vels = np.rad2deg(np.asarray(kdot, dtype=float) * (float(dc) / 1000.0))
        seq_idx = np.arange(vels.size, dtype=int)
        seq_frac = float(sp.get('bilateral_sequential_left_frac', getattr(self, 'bilateral_sequential_left_frac', 0.5)))

        if block_seq:
            if mirror_bm:
                raise ValueError(
                    "bilateral_mirror_motor cannot be used with block_sequence; set direction per block instead."
                )
            left_v = float(sp.get('left_stim_voltage', getattr(self, 'left_stim_voltage', 5.0)))
            right_v = float(sp.get('right_stim_voltage', getattr(self, 'right_stim_voltage', 5.0)))
            reset_ramp = float(
                sp.get('block_reset_ramp_duration_s', getattr(self, 'block_reset_ramp_duration_s', 2.0))
            )
            if not np.isfinite(reset_ramp) or reset_ramp < 0:
                raise ValueError(
                    f"block_reset_ramp_duration_s must be finite and >= 0; got {reset_ramp!r}."
                )
            self._validate_block_sequence_voltages(block_seq, left_v, right_v)
            dev = sp.get('device_name', None) or getattr(self, 'device_name', None)
            all_out = []
            last_deg = 0.0
            global_step = 0
            step_order_blocks = []
            vels_mag_base = np.abs(vels)
            for bi, block in enumerate(block_seq):
                last_deg = self._run_neutral_reset_segment(last_deg, reset_ramp, dev)
                dir_idx = self._lateral_index_for_block_direction(block['direction'])
                ms = self.motor_command_sign_for_bend_toward_index(dir_idx)
                theta0_block = abs(theta0_raw) * ms
                order = self._shuffled_step_order(vels_mag_base.size, randomize, random_seed, block_index=bi)
                step_order_blocks.append([int(x) for x in order])
                vels_mag = vels_mag_base[order]
                block_out = self._run_isovelocity_steps(
                    theta0_block,
                    vels_mag,
                    pre_hold_s=float(sp['pre_hold_s']),
                    iso_duration_s=float(sp['iso_duration_s']),
                    post_baseline_s=float(sp.get('post_baseline_s', 1.0) or 0.0),
                    stim_onset_s=sp.get('stim_onset_s'),
                    settle_before_stim_s=float(sp['settle_before_stim_s']),
                    pre_iso_stim_duration_s=float(sp.get('pre_iso_stim_duration_s', 0.0)),
                    stim_duration_s=sp.get('stim_duration_s'),
                    inter_step_interval_s=float(sp.get('inter_step_interval_s', 0.0) or 0.0),
                    reset_between_steps=bool(sp.get('reset_between_steps', False)),
                    is_stim=bool(sp.get('is_stim', False)),
                    stim_pulse_rate=sp.get('stim_pulse_rate', None),
                    stim_voltage=float(sp.get('stim_voltage', 5.0)),
                    device_name=dev,
                    block_direction=block['direction'],
                    stim_sides=block['stim_sides'],
                    left_stim_voltage=left_v,
                    right_stim_voltage=right_v,
                    block_index=bi,
                    bilateral_mirror_motor=False,
                )
                for j, e in enumerate(block_out):
                    e['starting_strain'] = float(starting_strain)
                    e['starting_strain_mode'] = starting_strain_mode
                    e['curvature_1_per_m'] = k0
                    e['sequence_index'] = int(order[j])
                    e['step_index'] = global_step
                    e['trial_index'] = global_step
                    e['cycle_index'] = global_step
                    global_step += 1
                last_deg = float(block_out[-1]['angle_cmd'][-1])
                all_out.extend(block_out)
            uidx_meta = None
            self.h5_protocol_metadata.update({
                'recruitment': rec_iso,
                'bilateral_mirror_motor': False,
                'bilateral_sequential_left_frac': seq_frac,
                'block_sequence': block_seq,
                'randomize_step_order': bool(randomize),
                'step_order': step_order_blocks,
                'left_stim_voltage': left_v,
                'right_stim_voltage': right_v,
                'block_reset_ramp_duration_s': reset_ramp,
                'specimen_side_index_left': int(self.specimen_side_index_left),
                'specimen_side_index_right': int(self.specimen_side_index_right),
                'specimen_lateral_index_on_positive_motor_side': int(
                    getattr(self, 'specimen_lateral_index_on_positive_motor_side', -1)
                ),
                'motor_positive_bend_toward_lateral_index': int(self.motor_positive_bend_lateral_index()),
                'unilateral_posture_lateral_index': uidx_meta,
                'n_trials': int(len(all_out)),
                'guard_triggered': bool(any(e.get('guard_triggered') for e in all_out)),
                'guard_angle_deg': next(
                    (float(e['guard_angle_deg']) for e in all_out if e.get('guard_triggered')),
                    float('nan'),
                ),
            })
            self.isovelocity_results = all_out
            self.test_type = 'isovelocity'
            self.trial_records = list(all_out)
            # Return to resting length: per-block resets fire only BEFORE each block, so drive the
            # motor to angle = 0° once after the final block before completing (1.7).
            self.return_to_resting_length(device_name=dev)
            return all_out

        # Single (non-block) path: shuffle the velocity-step order once if requested (1.2).
        order = self._shuffled_step_order(vels.size, randomize, random_seed)
        seq_idx = order
        vels = vels[order]
        theta0 = theta0_raw
        rec_iso = self._recruitment_with_bilateral_mirror_motor(rec_iso, mirror_bm)
        uidx0 = self.recruitment_unilateral_lateral_index(rec_iso)
        if uidx0 is not None:
            theta0 = theta0 * self.motor_command_sign_for_bend_toward_index(uidx0)
        out = self._run_isovelocity_steps(
            theta0,
            vels,
            pre_hold_s=float(sp['pre_hold_s']),
            iso_duration_s=float(sp['iso_duration_s']),
            post_baseline_s=float(sp.get('post_baseline_s', 1.0) or 0.0),
            stim_onset_s=sp.get('stim_onset_s'),
            settle_before_stim_s=float(sp['settle_before_stim_s']),
            pre_iso_stim_duration_s=float(sp.get('pre_iso_stim_duration_s', 0.0)),
            stim_duration_s=sp.get('stim_duration_s'),
            inter_step_interval_s=float(sp.get('inter_step_interval_s', 0.0) or 0.0),
            reset_between_steps=bool(sp.get('reset_between_steps', False)),
            is_stim=bool(sp.get('is_stim', False)),
            stim_pulse_rate=sp.get('stim_pulse_rate', None),
            stim_voltage=float(sp.get('stim_voltage', 5.0)),
            device_name=sp.get('device_name', None),
            recruitment=rec_iso,
            bilateral_mirror_motor=mirror_bm,
            bilateral_sequential_left_frac=seq_frac,
        )
        for e in out:
            e['starting_strain'] = float(starting_strain)
            e['starting_strain_mode'] = starting_strain_mode
            e['curvature_1_per_m'] = k0
        for i, e in enumerate(out):
            e['sequence_index'] = int(seq_idx[i])
        uidx_meta = self.recruitment_unilateral_lateral_index(rec_iso)
        self.h5_protocol_metadata.update({
            'recruitment': rec_iso,
            'bilateral_mirror_motor': mirror_bm,
            'bilateral_sequential_left_frac': seq_frac,
            'randomize_step_order': bool(randomize),
            'step_order': [int(x) for x in seq_idx],
            'specimen_side_index_left': int(self.specimen_side_index_left),
            'specimen_side_index_right': int(self.specimen_side_index_right),
            'specimen_lateral_index_on_positive_motor_side': int(
                getattr(self, 'specimen_lateral_index_on_positive_motor_side', -1)
            ),
            'motor_positive_bend_toward_lateral_index': int(self.motor_positive_bend_lateral_index()),
            'unilateral_posture_lateral_index': uidx_meta,
            'guard_triggered': bool(any(e.get('guard_triggered') for e in out)),
            'guard_angle_deg': next(
                (float(e['guard_angle_deg']) for e in out if e.get('guard_triggered')), float('nan')
            ),
        })
        # Return to resting length: drive the motor to angle = 0° before the trial completes (1.7).
        self.return_to_resting_length(device_name=sp.get('device_name', None))
        return out

    def test_isovelocity(
        self,
        min_vel,
        max_vel,
        starting_strain,
        num_steps,
        *,
        starting_strain_mode='strain',
        iso_duration_s=0.2,
        pre_hold_s=0.3,
        stim_params=None,
        **kwargs,
    ):
        """Backward-compatible alias for :meth:`isovelocity`."""
        return self.isovelocity(
            min_vel,
            max_vel,
            starting_strain,
            num_steps,
            starting_strain_mode=starting_strain_mode,
            iso_duration_s=iso_duration_s,
            pre_hold_s=pre_hold_s,
            stim_params=stim_params,
            **kwargs,
        )

    def run_isometric_series(self, *args, **kwargs):
        """Alias for :meth:`run_force_length_series`."""
        return self.run_force_length_series(*args, **kwargs)

    def _uniform_cycles_from_duration(self, duration, f0, amp_deg):
        """Build period/freq/amp arrays so ``make_cycles_dynamic`` runs uniform f and amplitude for ~``duration`` s."""
        dur = float(duration)
        fq = float(f0)
        amp = float(amp_deg)
        if not np.isfinite(dur) or dur <= 0:
            raise ValueError(
                f"_uniform_cycles_from_duration: duration must be finite and > 0 s; got {duration!r}."
            )
        if not np.isfinite(fq) or fq <= 0:
            raise ValueError(
                f"_uniform_cycles_from_duration: frequency f0 must be finite and > 0 Hz; got {f0!r}. "
                "Zero or NaN Hz yields Inf/NaN periods and breaks the time axis."
            )
        if not np.isfinite(amp):
            raise ValueError(f"_uniform_cycles_from_duration: amp_deg must be finite; got {amp_deg!r}.")
        n = max(1, int(round(dur * fq)))
        p = 1.0 / fq
        period_by_cycle = np.full(n, p, dtype=float)
        freq_by_cycle = np.full(n, fq, dtype=float)
        amp_by_cycle = np.full(n, amp, dtype=float)
        return period_by_cycle, freq_by_cycle, amp_by_cycle

    def make_cycles_frequency_step(self, all_freqs, all_curves, duration, waitbefore,
                                   nominal_frequency=None, nominal_curvature=None):
        """Returns (angle, anglevel, tnorm, freq, t); ``freq`` is NaN outside the motion window."""
        dclamp = getattr(self, 'dclamp', None)
        if dclamp is None:
            raise ValueError("frequency_step requires self.dclamp (set via organize_cycles).")
        if duration is None:
            raise ValueError(
                "frequency_step requires motion duration (seconds). "
                "Set bender.duration or use update_metadata(duration=...)."
            )
        fq_src = nominal_frequency if nominal_frequency is not None else all_freqs
        cq_src = nominal_curvature if nominal_curvature is not None else all_curves
        f0, f1 = _normalize_start_end(fq_src)
        c0, c1 = _normalize_start_end(cq_src)
        if abs(f1 - f0) > 1e-9 or abs(c1 - c0) > 1e-9:
            raise ValueError(
                "frequency_step expects a single Hz and 1/m setpoint (float or equal [start, end])."
            )
        if not np.isfinite(f0) or f0 <= 0:
            raise ValueError(
                "make_cycles_frequency_step needs a positive finite frequency (Hz); "
                f"got f={f0!r} from all_freqs/nominal_frequency {fq_src!r}."
            )
        amp_deg = np.rad2deg(c0 * dclamp / 1000.0)
        period_by_cycle, freq_by_cycle, amp_by_cycle = self._uniform_cycles_from_duration(duration, f0, amp_deg)
        # frequency_step never runs organize_cycles, so all_degs/all_freqs may be unset; read with
        # getattr so saving them for make_cycles_dynamic does not raise AttributeError (3.2B).
        saved_degs, saved_freqs = getattr(self, 'all_degs', None), getattr(self, 'all_freqs', None)
        self.all_degs = np.array([amp_deg, amp_deg], dtype=float)
        self.all_freqs = np.array([f0, f0], dtype=float)
        try:
            angle, anglevel, tnorm, t = self.make_cycles_dynamic(
                period_by_cycle, freq_by_cycle, amp_by_cycle, record_protocol=False)
        finally:
            self.all_degs, self.all_freqs = saved_degs, saved_freqs
        movedur = float(np.sum(period_by_cycle))
        freq = np.full_like(t, np.nan, dtype=float)
        mask = (t >= 0) & (t < movedur)
        freq[mask] = f0
        self._protocol_log('frequency_step', f0, f0, c0, c0, motion_duration_s=float(duration))
        return angle, anglevel, tnorm, freq, t

    def make_cycles_curvature_step(self, all_freqs, all_curves, duration, waitbefore,
                                   nominal_frequency=None, nominal_curvature=None):
        """Returns (angle, anglevel, tnorm, freq, t); ``freq`` is NaN outside the motion window."""
        dclamp = getattr(self, 'dclamp', None)
        if dclamp is None:
            raise ValueError("curvature_step requires self.dclamp (set via organize_cycles).")
        if duration is None:
            raise ValueError(
                "curvature_step requires motion duration (seconds). "
                "Set bender.duration or use update_metadata(duration=...)."
            )
        fq_src = nominal_frequency if nominal_frequency is not None else all_freqs
        cq_src = nominal_curvature if nominal_curvature is not None else all_curves
        f0, f1 = _normalize_start_end(fq_src)
        c0, c1 = _normalize_start_end(cq_src)
        if abs(f1 - f0) > 1e-9 or abs(c1 - c0) > 1e-9:
            raise ValueError(
                "curvature_step expects a single Hz and 1/m setpoint (float or equal [start, end])."
            )
        if not np.isfinite(f0) or f0 <= 0:
            raise ValueError(
                "make_cycles_curvature_step needs a positive finite bending frequency (Hz); "
                f"got f={f0!r} from all_freqs/nominal_frequency {fq_src!r}. "
                "Zero or NaN Hz makes cycle periods invalid and produces an empty or NaN time array."
            )
        amp_deg = np.rad2deg(c0 * dclamp / 1000.0)
        period_by_cycle, freq_by_cycle, amp_by_cycle = self._uniform_cycles_from_duration(duration, f0, amp_deg)
        saved_degs, saved_freqs = self.all_degs, self.all_freqs
        self.all_degs = np.array([amp_deg, amp_deg], dtype=float)
        self.all_freqs = np.array([f0, f0], dtype=float)
        try:
            angle, anglevel, tnorm, t = self.make_cycles_dynamic(
                period_by_cycle, freq_by_cycle, amp_by_cycle, record_protocol=False)
        finally:
            self.all_degs, self.all_freqs = saved_degs, saved_freqs
        movedur = float(np.sum(period_by_cycle))
        freq = np.full_like(t, np.nan, dtype=float)
        mask = (t >= 0) & (t < movedur)
        freq[mask] = f0
        self._protocol_log('curvature_step', f0, f0, c0, c0, motion_duration_s=float(duration))
        return angle, anglevel, tnorm, freq, t

    def make_cycles_frequency_sweep(self, all_freqs, all_curves, amplitude_frequency_exponent, duration, waitbefore,
                                    nominal_frequency=None, nominal_curvature=None):
        """Log-frequency sweep with curvature-based amplitude scaling (exponent matches legacy sweep).

        NOTE: ``amplitude_frequency_exponent`` (α) scales *commanded amplitude with instantaneous
        frequency*: θ ∝ f^α relative to the sweep start frequency. It is unrelated to
        ``self.velocity_exponent``, which only affects :func:`_ramp_progress` when ramp_mode
        is ``'exponential'`` (time-warp of generic ramps, not used in this sweep implementation).

        Returns:
            tuple: (angle, anglevel, tnorm, freq, t) — ``freq`` is NaN outside the motion window.
        """
        daq_ai_hz = self.daq_ai_sample_rate_hz
        dclamp = getattr(self, 'dclamp', None)
        if dclamp is None:
            raise ValueError("make_cycles_frequency_sweep requires self.dclamp (set via organize_cycles).")
        if duration is None:
            raise ValueError(
                "frequency_sweep requires motion duration (seconds). "
                "Set bender.duration or use update_metadata(duration=...)."
            )
        if amplitude_frequency_exponent is None:
            raise ValueError(
                "frequency_sweep requires amplitude_frequency_exponent. "
                "Set bender.amplitude_frequency_exponent or pass it via update_metadata(...)."
            )
        fq_src = nominal_frequency if nominal_frequency is not None else all_freqs
        cq_src = nominal_curvature if nominal_curvature is not None else all_curves
        f0, f1 = _normalize_start_end(fq_src)
        c0, c1 = _normalize_start_end(cq_src)

        total_duration = duration + waitbefore + self.waitafter
        t = np.arange(0, total_duration, 1.0 / daq_ai_hz)
        t -= waitbefore
        movedur = float(duration)
        startfreq, endfreq = f0, f1
        lnk = 1.0 / movedur * (np.log(endfreq) - np.log(startfreq))
        freq = startfreq * np.exp(t * lnk)
        tnorm = 2 * np.pi * startfreq * (np.exp(t * lnk) - 1) / lnk
        tnorm[t < 0] = -1
        tnorm[t > movedur] = np.ceil(np.max(tnorm))
        A0 = startfreq ** amplitude_frequency_exponent
        amplitude = np.rad2deg(c0 * dclamp / 1000.0)
        angle = amplitude / A0 * np.power(freq, amplitude_frequency_exponent) * np.sin(tnorm)
        anglevel = amplitude / A0 * np.exp(amplitude_frequency_exponent * t * lnk) * lnk * (
            amplitude_frequency_exponent * np.sin(tnorm) + 2 * np.pi / lnk * freq * np.cos(tnorm))
        freq[t < 0] = np.nan
        freq[t > movedur] = np.nan
        angle[t < 0] = 0
        angle[t > movedur] = 0
        anglevel[t < 0] = 0
        anglevel[t > movedur] = 0
        isramp = (t >= movedur) & (t < movedur + 0.5)
        k_index = np.argmax(t >= (waitbefore + movedur))
        pend = angle[k_index]
        velend = (0 - pend) / 0.5
        np.place(anglevel, isramp, velend)
        np.place(angle, isramp, pend + (t[isramp] - t[k_index]) * velend)

        self._protocol_log(
            'frequency_sweep', f0, f1, c0, c1,
            motion_duration_s=movedur,
            amplitude_frequency_exponent=float(amplitude_frequency_exponent),
        )

        self.t = t
        self.tnorm = tnorm
        # Order is (angle, anglevel, tnorm, freq, t) — freq has NaNs during waitbefore/waitafter.
        return angle, anglevel, tnorm, freq, t

    def make_cycles_dynamic(self, period_by_cycle, freq_by_cycle, amp_by_cycle, record_protocol=True):
        """
        Generates signals for dynamic tests with a custom sequence of frequencies and amplitudes.

        Includes amplitude ramps during step changes and start/end ramps.

        Args:
            period_by_cycle (list/array): Duration of each cycle in seconds.
            freq_by_cycle (list/array): Frequency for each cycle (Hz).
            amp_by_cycle (list/array): Amplitude for each cycle (degrees).
            record_protocol: If False, skip MasterLogger protocol entry (used when a wrapper logs).

        Returns:
            tuple: (angle, anglevel, tnorm, t) NumPy arrays.
        """

        # 1. KILL THE OLD CLOCK (This is the most common 'hang' source)
        self.t = None 
        self.tnorm = None
        self.angle = None
        
        # Access parameters from self (ensure they are set by the __init__ or another method)
        waitbefore = self.waitbefore
        waitafter = self.waitafter
        rampdur = self.rampdur
        amp_step_vel = self.amp_step_vel
        daq_ai_hz = self.daq_ai_sample_rate_hz
        all_degs = self.all_degs # Used for start/end ramps
        all_freqs = self.all_freqs # Used for start/end ramps

        if not np.isfinite(amp_step_vel) or float(amp_step_vel) == 0.0:
            raise ValueError(
                "make_cycles_dynamic: amp_step_vel must be finite and non-zero (deg/s scale for inter-cycle ramps)."
            )
        period_by_cycle = np.asarray(period_by_cycle, dtype=float).reshape(-1)
        freq_by_cycle = np.asarray(freq_by_cycle, dtype=float).reshape(-1)
        amp_by_cycle = np.asarray(amp_by_cycle, dtype=float).reshape(-1)
        if period_by_cycle.size != freq_by_cycle.size or period_by_cycle.size != amp_by_cycle.size:
            raise ValueError(
                "make_cycles_dynamic: period_by_cycle, freq_by_cycle, and amp_by_cycle must have the same length."
            )
        if not np.all(np.isfinite(period_by_cycle)) or not np.all(np.isfinite(freq_by_cycle)) or not np.all(np.isfinite(amp_by_cycle)):
            raise ValueError("make_cycles_dynamic: period_by_cycle, freq_by_cycle, and amp_by_cycle must be finite (no NaN/Inf).")
        if np.any(freq_by_cycle <= 0):
            raise ValueError("make_cycles_dynamic: every entry in freq_by_cycle must be > 0 Hz.")
        if np.any(period_by_cycle <= 0):
            raise ValueError("make_cycles_dynamic: every entry in period_by_cycle must be > 0 s.")
        all_freqs = np.asarray(all_freqs, dtype=float).reshape(-1)
        all_degs = np.asarray(all_degs, dtype=float).reshape(-1)
        if all_freqs.size < 1 or all_degs.size < 1:
            raise ValueError(
                "make_cycles_dynamic: self.all_freqs and self.all_degs must be non-empty (set via organize_cycles)."
            )
        if not np.isfinite(all_freqs[0]) or all_freqs[0] <= 0 or not np.isfinite(all_freqs[-1]) or all_freqs[-1] <= 0:
            raise ValueError(
                "make_cycles_dynamic: self.all_freqs[0] and [-1] must be finite and > 0 Hz (start/end ramp timing)."
            )

        # Calculate timings and durations
        movedur = float(np.sum(period_by_cycle))
        if not np.isfinite(movedur) or movedur <= 0:
            raise ValueError(
                f"make_cycles_dynamic: total motion duration from period_by_cycle is invalid ({movedur!r})."
            )
        totaldur = float(waitbefore + movedur + waitafter)
        if not np.isfinite(totaldur) or totaldur <= 0:
            raise ValueError(
                f"make_cycles_dynamic: total timeline waitbefore+movedur+waitafter is invalid ({totaldur!r})."
            )
        if not np.isfinite(daq_ai_hz) or daq_ai_hz <= 0:
            raise ValueError(f"make_cycles_dynamic: daq_ai_sample_rate_hz must be finite and > 0; got {daq_ai_hz!r}.")
        t = np.arange(0, totaldur, 1.0 / daq_ai_hz)
        if t.size < 2:
            raise ValueError(
                "make_cycles_dynamic: generated fewer than 2 time samples; increase duration or wait intervals, "
                "or check daq_ai_sample_rate_hz."
            )
        if not np.all(np.isfinite(t)):
            raise ValueError("make_cycles_dynamic: internal error — time axis t contains non-finite values.")
        t = np.asarray(t, dtype=float)
        t -= waitbefore # Shift time so the movement starts at t=0

        # Generate the base signals
        freq = np.zeros_like(t)
        amp = np.zeros_like(t)
        tnorm = np.zeros_like(t)

        cyclestart = np.cumsum(period_by_cycle)
        cyclestart = np.insert(cyclestart, 0, 0)

        for c, (cycstart1, f1, a1) in enumerate(zip(cyclestart, freq_by_cycle, amp_by_cycle)):
            cycend1 = cycstart1 + 1/f1
            iscyc = (t >= cycstart1) & (t < cycend1)
            freq[iscyc] = f1
            amp[iscyc] = a1
            # Use boolean assignment instead of np.place
            tnorm[iscyc] = (t[iscyc] - cycstart1) * f1 + c

        # Add linear ramps between cycles 
        for c, (cycstart1, a1, a2) in enumerate(zip(cyclestart[1:], amp_by_cycle[:-1], amp_by_cycle[1:])):
            amp_step_dur2 = (a2 - a1) / amp_step_vel / 2
            isstep = (t >= cycstart1 - amp_step_dur2) & (t < cycstart1 + amp_step_dur2)
            
            if np.any(isstep):
                amp_ramp = np.linspace(a1, a2, np.sum(isstep))
                amp[isstep] = amp_ramp # Use boolean assignment

        angle = amp * np.sin(2 * np.pi * tnorm)

        # Ensure signal is zero during wait periods
        angle[t < 0] = 0
        angle[t > movedur] = 0

        # Ramp to the start and end amplitudes (Original logic using boolean assignment)
        rampvel1 = all_degs[0] / rampdur
        tendramp1 = 0.25 / all_freqs[0]
        tstartramp1 = tendramp1 - rampdur
        rampvel2 = all_degs[-1] / rampdur
        tstartramp2 = movedur - 0.25 / all_freqs[-1]
        tendramp2 = tstartramp2 + rampdur

        if tstartramp1 > 0:
            pass # actual movement is slower than the ramp, so we won't bother adding the ramp
        else:
            # Use boolean assignment
            mask1 = (t >= tstartramp1) & (t < tendramp1)
            mask2 = (t >= tstartramp2) & (t < tendramp2)
            angle[mask1] = (t[mask1] - tstartramp1) * rampvel1
            angle[mask2] = (t[mask2] - tstartramp2 - rampdur) * rampvel2

        # Calculate angular velocity
        anglevel = np.zeros_like(angle)

        # DOUBLE CHECK: Avoid calculating velocity for the very first and last point
        anglevel[1:-1] = (angle[2:] - angle[:-2]) * (daq_ai_hz / 2.0)

        self.t = t
        self.tnorm = tnorm

        if record_protocol:
            lf0, lf1 = float(freq_by_cycle[0]), float(freq_by_cycle[-1])
            dc = getattr(self, 'dclamp', None)
            dyn_extra = {'dynamic_movedur_s': float(movedur)}
            if dc is not None:
                lc0 = float(np.deg2rad(amp_by_cycle[0]) * 1000.0 / dc)
                lc1 = float(np.deg2rad(amp_by_cycle[-1]) * 1000.0 / dc)
                self._protocol_log('dynamic', lf0, lf1, lc0, lc1, **dyn_extra)
            else:
                self.master_logger.record(
                    protocol='dynamic',
                    frequency_hz=_scalar_or_pair(lf0, lf1),
                    **dyn_extra,
                )

        return angle, anglevel, tnorm, t

    def make_cycles_step_change(self, frequencies, curves, cycles_per_step, dclamp=None,
                                amp_step_vel=None, record_protocol=True):
        """
        Step-change protocol (``step_change.ipynb``): each block has one (frequency, curvature)
        pair repeated ``cycles_per_step[i]`` times. Builds per-cycle arrays and delegates to
        :meth:`make_cycles_dynamic`.

        Arrays ``frequencies``, ``curves``, and ``cycles_per_step`` must have the same length
        (one row per step in the schedule).

        Also assigns ``self.freq_by_cycle``, ``self.amp_by_cycle``, and ``self.period_by_cycle``
        to the expanded per-cycle vectors for stimulation bookkeeping.
        """
        dc = float(dclamp) if dclamp is not None else getattr(self, 'dclamp', None)
        if dc is None:
            raise ValueError("make_cycles_step_change requires dclamp (argument or self.dclamp).")
        # Wrap bare scalars in a 1-element array before the length check so a single-step schedule
        # (one frequency, one curvature, one cycle count) always has matching lengths (3.4B).
        freqs = np.atleast_1d(np.asarray(frequencies, dtype=float)).reshape(-1)
        curv = np.atleast_1d(np.asarray(curves, dtype=float)).reshape(-1)
        cps = np.atleast_1d(np.asarray(cycles_per_step, dtype=int)).reshape(-1)
        if not (len(freqs) == len(curv) == len(cps)):
            raise ValueError("frequencies, curves, and cycles_per_step must have the same length.")
        allfreqs = np.concatenate([np.full((int(c),), f, dtype=float) for c, f in zip(cps, freqs)])
        allcurves = np.concatenate([np.full((int(c),), k, dtype=float) for c, k in zip(cps, curv)])
        period_by_cycle = 1.0 / allfreqs
        freq_by_cycle = allfreqs
        amp_by_cycle = np.rad2deg(allcurves * (dc / 1000.0))

        self.period_by_cycle = period_by_cycle
        self.freq_by_cycle = freq_by_cycle
        self.amp_by_cycle = amp_by_cycle

        # step_change never runs organize_cycles, so all_degs/all_freqs may be unset; read with
        # getattr so saving them for make_cycles_dynamic does not raise AttributeError (3.4C/3.2B).
        saved_degs = getattr(self, 'all_degs', None)
        saved_freqs = getattr(self, 'all_freqs', None)
        saved_v = getattr(self, 'amp_step_vel', None)
        try:
            self.all_degs = np.array([amp_by_cycle[0], amp_by_cycle[-1]], dtype=float)
            self.all_freqs = np.array([freq_by_cycle[0], freq_by_cycle[-1]], dtype=float)
            if amp_step_vel is not None:
                self.amp_step_vel = float(amp_step_vel)
            angle, anglevel, tnorm, t = self.make_cycles_dynamic(
                period_by_cycle, freq_by_cycle, amp_by_cycle, record_protocol=False)
        finally:
            self.all_degs = saved_degs
            self.all_freqs = saved_freqs
            self.amp_step_vel = saved_v

        movedur = float(np.sum(period_by_cycle))
        if record_protocol:
            f_lo, f_hi = float(np.min(allfreqs)), float(np.max(allfreqs))
            c_lo, c_hi = float(np.min(allcurves)), float(np.max(allcurves))
            # Log the per-step bend amplitude(s) in curvature units (1/m): scalar for a single-step
            # schedule, else the list of per-step curvatures (3.4A).
            step_amp = float(curv[0]) if curv.size == 1 else [float(x) for x in curv]
            self._protocol_log(
                'step_change', f_lo, f_hi, c_lo, c_hi,
                motion_duration_s=movedur,
                step_change_blocks=int(len(freqs)),
                step_change_total_cycles=int(len(allfreqs)),
                dclamp_mm=float(dc),
                step_amplitude_inv_m=step_amp,
            )
        return angle, anglevel, tnorm, t, movedur

    def make_stimuli(self, is_stim=None, phase_by_cycle=None, stim_pulse_rate=None, 
                      prestim_time=None, poststim_time=None, prepoststim_dur=None, 
                      prepoststim_sep=None, stimburstdur=None, duty_by_cycle=None, 
                      freq_by_cycle=None, movedur=None, t_basis=None, tnorm_basis=None):
        
        # 1. Reset ONLY the labels/lists (These are the 'ghost' sources)
        self.Lonoff = []
        self.Ronoff = []
        Lonoff = []
        Ronoff = []

        # 1. Quick-fill defaults using getattr
        is_stim = is_stim if is_stim is not None else getattr(self, 'is_stim', False)
        
        # EARLY EXIT: If no stim, don't even create the empty arrays yet
   # EARLY EXIT: If no stim, create empty arrays that MATCH the time length
        if not is_stim:
            t = t_basis if t_basis is not None else self.t
            self.S1stimcmd = np.zeros_like(t)
            self.S2stimcmd = np.zeros_like(t)
            self.Lonoff, self.Ronoff = [], []
            return self.S1stimcmd, self.S2stimcmd
        
        # 2. Grab the 'Self' versions if arguments are None
        phase_by_cycle = phase_by_cycle if phase_by_cycle is not None else getattr(self, 'phase_by_cycle', [])
        duty_by_cycle  = duty_by_cycle  if duty_by_cycle  is not None else getattr(self, 'duty_by_cycle', [])
        freq_by_cycle  = freq_by_cycle  if freq_by_cycle  is not None else getattr(self, 'freq_by_cycle', [])
        stimburstdur   = stimburstdur   if stimburstdur   is not None else getattr(self, 'stimburstdur', [])
        stim_pulse_rate = stim_pulse_rate if stim_pulse_rate is not None else getattr(self, 'stim_pulse_rate', 75)
        prepoststim_dur = prepoststim_dur if prepoststim_dur is not None else getattr(self, 'prepoststim_dur', 0.06)
        prestim_time    = prestim_time    if prestim_time    is not None else getattr(self, 'prestim_time', -2.0)
        prepoststim_sep = prepoststim_sep if prepoststim_sep is not None else getattr(self, 'prepoststim_sep', 1.0)
        poststim_time   = poststim_time   if poststim_time   is not None else getattr(self, 'poststim_time', 2.0)
        movedur         = movedur         if movedur         is not None else getattr(self, 'movedur', 0.0)


        # 3. HEAVY LIFTING
        # Use the passed-in basis if available, otherwise fallback to self
        t = t_basis if t_basis is not None else self.t
        tnorm = tnorm_basis if tnorm_basis is not None else self.tnorm

        S1stimcmd = np.zeros_like(t)
        S2stimcmd = np.zeros_like(t)
        Lonoff, Ronoff = [], []

        # Route left-phase and right-phase bursts to the AO channel that serves each
        # physical side, as declared in the hardware config (S1side / S2side).
        s1_is_left = str(getattr(self, 'S1side', 'left')).strip().lower() == 'left'

        def _assign_left(mask, pw):
            if s1_is_left:
                S1stimcmd[mask] = pw[mask]
            else:
                S2stimcmd[mask] = pw[mask]

        def _assign_right(mask, pw):
            if s1_is_left:
                S2stimcmd[mask] = pw[mask]
            else:
                S1stimcmd[mask] = pw[mask]

        # Calculate the pulse wave once
        # Using 5.0 for the stimulator 'On' voltage; pulse high time from pulse_width_ms.
        _pw_frac = self._pulse_high_fraction(stim_pulse_rate)
        pulse_wave = (np.mod(t * stim_pulse_rate, 1) < _pw_frac).astype(float) * 5.0
        bendphase = tnorm - 0.25

        # 4. Optimized Cycle Loop
        for c, (dur1, duty1, f1, p1) in enumerate(zip(stimburstdur, duty_by_cycle, freq_by_cycle, phase_by_cycle)):
            if dur1 <= 0: continue

            # Left-phase burst
            t_left = c + p1
            m_left = (bendphase >= t_left) & (bendphase < t_left + duty1)
            if np.any(m_left):
                _assign_left(m_left, pulse_wave)
                t_sub = t[m_left]
                Lonoff.append([t_sub[0], t_sub[-1]])

            # Right-phase burst (half-cycle offset)
            t_right = c + 0.5 + p1
            m_right = (bendphase >= t_right) & (bendphase < t_right + duty1)
            if np.any(m_right):
                _assign_right(m_right, pulse_wave)
                t_sub2 = t[m_right]
                Ronoff.append([t_sub2[0], t_sub2[-1]])

        # 5. Pre-stimulation (bilateral sequential bursts before bending at t=0)
        # S1 channel fires at prestim_time; S2 channel fires prepoststim_sep later; both last prepoststim_dur.
        pre_s1_start = prestim_time
        pre_s1_end = prestim_time + prepoststim_dur
        m_pre_s1 = (t >= pre_s1_start) & (t < pre_s1_end)
        if np.any(m_pre_s1):
            S1stimcmd[m_pre_s1] = pulse_wave[m_pre_s1]
            if s1_is_left:
                Lonoff.append([pre_s1_start, pre_s1_end])
            else:
                Ronoff.append([pre_s1_start, pre_s1_end])

        pre_s2_start = prestim_time + prepoststim_sep
        pre_s2_end = pre_s2_start + prepoststim_dur
        m_pre_s2 = (t >= pre_s2_start) & (t < pre_s2_end)
        if np.any(m_pre_s2):
            S2stimcmd[m_pre_s2] = pulse_wave[m_pre_s2]
            if s1_is_left:
                Ronoff.append([pre_s2_start, pre_s2_end])
            else:
                Lonoff.append([pre_s2_start, pre_s2_end])

        # 5b. Post-stimulation (bilateral sequential bursts after motion ends)
        # S1 channel fires at movedur + poststim_time; S2 channel fires prepoststim_sep later; both last prepoststim_dur.
        post_s1_start = movedur + poststim_time
        post_s1_end = post_s1_start + prepoststim_dur
        m_post_s1 = (t >= post_s1_start) & (t < post_s1_end)
        if np.any(m_post_s1):
            S1stimcmd[m_post_s1] = pulse_wave[m_post_s1]
            if s1_is_left:
                Lonoff.append([post_s1_start, post_s1_end])
            else:
                Ronoff.append([post_s1_start, post_s1_end])

        post_s2_start = movedur + poststim_time + prepoststim_sep
        post_s2_end = post_s2_start + prepoststim_dur
        m_post_s2 = (t >= post_s2_start) & (t < post_s2_end)
        if np.any(m_post_s2):
            S2stimcmd[m_post_s2] = pulse_wave[m_post_s2]
            if s1_is_left:
                Ronoff.append([post_s2_start, post_s2_end])
            else:
                Lonoff.append([post_s2_start, post_s2_end])

        # 6. Save and Return
        self.Lonoff = Lonoff
        self.Ronoff = Ronoff
        self.S1stimcmd = S1stimcmd
        self.S2stimcmd = S2stimcmd
        
        return S1stimcmd, S2stimcmd

    def set_physics(self, clamp_offset, front_h, front_w, back_h, back_w, spec_length, mode="lateral"):
        """
        Calculates the Moment of Inertia (MOI) for the system and specimen.
        """
        # --- THE TRANSLATOR (Fish -> Machine) ---
        # Fixed: Using the correct argument names (front_w, front_h, spec_length)
        if mode == "lateral":
            w_math, h_math, d_math = front_w, front_h, spec_length
        elif mode == "dorsoventral":
            w_math, h_math, d_math = front_h, front_w, spec_length
        elif mode == "torsional":
            w_math, h_math, d_math = front_w, spec_length, front_h 
        else:
            raise ValueError("Invalid mode. Choose 'lateral', 'dorsoventral', or 'torsional'.")

        # --- CONSTANTS ---
        RHO_PLA, RHO_STEEL, RHO_NEO, RHO_SPECIMEN = 0.001116, 0.008, 0.0075, 0.001 
        CLAMP_DIM = np.array([50.0, 100.0, 20.0]) # W, H, D
        M_BOLT = 12.0 # Total hardware per clamp

        # --- 1. SHAFT BASELINE (12" Steel) ---
        m_shaft = (np.pi * (9.525/2)**2 * 304.8) * RHO_STEEL
        i_shaft = 0.5 * m_shaft * (9.525/2)**2

        # --- 2. ROTATING CLAMP PAIR ---
        r_cm = clamp_offset + (CLAMP_DIM[2] / 2)
        m_unit = (np.prod(CLAMP_DIM) * RHO_PLA) + (np.pi*5**2*3*RHO_NEO) + M_BOLT
        i_rotating_clamps = 2 * (m_unit * r_cm**2)

        # --- 3. SPECIMEN MOI (The Fish) ---
        # Note: Using d_math, h_math, w_math from the translator above
        i_spec, m_spec = self.calculate_moi_specimen(
            rho_eff=RHO_SPECIMEN, 
            obj_depth_length=d_math,
            front_h_semi=h_math/2, 
            back_h_semi=h_math/2,
            front_w_semi=w_math/2, 
            back_w_semi=w_math/2,
            num_samples=50,      
            axis_offset_x=0.0,   
            axis_offset_z=0.0    
        )

        # --- 4. FINALIZE ---
        self.i_total_system = i_shaft + i_rotating_clamps + i_spec
        self.total_mass = m_shaft + (m_unit * 2) + m_spec
        
        # --- 5. REPORT ---
        print("\n" + "="*50)
        print(f"{'PHYSICS CONFIGURATION REPORT':^50}")
        print("="*50)
        print(f"{'Mode':<25} | {mode:<21}")
        print(f"{'Total Rotating MOI':<25} | {self.i_total_system:<12.2f} | g*mm²")
        print(f"{'Lever Arm (r)':<25} | {r_cm:<12.2f} | mm")
        print(f"{'Specimen Mass':<25} | {m_spec:<12.2f} | g")
        print("-" * 50)
        print(f"{'TOTAL SYSTEM MASS':<25} | {self.total_mass:<12.2f} | g")
        print("="*50 + "\n")

    def get_corrected_torque(self, raw_data_dict):
        """
        Subtracts the inertial load from the specific sensor channel 
        that matches the motor's rotation.
        """
        # Grab the raw data from the channel the user specified
        axis = str(getattr(self, 'primary_bending_axis', getattr(self, 'bending_axis_sensor', 'zTorque'))).strip()
        axis_l = axis.lower()
        key = {'x': 'xTorque', 'xtorque': 'xTorque', 'y': 'yTorque', 'ytorque': 'yTorque', 'z': 'zTorque', 'ztorque': 'zTorque'}.get(axis_l, axis)
        raw_torque = raw_data_dict[key]
        
        # alpha = angular acceleration
        alpha = np.gradient(self.anglevel) * self.daq_ai_sample_rate_hz
        
        # The Correction
        true_torque = raw_torque - (self.i_total_system * alpha)
        return true_torque

    def summary(self):
        # 1. Build the list of lines
        lines = [
            "="*50,
            f"{'BENDER SYSTEM SUMMARY':^50}",
            "="*50,
            f"Config:      {self.config_name}.py",
            f"Device:      {self.device_name}",
            f"Motor Port:  {self.motor_port}",
            f"Direction:   POSITIVE = {self.positive_motor_direction.upper()}",
            "-" * 50,
            f"Cal File:    {self.forcetorque_calibration_file}",
            f"Sample Rate: {self.daq_ai_sample_rate_hz} Hz",
            f"Ramp:        {self.rampdur} s",
            "="*50
        ]
        
        # 2. Join them into one big block of text
        report_text = "\n".join(lines)
        
        # 3. Print it now AND send it back to the notebook
        print(report_text)
        return report_text
    
    def update_metadata(self, **kwargs):
        """Saves any passed-in variables directly to the bender object."""
        # Optional high-level payload: protocol_params dict.
        # This lets notebook/GUI pass one object and keep UI fields simple.
        protocol_params = kwargs.pop('protocol_params', None)
        if isinstance(protocol_params, dict):
            pp = dict(protocol_params)
            # Common aliases for protocol name
            if 'test_type' in pp and 'test_type' not in kwargs:
                kwargs['test_type'] = pp.pop('test_type')
            if 'type' in pp and 'test_type' not in kwargs:
                kwargs['test_type'] = pp.pop('type')

            tt = str(kwargs.get('test_type', getattr(self, 'test_type', ''))).lower()
            if tt == 'isometric':
                alias_map = {
                    'initial': 'isometric_initial',
                    'final': 'isometric_final',
                    'num_steps': 'isometric_num_steps',
                    'mode': 'isometric_mode',
                    'randomize': 'isometric_randomize',
                    'random_seed': 'isometric_random_seed',
                    'time_between_steps_s': 'isometric_inter_step_interval_s',
                    'inter_step_interval_s': 'isometric_inter_step_interval_s',
                    'stim_params': 'isometric_stim_params',
                    'recruitment': 'recruitment',
                    'lateral_mode': 'lateral_mode',
                    'bilateral_mirror_motor': 'bilateral_mirror_motor',
                    'bilateral_sequential_left_frac': 'bilateral_sequential_left_frac',
                }
            elif tt == 'isovelocity':
                alias_map = {
                    'min_vel': 'isovelocity_min_vel',
                    'max_vel': 'isovelocity_max_vel',
                    'starting_strain': 'isovelocity_starting_strain',
                    'starting_strain_mode': 'isovelocity_starting_strain_mode',
                    'num_steps': 'isovelocity_num_steps',
                    'randomize': 'isovelocity_randomize',
                    'random_seed': 'isovelocity_random_seed',
                    'iso_duration_s': 'isovelocity_iso_duration_s',
                    'pre_hold_s': 'isovelocity_pre_hold_s',
                    'stim_params': 'isovelocity_stim_params',
                    'recruitment': 'recruitment',
                    'lateral_mode': 'lateral_mode',
                    'bilateral_mirror_motor': 'bilateral_mirror_motor',
                    'bilateral_sequential_left_frac': 'bilateral_sequential_left_frac',
                }
            elif tt == 'calibration':
                alias_map = {
                    'base_test_type': 'calibration_base_test_type',
                }
            else:
                alias_map = {}

            for k, v in pp.items():
                dest = alias_map.get(k, k)
                if dest not in kwargs:
                    kwargs[dest] = v

        for key, value in kwargs.items():
            setattr(self, key, value)
            print(f"  Stored: {key} = {value}")

        # Apply flag/alias normalization so notebooks and GUI only provide simple inputs.
        if 'use_sono' in kwargs:
            self.use_sono = bool(self.use_sono)
            self._apply_use_sono()
            print(f"  Auto-applied use_sono={self.use_sono}; channels={self.input_channel_names}")

        if 'primary_bending_axis' in kwargs or 'bending_axis_sensor' in kwargs:
            norm_axis = self._normalize_primary_bending_axis()
            print(f"  Auto-normalized primary_bending_axis -> {norm_axis}")

        if any(k in kwargs for k in (
            'isometric_initial', 'isometric_final', 'isometric_num_steps', 'initial', 'final', 'num_steps',
            'isovelocity_min_vel', 'isovelocity_max_vel', 'isovelocity_starting_strain', 'min_vel', 'max_vel', 'starting_strain'
        )):
            self._normalize_dispatch_aliases()

        # Optional auto-build of specimen inertia model from profiled stations (preferred).
        profile_keys = {
            'specimen_profile_stations', 'specimen_profile_length_mm', 'specimen_profile_density_g_per_mm3',
        }
        if profile_keys.intersection(kwargs.keys()):
            auto_flag = bool(getattr(self, 'use_frustum_inertial_model', True))
            if auto_flag:
                self.set_profiled_specimen_inertial_model(
                    stations=getattr(self, 'specimen_profile_stations'),
                    length_mm=getattr(self, 'specimen_profile_length_mm'),
                    density_g_per_mm3=getattr(self, 'specimen_profile_density_g_per_mm3'),
                    clamp_offset_mm=float(getattr(self, 'specimen_profile_clamp_offset_mm', 20.0)),
                    num_samples=int(getattr(self, 'specimen_profile_num_samples', 120)),
                )
                print("  Auto-built profiled specimen inertial model from setup inputs.")

        # Backward-compatible frustum auto-build.
        frustum_keys = {
            'frustum_height_mm', 'frustum_width_mm', 'frustum_length_mm',
            'frustum_density_g_per_mm3',
        }
        if frustum_keys.intersection(kwargs.keys()) and not profile_keys.intersection(kwargs.keys()):
            auto_flag = bool(getattr(self, 'use_frustum_inertial_model', True))
            if auto_flag:
                self.set_frustum_inertial_model(
                    height_mm=getattr(self, 'frustum_height_mm'),
                    width_mm=getattr(self, 'frustum_width_mm'),
                    length_mm=getattr(self, 'frustum_length_mm'),
                    density_g_per_mm3=getattr(self, 'frustum_density_g_per_mm3'),
                    tip_scale=float(getattr(self, 'frustum_tip_scale', 0.0)),
                    clamp_offset_mm=float(getattr(self, 'frustum_clamp_offset_mm', 20.0)),
                    num_samples=int(getattr(self, 'frustum_num_samples', 100)),
                )
                print("  Auto-built frustum inertial model from setup inputs.")

        # UX helper: if user provides curve input values + units, auto-populate all_curves.
        # This keeps notebook inputs simple (values + mode) while centralizing conversion logic.
        if 'all_amps' in kwargs and 'curve_input_values' not in kwargs:
            self.curve_input_values = kwargs.get('all_amps')
        if 'all_amps_mode' in kwargs and 'curve_input_mode' not in kwargs:
            self.curve_input_mode = kwargs.get('all_amps_mode')

        if 'all_curves' not in kwargs and (
            'curve_input_values' in kwargs or 'curve_input_mode' in kwargs
            or 'all_amps' in kwargs or 'all_amps_mode' in kwargs
        ):
            vals = getattr(self, 'curve_input_values', None)
            mode = getattr(self, 'curve_input_mode', None)
            if vals is not None and mode is not None:
                conv = self.get_all_amps(vals, mode=mode)
                self.all_curves = np.asarray(conv['curvature_1_per_m'], dtype=float).tolist()
                print(f"  Auto-converted all_amps ({mode}) -> all_curves = {self.all_curves}")

    def get_dispatch_schema(self):
        """Return GUI-friendly schema for run_experiment dispatcher fields."""
        return {
            'test_types': [
                'dynamic', 'frequency_sweep', 'frequency_step', 'curvature_step',
                'step_change', 'isometric', 'isovelocity', 'calibration',
            ],
            'common_optional': [
                'use_sono', 'use_inertial_calibration', 'inertial_calibration_file',
                'post_trial_notes', 'primary_bending_axis', 'protocol_params',
                'all_amps', 'all_amps_mode', 'curve_input_values', 'curve_input_mode',
                'use_frustum_inertial_model', 'use_theoretical_inertial_correction',
                'specimen_profile_stations', 'specimen_profile_length_mm',
                'specimen_profile_density_g_per_mm3', 'specimen_profile_clamp_offset_mm',
                'specimen_profile_num_samples',
                'frustum_height_mm', 'frustum_width_mm', 'frustum_length_mm',
                'frustum_density_g_per_mm3', 'frustum_tip_scale',
                'frustum_clamp_offset_mm', 'frustum_num_samples',
            ],
            'isometric_required': [
                'isometric_initial', 'isometric_final', 'isometric_num_steps',
                'block_sequence', 'left_stim_voltage', 'right_stim_voltage',
                'block_reset_ramp_duration_s',
            ],
            'isometric_optional': [
                'isometric_mode', 'randomize_step_order', 'isometric_random_seed',
                'rest_between_steps_s', 'reset_between_steps',
                'isometric_stim_params',
            ],
            'isovelocity_required': [
                'isovelocity_min_vel', 'isovelocity_max_vel',
                'isovelocity_starting_strain', 'isovelocity_starting_strain_mode',
                'isovelocity_velocity_mode', 'isovelocity_num_steps',
                'block_sequence', 'left_stim_voltage', 'right_stim_voltage',
                'block_reset_ramp_duration_s',
            ],
            'isovelocity_optional': [
                'randomize_step_order',
                'isovelocity_random_seed', 'isovelocity_iso_duration_s',
                'isovelocity_pre_hold_s', 'rest_between_steps_s', 'reset_between_steps',
                'isovelocity_stim_params',
            ],
            'calibration_required': ['calibration_base_test_type'],
            'calibration_optional': ['inertial_calibration_file'],
            'legacy_aliases_mirrored_to_canonical': {
                'isometric': {
                    'initial': 'isometric_initial',
                    'final': 'isometric_final',
                    'num_steps': 'isometric_num_steps',
                    'mode': 'isometric_mode',
                    'randomize': 'isometric_randomize',
                    'random_seed': 'isometric_random_seed',
                },
                'isovelocity': {
                    'min_vel': 'isovelocity_min_vel',
                    'max_vel': 'isovelocity_max_vel',
                    'starting_strain': 'isovelocity_starting_strain',
                    'num_steps': 'isovelocity_num_steps',
                    'randomize': 'isovelocity_randomize',
                    'random_seed': 'isovelocity_random_seed',
                },
            },
        }

    def _announce_pre_run_max_rotation(self, test_type=None):
        """Print peak commanded |angle| (deg) before protocol motor motion (from validation cache when set)."""
        from shared.utilities import compute_max_rotation_deg, format_max_rotation_message

        tt = str(test_type or getattr(self, 'test_type', 'dynamic'))
        deg = getattr(self, 'max_commanded_rotation_deg', None)
        if deg is None:
            try:
                deg = float(compute_max_rotation_deg({'bender': self, 'test_type': tt}))
                self.max_commanded_rotation_deg = deg
            except Exception as exc:
                logging.warning('Pre-run max rotation not computed for %s: %s', tt, exc)
                return
        print(format_max_rotation_message(float(deg)))

    def validate_dispatch_setup(self, test_type=None):
        """
        Validate whether required dispatcher fields are present for `test_type`.
        Returns dict: {'ok': bool, 'missing': [..], 'test_type': str, 'max_rotation_deg': float|None}.
        """
        tt = str(test_type or getattr(self, 'test_type', 'dynamic'))
        if tt in ('isometric', 'isovelocity'):
            self._normalize_dispatch_aliases()
        missing = []
        if tt == 'isometric':
            for k in ['isometric_initial', 'isometric_final', 'isometric_num_steps']:
                if getattr(self, k, None) is None:
                    missing.append(k)
            ns = getattr(self, 'isometric_num_steps', None)
            if ns is not None:
                try:
                    if int(ns) < 1:
                        missing.append('isometric_num_steps (must be ≥ 1)')
                except (TypeError, ValueError):
                    missing.append('isometric_num_steps (invalid)')
            if not self._clamp_spacing_mm_valid():
                missing.append('dclamp (mm)')
            _iso_mode = str(
                getattr(self, 'isometric_mode', None) or getattr(self, 'mode', None) or 'strain'
            ).lower()
            if _iso_mode in ('strain', 'strain_rate', 'strain_pct', 'strain_pct_rate'):
                if not self._xsec_width_mm_valid():
                    missing.append('xsec_width (mm)')
        elif tt == 'isovelocity':
            for k in [
                'isovelocity_min_vel', 'isovelocity_max_vel',
                'isovelocity_starting_strain', 'isovelocity_num_steps',
            ]:
                if getattr(self, k, None) is None:
                    missing.append(k)
            ns = getattr(self, 'isovelocity_num_steps', None)
            if ns is not None:
                try:
                    if int(ns) < 1:
                        missing.append('isovelocity_num_steps (must be ≥ 1)')
                except (TypeError, ValueError):
                    missing.append('isovelocity_num_steps (invalid)')
            if not self._clamp_spacing_mm_valid():
                missing.append('dclamp (mm)')
            _iv_mode = str(getattr(self, 'isovelocity_starting_strain_mode', None) or 'strain').lower()
            if _iv_mode in ('strain', 'strain_rate', 'strain_pct', 'strain_pct_rate'):
                if not self._xsec_width_mm_valid():
                    missing.append('xsec_width (mm)')
            _iv_vel_mode = str(getattr(self, 'isovelocity_velocity_mode', None) or 'angle_vel').lower()
            if _iv_vel_mode in ('strain_rate', 'strain_pct_rate'):
                if not self._xsec_width_mm_valid():
                    missing.append('xsec_width (mm) (required for strain-rate velocity mode)')
        elif tt == 'calibration':
            if getattr(self, 'calibration_base_test_type', None) is None:
                missing.append('calibration_base_test_type')
        elif tt in ('dynamic', 'frequency_sweep', 'frequency_step', 'curvature_step'):
            def _seq_missing(seq) -> bool:
                if seq is None:
                    return True
                try:
                    if np.asarray(seq).size == 0:
                        return True
                except Exception:
                    try:
                        return len(seq) == 0
                    except Exception:
                        return True
                return False

            af = getattr(self, 'all_freqs', None)
            ac = getattr(self, 'all_curves', None)
            if _seq_missing(af):
                missing.append('all_freqs (Hz)')
            else:
                try:
                    fa = np.asarray(af, dtype=float).ravel()
                    if fa.size == 0 or not np.all(np.isfinite(fa)) or np.any(fa <= 0):
                        missing.append('all_freqs (finite values > 0 Hz)')
                except Exception:
                    missing.append('all_freqs (valid numeric list)')
            # Amplitudes are valid whether supplied as all_curves (κ) or as all_amps + mode; the
            # motion path converts all_amps→all_curves via _ensure_all_curves_for_run at run time,
            # so requiring all_curves here wrongly failed sweeps that only set all_amps. Stim fields
            # are never required (stim is optional), so stim-OFF runs are unaffected (3.1).
            aa = getattr(self, 'all_amps', None)
            aa_mode = getattr(self, 'all_amps_mode', None) or getattr(self, 'curve_input_mode', None)
            if _seq_missing(ac) and (_seq_missing(aa) or aa_mode is None):
                missing.append('all_curves / amplitudes')
            if tt != 'dynamic':
                du = getattr(self, 'duration', None)
                try:
                    duf = float(du) if du is not None else None
                except (TypeError, ValueError):
                    duf = None
                if duf is None or not np.isfinite(duf) or duf <= 0:
                    missing.append('duration (s)')
            if tt == 'dynamic':
                cps = getattr(self, 'cycles_per_step', None)
                try:
                    cpsi = int(cps) if cps is not None else None
                except (TypeError, ValueError):
                    cpsi = None
                if cpsi is None or cpsi <= 0:
                    missing.append('cycles_per_step')
                nec = getattr(self, 'n_end_cycles', None)
                try:
                    neci = int(nec) if nec is not None else None
                except (TypeError, ValueError):
                    neci = None
                if neci is None or neci < 0:
                    missing.append('n_end_cycles')
            if not self._clamp_spacing_mm_valid():
                missing.append('dclamp (mm)')
            if not self._xsec_width_mm_valid():
                missing.append('xsec_width (mm)')
        elif tt == 'step_change':
            for attr, label in (
                ('step_change_frequencies', 'step_change_frequencies'),
                ('step_change_curves', 'step_change_curves'),
                ('step_change_cycles_per_step', 'step_change_cycles_per_step'),
            ):
                seq = getattr(self, attr, None)
                if seq is None:
                    missing.append(label)
                    continue
                try:
                    if np.asarray(seq).size == 0:
                        missing.append(label)
                except Exception:
                    try:
                        if len(seq) == 0:
                            missing.append(label)
                    except Exception:
                        missing.append(label)
            if not self._clamp_spacing_mm_valid():
                missing.append('dclamp (mm)')
            if not self._xsec_width_mm_valid():
                missing.append('xsec_width (mm)')
        ok = len(missing) == 0
        max_rotation_deg = None
        if ok:
            try:
                from shared.utilities import compute_max_rotation_deg

                max_rotation_deg = float(compute_max_rotation_deg({'bender': self, 'test_type': tt}))
                self.max_commanded_rotation_deg = max_rotation_deg
            except Exception as exc:
                logging.warning('validate_dispatch_setup: max rotation not computed for %s: %s', tt, exc)
                self.max_commanded_rotation_deg = None
        else:
            self.max_commanded_rotation_deg = None
        return {
            'ok': ok,
            'missing': missing,
            'test_type': tt,
            'max_rotation_deg': max_rotation_deg,
        }

    def make_cycle_tags(self):
        aidata = getattr(self, 'aidata', None)
        if aidata is None:
            raise ValueError('make_cycle_tags requires aidata')
        arr = np.asarray(aidata)
        if arr.ndim == 2:
            n_ai = int(arr.shape[1])
        elif arr.ndim == 1:
            n_ai = int(arr.size)
        else:
            raise ValueError(f'make_cycle_tags: aidata must be 1D or 2D; got ndim={arr.ndim}')
        t_arr = np.asarray(getattr(self, 't', np.array([])), dtype=float).reshape(-1)
        total_pts = int(t_arr.size) if t_arr.size > 0 else n_ai
        if t_arr.size > 0 and n_ai > 0 and t_arr.size != n_ai:
            import warnings
            warnings.warn(
                f'make_cycle_tags: len(t)={t_arr.size} != aidata samples={n_ai}; using len(t).',
                UserWarning,
                stacklevel=2,
            )
            total_pts = int(t_arr.size)
        freq_by_cycle = getattr(self, 'freq_by_cycle', None)
        if freq_by_cycle is None:
            raise AttributeError('make_cycle_tags requires freq_by_cycle')
        cycle_tag = np.full(total_pts, -1, dtype=int)  # -1 = pre/post motion (not a numbered cycle)
        
        # 2. Convert Pre-Stim Time to Points
        pre_time = abs(getattr(self, 'prestim_time', 0)) 
        pre_pts = int(pre_time * self.daq_ai_sample_rate_hz)
        
        # 3. Tag Active Cycles (Starting at 0 to match 22-element metadata)
        # We start 'current_pos' after the pre_pts (which remain -1)
        current_pos = pre_pts
        
        for i, freq in enumerate(freq_by_cycle):
            cycle_num = i  # 0, 1, 2... 21
            pts = int(round(self.daq_ai_sample_rate_hz / freq))
            end_pos = current_pos + pts
            
            # Safety check: don't overshoot
            if end_pos > total_pts: 
                end_pos = total_pts
                
            cycle_tag[current_pos:end_pos] = cycle_num
            current_pos = end_pos
            
            if current_pos >= total_pts: 
                break
                
        # Store it for the H5 saver
        self.cycle_index_history = cycle_tag