import numpy as np
import importlib
from scipy import interpolate
from datetime import datetime
from copy import copy
import time
import re
import os
import xml.etree.ElementTree as ElementTree
import json
print(f"DEBUG: Loading bender_functions.py from: {os.path.abspath(__file__)}")
import logging

try:
    import nidaqmx.constants as daq
    from nidaqmx import Task
    from nidaqmx.stream_writers import AnalogMultiChannelWriter, DigitalSingleChannelWriter
    from nidaqmx.stream_readers import AnalogMultiChannelReader, CounterReader
    from nidaqmx.errors import DaqError
    from nidaqmx.constants import TerminalConfiguration
except ImportError:
    logging.warning('No DAQmx available')
    
import xml.etree.ElementTree as ElementTree

import h5py


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
    """Map normalized time u in [0,1] to progress in [0,1]."""
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
        self.daq_ai_sample_rate_hz = cfg.daq_ai_sample_rate_hz
        self.daq_ao_do_sample_rate_hz = cfg.daq_ao_do_sample_rate_hz
        self.motor_gear_ratio = cfg.motor_gear_ratio
        self.motor_full_steps_per_rev = cfg.motor_full_steps_per_rev
        self.encoder_pulses_per_rev = cfg.encoder_pulses_per_rev
          
          
        # Grab the sensor lists
        self.input_channels = cfg.input_channels
        self.input_channel_names = cfg.input_channel_names

        # 3. Assign Calibration/Directionality from cfg
        self.forcetorque_calibration_file = cfg.forcetorque_calibration_file
        self.positive_motor_direction = cfg.positive_motor_direction
        self.S1side = cfg.S1side
        self.S2side = cfg.S2side
        self.motor_axis = cfg.motor_axis
        self.bending_axis_sensor = cfg.bending_axis_sensor

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
        self.velocity_exponent = 1.0
        self.ramp_mode_default = getattr(cfg, 'ramp_mode_default', 'linear')
        self.master_logger = MasterLogger()
        self.h5_protocol_metadata = {}

        # Standard 2D shapes for NI-DAQmx (Channels, Samples)
        self.stimcmdhi = np.zeros((2, 2))
        self.dig = np.zeros((1, 2), dtype='uint32')

        print(f"Bender initialized using: {config_module_name}.py")

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

    def run_experiment(self, test_type = "dynamic"):


       # --- THE CLEANING CREW ---
        self.aidata = None
        self.forcetorque = None
        self.sono_left_mm = None
        self.sono_right_mm = None
        self.angledata = None
        self.master_logger.clear()

        # Now, if the DAQ fails, you won't accidentally save old data!
        # Set input and output names/channels
        self.set_input_channels(input_channels=self.input_channels, input_channel_names=self.input_channel_names)
        self.set_stim_channels(*self.stim_channels) 

        # self.calibration is assumed to be loaded in the notebook using bender.loadCalibration(...)

        # General parameters (duration/sample_rate passed from notebook)
        sample_rate = self.daq_ai_sample_rate_hz

        # Define angle/anglevel/tnorm variables for scope consistency
        angle, anglevel, tnorm = None, None, None
        duration = None #Initialize, but will be set below

     # THIS LEVEL IS ABOUT CREATING MOTOR ANGLES
        if test_type == 'dynamic':
            if self.period_by_cycle is None:
                raise AttributeError("Dynamic test requires 'period_by_cycle' to be set via organize_cycles first.")

            duration = np.sum(self.period_by_cycle)

            angle, anglevel, tnorm, t = self.make_cycles_dynamic(
                self.period_by_cycle,
                self.freq_by_cycle,
                self.amp_by_cycle,
            )

        elif test_type == 'frequency_sweep':
            if self.duration is None:
                raise AttributeError("frequency_sweep requires 'self.duration' to be set in the notebook first.")

            duration = self.duration + self.waitbefore + self.waitafter

            angle, anglevel, tnorm, t, sweep_freq = self.make_cycles_frequency_sweep(
                self.all_freqs,
                self.all_curves,
                self.amplitude_frequency_exponent,
                self.duration,
                self.waitbefore,
                nominal_frequency=getattr(self, 'sweep_nominal_frequency', None),
                nominal_curvature=getattr(self, 'sweep_nominal_curvature', None),
            )
            self.sweep_instantaneous_freq = sweep_freq

        elif test_type == 'frequency_step':
            if self.duration is None:
                raise AttributeError("frequency_step requires 'self.duration' to be set in the notebook first.")

            duration = self.duration + self.waitbefore + self.waitafter

            angle, anglevel, tnorm, t, _ = self.make_cycles_frequency_step(
                self.all_freqs,
                self.all_curves,
                self.duration,
                self.waitbefore,
                nominal_frequency=getattr(self, 'step_nominal_frequency', None),
                nominal_curvature=getattr(self, 'step_nominal_curvature', None),
            )

        elif test_type == 'curvature_step':
            if self.duration is None:
                raise AttributeError("curvature_step requires 'self.duration' to be set in the notebook first.")

            duration = self.duration + self.waitbefore + self.waitafter

            angle, anglevel, tnorm, t, _ = self.make_cycles_curvature_step(
                self.all_freqs,
                self.all_curves,
                self.duration,
                self.waitbefore,
                nominal_frequency=getattr(self, 'step_nominal_frequency', None),
                nominal_curvature=getattr(self, 'step_nominal_curvature', None),
            )

        elif test_type == 'step_change':
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
                f"'frequency_step', 'curvature_step', or 'step_change')"
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

        self.master_logger.record(test_type=test_type)
        self.h5_protocol_metadata = self.master_logger.as_dict()

        # Create motor stepper pulses based on the generated angle/anglevel signals (MOTION ONLY)
        self.make_motor_stepper_pulses(
                        daq_ao_do_sample_rate_hz=self.daq_ao_do_sample_rate_hz,
                        motor_gear_ratio=self.motor_gear_ratio,
                        motor_full_steps_per_rev=self.motor_full_steps_per_rev
                    )
                            

        # Print file save location
        filename = self.increment_file_name(f'experiment_data_{test_type}_000.h5')
        print(f"Data will be saved to: {filename}")

        # Run the experiment using 'self'
        self.aidata = self.run(device_name=self.device_name)

        # --- 1. Process Force/Torque (Always assumed to exist) ---
        # Find exactly where the 6 SG channels are, regardless of order
        sg_names = ['xForce', 'yForce', 'zForce', 'xTorque', 'yTorque', 'zTorque']
        forcetorque_indices = [self.input_channel_names.index(n) for n in sg_names if n in self.input_channel_names]

        if len(forcetorque_indices) == 6:
            self.forcetorque = self.apply_calibration_forcetorque(self.aidata[forcetorque_indices, :])
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

   
    def get_data_by_name(self, name):
        """Returns the data row for a specific channel name."""
        try:
            idx = self.input_channel_names.index(name)
            return self.aidata[idx, :]
        except (ValueError, AttributeError):
            return None



    def organize_cycles(self, all_curves, all_freqs, randomize, cycles_per_step, n_end_cycles, dclamp, xsec_width, stim_cycles_in_step, all_stimduties, all_stimphases, stim_pulse_rate):
        start = time.time()
        self.dclamp = float(dclamp)
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
        all_strains_arr = (xsec_width / 2 / 1000) * all_curves_arr
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

        if np.any(np.array(stim_cycles_in_step) >= cycles_per_step):
            raise IndexError("stim_cycles_in_step have to be less than cycles_in_step")

        c = np.arange(0, cycles_per_step)
        is_stim_cycle = np.isin(c, stim_cycles_in_step)
        is_stim_cycle = np.tile(is_stim_cycle, len(all_freqs_arr))
        is_stim_cycle = np.concatenate((is_stim_cycle, [False] * n_end_cycles))

        stimburstdur = duty_by_cycle / freq_by_cycle
        stimburstdur = np.floor(stimburstdur * stim_pulse_rate * 2) / (stim_pulse_rate * 2)
        stimburstdur[is_stim_cycle == False] = 0

    # --- 6. STORE RESULTS (For BOTH Motor Control & H5 Metadata) ---
        self.period_by_cycle = 1.0 / freq_by_cycle
        self.freq_by_cycle   = freq_by_cycle
        self.amp_by_cycle    = amp_by_cycle
        
        # Use the original names your motor controller expects
        self.all_degs        = all_degs_arr
        self.all_strains     = all_strains_arr
        self.all_strainrates = all_strainrates_arr
        self.duty_by_cycle  = all_stimduties_arr
        self.phase_by_cycle = all_stimphases_arr
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
                

    def record_motor_signal(self, t, angle, anglevel, tnorm=None):
        self.daq_ai_sample_rate_hz = round(1.0 / (t[1] - t[0]))

        self.t = t
        self.angle = angle
        self.anglevel = anglevel

        if tnorm is None:
            self.tnorm = copy(t)
        else:
            self.tnorm = tnorm
    
    def record_stim_signal(self, S1stimcmd, S2stimcmd):
        self.S1stimcmd = S1stimcmd
        self.S2stimcmd = S2stimcmd

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

        stepnum = np.floor(poshi / stepsize)
        dstep = np.diff(stepnum)
        motorstep = np.concatenate((np.array([0], dtype='uint8'), (dstep != 0).astype('uint8')))
        motordirection = (velhi <= 0).astype('uint8')

        # Change enable back to ones_like (High = 5V)
        motorenable = np.ones_like(motordirection, dtype='uint8')

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

    def run(self, device_name):
        input_channels = ['/'.join((device_name, c1)) for c1 in self.input_channels]
        S1stim_chan = '/'.join((device_name, self.S1stim_chan))
        S2stim_chan = '/'.join((device_name, self.S2stim_chan))
        motor_port = '/'.join((device_name, self.motor_port))
        encoder_chan = '/'.join((device_name, self.encoder_chan))

        with Task() as analog_in, Task() as analog_out, \
                Task() as digital_out, Task() as angle_in:
                       # set up the input channels
            for c1, name1 in zip(input_channels, self.input_channel_names):
                # Check for 'sono' to set RSE mode, otherwise use Differential
                if 'sono' in name1.lower():
                    t_config = TerminalConfiguration.RSE # Need to reconfigure for Sonometrics DAC output channels
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

            # start everthing
            # make sure to start the output first, because it'll wait until the input starts
            digital_out.start()
            analog_out.start()
            angle_in.start()
            analog_in.start()
            
            # wait until we're done, record the time
            analog_in.wait_until_done(self.t[-1]+10)
            self.endTime = datetime.now()
            
            # and read the data
            reader.read_many_sample(self.aidata)
            angle_reader.read_many_sample_double(self.angledata)


            self.angle_measured = self.angledata  # This makes bender.angle_measured available

        return(self.aidata)

#
    def _uniform_cycles_from_duration(self, duration, f0, amp_deg):
        """Build period/freq/amp arrays so ``make_cycles_dynamic`` runs uniform f and amplitude for ~``duration`` s."""
        n = max(1, int(round(float(duration) * float(f0))))
        p = 1.0 / float(f0)
        period_by_cycle = np.full(n, p, dtype=float)
        freq_by_cycle = np.full(n, float(f0), dtype=float)
        amp_by_cycle = np.full(n, float(amp_deg), dtype=float)
        return period_by_cycle, freq_by_cycle, amp_by_cycle

    def make_cycles_frequency_step(self, all_freqs, all_curves, duration, waitbefore,
                                   nominal_frequency=None, nominal_curvature=None):
        dclamp = getattr(self, 'dclamp', None)
        if dclamp is None:
            raise ValueError("frequency_step requires self.dclamp (set via organize_cycles).")
        fq_src = nominal_frequency if nominal_frequency is not None else all_freqs
        cq_src = nominal_curvature if nominal_curvature is not None else all_curves
        f0, f1 = _normalize_start_end(fq_src)
        c0, c1 = _normalize_start_end(cq_src)
        if abs(f1 - f0) > 1e-9 or abs(c1 - c0) > 1e-9:
            raise ValueError(
                "frequency_step expects a single Hz and 1/m setpoint (float or equal [start, end])."
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
        self._protocol_log('frequency_step', f0, f0, c0, c0, motion_duration_s=float(duration))
        return angle, anglevel, tnorm, freq, t

    def make_cycles_curvature_step(self, all_freqs, all_curves, duration, waitbefore,
                                   nominal_frequency=None, nominal_curvature=None):
        dclamp = getattr(self, 'dclamp', None)
        if dclamp is None:
            raise ValueError("curvature_step requires self.dclamp (set via organize_cycles).")
        fq_src = nominal_frequency if nominal_frequency is not None else all_freqs
        cq_src = nominal_curvature if nominal_curvature is not None else all_curves
        f0, f1 = _normalize_start_end(fq_src)
        c0, c1 = _normalize_start_end(cq_src)
        if abs(f1 - f0) > 1e-9 or abs(c1 - c0) > 1e-9:
            raise ValueError(
                "curvature_step expects a single Hz and 1/m setpoint (float or equal [start, end])."
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
        """Log-frequency sweep with curvature-based amplitude scaling (exponent matches legacy sweep)."""
        daq_ai_hz = self.daq_ai_sample_rate_hz
        dclamp = getattr(self, 'dclamp', None)
        if dclamp is None:
            raise ValueError("make_cycles_frequency_sweep requires self.dclamp (set via organize_cycles).")
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

        # Calculate timings and durations
        movedur = np.sum(period_by_cycle)
        totaldur = waitbefore + movedur + waitafter
        t = np.arange(0, totaldur, 1.0 / daq_ai_hz)
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
        freqs = np.asarray(frequencies, dtype=float).reshape(-1)
        curv = np.asarray(curves, dtype=float).reshape(-1)
        cps = np.asarray(cycles_per_step, dtype=int).reshape(-1)
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

        saved_degs = self.all_degs
        saved_freqs = self.all_freqs
        saved_v = self.amp_step_vel
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
            self._protocol_log(
                'step_change', f_lo, f_hi, c_lo, c_hi,
                motion_duration_s=movedur,
                step_change_blocks=int(len(freqs)),
                step_change_total_cycles=int(len(allfreqs)),
                dclamp_mm=float(dc),
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


        # 3. HEAVY LIFTING
        # Use the passed-in basis if available, otherwise fallback to self
        t = t_basis if t_basis is not None else self.t
        tnorm = tnorm_basis if tnorm_basis is not None else self.tnorm
        
        S1stimcmd = np.zeros_like(t)
        S2stimcmd = np.zeros_like(t)
        Lonoff, Ronoff = [], []
        
        # Calculate the pulse wave once
        # Using 5.0 for the stimulator 'On' voltage
        pulse_wave = (np.mod(t * stim_pulse_rate, 1) <= 0.5).astype(float) * 5.0
        bendphase = tnorm - 0.25
    
        # 4. Optimized Cycle Loop
        for c, (dur1, duty1, f1, p1) in enumerate(zip(stimburstdur, duty_by_cycle, freq_by_cycle, phase_by_cycle)):
            if dur1 <= 0: continue
            
            # Left Side (S1)
            t_s1 = c + p1
            m1 = (bendphase >= t_s1) & (bendphase < t_s1 + duty1)
            if np.any(m1):
                S1stimcmd[m1] = pulse_wave[m1]
                t_sub = t[m1]
                Lonoff.append([t_sub[0], t_sub[-1]])

            # Right Side (S2)
            t_s2 = c + 0.5 + p1
            m2 = (bendphase >= t_s2) & (bendphase < t_s2 + duty1)
            if np.any(m2):
                S2stimcmd[m2] = pulse_wave[m2]
                t_sub2 = t[m2]
                Ronoff.append([t_sub2[0], t_sub2[-1]])

        # 5. Pre-stimulation (The "Start" burst)
        m_pre = (t >= prestim_time) & (t < (prestim_time + prepoststim_dur))
        if np.any(m_pre):
            S1stimcmd[m_pre] = pulse_wave[m_pre]
            Lonoff.append([prestim_time, prestim_time + prepoststim_dur])

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
        raw_torque = raw_data_dict[self.bending_axis_sensor]
        
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
        for key, value in kwargs.items():
            setattr(self, key, value)
            print(f"  Stored: {key} = {value}")

    def make_cycle_tags(self):
   
        # 1. Total samples from your data matrix [samples, channels]
        total_pts = self.aidata.shape[0]
        cycle_tag = np.full(total_pts, -1, dtype=int) # Initialize all as -1 (Pre/Post)
        
        # 2. Convert Pre-Stim Time to Points
        pre_time = abs(getattr(self, 'prestim_time', 0)) 
        pre_pts = int(pre_time * self.daq_ai_sample_rate_hz)
        
        # 3. Tag Active Cycles (Starting at 0 to match 22-element metadata)
        # We start 'current_pos' after the pre_pts (which remain -1)
        current_pos = pre_pts
        
        # Using 'freq_by_cycle' which you already organized
        for i, freq in enumerate(self.freq_by_cycle):
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