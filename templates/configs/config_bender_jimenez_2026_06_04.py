# -*- coding: utf-8 -*-
from config_bender_jimenez_2026_06_04 import *

S1side = 'left'
S2side = 'right'
SG_chan = ['ai0', 'ai1', 'ai2', 'ai3', 'ai4', 'ai5']
SG_name = ['xForce', 'yForce', 'zForce', 'xTorque', 'yTorque', 'zTorque']
amp_step_vel = 10
bending_axis_sensor = 'z'
bending_axis_specimen = 'lateral'
daq_ai_sample_rate_hz = 1000.0
daq_ao_do_sample_rate_hz = 60000.0
device_name = 'Dev1'
encoder_chan = 'ctr0'
encoder_pulses_per_rev = 10000
forcetorque_calibration_file = 'FT56491.cal'
motor_axis = 'z'
motor_full_steps_per_rev = 1600
motor_gear_ratio = 5
motor_port = 'port0'
positive_motor_direction = 'left'
poststim_time = 2.0
prepoststim_dur = 0.06
prepoststim_sep = 1.0
prestim_time = -2.0
primary_bending_axis = 'zTorque'
ramp_mode_default = 'linear'
rampdur = 0.25
sono_cal_left = [1.1, 4.5, 11.8, 47.0]
sono_cal_right = [1.1, 2.3, 3.4, 4.5, 11.8, 23.5, 35.3, 47.0]
sono_channel = ['ai6']
sono_internal_samplefreq = 242.2
sono_name = ['sono_right']
specimen_lateral_index_on_positive_motor_side = -1
stim_channels = ['ao0', 'ao1']
stim_monitor_chan = ['ai7']
stim_monitor_name = ['stim_monitor']
use_sono = True
waitafter = 4.0
waitbefore = 3.0

# Combined AI lists after overrides (stim monitor appended when listed)
input_channels = SG_chan + (sono_channel if use_sono else []) + list(stim_monitor_chan)
input_channel_names = SG_name + (sono_name if use_sono else []) + list(stim_monitor_name)
