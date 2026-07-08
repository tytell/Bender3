# jimenez_bender_config_A.py
# Test fixture config — standalone, no inheritance. Used by the test suite.

forcetorque_calibration_file = 'FT56491.cal'

# Apparatus-inertia calibration artifact (JSON from fit_apparatus_inertia.py). Machine-specific,
# same convention as forcetorque_calibration_file: put the file on this machine, then give an
# absolute path OR a bare filename resolved against the launch dir / repo root. A missing file is
# skipped silently (no auto-load); the GUI loader can still set one per session. This bare filename
# resolves against the repo root, where the committed calibration artifact lives, so Bender
# auto-loads it on init (mirrors how forcetorque_calibration_file loads FT56491.cal).
apparatus_inertia_calibration_file = '2026-07-07_apparatus_inertia_calibration.json'

positive_motor_direction = "right"

specimen_lateral_index_on_positive_motor_side = -1

apparatus_id = 'bender'

motor_axis = "z"
bending_axis_sensor = "z"
primary_bending_axis = "zTorque"
bending_axis_specimen = "lateral"

S1side = 'left'
S2side = 'right'

daq_ai_sample_rate_hz = 1000.0
daq_ao_do_sample_rate_hz = 60000.0
motor_full_steps_per_rev = 1600
motor_gear_ratio = 5

stim_channels = ["ao0", "ao1"]
motor_port = "port0"
encoder_chan = "ctr0"
encoder_pulses_per_rev = 10000
device_name = "Dev1"

SG_chan = ['ai0', 'ai1', 'ai2', 'ai3', 'ai4', 'ai5']
SG_name = ['xForce', 'yForce', 'zForce', 'xTorque', 'yTorque', 'zTorque']

stim_monitor_chan = ['ai7']
stim_monitor_name = ['stim_monitor']

use_sono = True
sono_channel = ["ai6"]
sono_name = ["sono_right"]
sono_internal_samplefreq = 241
sono_cal_left = [1.1, 4.5, 11.8, 47]
sono_cal_right = [1.1, 4.5, 11.8, 47]

input_channels = SG_chan + (sono_channel if use_sono else []) + list(stim_monitor_chan)
input_channel_names = SG_name + (sono_name if use_sono else []) + list(stim_monitor_name)

rampdur = 0.25
waitbefore = 3.0
waitafter = 4.0
prepoststim_dur = 0.06
prepoststim_sep = 1.0
prestim_time = -2.0
poststim_time = 2.0
amp_step_vel = 10
ramp_mode_default = 'linear'
