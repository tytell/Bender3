# jimenez_bender_config_A.py

# --- Calibration & Directionality ---
# Load sensor calibration file. Make sure max range/sensitivity match animal size!
forcetorque_calibration_file = 'FT56491.cal' 

# Does positive angle command make bender go left or right? (Depends on mounting/settings)
positive_motor_direction = "right"

# Specimen lateral axis (one number fixes both sides): signed index for the side named above.
# The opposite anatomical side gets the negated index. Non-zero (typically ±1). With
# positive_motor_direction="left" and -1 here, specimen LEFT = -1 and RIGHT = +1.
specimen_lateral_index_on_positive_motor_side = -1

# User Configuration
motor_axis = "z"           # Motor physically rotates along global 'Y' axis
bending_axis_sensor = "z" # Sensor's 'X' is actually the motor's rotation
primary_bending_axis = "zTorque"  # Preferred torque axis for QC plots/correction: xTorque|yTorque|zTorque
bending_axis_specimen = "lateral" # "dorsoventral", "lateral", or "anteroposterior"

S1side = 'left' # Double check stimulator channel 1 side!
S2side = 'right'

# --- ASSIGN DAQ and Motor Parameters ---
daq_ai_sample_rate_hz = 1000.0   # DAQ AI + encoder sample clock (Hz)
daq_ao_do_sample_rate_hz = 60000.0  # DAQ AO stim + DO motor stream (Hz)
motor_full_steps_per_rev = 1600  # Motor steps per revolution (e.g., 1/8 microstepping)
motor_gear_ratio = 5  # e.g., 5:1 — motor revolutions per one output revolution (gearbox)
# --- ASSIGN DAQ Hardware Ports ---
stim_channels = ["ao0", "ao1"]
motor_port = "port0"
encoder_chan = "ctr0" 
device_name = "Dev1"

# Add strain gauge input channels (if applicable) for six-axis force transducer (e.g., ATI Nano40). Make sure to assign correct channels and names based on your specific setup!
SG_chan = ['ai0', 'ai1', 'ai2', 'ai3', 'ai4', 'ai5']
SG_name = ['xForce', 'yForce', 'zForce', 'xTorque', 'yTorque', 'zTorque']

# Add stim monitor channel (if applicable) from S88 stimulator. Make sure to assign correct channel and name based on your specific setup!
stim_monitor_chan = ['ai7']
stim_monitor_name = ['stim_monitor']

# Add sonomicrometry channels from Sonometrics DS3 (if applicable)
use_sono = True
sono_channel = ["ai6"] # If using sonomicrometry, assign output channels for sonomicrometer excitation
sono_name = ["sono_right"]
sono_internal_samplefreq = 241 # Internal sample rate of the sonomicrometry system (e.g., 981 or 251 Hz for Sonometrics DS3)
# --- Sonometer Calibration (Linear: Volts to mm) ---
# Format: [Low_Volts, High_Volts, Low_mm, High_mm]
sono_cal_left = [1.1, 4.5, 11.8, 47] 
sono_cal_right = [1.1, 4.5, 11.8, 47]

# Combine all input channels and names into lists for Bender configuration
# Comment out any channels that you don't plant o use in the two lines below!! i.e., if not using sonomicrometry, comment out sono_channel and sono_name lines below.
input_channels = SG_chan + (sono_channel if use_sono else []) + list(stim_monitor_chan)
input_channel_names = SG_name + (sono_name if use_sono else []) + list(stim_monitor_name)

# --- Advanced / Stimulation Timing ---
amp_step_vel = 10 
encoder_pulses_per_rev = 10000  # E6 optical encoder (NI-DAQ: pulses_per_rev)

# Protocol ramps: 'linear' or 'exponential' (for code that calls _ramp_progress). This is
# separate from amplitude_frequency_exponent on the sweep — see NOTE in bender_run / bender_functions.
ramp_mode_default = 'linear'

# --- Time Buffers ---
waitbefore = 3.0 # Seconds to wait before bending after stimulation starts
waitafter = 4.0  # Seconds to wait after bending
rampdur = 0.25   # Seconds to ramp on/off motor motion

# duty of 0.3 at 5 Hz (Isometric tests)
prepoststim_dur = 0.06
prepoststim_sep = 1.0           # Time between left and right bursts
prestim_time = -2               # Time prestim left burst starts
poststim_time = 2               # Time *after* end of bending


# This dictionary is used by the HDF5 saver to label your data for R analysis.
units = {
    # Hardware & Sampling
    'daq_ai_sample_rate_hz': 'Hz',
    'daq_ao_do_sample_rate_hz': 'Hz',
    'motor_full_steps_per_rev': 'steps/rev',
    'motor_gear_ratio': 'multiplier',
    'encoder_pulses_per_rev': 'pulses/rev',
    
    # Timing & Buffers
    'waitbefore': 'seconds',
    'waitafter': 'seconds',
    'rampdur': 'seconds',
    'duration': 'seconds',
    
    # Physics & Results
    'xForce': 'Newtons',
    'yForce': 'Newtons',
    'zForce': 'Newtons',
    'xTorque': 'N-m',
    'yTorque': 'N-m',
    'zTorque': 'N-m',
    'angle_measured': 'degrees',
    'planned_angles': 'degrees',
    
    # Specimen & Setup
    'body_thickness': 'meters',
    'xsec_width': 'mm',
    'desired_curves': '1/m',
    'desired_strain_pct': 'percent',
    
    # Sonomics & Stimulation
    'sono_left': 'mm',
    'sono_right': 'mm',
    'S1volts': 'Volts',
    'S2volts': 'Volts',
    'stim_pulse_rate': 'Hz',
    'all_stimduties': 'fraction_of_cycle',
    'all_stimphases': 'fraction_of_cycle',
    'velocity_exponent': 'dimensionless',
    'protocol': 'string',
    'test_type': 'string',
    'frequency_hz': 'Hz',
    'curvature_1_per_m': '1/m',
    'motion_duration_s': 'seconds',
    'dynamic_movedur_s': 'seconds',
    'amplitude_frequency_exponent': 'dimensionless',
    'step_change_blocks': 'count',
    'step_change_total_cycles': 'count',
    'dclamp_mm': 'mm',
    'specimen_lateral_index_on_positive_motor_side': 'dimensionless',
}

unit_rules = {
    'width': 'mm',
    'height': 'mm',
    'thickness': 'mm',
    'length': 'mm',
    'depth': 'mm',
    'mass': 'grams',
    'weight': 'grams',
    'density': 'g/cm^3',
}