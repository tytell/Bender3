# jimenez_bender_config_A.py

# --- Calibration & Directionality ---
# Load sensor calibration file. Make sure max range/sensitivity match animal size!
calibration_file = 'FT56491.cal' 

# Does positive angle command make bender go left or right? (Depends on mounting/settings)
positive_motor_direction = "left"     

# User Configuration
motor_axis = "y"           # Motor physically rotates along global 'Y' axis
bending_axis_sensor = "x" # Sensor's 'X' is actually the motor's rotation
bending_axis_specimen = "dorsoventral" # "dorsoventral", "lateral", or "anteroposterior"

S1side = 'left' # Double check stimulator channel 1 side!
S2side = 'right'

# --- DAQ and Motor Parameters ---
samplefreq = 1000.0   # DAQ sample frequency
outputfreq = 100000.0 # DAQ output frequency
stepsperrev = 1600    # Motor steps per revolution (e.g., 1/8 microstepping)

# --- DAQ Hardware Ports ---
stim_channels = ["ao0", "ao1"]
motor_port = "port0"
encoder_chan = "ctr0" 
device_name = "Dev1"


# --- Advanced / Stimulation Timing ---
amp_step_vel = 10 
encoder_counts_per_rev = 10000 # E5 optical encoder (1000 PPR)

# --- Time Buffers ---
waitbefore = 3.0 # Seconds to wait before bending after stimulation starts
waitafter = 4.0  # Seconds to wait after bending
rampdur = 0.25   # Seconds to ramp on/off motor motion

# duty of 0.3 at 5 Hz (Isometric tests)
prepoststim_dur = 0.3 / 5       
prepoststim_sep = 1             # Time between left and right bursts
prestim_time = -2               # Time prestim left burst starts
poststim_time = 2               # Time *after* end of bending
