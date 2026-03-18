# jimenez_bender_config_A.py

# --- Calibration & Directionality ---
# Load sensor calibration file. Make sure max range/sensitivity match animal size!
forcetorque_calibration_file = 'FT56491.cal' 

# Does positive angle command make bender go left or right? (Depends on mounting/settings)
positive_motor_direction = "left"     

# User Configuration
motor_axis = "y"           # Motor physically rotates along global 'Y' axis
bending_axis_sensor = "x" # Sensor's 'X' is actually the motor's rotation
bending_axis_specimen = "dorsoventral" # "dorsoventral", "lateral", or "anteroposterior"

S1side = 'left' # Double check stimulator channel 1 side!
S2side = 'right'

# --- ASSIGN DAQ and Motor Parameters ---
samplefreq = 1000.0   # DAQ sample frequency
outputfreq = 60000.0 # DAQ output frequency
stepsperrev = 1600    # Motor steps per revolution (e.g., 1/8 microstepping)
gear_ratio = 5 # Gear ratio of the motor (e.g., 5:1 means 5 motor revolutions = 1 output revolution). Depends on gear box!
# --- ASSIGN DAQ Hardware Ports ---
stim_channels = ["ao0", "ao1"]
motor_port = "port0"
encoder_chan = "ctr0" 
device_name = "Dev1"

# Add strain gauge input channels (if applicable) for six-axis force transducer (e.g., ATI Nano40). Make sure to assign correct channels and names based on your specific setup!
SG_chan = ['ai0', 'ai1', 'ai2', 'ai3', 'ai4', 'ai5']
SG_name = ['xForce', 'yForce', 'zForce', 'xTorque', 'yTorque', 'zTorque']

# Add stim monitor channel (if applicable) from S88 stimulator. Make sure to assign correct channel and name based on your specific setup!
stim_monitor_chan = ['ai8']
stim_monitor_name = ['stim_monitor']

# Add sonomicrometry channels from Sonometrics DS3 (if applicable)
sono_channel = ["ai6", "ai7"] # If using sonomicrometry, assign output channels for sonomicrometer excitation
sono_name = ["sono_left", "sono_right"]
sono_internal_samplefreq = 241 # Internal sample rate of the sonomicrometry system (e.g., 981 or 251 Hz for Sonometrics DS3)
# --- Sonometer Calibration (Linear: Volts to mm) ---
# Format: [Low_Volts, High_Volts, Low_mm, High_mm]
sono_cal_left = [1.1, 4.5, 11.8, 47] 
sono_cal_right = [1.1, 4.5, 11.8, 47]

# Combine all input channels and names into lists for Bender configuration
# Comment out any channels that you don't plant o use in the two lines below!! i.e., if not using sonomicrometry, comment out sono_channel and sono_name lines below.
input_channels = SG_chan + sono_channel #+  stim_monitor_chan 
input_channel_names = SG_name + sono_name #+ stim_monitor_name

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
