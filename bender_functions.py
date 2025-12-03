import nidaqmx
import numpy as np
import plotly
import plotly.graph_objects as go
from scipy import interpolate
from datetime import datetime
import time
from copy import copy

import re
import os

import logging

try:
    import nidaqmx.constants as daq
    from nidaqmx import Task
    from nidaqmx.stream_writers import AnalogMultiChannelWriter, DigitalSingleChannelWriter
    from nidaqmx.stream_readers import AnalogMultiChannelReader, CounterReader
    from nidaqmx.errors import DaqError
except ImportError:
    logging.warning('No DAQmx available')
    
import xml.etree.ElementTree as ElementTree

import h5py

class Bender:
    def __init__(self):
        self.S1actcmd = None
        self.S2actcmd = None

    def set_bending_signal(self, t, angle, anglevel, tnorm=None):
        self.samplefreq = 1.0 / (t[1] - t[0])

        self.t = t
        self.angle = angle
        self.anglevel = anglevel

        if tnorm is None:
            self.tnorm = copy(t)
        else:
            self.tnorm = tnorm
    
    def set_activation(self, S1actcmd, S2actcmd):
        self.S1actcmd = S1actcmd
        self.S2actcmd = S2actcmd

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
            raise IOError('Bad calibration file')

        self.calibration = np.array(mat).T

    def applyCalibration(self, rawdata):
        self.forcetorque = np.dot(rawdata[:6, :].T, self.calibration)
        self.forcetorque = self.forcetorque.T

        return self.forcetorque

    def make_motor_stepper_pulses(self, outsampfreq,
                                scale=6.0,
                                stepsperrev=6400.0):

        self.outputfreq = outsampfreq
        self.scale = scale

        tout = np.arange(self.t[0], self.t[-1], 1.0/outsampfreq)

        poshi = interpolate.interp1d(self.t, self.angle, kind='linear', assume_sorted=True, bounds_error=False,
                                    fill_value=0.0)(tout)
        velhi = interpolate.interp1d(self.t, self.anglevel, kind='linear', assume_sorted=True, bounds_error=False,
                                    fill_value=0.0)(tout)

        poshi *= scale
        velhi *= scale

        if outsampfreq == 0 or stepsperrev == 0:
            raise ValueError('Problems with parameters')

        stepsize = 360.0 / stepsperrev
        maxspeed = stepsize * outsampfreq / 2

        if np.any(np.abs(self.anglevel) > maxspeed):
            raise ValueError('Motion is too fast!')

        stepnum = np.floor(poshi / stepsize)
        dstep = np.diff(stepnum)
        motorstep = np.concatenate((np.array([0], dtype='uint8'), (dstep != 0).astype('uint8')))
        motordirection = (velhi <= 0).astype('uint8')

        motorenable = np.ones_like(motordirection, dtype='uint8')
        motorenable[-5:] = 0

        dig = np.packbits(np.column_stack((np.zeros((len(motorstep), 5), dtype=np.uint8),
                                            motorenable,
                                            motorstep,
                                            motordirection)))
        # np.packbits always returns a uint8, so we need to convert to a uint32
        dig = dig.astype('uint32')

        self.tout = tout
        self.dig = dig

        if self.S1actcmd is None:
            self.actcmdhi = np.zeros((2, len(tout)))
        else:
            S1actcmdhi = interpolate.interp1d(self.t, self.S1actcmd, kind='linear', assume_sorted=True, bounds_error=False,
                                        fill_value=0.0)(tout)
            S2actcmdhi = interpolate.interp1d(self.t, self.S2actcmd, kind='linear', assume_sorted=True, bounds_error=False,
                                        fill_value=0.0)(tout)
            self.actcmdhi = np.vstack((S1actcmdhi, S2actcmdhi))

        return tout, dig, motorstep, motordirection

    def set_input_channels(self, inchannels, inchannel_names):
        self.inchannels = inchannels
        self.inchannel_names = inchannel_names

    def set_activation_channels(self, S1activation_chan, S2activation_chan):
        self.S1activation_chan = S1activation_chan
        self.S2activation_chan = S2activation_chan

    def set_motor_channel(self, motor_control_chan):
        self.motor_control_chan = motor_control_chan

    def set_encoder_channel(self, encoder_chan, 
                            counts_per_rev=10000):
        self.encoder_chan = encoder_chan
        self.encoder_counts_per_rev = counts_per_rev

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
                                num_samples, axis_offset_x, axis_offset_z):
            #Crude numerical calculation (voxelization) for the tapered specimen (H=Y, W=X, D=Z).
            #Uses summation I = Sum(m*r^2)

            # x_coords represents the Width axis (left/right, X)
            # z_coords represents the Depth axis (front/back, Z)
        x_coords = np.linspace(-max(front_w_semi, back_w_semi), max(front_w_semi, back_w_semi), num_samples)
        z_coords = np.linspace(0, obj_depth_length, num_samples) 
        
        # Calculate step sizes for voxel volume (dV = dx * dz * H_current)
        dx_step = (x_coords[-1] - x_coords[0]) / (num_samples - 1)
        dz_step = (z_coords[-1] - z_coords[0]) / (num_samples - 1)
        
        total_moi = 0.0
        mass_sum = 0.0
        
        for xi in x_coords:
            for zi in z_coords:
                f = zi / obj_depth_length 
                current_Rx = front_w_semi * (1 - f) + back_w_semi * f # Current half-width (X)
                current_Ry = front_h_semi * (1 - f) + back_h_semi * f # Current half-height (Y)
                
                current_H = current_Ry * 2 

                is_inside = (xi**2 / current_Rx**2) <= 1
                
                if is_inside:
                    mass_per_voxel = rho_eff * dx_step * dz_step * current_H

                    # Distance squared (r^2) to the global rotation axis (in the XZ plane)
                    distance_sq = (xi - axis_offset_x)**2 + (zi - axis_offset_z)**2
                    
                    total_moi += mass_per_voxel * distance_sq
                    mass_sum += mass_per_voxel
        return total_moi, mass_sum
    
    def run(self, device_name):
        inchannels = ['/'.join((device_name, c1)) for c1 in self.inchannels]
        S1activation_chan = '/'.join((device_name, self.S1activation_chan))
        S2activation_chan = '/'.join((device_name, self.S2activation_chan))
        motor_control_chan = '/'.join((device_name, self.motor_control_chan))
        encoder_chan = '/'.join((device_name, self.encoder_chan))

        with Task() as analog_in, Task() as analog_out, \
                Task() as digital_out, Task() as angle_in:
            # set up the input channels
            for c1, name1 in zip(inchannels, self.inchannel_names):
                analog_in.ai_channels.add_ai_voltage_chan(c1, name1)

            # set up the input sample frequency
            # just records as many samples as are in the output
            analog_in.timing.cfg_samp_clk_timing(self.samplefreq,
                                                sample_mode=daq.AcquisitionType.FINITE,
                                                samps_per_chan=len(self.t))

            # set up the encoder channel
            angle_in.ci_channels.add_ci_ang_encoder_chan(encoder_chan, 'encoder',
                                    pulses_per_rev=self.encoder_counts_per_rev)
            angle_in.timing.cfg_samp_clk_timing(self.samplefreq,
                                                source="ai/SampleClock",
                                                sample_mode=daq.AcquisitionType.FINITE,
                                                samps_per_chan=len(self.t))

            # set up the analog output channels
            analog_out.ao_channels.add_ao_voltage_chan(S1activation_chan, 'S1act')
            analog_out.ao_channels.add_ao_voltage_chan(S2activation_chan, 'S2act')
            # it will run much faster than the input channels, because the digital output is linked
            # to it, and it needs to run fast so that the pulses 
            # are output fast enough for smooth motion
            analog_out.timing.cfg_samp_clk_timing(self.outputfreq,
                                                sample_mode=daq.AcquisitionType.FINITE,
                                                samps_per_chan=len(self.tout))    

            # set it to start when the analog input starts
            analog_out.triggers.start_trigger.cfg_dig_edge_start_trig("ai/StartTrigger",
                                                    trigger_edge=daq.Edge.RISING)

            # set up the digital output channel
            digital_out.do_channels.add_do_chan(motor_control_chan, 'motor')
            # use the analog output clock for digital output timing
            digital_out.timing.cfg_samp_clk_timing(self.outputfreq, 
                                                source = "ao/SampleClock",
                                                sample_mode=daq.AcquisitionType.FINITE,
                                                samps_per_chan=len(self.tout))
            digital_out.triggers.start_trigger.cfg_dig_edge_start_trig("ai/StartTrigger",
                                                    trigger_edge=daq.Edge.RISING)

            # set up to read the input
            reader = AnalogMultiChannelReader(analog_in.in_stream)
            self.aidata = np.zeros((len(self.inchannels), len(self.t)), dtype=np.float64)
            
            angle_reader = CounterReader(angle_in.in_stream)
            self.angledata = np.zeros((len(self.t),), dtype=np.float64)

            # write the output
            analog_writer = AnalogMultiChannelWriter(analog_out.out_stream, 
                                                    auto_start=False)
            analog_writer.write_many_sample(self.actcmdhi)
            
            digital_writer = DigitalSingleChannelWriter(digital_out.out_stream,
                                                        auto_start=False)
            nwritten = digital_writer.write_many_sample_port_uint32(self.dig)
            # print(f"{nwritten=}")

            # start everthing
            # make sure to start the output first, because it'll wait until the 
            # input starts
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
        return(self.aidata)

    
## TRYING SOMETHING OUT HERE...




# --- New Execution Logic Starts Here (Put this at the VERY bottom of bender_functions.py) ---

def run_experiment_configuration(test_type, **kwargs):
    """Sets up the specific parameters for a given test type using kwargs."""
    bender_instance = Bender()
    device_name = "Dev1"
    
    # ... (Keep all channel and calibration setup code here) ...
    bender_instance.set_input_channels(inchannels=['ai0', ...], inchannel_names=['Fx', ...])
    bender_instance.set_activation_channels('ao0', 'ao1')
    bender_instance.calibration = np.identity(6) 

    # General parameters
    duration = kwargs.get('duration', 2.0)
    sample_rate = kwargs.get('sample_rate', 100.0) 
    t = np.linspace(0, duration, int(duration * sample_rate)) 
    bender_instance.t = t # Set the time base for the Bender instance

    if test_type == 'static':
        # ... (Looping logic for static tests as in the previous answer) ...
        # NOTE: The looping logic needs to create a new bender instance inside the loop
        pass # Placeholder - keep the loop structure from the previous answer

    elif test_type == 'frequency':
        print(f"Configuring frequency sweep...")
        # Extract specific parameters from kwargs
        startfreq = kwargs['startfreq']
        endfreq = kwargs['endfreq']
        amplitude = kwargs['amplitude']
        exponent = kwargs['exponent']
        waitbefore = kwargs.get('waitbefore', 0.0) # You may need this parameter

        # Generate signals using the new helper method
        angle, anglevel, tnorm, freq = bender_instance._generate_frequency_sweep_signals(
            startfreq, endfreq, amplitude, exponent, duration, waitbefore
        )

    elif test_type == 'custom_cycle':
        print(f"Configuring custom cycle dynamic test...")
        # This will require loading or defining the cycle parameters (freq_by_cycle, etc.)
        # These lists need to be passed via command line or loaded from a config file.
        # For this example, we assume they are passed as lists via the command line args
        period_by_cycle = kwargs['periods']
        freq_by_cycle = kwargs['frequencies']
        amp_by_cycle = kwargs['amplitudes']

        angle, anglevel, tnorm, freq = bender_instance._generate_custom_cycle_signals(
            period_by_cycle, freq_by_cycle, amp_by_cycle
        )
        
    else:
        raise ValueError(f"Unknown test type: {test_type}")

    # Set the generated signals for the Bender instance to use in .run()
    # (If using the static loop, this part needs to be inside that loop)
    bender_instance.set_bending_signal(t, angle, anglevel)
    bender_instance.make_motor_stepper_pulses(outsampfreq=500)
    
    filename = bender_instance.increment_file_name(f'experiment_data_{test_type}_000.h5')
    print(f"Data will be saved to: {filename}")
    
    # bender_instance.run(device_name="Dev1", operation_type=test_type)
    print("Run command is currently commented out for testing signal generation.")


#
  def _generate_frequency_sweep_signals(self, startfreq, endfreq, amplitude, amplitude_frequency_exponent, duration, waitbefore):
        """Generates the log sweep angle and anglevel signals based on input params."""
        t = self.t # Use the time base already set by set_bending_signal
        samplefreq = self.samplefreq

        lnk = 1.0/duration * (np.log(endfreq) - np.log(startfreq))

        freq = startfreq * np.exp(t * lnk)

        tnorm = 2*np.pi*startfreq * (np.exp(t * lnk) - 1) / lnk

        tnorm[t < 0] = -1
        tnorm[t > duration] = np.ceil(np.max(tnorm))

        A0 = startfreq ** amplitude_frequency_exponent

        angle = amplitude / A0 * np.power(freq, amplitude_frequency_exponent) * np.sin(tnorm)
        anglevel = amplitude / A0 * np.exp(amplitude_frequency_exponent * t * lnk) * lnk * \
            (amplitude_frequency_exponent * np.sin(tnorm) + 2*np.pi/lnk * freq * np.cos(tnorm))

        freq[t < 0] = np.nan
        freq[t > duration] = np.nan

        angle[t < 0] = 0
        angle[t > duration] = 0

        anglevel[t < 0] = 0
        anglevel[t > duration] = 0

        isramp = (t >= duration) & (t < duration+0.5)
        # Ensure k calculation handles array indexing correctly
        k_index = np.argmax(t >= (waitbefore + duration))
        
        pend = angle[k_index]
        velend = (0 - pend) / 0.5

        ramp = pend + (t[isramp] - t[k_index])*velend

        np.place(anglevel, isramp, velend)
        np.place(angle, isramp, ramp)

        return angle, anglevel, tnorm, freq

    def _generate_custom_cycle_signals(self, period_by_cycle, freq_by_cycle, amp_by_cycle):
        """Generates signals for dynamic tests with a custom sequence of frequencies and amplitudes."""
        t = self.t
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

            np.place(tnorm, iscyc, (t[iscyc] - cycstart1) * f1 + c)
        
        # Calculate angle and anglevel from tnorm, freq, and amp here 
        # (Assuming you have this logic missing from your snippet, you need an actual motion calculation)
        # Placeholder calculation for demonstration:
        angle = amp * np.sin(2 * np.pi * tnorm) 
        anglevel = np.zeros_like(t) # This needs correct differentiation

        return angle, anglevel, tnorm, freq

# The main execution block starts here, at the bottom of the file
if __name__ == "__main__":
    # The argparse logic and the try/except block to run the program
    import argparse
    # ... (rest of the argparse code from the previous answer) ...
    try:
        run_experiment_configuration(args.test_type, **params)
    except Exception as e:
        logging.error(f"An error occurred: {e}")