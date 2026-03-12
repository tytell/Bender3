import numpy as np
from scipy import interpolate
from datetime import datetime
from copy import copy
import time
import re
import os
print(f"DEBUG: Loading bender_functions.py from: {os.path.abspath(__file__)}")
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
        self.S1stimcmd = None
        self.S2stimcmd = None
#Sets up the specific parameters for a given test type using data stored in self and general parameters passed via kwargs.

    def run_experiment(self, device_name = "Dev1", test_type = "dynamic"):

        # Replace ... with actual channel names as needed (Operate on self)
        self.set_input_channels(inchannels=['ai0'], inchannel_names=['Fx'])
        self.set_stim_channels('ao0', 'ao1')
        # self.calibration is assumed to be loaded in the notebook using bender.loadCalibration(...)

        # General parameters (duration/sample_rate passed from notebook)
        sample_rate = self.samplefreq

        # Define angle/anglevel/tnorm variables for scope consistency
        angle, anglevel, tnorm = None, None, None
        duration = None #Initialize, but will be set below

         # THIS LEVEL IS ABOUT CREATING MOTOR ANGLES
        if test_type == 'dynamic':            
            # Add safety check before using the attribute
            if self.period_by_cycle is None:
                raise AttributeError("Dynamic test requires 'period_by_cycle' to be set via organize_cycles first.")
            # Calculate duration. Needed for stimulation timing
            duration = np.sum(self.period_by_cycle)

            # Generate signals using the method on 'self', passing required inputs from self.attributes
            angle, anglevel, tnorm, _ = self.make_dynamic_cycles(
                self.period_by_cycle, 
                self.freq_by_cycle, 
                self.amp_by_cycle)
            
            # CHECK INDENT AND STRUCTURE HERE
        elif test_type == 'sweep':
            # Check a SWEEP attribute instead
            if self.duration is None:
                raise AttributeError("Sweep test requires 'self.duration' to be set in the notebook first.")
            duration = self.duration + self.waitbefore + self.waitafter

            angle, anglevel, tnorm, _ = self.make_frequency_sweep(
                self.all_freqs, 
                self.all_curves, 
                self.amplitude_frequency_exponent, 
                self.waitbefore)
        else:
            raise ValueError(f"Unknown test type: {test_type}")
    
        # NOW ACTUALLY MAKE THE SIGNAL FOR THE MOTOR (BASED ON PREVIOUS DYNAMIC, SWEEP, ETC)
        # Set the generated signals for     the Bender instance ('self') to use in .run()
        self.set_bending_signal(self.t, angle, anglevel, tnorm) # Pass tnorm here as well


        # NOW LET'S ADD ELECTRICAL STIMULI!
        # Access stim parameters stored in 'self' (set by organize_cycles)
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

        # FINISH Generate electrical stimuli
        S1stimcmd, S2stimcmd = self._make_stimuli(**stimulation_params)
            
        # Set the generated stim signals in 'self'
        self.set_stim(S1stimcmd, S2stimcmd)
        self.make_motor_stepper_pulses(outputfreq=self.outputfreq)
        filename = self.increment_file_name(f'experiment_data_{test_type}_000.h5')
        print(f"Data will be saved to: {filename}")
            
        # Run the experiment using 'self'
        self.run(device_name=device_name)



    def organize_cycles(self, all_curves, all_freqs, randomize, cycles_per_step, n_end_cycles, dclamp, xsec_width, stim_cycles_in_step, all_stimduties, all_stimphases, stim_pulse_rate):
        start = time.time()
        # Build combinations without modifying input lists
        combos = []
        for c1 in all_curves:
            for f1 in all_freqs:
                for d1 in all_stimduties:
                    for p1 in all_stimphases:
                        combos.append((c1, f1, d1, p1))
        # Unpack combos into separate arrays
        all_curves_arr, all_freqs_arr, all_stimduties_arr, all_stimphases_arr = map(np.array, zip(*combos))

        # Randomize if needed
        if randomize:
            order = np.arange(len(all_freqs_arr))
            np.random.shuffle(order)
            all_curves_arr = all_curves_arr[order]
            all_freqs_arr = all_freqs_arr[order]
            all_stimduties_arr = all_stimduties_arr[order]
            all_stimphases_arr = all_stimphases_arr[order]

        # Calculate amplitudes, strains, and strain rates
        all_degs = np.rad2deg(all_curves_arr * (dclamp/1000))
        allstrains = xsec_width/2/1000 * all_curves_arr
        allstrainrates = 2*np.pi * allstrains * all_freqs_arr

        # Create frequency, amplitude, and period arrays by cycle
        freq_by_cycle = np.repeat(all_freqs_arr, cycles_per_step)
        amp_by_cycle = np.repeat(all_degs, cycles_per_step)
        duty_by_cycle = np.repeat(all_stimduties_arr, cycles_per_step)
        phase_by_cycle = np.repeat(all_stimphases_arr, cycles_per_step)

        # Add end cycles
        freq_by_cycle = np.concatenate((freq_by_cycle, [all_freqs_arr[-1]] * n_end_cycles))
        amp_by_cycle = np.concatenate((amp_by_cycle, [all_degs[-1]] * n_end_cycles))
        duty_by_cycle = np.concatenate((duty_by_cycle, [all_stimduties_arr[-1]] * n_end_cycles))
        phase_by_cycle = np.concatenate((phase_by_cycle, [all_stimphases_arr[-1]] * n_end_cycles))

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

        # Store results
        self.period_by_cycle = period_by_cycle
        self.freq_by_cycle = freq_by_cycle
        self.amp_by_cycle = amp_by_cycle
        self.duty_by_cycle = duty_by_cycle
        self.phase_by_cycle = phase_by_cycle
        self.all_freqs = all_freqs_arr
        self.all_curves = all_curves_arr
        self.all_degs = all_degs
        self.stimburstdur = stimburstdur
        self.allstrains = allstrains
        self.allstrainrates = allstrainrates
        print("organize_cycles took", time.time() - start, "seconds")
# ...existing code...
    

    def set_bending_signal(self, t, angle, anglevel, tnorm=None):
        self.samplefreq = 1.0 / (t[1] - t[0])

        self.t = t
        self.angle = angle
        self.anglevel = anglevel

        if tnorm is None:
            self.tnorm = copy(t)
        else:
            self.tnorm = tnorm
    
    def set_stim(self, S1stimcmd, S2stimcmd):
        self.S1stimcmd = S1stimcmd
        self.S2stimcmd = S2stimcmd

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

    def make_motor_stepper_pulses(self, outputfreq = 1000,
                                scale=6.0,
                                stepsperrev=6400.0):

        self.outputfreq = outputfreq
        self.scale = scale

        tout = np.arange(self.t[0], self.t[-1], 1.0/outputfreq)

        poshi = interpolate.interp1d(self.t, self.angle, kind='linear', assume_sorted=True, bounds_error=False,
                                    fill_value=0.0)(tout)
        velhi = interpolate.interp1d(self.t, self.anglevel, kind='linear', assume_sorted=True, bounds_error=False,
                                    fill_value=0.0)(tout)

        poshi *= scale
        velhi *= scale


        if outputfreq == 0 or stepsperrev == 0:
            raise ValueError('Problems with parameters')

        stepsize = 360.0 / stepsperrev
        maxspeed = stepsize * outputfreq / 2

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

        if self.S1stimcmd is None:
            self.stimcmdhi = np.zeros((2, len(tout)))
        else:
            S1stimcmdhi = interpolate.interp1d(self.t, self.S1stimcmd, kind='linear', assume_sorted=True, bounds_error=False,
                                        fill_value=0.0)(tout)
            S2stimcmdhi = interpolate.interp1d(self.t, self.S2stimcmd, kind='linear', assume_sorted=True, bounds_error=False,
                                        fill_value=0.0)(tout)
            self.stimcmdhi = np.vstack((S1stimcmdhi, S2stimcmdhi))

        return tout, dig, motorstep, motordirection

    def set_input_channels(self, inchannels, inchannel_names):
        self.inchannels = inchannels
        self.inchannel_names = inchannel_names

    def set_stim_channels(self, S1stim_chan, S2stim_chan):
        self.S1stim_chan = S1stim_chan
        self.S2stim_chan = S2stim_chan

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
        S1stim_chan = '/'.join((device_name, self.S1stim_chan))
        S2stim_chan = '/'.join((device_name, self.S2stim_chan))
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
            analog_out.ao_channels.add_ao_voltage_chan(S1stim_chan, 'S1stim')
            analog_out.ao_channels.add_ao_voltage_chan(S2stim_chan, 'S2stim')
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
            analog_writer.write_many_sample(self.stimcmdhi)
            
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
# --- New Execution Logic Starts Here---








#
    def make_frequency_sweep(self, all_freqs, all_curves, amplitude_frequency_exponent, duration, waitbefore):
        # Generates the log sweep angle and anglevel signals based on input params
        samplefreq = self.samplefreq
        
        # Define start and end frequencies for the sweep based on all_freqs (first and last values)
        startfreq = all_freqs[0]
        endfreq = all_freqs[-1]


        # Calculate the total duration needed for the sweep + wait times
        total_duration = duration + waitbefore + self.waitafter
        t = np.arange(0, total_duration, 1.0 / samplefreq)
        t -= waitbefore # Shift time so the movement starts at t=0

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

    # Assuming this method exists inside a class definition
    def make_dynamic_cycles(self, period_by_cycle, freq_by_cycle, amp_by_cycle):
        """
        Generates signals for dynamic tests with a custom sequence of frequencies and amplitudes.
        
        Includes amplitude ramps during step changes and start/end ramps.
        
        Args:
            period_by_cycle (list/array): Duration of each cycle in seconds.
            freq_by_cycle (list/array): Frequency for each cycle (Hz).
            amp_by_cycle (list/array): Amplitude for each cycle.

        Returns:
            tuple: (angle, anglevel, tnorm, freq) NumPy arrays.
        """
        
        # Access parameters from self (ensure they are set by the __init__ or another method)
        waitbefore = self.waitbefore
        waitafter = self.waitafter
        rampdur = self.rampdur
        amp_step_vel = self.amp_step_vel
        samplefreq = self.samplefreq # <--- Use the value from the instance
        all_degs = self.all_degs # Used for start/end ramps
        all_freqs = self.all_freqs # Used for start/end ramps

        # Calculate timings and durations
        movedur = np.sum(period_by_cycle)
        totaldur = waitbefore + movedur + waitafter
        t = np.arange(0, totaldur, 1.0 / samplefreq)
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
        anglevel[1:-1] = (angle[2:] - angle[:-2]) * (samplefreq / 2.0)
        
        self.t = t
        self.tnorm = tnorm

        # Return all generated signals
        return angle, anglevel, tnorm, freq

    # Lots of unused arguments, need to check on that later.        
    def _make_stimuli(self, 
                      is_stim=None, 
                      phase_by_cycle=None, 
                      stim_pulse_rate=None, 
                      prestim_time=None, 
                      poststim_time=None, 
                      prepoststim_dur=None, 
                      prepoststim_sep=None, 
                      stimburstdur=None, 
                      duty_by_cycle=None, 
                      freq_by_cycle=None, 
                      movedur=None):

 # --- AUTO-FILL INGREDIENTS ---
        # This checks if the user provided an input. If not, it pulls from 'self'.
        is_stim         = is_stim         if is_stim is not None         else self.is_stim
        phase_by_cycle  = phase_by_cycle  if phase_by_cycle is not None  else self.phase_by_cycle
        stim_pulse_rate = stim_pulse_rate if stim_pulse_rate is not None else self.stim_pulse_rate
        prestim_time    = prestim_time    if prestim_time is not None    else self.prestim_time
        poststim_time   = poststim_time   if poststim_time is not None   else self.poststim_time
        prepoststim_dur = prepoststim_dur if prepoststim_dur is not None else self.prepoststim_dur
        prepoststim_sep = prepoststim_sep if prepoststim_sep is not None else self.prepoststim_sep
        stimburstdur    = stimburstdur    if stimburstdur is not None    else self.stimburstdur
        duty_by_cycle   = duty_by_cycle   if duty_by_cycle is not None   else self.duty_by_cycle
        freq_by_cycle   = freq_by_cycle   if freq_by_cycle is not None   else self.freq_by_cycle
        
        # Logic for movedur: use input, or calculate from periods, or use self.duration
        if movedur is None:
            movedur = np.sum(self.period_by_cycle) if hasattr(self, 'period_by_cycle') else self.duration
            
        t = self.t
        tnorm = self.tnorm
        S1stimcmd = np.zeros_like(t)
        S2stimcmd = np.zeros_like(t)
        Lonoff = []
        Ronoff = []
        bendphase = np.zeros_like(tnorm)
        prepostburst = (np.mod(t * stim_pulse_rate, 1) <= 0.5).astype(float)
        prepostburst *= 5.0

        if is_stim:
            pulsedur = 0.01         
            burst = (np.mod(t * stim_pulse_rate, 1) <= 0.5).astype(float)
            burst *= 5.0
            prepostburst = (np.mod(t * stim_pulse_rate, 1) <= 0.5).astype(float)
            prepostburst *= 5.0
            bendphase = tnorm - 0.25

        for c, (dur1, duty1, f1, p1) in enumerate(zip(stimburstdur, duty_by_cycle, freq_by_cycle, phase_by_cycle)):
            if dur1 == 0:
                continue
            
            # We add the safety check recommended previously
            target_time = c + p1
            if not np.any(bendphase >= target_time):
                continue # Skip cycles that fall outside the time array bounds

            k = np.argmax(bendphase >= target_time)
            tstart = t[k]
            tend = tstart + dur1

            if np.any(bendphase >= target_time):
                Lonoff.append([tstart, tend])
            if np.any(bendphase >= c + 0.5 + p1):
                Ronoff.append(np.array([tstart, tend]) + 0.5 / f1)

            # FIX 1: Replace np.place with boolean assignment for S1stimcmd
            mask1 = (bendphase >= c + p1) & (bendphase < c + p1 + duty1)
            S1stimcmd[mask1] = burst[mask1] # Assign only where the mask is True

            # FIX 2: Replace np.place with boolean assignment for S2stimcmd
            mask2 = (bendphase >= c + 0.5 + p1) & (bendphase < c + 0.5 + p1 + duty1)
            S2stimcmd[mask2] = burst[mask2] # Assign only where the mask is True


        # --- Pre/Post stimulation logic ---
        tstart = prestim_time
        tend = tstart + prepoststim_dur
        Lonoff.append([tstart, tend])
        
        # FIX 3: Replace np.place with boolean assignment for pre/post stim
        mask3 = (t >= tstart) & (t < tend)
        S1stimcmd[mask3] = prepostburst[mask3] # Assign only where the mask is True

        # Store these in the instance if you need to access them later for analysis
        self.Lonoff = np.array(Lonoff)
        self.Ronoff = np.array(Ronoff)

        return S1stimcmd, S2stimcmd
