import numpy as np
from scipy import interpolate
from datetime import datetime
import time
from copy import copy

import re
import os

import nidaqmx.constants as daq
from nidaqmx import Task
from nidaqmx.stream_writers import AnalogMultiChannelWriter, DigitalSingleChannelWriter
from nidaqmx.stream_readers import AnalogMultiChannelReader, CounterReader
from nidaqmx.errors import DaqError

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
        m = re.search('(\d+)\.h5', filename)
        if m is None:
            basename, ext = os.path.splitext(filename)
            num = 0
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
    
