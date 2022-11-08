import numpy as np
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
    from nidaqmx.stream_writers import AnalogMultiChannelWriter
    from nidaqmx.stream_readers import AnalogMultiChannelReader
    from nidaqmx.errors import DaqError
except ImportError:
    logging.warning('No DAQmx available')
    
import xml.etree.ElementTree as ElementTree

import h5py

class Ergometer:
    # TODO: Check the scales
    def __init__(self, lengthscale=0.1, # mm/V
                       forcescale=0.1, # N/V
                       stimvoltscale=20,  # output volt stimulus per input volt
                       stimampscale=0.1,  # GUESS! mA per V
                       stimtype = 'voltage'   # or current
                       ):
        self.lengthscale = lengthscale
        self.forcescale = forcescale
        self.stimvoltscale = stimvoltscale
        self.stimampscale = stimampscale
        self.stimtype = stimtype

    def set_length_signal(self, t, length, tnorm=None):
        self.samplefreq = 1.0 / (t[1] - t[0])

        self.t = t
        self.length = length / self.lengthscale

        if tnorm is None:
            self.tnorm = copy(t)
        else:
            self.tnorm = tnorm
    
    def set_activation(self, actcmd):
        if self.stimtype == 'voltage':
            self.actcmd = actcmd #/ self.stimvoltscale
        elif self.stimtype == 'current':
            self.actcmd = actcmd #/ self.stimampscale

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

    def run(self, device_name, lengthinchan, forceinchan, stimvoltinchan, stimampinchan,
                    lengthoutchan, stimvoltoutchan):
        inchannels = ['/'.join((device_name, c1)) for c1 in [lengthinchan, forceinchan, stimvoltinchan, stimampinchan]]
        inchannelnames = ['length', 'force', 'stimvolts', 'stimamps']
        outchannels = ['/'.join((device_name, c1)) for c1 in [lengthoutchan, stimvoltoutchan]]
        outhchannelnames = ['length', 'stimulus']

        with Task() as analog_in, Task() as analog_out:
            # set up the input channels
            for c1, name1 in zip(inchannels, inchannelnames):
                analog_in.ai_channels.add_ai_voltage_chan(c1, name1)

            # set up the input sample frequency
            # just records as many samples as are in the output
            analog_in.timing.cfg_samp_clk_timing(self.samplefreq,
                                                sample_mode=daq.AcquisitionType.FINITE,
                                                samps_per_chan=len(self.t))

            # set up the analog output channels
            for c1, name1 in zip(outchannels, outhchannelnames):
                analog_out.ao_channels.add_ao_voltage_chan(c1, name1)

            # it will run much faster than the input channels, because the digital output is linked
            # to it, and it needs to run fast so that the pulses 
            # are output fast enough for smooth motion
            analog_out.timing.cfg_samp_clk_timing(source='/'.join((device_name, 'ai', 'SampleClock')),
                                                rate = self.samplefreq,
                                                sample_mode=daq.AcquisitionType.FINITE,
                                                samps_per_chan=len(self.t))    

            # set it to start when the analog input starts
            # analog_out.triggers.start_trigger.cfg_dig_edge_start_trig("ai/StartTrigger",
            #                                         trigger_edge=daq.Edge.RISING)

            # set up to read the input
            reader = AnalogMultiChannelReader(analog_in.in_stream)
            self.aidata = np.zeros((len(inchannels), len(self.t)), dtype=np.float64)
            
            # write the output
            outdata = np.vstack((self.length, self.actcmd))
            analog_writer = AnalogMultiChannelWriter(analog_out.out_stream, 
                                                    auto_start=False)
            analog_writer.write_many_sample(outdata)
            
            # start everthing
            # make sure to start the output first, because it'll wait until the 
            # input starts
            analog_out.start()
            analog_in.start()
            
            # wait until we're done, record the time
            analog_in.wait_until_done(self.t[-1]+10)
            self.endTime = datetime.now()
            
            # and read the data
            reader.read_many_sample(self.aidata)
    
        return(self.aidata)

    def scale_data(self, aidata):
        length = aidata[0, :] * self.lengthscale
        force = aidata[1, :] * self.forcescale
        stimvolt = aidata[2, :] * self.stimvoltscale
        stimamp = aidata[3, :] * self.stimampscale

        return((length, force, stimvolt, stimamp))
    
