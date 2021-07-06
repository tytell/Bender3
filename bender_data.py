import os
from collections import defaultdict
import numpy as np
import h5py
import re

import plotly.graph_objects as go
from plotly.subplots import make_subplots 

import pandas as pd

def signed_max(a):
    m1 = np.max(a)
    m2 = np.min(a)

    if abs(m2) > abs(m1):
        return m2
    else:
        return m1

class BenderData:
    def __init__(self, filename):
        self.filename = filename

        with h5py.File(filename, 'r') as h5file:
            self.yforce = np.array(h5file['/Calibrated/yForce'])
            self.xtorque = np.array(h5file['/Calibrated/xTorque'])
            self.ztorque = np.array(h5file['/Calibrated/zTorque'])
            self.t = np.array(h5file['/NominalStimulus/t'])
            self.tnorm = np.array(h5file['/NominalStimulus/tnorm'])
            self.angle = np.array(h5file['/Calibrated/Encoder'])
            self.angle_cmd = np.array(h5file['/NominalStimulus/Position'])

            if 'Frequencies' in h5file['/NominalStimulus'].attrs:
                self.frequencies = h5file['/NominalStimulus'].attrs['Frequencies']
                self.curvatures = h5file['/NominalStimulus'].attrs['Curvatures']
                self.amplitudes = h5file['/NominalStimulus'].attrs['Amplitudes']
            else:
                self.frequencies = [h5file['/NominalStimulus'].attrs['Frequencies']]
                self.curvatures = [h5file['/NominalStimulus'].attrs['Curvature']]
                self.amplitudes = [h5file['/NominalStimulus'].attrs['Amplitude']]

            if 'Cycles' in h5file['/NominalStimulus'].attrs:
                self.ncycles = h5file['/NominalStimulus'].attrs['Cycles']
                self.movedur = self.ncycles / self.frequencies[0]
                self.freqbycycle = [self.frequencies[0]] * self.ncycles
            else:
                self.freqbycycle = h5file['/NominalStimulus'].attrs['FrequencyByCycle']
                movedur = np.sum(1.0 / self.freqbycycle)
                self.ncycles = len(self.freqbycycle)
                self.movedur = movedur

            self.Lonoff = np.empty((0,0))
            self.Ronoff = np.empty((0,0))
            self.is_active = h5file['/NominalStimulus'].attrs['ActivationOn']      
            self.activation_duty = 0
            self.activation_phase = np.nan
            self.start_cycle = None      
            if self.is_active:
                self.activation_duty = h5file['/NominalStimulus'].attrs['ActivationDuty']
                if 'ActivationStartPhase' in h5file['/NominalStimulus'].attrs:
                    self.activation_phase = h5file['/NominalStimulus'].attrs['ActivationStartPhase']
                else:
                    self.activation_phase = h5file['/NominalStimulus'].attrs['ActivationPhase']

                if 'ActivationStartCycle' in h5file['/NominalStimulus'].attrs:
                    self.start_cycle = h5file['/NominalStimulus'].attrs['ActivationStartCycle']
                    self.active_cycles = np.full((self.ncycles,), True)
                    self.active_cycles[:self.start_cycle] = False
                else:
                    self.start_cycle = None
                    self.active_cycles = h5file['/NominalStimulus'].attrs['IsActiveByCycle']

                if 'Lonoff' in h5file['/NominalStimulus']:
                    self.Lonoff = np.array(h5file['/NominalStimulus/Lonoff'])
                    self.Ronoff = np.array(h5file['/NominalStimulus/Ronoff'])
                elif self.start_cycle is not None:
                    self.generate_activation()

            self.bodytorque = self.xtorque

            self.remove_bias()
            self.get_isometric_torques()

    def generate_activation(self):
        if not self.is_active:
            self.Lonoff = np.empty((0,0))
            self.Ronoff = np.empty((0,0))
            return
            
        self.Lonoff = []
        self.Ronoff = []

        bendphase = self.tnorm - 0.25

        # list of cycles when we'll have activation
        actcycles = list(range(self.start_cycle, self.ncycles))

        # also do one cycle before the bending starts
        actcycles.insert(0, -3)
        # and one after
        actcycles.append(self.ncycles+1)

        for c in actcycles:
            k = np.argmax(bendphase >= c + self.activation_phase)
            tstart = self.t[k]
            tend = tstart + self.activation_duty / self.freq

            if any(bendphase >= c + self.activation_phase):
                self.Lonoff.append([tstart, tend])
            if any(bendphase >= c + self.activation_phase + 0.5):
                self.Ronoff.append(np.array([tstart, tend]) + 0.5 / self.freq)

        self.Lonoff = np.array(self.Lonoff)
        self.Ronoff = np.array(self.Ronoff)    

    def remove_bias(self):
        sigs = [self.xtorque, self.ztorque, self.yforce]

        if len(self.Lonoff) > 0:
            pre = self.Lonoff[0,0] - 0.5
        else:
            pre = -0.5

        for sig1 in sigs:
            s0 = np.mean(sig1[self.t < pre])
            sig1 -= s0

    def get_isometric_torques(self):
        self.iso_torque = np.zeros((2,2))
        if self.Lonoff is not None and self.Lonoff.shape[0] >= 2:
            self.iso_torque[0,0] = signed_max(self.bodytorque[(self.t >= self.Lonoff[0,0]) & (self.t < self.Lonoff[0,1])])
            self.iso_torque[0,1] = signed_max(self.bodytorque[(self.t >= self.Ronoff[0,0]) & (self.t < self.Ronoff[0,1])])
            self.iso_torque[1,0] = signed_max(self.bodytorque[(self.t >= self.Lonoff[-1,0]) & (self.t < self.Lonoff[-1,1])])
            self.iso_torque[1,1] = signed_max(self.bodytorque[(self.t >= self.Ronoff[-1,0]) & (self.t < self.Ronoff[-1,1])])
    
    def labeller(self, fmt):
        if fmt is None:
            return None

        fn = os.path.basename(self.filename)

        labels = defaultdict(str, [('fn', fn),
                  ('ph', '{:.2f}'.format(self.activation_phase)),
                  ('dc', '{:.2f}'.format(self.activation_duty))])

        if len(self.curvatures) == 1:
            labels['c'] = '{:.1f}'.format(self.curvatures[0])
        else:
            labels['c'] = '{:.1f}-{:.1f}'.format(min(self.curvatures), max(self.curvatures))

        if len(self.frequencies) == 1:
            labels['f'] = '{:.2f}'.format(self.frequencies[0])
        else:
            labels['f'] = '{:.1f}-{:.1f}'.format(min(self.frequencies), max(self.frequencies))

        m = re.search('(\d+)\.h5', fn)
        if m is not None:
            labels['num'] = m.group(1)

        return(fmt.format_map(labels))
        
    def get_data(self):
        fn = os.path.basename(self.filename)

        m = re.search('(\d+)\.h5', fn)
        if m is not None:
            trial = m.group(1)
            trial = int(trial)
            i = trial
        else:
            trial = None
            i = 0
        
        df = pd.DataFrame({'filename': fn,
                           'trial': trial,
                           'mincurvature': min(self.curvatures),
                           'maxcurvature': max(self.curvatures),
                           'minfreq': min(self.frequencies),
                           'maxfreq': max(self.frequencies),
                           'minamplitude': min(self.amplitudes),
                           'maxamplitude': max(self.amplitudes),
                           'is_active': self.is_active,
                           'act_phase': self.activation_phase,
                           'act_duty': self.activation_duty,
                           'Ltorque_before': self.iso_torque[0,0],
                           'Rtorque_before': self.iso_torque[0,1],
                           'Ltorque_after': self.iso_torque[1,0],
                           'Rtorque_after': self.iso_torque[1,1]}, index=[i])
        
        return(df)
        
    def plot_time_series(self, label=None):
        fig = make_subplots(specs=[[{"secondary_y": True}]])

        fig.add_trace(
            go.Scatter(x = self.t, y = self.angle, mode="lines", name="angle"),
            secondary_y=True)
        fig.add_trace(
            go.Scatter(x = self.t, y = self.angle_cmd, mode="lines", name="angle_cmd"),
            secondary_y=True)
        fig.add_trace(
            go.Scatter(x = self.t, y = self.xtorque, mode="lines", name="stim"),
            secondary_y=False)

        for onoff in self.Lonoff:
            fig.add_vrect(x0 = onoff[0], x1=onoff[1], fillcolor="black", opacity=0.25, line_width=0,
                                row="all", col="all")

        for onoff in self.Ronoff:
            fig.add_vrect(x0 = onoff[0], x1=onoff[1], opacity=0.7, line_width=1,
                            row="all", col="all")

        angrng = np.max(np.abs(self.angle))
        torquerng = np.max(np.abs(self.xtorque))

        fig.update_yaxes(title_text = "angle (deg)", secondary_y=True,
                range=[-angrng, angrng], showgrid=False)
        fig.update_yaxes(title_text = "torque (Nm)", secondary_y=False,
                range=[-torquerng, torquerng])

        if len(self.Lonoff > 0):
            rng = [self.Lonoff[0,0]-0.5, self.Ronoff[-1,1]+0.5]
        else:
            rng = [-0.5, self.movedur + 0.5]

        fig.update_xaxes(title_text = "time (s)",
                        range = rng)
        if label is not None:
            title = self.labeller(label)
            fig.update_layout(title=title)

        return(fig)

    def get_is_active_passive(self):
        periods = 1.0 / self.freqbycycle
        tcycle = np.cumsum(periods)
        tcycle = np.insert(tcycle, 0,0)
        
        is_active = np.full_like(self.t, False, dtype='bool')
        is_passive = np.full_like(self.t, False, dtype='bool')
        for tc1,ac1,p1 in zip(tcycle, self.active_cycles, periods):
            iscyc = (self.t >= tc1) & (self.t < tc1+p1)
            if ac1:
                is_active[iscyc] = True
                is_passive[iscyc] = False
            else:
                is_active[iscyc] = False
                is_passive[iscyc] = True
        return is_active, is_passive

    def plot_passive_and_active_loop(self, fig=None, row=1, title=None, 
                        passive_name="passive", active_name="active"):
        if fig is None:
            fig = make_subplots(rows=1, cols=2, shared_xaxes=True, shared_yaxes=True)

        is_active, is_passive = self.get_is_active_passive()

        # tpassive = [0.25/self.freq, self.Lonoff[1,0]]
        # tactive = [self.Lonoff[1,0], (self.ncycles-0.25) / self.freq]

        # is_passive = (self.t >= tpassive[0]) & (self.t < tpassive[1])
        # is_active = (self.t >= tactive[0]) & (self.t < tactive[1])

        passive_name = self.labeller(passive_name)
        active_name = self.labeller(active_name)
        
        fig.add_trace(go.Scatter(x = self.angle[is_passive], y = self.xtorque[is_passive],
                                name=passive_name, showlegend=(passive_name is not None)),
                      row = row, col = 1)
        fig.add_trace(go.Scatter(x = self.angle[is_active], y = self.xtorque[is_active],
                                name=active_name, showlegend=(active_name is not None)),
                      row = row, col = 2)
        fig.add_hline(y = self.iso_torque[0,0])
        fig.add_hline(y = self.iso_torque[0,1])

        title = self.labeller(title)
        fig.update_layout(title=title)

        return(fig)
                      
    def plot_active_loop(self, fig=None, index=1, rows=1, cols=1, title=None, 
                        passive_name="passive", active_name="active",
                        show_pre_post=True):
        if fig is None:
            fig = make_subplots(rows=rows, cols=cols, shared_xaxes=True, shared_yaxes=True)

        r = int(np.floor(index / cols))
        c = int(index - (r*cols))

        print(f"{r=}, {c=}")

        is_active, _ = self.get_is_active_passive()
        
        active_name = self.labeller(active_name)
        
        fig.add_trace(go.Scatter(x = self.angle[is_active], y = self.xtorque[is_active],
                                name=active_name, showlegend=(active_name is not None)),
                      row = r+1, col = c+1)
        if show_pre_post:
            fig.add_hline(y = self.iso_torque[0,0])
            fig.add_hline(y = self.iso_torque[0,1])

        title = self.labeller(title)
        fig.update_layout(title=title)

        return(fig)

