'''
growthrate.py
----

Class to obtain instability growth rates 
with optimized MWFFT and interactive start-end 
choosing

Date: 21/02/23
Author: edelafue
'''

import os
import glob
import time
import csv

import numpy as np
import h5py
import matplotlib.pyplot as plt
import scipy.signal
from scipy.constants import c

def onclick(event):
    global _grstart
    _grstart = event.xdata
    plt.close()

class GrowthRate():

    def __init__(self, file='bunchmonitor.h5', nperseg=65, noverlap=50, 
                 turn_min=1 , turn_max=None, interactive=True, grstart=None,
                 machine='SPS', grmin=1e-5, plane='V', pad=3, source='BQ.KICKED',
                 **kwargs):

        self.file = file
        self.nperseg = nperseg
        self.noverlap = noverlap
        self.turn_min = turn_min
        self.turn_max = turn_max
        self.interactive = interactive
        self.grstart = grstart
        self.grmin = grmin #gr limit to consider beam stable
        self.beam_unstable = True
        self.pad = pad
        self.plane = plane
        self.source = source

        self.qpv = None 
        self.qph = None 
        
        self._machine = machine
        if self._machine == 'LHC':
            self.BUCKET_MAX = 3564
            RING_CIRCUMFERENCE = 26658.883   #[m]
            GAMMA_R = 7461                   # flat top
                        
        elif self._machine =='SPS':
            self.BUCKET_MAX = 920 #924 is the correct number
            RING_CIRCUMFERENCE = 6895        #[m]
            #GAMMA_R = 251    # flat top 236 GeV
            GAMMA_R = 27.7   # flat bottom value 26 GeV

        elif self._machine == 'PS':
            self.BUCKET_MAX = 21
            RING_CIRCUMFERENCE = 628
            GAMMA_R = 27.7366 #28.7185  # p=26GeV

        BETA_R = np.sqrt(1 - (1/GAMMA_R**2))
        self.T_1_TURN = RING_CIRCUMFERENCE/(c*BETA_R)
        self.REV_FREQ = 1/self.T_1_TURN

        self.pos = None 
        self.turn = None 
        self.ms = None 

        if file.endswith('.h5'):
            self.get_data()

        elif file.endswith('.parquet'):
            self.get_data_parquet()
            if (self.qpv is None) or (self.qph is None): 
                self.get_QP_knob()

        if self.turn_max is None:
            self.turn_max = len(self.pos)
            
        for key, val in kwargs.items():
            setattr(self, key, val)

    def get_data(self):

        if self.plane == 'H':
            plane = 'x'
        else:
            plane = 'y'
            
        f = h5py.File(self.file)
        self.data = f['Bunch']
        self.pos = np.array(f['Bunch']['mean_'+plane])
        self.turns = np.arange(len(self.pos))
        self.ms = self.turns*self.T_1_TURN

    def get_data_parquet(self):
        import datascout as ds
        d = ds.parquet_to_dict(self.file)
        self.data = d

        # Define key for source of data
        if self.source == 'BQ.KICKED':
                key = 'SPS.BQ.KICKED/ContinuousAcquisition'
        elif self.source == 'BQ.QC':
                key = 'SPS.BQ.QC/Acquisition'
        # Read dict data
        try:
            self.pos = d[key]['value']['rawData'+self.plane]
        except:
            print(f'{self.file} contains no valid data')
            self.beam_unstable = False
            return
        self.turns = np.arange(len(self.pos))
        self.ms = self.turns*self.T_1_TURN

    def get_QP_knob(self, cycle_offset=4015, cycle_start=100):
        X = self.data['SPSBEAM/QPV']['value']['JAPC_FUNCTION']['X'] - cycle_offset
        Y = self.data['SPSBEAM/QPV']['value']['JAPC_FUNCTION']['Y']
        self.qpv = float(Y[X > cycle_start][0])

        X = self.data['SPSBEAM/QPH']['value']['JAPC_FUNCTION']['X'] - cycle_offset
        Y = self.data['SPSBEAM/QPH']['value']['JAPC_FUNCTION']['Y']
        self.qph = float(Y[X > cycle_start][0])

        return self.qph, self.qpv

    def get_grstart(self):

        fig, ax = plt.subplots(1,1, num=1)
        fig.tight_layout()
        fig.set_size_inches(6,6)
        ax.plot(self.turns, self.pos, color='b')

        if self.turn_min is not None:
            ax.set_xlim(xmin = self.turn_min)
        if self.turn_max is not None:    
            ax.set_xlim(xmax = self.turn_max)
        if self.turn_max is not None and self.turn_min is not None:
            ax.set_ylim(ymin=self.pos[int(self.turn_min):int(self.turn_max)].min(), ymax=self.pos[int(self.turn_min):int(self.turn_max)].max())
            

        ax.set_ylabel('Vertical position', color = 'b')
        ax.set_xlabel('Turns [-]')
        ax.tick_params(axis='y', colors='b')
        ax.set_title('Click on instability start', fontweight='bold', color='red')
        fig.tight_layout()

        #get instability start 
        fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show()
        self.grstart = _grstart

    def MWFFT_analysis(self):
        print('Performing MWFFT analysis')

        if self.pos is None:
            return 

        if self.grstart is None and self.interactive:
            self.get_grstart()

        elif self.grstart is None:
            raise Exception('`grstart` attribute not provided. Set `interactive` to `True`')

        xstart = np.where(self.turns >= self.grstart)[0][0]
        xpos = np.where(self.turns >= self.grstart)[0]
        xend = np.argmax(self.pos[xpos])+xstart+1
        
        if xend > self.turn_max: #to avoid dump signal
            print('\033[93m'+f'[!] Found max. at {xend} > turn_max...'+'\033[0m')
            pos = self.pos[xpos]
            xend = np.argmax(pos[:self.turn_max])+xstart+1
            
        self.beam_unstable = True

        self.xstart = xstart
        self.xpos = xpos
        self.xend = xend

        if (xend-xstart) < self.nperseg:
            self.slope = 0.
            self.error = 0.
            self.beam_unstable = False
            print('\033[91m'+'[!!] Beam is not unstable, slope set to 0.0'+'\033[0m')
            return

        #perform MWFFT (scipy.signal)
        pad = self.pad
        pos_pad = self.pos[xstart-pad:xend+pad]
        MWFFT_all = scipy.signal.stft(self.pos, fs=1.0, window='hann', nperseg=self.nperseg, noverlap=self.noverlap, nfft=None, detrend=False, return_onesided=True, boundary='zeros', padded=True, axis=-1)
        MWFFT = scipy.signal.stft(pos_pad, fs=1.0, window='hann', nperseg=self.nperseg, noverlap=self.noverlap, nfft=None, detrend=False, return_onesided=True, boundary='zeros', padded=True, axis=- 1)
        self.conversion = int(self.nperseg - self.noverlap)
        self.points = MWFFT_all[1]
        self.MWFFT = MWFFT_all
        
        #get the log value of the max amplitude per window
        log_peak_amplitude = []
        for i in range(len(MWFFT[1])):
            log_peak_amplitude.append(np.log(np.max(np.abs(MWFFT[2][:,i]))))

        log_peak_amplitude_all = []
        for i in range(len(self.points)):
            log_peak_amplitude_all.append(np.log(np.max(np.abs(MWFFT_all[2][:,i]))))

        self.gry = np.array(log_peak_amplitude)[pad:-pad]
        self.grx = MWFFT[1][pad:-pad]+xstart
        self.log_peak_amplitude_all = log_peak_amplitude_all

        # linear fit
        try:
            coef, cov = np.polyfit(self.grx, self.gry, 1, cov=True)
            self.poly1d_fn = np.poly1d(coef)
            self.error = np.sqrt(np.diag(cov)[0])
            self.slope = float(coef[0])
            
            if self.slope < 0:
                print('\033[91m'+f'[!!] gr.slope is negative due to a bad fit and will be set to 0.0'+'\033[0m')
                self.error = 0.
                self.slope = 0.
                self.beam_unstable = False
            elif self.slope < self.grmin:
                print('\033[91m'+f'[!!] gr.slope is too small and will be set to 0.0'+'\033[0m')
                self.error = 0.
                self.slope = 0.
                self.beam_unstable = False
            else:
                print(f'gr.slope = {self.slope*1e3}')

        except:
            print('\033[91m'+'[!!] Linear Fit could not be performed, slope set to 0.0'+'\033[0m')
            self.error = 0.
            self.slope = 0.
            self.beam_unstable = False

    def plot(self):

        xstart = self.xstart
        xend = self.xend

        fig, (ax1, ax2) = plt.subplots(2, 1, num=1)

        ax1.axvline(self.grstart,c='r')
        ax1.plot(self.turns, self.pos, color='b', alpha = 0.5)
        ax1.plot(self.turns[xstart:xend], self.pos[xstart:xend], color='b') 

        if self.turn_min is not None:
            ax1.set_xlim(xmin = self.turn_min)
        if self.turn_max is not None:    
            ax1.set_xlim(xmax = self.turn_max)

        ax1.set_ylim(-np.max(self.pos[xstart:xend])*1.2, np.max(self.pos[xstart:xend])*1.2)
        ax1.set_ylabel('Vertical position', color = 'b')
        ax1.set_xlabel('Turns [-]')
        ax1.tick_params(axis='y', colors='b')

        if self.qpv is not None and self.qph is not None:
            ax1.set_title(f'Growth rate for $Q\'_v$ = {round(self.qpv,2)}, $Q\'_h$ = {round(self.qph,2)}')

        if self.beam_unstable:

            ax2.plot(self.points, self.log_peak_amplitude_all,'.-b', label='MWFFT')
            ax2.plot(self.grx, self.gry, '.-r', label='points for fit')
            ax2.plot(self.grx, self.poly1d_fn(self.grx), '--k', lw=2.0, label='fit masked data')

            ax2.text(0.1,0.1," $\u03C4^{-1}$ = %.2E +/- %.E" %(self.slope, self.error),va='top',transform=ax2.transAxes, color='black',\
                     bbox=dict(boxstyle='round',facecolor='white',alpha=0.9), fontsize=11)
            ax2.set_xlabel('Turns [-]')
            ax2.set_ylabel('log(Amplitude of Peak in FFT)', color = 'red')
            ax2.tick_params(axis='y', colors='red')

            if self.turn_min is not None:
                ax2.set_xlim(xmin = self.turn_min)
            if self.turn_max is not None:    
                ax2.set_xlim(xmax = self.turn_max)

            ax2.legend()

        else:
            ax2.text(0.1,0.4,"Beam is not unstable",va='top',transform=ax2.transAxes, color='black',\
                     bbox=dict(boxstyle='round',facecolor='white',alpha=0.9), fontsize=11)


        fig.suptitle('Growth rate MWFFT analysis', fontweight='bold')
        fig.set_size_inches(6,10)
        fig.tight_layout()

        plt.show()

    def save(self, fname='growthrate.csv', header=None, **kwargs):

        # check for nan [TODO]

        # set header and rows
        for key, val in kwargs.items():
            setattr(self, key, val)

        if header is None:
            header = ['slope', 'error']
            row = [self.slope, self.error]

            for key in kwargs.keys():
                header.append(key)
                row.append(getattr(self, key))
        else:
            row = []
            for key in header:
                row.append(getattr(self, key))

        # write csv
        if not os.path.exists(fname):
            with open(fname, 'w') as f:
                writer = csv.writer(f, delimiter=' ')
                writer.writerow(header)

        with open(fname, 'a') as f:  
            writer=csv.writer(f, delimiter=' ')
            writer.writerow(row)


