'''
Chromaticity estimation from
SPS BQHT headtail monitor data

`.h5` files are stored in:

/nfs/cs-ccr-bqhtnfs/sps_data/SPS.BQHT/yyyy_mm_dd

August'24
'''

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import csv
import os

class HTchromaticity():

    def __init__(self, file=None, nturns=300, plane='V', knob=None, parquet=None,
                 bl=3, monitor=False, save_png=False, dpi=100, **kwargs):

        self.file = file
        self.nturns = nturns
        self.plane = plane
        self.bl = bl #TODO use measurement 
        self.knob = knob
        self.parquet = parquet
        self.monitor = monitor
        self.dpi = dpi
        self.save_png = save_png

        # plotting
        try:
            import scienceplots
            plt.style.use(['science','no-latex', 'ieee'])
        except:
            pass
        matplotlib.rcParams.update({'font.size': 10})

        # monitor
        if self.monitor:
            self.fig, self.axs = plt.subplots(1,2, dpi=self.dpi, figsize=(6.25, 3.25), num=1)

        # read data
        self.delta = None
        self.len = None
        self.xp, self.yp = None, None 
        self.qp = None

        if type(self.knob) is str:
            self.knob = float(self.knob)

        if self.parquet is not None:
            self.get_knob_from_parquet()

        self.get_data()

    def get_knob_from_parquet(self, cycle_offset=5515, cycle_start=100):
        import datascout as ds
        d = ds.parquet_to_dict(self.parquet)
        try:
            self.pos = d['SPS.BQ.KICKED/ContinuousAcquisition']['value']['rawData'+self.plane]
        except:
            print(f'{self.file} contains no valid data')
            self.beam_unstable = False
            return

        X = d['SPSBEAM/QPV']['value']['JAPC_FUNCTION']['X'] - cycle_offset
        Y = d['SPSBEAM/QPV']['value']['JAPC_FUNCTION']['Y']
        self.qpv = float(Y[X > cycle_start][0])

        X = d['SPSBEAM/QPH']['value']['JAPC_FUNCTION']['X'] - cycle_offset
        Y = d['SPSBEAM/QPH']['value']['JAPC_FUNCTION']['Y']
        self.qph = float(Y[X > cycle_start][0])

        if self.plane == 'V':
            self.knob = self.qpv
        elif self.plane == 'H':
            self.knob = self.qph

    def get_data(self):

        f = h5py.File(self.file)

        # Read and reshape
        if self.plane == 'H':
            self.delta = f['horizontal']['delta']
        elif self.plane == 'V':
            self.delta = f['vertical']['delta']
        else:
            raise Exception('Specify a valid plane: "H" or "V"')

        self.len = int(len(self.delta)/self.nturns) #should be 250
        self.delta = np.reshape(self.delta, [self.nturns, self.len])

        # Remove baseline
        for i in range(self.nturns):
            self.delta[i, :] = self.delta[i, :] - np.mean(self.delta, axis=0) 

    def calc_chromaticity(self):
        
        # 2D fourier transform
        ft = np.fft.ifftshift(self.delta)
        ft = np.fft.fft2(self.delta)
        ft = np.fft.fftshift(ft)
        self.ft = np.abs(ft)

        # Compute phase offset
        ftleft = np.abs(self.ft[:(self.delta.shape[0]//2-5), :])
        ftyfreq = np.fft.fftshift(np.fft.fftfreq(self.len, 25e-9/self.len))
        self.xp, self.yp = np.unravel_index(np.argmax(ftleft), ftleft.shape)
        
        # normalize to chromaticity
        self.qp_freq = ftyfreq[self.yp]
        frev = 43.4e3 # 26 GeV [Hz]
        #slip_factor = 0.62e-3    # 26GeV q26
        slip_factor = 1.6e-3      # 100 GeV q26
        tune = 26.18
        self.qp = -self.qp_freq*slip_factor/frev/tune

        # old dirty way
        #phase = -(self.yp - self.ft.shape[1]//2) / (self.ft.shape[1]//2)
        #self.qp = phase*3       

    def plot(self, xlim=None, ylim=None):
        '''
        Params
        ------
        xlim: list, opt
            x limits of the HT data in turns
        ylim: list, opt
            y limits of the HT data in time [ns]
        '''

        if self.monitor:
            for ax in self.axs:
                ax.cla()
        else:
            self.fig, self.axs = plt.subplots(1,2, dpi=self.dpi, figsize=(6.25, 3.25), num=1)

        t = np.linspace(0, 25, self.len) #s --> 25 ns discussed with Tom
        n = np.arange(self.nturns) #s
        fx = np.fft.fftfreq(len(n), d=1/43e3) #sps rev freq
        fy = np.fft.fftfreq(len(t), d=(t[2]-t[1]))

        # z-t HT signal
        extent_t = [0, self.nturns, t[0], t[-1]]
        ht = self.axs[0].imshow(self.delta.transpose(), cmap='rainbow', origin='lower', interpolation='spline36', extent=extent_t, aspect='auto')
        
        # 2D FFT signal
        extent_f = [-fx.max()/1e3, fx.max()/1e3, -fy.max()/1e9, fy.max()/1e9]
        ft = self.axs[1].imshow(self.ft.transpose(), cmap='rainbow', origin='lower', interpolation='spline36', extent=extent_f, aspect='auto')

        self.axs[0].set_title(f'HT monitor - Bunch $\Delta_{self.plane}$')
        self.axs[0].set_ylabel('Time [ns]')
        self.axs[0].set_xlabel('Turns')
        self.axs[0].set_ylim(12.5-self.bl,12.5+self.bl) #in [ns]
        if ylim is not None:
            self.axs[0].set_ylim(ylim[0],ylim[1])
        if xlim is not None:
            self.axs[0].set_xlim(xlim[0],xlim[1])
        
        self.axs[0].text(0.65, 0.95, f'Calc: $Q\'_{self.plane}={round(self.qp,2)}$', ha='center',va='center', transform=self.axs[0].transAxes)

        self.axs[1].set_title('2D FFT Magnitude')
        self.axs[1].set_yticklabels([])
        self.axs[1].set_xticklabels([])
        #axs[1].set_ylabel('Frequency [GHz]')
        #axs[1].set_xlabel('Frequency [kHz]')
        #axs[1].set_ylim(-2,2)
        #axs[1].set_xlim(-10,10)

        self.fig.colorbar(ht, ax=self.axs[0], label=f'$\Delta_{self.plane}$ amplitude')
        self.fig.colorbar(ft, ax=self.axs[1], label='FT amplitude')

        if self.knob is not None:
            fname = self.file.split('/')[-1]
            self.fig.suptitle(f'{fname} | QP{self.plane} knob={str(self.knob)} | QP{self.plane} calc={round(self.qp,2)}', fontsize=6)
        else:
            fname = self.file.split('/')[-1]
            self.fig.suptitle(f'{fname} | QP{self.plane} calc={round(self.qp,2)}', fontsize=6)

        self.fig.tight_layout()

        if self.save_png:
            self.fig.savefig(f'{self.file}.png')

        if self.monitor:
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
        else:
            plt.show()

    def save(self, fname='HTchromaticity.csv', **kwargs):
        # **kwargs = {'key-to-save' : val-to-save}

        fname = self.file.split('/')[-1]
        
        # set header and rows
        for key, val in kwargs.items():
            setattr(self, key, val)

        if self.knob is not None:
            header = ['file', 'plane', 'QP knob', 'QP calc', 'error']
            row = [fname, self.plane, self.knob, self.qp, np.abs(self.qp-self.knob)/self.knob ]

            for key in kwargs.keys():
                header.append(key)
                row.append(getattr(self, key))
        else:
            header = ['file', 'plane', 'QP calc']
            row = [fname, self.plane, self.qp]
            
            for key in kwargs.keys():
                header.append(key)
                row.append(getattr(self, key))

        # write csv
        if not os.path.exists(fname):
            with open(fname, 'w') as f:
                writer = csv.writer(f, delimiter=' ')
                writer.writerow(header)

        with open(fname, 'a') as f:  
            writer=csv.writer(f, delimiter=' ')
            writer.writerow(row)
