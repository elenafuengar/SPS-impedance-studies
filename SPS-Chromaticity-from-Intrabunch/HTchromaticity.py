'''
Chromaticity estimation from
SPS BQHT headtail monitor data

`.h5` files are stored in:

/nfs/cs-ccr-bqhtnfs/sps_data/SPS.BQHT/yyyy_mm_dd

Mar 24
'''

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import csv
import os

class HTchromaticity():

    def __init__(self, file=None, nturns=300, plane='V', knob=None,
                 bl=3, monitor=False, save_png=False, dpi=100, **kwargs):

        self.file = file
        self.nturns = nturns
        self.plane = plane
        self.bl = bl #TODO use measurement --> 25 ns discussed with Tom
        self.knob = knob
        self.monitor = monitor
        self.dpi = dpi
        self.save_png = save_png

        # plotting
        try:
            import scienceplots
            plt.style.use(['science','no-latex', 'ieee'])
        except:
            pass
        matplotlib.rcParams.update({'font.size': 8})

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

        self.get_data()

    def get_data(self):

        f = h5py.File(self.file)

        # Read and reshape
        if self.plane == 'H':
            self.delta = f['horizontal']['delta']
        elif self.plane == 'V':
            self.delta = f['vertical']['delta']
        else:
            raise Exception('Specify a valid plane: "H" or "V"')

        self.len = int(len(self.delta)/self.nturns)
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
        self.xp, self.yp = np.unravel_index(np.argmax(ftleft), ftleft.shape)
        phase = -(self.yp - self.ft.shape[1]//2) / (self.ft.shape[1]//2)
        self.qp = phase*3  #TODO normalize to chromaticity

    def plot(self):

        if self.monitor:
            for ax in self.axs:
                ax.cla()
        else:
            self.fig, self.axs = plt.subplots(1,2, dpi=self.dpi, figsize=(6.25, 3.25), num=1)

        t = np.linspace(-self.bl*10, self.bl*10, self.nturns)*1e-9 #s
        n = np.arange(self.nturns) #s
        fx = np.fft.fftfreq(len(n), d=1/43e3) #sps rev freq
        fy = np.fft.fftfreq(len(t), d=(t[2]-t[1]))

        # z-t HT signal
        extent_t = [0, self.nturns, t[0]*1e9, t[-1]*1e9]
        ht = self.axs[0].imshow(self.delta.transpose(), cmap='rainbow', origin='lower', interpolation='spline36', extent=extent_t, aspect='auto')
        
        # 2D FFT signal
        extent_f = [-fx.max()/1e3, fx.max()/1e3, -fy.max()/1e9, fy.max()/1e9]
        ft = self.axs[1].imshow(self.ft.transpose(), cmap='rainbow', origin='lower', interpolation='spline36', extent=extent_f, aspect='auto')

        self.axs[0].set_title(f'HT monitor - Bunch $\Delta_{self.plane}$')
        self.axs[0].set_ylabel('Time [ns]')
        self.axs[0].set_xlabel('Turns')
        self.axs[0].set_ylim(-3*self.bl,3*self.bl)
        #self.axs[0].set_xlim(100,200)
        self.axs[0].text(0.65, 0.95, f'$Q\'_{self.plane}={round(self.qp,2)}$', ha='center',va='center', transform=self.axs[0].transAxes)

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