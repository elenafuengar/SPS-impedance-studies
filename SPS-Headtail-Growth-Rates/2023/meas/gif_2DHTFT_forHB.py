'''BQHT headtail monitor file analysis

Jul 23
'''

import glob
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.constants import c

#path = '/nfs/cs-ccr-bqhtnfs/sps_data/SPS.BQHT/2023_07_04/'
#files = sorted(glob.glob(path+'*.h5'))

### Plot 2D FFT
case = 'QPH'
stored_turns = 300

# QPH files
if case == 'QPH':
    files = {
    '-1.4' : 'SPS.BQHT_MD2_20230824_115423.h5', #QPH=-1.4
    '-1.3' : 'SPS.BQHT_MD2_20230824_115135.h5', #QPH=-1.3
    '-0.9' : 'SPS.BQHT_MD2_20230824_112211.h5', #QPH=-1.0
    '-0.8' : 'SPS.BQHT_MD2_20230824_111017.h5', #QPH=-0.8
    '-0.7' : 'SPS.BQHT_MD2_20230824_110811.h5', #QPH=-0.7
    '-0.6' : 'SPS.BQHT_MD2_20230824_110359.h5', #QPH=-0.6
    '-0.55' : 'SPS.BQHT_MD2_20230824_105823.h5', #QPH=-0.55
    '-0.4' : 'SPS.BQHT_MD2_20230824_105247.h5', #QPH=-0.4
    '-0.35' : 'SPS.BQHT_MD2_20230824_104917.h5', #QPH=-0.35
    '-0.2' : 'SPS.BQHT_MD2_20230824_103229.h5', #QPH=-0.2 
    '-0.05' : 'SPS.BQHT_MD2_20230704_113250.h5', #QPH=0.05
    }

#QPV files
if case == 'QPV':
    files =  {
        '-1.32' : 'SPS.BQHT_MD2_20230703_182829.h5',
        '-1.12' : 'SPS.BQHT_MD2_20230703_182452.h5',
        '-1.02' : 'SPS.BQHT_MD2_20230703_170914.h5',
        '-0.92' : 'SPS.BQHT_MD2_20230703_181726.h5',
        '-0.82' : 'SPS.BQHT_MD2_20230703_170608.h5',
        '-0.72' : 'SPS.BQHT_MD2_20230703_181342.h5',
        '-0.52' : 'SPS.BQHT_MD2_20230703_181151.h5',
        '-0.0' : 'SPS.BQHT_MD2_20230703_174855.h5',
    }

for k, key in enumerate(sorted(files.keys())):

    file = files[key]

    with h5.File(file, 'r') as hf:

        delta_H = np.array(hf.get('horizontal/delta/'))
        delta_V = np.array(hf.get('vertical/delta/'))
        len_z_H = int(len(delta_H)/stored_turns)
        len_z_V = int(len(delta_V)/stored_turns)

        #Reshape
        delta_H = np.reshape(delta_H, [stored_turns, len_z_H])
        delta_V = np.reshape(delta_V, [stored_turns, len_z_V])

    for i in range(stored_turns):
        # Remove baseline
        delta_H[i, :] = delta_H[i, :] - np.mean(delta_H, axis=0) 
        delta_V[i, :]  = delta_V[i, :] - np.mean(delta_V, axis=0) 

    bl = 3 #ns

    # 2D fourier analysis
    if case == 'QPH':
        delta = delta_H
    elif case == 'QPV':
        delta = delta_V

    import scienceplots
    plt.style.use(['science','no-latex', 'ieee'])
    matplotlib.rcParams.update({'font.size': 8})
    fig, axs = plt.subplots(1,2, dpi=300, figsize=(3.25, 1.65), width_ratios=[2, 1.2])

    ft = np.fft.ifftshift(delta)
    ft = np.fft.fft2(delta)
    ft = np.fft.fftshift(ft)
    magnitude2d = np.abs(ft)

    # Compute phase
    ftleft = np.abs(ft[:(delta.shape[0]//2-5), :])
    xphase, yphase = np.unravel_index(np.argmax(ftleft), ftleft.shape)

    slice_i = 0
    slice_f = 500
    t = np.linspace(-bl*14, bl*10, 500)*1e-9 #s

    if key == '-1.3' or key == '-1.4' or key == '-1.32' or key == '-1.12':
        t = np.linspace(-bl*15, bl*10, 500)*1e-9 #s

    n = np.arange(stored_turns) #s

    fx = np.fft.fftfreq(len(n), d=1/43e3) #sps rev freq
    fy = np.fft.fftfreq(len(t), d=(t[2]-t[1]))

    extent_t = [-200, 100, t[0]*1e9, t[-1]*1e9]
    if key =='-0.72':
        extent_t = [-230, 70, t[0]*1e9, t[-1]*1e9]
    axs[0].imshow(delta[:, slice_i:slice_f].transpose(), cmap='rainbow', origin='lower', interpolation='spline36', extent=extent_t, aspect='auto')
    extent_f = [-fx.max()/1e3/5.4, fx.max()/1e3/5.4, -1.2, 1.2]
    if case == 'QPV':
        extent_f = [-fx.max()/1e3/7.4, fx.max()/1e3/7.4, -1.2, 1.2]
    axs[1].imshow(magnitude2d[:, slice_i:slice_f].transpose(), vmax=np.max(magnitude2d)/2, cmap='rainbow', origin='lower', interpolation='spline36', extent=extent_f, aspect='auto')
    if key == '-0.82' or key == '-1.02':
        axs[1].imshow(magnitude2d[:, slice_i:slice_f].transpose(), vmax=np.max(magnitude2d)/10, cmap='rainbow', origin='lower', interpolation='spline36', extent=extent_f, aspect='auto')
    ft = np.fft.ifftshift(delta)
    ft = np.fft.fft2(delta)
    ft = np.fft.fftshift(ft)
    magnitude2d = np.abs(ft)

    if case == 'QPH':
        axs[0].set_title(f'HT monitor | $Q\'_H={key}$')
    elif case == 'QPV':
        axs[0].set_title(f'HT monitor | $Q\'_V={key}$')
        
    axs[0].set_ylabel('z position [ns]', labelpad=0.05)
    axs[0].set_xlabel('Turns')
    axs[0].set_ylim(-bl, bl)
    axs[0].set_xlim(0,70)
    #axs[0].text(0.65, 0.90, f'$Q\'_H={key}$', ha='center',va='center', transform=axs[0].transAxes)

    axs[1].set_title('2D FFT Magnitude')
    axs[1].set_ylabel('$\Delta\psi_{HT}$ [a.u.]', labelpad=0.05 )
    axs[1].set_xlabel(r'Betatron phase')
    axs[1].set_ylim(-1,1)
    axs[1].set_xlim(-2,2)

    fig.tight_layout(pad=0.3)
    fig.savefig(f'gif2DHTFT/{str(k).zfill(2)}_2DHTFT_{case}{key}.png')
    #plt.show()
