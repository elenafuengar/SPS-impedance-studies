import sys, glob
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
from growthrate import GrowthRate
import h5py as h5 

End_turns = 200
Start_turns = 100
Freq_turns = 1
stored_turns = 300

parquet = '2023.07.03.16.56.17.560516.parquet' #'2023.07.04.11.34.41.701026.parquet'
htfile = 'SPS.BQHT_MD2_20230704_113250.h5'
sim = '../pyheadtail/q26_qpv_N2.8e10_qph0.2/bunchmonitor_qpv-0.1_q26.h5'

#grsim = GrowthRate(file=sim, nperseg=65, noverlap=50, plane='V')
#grsim.MWFFT_analysis()

gr = GrowthRate(file=parquet, nperseg=65, noverlap=50, plane='V', grstart=45000, turn_min=0, turn_max=88000)
gr.MWFFT_analysis()

import scienceplots
plt.style.use(['science','no-latex', 'ieee'])
matplotlib.rcParams.update({'font.size': 8})
plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.rainbow(np.linspace(0,1, int( (End_turns - Start_turns)/ Freq_turns) )))

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, dpi=300, figsize=(3.25, 3.25))

ax1.plot(gr.intms, gr.inty, color='darkorange')
ax1.set_xlabel('Time [ms]')
ax1.set_ylabel('Intensity [p/b]')
ax1.set_xlim(0,5000)
ax1.text(0.65, 0.90, f'$Q\'_V=-0.2$', ha='center',va='center', transform=ax1.transAxes)

ax2.plot(gr.ms*1e3, gr.pos, color='b')
ax2.set_xlabel('Time [ms]')
ax2.set_ylabel('Bunch $\hat{y}$')
ax2.set_xlim(1015,1135)

with h5.File(htfile, 'r') as hf:

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

z = np.linspace(-0.2*4*6, 0.2*4*5, len(delta_H[i,:]))
for i in range(Start_turns,End_turns,Freq_turns):
    ax3.plot(z, delta_H[i,:]/10000)

ax3.set_xlabel('z [m]')
ax3.set_ylabel('Bunch $\Delta_y$')
ax3.set_xlim(-0.2*4, 0.2*4)
ax3.text(0.7, 0.92, '100 turns', ha='center',va='center', transform=ax3.transAxes)

ax4.plot(gr.points, gr.log_peak_amplitude_all,'.-b', label='MWFFT')
ax4.plot(gr.grx, gr.gry, '.-r', label='points for fit')
ax4.plot(gr.grx, gr.poly1d_fn(gr.grx), '--k', lw=2.0, label='fit masked data')
ax4.text(0.5, 0.70, r'$\tau^{-1}$', ha='center',va='center', transform=ax4.transAxes, fontsize=12)
ax4.set_xlabel('Turns')
ax4.set_ylabel('$\log_{10}$ of MWFT Amplitude')
ax4.set_xlim(1015*43.35,1135*43.35)
ax4.set_ylim(8,22)

fig.tight_layout()
fig.savefig('GR_forHB.pdf')
fig.savefig('GR_forHB.png')
fig.savefig('GR_forHB.svg')
plt.show()