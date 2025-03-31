import sys, glob

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


fig, ax = plt.subplots()
#------------------------------
# measurement data

df = pd.read_csv('../../2022/gr-20220823.csv', sep=' ').sort_values(by='QPV')
qpvs=np.unique(df['QPV'])
for qpv in qpvs:
    mean = df.loc[df['QPV'] == qpv, 'slope'].mean()*1e3
    error = (df.loc[df['QPV'] == qpv, 'slope']*1e3).std(ddof=0)
    markers, caps, bars = ax.errorbar(qpv-0.12, mean, yerr=error, c='r', mfc='w', fmt = 'o', ecolor='k', elinewidth = 1, capsize=5)
markers, caps, bars = ax.errorbar(qpv-0.12, mean, yerr=error, c='r', mfc='w', fmt = 'o', ecolor='k', elinewidth = 1, capsize=5, label='measurements')

#------------------------------
# simulations data
#pyheadtail data: slope error QPH QPV 

# resonator data
colors = ['#0C5DA5', '#00B945', '#FF9500', '#FF2C00', '#845B97', '#474747', '#9e9e9e']

dfHT = pd.read_csv('qpvscan.csv', sep=' ').sort_values(by='QPV') 
qpvvals = dfHT['QPV']
grvals = dfHT['slope']*1e3
errors = dfHT['error']*1e3

ax.plot(qpvvals, grvals, 'o-', color=colors[3], mec=colors[3], label='pyHT q26 + resonator')
ax.fill_between(qpvvals, grvals-errors, grvals+errors, color=colors[3], alpha=0.3)
ax.set_xlabel('Chromaticity $Q\'_v$ value')
ax.set_ylabel('Growthrate [1e3/nturns]')
ax.invert_xaxis()

axx = ax.twiny()
axx.plot(-qpvvals*1.6, grvals, color=colors[3], ls='None')
axx.set_xlabel(r'Chromaticity frequency $f_{\xi}$ [GHz]', fontweight='bold')
axx.set_xlim(0, 2*1.6)

#add nominal q26 non linear
dfHT = pd.read_csv('../qpvscan-nonlinearq26/qpvscan.csv', sep=' ').sort_values(by='QPV') 
qpvvals = dfHT['QPV']
grvals = dfHT['slope']*1e3
errors = dfHT['error']*1e3

ax.plot(qpvvals, grvals, 'o-k', mec='k', label='pyHT q26 wake')
ax.fill_between(qpvvals, grvals-errors, grvals+errors, color='k', alpha=0.3)

ax.set_xlabel('Chromaticity $Q\'_v$ value', fontweight='bold')
ax.set_ylabel('Growthrate [1e3/nturns]', fontweight='bold')
ax.set_title('Growthrate vs chromaticity Q26 \n Added resonator Rs=2e7 Ohm/m, Q=100, f=2.6GHz')
ax.set_xlim(0,-2)
ax.set_ylim((0,10))
ax.legend()

fig.savefig('qpvscan-meas.png')
fig.savefig('qpvscan-meas.svg')
fig.tight_layout()
plt.show()