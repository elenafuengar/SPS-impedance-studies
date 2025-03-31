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
dfHT = pd.read_csv('qpvscan.csv', sep=' ').sort_values(by='QPV') 
qpvvals = dfHT['QPV']
grvals = dfHT['slope']*1e3
errors = dfHT['error']*1e3

ax.plot(qpvvals, grvals, 'o-g', mec='g', label='pyHT $Q\'\'_v=13$, $Q\'\'\'_v=1e6$')
ax.fill_between(qpvvals, grvals-errors, grvals+errors, color='g', alpha=0.3)
ax.set_xlabel('Chromaticity $Q\'_v$ value')
ax.set_ylabel('Growthrate [1e3/nturns]')
ax.set_title('Growthrate vs chromaticity simulations')
ax.invert_xaxis()

axx = ax.twiny()
axx.plot(-qpvvals*1.6, grvals, '-k')
axx.set_xlabel(r'Chromaticity frequency $f_{\xi}$ [GHz]')
axx.set_xlim(0, 2*1.6)

#add nominal q26 non linear
dfHT = pd.read_csv('../qpvscan-nonlinearq26/qpvscan.csv', sep=' ').sort_values(by='QPV') 
qpvvals = dfHT['QPV']
grvals = dfHT['slope']*1e3
errors = dfHT['error']*1e3

ax.plot(qpvvals, grvals, 'o-k', mec='k', label='pyHT $Q\'\'_v=130$, $Q\'\'\'_v=2e5$')
ax.fill_between(qpvvals, grvals-errors, grvals+errors, color='k', alpha=0.3)

ax.set_xlabel('Chromaticity $Q\'_v$ value')
ax.set_ylabel('Growthrate [1e3/nturns]')
ax.set_title('Growthrate vs chromaticity | non-linear Q\', Q26', fontweight='bold')
ax.set_xlim(0,-2)
ax.set_ylim((0,10))
ax.legend()

fig.savefig('qpvscan-meas.png')
fig.savefig('qpvscan-meas.svg')
fig.tight_layout()
plt.show()
