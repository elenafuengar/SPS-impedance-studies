import sys, glob

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

import scienceplots
plt.style.use(['science','no-latex', 'ieee'])
matplotlib.rcParams.update({'font.size': 8})

#####   q26-H   #####
fig, ax = plt.subplots(dpi=300, figsize=(3.25, 3.25))
cmap = plt.cm.coolwarm
red = cmap(1.)
blue = cmap(0.)

#------------------------------

# measurement data

df = pd.read_csv('../../2022/gr-20220823.csv', sep=' ').sort_values(by='QPV')
qpvs=np.unique(df['QPV'])
for qpv in qpvs:
    mean = df.loc[df['QPV'] == qpv, 'slope'].mean()*1e3
    error = (df.loc[df['QPV'] == qpv, 'slope']*1e3).std(ddof=0)
    markers, caps, bars = ax.errorbar(qpv-0.12, mean, yerr=error, c='grey', mfc='w', ms=3., fmt = 'o', ecolor='k', elinewidth=0.6, capsize=3, markeredgewidth=0.6)
m1 = ax.errorbar(qpv-0.12, mean, yerr=error, c='grey', mfc='w', fmt = 'o', ecolor='k',  ms=3., elinewidth=0.6, capsize=3, markeredgewidth=0.6, label='meas. 30/08/22')

df = pd.read_csv('../../2023/meas/offline_MD_04072023_d1000_V.csv', sep=' ').sort_values(by='qpv')
qpvs=np.unique(df['qpv'])
for qpv in qpvs:
    mean = df.loc[df['qpv'] == qpv, 'slope'].mean()*1e3*3/2.5
    error = (df.loc[df['qpv'] == qpv, 'slope']*1e3).std(ddof=0)
    markers, caps, bars = ax.errorbar(qpv, mean, yerr=error, c='r', mfc='w',  ms=3.5, fmt = 'o', ecolor='k', elinewidth=0.6, markeredgewidth=0.6, capsize=3)
m2 = ax.errorbar(qpv, mean, yerr=error, c='r', mfc='w', fmt = 'o', ecolor='k',  ms=3.5, elinewidth=0.6, markeredgewidth=0.6, capsize=3, label='meas. 04/07/23')

df = pd.read_csv('../../2023/meas/offline_MD_03072023_d1000_V.csv', sep=' ').sort_values(by='qpv')
qpvs=np.unique(df['qpv'])
for qpv in qpvs:
    mean = df.loc[df['qpv'] == qpv, 'slope'].mean()*1e3*3/2.5
    error = (df.loc[df['qpv'] == qpv, 'slope']*1e3).std(ddof=0)
    markers, caps, bars = ax.errorbar(qpv, mean, yerr=error, c='b', mfc='w', fmt = 'o',  ms=3.5, ecolor='k', elinewidth=0.6, markeredgewidth=0.6, capsize=3)
m3 = ax.errorbar(qpv, mean, yerr=error, c='b', mfc='w', fmt = 'o', ecolor='k',  ms=3.5, elinewidth=0.6, markeredgewidth=0.6, capsize=3, label='meas. 03/07/23')

#------------------------------
# simulations data 2023 
qpshift = 1.6 #43.4e-6*26/(0.6230630449363447e-3)

#pyheadtail data: slope error QPH QPV 
dfHT = pd.read_csv('q26_qpv_linear/qpvscan.csv', sep=' ').sort_values(by='QPV') 
qpvvals = dfHT['QPV']
grvals = dfHT['slope']*1e3
errors = dfHT['error']*1e3*10
l2, = ax.plot(qpvvals, grvals, ':', color='grey', marker='s', ms=1.5, lw=0.8, label='pyHT q26 linear')
ax.fill_between(qpvvals, grvals-errors, grvals+errors, color='k', alpha=0.15)

#pyheadtail data: slope error QPH QPV 
dfHT = pd.read_csv('q26_qpv_N2.8e10_qph0.2/qpvscan.csv', sep=' ').sort_values(by='QPV') 
qpvvals = dfHT['QPV']
grvals = dfHT['slope']*1e3
errors = dfHT['error']*1e3*10
l1, = ax.plot(qpvvals, grvals, '-k', marker='d', ms=1.5, lw=0.8, label='pyHT q26 non-linear')
fb = ax.fill_between(qpvvals, grvals-errors, grvals+errors, color='k', alpha=0.15)



#-------------------------------
qpvmax = 2.0
#ax.set_xlabel('Chromaticity $Q\'_V$ [knob]')
ax.set_xlabel(r"Chromaticity $Q_V'/ Q_V$")
ax.set_ylabel(r'1000 $\cdot$ Growth rate $\tau^{-1}$ [1/turns]')
ax.set_xlim(0.,-qpvmax)
ax.set_ylim(ymin=0, ymax=10)


axx = ax.twiny()
axx.plot(-qpvvals*1.6, grvals, ls='None')
axx.set_xlabel(r'Chromatic frequency $f_{\xi}$ [GHz]')
axx.set_xlim(0., qpvmax*1.6)

fig.tight_layout(rect=[0, 0.12, 1, 1])

from matplotlib.legend_handler import HandlerLine2D
from matplotlib.legend_handler import HandlerTuple

fig.legend(handles=[(fb, l1), (fb, l2), m1, (m2,m3)], 
           labels=['pyHT non-linear $Q\'$', 'pyHT linear $Q\'$', 'Meas. 08/2022', 'Meas. 07/2023'],  
           bbox_to_anchor=(0.5, 0.0), fontsize=8, loc='lower center', ncol=2,
           handler_map={(m2, m3): HandlerTuple(ndivide=None), l1:HandlerLine2D(marker_pad = 0)} 
            )

plt.show()

fig.savefig('q26-qpvscan-forHB.png')
fig.savefig('q26-qpvscan-forHB.svg')
fig.savefig('q26-qpvscan-forHB.pdf')



