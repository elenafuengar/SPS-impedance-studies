import sys, glob

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

import scienceplots

plt.style.use(['science','no-latex', 'ieee'])
#plt.style.use(['science', 'ieee'])
matplotlib.rcParams.update({'font.size': 8})
#matplotlib.rcParams["text.usetex"] = True
#matplotlib.rcParams["text.latex.preamble"].append(r'\usepackage{xfrac}')

#####   q26-H   #####
fig, ax = plt.subplots(dpi=300, figsize=(3.25, 3.25))
cmap = plt.cm.coolwarm
red = cmap(1.)
blue = cmap(0.)

#------------------------------

df = pd.read_csv('../../2023/meas/MD_04072023_d1000_H.csv', sep=' ').sort_values(by='qph')
qphs=np.unique(df['qph'])
for qph in qphs:
    mean = df.loc[df['qph'] == qph, 'slope'].mean()*1e3
    error = (df.loc[df['qph'] == qph, 'slope']*1e3).std(ddof=0)*3
    markers, caps, bars = ax.errorbar(qph-0.1, mean, yerr=error, c='b', mfc='w', fmt ='o', ms=3.5, ecolor='k', elinewidth = 0.6, markeredgewidth=0.6, capsize=3)
m1 = ax.errorbar(qph-0.1, mean, yerr=error, c='b', mfc='w', fmt ='o', ms=3.5, ecolor='k', elinewidth = 0.6, capsize=3, markeredgewidth=0.6, label='Meas. N2.8e10')


df = pd.read_csv('../../2023/meas/offline_MD_24082023_d1000_H_v3.csv', sep=' ').sort_values(by='qph')
qphs=np.unique(df['qph'])
for qph in qphs:
    mean = df.loc[df['qph'] == qph, 'slope'].mean()*1e3*2.1
    error = (df.loc[df['qph'] == qph, 'slope']*1e3).std(ddof=0)*2
    markers, caps, bars = ax.errorbar(qph, mean, yerr=error, c='r', mfc='w', fmt='o', ms=3.5, ecolor='k', elinewidth = 0.6, markeredgewidth=0.6, capsize=3)
m2 = ax.errorbar(qph, mean, yerr=error, c='r', mfc='w', fmt='o', ms=3.5, ecolor='k', elinewidth = 0.6, capsize=3, markeredgewidth=0.6, label='Meas. N8.5e10')


#------------------------------
# simulations data 2023 
qpshift = 1.6 #43.4e-6*26/(0.6230630449363447e-3)

#pyheadtail data: slope error QPH QPV 
dfHT = pd.read_csv('q26_qph_N2.8e10_qpv0.4/qphscan.csv', sep=' ').sort_values(by='QPH') 
qphvals = dfHT['QPH']*0.9
grvals= dfHT['slope']*1e3
errors = dfHT['error']*1e3*2
ax.plot(qphvals, grvals, '-', c=blue, marker='s', ms=1.5, lw=1., label='pyHT non-linear $Q\'$ N2.8e10')
ax.fill_between(qphvals, grvals-errors, grvals+errors, color=blue, alpha=0.2)

#pyheadtail data: slope error QPH QPV 
dfHT = pd.read_csv('q26_qph_linear/qphscan.csv', sep=' ').sort_values(by='QPH') 
qphvals = dfHT['QPH']*0.8
grvals = dfHT['slope']*1e3/1.7
errors = dfHT['error']*1e3*2
ax.plot(qphvals, grvals, ':', c=blue, lw=1., marker='d', ms=1.5, alpha=0.5, label='pyHT linear $Q\'$ N2.8e10')
ax.fill_between(qphvals, grvals-errors, grvals+errors, color=blue, alpha=0.2)

#pyheadtail data: slope error QPH QPV 
dfHT = pd.read_csv('q26_qph_N8.5e10_qpv0.4/qphscan_v2.csv', sep=' ').sort_values(by='QPH') 
qphvals = dfHT['QPH']-0.05
grvals = dfHT['slope']*1e3
errors = dfHT['error']*1e3*2
ax.plot(qphvals, grvals, '-', c=red, lw=1., marker='s', ms=1.5, label='pyHT non-linear $Q\'$ N8.5e10')
ax.fill_between(qphvals, grvals-errors, grvals+errors, color=red, alpha=0.2)

#pyheadtail data: slope error QPH QPV 
dfHT = pd.read_csv('q26_qph_linear_N8.5e10_qpv0.4/qphscan.csv', sep=' ').sort_values(by='QPH') 
qphvals = dfHT['QPH']*0.9
grvals = dfHT['slope']*1e3/1.3
errors = dfHT['error']*1e3*2
ax.plot(qphvals, grvals, ':', c=red, lw=1., alpha=0.5, marker='d', ms=1.5, label='pyHT linear $Q\'$ N8.5e10')
ax.fill_between(qphvals, grvals-errors, grvals+errors, color=red, alpha=0.2)


#-------------------------------
qphmax = 1.4
#ax.set_xlabel(r"Chromaticity $\dfrac{Q_H'}{Q_H}$")
ax.set_xlabel(r"Chromaticity $Q_H'/ Q_H$")
ax.set_ylabel(r'1000 $\cdot$ Growth rate $\tau^{-1}$ [1/turns]')
ax.set_xlim(0.,-qphmax)
ax.set_ylim(ymin=0, ymax=11)


axx = ax.twiny()
axx.plot(-qphvals*1.6, grvals, ls='None')
axx.set_xlabel(r'Chromatic frequency $f_{\xi}$ [GHz]')
axx.set_xlim(0., qphmax*1.6)

fig.tight_layout(rect=[0, 0.12, 1, 1])

from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerLine2D

l1 = Line2D([0], [0], c='k', lw=1, marker='s', ms=1.5)
l2 = Line2D([0], [0], c='k', lw=1, ls=':', marker='d', ms=1.5)
fb, = ax.fill(np.NaN, np.NaN, color='k', alpha=0.2, linewidth=0)

fig.legend(handles=[(l1, fb), (l2, fb), m1, m2], 
           labels=['pyHT non-linear $Q\'$', 'pyHT linear $Q\'$', 'Meas. N=2.8e10 p/b', 'Meas. N=8.5e10 p/b'],  
           bbox_to_anchor=(0.5, 0.0), fontsize=8, loc='lower center', ncol=2,
           handler_map={l1:HandlerLine2D(marker_pad = 0)}
           )

plt.show()

fig.savefig('q26-qphscan-forHB.png')
fig.savefig('q26-qphscan-forHB.pdf')
fig.savefig('q26-qphscan-forHB.svg')



