import sys, glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

sys.path.append('../')
from growthrate import GrowthRate

h5list = glob.glob('*.h5')
qpv = 0.4
qphvals, grvals = np.array([]), np.array([])

#analysis
for file in h5list:

    print(f'Analyzing file {file.split("/")[-1]}')
    #get params
    qph = float(file.split('qph')[1].split('_')[0])
    qphvals = np.append(qphvals, qph)

    #analysis
    gr = GrowthRate(file=file, nperseg=65, plane='H', noverlap=40, grstart=4000)
    gr.MWFFT_analysis()
    #gr.plot()
    gr.save(fname='qphscan.csv', QPH=qph, QPV=qpv)
    grvals = np.append(grvals, gr.slope*1e3)

#sort
ind = np.argsort(qphvals)
qphvals = qphvals[ind]
grvals = grvals[ind]

# plot simulation
# ---------------

fig, ax = plt.subplots()

ax.plot(qphvals, grvals, 'ok', mec='k', label='pyHEADTAIL')
ax.set_xlabel('Chromaticity $Q\'_v$ value')
ax.set_ylabel('Growthrate [1e3/nturns]')
ax.set_title('Growthrate vs chromaticity simulations \n Added resonator Q=100, f=2.2GHz, R=2e6 Ohm/m ')
ax.invert_xaxis()

axx = ax.twiny()
axx.plot(-qphvals*1.6, grvals, '--k')
axx.set_xlabel(r'Chromaticity frequency $f_{\xi}$ [GHz]')

fig.tight_layout()
plt.show()

# plot with measurements
# -----------------------

fig, ax = plt.subplots()

# measurement data
df = pd.read_csv('../../../2023/meas/MD_04072023_d1000_H.csv', sep=' ').sort_values(by='qph')
qphs=np.unique(df['qph'])
for qph in qphs:
    mean = df.loc[df['qph'] == qph, 'slope'].mean()*1e3*8.5/3
    error = (df.loc[df['qph'] == qph, 'slope']*1e3).std(ddof=0)*2
    markers, caps, bars = ax.errorbar(qph-0.12, mean, yerr=error, c='b', mfc='w', fmt = 'o', ecolor='k', elinewidth = 1, capsize=5)
markers, caps, bars = ax.errorbar(qph-0.12, mean, yerr=error, c='b', mfc='w', fmt = 'o', ecolor='k', elinewidth = 1, capsize=5, label='online N3e10')

df = pd.read_csv('../../../2023/meas/offline_MD_04072023_d1000_H.csv', sep=' ').sort_values(by='qph')
qphs=np.unique(df['qph'])
for qph in qphs:
    mean = df.loc[df['qph'] == qph, 'slope'].mean()*1e3*8.5/3
    error = (df.loc[df['qph'] == qph, 'slope']*1e3).std(ddof=0)*2
    markers, caps, bars = ax.errorbar(qph-0.12, mean, yerr=error, c='b', mfc='w', fmt = 'o', ecolor='k', elinewidth = 1, capsize=5, alpha=0.5)
markers, caps, bars = ax.errorbar(qph-0.12, mean, yerr=error, c='b', mfc='w', fmt = 'o', ecolor='k', elinewidth = 1, capsize=5, alpha=0.5, label='offline N3e10')

# measurement data
df = pd.read_csv('../../../2023/meas/q26_24082023_d1000_H.csv', sep=' ').sort_values(by='qph')
qphs=np.unique(df['qph'])
for qph in qphs:
    mean = df.loc[df['qph'] == qph, 'slope'].mean()*1e3*2
    error = (df.loc[df['qph'] == qph, 'slope']*1e3).std(ddof=0)*2
    markers, caps, bars = ax.errorbar(qph, mean, yerr=error, c='r', mfc='w', fmt = 'o', ecolor='k', elinewidth = 1, capsize=5)
markers, caps, bars = ax.errorbar(qph, mean, yerr=error, c='r', mfc='w', fmt = 'o', ecolor='k', elinewidth = 1, capsize=5, label='online N8e10')


#--------------

#pyheadtail data: slope error QPH QPH 
dfHT = pd.read_csv('qphscan.csv', sep=' ').sort_values(by='QPH') 
qphvals = dfHT['QPH']
grvals = dfHT['slope']*1e3
errors = dfHT['error']*2*1e3

ax.plot(qphvals, grvals, 'o-k', mec='k', label='pyHT non-linear')
ax.fill_between(qphvals, grvals-errors, grvals+errors, color='k', alpha=0.2)
ax.set_xlabel('Chromaticity $Q\'_h$ value')
ax.set_ylabel('Growthrate [1e3/nturns]')
ax.set_title('Growthrate vs chromaticity simulations')
ax.invert_xaxis()

axx = ax.twiny()
axx.plot(-qphvals*1.6, grvals, '-k')
axx.set_xlabel(r'Chromaticity frequency $f_{\xi}$ [GHz]')
axx.set_xlim(0, 2*1.6)

ax.set_xlabel('Chromaticity $Q\'_v$ value')
ax.set_ylabel('Growthrate [1e3/nturns]')
#ax.set_title('Growthrate vs chromaticity | non-linear Q\', Q26 \n Added resonator Q=100, f=2.2GHz, R=2e6 Ohm/m ', fontweight='bold')
ax.set_xlim(0,-2)
#ax.set_ylim((0,3))
ax.legend()

fig.savefig('qphscan-meas.png')
fig.savefig('qphscan-meas.svg')
fig.tight_layout()
plt.show()
