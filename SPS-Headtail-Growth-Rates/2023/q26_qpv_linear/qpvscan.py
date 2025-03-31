import sys, glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

sys.path.append('../')
from growthrate import GrowthRate

h5list = glob.glob('*.h5')
qph = 0.1
qpvvals, grvals = np.array([]), np.array([])

#analysis
for file in h5list:

    print(f'Analyzing file {file.split("/")[-1]}')
    #get params
    qpv = float(file.split('qpv')[1].split('_')[0])
    qpvvals = np.append(qpvvals, qpv)

    #analysis
    gr = GrowthRate(file=file, nperseg=65, noverlap=50, grstart=3000)
    gr.MWFFT_analysis()
    gr.plot()
    gr.save(fname='qpvscan.csv', QPH=qph, QPV=qpv)
    grvals = np.append(grvals, gr.slope*1e3)

#sort
ind = np.argsort(qpvvals)
qpvvals = qpvvals[ind]
grvals = grvals[ind]

#plot simulation
fig, ax = plt.subplots()

ax.plot(qpvvals, grvals, 'ok', mec='k', label='pyHEADTAIL')
ax.set_xlabel('Chromaticity $Q\'_v$ value')
ax.set_ylabel('Growthrate [1e3/nturns]')
ax.set_title('Growthrate vs chromaticity simulations \n Added resonator Q=100, f=2.2GHz, R=2e6 Ohm/m ')
ax.invert_xaxis()

axx = ax.twiny()
axx.plot(-qpvvals*1.6, grvals, '--k')
axx.set_xlabel(r'Chromaticity frequency $f_{\xi}$ [GHz]')

fig.tight_layout()
plt.show()


fig, ax = plt.subplots()

#read measurement results and plot .csv
df = pd.read_csv('../../../2023/meas/offline_MD_04072023_d1000_V.csv', sep=' ').sort_values(by='qpv')
qpvs=np.unique(df['qpv'])
for qpv in qpvs:
    mean = df.loc[df['qpv'] == qpv, 'slope'].mean()*1e3
    error = (df.loc[df['qpv'] == qpv, 'slope']*1e3).std(ddof=0)
    markers, caps, bars = ax.errorbar(qpv, mean, yerr=error, c='r', mfc='w', fmt = 'o', ecolor='k', elinewidth = 1, capsize=5)
markers, caps, bars = ax.errorbar(qpv, mean, yerr=error, c='r', mfc='w', fmt = 'o', ecolor='k', elinewidth = 1, capsize=5, label='measurements 04/07/23')

'''
#read measurement results and plot .csv
df = pd.read_csv('../../../2023/meas/MD_03072023_d1000.csv', sep=' ').sort_values(by='qpv')
qpvs=np.unique(df['qpv'])
for qpv in qpvs:
    mean = df.loc[df['qpv'] == qpv, 'slope'].mean()*1e3
    error = (df.loc[df['qpv'] == qpv, 'slope']*1e3).std(ddof=0)
    markers, caps, bars = ax.errorbar(qpv, mean, yerr=error, c='b', mfc='w', fmt = 'o', ecolor='k', elinewidth = 1, capsize=5)
markers, caps, bars = ax.errorbar(qpv, mean, yerr=error, c='b', mfc='w', fmt = 'o', ecolor='k', elinewidth = 1, capsize=5, label='measurements 03/07/23')
'''

#pyheadtail data: slope error QPH qpv 
dfHT = pd.read_csv('qpvscan.csv', sep=' ').sort_values(by='QPV') 
qpvvals = dfHT['QPV']
grvals = dfHT['slope']*1e3
errors = dfHT['error']

ax.plot(qpvvals, grvals, 'o-k', mec='k', label='pyHEADTAIL')
ax.fill_between(qpvvals, grvals-errors, grvals+errors, color='k', alpha=0.6)
ax.set_xlabel('Chromaticity $Q\'_v$ value')
ax.set_ylabel('Growthrate [1e3/nturns]')
ax.set_title('Growthrate vs chromaticity simulations')
ax.invert_xaxis()

axx = ax.twiny()
axx.plot(-qpvvals*1.6, grvals, '-k')
axx.set_xlabel(r'Chromaticity frequency $f_{\xi}$ [GHz]')
axx.set_xlim(0, 2*1.6)

ax.set_xlabel('Chromaticity $Q\'_v$ value')
ax.set_ylabel('Growthrate [1e3/nturns]')
ax.set_xlim(0,-2)
ax.set_ylim((0,8))
ax.legend()

fig.savefig('qpvscan-meas.png')
fig.savefig('qpvscan-meas.svg')
fig.tight_layout()
plt.show()
