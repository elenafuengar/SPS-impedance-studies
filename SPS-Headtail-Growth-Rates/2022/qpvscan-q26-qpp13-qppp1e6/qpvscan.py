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
    gr = GrowthRate(file=file, nperseg=65, noverlap=40, grstart=3500)
    gr.MWFFT_analysis()
    gr.plot()
    gr.save(fname='qpvscan.csv', QPH=qph, QPV=qpv)
    grvals = np.append(grvals, gr.slope*1e3)

#sort
ind = np.argsort(qpvvals)
qpvvals = qpvvals[ind]
grvals = grvals[ind]

#plot
fig, ax = plt.subplots()

ax.plot(qpvvals, grvals, 'ok', mec='k', label='pyHEADTAIL')
ax.set_xlabel('Chromaticity $Q\'_v$ value')
ax.set_ylabel('Growthrate [1e3/nturns]')
ax.set_title('Growthrate vs chromaticity Q26 \n $Q\'\'_v = 13$, $Q\'\'\'_v = {10}^{6}$')
ax.invert_xaxis()

axx = ax.twiny()
axx.plot(-qpvvals*1.6, grvals, '--k')
axx.set_xlabel(r'Chromaticity frequency $f_{\xi}$ [GHz]')

fig.tight_layout()
plt.show()

