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
    gr = GrowthRate(file=file, nperseg=65, noverlap=40, grstart=2000)
    gr.MWFFT_analysis()
    #gr.plot()
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
ax.set_title('Growthrate vs chromaticity simulations')
ax.invert_xaxis()

axx = ax.twiny()
axx.plot(-qpvvals*1.6, grvals, '--k')
axx.set_xlabel(r'Chromaticity frequency $f_{\xi}$ [GHz]')

fig.tight_layout()
plt.show()

#read measurement results .csv
#0,time,QPH,QPV,grstart,grend,slope,intensity,error,nturns,acqPeriod,current_start,current_end,current_mapping,QPV_error

#df1 = pd.read_csv('../../2022/20220823_with_errors_new.csv').sort_values(by='QPV') 
#df2 = pd.read_csv('../../2022/20220830_with_errors_new.csv').sort_values(by='QPV') 

df1 = pd.read_csv('../../2022/growthrate_220823.csv').sort_values(by='QPV') 
df2 = pd.read_csv('../../2022/growthrate_220830.csv').sort_values(by='QPV') 
df3 = pd.read_csv('../../2022/growthrate_220906.csv').sort_values(by='QPV') 

#df1 = pd.read_csv('../../2022/20220823xls.csv').sort_values(by='QPV') 
#df2 = pd.read_csv('../../2022/20220830xls.csv').sort_values(by='QPV') 

fig, ax = plt.subplots()

ax.plot(qpvvals, grvals, '+-k', label='pyHEADTAIL')
ax.plot(df1['QPV']-0.12, df1['slope']*1e3, 'or', label='meas. 23/08/22', alpha=0.6)
ax.plot(df2['QPV']-0.12, df2['slope']*1e3, 'ob', label='meas. 30/08/22', alpha=0.6)
ax.plot(df3['QPV']-0.12, df3['slope']*1e3*10, 'og', label='meas. 06/09/22', alpha=0.6)

#ax.plot(df1['QPV']-0.12, df1['GR']*1e3, 'o--r', mec='k', ls=1.5, label='meas. 23/08/22', alpha=0.6)
#ax.plot(df2['QPV']-0.12, df2['GR']*1e3, 'o--b', mec='k', ls=1.5, label='meas. 30/08/22', alpha=0.6)


ax.set_xlabel('Chromaticity $Q\'_v$ value')
ax.set_ylabel('Growthrate [1e3/nturns]')
ax.set_title('Growthrate vs chromaticity | non-linear Q\', Q26', fontweight='bold')
ax.set_ylim((0,10))
ax.invert_xaxis()
ax.legend()

fig.tight_layout()
plt.show()
