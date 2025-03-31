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
    gr = GrowthRate(file=file, nperseg=65, noverlap=40, grstart=3000)
    gr.MWFFT_analysis()
    #gr.plot()
    gr.save(fname='qpvscan.csv', QPH=qph, QPV=qpv)
    grvals = np.append(grvals, gr.slope*1e3)
