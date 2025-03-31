import glob
import time
import pandas as pd
from HTchromaticity import HTchromaticity
from datetime import date

today = '2024_03_12'  #date.today().isoformat()
datapath = f'/nfs/cs-ccr-bqhtnfs/sps_data/SPS.BQHT/{today}/'

#monitor
while True:

    flist = sorted(glob.glob(datapath+'*.h5'))
    file = flist[-1]
    print(f'Analyzing file {file.split("/")[-1]}')

    #analysis
    ht = HTchromaticity(file=file, nturns=500, plane='H', knob=None, monitor=False, bl=3, dpi=150)
    ht.calc_chromaticity()
    ht.plot()

    inp = input('Save? y/n:')
    if inp == 'y': 
        ht.save(fname=f'HT_{today}.csv')

    #time.sleep(60) #supercycle time


