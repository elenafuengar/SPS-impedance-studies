import sys, glob
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from growthrate import GrowthRate

user = 'SPS.USER.MD2'
date = '2023-08-24'
datapath = f'/user/spsscrub/2023/sps_beam_monitoring/sps_beam_monitoring/data/GR/{user}/{date}/'
datapath = f'/user/spsscrub/2024/sps_beam_monitoring/sps_beam_monitoring/general/GR/data/GR/{user}/{date}/'

#analysis
while True:

    parquetlist = sorted(glob.glob(datapath+'*.parquet'))
    file = parquetlist[-1]
    print(f'Analyzing file {file.split("/")[-1]}')

    #analysis
    gr = GrowthRate(file=file, nperseg=65, noverlap=50, plane='V', 
                    pad=3, turn_min=43000, turn_max=80000)
    print(f'QPH = {gr.qph}, QPV = {gr.qpv}')
    gr.MWFFT_analysis()
    gr.plot()
    inp = input('Save? y/n:')
    if inp == 'y': 
        gr.save(fname='q26_24082023_d1000_H_11e10.csv', header=['slope', 'error', 'qpv', 'qph'])

    print('waiting for cycle...')
    time.sleep(20) #supercycle time
