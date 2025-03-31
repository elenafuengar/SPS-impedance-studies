import glob
import time
import pandas as pd
from HTchromaticity import HTchromaticity
from datetime import date

path = '/home/edelafue/cernbox/measurements/HTModeZeroGrow/2023/meas/'
#flist = glob.glob(datapath+'*.h5')

# QPH files
files = {
'-1.4' : 'SPS.BQHT_MD2_20230824_115423.h5', #QPH=-1.4
'-1.3' : 'SPS.BQHT_MD2_20230824_115135.h5', #QPH=-1.3
'-0.9' : 'SPS.BQHT_MD2_20230824_112211.h5', #QPH=-1.0
'-0.8' : 'SPS.BQHT_MD2_20230824_111017.h5', #QPH=-0.8
'-0.7' : 'SPS.BQHT_MD2_20230824_110811.h5', #QPH=-0.7
'-0.6' : 'SPS.BQHT_MD2_20230824_110359.h5', #QPH=-0.6
'-0.55' : 'SPS.BQHT_MD2_20230824_105823.h5', #QPH=-0.55
'-0.4' : 'SPS.BQHT_MD2_20230824_105247.h5', #QPH=-0.4
'-0.35' : 'SPS.BQHT_MD2_20230824_104917.h5', #QPH=-0.35
'-0.2' : 'SPS.BQHT_MD2_20230824_103229.h5', #QPH=-0.2 
'-0.05' : 'SPS.BQHT_MD2_20230704_113250.h5', #QPH=0.05
}

#analysis
for knob in files:    

    file = path + files[knob]

    print(f'Analyzing file {file.split("/")[-1]}')

    #analysis
    ht = HTchromaticity(file=file, nturns=300, plane='H', knob=knob,  bl=3)
    ht.calc_chromaticity()
    ht.plot()
    
    inp = input('Save? y/n:')
    if inp == 'y': 
        ht.save(fname=f'HT_offline.csv')



