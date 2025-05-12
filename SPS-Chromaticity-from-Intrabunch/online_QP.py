import glob
from HTchromaticity import HTchromaticity
from datetime import date

today = '2024_08_02'  
datapath = f'/nfs/cs-ccr-bqhtnfs/sps_data/SPS.BQHT/{today}/'

#For chroma knob:
parquetpath = f'/user/spsscrub/2024/sps_beam_monitoring/sps_beam_monitoring/general/GR/data/GR/SPS.USER.MD2/{date.today().isoformat()}/'

#monitor
while True:

    # For HT
    flist = sorted(glob.glob(datapath+'*.h5'))
    file = flist[-1]
    print(f'Analyzing file {file.split("/")[-1]}')

    # For chroma knob
    parquetlist = sorted(glob.glob(parquetpath+'*.parquet'))
    parquet = parquetlist[-1]
    print(f'Obtaining knob from {parquet.split("/")[-1]}')

    ht = HTchromaticity(file=file, nturns=500, plane='V', 
                        parquet=parquet, bl=3.5, dpi=200, monitor=False)
    ht.calc_chromaticity()
    ht.plot(ylim=(0,25))

    inp = input('Save? y/n:')
    if inp == 'y': 
        ht.save(fname=f'HT_{today}.csv')

    #time.sleep(40) #supercycle time


