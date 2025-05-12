'''
Record SPS LSA knobs QPV & QPH,
BBQcont position data 
and BCT intensity data

Aug'24
'''

from pyjapcscout import PyJapcScout
from helper_functions import callback_core
import PlottingClassesSPS as pc

from datetime import date

def outer_callback(plot_func_list):
    def callback(data, h):
        callback_core(data, h, plot_func_list)
    return callback

# Set up monitor
selector = 'SPS.USER.MD3' #MD2 q26, MD5 q22
BBQ_device = 'SPS.BQ.QC/Acquisition'
BBQCONT_device = 'SPS.BQ.KICKED/ContinuousAcquisition'
devices = [BBQ_device, BBQCONT_device,'SPS.BCTDC.51456/Acquisition', 'SPSBEAM/QPV', 'SPSBEAM/QPH']

# start PyJapcScout and so incaify Python instance
myPyJapc = PyJapcScout(incaAcceleratorName='SPS')

myMonitor = myPyJapc.PyJapcScoutMonitor(selector, devices)

# Saving data configuration
date=str(date.today())
myMonitor.saveDataPath = f'./data/GR/{selector}/{date}/'
myMonitor.saveData = True
myMonitor.saveDataFormat = 'parquet'

# Start acquisition
myMonitor.startMonitor()


'''
if 0:
    ## for controlling data acquisition:
    #
    myMonitor.saveData = True
    # myMonitor.saveData = False

    ## to stop the monitor
    #myMonitor.stopMonitor()
'''