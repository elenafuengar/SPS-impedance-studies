## SPS HeadTail mode 0
> Scripts needed to run a PyHEADTAIL simulation of the HeadTail mode 0 instability and post-process the results

- sps_singlebunchQ26.py: PyHEADTAIL simulation script defining the accelerator optics and tracking
- SPSoctupolesNewConfiguration.py: class to add Octupole effects to high-order chromaticity
- growthrate.py: script containing the class GrowthRate that post-process the PyHeadtail data and give the instability GrowthRate
- qpvscan.py: Example of post-processing script using Growthrate class and PyHeadtail Bunch Monitor .h5 output files
- SPS_Q26.wake: ASCII table containin the SPS Wake model for q26 optics, generated from: https://gitlab.cern.ch/IRIS/SPS_IW_model/SPS_IW_merged_SingleMulti_bunch_model

### Quick PyHeadtail installation guide:
At CERN, we recommend installing a Miniconda environment on `/afs/cern.ch/work/` from LxPlus and submit the simulation to `HTCondor`
```
# get, install and activate miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh 
source miniconda3/bin/activate

# create dev python environment python>3.9
conda create --name pyHT python=3.9
conda activate pyHT

# Install useful packages [optional]:
pip install ipython matplotlib scipy numpy pandas h5py
```
For using PyHEADTAIL without modifying the source code, we recommend to install the latest version via PyPI:
```
pip install PyHEADTAIL
```
To run PyHEADTAIL scripts, `ipython` run is recommended:
```
(pyHT)user/path$ ipython
(...)
In [1] run pyheadtail_script.py
```
