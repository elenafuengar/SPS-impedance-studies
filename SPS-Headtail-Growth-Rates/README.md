## SPS HeadTail mode 0
> Scripts needed to run a PyHEADTAIL simulation of the HeadTail mode 0 instability and post-process the results

- sps_singlebunchQ26.py: PyHEADTAIL simulation script defining the accelerator optics and tracking
- SPSoctupolesNewConfiguration.py: class to add Octupole effects to high-order chromaticity
- growthrate.py: script containing the class GrowthRate that post-process the PyHeadtail data and give the instability GrowthRate
- qpvscan.py: Example of post-processing script using Growthrate class and PyHeadtail Bunch Monitor .h5 output files
- SPS_Q26.wake: ASCII table containin the SPS Wake model for q26 optics, generated from: https://gitlab.cern.ch/IRIS/SPS_IW_model/SPS_IW_merged_SingleMulti_bunch_model

### SPS parameters
Simulations use single bunch low-intensity beam with Q26 optics. 
SPS machine parameters for Q20 and Q26 optics can be found in: 
* [H. Bartosik PhD Thesis](http://cds.cern.ch/record/1644761/files/CERN-THESIS-2013-257.pdf?version=1)

![](https://codimd.web.cern.ch/uploads/upload_865d76bd94f7bb6626ab09fc454c878f.png)

![](https://codimd.web.cern.ch/uploads/upload_1a18ca53b55a113b0be0116e7b4c4a86.png)

### Quick PyHeadtail installation guide:
This guide works for single-bunch PyHEADTAIL simulations. For multi-bunch simulations, we recommend using `spack` and run on `SLURM` HPC cluster.
Installing a Miniconda environment on `/afs/cern.ch/work/` from LxPlus and submit the simulation to `HTCondor`
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
