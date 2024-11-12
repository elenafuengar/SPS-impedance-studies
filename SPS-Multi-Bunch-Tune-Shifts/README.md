## SPS multi-bunch tune shift and stability studies
> Scripts needed to run a PyHEADTAIL simulation for tune shift and stability studies

- sps_multi_bunch_4x72bunches.py: PyHEADTAIL simulation script defining the accelerator optics and tracking
- SPSoctupolesNewConfiguration.py: class to add Octupole effects to high-order chromaticity
- SPS_wake_model_with_EF_BPH_BPV_kicker_wall_steps_ZS_PostLS2_Q20.txt: table containin the SPS Wake model for q20 optics, generated from: https://gitlab.cern.ch/IRIS/SPS_IW_model/SPS_IW_merged_SingleMulti_bunch_model
- submit_job_sps_multi_bunch_4x72bunches: example of slurm submition script

### SPS parameters
Simulations use multi-bunch high-intensity beam with Q20 optics. 
SPS machine parameters for Q20 and Q26 optics can be found in: 
* [H. Bartosik PhD Thesis](http://cds.cern.ch/record/1644761/files/CERN-THESIS-2013-257.pdf?version=1)

![](https://codimd.web.cern.ch/uploads/upload_865d76bd94f7bb6626ab09fc454c878f.png)

![](https://codimd.web.cern.ch/uploads/upload_1a18ca53b55a113b0be0116e7b4c4a86.png)

### Installation guide:
For multi-bunch simulations, we recommend using `spack` and run on `SLURM` HPC cluster.
The branch used for multi-bunch PyHEADTAIL simulations is `release/v1.17.1`.
This [guide](https://codimd.web.cern.ch/t8xXDI9CSpquoZM3OpKi5w?view) by L. Giacomel explains how to install spack and PyHEADTAIL on `SLURM` HPC cluster.
