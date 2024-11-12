## SPS HeadTail mode 0
> Scripts needed to run a PyHEADTAIL simulation of the HeadTail mode 0 instability and post-process the results

- sps_singlebunchQ26.py: PyHEADTAIL simulation script defining the accelerator optics and tracking
- SPSoctupolesNewConfiguration.py: class to add Octupole effects to high-order chromaticity
- growthrate.py: script containing the class GrowthRate that post-process the PyHeadtail data and give the instability GrowthRate
- qpvscan.py: Example of post-processing script using Growthrate class and PyHeadtail Bunch Monitor .h5 output files
- SPS_Q26.wake: ASCII table containin the SPS Wake model for q26 optics, generated from: https://gitlab.cern.ch/IRIS/SPS_IW_model/SPS_IW_merged_SingleMulti_bunch_model

