#!/bin/bash


# define global shortcuts
export EOS_MGM_URL=root://eosuser.cern.ch
export myfolder=/afs/cern.ch/work/m/miguelgo/private/pyhead/match_10_07_2024 #folder where the program is
export myfilename="sps_singlebunchQ26" #name of the program without extension
cd $myfolder # move htcondor session to the folder where we are working
# copy file from eos
# cp $mycernboxfolder/$myfilename.py $myfolder/.

# change permission to open the file
chmod +x $myfolder/$myfilename.py

# launch program, remember to put the correct location of the python installation
/afs/cern.ch/work/m/miguelgo/miniconda3/bin/python $myfolder/$myfilename.py 

# User notes:
# -------------
# Results will be in the Export/ folder only if the post-processing step General1D>ASCIIexport>Export All 1D results is checked in the .cst file
# numthreads should be equal to the number of cpu requested in the sub file
# Before submitting, change permission on .sh and .sub files with chmod -x
# and try to run the .sh file by doing: . cst.sh
# Before submitting, change to the least used batch: myschedd bump
# To submit the file to htcondor do: condor_submit cst.sub
# To check the status of the htcondot job do: condor_q or: cat $filename/Result/progress.log
# To monitor the status of the job: watch condor_q or condor_tail JOBID (prints the .out file)
# To cancel a job do: condor_rm JOBID
# To access jobs local files while running: condor_ssh_to_job JOBID
# To know the total running time when it finishes: condor_history JOBID -limit 1 -af:h RemoteWallClockTime