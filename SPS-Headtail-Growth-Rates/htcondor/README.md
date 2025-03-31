## HT Condor
> Example files to submit a PyHEADTAIL job to HTCONDOR from LXPlus

---

### User notes:

- `numthreads` should be equal to the number of cpu requested in the sub file
- If you find permission problems for submitting, change permission on .sh and .sub files with chmod -x and try to run the .sh file on you LXPlus node by doing: . simulation.sh
- Before submitting, change to the least used batch: `myschedd bump`
- To submit the file to htcondor do: `condor_submit simulation.sub`
- To check the status of the htcondot job do: condor_q 
- To monitor the status of the job: `watch condor_q` or `watch condor_tail JOBID` (prints the .out file). Exit by pressing `q`
- To cancel a job do: `condor_rm JOBID`
- To access jobs local files while running: `condor_ssh_to_job JOBID`
- To know the total running time when it finishes: `condor_history JOBID -limit 1 -af:h RemoteWallClockTime`

