#!/bin/bash
JOBONE=$(qsub -pe smp 8 -R y -binding linear:8 -l h_rt=48:00:00 -l h_vmem=2G -l os=RedHat7 /YOUR_PATH/RCTD/R_scripts/pre_sample.sh)
jobid=$(echo $JOBONE | cut -d " " -f 3)
qsub -hold_jid $jobid -t 1:NUM_FOLDS -pe smp 8 -R y -binding linear:8 -l h_rt=48:00:00 -l h_vmem=2G -l os=RedHat7 /YOUR_PATH/RCTD/R_scripts/doub_sample.sh
