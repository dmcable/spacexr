#!/bin/bash
JOBONE=$(qsub -pe smp 8 -R y -binding linear:8 -l h_rt=48:00:00 -l h_vmem=2G -l os=RedHat7 /YOUR_PATH/RCTD/pre.sh)
jobid=$(echo $JOBONE | cut -d " " -f 3)
qsub -hold_jid $jobid -t 1:40 -pe smp 8 -R y -binding linear:8 -l h_rt=48:00:00 -l h_vmem=2G -l os=RedHat7 /broad/thechenlab/Dylan/slideseq/Cell\ Demixing/RCTD/doub.sh
