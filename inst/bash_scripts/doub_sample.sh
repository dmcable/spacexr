#!/bin/sh
source /broad/software/scripts/useuse
reuse -q R-3.5
cd /YOUR_PATH/RCTD/
set -e
mkdir -p logs
mkdir -p logs/callDoublets/
Rscript R_scripts/callDoublets.R ${SGE_TASK_ID} > logs/callDoublets/chunk-${SGE_TASK_ID}.txt
