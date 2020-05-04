#!/bin/sh
source /broad/software/scripts/useuse
reuse -q R-3.5
cd /YOUR_PATH/RCTD/
set -e
mkdir -p logs
Rscript R_scripts/fitBulk.R > logs/fitBulk.txt
Rscript R_scripts/chooseSigma.R > logs/chooseSigma.txt
