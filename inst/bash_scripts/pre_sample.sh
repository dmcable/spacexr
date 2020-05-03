#!/bin/sh
source /broad/software/scripts/useuse
reuse -q R-3.5
cd /YOUR_PATH/RCTD/
set -e
mkdir -p logs
Rscript fitBulk.R > logs/fitBulk.txt
Rscript chooseSigma.R > logs/chooseSigma.txt
