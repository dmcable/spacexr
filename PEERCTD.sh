#!/bin/sh
mkdir logs
Rscript prepareBulkData.R > logs/prepareBulkData.txt
python fitBulkWeights.py > logs/fitBulkWeights.txt
Rscript callBeads.R > logs/callBeads.txt
