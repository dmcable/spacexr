#!/bin/sh
set -e
mkdir -p Data
mkdir -p Data/Slideseq
mkdir -p Data/Reference
mkdir -p logs
Rscript prepareBulkData.R > logs/prepareBulkData.txt
python fitBulkWeights.py > logs/fitBulkWeights.txt
Rscript callBeads.R > logs/callBeads.txt
