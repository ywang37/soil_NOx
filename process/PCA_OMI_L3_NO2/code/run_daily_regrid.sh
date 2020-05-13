#!/bin/sh
#$ -N regrid
#$ -q ARROMA
#$ -cwd
#$ -pe smp 6
#$ -o ./$JOB_ID.out
#$ -e ./$JOB_ID.err

ulimit -s unlimited
export KMP_STACKSIZE=209715200
cd /Dedicated/jwang-data/ywang/soil_NOx/process/PCA/code/tmp
time python /Dedicated/jwang-data/ywang/soil_NOx/process/PCA/code/main_daily_regrid.py startDate endDate
