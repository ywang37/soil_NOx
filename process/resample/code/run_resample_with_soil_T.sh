#!/bin/sh
#$ -N soil_NOx
#$ -q ARROMA
#$ -cwd
#$ -pe smp 6
#$ -o ./$JOB_ID.out
#$ -e ./$JOB_ID.err

ulimit -s unlimited
export KMP_STACKSIZE=209715200
cd /Dedicated/jwang-data/ywang/soil_NOx/process/resample/code/tmp
time python /Dedicated/jwang-data/ywang/soil_NOx/process/resample/code/main_resample_with_soil_T.py startDate endDate
