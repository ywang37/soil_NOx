#!/bin/sh
#$ -N soil_NOx
#$ -q INFORMATICS
#$ -cwd
#$ -pe smp 6
#$ -o ./$JOB_ID.out
#$ -e ./$JOB_ID.err

ulimit -s unlimited
export KMP_STACKSIZE=209715200
cd /Dedicated/jwang-data/ywang/soil_NOx/MERRA2_plot/process/resample/code/tmp
time python /Dedicated/jwang-data/ywang/soil_NOx/MERRA2_plot/process/resample/code/main_resample.py startDate endDate
