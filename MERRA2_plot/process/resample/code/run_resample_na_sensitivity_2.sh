#!/bin/sh
#$ -N na_soil_NOx
#$ -q ARROMA
##$ -l h=argon-itf-bx47-37
#$ -cwd
#$ -pe smp 12
#$ -o ./$JOB_ID.out
#$ -e ./$JOB_ID.err

ulimit -s unlimited
export KMP_STACKSIZE=209715200
cd /Dedicated/jwang-data/ywang/soil_NOx/MERRA2_plot/process/resample/code/tmp_na_sensitivity_2
time python /Dedicated/jwang-data/ywang/soil_NOx/MERRA2_plot/process/resample/code/main_resample_na_sensitivity_2.py startDate endDate
