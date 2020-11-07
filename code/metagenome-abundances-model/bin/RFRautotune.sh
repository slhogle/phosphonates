#!/bin/bash -l
#SBATCH -J RFRautotune
#SBATCH -o RFRautotune.log
#SBATCH -A project_2001175
#SBATCH -p small
#SBATCH -t 0-48:00:00
#SBATCH -n 20
#SBATCH --mem-per-cpu=4000

echo "Batch job started at `date`" 

T1=`date +%s`

module load r-env

Rscript --vanilla RFRautotune.R

T2=`date +%s`

TDIFF=`echo 'scale=3;('$T2'-'$T1')/60' | bc`

echo "Total runtime: ${TDIFF} mins"