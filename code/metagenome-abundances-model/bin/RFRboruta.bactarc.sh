#!/bin/bash -l
#SBATCH --job-name RFRboruta
#SBATCH --output RFRboruta.log
#SBATCH --account project_2001175
#SBATCH --partition small
#SBATCH --time 0-48:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 20
#SBATCH --cpus-per-task 1
#SBATCH --mem-per-cpu=6000

echo "Batch job started at `date`" 

T1=`date +%s`

module load r-env

Rscript --vanilla RFRboruta.bactarc.R

T2=`date +%s`

TDIFF=`echo 'scale=3;('$T2'-'$T1')/60' | bc`

echo "Total runtime: ${TDIFF} mins"