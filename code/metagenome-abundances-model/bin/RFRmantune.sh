#!/bin/bash -l
#SBATCH -J RFRtune
#SBATCH -o RFRtune.log
#SBATCH -A project_2001175
#SBATCH -p small
#SBATCH -t 0-36:00:00
#SBATCH -n 20
#SBATCH --mem-per-cpu=3000

module load r-env

Rscript --vanilla RFRmantune.R