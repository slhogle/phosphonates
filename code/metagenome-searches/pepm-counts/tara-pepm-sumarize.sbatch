#!/usr/bin/env bash
#SBATCH -J marker_summarize
#SBATCH --time 0-04:00:00
#SBATCH -o ../logs/reports/%j__tara_ST%a.log
#SBATCH -p sched_mit_chisholm
#SBATCH --mem 16GB
#SBATCH -n 1

MYLIB=tara_S`printf %04d $SLURM_ARRAY_TASK_ID`

MYMETA="tara"

T1=`date +%s`
echo $MYLIB, $SLURM_JOB_ID

python pepm_parser.py ${MYLIB} ${MYMETA}

T2=`date +%s`

TDIFF=`echo 'scale=3;('$T2'-'$T1')/60' | bc`
echo "${MYLIB} runtime: ${TDIFF} mins"
