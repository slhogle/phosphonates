#!/usr/bin/env bash
#SBATCH -J diamond
#SBATCH --time 0-15:00:00
#SBATCH -o ../logs/tara/%j__ST%a.log
#SBATCH -p sched_mit_chisholm
#SBATCH --mem 200GB
#SBATCH -n 20

T1=`date +%s`

module load engaging/zlib/1.2.8

MYLIB=tara_S`printf %04d $SLURM_ARRAY_TASK_ID`
MYMETA="tara"
THREADS=20

READDIR="/nfs/chisholmlab001/chisholmlab/tara-shogle/qc-reads"
OUTDIR="../raw_output/${MYMETA}"
UUID=$SLURM_JOB_ID-$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 16 | head -n 1)

mkdir ${OUTDIR}/${UUID}
mkdir ${OUTDIR}/${UUID}/dmd

echo $MYLIB, $SLURM_JOB_ID

rename.sh in=${READDIR}/${MYLIB}_overlappedandsingletons_qctrimmed.fa.gz out=${OUTDIR}/${UUID}/${MYLIB}-SGL.fa prefix=SGL

rename.sh in=${READDIR}/${MYLIB}_nonoverlapped_qctrimmed_1.fa.gz in2=${READDIR}/${MYLIB}_nonoverlapped_qctrimmed_2.fa.gz out=${OUTDIR}/${UUID}/${MYLIB}-PRD.fa prefix=PRD

cat ${OUTDIR}/${UUID}/${MYLIB}-SGL.fa ${OUTDIR}/${UUID}/${MYLIB}-PRD.fa > ${OUTDIR}/${UUID}/${MYLIB}-ALL.fa
rm ${OUTDIR}/${UUID}/${MYLIB}-SGL.fa
rm ${OUTDIR}/${UUID}/${MYLIB}-PRD.fa

diamond blastx \
-p ${THREADS} \
-f 6 qseqid sseqid pident evalue bitscore qcovhsp qlen slen length qstart qend qframe \
-c 1 \
-b 20 \
--top 5 \
-e 1e-5 \
--id 55 \
-t ${OUTDIR}/${UUID}/dmd \
-d ../db/PEPM.dmnd \
-q ${OUTDIR}/${UUID}/${MYLIB}-ALL.fa \
-o ${OUTDIR}/${UUID}/${MYLIB}.m8 \
--sensitive

rg -N "^SGL" ${OUTDIR}/${UUID}/${MYLIB}.m8 > ${OUTDIR}/${UUID}/${MYLIB}-SGL.m8

rg -N "^PRD" ${OUTDIR}/${UUID}/${MYLIB}.m8 | sed 's/\s[1-2]:$//g' | sort -k1,1 -k5,5rn -k3,3g > ${OUTDIR}/${UUID}/${MYLIB}-PRD.m8

cat ${OUTDIR}/${UUID}/${MYLIB}-SGL.m8 ${OUTDIR}/${UUID}/${MYLIB}-PRD.m8 > ${OUTDIR}/${MYLIB}.m8


if [ -s ${OUTDIR}/${MYLIB}.m8 ]
then
	rm -r ${OUTDIR}/${UUID}
fi

T2=`date +%s`

TDIFF=`echo 'scale=3;('$T2'-'$T1')/60' | bc`
echo "${MYLIB} runtime: ${TDIFF} mins"
