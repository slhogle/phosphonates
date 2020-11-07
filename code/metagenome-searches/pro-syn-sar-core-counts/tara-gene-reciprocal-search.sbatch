#!/usr/bin/env bash
#SBATCH -J diamond
#SBATCH --time 0-10:00:00
#SBATCH -o /nobackup1/shogle/projects/pro-syn-sar-core-counts/logs/tara/%j__ST%a.log
#SBATCH -p newnodes
#SBATCH --mem 140GB
#SBATCH -n 20

module load engaging/zlib/1.2.8

MYLIB=tara_S`printf %04d $SLURM_ARRAY_TASK_ID`

T1=`date +%s`

echo $MYLIB, $SLURM_JOB_ID

############################################################################################################################
############################################### PATHS AND GLOBAL VARIABLES #################################################
############################################################################################################################
MYMETA="tara"
THREADS="20"

MARKERDIR="/nobackup1/shogle/seq_databases/diamond/MARMICRODB/seq_families/pro-syn-pel-kaks10"
DIMNDMARMICRODB="/nobackup1/shogle/seq_databases/diamond/MARMICRODB/MARMICRODB.dmnd"
READDIR="/nfs/chisholmlab001/chisholmlab/tara-shogle/qc-reads"
OUTDIR="/nobackup1/shogle/projects/pro-syn-sar-core-counts/raw_output/${MYMETA}"

## Directory structure for search output
MYTMP="${OUTDIR}/${MYLIB}"
mkdir ${MYTMP}
MYDMDTMP="${MYTMP}/dmd-tmp"
mkdir ${MYDMDTMP}
MY1STSRCH="${MYTMP}/SRCH1"
mkdir ${MY1STSRCH}
MY2NDSRCH="${MYTMP}/SRCH2"
mkdir ${MY2NDSRCH}

############################################################################################################################
################################# 1ST SEARCH: READS VS ALL IDENTFIED RIBOSOMAL PROTEINS ####################################
############################################################################################################################

rename.sh in=${READDIR}/${MYLIB}_overlappedandsingletons_qctrimmed.fa.gz out=${MYTMP}/${MYLIB}-SGL.fa prefix=SGL 

rename.sh in1=${READDIR}/${MYLIB}_nonoverlapped_qctrimmed_1.fa.gz in2=${READDIR}/${MYLIB}_nonoverlapped_qctrimmed_2.fa.gz out=${MYTMP}/${MYLIB}-PRD.fa prefix=PRD 

cat ${MYTMP}/${MYLIB}-SGL.fa ${MYTMP}/${MYLIB}-PRD.fa > ${MY1STSRCH}/${MYLIB}-ALL.fa

diamond blastx \
	-p ${THREADS} \
	-f 6 qseqid sseqid bitscore evalue pident qcovhsp slen qlen length \
	-c 1 \
	-b 20 \
	--top 5 \
	-e 0.00001 \
	-t ${MYDMDTMP} \
	-d ${MARKERDIR}/ALLMARKERS/ALLMARKERS.dmnd \
	-q ${MY1STSRCH}/${MYLIB}-ALL.fa \
	-o ${MY1STSRCH}/${MYLIB}.m6

cut -f1 ${MY1STSRCH}/${MYLIB}.m6 | sort -u > ${MY1STSRCH}/${MYLIB}-ALL.names

rg -N "^SGL" ${MY1STSRCH}/${MYLIB}-ALL.names > ${MY1STSRCH}/${MYLIB}-SGL.names

rg -N "^PRD" ${MY1STSRCH}/${MYLIB}-ALL.names | sed 's/ [1-2]:$//g' | sort -u | sed 's/$/ 1:/g' > ${MY1STSRCH}/${MYLIB}-PRD.names

filterbyname.sh \
	include=t \
	in=${MYTMP}/${MYLIB}-SGL.fa \
	names=${MY1STSRCH}/${MYLIB}-SGL.names \
	out=${MY1STSRCH}/${MYLIB}-SGL-1STSEARCH.fa

filterbyname.sh \
	include=t \
	in=${MYTMP}/${MYLIB}-PRD.fa \
	names=${MY1STSRCH}/${MYLIB}-PRD.names \
	out=${MY1STSRCH}/${MYLIB}-PRD-1STSEARCH.fa

cat ${MY1STSRCH}/${MYLIB}-SGL-1STSEARCH.fa ${MY1STSRCH}/${MYLIB}-PRD-1STSEARCH.fa > ${MY2NDSRCH}/${MYLIB}-ALL-2NDSEARCH.fa

rm -r ${MY1STSRCH}

############################################################################################################################
######################## 2ND SEARCH: DIAMOND RECIPROCAL SEARCH TO FILTER READS HITTING ALL MARKERS #########################
############################################################################################################################

diamond blastx \
	-p ${THREADS} \
	-f 6 qseqid sseqid bitscore evalue pident qcovhsp slen qlen length \
	--top 5 \
	-c 1 \
	-b 20 \
	-e 0.00001 \
	-t ${MYDMDTMP} \
	-d ${DIMNDMARMICRODB} \
	-q ${MY2NDSRCH}/${MYLIB}-ALL-2NDSEARCH.fa \
	-o ${MY2NDSRCH}/${MYLIB}.m6

rg -N "^SGL" ${MY2NDSRCH}/${MYLIB}.m6 > ${MY2NDSRCH}/${MYLIB}-SGL.m6

rg -N "^PRD" ${MY2NDSRCH}/${MYLIB}.m6 | sed 's/\s[1-2]:$//g' | sort -k1,1 -k3,3rn -k4,4g > ${MY2NDSRCH}/${MYLIB}-PRD.m6

cat ${MY2NDSRCH}/${MYLIB}-SGL.m6 ${MY2NDSRCH}/${MYLIB}-PRD.m6 > ${OUTDIR}/${MYLIB}-ALL-RECIPR.m6

rm -r ${MYTMP}

T2=`date +%s`

TDIFF=`echo 'scale=3;('$T2'-'$T1')/60' | bc`
echo "${MYLIB} runtime: ${TDIFF} mins"
