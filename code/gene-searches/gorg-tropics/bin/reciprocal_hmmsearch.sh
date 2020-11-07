#!/usr/bin/env bash

# Searches HMM models from:
# https://github.com/slhogle/phosphonates/blob/master/data/phosphonate_biosynthesis/HMM_models
# https://github.com/slhogle/phosphonates/blob/master/data/phosphonate_catabolism/HMM_models
# https://github.com/slhogle/phosphonates/blob/master/data/core_gene_families/HMM_models 
# against a collection of proteins (eg. GORG-tropics). Pulls best hits, then searches them 
# against the whole PFAM and TIGRFAM to make sure they are reciprocal best hits

### ARGUMENTS
### ${1} = PHNK, PEPM, MPNS, etc, etc...

MYFAM=${1}

HMMDIR="/nobackup1/shogle/projects/phosphonates/hmms/phn-utilization"
#phn-utilization
#phn-synthesis

HMMDBDIR="/nobackup1/shogle/seq_databases/hmms"
SEARCHPATH="/nobackup1/shogle/seq_databases/raw_sequence_files/GORG-TROPICS"

echo "performing initial hmmsearch..."

hmmsearch \
  --noali \
  -o "${MYFAM}-01.tmp" \
  --tblout "${MYFAM}-01.hmmsearch" \
  --cpu 20 \
  --cut_ga \
  "${HMMDIR}/${MYFAM}/${MYFAM}.hmm" "${SEARCHPATH}/GORG-TROPICS-genes.faa"

rm "${MYFAM}-01.tmp"

cat "${SEARCHPATH}/GORG-TROPICS-genes.faa" | \
seqkit grep -n -f <(cat ${MYFAM}-01.hmmsearch | bioawk '{ print $1 }' | sort -u | grep -v "^#"
) > "${MYFAM}-01.faa"

# pipe in file
# exclude lines starting with #
# replace one or more whitespaces with tabs
# select useful columns with bioawk
# sort first on query, then bitscore, then evalue - 1=query, 3=pfam_desci, 4=pfam_id, 5=evalue, 6=bitscore

cat "${MYFAM}-01.hmmsearch" |
grep -v "^#" |
sed -r 's/\s+/\t/g' |
bioawk -t '{ print $1,$4,$6,$5,$9,$8 }' |
sort -k1,1 -k3,3rn -k4,4g > ${MYFAM}-01-formatted.tsv

echo "performing reciprocal search for ${MYFAM}..."

hmmsearch \
  --noali \
  -o "${MYFAM}-02.tmp" \
  --cut_ga \
  --tblout "${MYFAM}-02.hmmsearch" \
  --cpu 20 \
  "${HMMDBDIR}/database.hmm" "${MYFAM}-01.faa"

rm "${MYFAM}-02.tmp"

# pipe in file
# exclude lines starting with #
# replace one or more whitespaces with tabs
# select useful columns with bioawk
# sort first on query, then bitscore, then evalue - 1=query, 3=pfam_desci, 4=pfam_id, 5=evalue, 6=bitscore

cat "${MYFAM}-02.hmmsearch" |
grep -v "^#" |
sed -r 's/\s+/\t/g' |
bioawk -t '{ print $1,$4,$6,$5,$9,$8 }' |
sort -k1,1 -k3,3rn -k4,4g > ${MYFAM}-02-formatted.tsv

cat ${MYFAM}-02-formatted.tsv |
sort -u -k1,1 --merge |
rg -F -N -f ${HMMDIR}/${MYFAM}/${MYFAM}.list > ${MYFAM}-reciprocalbesthits.tsv

#awk '{ if ($5 < 1e-30) { print } }' # only print hits if evalue is less than 1e-30 (for removing weak hits)

cat "${SEARCHPATH}/GORG-TROPICS-genes.faa" |
seqkit grep -n -f <(awk 'FS="\t" { print $1 } ' ${MYFAM}-reciprocalbesthits.tsv | sort -u) > ${MYFAM}-reciprocalbesthits.faa

rm "${MYFAM}-01.hmmsearch"
rm "${MYFAM}-01.faa"
rm "${MYFAM}-02.hmmsearch"