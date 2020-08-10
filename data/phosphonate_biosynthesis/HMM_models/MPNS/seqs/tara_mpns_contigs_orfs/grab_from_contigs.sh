#!/usr/bin/env bash

BASE="/nobackup1/shogle/projects/phosphonates/MPNS_custom_hmm"
FAM="mpns"

while read LINE; do ASSBASE=$(echo ${LINE} | cut -c1-5); zcat ${BASE}/tara_assemblies/${ASSBASE}*.fasta.gz | seqkit grep -r -p ${LINE} | seqkit seq -i | orfm >> "all_${FAM}_orfs.faa"; done < "tara-${FAM}.ids"

makeblastdb -in "all_${FAM}_orfs.faa" -dbtype prot -out "all_${FAM}_orfs"
blastp -db "all_${FAM}_orfs" -query "${BASE}/ncbi-${FAM}.faa" -evalue 1e-20 -outfmt 6 -out tara_${FAM}_orfs-v-ncbi_${FAM}.tsv

cut -f2 tara_${FAM}_orfs-v-ncbi_${FAM}.tsv | sort -u > tara_${FAM}_verified_orfs_ids.txt

cat "all_${FAM}_orfs.faa" | seqkit grep -n -f "tara_${FAM}_verified_orfs_ids.txt" | seqkit seq -i --id-regexp "^ENA\|.+\|(.+)" > tara_${FAM}.faa

rm all_${FAM}_orfs.*
rm tara_${FAM}_orfs-v-ncbi_${FAM}.tsv
rm tara_${FAM}_verified_orfs_ids.txt
