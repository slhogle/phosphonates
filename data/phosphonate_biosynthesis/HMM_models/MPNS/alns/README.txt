# commands to produce hmms

einsi <(cat ../seqs/ncbi-mpns.faa ../seqs/tara_mpns_contigs_orfs/tara_mpns.faa) > mpns-einsi.faa.aln
trimal -in mpns-einsi.faa.aln -out mpns-einsi.faa.gappyout.aln -gappyout

## stupid trimal can't select regions of alignment...
read_fasta -i mpns-einsi.faa.gappyout.aln | slice_align -b 317 -e 792 | write_fasta -x > mpns-einsi.faa.gappyout.select317-792.aln

seqmagick convert --input-format fasta --output-format stockholm mpns-einsi.faa.gappyout.select317-792.aln mpns-einsi.faa.gappyout.select317-792.sto

hmmbuild -n mpns mpns.hmm mpns-einsi.faa.gappyout.select317-792.sto
