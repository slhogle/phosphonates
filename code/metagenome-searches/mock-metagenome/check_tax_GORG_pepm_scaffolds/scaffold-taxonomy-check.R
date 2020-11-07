library(tidyverse)
library(here)

seqs <- read_tsv(here("data", "PEPMscaffs-v-MARMICRODB.m8"), 
                  col_names=c("qseqid", "tseqid", "pident", 
                              "evalue", "bitscore", "qcov", "tcov", 
                              "qlen", "tlen", "alnlen"),  col_types="ccnnnnnnnn")

pepm <- read_tsv(here("data", "pepm-ver-summary-mod.tsv"), 
                 col_names=c("tseqid", "group", "seqtype"),  col_types="ccc") 

seqsTOP <- seqs %>% 
  group_by(qseqid) %>% 
  arrange(desc(bitscore)) %>% 
  filter(bitscore == max(bitscore)) %>%
  ungroup()  

seqsTOPchop <- seqsTOP %>%
  separate(qseqid, c("genome", "scaffold", "gene"), 3,
           sep="_", extra="drop", remove=FALSE) %>%
  separate(tseqid, c("remainder", "taxid"), 2,
           sep="_(?=[:digit:]+$)", extra="merge", remove=FALSE) %>% select(-remainder) %>%
    arrange(scaffold, gene) %>%
    mutate(taxid=as.character(taxid)) %>%
    left_join(., pepm)

write_tsv(seqsTOPchop %>% select(taxid), here("data", "taxids.txt"))

taxids_summary <- read_tsv(here("data", "taxidsids-lineage.tsv"), 
                 col_names=c("taxid", "taxsummary", "taxlineage"),  col_types="ccc")

seqsTOPchop <- left_join(seqsTOPchop, taxids_summary)

taxsum <- seqsTOPchop %>%
    group_by(genome, scaffold) %>% 
    count(taxsummary) %>%
    mutate(seqsum=sum(n)) %>%
    mutate(mostabundanttax_perc=max(n/seqsum)) %>%
    filter(n==max(n))

write_tsv(taxsum, here("data", "taxsummary.tsv"))
