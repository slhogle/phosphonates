library(tidyverse)
library(here)
library(fs)

readsmrg <- read_tsv(here("data", "readsim-more-sensitive.merged.m8"), 
                     col_names=c("qseqid", "sseqid", "pident", 
                                 "evalue", "bitscore", "qcovhsp", "qlen", 
                                 "slen", "alnlen", "qstart",
                                 "qend", "qframe"),  col_types="ccnnnnnnnnnn")

readssng <- read_tsv(here("data", "readsim-more-sensitive.single.m8"), 
                     col_names=c("qseqid", "sseqid", "pident", 
                                 "evalue", "bitscore", "qcovhsp", "qlen", 
                                 "slen", "alnlen", "qstart",
                                 "qend", "qframe"),  col_types="ccnnnnnnnnnn")

reads <- bind_rows(readsmrg, readssng)

#pepm <- read_tsv(here("data", "pepm-summary.tsv"), 
#                 col_names=c("sseqid", "group", "seqtype"),  col_types="ccc") 

#### changed tax affiliation of TMED122, TMED191, TMED37,
#### ARS1441, and NAT185 to pelagibacterales

pepm <- read_tsv(here("data", "pepm-ver-summary-alttax.tsv"), 
                 col_names=c("sseqid", "group", "seqtype"),  col_types="ccc") 

gorgtaxonomy <- read_tsv(here("data", "GORG-taxonomy.tsv"), 
                 col_names=c("genome", "qseqtax"),  col_types="cc")

simulation_summary <- tibble::tribble(
  ~qseqlen,      ~source, ~qseqtype, ~readcount,  ~totalbp,
   "mergedreads", "backgroundcommunity",   "other",    1523820, 354372923,
   "singleread", "backgroundcommunity",   "other",     243784,  36567600,
   "mergedreads",     "prochlorococcus",       "PEPM",       2751,    638939,
   "singleread",     "prochlorococcus",       "PEPM",        420,     63000,
   "mergedreads",    "pelagibacterales",       "PEPM",      10960,   2546749,
   "singleread",    "pelagibacterales",       "PEPM",       1868,    280200,
   "mergedreads",    "pelagibacterales",  "falsePEPM",      11501,   2674579,
   "singleread",    "pelagibacterales",  "falsePEPM",       1838,    275700,
   "mergedreads",           "othertaxa",       "PEPM",       1648,    382832,
   "singleread",           "othertaxa",       "PEPM",        248,     37200,
   "mergedreads",           "othertaxa",  "falsePEPM",       4519,   1052100,
   "singleread",           "othertaxa",  "falsePEPM",        686,    102900
  )

reads1 <- left_join(reads, pepm) %>%
  rename(sseqtax=group) %>% select(-seqtype) %>%
  separate(qseqid, c("readtax", "readsource", "pepmverified", "readoverlap", "remainder"), 5,
           sep="\\|", extra="drop", remove=FALSE) %>% 
  mutate(qseqsource=ifelse(readsource=="GRG", "GORG", ifelse(readsource=="GEO", "geotraces_assemblies", ifelse(readsource=="RGC", "TARA-OM-RGC", readsource)))) %>%
  mutate(qseqtype=ifelse(pepmverified=="TRU", "PEPM", ifelse(pepmverified=="FLS", "falsePEPM", ifelse(pepmverified=="NAN", "other", pepmverified)))) %>%
  mutate(qseqlen=ifelse(readoverlap=="MRG", "mergedreads", ifelse(readoverlap=="SNG", "singleread", readoverlap))) %>%
  separate(remainder, c("junk", "seqnum", "strand", "start", "end", "genome", "scaffold", "gene"), 8,
           sep="_", extra="drop", remove=TRUE) %>% select(-junk, -readtax, -readsource, -pepmverified, -readoverlap) %>%
    left_join(., gorgtaxonomy) %>%
    mutate(qseqtax=ifelse(is.na(qseqtax), "prochlorococcus", qseqtax))

reads1TOP <- reads1 %>% 
  group_by(qseqid) %>% 
  arrange(evalue) %>% 
  filter(evalue == min(evalue)) %>%
  ungroup()  

###########################################  
########### homology sens/spec ############
###########################################

#### all groups combined ##################

#### distribution no filtering all combined
p1 <- reads1TOP %>% 
  group_by(qseqid) %>% 
  slice(1) %>%
  ungroup() %>%
  mutate(log10evalue=-log10(evalue)) %>%
  gather(param, value, log10evalue, pident) %>%
  mutate(qseqtype=factor(qseqtype, levels=c("PEPM", "falsePEPM", "other"))) %>%
  ggplot() +
  geom_histogram(aes(x=value), bins=50) +
  labs(x="", y="Read Count") + 
  facet_grid(qseqtype ~ param, scales="free") +
  theme_bw()

ggsave(p1, here("figs", "homology_histogram.png"), height=5, width=5)

#### calculate sensitivity and specificity
#### see here: https://en.wikipedia.org/wiki/Sensitivity_and_specificity
#### see here: https://en.wikipedia.org/wiki/F1_score
#### TP = True positive, PEPM reads correctly identified
#### TN = True Negative, other+falsePEPM reads correctly rejected
#### FP = False positive, other+falsePEPM reads incorrectly identified
#### FN = False negative, PEPM reads incorrectly rejected
#### sensitivity = recall = TP / (TP + FN)
#### specificity = TN / (TN + FP)
#### precision = TP / (TP + FP)
#### F1 = 2 * (precision * recall) / (precision + recall)

dat <- data.frame()           
for (e in c(1, 1e-5, 1e-10, 1e-15, 1e-20)){
 for (p in c(c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60)
 )){
   t <- reads1TOP %>%
     filter(evalue <= e) %>%
     filter(pident >= p) %>%
     group_by(qseqid) %>% 
     slice(1) %>%
     ungroup() %>%
     group_by(qseqtype, qseqlen) %>%
     count() %>% rename(hits=n) %>%
     full_join(., simulation_summary %>%
                 group_by(qseqlen, qseqtype) %>%
                 summarize(reads=sum(readcount))) %>%
     replace_na(list(hits = 0, reads = 0)) %>%
     gather(obstype, count, reads, hits) %>%
     group_by(qseqlen) %>%
     unite(newvar, qseqtype, obstype) %>%
     spread(newvar, count, fill=0) %>%
     mutate(TP=PEPM_hits,
            FP=other_hits+falsePEPM_hits,
            TN=(other_reads-other_hits)+(falsePEPM_reads-falsePEPM_hits),
            FN=PEPM_reads-PEPM_hits
     ) %>%
     mutate(sensitivity=TP/(TP+FN),
            specificity=TN/(TN+FP),
            precision=TP/(TP+FP),
            f1=2*(precision*sensitivity)/(precision+sensitivity)) %>%
     select(qseqlen, sensitivity, specificity, precision, f1) %>% 
     mutate(evalue=e, pident=p)
   
   dat <- bind_rows(dat, t)
 }
}

#### plot of sensitivity and specificticy and f1 score
p2 <- dat %>%
  gather(p, v, sensitivity, precision, f1) %>%
  ggplot() +
  geom_line(aes(x=pident, y=v, group=as.factor(evalue), color=as.factor(evalue))) +
  geom_point(data=dat %>% gather(p, v, sensitivity, precision, f1) %>% filter(evalue==1e-5), 
             aes(x=pident, y=v), color="black") +
  geom_line(data=dat %>% gather(p, v, sensitivity, precision, f1) %>% filter(evalue==1e-5), 
             aes(x=pident, y=v), color="black") +
  geom_vline(xintercept=55) +
  labs(y="", x="Percent Id.", color="E value") +
  facet_wrap(qseqlen ~ p, scales="free") +
  theme_bw() 

ggsave(here::here("figs", "sens_precis.svg"), plot=p2, device="svg", units="cm", height=10, width=17)

#### plot of histograms with homology filtering
p3 <- reads1TOP %>% 
  mutate(filtered=ifelse(evalue<=1e-5 & pident >=55, "no", "yes")) %>%
  group_by(qseqid) %>% 
  slice(1) %>%
  ungroup() %>%
  mutate(evalue=-log10(evalue)) %>%
  mutate(qseqtype=factor(qseqtype, levels=c("PEPM", "falsePEPM", "other"))) %>%
  gather(param, value, evalue, pident) %>%
  ggplot() +
  geom_histogram(aes(x=value, fill=filtered), bins=50) +
  facet_grid(qseqtype ~ param, scales="free") +
  labs(x="", y="Read Count", fill="") + 
  theme_bw()

############################################ 
library(patchwork)

p4 <- (p1 / p3) + plot_annotation(tag_levels = 'A')
ggsave(here::here("figs", "homology_histograms.svg"), plot=p4, device="svg", units="cm", height=17, width=14)


#############################################  
############ taxonomy sens/spec #############
#############################################

taxonomyprecision1 <- reads1TOP %>% 
  filter(evalue<=1e-5 & pident >=55) %>%
  mutate(prorecall=ifelse(qseqtype=="PEPM" & qseqtax=="prochlorococcus" & sseqtax=="prochlorococcus", "TP", NA)) %>%
  mutate(prorecall=ifelse(qseqtype=="PEPM" & qseqtax!="prochlorococcus" & sseqtax=="prochlorococcus", "FP", prorecall)) %>%
  mutate(prorecall=ifelse(qseqtype=="PEPM" & qseqtax=="prochlorococcus" & sseqtax!="prochlorococcus", "FN", prorecall)) %>%
  mutate(prorecall=ifelse(qseqtype=="PEPM" & qseqtax!="prochlorococcus" & sseqtax!="prochlorococcus", "TN", prorecall)) %>%
  
  mutate(pelrecall=ifelse(qseqtype=="PEPM" & qseqtax=="pelagibacterales" & sseqtax=="pelagibacterales", "TP", NA)) %>%
  mutate(pelrecall=ifelse(qseqtype=="PEPM" & qseqtax!="pelagibacterales" & sseqtax=="pelagibacterales", "FP", pelrecall)) %>%
  mutate(pelrecall=ifelse(qseqtype=="PEPM" & qseqtax=="pelagibacterales" & sseqtax!="pelagibacterales", "FN", pelrecall)) %>%
  mutate(pelrecall=ifelse(qseqtype=="PEPM" & qseqtax!="pelagibacterales" & sseqtax!="pelagibacterales", "TN", pelrecall))
  
prorecall <- taxonomyprecision1 %>% 
  count(prorecall) %>%
  spread(prorecall, n) %>%
  mutate(sensitivity=TP/(TP+FN)) %>%
  mutate(specificity=TN/(TN+FP)) %>%
  mutate(precision=TP/(TP+FP)) %>%
  select(-`<NA>`) %>%
  mutate(group="prochlorococcus")

pelrecall <- taxonomyprecision1 %>% 
  count(pelrecall) %>%
  spread(pelrecall, n) %>%
  mutate(sensitivity=TP/(TP+FN)) %>%
  mutate(specificity=TN/(TN+FP)) %>%
  mutate(precision=TP/(TP+FP)) %>%
  select(-`<NA>`) %>%
  mutate(group="pelagibacterales")

totalrecall <- bind_rows(prorecall, pelrecall)

write_tsv(totalrecall, here("tables", "tax-recall-alternate.tsv"))
#############################################

