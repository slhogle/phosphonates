library(tidyverse)
library(gggenes)
library(svglite)
library(here)

all_coord <- read_csv(here("data", "phn-syn-prok_coords.csv")) %>%
  mutate(strand=ifelse(strand=="-", -1, 
                       ifelse(strand=="+", 1, strand)))

# this is a list of prokaryotic PEPM gene_ids that met the criteria of having the 
# essential, conserved catalytic residue EDK.....N and aligned well

verified_pepm_ids <- readr::read_tsv(here("data","prok_verified_pepm_seqs.txt"), col_names = FALSE) %>% pull()

verified_mpns <- readr::read_tsv(here("data","prok_verified_mpns_seqs.tsv"), col_names = c("kaijudb_id", "class"))

# prok genomes with aln-verified pepm sequence
verified_pepm_genomes_proks <- all_coord %>% 
  mutate(kaijudb_id=ifelse(str_detect(locus_id, "^specI-"), locus_id, kaijudb_id)) %>% 
  filter(kaijudb_id %in% verified_pepm_ids) %>%
  pull(genome)

ppda <- read_tsv(here("data","prok-PPDA-reciprocalbesthits.tsv"), 
                 col_names=c("kaijudb_id", "protein_domain", "score", 
                             "evalue", "dom_score", "dom_evalue")) %>%
  extract(kaijudb_id, remove = FALSE, c("genome", "b", "c"), "(\\S+)__(\\d+)_(\\d+)") %>%
  unite(locus_id, genome, b, remove = FALSE) %>%
  mutate(locus_id=ifelse(locus_id=="NA_NA", kaijudb_id, locus_id)) %>%
  mutate(genome=ifelse(is.na(genome), str_match(kaijudb_id, "(specI-\\d+)"), genome)) %>%
  select(-b, -c)


all_coord %>% 
  filter(genome %in% verified_pepm_genomes_proks) %>%
  filter(gene_symbol=="PPDA") %>% nrow()

all_coord %>% 
  filter(genome %in% verified_pepm_genomes_proks) %>%
  filter(gene_symbol=="PPDA") %>% distinct(genome) %>% nrow()

nrow(ppda)

# filters from all coordinates only the PEPM sequences that are alignment
# verified and the genes that are within 5000 bp of the start and end
# of the verified PEPM sequence. Matches alignment verified MPNS family seqs
# to those neighborhoods

verified_pepm_neighborhoods <- all_coord %>% 
  mutate(kaijudb_id=ifelse(str_detect(locus_id, "^specI-"), locus_id, kaijudb_id)) %>% 
  left_join(., verified_mpns) %>% 
  mutate(gene_symbol=ifelse(gene_symbol=="MPNS", class, gene_symbol)) %>%
  select(-class) %>%
  mutate(gene_symbol=ifelse(!(kaijudb_id %in% verified_pepm_ids) & 
                              gene_symbol=="PEPM", "PEPM_like", gene_symbol)) %>%
  group_by(genome, scaffold) %>%
  mutate(hi=ifelse(kaijudb_id %in% verified_pepm_ids, end+5000, NA)) %>%
  mutate(lo=ifelse(kaijudb_id %in% verified_pepm_ids, start-5000, NA)) %>%
  fill(hi, lo, .direction = "up") %>%
  fill(hi, lo, .direction = "down") %>%
  filter(start < hi & start > lo) %>%
  ungroup()

# get counts of the core phosphonate biosynthetic genes from the alignment
# verified PEPM neighborhoods. Counts are per scaffold per genome. Some 
# MAGs appear to have multiple scaffolds with a verified PEPM and PPDA.
# I think this is more likely a product of metagenomic binning errors
# rather than true duplication events

verified_pepm_neighborhoods_counts <- verified_pepm_neighborhoods %>% 
  group_by(genome, scaffold) %>% #
  filter(!(is.na(gene_symbol))) %>%
  count(gene_symbol) %>% 
  spread(gene_symbol, n) %>%
  replace_na(list(PEPM=0, PEPM_like=0, PPDA=0, PPDH=0, MPNS=0, HEPDI=0, HEPDII=0)) %>%
  summarize(PEPM=sum(PEPM), PEPM_like=sum(PEPM_like), MPNS=sum(MPNS), HEPDI=sum(HEPDI), HEPDII=sum(HEPDII), PPDA=sum(PPDA), PPDH=sum(PPDH)) %>%
  ungroup() %>%
  left_join(., verified_pepm_neighborhoods %>% select(genome, origin) %>% distinct(genome, origin))

# There are 1062 alignment verified PEPM sequences
length(verified_pepm_ids)

verified_pepm_neighborhoods_counts %>% summarize(sum=sum(PEPM))

# There are 1044 (339+705) scaffolds/genomes with one alignment verified PEPM sequence

# There are 8 (2+2+2+2) genomes with multiple verified PEPM sequences (accounting
# for 18 PEPMs) on the same contig/scaffold. Most are from complete progenomes

# There are 707 (705+2) scaffolds/genomes with one alignment verified PPDA sequence

# There are 4 genomes (accounting for 8 PPDA) with multiple verified PPDA sequences on the same
# contig/scaffold. Most are from complete progenomes

# There are 341 (339+2) genomes with no PPDA

verified_pepm_neighborhoods_counts %>% count(PEPM, PPDA)

# accounting for all PEPM sequences
339+705+(2*2)+(2*2)+(2*2)+(3*2)

# All the multiple PPDA genomes correspond to a genome with multiple
# alignment verified PEPM sequences

verified_pepm_neighborhoods_counts %>% filter(PPDA >= 1 & PEPM > 1) 

# extract the 6 scaffolds with multiple verified PEPM hits and one or more PPDA sequences
multiscaffs <- verified_pepm_neighborhoods_counts %>% filter(PPDA >= 1 & PEPM > 1) %>% pull(scaffold)

# extract the 705 scaffolds/genomes with 1 PEPM and 1 PPDA
singlescaffs <- verified_pepm_neighborhoods_counts %>% filter(PPDA == 1 & PEPM == 1) %>% pull(scaffold)

# 719/1062 PEPMs and 715/852 PPDAs are co-located in neighborhoods

verified_pepm_neighborhoods_counts %>% 
  count(PEPM, PPDA) %>% 
  mutate(num_PEPM=PEPM*n, num_PPDA=PPDA*n)

verified_pepm_neighborhoods_counts %>% 
  count(PEPM, PPDA) %>% 
  mutate(num_PEPM=PEPM*n, num_PPDA=PPDA*n) %>%
  summarize(PDPM_sum=sum(num_PEPM), PPDA_sum=sum(num_PPDA))

# PEPM with PPDA
705+4+4+6

# PPDA 
705+2+4+4

verified_pepm_neighborhoods_counts %>% 
  count(PEPM, PPDA) %>% 
  mutate(num_PEPM=PEPM*n, num_PPDA=PPDA*n) %>%
  summarize(PEPM_sum=sum(num_PEPM), PPDA_sum=sum(num_PPDA))

# good these 2 show the same
verified_pepm_neighborhoods_counts %>% 
  filter(PPDA==0) %>% count(PEPM, PPDA)

verified_pepm_neighborhoods_counts %>% 
  filter(!(scaffold %in% c(multiscaffs, singlescaffs))) %>% count(PEPM, PPDA)

verified_pepm_neighborhoods %>%
  filter(!(scaffold %in% c(multiscaffs, singlescaffs)))

# Now let's look if any of the 343 PEPMs have one of the unaccounted
# for 137 PPDA sequences on a different contig/scaffold somewhere 
# else in the genome 

1062-719
852-715

x <- verified_pepm_neighborhoods_counts %>% 
  filter(!(scaffold %in% c(multiscaffs, singlescaffs))) %>%
  pull(genome)

y <- verified_pepm_neighborhoods %>% filter(gene_symbol=="PPDA") %>% pull(kaijudb_id)

z <- ppda %>% 
  filter((genome %in% x) & !(kaijudb_id %in% y)) %>% 
  mutate(gene_symbol="PPDAorphan")%>% 
  select(-score, -dom_score, -dom_evalue)

verified_pepm_orphans_counts <- verified_pepm_neighborhoods %>%
  filter(!(scaffold %in% c(multiscaffs, singlescaffs))) %>%
  bind_rows(., z) %>%
  group_by(genome) %>%
  filter(!(is.na(gene_symbol))) %>%
  count(gene_symbol) %>% 
  spread(gene_symbol, n) %>%
  replace_na(list(PEPM=0, PEPM_like=0, PPDAorphan=0, PPDH=0, MPNS=0, HEPDII=0)) %>%
  summarize(PEPM=sum(PEPM), PEPM_like=sum(PEPM_like), MPNS=sum(MPNS), 
            HEPDII=sum(HEPDII), PPDAorphan=sum(PPDAorphan), PPDH=sum(PPDH)) %>%
  ungroup() %>%
  left_join(., verified_pepm_neighborhoods %>% select(genome, origin) %>% distinct(genome, origin))

# 68 PPDAs are located on other scaffolds in genomes with an
# alignment verified PEPM sequence and without a co-located 
# PPDA sequence. 

verified_pepm_orphans_counts %>%
  count(PEPM, PPDAorphan) %>%
  mutate(num_PEPM=PEPM*n, num_PPDAorphan=PPDAorphan*n)

# Get genome that has 2 PPDAorphans, UBA3669
verified_pepm_orphans_counts %>% filter(PPDAorphan==2 & PEPM==1)

verified_pepm_orphans_counts1 <- semi_join(verified_pepm_neighborhoods_counts %>% filter(PPDA==0),
                                           verified_pepm_orphans_counts %>% filter(PPDAorphan >= 1), by="genome") %>%
  mutate(PPDAorphan=ifelse(genome=="UBA3669", 2, 1)) 

verified_pepm_orphans_counts1 %>% filter(genome=="UBA4338")
verified_pepm_neighborhoods_counts %>% filter(genome=="UBA4338")

final_counts <- left_join(verified_pepm_neighborhoods_counts,
                          verified_pepm_orphans_counts1 %>% select(genome, PPDAorphan, origin),
                          by=c("genome", "origin")) %>% 
  # to account for the annoyingness of UBA4338 already having a 
  # scaffold with 1 PEPM and 1 PPDA. It's the only genome with this...
  mutate(PPDAorphan=ifelse(scaffold=="gnl|ACE|UBA4338_6", NA, PPDAorphan)) %>%
  replace_na(list(PPDAorphan=0))

final_counts %>% count(PEPM, PPDA, PPDAorphan)

final_counts %>% count(PEPM, PPDA, PPDAorphan) %>%
  summarize(sumPEPM=sum(PEPM*n), sumPPDA=sum(PPDA*n), sumPPDAorphan=sum(PPDAorphan*n))

#PEPM --> 1062 all accounted for

#PPDA --> 715 adjacent PEPM accounted for

#PPDAorphan --> additional 68 accounted for

# Total number of PPDA's accounted for --> 783

final_counts %>% count(PEPM, PPDA, PPDAorphan, MPNS, HEPDI, HEPDII)

missingMPNS <- anti_join(verified_mpns, verified_pepm_neighborhoods) %>%
  extract(kaijudb_id, remove = FALSE, c("genome", "b", "c"), "(\\S+)__(\\d+)_(\\d+)")%>% mutate(genome=ifelse(is.na(genome), str_match(kaijudb_id, "(specI-\\d+)"), genome)) %>%
  select(-b, -c, -kaijudb_id) %>%
  group_by(genome) %>%
  count(class) %>% 
  spread(class, n)



lost_ids <- tribble(~origin, ~genome, "gbkmags", "CPC51", "gbkmags", "EAC1785", "gbkmags", "EAC1786", "gbkmags", "EAC654", "gbkmags", "MED821", "gbkmags", "NAT220", "gbkmags", "NAT226", "gbkmags", "REDSEA-S12_B4", "gbkmags", "RS374", "gbkmags", "RS457", "gbkmags", "SAT157", "gbkmags", "SAT42", "gbkmags", "SP115", "gbkmags", "SP232", "gbkmags", "TMED256", "gbkmags", "UBA1136", "gbkmags", "UBA1226", "gbkmags", "UBA1326", "gbkmags", "UBA1361", "gbkmags", "UBA1706", "gbkmags", "UBA2338", "gbkmags", "UBA2588", "gbkmags", "UBA2837", "gbkmags", "UBA2845", "gbkmags", "UBA2943", "gbkmags", "UBA2962", "gbkmags", "UBA2973", "gbkmags", "UBA2976", "gbkmags", "UBA3211", "gbkmags", "UBA3814", "gbkmags", "UBA4420", "gbkmags", "UBA4468", "gbkmags", "UBA4479", "gbkmags", "UBA4587", "gbkmags", "UBA4697", "gbkmags", "UBA4723", "gbkmags", "UBA4898", "gbkmags", "UBA4914", "gbkmags", "UBA5263", "gbkmags", "UBA5479", "gbkmags", "UBA6407", "gbkmags", "UBA6605", "gbkmags", "UBA7162", "gbkmags", "UBA7503", "gbkmags", "UBA830", "geba", "DSM2950", "marref", "MMP02469388", "marref", "MMP02469391", "marref", "MMP03009086", "marref", "MMP03699827", "marref", "MMP05282066", "marref", "MMP08161523", "pro", "AG-363-B05", "progenomes", "specI-1121114", "progenomes", "specI-1137268", "progenomes", "specI-1157637", "progenomes", "specI-1169161", "progenomes", "specI-1232451", "progenomes", "specI-1380770", "progenomes", "specI-1415559", "progenomes", "specI-1415626", "progenomes", "specI-1464079", "progenomes", "specI-146922", "progenomes", "specI-1609104", "progenomes", "specI-189753", "progenomes", "specI-397287", "progenomes", "specI-443598", "progenomes", "specI-452637", "progenomes", "specI-66897", "progenomes", "specI-67287", "progenomes", "specI-68570", "progenomes", "specI-742733", "progenomes", "specI-796334", "progenomes", "specI-882085", "simonshets", "AG-410-A08")

crappers <- c("specI-1463830", "specI-1463894", "specI-1463906", "specI-1463907", "specI-1463914", "specI-1463923", "specI-1463931", "specI-281186", "specI-67352", "specI-80858")

a <- verified_mpns %>%
  extract(kaijudb_id, remove = FALSE, c("genome", "b", "c"), "(\\S+)__(\\d+)_(\\d+)")%>% mutate(genome=ifelse(is.na(genome), str_match(kaijudb_id, "(specI-\\d+)"), genome)) %>% 
  rename(gene_symbol=class) %>%
  select(-kaijudb_id, -b, -c)

b <- all_coord %>% 
  mutate(kaijudb_id=ifelse(str_detect(locus_id, "^specI-"), locus_id, kaijudb_id)) %>% 
  mutate(gene_symbol=ifelse(!(kaijudb_id %in% verified_pepm_ids) & 
                              gene_symbol=="PEPM", "PEPM_like", gene_symbol)) %>%
  #mutate(gene_symbol=ifelse(gene_symbol=="PPDA", "PPDA_coloco", gene_symbol)) %>%
  filter(gene_symbol %in% c("PEPM", "PEPM_like", "PPDH")) %>% #"PPDA_coloco", 
  select(genome, gene_symbol, origin)

c <- ppda %>% mutate(gene_symbol="PPDA") %>% select(genome, gene_symbol)

d <- bind_rows(a, b, c) %>% 
  group_by(genome) %>% #
  count(gene_symbol) %>% 
  spread(gene_symbol, n) %>%
  replace_na(list(PEPM=0, PEPM_like=0, PPDA=0, PPDH=0, MPNS=0, HEPDI=0, HEPDII=0)) %>% #PPDA_coloco=0, 
  filter(PEPM>=1 | PPDA>=1 | MPNS >=1 | HEPDI >=1 | HEPDII >=1) %>%
  ungroup() %>% arrange(desc(PEPM_like)) %>% 
  left_join(., final_counts %>% select(genome, origin)) %>%
  left_join(., lost_ids, by="genome") %>%
  mutate(origin=ifelse(is.na(origin.x), origin.y, origin.x)) %>%
  select(-origin.x, -origin.y)

bind_rows(a, b, c) %>% count(gene_symbol)
d %>% filter(genome %in% crappers)

d_binary <- bind_rows(a, b, c) %>% 
  group_by(genome) %>% #
  count(gene_symbol) %>% 
  spread(gene_symbol, n) %>%
  replace_na(list(PEPM=0, PEPM_like=0, PPDA=0, PPDH=0, MPNS=0, HEPDI=0, HEPDII=0)) %>%
  filter(PEPM>=1 | PPDA>=1 | MPNS >=1 | HEPDI >=1 | HEPDII >=1) %>%
  mutate_all(funs(ifelse(. > 0, 1, 0))) %>% 
  left_join(., final_counts %>% select(genome, origin))  %>%
  left_join(., lost_ids, by="genome") %>%
  mutate(origin=ifelse(is.na(origin.x), origin.y, origin.x)) %>%
  select(-origin.x, -origin.y)

d_reduced <- bind_rows(a, b, c) %>% 
  mutate(gene_symbol=ifelse(gene_symbol %in% c("MPNS", "HEPDI", "HEPDII"), "MPNS", gene_symbol)) %>%
  group_by(genome) %>% #
  count(gene_symbol) %>% 
  spread(gene_symbol, n) %>%
  replace_na(list(PEPM=0, PEPM_like=0, PPDA=0, PPDH=0, MPNS=0)) %>% 
  filter(PEPM>=1 | PPDA>=1 | MPNS >=1 ) %>% select(-PEPM_like) %>% 
  left_join(., final_counts %>% select(genome, origin))  %>%
  left_join(., lost_ids, by="genome") %>%
  mutate(origin=ifelse(is.na(origin.x), origin.y, origin.x)) %>%
  select(-origin.x, -origin.y)

d_reduced_binary <- bind_rows(a, b, c) %>% 
  mutate(gene_symbol=ifelse(gene_symbol %in% c("MPNS", "HEPDI", "HEPDII"), "MPNS", gene_symbol)) %>%
  group_by(genome) %>% #
  count(gene_symbol) %>% 
  spread(gene_symbol, n) %>%
  replace_na(list(PEPM=0, PEPM_like=0, PPDA=0, PPDH=0, MPNS=0)) %>% 
  filter(PEPM>=1 | PPDA>=1 | MPNS >=1 ) %>% select(-PEPM_like) %>%
  mutate_all(funs(ifelse(. > 0, 1, 0))) %>% 
  left_join(., final_counts %>% select(genome, origin))  %>%
  left_join(., lost_ids, by="genome") %>%
  mutate(origin=ifelse(is.na(origin.x), origin.y, origin.x)) %>%
  select(-origin.x, -origin.y)

write_tsv(d, here("tables", "genome_counts.tsv"))
write_tsv(d_binary, here("tables", "genome_counts.binary.tsv"))
write_tsv(d_reduced, here("tables", "genome_counts.reduced.tsv"))
write_tsv(d_reduced_binary, here("tables", "genome_counts.reduced.binary.tsv"))

##################### MAKING GENE DIAGRAMS ##########################
jointgenomes <- tibble::tribble(
  ~genome, ~group,         ~clade, ~order,
  "AG-325-E08",  "pel",        "Ia_3a",      1,
  "AG-422-C23",  "pel",        "Ia_3a",      2,
  "AG-447-N18",  "pel",        "Ia_3a",      3,
  "AG-345-O09",  "pel",        "Ia_3a",      4,
  "HTCC7217",  "pel",        "Ia_3b",      5,
  "AG-319-F18",  "pel",        "Ia_3b",      6,
  "RS40",  "pel",         "Ib_1",      7,
  "AG-410-D23",  "pel",         "Ib_1",      8,
  "AG-365-D07",  "pel",         "Ib_2",      9,
  "AG-426-G09",  "pel",         "Ib_2",     10,
  "TMED64",  "pel",          "IIb",     11,
  "MED605",  "pel",          "IIb",     12,
  "AG-345-F21",  "pel",        "IIa_1",     13,
  "AG-390-C22",  "pel",        "IIa_1",     14,
  "TMED216",  "pel",        "Iia_1",     15,
  "MED1153",  "pel",        "IIa_1",     16,
  "AG-430-I06",  "pel",        "IIa_2",     17,
  "AG-430-L02",  "pel",        "IIa_2",     18,
  "NAT182",  "pel",           "5a",     19,
  "SP38",  "pel",      "unknown",     20,
  "TMED65",  "pel",      "unknown",     21,
  "SB",  "pro",         "HLII",     22,
  "AAA795-I06",  "pro",         "HLII",     23,
  "AG-347-I06",  "pro",         "HLII",     24,
  "AG-347-I22",  "pro",         "HLII",     25,
  "AG-418-O03",  "pro",         "HLII",     26,
  "AG-418-P13",  "pro",         "HLII",     27,
  "AG-436-D21",  "pro",         "HLII",     28,
  "AG-402-K10",  "pro",         "HLVI",     29,
  "AG-469-F22",  "pro",          "HLI",     30,
  "AG-670-L08",  "pro",          "HLI",     31,
  "AG-670-O11",  "pro",          "HLI",     32,
  "AG-436-C13",  "pro",         "LLII",     33,
  "AG-363-C20",  "pro",       "LLVIIa",     34,
  "AG-409-M05",  "pro",       "LLVIIe",     35,
  "AG-676-A21",  "syn", "5.1BV.VI.VII",     36,
  "AG-670-B23",  "syn",       "5.1AII",     37
)


pro <- verified_pepm_neighborhoods %>% 
  filter(origin %in% c("prochlorococcus", "synechococcus")) %>% 
  mutate(gene_symbol=ifelse(is.na(gene_symbol), "Other", gene_symbol)) %>%
  mutate(start1=ifelse(strand==-1, end, start)) %>% 
  mutate(end1=ifelse(strand==-1, start, end)) %>%
  left_join(., jointgenomes) %>%
  arrange(order) %>%
  mutate(genome = factor(genome))

## make a dummy dataframe to force centering all coordinates on PEPM gene
dummies <- make_alignment_dummies(
  pro,
  aes(xmin = start1, xmax = end1, y = genome, id = gene_symbol),
  on = "PEPM"
)

ggplot(pro) +
  geom_gene_arrow(aes(y=genome, xmin = start1, xmax = end1, fill = gene_symbol),
                  arrowhead_height = unit(2.25, "mm"), 
                  arrowhead_width = unit(1.25, "mm"),
                  arrow_body_height = unit(2.25, "mm")) +
  geom_blank(data = dummies, aes(y=genome, xmin = start1, xmax = end1)) +
  facet_wrap(genome~.,  scales = "free" ,ncol = 1) +
  scale_fill_manual(values=c("#66c2a5", "#fc8d62", "#FFFFFF", "#8da0cb", "#e78ac3", "#a6d854")) +
  ylab("Genome") +
  xlab("Contig coordinate [bp]") +
  theme_genes()

ggsave(here("figs", "progenes.svg"), width=5, height=6.46, units="in")


pel <- verified_pepm_neighborhoods %>% filter(origin=="pelagibacterales") %>% 
  mutate(genome=ifelse(scaffold=="MED605_11", "MED605_1", genome)) %>%
  mutate(gene_symbol=ifelse(is.na(gene_symbol), "Other", gene_symbol)) %>%
  mutate(start1=ifelse(strand==-1, end, start)) %>% 
  mutate(end1=ifelse(strand==-1, start, end))

## make a dummy dataframe to force centering all coordinates on PEPM gene
dummies <- make_alignment_dummies(
  pel,
  aes(xmin = start1, xmax = end1, y = scaffold, id = gene_symbol),
  on = "PEPM"
)

ggplot(pel) +
  geom_gene_arrow(aes(y=scaffold, xmin = start1, xmax = end1, fill = gene_symbol),
                  arrowhead_height = unit(2.25, "mm"), 
                  arrowhead_width = unit(1.25, "mm"),
                  arrow_body_height = unit(2.25, "mm")) +
  geom_blank(data = dummies, aes(y=scaffold, xmin = start1, xmax = end1)) +
  facet_wrap(scaffold~. , scales = "free", ncol = 1) +
  scale_fill_manual(values=c("#66c2a5", "#fc8d62", "#FFFFFF", "#8da0cb", "#e78ac3", "#a6d854")) +
  ylab("Genome") +
  xlab("Contig coordinate [bp]") +
  theme_genes()

ggsave(here("figs", "pelgenes.svg"), width=5, height=8.46, units="in")