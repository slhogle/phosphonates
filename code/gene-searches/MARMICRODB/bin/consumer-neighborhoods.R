library(tidyverse)

### load phn-util coordinates
all_coord <- read_csv("prok_coords.csv") %>%
  mutate(strand=ifelse(strand=="-", -1, 
                       ifelse(strand=="+", 1, strand))) %>%
  mutate(kaijudb_id=ifelse(str_detect(locus_id, "^specI-"), locus_id, kaijudb_id))

### CPLyase neighborhoods
CPgenomes <- all_coord %>% filter(gene_symbol=="phnK") %>% distinct(genome) %>% pull()

verified_CP_neighborhoods <- all_coord %>% 
  filter(genome %in% CPgenomes) %>%
  group_by(genome, scaffold) %>%
  mutate(hi=ifelse(gene_symbol=="phnK", end+8000, NA)) %>%
  mutate(lo=ifelse(gene_symbol=="phnK", start-8000, NA)) %>%
  fill(hi, lo, .direction = "up") %>%
  fill(hi, lo, .direction = "down") %>%
  filter(start < hi & start > lo) %>%
  ungroup()

CPbinary <- verified_CP_neighborhoods %>% 
  group_by(genome) %>% #
  count(gene_symbol) %>% 
  spread(gene_symbol, n, fill = 0) %>%
  select(genome, phnC, phnD, phnE, phnF, phnG, phnH, phnI, phnJ, phnK, phnL,
         phnM, phnN, phnO, phnP, phn01, phn02) %>% 
  mutate_all(funs(ifelse(. > 0, 1, 0)))

write_tsv(CPbinary, "CPlyase_counts.binary.tsv")

### Phosphite neighborhoods
verified_phosphite_neighborhoods <- all_coord %>% 
  group_by(genome, scaffold) %>%
  mutate(hi=ifelse(gene_symbol=="phnz" | gene_symbol=="ptxD", end+5000, NA)) %>%
  mutate(lo=ifelse(gene_symbol=="phnz" | gene_symbol=="ptxD", start-5000, NA)) %>%
  fill(hi, lo, .direction = "up") %>%
  fill(hi, lo, .direction = "down") %>%
  filter(start < hi & start > lo) %>%
  ungroup()

phosphitebinary <- verified_phosphite_neighborhoods %>%
  group_by(genome) %>% #
  count(gene_symbol) %>% 
  spread(gene_symbol, n, fill = 0) %>%
  select(genome, phnC, `phnD-like`, phnE, ptxD, phyH, phnZ) %>%
  mutate_all(funs(ifelse(. > 0, 1, 0))) %>%
  mutate(mysum=(phnC + `phnD-like` + phnE + ptxD + phyH + phnZ)) %>%
  filter(mysum > 3) %>% select(-mysum)

write_tsv(phosphitebinary, "phosphite_counts.binary.tsv")


### Weird 2AEp GenPropt0713 and GenPropt0238 thing

verified_07130238_neighborhoods <- all_coord %>% 
  #filter(genome %in% CPgenomes) %>%
  group_by(genome, scaffold) %>%
  mutate(hi=ifelse(gene_symbol=="phnW" | gene_symbol=="ptxA" | gene_symbol=="phnX", end+5000, NA)) %>%
  mutate(lo=ifelse(gene_symbol=="phnW" | gene_symbol=="ptxA" | gene_symbol=="phnX", start-5000, NA)) %>%
  fill(hi, lo, .direction = "up") %>%
  fill(hi, lo, .direction = "down") %>%
  filter(start < hi & start > lo) %>%
  ungroup()

binary07130238 <- verified_0713_neighborhoods %>%
  group_by(genome) %>% #
  count(gene_symbol) %>% 
  spread(gene_symbol, n, fill = 0) %>%
  select(genome, phnX, phnW, phnA, phnY, lysR, phnS2, phnT2, phnU2, phnS1=phnS, phnT1=phnT, phnU1=phnU) %>%
  mutate(phnS=phnS1+phnS2, phnT=phnT1+phnT2, phnU=phnU1+phnU2) %>% select(-phnT1, -phnT2, -phnS1, -phnS2, -phnU1, -phnU2) %>%
  mutate_all(funs(ifelse(. > 0, 1, 0))) %>%
  mutate(mysum1=(phnW + phnA + phnY)) %>%
  mutate(mysum2=(phnW + phnX + phnS + phnT + phnU)) %>%
  filter(mysum1 > 2 | mysum2 >= 3) %>% select(-mysum1, -mysum2) 

write_tsv(binary07130238, "2AEP_counts.binary.tsv")

### HpnWXZ-GenProp0736 

verified_0736_neighborhoods <- all_coord %>% 
  group_by(genome, scaffold) %>%
  mutate(hi=ifelse(gene_symbol=="phnZ", end+5000, NA)) %>%
  mutate(lo=ifelse(gene_symbol=="phnZ", start-5000, NA)) %>%
  fill(hi, lo, .direction = "up") %>%
  fill(hi, lo, .direction = "down") %>%
  filter(start < hi & start > lo) %>%
  ungroup()

binary0736 <- verified_0736_neighborhoods %>%
  group_by(genome) %>% #
  count(gene_symbol) %>% 
  spread(gene_symbol, n, fill = 0) %>%
  select(genome, gntR, phnZ, hpnW, hpnZ) %>%
  mutate_all(funs(ifelse(. > 0, 1, 0))) %>%
  mutate(mysum=sum(gntR, phnZ, hpnW, hpnZ)) %>%
  filter(mysum > 2) %>% select(-mysum) 

write_tsv(binary0736, "Hpn_counts.binary.tsv")
