library(tidyverse)
library(here)

mmdb <- read_tsv(here("tables", "MARMICRODB_compiled.tsv"))

phn <- read_tsv(here("tables", "genome_counts.binary.tsv")) %>%
  select(-origin)

phn_a <- left_join(phn, mmdb, by = "genome")

write_tsv(phn_a %>% select(lineage_assignment), here("tables", "import2krona.txt"))

grab <- read_tsv("/data/veriried-pepm-selfsim.m8", col_names = FALSE)

grab1 <- grab %>%
  filter(X3 > 0.98) %>%
  filter(X1 != X2) %>%
  separate(X1, c("genome1", "gene1"), "_+", remove = FALSE) %>%
  separate(X2, c("genome2", "gene2"), "_+", remove = FALSE)

grab2 <- left_join(grab1, mmdb %>% rename(genome1 = genome)) %>%
  select(X1, genome1, X2, genome2, X3, X4, full_name, lineage_assignment) %>%
  rename(full_name1 = full_name, lineage_assignemnt1 = lineage_assignment)

grab3 <- left_join(grab2, mmdb %>% rename(genome2 = genome)) %>%
  select(X1, genome1, full_name1, lineage_assignemnt1, X2, genome2, X3, X4, full_name, lineage_assignment) %>%
  rename(full_name2 = full_name, lineage_assignemnt2 = lineage_assignment) %>%
  select(genome1, genome2, lineage_assignemnt1, lineage_assignemnt2, X3, X4)
