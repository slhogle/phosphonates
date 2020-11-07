library(tidyverse)
library(here)

pepmsummary <- tibble::tribble(
                 ~group, ~num_pepm,   ~sum_pepm_len,  ~min_pepm_len,    ~avg_pepm_len,  ~max_pepm_len,
  "pelagibacterales",  22, 28611,  834, 1300.5, 1674,
   "prochlorococcus",  14, 18039, 1014, 1288.5, 1338,
     "synechococcus",   2,  2898, 1278,   1449, 1620
  )

pel_ex <- read_tsv(here("data", "pelagibacterales-excluded-genomes.tsv"), 
              col_names=c("genome", "type", "fails"),
              col_types="ccc")

pro_ex <- read_tsv(here("data", "prochlorococcus-excluded-genomes.tsv"), 
                   col_names=c("genome", "type", "fails"),
                   col_types="ccc")

pel <- read_tsv(here("data", "pelagibacterales-qa-2.tsv"), 
                   col_names=TRUE, 
                col_types="ccnnnnnnnnnnnnnnnnnnnnnnnnnnn") %>%
  select(genome=`Bin Id`, comp=Completeness, size=`Genome size (bp)`) %>%
  mutate(group="pelagibacterales")

pro <- read_tsv(here("data", "prochlorococcus-qa-2.tsv"), 
                col_names=TRUE, 
                col_types="ccnnnnnnnnnnnnnnnnnnnnnnnnnnn") %>%
  select(genome=`Bin Id`, comp=Completeness, size=`Genome size (bp)`) %>%
  mutate(group="prochlorococcus")

filtered <- bind_rows(pel_ex, pro_ex) %>% pull(genome)

result <- bind_rows(pel, pro)%>% 
  filter(!(genome %in% filtered)) %>%
  mutate(estsize=size/(comp/100)) %>%
  group_by(group) %>%
  summarize(sum_size=sum(size),
            sum_estsize=sum(estsize)) %>%
  left_join(., pepmsummary) %>%
  mutate(sum_pepm_estlen=sum_pepm_len/sum_size*sum_estsize) %>%
  mutate(est_num_pepm=sum_pepm_estlen/avg_pepm_len) %>%
  left_join(., bind_rows(pel, pro) %>% 
              filter(!(genome %in% filtered)) %>% 
              count(group, group) %>% 
              rename(num_genomes=n)) %>%
  mutate(revised_prevalence=est_num_pepm/num_genomes*100,
         prevalence=num_pepm/num_genomes*100)

##### RESULTS #####
# pelagibacter prevalence = 11.7; revised prevalence = 14.5
# prochlorococcus prevalence = 2.3; revised prevalence = 2.8

