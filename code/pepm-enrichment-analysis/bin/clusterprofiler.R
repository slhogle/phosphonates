library(httr)
library(tidyverse)
library(clusterProfiler)
library(KEGGREST)
library(DOSE)
library(here)

#### only run once... slow
#### actually not really worth it
#### but I learned how to use textConnection... 

path2ko <- read.delim(textConnection(content(GET("http://rest.kegg.jp/link/ko/path"), "text")), sep = "\t", header=FALSE) %>% rename(path=V1, ko=V2)

pathlist <- read.delim(textConnection(content(GET("http://rest.kegg.jp/list/path"), "text")), sep = "\t", header=FALSE) %>% rename(path=V1, desc=V2)

module2ko <- read.delim(textConnection(content(GET("http://rest.kegg.jp/link/ko/module"), "text")), sep = "\t", header=FALSE) %>% rename(module=V1, ko=V2)

modulelist <- read.delim(textConnection(content(GET("http://rest.kegg.jp/list/module"), "text")), sep = "\t", header=FALSE) %>% rename(module=V1, desc=V2)

################ UNIVERSAL ENRICHMENT ANALYSIS
### see https://en.wikipedia.org/wiki/Hypergeometric_distribution#Hypergeometric_test
### https://stats.stackexchange.com/questions/300309/how-to-use-hyper-geometric-test
### http://guangchuangyu.github.io/2015/05/use-clusterprofiler-as-an-universal-enrichment-analysis-tool/
### https://bioconductor.org/packages/devel/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html#universal-enrichment-analysis

path2ko <- readr::read_tsv(here("data", "kegg", "path2ko.tsv"), col_names = c("PATH", "KO"))
pathdesc <- readr::read_tsv(here("data", "kegg", "path_desc.tsv"), col_names = c("PATH", "DESC"))

marmicrodb <- readr::read_tsv(here("data", "MARMICRODB_catalog.tsv"), col_names = TRUE)

#### for reading in coordinates
READcoords <- function(pattern = pattern, path = path){
  map(list.files(path = path, pattern = pattern, full.names = TRUE), read_csv)
}

safely_READcoords <- safely(READcoords)

coords <- flatten_dfr(safely_READcoords(".csv", "./data/coords")) %>%
  extract(locus_id, remove = FALSE, c("genome", "b"), "(\\S+)_(\\d+)") %>% select(-b)


#### for reading in eggnogg mapper results
READeggnogg <- function(pattern = pattern, path = path){
  map(list.files(path = path, pattern = pattern, full.names = TRUE), read_tsv)
}

safely_READeggnogg <- safely(READeggnogg)

noggs <- flatten_dfr(safely_READeggnogg(".annotations", "./data/eggnoggmapped")) %>%
  extract(query_name, remove = FALSE, c("genome", "b"), "(\\S+)_(\\d+)") %>% select(-b) %>%
  select(genome, ID=query_name, KO=KEGG_KOs, desc=`eggNOG annot`)

##### Additional formatting
noggs1 <- noggs %>%
  separate_rows(KO) %>%
  mutate(dum="ko:") %>%
  unite(KO, dum, KO, sep='', remove=TRUE) %>% 
  mutate(KO=ifelse(KO=="ko:NA", NA, KO))

##### join with kegg descriptions and pathways
GENOME_KO_PATH <- left_join(noggs1, path2ko) %>% 
  filter(str_detect(PATH, "map") | is.na(PATH)) %>% 
  left_join(., pathdesc)

##### calculating enrichments for individual genomes
genomes <- coords %>% select(genome) %>% distinct() %>% pull()
dat <- data.frame()      

for (g in genomes){
  PATH <- GENOME_KO_PATH %>% filter(genome==g)
  GENES <- coords %>% 
    filter(genome==g) %>% 
    pull(locus_id)
  
  path2name = PATH[, c("PATH", "DESC")] %>% drop_na(PATH)
  path2gene = PATH[, c("PATH", "ID")]
  
  x = enricher(GENES, TERM2GENE=path2gene, 
             TERM2NAME=path2name, 
             pAdjustMethod = "BH") #, pvalueCutoff = 1, qvalueCutoff = 1

  dat <- bind_rows(dat, as.data.frame(x))
}

t <- dat %>% separate_rows(geneID, sep = "/") %>% 
    left_join(., noggs %>% rename(geneID=ID)) %>% 
    left_join(., coords %>% rename(geneID=locus_id)) %>%
    left_join(., marmicrodb)

### keep only the genes that are not PEPM, PPDA, PPDH, or MPNS
t1 <- t %>% filter(is.na(gene_symbol)) 

#### write table of KOs matching amino sugar metabolism
#### to upload to https://www.genome.jp/kegg/tool/map_pathway1.html
write_tsv(t1 %>% 
  filter(Description=="Amino sugar and nucleotide sugar metabolism") %>% 
  dplyr::select(KO) %>%
  separate_rows(KO), here("output", "map00520_KO.txt"), col_names=FALSE)

#### write a full table of the results with qvalues, enrichment score, etc...
write_tsv(t %>% arrange(geneID, genome) %>% select(-ft_type, -kaijudb_id, -binlocation, -marine, -source, -MARMICRODBtaxid),
          here("output", "enrichment_scores.tsv"))   
