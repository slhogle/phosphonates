library(tidyverse)
library(here)
library(ape)
library(treeman)
library(dendextend)
library(DECIPHER)

#### load trees
genome_bac <- read.tree(file=here("data", "gtdbtk.bac120.user_msa.trimgappyout.shared.tre"))
pepm_bac <- read.tree(file=here("data", "bacteria_PEPM.curated.gappyout.genomes.shared.tre"))

genome_arc <- read.tree(file=here("data", "gtdbtk.ar122.user_msa.trimgappyout.shared.tre"))
pepm_arc <- read.tree(file=here("data", "archaea_PEPM.curated.gappyout.genomes.shared.tre"))

random_bac <- rtree(length(genome_bac$tip.label), tip.label = genome_bac$tip.label, rooted=FALSE)
random_arc <- rtree(length(genome_arc$tip.label), tip.label = genome_arc$tip.label, rooted=FALSE)


ete3 compare -t randomtree.tre -r pruned_species.tre --unrooted


## calculate distances between trees for bacteria
dist.topo(genome_bac, genome_bac, method = "PH85")
dist.topo(genome_bac, pepm_bac, method = "PH85")
dist.topo(genome_bac, random_bac, method = "PH85")

dist.topo(genome_bac, genome_bac, method = "score")
dist.topo(genome_bac, pepm_bac, method = "score")
dist.topo(genome_bac, random_bac, method = "score")

association <- cbind(genome_bac$tip.label, pepm_bac$tip.label)
cophyloplot(genome_bac, pepm_bac, assoc = association,
            length.line = 5, space = 68, gap = 3, show.tip.label = FALSE)

a <- ReadDendrogram(here("data", "gtdbtk.bac120.user_msa.trimgappyout.shared.tre"))
b <- ReadDendrogram(here("data", "bacteria_PEPM.curated.gappyout.genomes.shared.tre"))

#### yuck...
tanglegram(a, b, 
           common_subtrees_color_lines=FALSE, 
           highlight_distinct_edges=FALSE, 
           highlight_branches_lwd=FALSE, 
           highlight_branches_col=FALSE, 
           remove_nodePar=TRUE,
           faster=TRUE)

## calculate distances between trees for archaea
dist.topo(genome_arc, genome_arc, method = "PH85")
dist.topo(genome_arc, pepm_arc, method = "PH85")
dist.topo(genome_arc, random_arc, method = "PH85")

dist.topo(genome_arc, genome_arc, method = "score")
dist.topo(genome_arc, pepm_arc, method = "score")
dist.topo(genome_arc, random_arc, method = "score")

a <- ReadDendrogram(here("data", "gtdbtk.ar122.user_msa.trimgappyout.shared.tre"))
b <- ReadDendrogram(here("data", "archaea_PEPM.curated.gappyout.genomes.shared.tre"))

tanglegram(a, b, 
           common_subtrees_color_lines=TRUE, 
           highlight_distinct_edges=TRUE, 
           highlight_branches_lwd=FALSE, 
           highlight_branches_col=FALSE, 
           remove_nodePar=FALSE,
           faster=TRUE)