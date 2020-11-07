library(tidyverse)
library(ggtree)
library(ape)
library(patchwork)
library(here)
library(phytools)
library(thacklr)

# analyze -------------------------------------------------------------------
islands_0 <- read_tsv(here("data", "pro-623_islands-novt.tsv"))
gene_coords_0 <- read_tsv(here("data", "pro-623_gene-coords.tsv"))
meta_0 <- read_tsv(here("data", "pro-623_meta-data.tsv"))
tree_0 <- read.tree(here("data", "pro-623_bac120-raxml.nwk"))
clade_colors <- read_rds(here("data", "pro-clade-colors.rds"))

slh_genes_0 <- read_tsv(here("data", "island_searches.tsv")) %>%
  rename(genome_id = genome) %>%
  mutate_at(vars(genome_id, gene_id), str_replace_all, "-", "_")

slh_genes_1 <-  slh_genes_0 %>% left_join(gene_coords_0) %>%
  filter(contig_id %in% islands_0$contig_id) %>%
  filter(gene_type != "siderophore")
slh_genes_1

islands_1 <- filter(islands_0, contig_id %in% slh_genes_1$contig_id)

SB <- tibble::tribble(
        ~genome_id, ~contig_id,  ~start,    ~end,
              "SB",     "SB_c",       1,  420066,
              "SB",     "SB_c",  420067, 1657595,
              "SB",     "SB_c", 1657596, 1668514
        )

ggplot(islands_1) +
  geom_segment(aes(x=start, xend=end, y=contig_id, yend=contig_id), color="grey50", size=1) +
  geom_point(aes(x=(start+end)/2, y=contig_id, color=gene_type), slh_genes_1) +
  theme_minimal()

# tree ----------------------------------------------------------------------
# order by tree
# root tree
LLIVs <- filter(meta_0, clade == "LLIV") %>% pull(genome_id)
tree_1 <- reroot_by_outgroup(tree_0, LLIVs)

# subset tree
genomes_to_drop <- setdiff(tree_1$tip.label, slh_genes_1$genome_id)
tree_1 <- drop.tip(tree_1,  genomes_to_drop)
tree_1 <- rotate(tree_1, 13)

trgg_1 <- fortify(tree_1) %>% left_join(meta_0, by=c(label="genome_id"))

gg_trgg <- ggtree(trgg_1) +
  geom_tiplab(aes(color=clade), align=T, offset=.02) +
  scale_x_continuous(expand=c(0, .5)) +
  scale_color_manual(values = clade_colors)

# islands and tree
tree_y <- trgg_1 %>% select(genome_id=label, y)

islands_2 <- islands_1 %>% left_join(tree_y)
slh_genes_2 <- slh_genes_1 %>% left_join(tree_y)

SB <- tibble::tribble(
        ~genome_id, ~contig_id,  ~start,    ~end, ~y,
              "SB",     "SB_1",       1, 1237529,  9,
              "SB",     "SB_2", 1237530, 1657596,  9,
              "SB",     "SB_3", 1657597, 1668516,  9
        )

gg_isl <- ggplot(islands_2) +
  geom_segment(aes(x=start, xend=end, y=y, yend=y, color=contig_id), size=0.5, SB) +
  geom_segment(aes(x=start, xend=end, y=y, yend=y), color="grey50", size=2) +
  geom_point(aes(x=(start+end)/2, y=y, color=gene_type), slh_genes_2, size=3) + 
  labs(x="Genomic Position", y="") +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y.left = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = "none",
        axis.line = element_line(colour = "black"))

gg_trgg + gg_isl + plot_layout(widths = c(2,8))

ggsave(here("figs", "pepm_islands.svg"))

       