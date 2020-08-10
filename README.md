# Supporting data and code for "Phosphonate production in marine microbes: exploring new sources and potential function"

Data and code here is provided “AS IS” WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED OR IMPLIED under the GNU General Public License v3.0. Feel free to use or remix as you see fit.

The `data` directory contains results and summary tables needed to support the conclusions of the paper. It is basically formatted output from the `code` directory and should be clearly organized and user friendly. The `code` directory has many different bash and R scripts in varying degrees of "user-friendliness" and is presented here mainly in the interest of transparency and not as a software tool for production use. The directories are organized by the type of analysis performed in the paper and *should be* modifiable to suite your computing environment and needs if you have the time and programming experience.

## Handy shortcuts

### Phosphonate biosynthesis
- [Summary of PepM genes](data/phosphonate_biosynthesis/MARMICRODB/pepm_genes.md) identified in [MARMICRODB](https://zenodo.org/record/3520509)
- [Summary of MpnS genes](data/phosphonate_biosynthesis/MARMICRODB/mpns_genes.md) identified in [MARMICRODB](https://zenodo.org/record/3520509)
- [HMM models](data/phosphonate_biosynthesis/HMM_models) for identifying phosphonate biosynthesis genes

### Phosphonate catabolism
- [Gene families](data/phosphonate_catabolism/HMM_models/phosphonate_utilization_families.md) used to identify 4 discrete phosphonate catabolism pathways in genomes
- [HMM models](data/phosphonate_catabolism/HMM_models) for identifying phosphonate catabolism genes

### Core genes families
- [Summary of core genes](data/core_gene_families/core_gene_families.tsv) used for metagenome normalization
- [HMM models](data/core_gene_families/HMM_models) used for identifying genes used in metagenome normalization
- [HMM alignments](data/core_gene_families/HMM_alignments) used to produce core gene HMMs

```
.
├── code
└── data
    ├── core_gene_families
    │   ├── HMM_alignments
    │   └── HMM_models
    ├── phosphonate_biosynthesis
    │   ├── GORG_tropics
    │   ├── HMM_models
    │   └── MARMICRODB
    └── phosphonate_catabolism
        ├── GORG_tropics
        ├── HMM_models
        └── MARMICRODB
```