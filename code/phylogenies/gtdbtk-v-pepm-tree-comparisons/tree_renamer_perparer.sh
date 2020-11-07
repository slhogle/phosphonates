#!/usr/bin/env bash

# the phyx tools for phylogenetic tree manipulation
# Brown, J. W., J. F. Walker, and S. A. Smith; Phyx: phylogenetic tools for unix. Bioinformatics 2017; 33 (12): 1886-1888. doi: 10.1093/bioinformatics/btx063

### instructions for getting tools. see github for dependencies
#wget https://github.com/FePhyFoFum/phyx/archive/v0.999.tar.gz
#tar xvf phyx-0.999.tar.gz
#cd phyx-0.999
#make
#make install

### for bacteria
pxlstr -i -t bacteria_PEPM.curated.gappyout.tre > bac_pepm_nodes.old.txt
pxlstr -i -t bacteria_PEPM.curated.gappyout.tre | awk -F'_' '{print $1}' > bac_pepm_nodes.new.txt
pxrlt -t bacteria_PEPM.curated.gappyout.tre -c bac_pepm_nodes.old.txt -n bac_pepm_nodes.new.txt > bacteria_PEPM.curated.gappyout.genomes.tre
comm -12 <(pxlstr -i -t gtdbtk.bac120.user_msa.trimgappyout.tre | sort) <(pxlstr -i -t bacteria_PEPM.curated.gappyout.genomes.tre | sort) > bac-shared-ids.txt
pxrmt -c -t bacteria_PEPM.curated.gappyout.genomes.tre -f bac-shared-ids.txt > bacteria_PEPM.curated.gappyout.genomes.shared.tre
pxrmt -c -t gtdbtk.bac120.user_msa.trimgappyout.tre -f bac-shared-ids.txt > gtdbtk.bac120.user_msa.trimgappyout.shared.tre 
rm bac_pepm_nodes.old.txt bac_pepm_nodes.new.txt bacteria_PEPM.curated.gappyout.genomes.tre bac-shared-ids.txt

### for archaea
pxlstr -i -t archaea_PEPM.curated.gappyout.tre > arc_pepm_nodes.old.txt
pxlstr -i -t archaea_PEPM.curated.gappyout.tre | awk -F'_' '{print $1}' > arc_pepm_nodes.new.txt
pxrlt -t archaea_PEPM.curated.gappyout.tre -c arc_pepm_nodes.old.txt -n arc_pepm_nods.new.txt > archaea_PEPM.curated.gappyout.genomes.tre
comm -12 <(pxlstr -i -t gtdbtk.ar122.user_msa.trimgappyout.tre | sort) <(pxlstr -i -t archaea_PEPM.curated.gappyout.genomes.tre | sort) > arc-shared-ids.txt
pxrmt -c -t archaea_PEPM.curated.gappyout.genomes.tre -f arc-shared-ids.txt > archaea_PEPM.curated.gappyout.genomes.shared.tre
pxrmt -c -t gtdbtk.ar122.user_msa.trimgappyout.tre -f arc-shared-ids.txt > gtdbtk.ar122.user_msa.trimgappyout.shared.tre
rm arc_pepm_nodes.old.txt arc_pepm_nodes.new.txt archaea_PEPM.curated.gappyout.genomes.tre arc-shared-ids.txt