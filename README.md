# Colonization and diversification history of Madagascan palms with new phylogenomic evidence from the genus *Orania* (Arecaceae)


This repository contains all scripts and supplementary material generated during the preparation of my Master's thesis.

The file XXX contains the submitted .pdf version of this Master's thesis.

The supplementary materials to the Part B and Part C of the thesis can be found under [/supplementary_material](./supplementary_material).

### PART B

The folders [/1_phylo_reconstruction](./1_phylo_reconstruction), [/2_phylo_dating](./2_phylo_dating), and [/3_biogeography](./3_biogeography) correspond to the "Part B" of the Master's thesis, which aimed to reconstruct a dated phylogeny of the POS clade and apply  ancestral range estimations to this phylogeny. Each folder contains a detailed description of the methods applied and how to run the scripts.

* [/1_phylo_reconstruction](./1_phylo_reconstruction): scripts related to the reconstruction of gene and species trees. 

* [/2_phylo_dating](./2_phylo_dating): scripts related to dating of the previously estimated species tree

* [/3_biogeography](./3_biogeography): scripts related to the ancestral range estimation based on the dated species tree

### PART C
The script `partc_tfsp_correlation.R` corresponds to the analysis in "Part C" of the thesis, which aimed to test the Time-for-Speciation hypothesis in Madagascan palms. The correlation between species richness and arrival time was tested for different datasets while sampling over all possible arrival time intervals.


### General
The repository needs to be cloned into the folder where you want to have the data generated. The folder containing all scripts should be renamed to "scripts".
All scripts should be run from the folder `../scripts`. For example: `./scripts/1_phylo_reconstruction/0.trimming.sh` to trimm reads. For details on how to run the scripts and script flags, see in the corresponding folders.