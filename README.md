# Colonization and diversification history of Madagascan palms with new phylogenomic evidence from the genus *Orania* (Arecaceae)


This repository contains all scripts and supplementary material generated during the preparation of my Master's thesis. The file `Master_Thesis_Maya_Schrodl.pdf` contains the submitted .pdf version of this Master's thesis. If you have any questions regarding the thesis or the scripts, feel free to contact me: [maya.schroedl@bios.au.dk](mailto:maya.schroedl@bios.au.dk)


This thesis mainly dealt with understanding colonization and diversification patterns in Madagascan palms using phylogenomics. It is structured in three parts: Part A, B, and C (see below for description of each part). Part A was just a general introduction and did therefore inlude no scripts.


The supplementary materials to the Part B and Part C of the thesis can be found under [/supplementary_material](./supplementary_material).


### PART B: Testing the 'Time-for-Speciation hypothesis' in Madagascan palms
This part aimed to understand whether differences in species richness between Madagascan palm groups could be explained only by the time they have been present on the island. 

An interval of possible arrival times and species richness for each group was constructed, based on available species-level phylogenies and the [World Checklist of Selected Plant Families (WCSP)](https://wcsp.science.kew.org/qsearch.do). Second, it was tested whether a correlation between species richness and arrival time could be found.

##### **`partc_tfsp_correlation.R`**
The script `partc_tfsp_correlation.R` test how "probable" a significant positive/negative correlation between species richness (ln) and arrival time is when sampling randomly on all possible arrival time intervals (on a uniform distribution). Therefore, for each species group, an arrival time was randomly sampled and the correlation was tested. This was repeated for 100,000 samples and then the proportion of positive/negative correlations was calculated. In addition to that, the script outputs a plot of species richness (ln) as a function of arrival time with all intervals plotted. Additionally, a plot with all slopes of the correlations over all repetitions is outputted.

### PART C: Phylogenomics & Biogeography of the genus *Orania*

The folders [/1_phylo_reconstruction](./1_phylo_reconstruction), [/2_phylo_dating](./2_phylo_dating), and [/3_biogeography](./3_biogeography) correspond to the "Part C" of the Master's thesis, which aimed to reconstruct a dated phylogeny of the POS clade and apply  ancestral range estimations to this phylogeny. Each folder contains a detailed description of the methods applied and how to run the scripts.

* [/1_phylo_reconstruction](./1_phylo_reconstruction): scripts related to the reconstruction of gene and species trees, based on targeted sequencing data

* [/2_phylo_dating](./2_phylo_dating): scripts related to dating of the previously estimated species tree

* [/3_biogeography](./3_biogeography): scripts related to the ancestral range estimation based on the dated species tree

### General
The folder [/general](./general) contains scripts that were useful for many phylogenetic analyses, like e.g. rooting and visualizing trees. Detailed description of the scripts can be found in this folder.

### How to run scripts
The repository needs to be cloned into the folder where you want to have the data generated. The folder containing all scripts should be renamed to "scripts".
All scripts should be run from the folder `../scripts`. For example: `./scripts/1_phylo_reconstruction/0.trimming.sh` to trimm reads. For details on how to run the scripts and script flags, see in the corresponding folders.