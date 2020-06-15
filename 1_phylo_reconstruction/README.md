# Reconstructing phylogenies based on targeted sequencing data

These scripts document the whole process of going from raw reads towards a species phylogeny.

First, you should create a new folder for the project, where the scripts will be run from. Create a directory named `1_phylogeny_reconstruction/input_reads_and_info` and place your read files in this directory. Each read file should contain all read sequences for one individual. Each read file should be named

<p style="text-align: center;">
`TAG-x_R*_y_z.fastq`
</p>

where x should be a number that corresponds to the individual
* should be 1 or 2, standing for either the reversed or forward read
y should be "unpaired" or "paired" (if unpaired, remove the `_R*_` part)
z should be "trimmed" or "" (empty), depending on whether trimming was alread done or not.

In addition to that, the folder should contain a fasta file named `target_file.fasta` which contains all the sequences targeted by the used baits for target capture.

The folder should also contain a file called `tags_indiv.txt` which contains a table that maps the "TAGs" (column: tag) to the "individual names" (i.e. species name underscore individual number) (column: indiv).

That's it. Now the tree building can begin. Make sure that within your working directory you have a folder named "scripts" where this repository was cloned in. Scripts should be run from you working directory.

All programs should be installed in a folder in the working directory called "programs"


### **`0_trimming.sh`**

This script trimms all samples which were not yet trimmed with Trimmomatic. You need to have installed [Trimmomatic v0.39](http://www.usadellab.org/cms/?page=trimmomatic). 

We chose to use the following settings: LEADING:30 TRAILING:30 SLIDINGWINDOW:4:30 MINLEN:36 AVGQUAL:30. For details see the [Trimmomatic website]((http://www.usadellab.org/cms/?page=trimmomatic)).

[FASTQC v0.11.9](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) needs to be installed for Quality control before and after trimming.

### **`1_hybpiper.sh`**

This script generates exon, intron, and supercontig sequences using [HybPiper] (https://github.com/mossmatters/HybPiper), which needs to be installed beforehand. For dependencies and methods see the [HybPiper Wiki](https://github.com/mossmatters/HybPiper/wiki). We used paired and unpaired trimmed reads as input.

##### **`1.0.1_gene_recovery_heatmap_ggplot.R`**

[HybPiper] (https://github.com/mossmatters/HybPiper) script that outputs a histogram representing the amount of recovery for each gene and individual.

### **`1.1_gene_select_75_75.R`**

In order to avoid gene tree errors due to fragmented sequences and to limit the amount of missing data, only this script selects the genes that had 75% of their target sequence length recovered in 75% of the individuals.

### **`1.2_coverage_maps.sh`**

This script computes a "coverage map" for each supercontig using [bwa](http://bio-bwa.sourceforge.net/) and [samtools](http://www.htslib.org/). A coverage map represents how many reads support each position over a sequence.

##### **`1.2.1_visualize_covmaps_stats.R`**

This script can be used to visualize the coverage maps produced with `1.2_coverage_maps.sh` and get summary statistics over all coverage maps.

### **`1.3_cov_trimm_contigs.sh`**

##### **`1.3.1_reorga_coverage_to_genes.py`**

##### **`1.3.2_cov_trimm_contigs.py`**

### **`2_alignment_script.sh`**

##### **`2.1_alignment_stats.R`**

### **`3_raxml_script.sh`**

##### **`3.1_collapse_low_brnlen.R`**

### **`4_astral_script.sh`**

### **`5_piecharts.sh`**

##### **`5.1_phypartstopiecharts.py`**

##### **`5.2_phypart_alternative_topo_viz.R`**
