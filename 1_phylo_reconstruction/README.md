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

All programs should be installed in a folder in the working directory called "programs".

Python v2.7 and R3.6.2 were used.


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

This script computes a "coverage map" for each supercontig using [bwa](http://bio-bwa.sourceforge.net/) and [samtools](http://www.htslib.org/). A coverage map represents how many reads support each position over a sequence. The output is one coverage map per TAG (i.e. individual).

##### **`1.2.1_visualize_covmaps_stats.R`**

This script can be used to visualize the coverage maps produced with `1.2_coverage_maps.sh` and get summary statistics over all coverage maps.

### **`1.3_cov_trimm_contigs.sh`**

In order to minimize the probability of sequencing errors, this script maskes the positions on a supercontig that present low read coverage. A threshold of read coverage (e.g. 2) can be given with the flag `-t` (e.g. -t 2). All positions each contig which are <= t, will be masked with a N and the read endings which contain a chain of N (NNNN) will be trimmed. A new folder will be created within the directory "coverage_maps", where the trimmed contigs will be placed. The folder name will be the same as the set threshold.

It is based on the two following scripts (two steps for coverage trimming):

##### **`1.3.1_reorga_coverage_to_genes.py`**
With the script `1.2_coverage_maps.sh`, one coverage map per TAG was generated. Now, these coverage maps were restructured into one file per gene, because the contig files were structured into genes. This restructering was necessary for the following script:

##### **`1.3.2_cov_trimm_contigs.py`**
This script does the actual trimming work and replaced each position under a certain threshold t with an "N" and trimms the endings. The threshold should be the first argument (e.g. `1.3.2_cov_trimm_contigs.py 2` for a threshold of 2) The input are coverage maps and contigs per gene, and the script outputs a masked and trimmed contig per gene.

### **`2_alignment_script.sh`**
This script aligns the contigs using [MAFFT](https://mafft.cbrc.jp/alignment/software/) v7.4 and cleans the alignments with [Gblocks](http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_documentation.html) v0.91b. Both programms should be previously installed. Gblocks was run with the option ("-b5=a"), meaning that "all gap positions can be selected. Positions with gaps are not treated differently from other positions" (cf. [Gblocks documentation](http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_documentation.html)). To select the right contigs (i.e. the right folder name, which is called by the threshold), one can specify the flag `-t`, indicating the threshold which was chosen for coverage trimming. The flag is optional and can be left out if no coverage trimming was done.

Alignments should be manually cleaned after running this script with e.g. [AliView](https://ormbunkar.se/aliview/) to remove poorly aligned, short sequence bits.

##### **`2.1_alignment_stats.R`**

Get general alignment statistics, like mean alignement length and percentage of missing data over all alignments.

### **`3_raxml_script.sh`**
Generate gene trees with [RAxML-ng](https://github.com/amkozlov/raxml-ng) v0.9 using bootstrapping for support. Also collapses branches with low bootstrap support (<10%) using [Newick Utilities](http://cegg.unige.ch/newick_utils) and branches with low branchlengths (< 0.00001) using the script `3.1_collapse_low_brnlen.R`. Rooted gene trees are outputted and also gene trees with their labels changed from TAGs to individual names.

RAxML-ng and Newick Utilities should be previously installed.

The flags for this script are:

* `-t`: indicating the threshold which was chosen for coverage trimming. The flag is optional and can be left out if no coverage trimming was done.
* `-o`: name of the outgroup species (here: only one species was used) in TAG format, for rooting
* `-l`: name of the outgroup species (here: only one species was used) in individual name format, for rooting on labeled trees

An example to run the script from the project working directory is: `./scripts/1_phylo_reconstruction/3_raxml_script.sh -t 2 -o TAG-47 -l Dypsis_mananjarensis`

##### **`3.1_collapse_low_brnlen.R`**

This short R script collapses clades that present a branchlengths < 0.00001 using the package "ape". This was done because it might be that 200 bootstrap replicates are not enough to assign a low bootstrap value for these very short branches.

### **`4_astral_script.sh`**
This script generates a coalescent tree of individuals and a species tree, where individuals are forced into monophyletic species, using [ASTRAL-III](https://github.com/smirarab/ASTRAL/). Support values were generated as quartet support and local posterior probabilities.
Moreover, to test the monophyly of individuals of the same species, in a second analysis, a monophyletic constraint was done on these individuals (option "-j", cf. [ASTRAL Constrained-search github](https://github.com/maryamrabiee/Constrained-search)). 

The flags for this script are:

* `-t`: indicating the threshold which was chosen for coverage trimming. The flag is optional and can be left out if no coverage trimming was done.
* `-o`: name of the outgroup species (here: only one species was used) in TAG format, for rooting
* `-l`: name of the outgroup species (here: only one species was used) in individual name format, for rooting on labeled trees


### **`5_piecharts.sh`**

In order to get gene tree conflict, this script used [PhyParts](https://bitbucket.org/blackrim/phyparts/src/master/) to generate statistics on gene tree conflict for each node, based on bipartitions. The script `5.1_phypartstopiecharts.py` (by [Matt Johnson](https://github.com/mossmatters/phyloscripts/tree/master/phypartspiecharts)) visualizes this gene tree conflict. For details, have a look at his well-documented github tutorial.

##### **`5.1_phypartstopiecharts.py`**
This is Matt Johnsons script to visualize gene tree conflict based on the PhyPart results. This script is taken from https://github.com/mossmatters/phyloscripts/tree/master/phypartspiecharts and was not modified. All credits go to Matt Johnson.

##### **`5.2_phypart_alternative_topo_viz.R`**
This R script visualizes the first alternative topology for each node with the results from [PhyParts](https://bitbucket.org/blackrim/phyparts/src/master/). A manual step before running this script is needed: the out.hist.alts of the PhyParts output needs to be manually transformed into table, and the included trees transformed into newick format.This new file should be named "out.hist.alts.sum.trees".
