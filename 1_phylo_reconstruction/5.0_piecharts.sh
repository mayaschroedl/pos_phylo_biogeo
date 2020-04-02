#!/bin/bash
###########################################################################
# Project: Orania Phylogeny MT
# Script: phyparts.sh
# --- Action: This script generates and vizualises phyparts on sp trees as 
# ----------- piecharts.
# --- Input: gene trees; species tree
# --- Output: species tree with piecharts
# Author: Maya Schroedl
# Date: 11/2019
###########################################################################

###################
#----ARGUMENTS----#
###################

#t: if you want trimmed contigs, based on coverage and threshold, enter t=threshold
#o: one outgroup
while getopts ":t:o:" opt; do
  case ${opt} in
    t)
      echo "-t was triggered, Threshold: $OPTARG"
      t=$OPTARG
      ;;
	  
	o)
	 echo "-o was triggered, Outgroup: $OPTARG"
	 outgroup=$OPTARG
	 ;;

    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac

done

###################################
#----DIRECTORIES AND VARIABLES----#
###################################

#working directories
GWD=$PWD #global working directory, with subprojects and scripts
WD="$PWD"/1_phylo_reconstruction #current working directory

if ! [ -z "$t" ]; #if there is a threshold definded, then work in the threshold directory
then
        dir_value=$t/
fi

mkdir -p $WD/3_gene_trees/$dir_value"rooted"
mkdir -p $WD/5_phypartstopiecharts/

######################
#----REROOT TREES----#
######################

#----REROOT GENE TREES----#
cd $WD/3_gene_trees/$dir_value

for i in RAxML_parsimonyTree.*; do nw_reroot $i $outgroup > rooted/${i/.tre}_rooted.tre; done 

cat $WD/3_gene_trees/rooted/RAxML_parsi* > $WD/3_gene_trees/"$dir_value"/rooted/RAxML_combined_rooted_gene.tre

#----REROOT SPECIES TREE----#
nw_reroot $WD/4_coalescent_trees/"$FILENAMES$t"_coalescent_indiv_diff.tre $outgroup > $WD/4_coalescent_trees/"$FILENAMES$t"_coalescent_indiv_diff_rooted.tre | nw_topology -


##################
#----PHYPARTS----#
##################

java -jar $WD/programs/phyparts/target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 1 -v -d $WD/3_gene_trees/rooted/raxml_parsimony_"$FILENAMES$t"_gene_rooted.tre -m $WD/4_coalescent_trees/"$FILENAMES$t"_coalescent_indiv_diff_rooted.tre -o $WD/5_phypartstopiecharts/"$FILENAMES$t"

###############################
#----PHYPARTS 2 PIE CHARTS----#
###############################
#phyparts to pie charts
display=$(shuf -i 100-200 -n 1)
export DISPLAY=:${display}
Xvfb :${display} -screen 0 1024x768x16 > /dev/null 2>&1 &
echo "export DISPLAY=:${display}" > ~/.xvfb

python $WD/scripts/phypartstopiecharts.py $WD/4_coalescent_trees/"$FILENAMES$t"_coalescent_indiv_diff_rooted.tre $WD/5_phypartstopiecharts/"$FILENAMES$t" 174 --svg_name $WD/5_phypartstopiecharts/"$FILENAMES$t"_phypartspiecharts.svg


