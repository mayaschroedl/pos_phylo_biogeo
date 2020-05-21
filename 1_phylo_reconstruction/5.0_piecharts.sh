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

#################################
#----DIRECTORIES & VARIABLES----#
#################################

#working directories
GWD=$PWD #global working directory, with subprojects and scripts
WD="$GWD"/1_phylo_reconstruction #current working directory

if ! [ -z "$t" ]; #if there is a threshold definded, then work in the threshold directory
then
        dir_value=$t/
fi

mkdir -p $WD/5_phypartstopiecharts/$dir_value


#remove spaces from species tree (because this might be a problem with phyparts)
sed -i 's/\s*$//' $WD/4_coalescent_trees/"$dir_value"coalescent_qs.tree_lab_rooted


##################
#----PHYPARTS----#
##################

cd  $WD/5_phypartstopiecharts/$dir_value

java -jar $GWD/programs/phyparts/target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -a 1 -v -d $WD/3_gene_trees/"$dir_value"4_collapsed/all_genes.raxml.support.coll_lab_rooted -m $WD/4_coalescent_trees/"$dir_value"coalescent_lpp.tree_lab_rooted

###############################
#----PHYPARTS 2 PIE CHARTS----#
###############################
#phyparts to pie charts
display=$(shuf -i 100-200 -n 1)
export DISPLAY=:${display}
Xvfb :${display} -screen 0 1024x768x16 > /dev/null 2>&1 &
echo "export DISPLAY=:${display}" > ~/.xvfb

#number of genes
gene_num() { wc -l < $WD/genelist_7575.txt ;} 
gene_n=$(gene_num)

#execute Matt Johnson's script
python $GWD/scripts/1_phylo_reconstruction/5.1_phypartstopiecharts.py $WD/4_coalescent_trees/"$dir_value"coalescent_lpp.tree_lab_rooted $WD/5_phypartstopiecharts/$dir_value $gene_n --svg_name $WD/5_phypartstopiecharts/$dir_value/phypartspiecharts_lab.svg 

cd $GWD