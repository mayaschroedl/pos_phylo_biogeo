#!/bin/bash
###########################################################################
# Project: Orania Phylogeny MT
# Script: alignment_script.sh
# --- Action: This script builts a coalescent tree using Astral.
# --- Input: genetrees from raxml; file with species names
# --- Output: coalescent tree
# Author: Maya Schroedl
# Date: 11/2019
###########################################################################


###################
#----ARGUMENTS----#
###################


#t:if you want trimmed contigs, based on coverage and threshold, enter t:threshold
#o: outgroup (TAG)
#l: labeled outgroup (individual name)
while getopts ":t:o:l:" opt; do
  case ${opt} in
    t)
      echo "-t was triggered, File: $OPTARG"
      t=$OPTARG
      ;;
o)
      echo "-o was triggered, Outgroup: $OPTARG"
      outgroup=$OPTARG
      ;;
	  
	  l)
      echo "-l was triggered, labeled Outgroup: $OPTARG"
      outgroup_lab=$OPTARG
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

mkdir -p $WD/4_coalescent_trees/$dir_value

################
#----ASTRAL----#
################

cd $WD/4_coalescent_trees/$dir_value
rm *

#ASTRAL tree with tips representing individuals (individuals of same species are not merged into one species tip)

#execute with branch support as local posterior probabilities (lpp) [default]
java -jar $GWD/programs/Astral/astral.5.6.3.jar -i $WD/3_gene_trees/"$dir_value"4_collapsed/all_genes.raxml.support.coll -o $WD/4_coalescent_trees/"$dir_value"coalescent_lpp.tree 2>$WD/4_coalescent_trees/"$dir_value"coalescent_lpp.log

#execute with branch support as all quartet supports (qs) [-t 8]
java -jar $GWD/programs/Astral/astral.5.6.3.jar -t 8 -i $WD/3_gene_trees/"$dir_value"4_collapsed/all_genes.raxml.support.coll -o $WD/4_coalescent_trees/"$dir_value"coalescent_qs.tree 2>$WD/4_coalescent_trees/"$dir_value"coalescent_qs.log

#bootstrap
#java -jar $WD/programs/Astral/astral.5.6.3.jar -i $WD/3_gene_trees/"$dir_value"4_collapsed/all_genes.raxml.support.coll  -b $WD/3_gene_trees/"$dir_value"2_bootstrap  -o $WD/4_coalescent_trees/"$dir_value"coalescent_bstrp.tree 2>$WD/4_coalescent_trees/"$dir_value"coalescent_bstrp.log

######################
#----LABEL TREES----#
######################
Rscript $GWD/scripts/general/change_tiplabels.R $WD/4_coalescent_trees/"$dir_value"coalescent_lpp.tree $WD/4_coalescent_trees/"$dir_value"coalescent_lpp.tree_lab $WD/renamed_reads/tags_indiv.txt

Rscript $GWD/scripts/general/change_tiplabels.R $WD/4_coalescent_trees/"$dir_value"coalescent_qs.tree $WD/4_coalescent_trees/"$dir_value"coalescent_qs.tree_lab $WD/renamed_reads/tags_indiv.txt

######################
#----REROOT TREES----#
######################

#for unlabeled
Rscript $GWD/scripts/general/root_tree.R $WD/4_coalescent_trees/"$dir_value"coalescent_lpp.tree $WD/4_coalescent_trees/"$dir_value"coalescent_lpp.tree_rooted $outgroup
 
Rscript $GWD/scripts/general/root_tree.R $WD/4_coalescent_trees/"$dir_value"coalescent_qs.tree $WD/4_coalescent_trees/"$dir_value"coalescent_qs.tree_rooted  $outgroup

#for labeled
Rscript $GWD/scripts/general/root_tree.R $WD/4_coalescent_trees/"$dir_value"coalescent_lpp.tree_lab $WD/4_coalescent_trees/"$dir_value"coalescent_lpp.tree_lab_rooted $outgroup_lab
 
Rscript $GWD/scripts/general/root_tree.R $WD/4_coalescent_trees/"$dir_value"coalescent_qs.tree_lab $WD/4_coalescent_trees/"$dir_value"coalescent_qs.tree_lab_rooted  $outgroup_lab



cd $GWD