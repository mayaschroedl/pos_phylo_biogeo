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

#t: if you want trimmed contigs, based on coverage and threshold, enter t=threshold
while getopts ":t:" opt; do
  case ${opt} in
    t)
      echo "-t was triggered, File: $OPTARG"
      t=$OPTARG
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

mkdir -p $WD/4_coalescent_trees/$dir_value

################
#----ASTRAL----#
################

cd $WD/4_coalescent_trees/$dir_value
rm *.tre
rm *.log

java -jar $WD/programs/Astral/astral.5.6.3.jar -t 1 -i $WD/3_gene_trees/"$dir_value"RAxML_combined_gene.tre -o "$dir_value"coalescent.tre 2>"$dir_value"coalescent.log -a $WD/renamed_reads/sp_TAGs.txt
#sp_TAGs.txt: doublesp file - which TAGs correspond to which species names.


#bootstrap
#java -jar $WD/programs/Astral/astral.5.6.3.jar -t 1 -i $WD/3_gene_trees/"$FILENAMES$t"/raxml_parsimony_"$FILENAMES$t"_gene.tre -b bs-files -o "$FILENAMES$t"_coalescent.tre 2>"$FILENAMES$t"_coalescent.log  -a $WD/4_coalescent_trees/"$FILENAMES$t"/$DOUBLESP 
#make bs-files

#new astral tree, where individuals are treated differently
cd $WD/4_coalescent_trees/"$FILENAMES$t"

java -jar $WD/programs/Astral/astral.5.6.3.jar -t 1 -i $WD/3_gene_trees/raxml_parsimony_gene.tre -o coalescent_indiv_diff.tre 2>coalescent_indiv_diff.log


cd $WD