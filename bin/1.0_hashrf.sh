#!/bin/bash
###########################################################################
# Project: Orania Phylogeny MT
# Script: hashrf.sh
# --- Action: This script evaluates the simiarity between each gene tree
# ------------ and the species tree.
# --- Input: gene trees; species tree
# --- Output: matrix of RF values
# Author: Maya Schroedl

###########################################################################
#install newest version of https://code.google.com/archive/p/hashrf/downloads

###################
#----ARGUMENTS----#
###################

#t:if you want trimmed contigs, based on coverage and threshold, enter t:threshold

while getopts ":t:" opt; do
    t)
      echo "-t was triggered, Threshold: $OPTARG"
      t=$OPTARG
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
WD_phylo="$GWD"/1_phylo_reconstruction #phylo reconstruction directory
WD="$GWD"/2_phylo_dating #dating directory

if ! [ -z "$t" ]; #if there is a threshold definded, then work in the threshold directory
then
        dir_value=$t/
fi

mkdir -p $WD/1_hashrf/$dir_value


##########################
#----FILE PREPARATION----#
##########################

#put astral tree and gene trees in one file 
cat $WD_phylo/4_coalescent_trees/"$dir_value"coalescent_lpp_rooted.tree $WD_phylo/3_gene_trees/"$dir_value"4_collapsed/all_genes_rooted.raxml.support.coll > $WD/1_hashrf/$dir_value"genetrees_astral_combi_rooted.tre"

################
#----HASHRF----#
################

#hashrf calculates the distance between every tree and makes a matrix of RF values
hashrf $WD/1_hashrf/$dir_value"genetrees_astral_comb_rooted.tre" 0 -o $WD/1_hashrf/$dir_value"RF_matrix_rooted.txt"

###########################################
#----SELECT GENES WITH LOWEST RF VALUE----#
###########################################
#RF_sort.R ?

#output: list of selected genes: $WD/1_hashrf/RF_selgenes.txt


###################################################
#----TAKE SELECTED ALIGNMENTS & ADD EMPTY SEQS----#
###################################################
#make a new folder with alignements corresponding to genes selected by RF
mkdir -p $WD_phylo/2_alignment/"$dir_value"RF_selected

while read gene;
do gene=$(echo "$gene" | tr -d '\r') #need to remove carriage return
cp $WD_phylo/2_alignment/"$dir_value"$gene"_aligned_gb.fasta" $WD_phylo/2_alignment/"$dir_value"RF_selected;
done < $WD/1_hashrf/RF_selgenes.txt

#---ADD EMPTY SEQUENCES FOR MISSING TAXA----#

#run 7.1_biopython_RF_sp_add.py on windows (or modify script)

#output: alignemnt with empty sequences added

################################
#----CONCATENATE ALIGNMENTS----#
################################

$WD_phylo/2_alignment/"$dir_value"RF_selected
cat *_sp_added.fasta|awk -v RS=">" -v FS="\n" -v OFS="\n" '{for(i=2; i<=NF; i++) {seq[$1] = seq[$1]$i}}; END {for(id in seq){print ">"id, seq[id]}}' > "RF_selected_concat_alignment.fasta"

