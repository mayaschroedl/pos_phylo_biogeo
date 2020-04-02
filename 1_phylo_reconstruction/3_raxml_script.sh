#!/bin/bash
###########################################################################
# Project: Orania Phylogeny MT
# Script: raxml_script.sh
# --- Action: This script creates gene trees with RaxML
# --- Input: aligned contigs
# --- Output: gene trees
# Author: Maya Schroedl
# Date: 11/2019
###########################################################################

###################
#----ARGUMENTS----#
###################

#t:if you want trimmed contigs, based on coverage and threshold, enter t:threshold
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

#make directory
mkdir -p $WD/3_gene_trees/$dir_value

#########################
#----MAKE GENE TREES----#
#########################
# with GTRGAMMA model

cd $WD/3_gene_trees/$dir_value
rm *
# while read gene; 
# do raxml -s $WD/2_alignment/$dir_value"$gene"_aligned_gb.fasta -m GTRGAMMA -n "$gene"_gene.tre -p 12345;
# done < $WD/1_hybpiper/genelist_7575.txt

#combine all raxml genetrees in one file
# cat $WD/3_gene_trees/"$dir_value"RAxML_parsi* > $WD/3_gene_trees/"$dir_value"RAxML_combined_gene.tre

# #raxml-ng + bootstrap
while read gene;
 do raxml-ng --check --msa $WD/2_alignment/$dir_value"$gene"_aligned_gb.fasta --model GTRGAMMA; 
 do raxml-ng --msa $WD/2_alignment/$dir_value"$gene"_aligned_gb.fasta --model GTRGAMMA --seed 2;
 do raxml-ng --bootstrap --msa $WD/2_alignment/$dir_value"$gene"_aligned_gb.fasta --model GTRGAMMA --seed 2 --bs-trees 200;
done < $WD/1_hybpiper/genelist_7575.txt

#check output
#collapse all branches with bootstrap <10
#module load bioinfo/newick-utils/1.6
#$WD/3_gene_trees/"$FILENAMES$t"/raxml_parsimony_"$FILENAMES$t"_gene.tre 'i & b<=10' o > $WD/3_gene_trees/"$FILENAMES$t"/raxml_parsimony_"$FILENAMES$t"_gene_bs10.tre 


cd $WD




