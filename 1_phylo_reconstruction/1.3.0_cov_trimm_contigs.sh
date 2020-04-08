#!/bin/bash
# -*- coding: utf-8 -*-

###########################################################################
# Project: Orania Phylogeny MT
# Script: trimm_contigs.sh
# --- Action: This script replaces the bases that have a coverage under 
# ----------- a certain threshold “t” with an “N” and then trimms the ends
# ----------- of the gene.
# --- Input:  contigs and coverage maps for each gene, for different
# ----------- species families ($ID)
# --- Output: trimmed contigs for each speciesl family ($ID)
# Author: Maya Schroedl
# Date: 11/2019
###########################################################################


####################
#----ARGUMENTS----#
####################

#t: trimming threshold
while getopts ":t:" opt; do
  case ${opt} in
	t)
      echo "-t was triggered, Threshold: $OPTARG"
      t=$OPTARG
      ;;

    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
  
  done

#####################
#----DIRECTORIES----#
#####################

#working directories
GWD=$PWD #global working directory, with subprojects and scripts
WD="$GWD"/1_phylo_reconstruction #current working directory

mkdir -p $WD/1.1_coverage_maps/per_gene
mkdir -p $WD/1.1_coverage_maps/trimmed_contigs/$t/contigs

#copy supercontigs to mkdir -p $WD/1.1_coverage_maps/trimmed_contigs/$t/contigs
rm $WD/1.1_coverage_maps/trimmed_contigs/$t/contigs/*
while read gene;
do cp $WD/1.0_hybpiper/supercontigs/$gene* $WD/1.1_coverage_maps/trimmed_contigs/$t/contigs
done < $WD/genelist_7575.txt

cd $WD/1.1_coverage_maps/trimmed_contigs/$t/contigs
#rename header of each  file (*$gene.fasta) from TAGX-geneX to TAGX
for file in *.fasta; do sed -i 's/-[^-]*//2g' $file; done #


##########################
#----FILE PREPARATION----#
##########################
#----PROBLEM----#
#We generated coverage maps per TAG. But we want the trimmed conitgs organized in genes. 
#----SOLUTION---#
#--> we have to reorganize the coverage maps per gene

cd $WD
python $GWD/scripts/1_phylo_reconstruction/1.3.1_reorga_coverage_to_genes.py $t #t=threshold;


###############
#----TRIMM----#
###############

cd $WD
python $GWD/scripts/1_phylo_reconstruction/1.3.2_cov_trimm_contigs.py $t #t=threshold;

rm -r $WD/1.1_coverage_maps/trimmed_contigs/$t/contigs

cd $GWD

