#!/bin/bash
###########################################################################
# Project: Orania Phylogeny MT
# Script: hashrf.sh
# --- Action: Concatenated the selected gene alignments for dating
# --- Input: selected genes
# --- Output: concatenated alignments
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
WD="$GWD"/2_phylo_dating #dating directory


mkdir -p $WD/1_alignment/concat
mkdir -p $WD/3_gene_shop/selected_genes

################################
#----CONCATENATE ALIGNMENTS----#
################################

cd $WD/1_alignment/

#concatenate all genes

pxcat -s $WD/1_alignment/*_gb_mod.fasta -o $WD/2_alignment/concat/all_genes_concat.fasta #install pxcat
#### AFTER CONCATENATING REMOVE THE INDIVIDUAL SEQUENCES WHICH WERE SELECTED FOR DATING (SEE astra_dropped_tags.txt)



########## select genes ###############
# with script gene_shopping_dat.R

### concatenate selected genes

concat() {
	filename=$1
	cp $WD/3_gene_shop/selected_genes/$filename.txt $WD/3_gene_shop/selected_genes/"$filename"_files.txt
	sed -i -e "s#^#/data_vol/maya/2_phylo_dating/1_alignment/##" $WD/3_gene_shop/selected_genes/"$filename"_files.txt
	sed -i -e "s/$/_aligned_gb_mod.fasta/" $WD/3_gene_shop/selected_genes/"$filename"_files.txt
	tr -d "\r" < $WD/3_gene_shop/selected_genes/"$filename"_files.txt  >  $filename.txt_ch  && mv  $filename.txt_ch $WD/3_gene_shop/selected_genes/"$filename"_files.txt 

	selected_genes_str=$(cat $WD/3_gene_shop/selected_genes/"$filename"_files.txt  | tr "\r" " ")
	
	echo $selected_genes_str
	
	pxcat -s $selected_genes_str -o $WD/1_alignment/concat/"$filename"_concat.fasta
}


rm $WD/1_alignment/concat/*
concat no_gdis
concat most_diff
concat clock_one
concat clock_three
concat clock_nine

#### AFTER CONCATENATING REMOVE THE INDIVIDUAL SEQUENCES WHICH WERE SELECTED FOR DATING (SEE astra_dropped_tags.txt)



