#!/bin/bash
###########################################################################
# Project: Orania Phylogeny MT
# Script: hashrf.sh
# --- Action: Concatenated the selected gene alignments for dating
# --- Input: selected genes
# --- Output: concatenated alignments
# Author: Maya Schroedl
# Date: 11/2019
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

if ! [ -z "$t" ]; #if there is a threshold definded, then work in the threshold directory
then
        dir_value=$t/
fi

mkdir -p $WD/3_all_genes/$dir_value


################################
#----CONCATENATE ALIGNMENTS----#
################################

cd $WD/2_alignment/"$dir_value"

python $GWD/scripts/1_phylo_dating/Gene_Sticher.py -in *_gb.fasta

cat *_gb.fasta|awk -v RS=">" -v FS="\n" -v OFS="\n" '{for(i=2; i<=NF; i++) {seq[$1] = seq[$1]$i}}; END {for(id in seq){print ">"id, seq[id]}}' > $WD/3_all_genes/$dir_value"all_genes_concat.fasta"

