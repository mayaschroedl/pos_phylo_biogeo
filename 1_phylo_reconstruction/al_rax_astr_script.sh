#!/bin/bash
###########################################################################
# Project: Orania Phylogeny MT
# Script: al_rax_astr_script
# --- Action: This script combines the alignment script, the raxml and the
# ----------- astral script.
# --- Input: assembled contigs from script "hybpiper_script.sh"
# --- Output: coalescent tree
# Author: Maya Schroedl
# Date: 11/2019
###########################################################################

#set working directory to current directory
WD=$PWD

#f: which files do you want to work on #foldernames seperated by commata
#t: if you want trimmed contigs, based on coverage and threshold, enter t:threshold
#d: filename of the doublespecies file in $WD
while getopts ":f:d:t:" opt; do
  case ${opt} in
    f)
      echo "-f was triggered, File: $OPTARG"
      FILENAMES=${OPTARG}
      ;;
    d)
      echo "-d was trigerred, File: $OPTARG"
      DOUBLESP=${OPTARG}
      ;;
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

#make directory
mkdir -p $WD/2_alignment
mkdir -p $WD/3_gene_trees
mkdir -p $WD/4_coalescent_trees

#make directory
IFS=" "
mkdir -p $WD/2_alignment/$FILENAMES$t
mkdir -p $WD/3_gene_trees/$FILENAMES$t
mkdir -p $WD/4_coalescent_trees/$FILENAMES$t

#align
$WD/scripts/2_alignment_script.sh -f $FILENAMES -t $t

#make gene trees
$WD/scripts/3_raxml_script.sh -f $FILENAMES -t $t

#make coalescent tree
$WD/scripts/4_astral_script.sh -f $FILENAMES -d $DOUBLESP -t $t

cd $WD