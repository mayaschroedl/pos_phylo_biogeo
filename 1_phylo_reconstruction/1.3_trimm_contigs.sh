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
WD="$PWD"/1_phylo_reconstruction #current working directory

mkdir -p $WD/1.1_coverage_maps/trimmed_contigs/$t

###############
#----TRIMM----#
###############

rm -rf $WD/1.1_coverage_maps/trimmed_contigs/$t/* #if already exists: remove

cd $WD
python $GWD/scripts/1.2.1_trimm_contigs.py $t #t=threshold;

rm -r $WD/1.1_coverage_maps/trimmed_contigs/$t/contigs

