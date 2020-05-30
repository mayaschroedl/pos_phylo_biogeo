#!/bin/bash
###########################################################################
# Project: Orania Phylogeny MT
# Script: alignment_script.sh
# --- Action: This script aligns the contigs using Mafft and Gblocks.
# ----------- We used the default parameters of Mafft and Gblocks.
# --- Input: assembled contigs from script "hybpiper_script.sh"
# --- Output: aligned contigs
# Author: Maya Schroedl

###########################################################################

#remove files not listed in list

for f in $(cat $WD/1_hybpiper/or,orscldyp,sclpod_not7575.txt) ; do 
  find . -type f -name "$f""_" -delete
done