#!/bin/bash
###########################################################################
# Project: Orania Phylogeny MT
# Script: alignment_script.sh
# --- Action: This script aligns the contigs using Mafft and Gblocks.
# ----------- We used the default parameters of Mafft and Gblocks.
# --- Input: assembled contigs from script "hybpiper_script.sh"
# --- Output: aligned contigs
# Author: Maya Schroedl (maya.schroedl@bios.au.dk)

###########################################################################

###################
#----ARGUMENTS----#
###################

#t:if you want trimmed contigs, based on coverage and threshold, enter t:threshold

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
  esac

done


#################################
#----DIRECTORIES & VARIABLES----#
#################################

#working directories
GWD=$PWD #global working directory, with subprojects and scripts
WD="$GWD"/1_phylo_reconstruction #current working directory


if [ -z "$t" ]; #if there is no threshold definded, then look at the normal files
then
        FILES=$WD/1.0_hybpiper/supercontigs
else #else, have a look at the files trimmed according to coverage maps and threshold t
        FILES=$WD/1.1_coverage_maps/trimmed_contigs/$t
        dir_value=$t/
fi

#make a new directory
mkdir -p $WD/2_alignment/$dir_value

rm -rf $WD/2_alignment/$dir_value* #if directory not empty: clear directory

##########################
#----FILE PREPARATION----#
##########################

#remove carriage return in genelist_7575_X.txt:
sed -i 's/\r$//g' $WD/genelist_7575.txt

#combine contigs into one file
while read gene; do
cat $FILES/"$gene"_supercontig.fasta >> $WD/2_alignment/$dir_value$gene"_combined.fasta"
done < $WD/genelist_7575.txt

###################
#----ALIGNMENT----#
###################

#### MAFFT ####
while read gene;
do mafft --auto $WD/2_alignment/$dir_value$gene"_combined.fasta" > $WD/2_alignment/$dir_value$gene"_aligned.fasta";
done < $WD/genelist_7575.txt

#### GBLOCKS ####
while read gene;
do $GWD/programs/Gblocks_0.91b/Gblocks $WD/2_alignment/$dir_value$gene"_aligned.fasta" -t=d -b5="a"; 
 done < $WD/genelist_7575.txt
#b5=a: all gap positions allowed

cd $WD/2_alignment/$dir_value
rename 's/.fasta-gb/_gb.fasta/' * #sudo apt-get install rename

### copy original alignements to new files which will be modified manually
for file in $WD/2_alignment/$dir_value*_gb.fasta; 
do cp "$file" ${file/.fasta}_mod.fasta;
done

#alignments were manually corrected using AliView v 1.26 because some sequences were misaligned.



rm *.htm
rm *combined.fasta


cd $GWD
