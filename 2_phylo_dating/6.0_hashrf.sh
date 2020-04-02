#!/bin/bash
###########################################################################
# Project: Orania Phylogeny MT
# Script: hashrf.sh
# --- Action: This script evaluates the simiarity between each gene tree
# ------------ and the species tree.
# --- Input: gene trees; species tree
# --- Output: matrix of RF values
# Author: Maya Schroedl
# Date: 11/2019
###########################################################################
#install newest version of https://code.google.com/archive/p/hashrf/downloads

WD=$PWD

#f: which files do you want to use ex: or,orscldyp,sclpod
#t: if you want trimmed contigs, based on coverage and threshold, enter t=threshold
while getopts ":f:t:" opt; do
  case ${opt} in
    f)
      echo "-f was triggered, File: $OPTARG"
      FILENAMES=${OPTARG}
      ;;
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

mkdir -p 6_hashrf

#make a genelist
cd $WD/3_gene_trees/$FILENAMES$t
rm "$FILENAMES$t"_RF_genelist.txt
echo astral > "$FILENAMES$t"_RF_genelist.txt
cat RAxML_parsimonyTree.* > "raxml_parsimony_"$FILENAMES$t"_gene.tre" && for i in RAxML_parsimonyTree.*; do ls $i | cut -d'_' -f 2| cut -d'.' -f 2 >> "$FILENAMES$t"_RF_genelist.txt ; done


cd $WD/6_hashrf

cat $WD/4_coalescent_trees/$FILENAMES$t/$FILENAMES$t"_coalescent_indiv_diff.tre" $WD/3_gene_trees/$FILENAMES$t/"raxml_parsimony_"$FILENAMES$t"_gene.tre" > $FILENAMES$t"_genetrees_astral_combi.tre"

hashrf $FILENAMES$t"_genetrees_astral_combi.tre" 0 -o $FILENAMES$t"_RF_matrix.txt"

cp $WD/3_gene_trees/"$FILENAMES$t"/"$FILENAMES$t"_RF_genelist.txt $WD/6_hashrf

#RF_sort.R ?

#make a new folder with genes selected by RF
mkdir -p $WD/2_alignment/$FILENAMES$t/RF_selected

while read name;
do
name=$(echo "$name" | tr -d '\r') #need to remove carriage return
cp $WD/2_alignment/$FILENAMES$t/$name"_aligned_gb_head.fasta" $WD/2_alignment/$FILENAMES$t/RF_selected;
done < $WD/6_hashrf/"$FILENAMES$t$t"_RF_selgenes.txt

####add empty sequences for taxa#####

#create species list
IFS=',' 
read -r -a FILENAMES$t$t_array <<< $FILENAMES$t$t

IFS=' '

rm $WD/specieslist.txt
for i in ${FILENAMES$t$t_array[@]};
do
cat $WD/../1_hybpiper/$i/namelist_$i".txt" >> $WD/specieslist.txt;
done

#run 7.1_biopython_RF_sp_add.py on windows (or modify script)


#concatenate alignments
cd $WD/2_alignment/$FILENAMES$t$t/RF_selected
cat *_sp_added.fasta|awk -v RS=">" -v FS="\n" -v OFS="\n" '{for(i=2; i<=NF; i++) {seq[$1] = seq[$1]$i}}; END {for(id in seq){print ">"id, seq[id]}}' > $FILENAMES$t$t"_RF_selected_concat.fasta"


#mv to $FILENAMES$t$t"_RF_selected_concat.fasta" windows 
#remove manually all numbers in astral tree (newick)
#do beauti file with astral monophyletic constraints (astral tree)
