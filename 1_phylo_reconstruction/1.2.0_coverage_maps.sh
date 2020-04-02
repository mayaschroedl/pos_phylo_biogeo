#!/bin/bash
###########################################################################
# Project: Orania Phylogeny MT
# Script: coverage_maps.sh
# --- Action: This script computes a "coverage map" for each supercontig. 
# ----------- A coverage map represents how many reads support each
# ----------- position over the supercontig. Maps reads to supercontigs.
# --- Input: supercontigs from script "hybpiper_script.sh"
# --- Output: coverage map for each supercontig (per gene & species)
# Author: Maya Schroedl
# Date:11/2019
###########################################################################

#####################
#----DIRECTORIES----#
#####################

#working directories
GWD=$PWD #global working directory, with subprojects and scripts
WD="$PWD"/1_phylo_reconstruction #current working directory

mkdir -p $WD/1.0_hybpiper/supercontigs_reorga #here we're going to store the supercontigs which are reorganized into TAGs
mkdir -p $WD/1.1_coverage_maps/ #coverage maps are going to be stored here

##########################
#----FILE PREPARATION----#
##########################

#----PROBLEM----#
#Hybpiper gives out supercontigs for each gene in one file (fasta headers: TAGX-geneX). 
#but bwa wants one file per TAG with supercontigs for each gene (fasta headers: gene)
#----SOLUTION---#
#--> we have to reorganize and transform our data

cd $WD/1.0_hybpiper/supercontigs_reorga

rm *

while read gene; #for each gene
    #split each multifastafile into one single fasta file [TAGX,geneX]
	do cat $WD/1.0_hybpiper/supercontigs/"$gene"_supercontig.fasta  | awk '{if (substr($0, 1, 1)==">") {ID=(substr($0,2) ".fasta")} print $0 > ID}'
	#rename header of each new splitted file (*$gene.fasta) from TAGX-geneX to geneX
	for file in *$gene.fasta; do sed -i "s/>.*/>$gene/" $file; done #
done < $WD/genelist_7575.txt

#combine all files for each TAG
while read name; #for each TAG: combine all gene supercontig files in one file
	do cat $WD/1.0_hybpiper/supercontigs_reorga/"$name"*.fasta > "$name"_supercontigs.fasta;
done < $WD/namelist.txt #namelist = list of all TAGs

#-cleanup-#
cd $WD/1.0_hybpiper/supercontigs_reorga
rm !("TAG"*"_supercontigs.fasta") #remove all files except the combined supercontigs files

###################################
#----MAP-READS-TO-SUPERCONTIGS----#
###################################
# Like in https://github.com/sidonieB/bioinfo-utils/blob/master/docs/advice/target_capture_data_analysis.md; 6. Assess target capture efficiency; Coverage of the recovered regions

cd $WD/1.0_hybpiper/supercontigs_reorga

#index supercontigs
for f in *.fasta; do bwa index $f; done

#map reads (.fastq) with bwa to supercontigs (*_supercontigs.fasta)
#map paired reads and unpaired reads (if there are) independently. we will merge them later
while read name; 
	do bwa mem -t 32 $WD/1.0_hybpiper/supercontigs_reorga/"$name"_supercontigs.fasta $WD/0_trimmed/"$name"_R1_paired_trimmed.fastq $WD/0_trimmed/"$name"_R2_paired_trimmed.fastq > $WD/1.1_coverage_maps/"$name"_supercontigs_BWA.sam; 
	if [ -f $WD/0_trimmed/"$name"_unpaired_trimmed.fastq ]; then bwa mem -t 32 $WD/1.0_hybpiper/supercontigs_reorga/"$name"_supercontigs.fasta $WD/0_trimmed/"$name"_unpaired_trimmed.fastq > $WD/1.1_coverage_maps/"$name"_unpaired_supercontigs_BWA.sam; fi
done < $WD/namelist.txt

#convert sam to bam
cd $WD/1.1_coverage_maps/
for f in *BWA.sam; do (samtools view -b $f -o ${f/.sam}.bam); done

#merge upaired and paired .bam files
while read name; 
do if [ -f $WD/1.1_coverage_maps/"$name"_unpaired_supercontigs_BWA.bam ]; 
then samtools merge $WD/1.1_coverage_maps/"$name"_supercontigs_BWA.sam $WD/1.1_coverage_maps/"$name"_supercontigs_BWA.sam $WD/1.1_coverage_maps/"$name"_unpaired_supercontigs_BWA.sam;
fi;
done < $WD/namelist.txt
# samtools merge output input1(paired) input2(unpaired)

#sort by position
for f in *BWA.bam; do (samtools sort $f -o ${f/.bam}_sorted.bam); done

#index bam files #index positions
for f in *_sorted.bam; do samtools index $f; done

#calculate coverage for each position (how many reads per base pair)
for f in *_sorted.bam; do samtools depth $f > ${f/_BWA_sorted.bam}"_sorted.coverage"; done;

#remove unnecessary files
rm *.sam
rm *.bam

###############################
#----VISUALIZATION + STATS----#
###############################
# check script "visualize_covmaps_stats.R

#-----------------------
#-----------------------

cd $WD
