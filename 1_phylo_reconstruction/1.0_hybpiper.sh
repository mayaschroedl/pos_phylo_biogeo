#!/bin/bash
###########################################################################
# Project: Orania Phylogeny MT
# Script: hybpiper_script.sh
# --- Action: This script maps reads to the bait target file using HybPiper
# --- Input: trimmed reads & target file
# --- Output: assembled supercontigs, exons, and introns
# Author: Maya Schroedl
# Date: 11/2019
###########################################################################


#################################
#----DIRECTORIES & VARIABLES----#
#################################

#working directories
GWD=$PWD #global working directory, with subprojects and scripts
WD="$PWD"/1_phylo_reconstruction #current working directory

#TARGET FILE
TARGET_FILE=$WD/renamed_reads/target_file.fasta


#make a new directories
mkdir -p $WD/1.0_hybpiper/

mkdir -p $WD/1.0_hybpiper/mapped #mapped reads
mkdir -p $WD/1.0_hybpiper/contigs_exon #exon contigs
mkdir -p $WD/1.0_hybpiper/contigs_intron #intron contigs
mkdir -p $WD/1.0_hybpiper/supercontigs #supercontigs


##########################
#----FILE PREPARATION----#
##########################

cd $WD

####make genelist.txt
#create a genelist from target file, if genelist does not exist
if [ ! -f genelist.txt ]; then
    grep  '^>'  $TARGET_FILE >  headers.txt
    cut -f 2- -d '-' headers.txt > genelist.txt
    rm headers.txt
fi


####make namelist####
#namelist = list with all TAG names
if [ -s $WD/namelist.txt ] ; then rm $WD/namelist.txt ; fi #if file exists and is not empty, delete the file

cd $WD/0_trimmed/
for f in *"_R1_paired_trimmed.fastq"; do (echo ${f/"_R1_paired_trimmed.fastq"}>> $WD/namelist.txt); done #add all TAGs to namelist (all trimmed)

##################
#----HYBPIPER----#
##################

####execute hybpiper####
cd $WD/1.0_hybpiper #results will be stored here

#map reads to target file
while read name; #for each TAG, execute hybpiper with default settings
if [ -f $WD/0_trimmed/"$name"_unpaired_trimmed.fastq ]; then unpaired_flag=(--unpaired $WD/0_trimmed/"$name"_unpaired_trimmed.fastq); fi  #if unpaired trimmed file exists, then use the --unpaired flag in Hybpiper
do python2 $GWD/programs/HybPiper/reads_first.py -b $TARGET_FILE -r $WD/0_trimmed/"$name""_R"*"_paired_trimmed.fastq" --prefix "$name" --bwa "${unpaired_flag[@]}" ; 
done < $WD/namelist.txt

#get sequence lengths
python2 $GWD/programs/HybPiper/get_seq_lengths.py $TARGET_FILE $WD/namelist.txt dna > $WD/1.0_hybpiper/seq_lengths.txt

#get statistics
python2 $GWD/programs/HybPiper/hybpiper_stats.py $WD/1.0_hybpiper/seq_lengths.txt $WD/namelist.txt > $WD/1.0_hybpiper/stats.txt


#----EXONS----#

#get exon contigs
python2 $GWD/programs/HybPiper/retrieve_sequences.py $TARGET_FILE . dna


#move exon contigs to separate folder
mv *.FNA $WD/1.0_hybpiper/contigs_exon

#rename
cd $WD/1.0_hybpiper/contigs_exon
for f in *.FNA; do mv $f ${f/.FNA}_contig.fasta; done


#----INTRONS-AND-SUPERCONTIGS----#
#"introns" = non-coding sequences
#supercontigs = exons + "introns"

cd $WD/1.0_hybpiper/

#intronerate to detect introns
while read name; #for each TAG, intronerate and output to mapped/$name
do python2 $GWD/programs/HybPiper/intronerate.py --prefix $WD/1.0_hybpiper/$name;
done < $WD/namelist.txt

# #retrieve introns
python2 $GWD/programs/HybPiper/retrieve_sequences.py $TARGET_FILE . intron

#move intron contigs to separate folder
mv *_introns.fasta $WD/1.0_hybpiper/contigs_intron



#----SUPERCONTIGS----#
cd $WD/1.0_hybpiper/

#retrieve supercontigs
python2 $GWD/programs/HybPiper/retrieve_sequences.py $TARGET_FILE . supercontig

#move supercontigs to separate folder
mv *_supercontig.fasta $WD/1.0_hybpiper/supercontigs;

#-----------------------
#-----------------------

#---cleanup---#
mv TAG* $WD/1.0_hybpiper/mapped

cd $WD

