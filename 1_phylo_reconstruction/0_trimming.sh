#!/bin/bash
###########################################################################
# Project: Orania Phylogeny MT
# Script: trimming.sh
# --- Action: This script trimms the reads
# --- Input: untrimmed reads
# --- Output: trimmed reads and fastqc files
# Author: Maya Schroedl
# Date: 03/2019
###########################################################################
  
#####################
#----DIRECTORIES----#
#####################

#working directories
GWD=$PWD #global working directory, with subprojects and scripts
WD="$PWD"/1_phylo_reconstruction #current working directory

#make a new directory 
mkdir -p  $WD/0_trimmed/fastqc

##########################
#----FILE PREPARATION----#
##########################

#copy trimmed reads in renamed_reads to 0_trimmed/
cp $WD/renamed_reads/*trimmed.fastq $WD/0_trimmed/


####make namelist for untrimmed reads####
#namelist = list with all TAG names
if [ -s $WD/renamed_reads/namelist_untrimmed.txt ] ; then rm $WD/renamed_reads/namelist_untrimmed.txt ; fi #if file exists and is not empty, delete the file

cd $WD/renamed_reads/
for f in *"_R1_paired.fastq"; do (echo ${f/"_R1_paired.fastq"}>> $WD/renamed_reads/namelist_untrimmed.txt); done #add TAGs to namelist



################################
#----TRIMMOMATIC AND FASTQC----#
################################

#----Trimmomatic----#
# with default settings for adapter trimming (Illumina Adapters) and filtering of reads
while read f; 
do 
	##TRIMMOMATIC #default: java -jar trimmomatic-0.39.jar PE input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
	 java -jar $GWD/programs/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 $WD/renamed_reads/"$f"_R1_paired.fastq $WD/renamed_reads/"$f"_R2_paired.fastq $WD/0_trimmed/"$f"_R1_paired_trimmed.fastq $WD/0_trimmed/"$f"_R1_unpaired.fastq $WD/0_trimmed/"$f"_R2_paired_trimmed.fastq $WD/0_trimmed/"$f"_R2_unpaired.fastq ILLUMINACLIP:$GWD/programs/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:1:true LEADING:30 TRAILING:30 SLIDINGWINDOW:4:30 MINLEN:36 AVGQUAL:30 2>$WD/0_trimmed/"$f"_trim.log
	cat $WD/0_trimmed/"$f"_R*_unpaired.fastq > $WD/0_trimmed/"$f"_unpaired_trimmed.fastq #put unpaired reads to one file (HybPiper needs only one unpaired file)
	rm $WD/0_trimmed/"$f"_R*_unpaired.fastq #remove unpaired files seperated into reads
 done < $WD/renamed_reads/namelist_untrimmed.txt


# #----FASTQC----#
# # for all trimmed reads for visualization of the read qualities

# cd $WD/0_trimmed/

# for f in *.fastq; #for all trimmed reads
# do fastqc $f --outdir  $WD/0_trimmed/fastqc; 
# done

# rm $WD/0_trimmed/fastqc/*.zip #remove all unnecessary zip files

# #### CHECK FASTQC OUTPUT ####

#-----------------------
#-----------------------

cd $WD
