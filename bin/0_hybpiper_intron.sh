#!/bin/bash
###########################################################################
# Project: Orania Phylogeny MT
# Script: hybpiper_script.sh
# --- Action: This script maps reads to the bait target file using HybPiper
# --- Input: trimmed reads & target file
# --- Output: assembled contigs
# Author: Maya Schroedl
# Date: 01/2020
###########################################################################

#set working directory to current directory
WD=$PWD

# f: get input folder name -f = foldername in /0_trimmed
# t: get targetfile name
while getopts ":f:t:i" opt; do
  case ${opt} in
    f)
      echo "-f was triggered, File: $OPTARG"
      ID=${OPTARG}
      ;;
    t)
      echo "-t was triggered, Targetfile: $OPTARG"
      TARGET_FILE=${OPTARG}
      ;;
    i)
      echo "-i was triggered, Intron: $OPTARG"
      intron_bol=${OPTARG}
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
#HybPiper begins
#####################
cd $WD

####make genelist.txt
#create a genelist from target file, if genelist does not exist
if [ ! -f genelist.txt ]; then
    grep  '^>'  $TARGET_FILE >  headers.txt
    cut -f 2- -d '-' headers.txt > genelist.txt
    rm headers.txt
fi

#make a new directory 
mkdir -p $WD/1_hybpiper/
mkdir -p $WD/1_hybpiper/$ID #for this hybpiper analysis

####make namelist####
if [ -s $WD/1_hybpiper/$ID/namelist_$ID.txt ] ; then rm $WD/1_hybpiper/$ID/namelist_$ID.txt ; fi #if file exists and is not empty, delete the file

cd $WD/0_trimmed/$ID/
for f in *_R1_paired.fastq; do (echo ${f/_R1_paired.fastq}>> $WD/1_hybpiper/$ID/namelist_$ID.txt); done


####execute hybpiper####
cd $WD/1_hybpiper/$ID #results will be stored here

if ![ "$(ls -A $WD/1_hybpiper/$ID)" ];
then

#map reads to target file
while read name; 
do python2 $WD/programs/HybPiper/reads_first.py -b $WD/$TARGET_FILE -r $WD/0_trimmed/$ID/"$name"_R*_paired.fastq --prefix "$name" --bwa; 
done < $WD/1_hybpiper/$ID/namelist_$ID.txt


#get sequence lengths
python2 $WD/programs/HybPiper/get_seq_lengths.py $WD/target_file_Heyduk_baits_nuclear_exons_concatenated.fasta $WD/1_hybpiper/$ID/namelist_$ID.txt dna > seq_lengths_$ID.txt

#get statistics
python2 $WD/programs/HybPiper/hybpiper_stats.py seq_lengths_$ID.txt namelist_$ID.txt > stats_$ID.txt

#get exon contigs
python2 $WD/programs/HybPiper/retrieve_sequences.py $WD/target_file_Heyduk_baits_nuclear_exons_concatenated.fasta . dna


#cleanup
mkdir -p $WD/1_hybpiper/$ID/mapped_$ID
mv TAG* $WD/1_hybpiper/$ID/mapped_$ID

mkdir -p $WD/1_hybpiper/$ID/contigs_exon_$ID
mv *.FNA $WD/1_hybpiper/$ID/contigs_exon_$ID

#rename
cd $WD/1_hybpiper/$ID/contigs_exon_$ID

for f in *.FNA; do mv $f ${f/.FNA}_contig_$ID.FNA; done;
fi


### INTRONS ###
# generate introns & supercontigs

if intron_bol==TRUE;
then
mkdir -p $WD/introns
cd $WD/1_hybpiper/$ID/mapped_"$ID"

while read name; 
do python2 $WD/programs/HybPiper/intronerate.py --prefix $WD/1_hybpiper/$ID/mapped_"$ID"/"$name";
done < $WD/1_hybpiper/$ID/namelist_$ID.txt

cd $WD/1_hybpiper/$ID/mapped_"$ID"

python2 $WD/programs/HybPiper/retrieve_sequences.py $WD/target_file_Heyduk_baits_nuclear_exons_concatenated.fasta . intron

python2 $WD/programs/HybPiper/retrieve_sequences.py $WD/target_file_Heyduk_baits_nuclear_exons_concatenated.fasta . supercontig

mkdir -p $WD/1_hybpiper/$ID/contigs_intron_$ID
mv *_introns.fasta $WD/1_hybpiper/$ID/contigs_intron_$ID

mkdir -p $WD/1_hybpiper/$ID/supercontigs_$ID
mv *_supercontig.fasta $WD/1_hybpiper/$ID/supercontigs_$ID;
fi



