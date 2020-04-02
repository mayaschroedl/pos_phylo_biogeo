#!/bin/bash
# -*- coding: utf-8 -*-
"""
###########################################################################
# Project: Orania Phylogeny MT
# Script: trimm_contigs.py
# --- Action: This script replaces the bases that have a coverage under 
# ----------- a certain threshold “t” with an “N” and then trimms the ends
# ----------- of the gene.
# --- Input:  contigs and coverage maps for each gene
# --- Output: trimmed contigs
# Author: Maya Schroedl
# Date: 11/2019
###########################################################################
"""

import os
import sys
import pandas #to manipulate dataframes
from Bio import SeqIO #to manipulate fastas
from pyfaidx import Fasta #to overwrite fastas
import shutil #to easily copy files and whole folders

#set working directory to current directory
WD=os.getcwd()

t=sys.argv[1] #threshold

genelist=pandas.read_csv(WD+"/1_hybpiper/genelist_7575.txt", names=["gene"])#read genelist containing the filtered genes with the script

if not os.path.exists(WD+"/1.1_coverage_maps/trimmed_contigs/"+str(t)+"/contigs"):
    shutil.copytree(WD+"/1_hybpiper/supercontigs" #contigs for each gene
                    , WD+"/1.1_coverage_maps/trimmed_contigs/"+str(t)+"/contigs")

for index,row in genelist.iterrows():
    gene=row[0]
    if os.path.exists(WD+"/1.1_coverage_maps/"+gene+"_supercontigs_sorted.coverage"):
        #if ... : because some genes are in the genelist, but they did not get mapped on the reffile and therefore do not have a contig
        cleanrecs=[]    
        mut_table=pandas.read_csv(WD+"/1.1_coverage_maps/"+gene+"_supercontigs_sorted.coverage", sep="\t",names=["sp","pos","cov"])
        records=Fasta(WD+"/1.1_coverage_maps/trimmed_contigs/"+str(t)+"/contigs/"+gene+"_supercontig".fasta", mutable=True)
        for record in records:
            mut_table.select=mut_table[mut_table["sp"]==record.name]
            for row in mut_table.select.itertuples(index=True, name='Pandas'):
                if (int(row[3]) < int(t)):
                    record[int(row[2])-1]="N"
                
            
        records2=list(SeqIO.parse(WD+"/1.1_coverage_maps/trimmed_contigs/"+str(t)+"/contigs/"+gene+"_supercontig.fasta","fasta"))
        for record in records2:
            record.seq=record.seq.strip("N")
            cleanrecs.append(record)
        
        with open(WD+"/1.1_coverage_maps/trimmed_contigs/"+str(t)+"/"+gene+"_covtrimmed.fasta", 'w') as f:
            SeqIO.write(cleanrecs, f, 'fasta')     
