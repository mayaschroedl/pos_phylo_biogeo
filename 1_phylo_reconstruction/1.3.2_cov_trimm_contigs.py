#!/bin/bash
# -*- coding: utf-8 -*-
"""
###########################################################################
# Project: Orania Phylogeny MT
# Script: trimm_contigs.py
# --- Action: This script replaces the bases that have a coverage under 
# ----------- a certain threshold t with an N and then trimms the ends
# ----------- of the gene.
# --- Input:  contigs and coverage maps for each gene
# --- Output: trimmed contigs
# Author: Maya Schroedl (maya.schroedl@bios.au.dk)

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

genelist=pandas.read_csv(WD+"/genelist_7575.txt", names=["gene"])#read genelist containing the filtered genes with the script

for index,row in genelist.iterrows(): #iterate over the genelist
    gene=row[0]; print(gene); print("Progress: "+str(index+1)+" gene(s) out of "+ str(len(genelist)))
    
    if os.path.exists(WD+"/1.1_coverage_maps/per_gene/"+gene+"_supercontigs_sorted.coverage"): 
        # it is possible that some genes are in the genelist, but they did not get 
        # mapped on the reference file and therefore do not have a contig. 
        # should have been filtered with script gene_selection_75:75.R; but just in case.
        
        cleanrecs=[]    # here we are going to store the trimmed sequences (records)
        covmap_gene=pandas.read_csv(WD+"/1.1_coverage_maps/per_gene/"+gene+
            "_supercontigs_sorted.coverage", sep="\t",names=["sp","pos","cov"]) #coverage map corresponding to this gene
        records=Fasta(WD+"/1.1_coverage_maps/trimmed_contigs/"+str(t)+
            "/contigs/"+gene+"_supercontig.fasta", mutable=True) #open and read the contigs for this gene
        
        for record in records: #for each species record
            #select the coverage map corresponding to this species (and gene)
            covmap_gene.select=covmap_gene[covmap_gene["sp"]==record.name] 
            
            for row in covmap_gene.select.itertuples(index=True, name='Pandas'):
                #for each row (each position) in coverage map, 
                #if the depth of this position is under the threshold
                if (int(row[3]) < int(t)): 
                    #replace the position [row(2)] in the record of this contig with "N" 
                    #(is automatically overwritten in the original contig file)
                    record[int(row[2])-1]="N"; 
					
                  
        masked_records=list(SeqIO.parse(WD+"/1.1_coverage_maps/trimmed_contigs/"+str(t)+
            "/contigs/"+gene+"_supercontig.fasta","fasta")) #reread the file which was just modified and masked with "N"s
       
        for record in masked_records: # for each contig
            #trimm the ends NNNNNNN
            record.seq=record.seq.strip("N") 
            cleanrecs.append(record) # put the cleaned record to cleanrecs
        
        #write the trimmed records for this gene (cleanrecs) into a new file
        with open(WD+"/1.1_coverage_maps/trimmed_contigs/"+str(t)+"/"+gene+
            "_supercontig.fasta", 'w') as f: 
            SeqIO.write(cleanrecs, f, 'fasta')
