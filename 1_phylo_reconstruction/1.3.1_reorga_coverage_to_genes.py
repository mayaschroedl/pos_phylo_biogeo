#!/bin/bash
# -*- coding: utf-8 -*-
"""
###########################################################################
# Project: Orania Phylogeny MT
# Script: reorga_coverage_to_genes.py
# --- Action: This script reorganizes the coverage files. From TAGs to genes.
# --- Input: coverage maps per TAG (one file per TAG)
# --- Output: coverage maps per gene (one file per gene)
# Author: Maya Schroedl (maya.schroedl@bios.au.dk)

###########################################################################
"""

import pandas as pd
import os
import sys

#set working directory to current directory
WD=os.getcwd()

#read specieslist (TAGlist)
splist=pd.read_csv(WD+"/namelist.txt",names=["tag"]) 

for index,row in splist.iterrows(): #iterate over TAGs
    tag=row[0] #take TAG
	#open coverage map corresponding to this TAG
    cov_file=pd.read_csv(WD+"/1.1_coverage_maps/"+tag+"_supercontigs_sorted.coverage",sep="\t",names=["gene","pos","cov"]) 
    
    for gene_i in cov_file["gene"].unique(): #iterate over genes in coverage map
		# we have a TAGX-geneX selected
        select=cov_file[cov_file["gene"]==gene_i] #which rows in coverage map file correspond to current gene gene_i
        
        #replace column "gene" by "tags" for selected TAGX-geneX combination #IGNORE WARNING MESSAGE
        select["gene"]=[tag]*len(select["gene"]) 
        select2=select.to_csv(sep="\t",index=False,header=False)# convert covmap for TAGX-geneX to csv
		
		#write to new covmap file, organized per gene (append to new file)
        f=open(WD+"/1.1_coverage_maps/per_gene/"+gene_i+"_supercontigs_sorted.coverage","a+")
        f.write(select2)
        f.close()
		
#result: one file for each gene containing coverage maps for each TAG
        




