# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 13:03:40 2019

@author: mayaa
"""

import os
import pandas
import numpy as np
import csv
from Bio import SeqIO
from pyfaidx import Fasta
import shutil

wd="D:/Programme/Google_Drive/Aarhus/Orania/Orania_project/2_Assembly_Alignment/Contigmaps_correction/"
dataset="or"
t=20

genelist=pandas.read_csv(wd+dataset+"_bwa_to_contigs/genelist_contigs.txt", names=["gene"])

if not os.path.exists(wd+dataset+"_trimmed_contigs/"+str(t)):
    os.makedirs(wd+dataset+"_trimmed_contigs/"+str(t))

if not os.path.exists(wd+dataset+"_trimmed_contigs/"+str(t)+"/or_contigs_exons"):
    shutil.copytree(wd+dataset+"_contigs_exons", wd+ dataset+"_trimmed_contigs/"+str(t)+"/or_contigs_exons")

for index,row in genelist.iterrows():
    gene=row[0]
    if os.path.exists(wd+dataset+"_bwa_to_contigs_genereorga/"+gene+"_contigs_sorted.coverage"):
        #if ... : because some genes are in the genelist, but they did not get mapped on the reffile and therefore do not have a contig
        cleanrecs=[]    
        mut_table=pandas.read_csv(wd+dataset+"_bwa_to_contigs_genereorga/"+gene+"_contigs_sorted.coverage", sep=" ",names=["sp","pos","cov"])
        records=Fasta(wd+dataset+"_trimmed_contigs/"+str(t)+"/"+dataset+"_contigs_exons/"+gene+".FNA", mutable=True)
        for record in records:
            mut_table.select=mut_table[mut_table["sp"]==record.name]
            for row in mut_table.select.itertuples(index=True, name='Pandas'):
                if int(row[3]) <= t:
                    record[int(row[2])-1]="N"
            
        records2=list(SeqIO.parse(wd+dataset+"_trimmed_contigs/"+str(t)+"/"+dataset+"_contigs_exons/"+gene+".FNA","fasta"))
        for record in records2:
            record.seq=record.seq.strip("N")
            cleanrecs.append(record)
        
        with open(wd+dataset+"_trimmed_contigs/"+str(t)+"/"+gene+"_"+"trimmed_"+str(t)+".fasta", 'w') as f:
            SeqIO.write(cleanrecs, f, 'fasta')     
