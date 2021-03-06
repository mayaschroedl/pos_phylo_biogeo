# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 15:43:30 2020

@author: mayaa
"""
#script to put empty sequence for species that have no seq (after alignment)
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
import glob

t=str(2)

WD="D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project"
if not os.path.exists(WD+ "/2_phylo_dating/1_alignment"):
	    os.makedirs(WD+ "/2_phylo_dating/1_alignment")
	    
phylo_dir = WD+"/1_phylo_reconstruction/2_alignment/"+t
os.setdir(phylo_dir)

matching_filenames=glob.glob("*_mod.fasta")

for filename in matching_filenames:
    with open(WD+"/1_phylo_reconstruction/namelist.txt","r") as totalsp:
        totalsp=totalsp.read().split('\n')
    
    records = list(SeqIO.parse(WD+"/1_phylo_reconstruction/2_alignment/"+t+"/"+filename, "fasta"))
    recordids=[]
    
    for r in records:
        recordids.append(r.id)
        
        # Python code to get difference of two lists
    diff=list(set(totalsp) - set(recordids))
    
    sequence=Seq("-"*len(records[2].seq),SingleLetterAlphabet) 
    
    for sp_diff in diff: # add sequence of "-" for species which are not represented in this gene alignment
        newrec=SeqRecord(sequence,id=sp_diff, description="")
        records.append(newrec)  
    
	
	
    SeqIO.write(records, WD+ "/2_phylo_dating/1_alignment/"+filename, "fasta")    
        
