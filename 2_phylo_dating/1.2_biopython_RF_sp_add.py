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

WD="D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/2_hyb_align_reconstruction/AU_server/supercontigs/"

for filename in os.listdir(WD+"2_alignment/or,orscldyp,sclpod/RF_selected/"):
    with open(WD+"specieslist.txt","r") as totalsp:
        totalsp=totalsp.read().split('\n')
    
    records = list(SeqIO.parse(WD+"2_alignment/or,orscldyp,sclpod/RF_selected/"+filename, "fasta"))
    recordids=[]
    
    for r in records:
        recordids.append(r.id)
        
        # Python code to get difference of two lists
    diff=list(set(totalsp) - set(recordids))
    
    
    sequence=Seq("-"*len(records[2].seq),single_letter_alphabet)
    
    for sp_diff in diff:
        newrec=SeqRecord(sequence,id=sp_diff, description="")
        records.append(newrec)  
    
    SeqIO.write(records, WD+"2_alignment/or,orscldyp,sclpod/RF_selected/"+filename.strip(".fasta")+"_sp_added.fasta", "fasta")    
        
    
   # with open(WD+"2_alignment/or,orscldyp,sclpod/RF_selected/"+filename+"_diff","w") as difffile:
    #    for item in diff:
    #        difffile.write("%s\n" % item)

    