# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 12:34:22 2019

@author: mayaa
"""
import pandas as pd

wd="D:/Programme/Google_Drive/Aarhus/Orania/Orania_project/2_Assembly_Alignment/Contigmaps_correction/"

## reorganize data to genes
splist=pd.read_csv(wd+"or_bwa_to_contigs/splist_thomas.txt",names=["tag"])

for index,row in splist.iterrows():
    tag=row[0]
    TAG=pd.read_csv(wd+"or_bwa_to_contigs/"+tag+"_contigs_sorted.coverage",sep="\t",names=["gene","pos","cov"])
    for ex in TAG["gene"].unique():
        select=TAG[TAG["gene"]==ex]
        select["gene"]=[tag]*len(select["gene"])
        select2=select.to_csv(sep="\t",index=False,header=False)
        f=open(wd+"or_bwa_to_contigs_genereorga2/"+ex+"_contigs_sorted.coverage","a+")
        f.write(select2)
        f.close()
         




