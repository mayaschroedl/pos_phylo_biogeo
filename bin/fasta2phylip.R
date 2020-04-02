install.packages("seqinr")
install.packages("ape")

library(seqinr)
library(ape)

setwd("D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/3_phylogeny_dating/GENOMEDK_Server/paml/or,orscldyp,sclpod_align_supercontig_75/")
fastaobject<-seqinr::read.fasta("E1007_aligned_gb_head.fasta")
ape::write.dna(fastaobject, "E1007_aligned_gb_head.phy", nbcol=1,colsep="", colw=1000000)
