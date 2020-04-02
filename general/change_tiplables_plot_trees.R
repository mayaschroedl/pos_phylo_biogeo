#install.packages("ape")
#install.packages("phylotools")
library(ape)
library(phylotools)

input="D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/3_phylogeny_dating/GENOMEDK_Server/beast.or,orscldyp,sclpod_RF/calib/S886_aligned_gb_head_sp_added_R_calib.tree"
tree=read.tree(input) #change to read.nexus
list=read.table("D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/sp_TAGs.txt")

sub=sub.taxa.label(tree, list)

plot.phylo(sub)
plot.phylo(tree)

sub$tip.label


write.tree(sub, paste0(substr(input, 1, nchar(input)-5),"_lab.tre"))

