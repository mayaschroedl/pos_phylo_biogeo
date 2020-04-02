library(ape)
library(phyclust)

rm(list=ls())

setwd("D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/3_phylogeny_dating/GENOMEDK_Server/beast.or,orscldyp,sclpod_RF/uncalib")

tree <- read.nexus("S886_aligned_gb_head_sp_added_cons.tree")
plot(tree)
add.scale.bar()

scale=67/get.rooted.tree.height(tree)

tree$edge.length=tree$edge.length*scale

plot(tree)
edgelabels(round(tree$edge.length,3),bg="white")
add.scale.bar()


setwd("D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/3_phylogeny_dating/GENOMEDK_Server/beast.or,orscldyp,sclpod_RF/calib")

######################################################################################
write.tree(tree, file="S886_aligned_gb_head_sp_added_R_calib.tree")

