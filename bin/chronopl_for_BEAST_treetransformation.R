#https://osf.io/7y59t/wiki/home/

rm(list=ls())
library(ape)
library(geiger)
library(phytools)

setwd("D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/3_phylogeny_dating/GENOMEDK_Server/beast.or,orscldyp,sclpod_RF")

tree <- read.nexus("S886_aligned_gb_head_sp_added_cons.tree")
plot(tree, cex = 0.5)


#Get nodes of interest for dating
# syntax = mrca(phylogeny)["taxon1", "taxon2"]

mrca(tree)["TAG-32","TAG-44"]

######################################################################################
#chronopl explanation

#If age.max = NULL (the default), it is assumed that age.min gives 
#exactly known ages. Otherwise, age.max and age.min must be of the 
#same length and give the intervals for each node. Some node may be 
#known exactly while the others are known within some bounds: the 
#values will be identical in both arguments for the former 
#(e.g., age.min = c(10, 5), age.max = c(10, 6), node = c(15, 18) 
#means that the age of node 15 is 10 units of time, and the age of 
#node 18 is between 5 and 6).'

#ex.
#chronopl(phy, lambda, age.min = c(1,2), age.max = c(1,2),
#         node = c(1,2))

######################################################################################
#Calibration set 1
#["Phrynomantis_microps","Hemisus_marmoratus"] = 284; 90
#["Breviceps_adspersus","Hemisus_marmoratus"] = 286; 60
#["Trichobatrachus_robustus","Cardioglossa_gracilis"] = 318; 52
#["Kassina_decorata","Hyperolius_montanus"] = 394; 52

#examine difference between lambda values

#lambda = 0 
calib1a=chronopl(tree, lambda=1, age.min = 42.8, 
                    age.max = 42.8, node = 30)

plot(calib1a, cex = 0.3)

setwd("D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/3_phylogeny_dating/GENOMEDK_Server/beast.or,orscldyp,sclpod_RF/calib")

######################################################################################
write.tree(calib1a, file="S886_aligned_gb_head_sp_added_cons_chronopl_starting.tree")

