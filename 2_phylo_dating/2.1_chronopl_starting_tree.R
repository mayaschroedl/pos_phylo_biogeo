#!/usr/bin/env Rscript
###########################################################################
# Project: Orania phylo MT
# Script: root_tree.R
# --- Action: Makes starting tree for BEAST.
# --- Input: rooted astral tree
# --- Output: rooted ultrametric starting tree (calibrated)
# Author: https://osf.io/7y59t/wiki/home/, adapted by Maya Schrödl
###########################################################################
rm(list=ls())

# Packages ---------------------------------------------------------------

if (!require('ape')) install.packages('ape'); library('ape')
if (!require('geiger')) install.packages('geiger'); library('geiger')
if (!require('phytools')) install.packages('phytools'); library('phytools')


# Working directory -------------------------------------------------------
wd = getwd()

# Input & output -------------------------------------------------------------
input = file.path("1_phylo_reconstruction","4_coalescent_trees","2","coalescent_lpp.tree_rooted")

output = file.path(wd, "2_phylo_dating", "3_babette","coalescent_lpp_chronopl.tree_rooted")
                   
# Astral tree -------------------------------------------------------------
astral_tree <- read.tree(input)
plot(astral_tree, cex = 0.5)

#Get nodes of interest for dating
calib_node = mrca(astral_tree)["TAG-32","TAG-47"] #we will calibrate the root

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


# because astral outputs no branchlengths for tips (NaN), we replace the NaNs with a very low branchlength to make the tree calibratable with chronopl
astral_tree$edge.length[which(astral_tree$edge.length=="NaN")] = rep(0.000001, length(astral_tree$edge.length[which(astral_tree$edge.length=="NaN")]))


# create starting tree
calib_starting= chronopl(astral_tree, lambda=1, age.min = 57, 
                    age.max = 57, node = calib_node) #we will later calibrate the root at the age "57" +- 95% HP (taken from BakerCouvreur2013a)

# check if worked
plot(calib_starting, cex = 0.3)
add.scale.bar()

# root starting tree
calib_starting_r <- root(calib_starting, "TAG-47", edgelabel = TRUE,resolve.root=TRUE)


######################################################################################
write.tree(calib_starting_r, file = output)

#Remove "Root" label
system(paste('sed -i "s/Root//g"', output))
