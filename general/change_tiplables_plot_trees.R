#!/usr/bin/env Rscript
###########################################################################
# Project: Orania Phylogeny MT
# Script: change_tiplables_plot_trees.R
# --- Action: Changes the tiplabels of a tree from TAGs to individuals (species)
# --- Input: tree with TAGS tiplables (newick) + dataframe indicating which
# ---------- TAG belongs to which individual
# --- Output: tree with individual tiplables (newick)
# Author: Maya Schroedl
# Date: 11/2019
###########################################################################

rm(list=ls())

# Libraries ---------------------------------------------------------------
if (!require('ape')) install.packages('ape'); library('ape')
if (!require('phylotools')) install.packages('phylotools'); library('phylotools')


# WD + Input --------------------------------------------------------------
wd=getwd()
t=2 # threshold
input_tree=file.path(wd, "4_coalescent_trees", t, "coalescent_lpp_rooted.tree")
tag_indiv_df=read.table(file.path(wd,"renamed_reads","tags_indiv.txt"),h=T) 

# Script -----------------------------------------------------------------
tree=read.tree(input_tree) #read the newick tree
lab_tree=sub.taxa.label(tree, tag_indiv_df) #replace tiplables

plot.phylo(lab_tree) #plot
plot.phylo(tree)

lab_tree$tip.label


# Output ------------------------------------------------------------------
write.tree(lab_tree, paste0(substr(input, 1, nchar(input)-6),"_lab.tree"))

