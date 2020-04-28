#!/usr/bin/env Rscript
###########################################################################
# Project: Orania phylo MT
# Script: root_tree.R
# --- Action: Roots a tree (with correct nodelables [support] and correct astral output). Because astral has a weird format
# --- Input: unrooted tree (astral or gene tree)
# --- Output: rooted tree in same format as input (astral) tree
# Author: Maya Schroedl
###########################################################################

rm(list=ls())

# Packages ---------------------------------------------------------------
if (!require('ape')) install.packages('ape'); library('ape')

# Input -------------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)

tree = args[1]
output_tree = args[2]
outgroup = args[3]

# Root tree ---------------------------------------------------------------
tr <- read.tree(tree)
tr <- root(tr, outgroup, edgelabel = TRUE,resolve.root=TRUE)
write.tree(tr,output_tree)

#Remove "Root" label
system(paste('sed -i "s/Root//g"', output_tree))

#When astral output with "q1, q2, & q3", the [q1, q2, q3] needs to be written as '[...]'
#
system(paste0('sed -i -e ', '"s/\\[/\'\\[/g" ', output_tree))
system(paste0('sed -i -e ', '"s/\\]/\\]\'/g" ', output_tree))
# and - needs to be replaced by ;
system(paste0('sed -i -e ', '"s/-q/;q/g" ', output_tree))

# Because this produces NaN as branchlengths for the tiplabels in astral trees, we want to remove the ":NaN"
system(paste('sed -i -e "s/:NaN//g"', output_tree))

