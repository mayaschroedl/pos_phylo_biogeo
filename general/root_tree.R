#!/usr/bin/env Rscript
###########################################################################
# Project: Orania phylo MT
# Script: root_tree.R
# --- Action: 
# --- Input: 
# --- Output:
# Author: Maya Schroedl
# Date: 11/2019
###########################################################################

rm(list=ls())

# Libraries ---------------------------------------------------------------
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

# Because this produces NaN as branchlengths for the tiplabels in astral trees, we want to remove the ":NaN"
system(paste('sed -i -e "s/:NaN//g"', output_tree))
