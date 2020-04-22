#!/usr/bin/env Rscript
###########################################################################
# Project: Orania Phylogeny MT
# Script: collapse_low_brnlen.R
# --- Action: Collapses clades that present low branchlengths
# --- Input: directory to tree (newick format), output file directory
# --- Output: tree with collapsed clades (brnlen < 1 e-5)
# Author: Maya Schroedl
# Date: 04/2020
###########################################################################

rm(list = ls()) #clear environment

# Libraries ---------------------------------------------------------------
if (!require('ape')) install.packages('ape'); library('ape')

# Arguments ---------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)

tree = read.tree(args[1])
output = args[2]

collapsed_tree = di2multi(tree, 0.00001)

write.tree(collapsed_tree, output)

