#!/usr/bin/env Rscript
###########################################################################
# Project: Orania Phylogeny MT
# Script: change_tiplabels.R
# --- Action: 
# --- Input: 
# --- Output: 
# Author: Maya Schroedl
# Date: 04/2020
###########################################################################

rm(list=ls())

# Packages ---------------------------------------------------------------
if (!require('ape')) install.packages('ape'); library('ape')
if (!require('phylotools')) install.packages('phylotools'); library('phylotools')

# Arguments ---------------------------------------------------------------
args = commandArgs(trailingOnly=TRUE)

input_tree = args[1]
output_tree = args[2]
tag_indiv_df = read.table(args[3],sep="\t",h=T)

# Function -----------------------------------------------------------------
change_tip_labs=function(input_tree,output_tree,tag_indiv_df){
  tree=read.tree(input_tree) #read the newick tree
  lab_tree=sub.taxa.label(tree, tag_indiv_df) #replace tiplables
  #output
  write.tree(lab_tree, output_tree)
  
  #When astral output with "q1, q2, & q3", the [q1, q2, q3] needs to be written as '[...]' 
  system(paste0('sed -i -e ', '"s/\\[/\'\\[/g" ', output_tree))
  system(paste0('sed -i -e ', '"s/\\]/\\]\'/g" ', output_tree))
  # and - needs to be replaced by ;
  system(paste0('sed -i -e ', '"s/-q/;q/g" ', output_tree))
  
  # Because this produces NaN as branchlengths for the tiplabels in astral trees, we want to remove the ":NaN"
  system(paste('sed -i -e "s/:NaN//g"', output_tree))
}

change_tip_labs(input_tree,output_tree,tag_indiv_df)


