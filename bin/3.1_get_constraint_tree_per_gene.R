#!/usr/bin/env Rscript
###########################################################################
# Project: Orania phylo MT
# Script: get_constraint_tree_per_gene.R
# --- Action: 
# --- Input: 
# --- Output: 
# Author: Maya Schroedl
###########################################################################

rm(list=ls())


# Packages ----------------------------------------------------------------
if (!require('ape')) install.packages('ape'); library('ape')

if (!require('data.table')) install.packages('data.table'); library('data.table')



# Working directory -------------------------------------------------------
wd = file.path(getwd(),"1_phylo_reconstruction")
t = 2 #coverage trimming threshold



# Function ----------------------------------------------------------------
function(gene,
         wd = file.path(getwd(),"1_phylo_reconstruction"), 
         t = 2, #coverage trimming threshold
        ){
#### Input #### 
  # Constraint tree
  const.tree = read.tree(file.path(wd, "input_reads_and_info", "raxml_constraint.tree"))
  
  # Alignment
  in.fasta = file.path(wd, "2_alignment", t, paste0(gene,"_aligned_gb.fasta" ))
  
#### Output ####
output.dir = file.path(wd, "3_gene_trees", t, "0_constraint_trees")
if (!dir.exists(output.dir)){dir.create(output.dir)}
  
output.tree = file.path(output.dir, paste0(gene,"_constraint.tree"))
  
##### Get fasta headers (alignment "species" (individual) names) ####
al.names = fread(paste0("sed '/^>/s/;/\t/g;/^[^>]/d;s/^>//' ",in.fasta))



}
# Alignment


