###########################################################################
# Project: Orania Phylogeny MT
# Script: alignment_stats.R
# --- Action: get statistics for cleaned alignments
# --- Input: 
# --- Output: 
# Author: Maya Schroedl (maya.schroedl@bios.au.dk)
###########################################################################


rm(list=ls())

# Packages ----------------------------------------------------------------
if (!require('seqinr')) install.packages('seqinr'); library('seqinr')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('rlist')) install.packages('rlist'); library('rlist')
if (!require('ape')) install.packages('ape'); library('ape')

#not in
'%!in%' <- function(x,y)!('%in%'(x,y))

# Working directory -------------------------------------------------------
gwd = getwd()
t = 2 #coverage trimming threshold
wd = file.path(getwd(),"1_phylo_reconstruction","2_alignment",t)

# Alignment length --------------------------------------------------------
al.len = function(gene){ #gene: gene name
  al = read.fasta(file.path(wd, paste0(gene, "_aligned_gb_mod.fasta"))) # alignment
  len = length(al[1][[1]])
  return(len)
}

#genelist
genelist = read.table(file.path(gwd, "1_phylo_reconstruction","genelist_7575.txt"))[,1]

al_len_all = sapply(genelist,al.len) # list of alignement length for all genes

#mean alignment length over all genes
mean(al_len_all)

#standard deviation
sd(al_len_all)

#range of alignment length over all genes
min(al_len_all)
max(al_len_all)


# Percentage of missing data ----------------------------------------------
# Get concatenated alignments
all_alignments_concat = read.fasta(file.path(gwd,"2_phylo_dating", "1_alignment", "concat", "all_genes_concat.fasta")) #get concatenated alignments for all individuals

unlisted = unlist(all_alignments_concat) # put all alignements into one row

#this is the percentage of missing data:
perc_mis_dat = (length(which(unlisted == "-")) + length(which(unlisted == "n")))/length(unlisted)
  
  
