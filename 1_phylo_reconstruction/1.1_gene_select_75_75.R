#!/usr/bin/env Rscript
###########################################################################
# Project: Orania Phylogeny MT
# Script: gene_select_75_75.R
# --- Action: How to filter loci with 75% of their length recovered in 75% of the TAGs
# --- Input: 
# --- Output: 
# Author: Maya Schroedl
# Date: 11/2019
###########################################################################

rm(list=ls())

# Libraries ---------------------------------------------------------------

if (!require('dplyr')) install.packages('dplyr'); library('dplyr')

####

# wd & arguments --------------------------------------------------------------
wd=file.path(getwd(),"1_phylo_reconstruction")

seq.ref.len=read.table(file.path(wd,"1.0_hybpiper","seq_lengths.txt")
                  ,header=T,sep="\t")

#first line is the reference length (length of the gene in the target file)
ref.len=seq.ref.len[1,]

#extract only recovered gene/sequence lengths (by HybPiper)
seq.len=seq.ref.len[-1,]

species=seq.len$Species #extract species (TAG) names

seq.len=seq.len[,!(names(seq.len)=="Species")]#remove species column
ref.len=ref.len[,!(names(ref.len)=="Species")]#remove species column

#calculate percentage of seq recovered for each seq in each TAG
percent.len = sweep(seq.len, 2, ref.len, "/")
percent.len.75=percent.len

#if percentage gene length recovered is >= 75% make value 1
#if not make value 0
percent.len.75[percent.len>=0.75] = 1
percent.len.75[percent.len<0.75] = 0


# Genes with 75% of species and 75% of their length reconstructed -------------

#calculate % TAGs with >= 75% of seq.len recovered for each gene
col = colSums(percent.len.75 != 0) / nrow(seq.len)

table(col>=0.75)


# Extract 75_75 genes -----------------------------------------------------

#Extract only 75_75 genes from original dataset
length7575<-seq.len[,col>=0.75]

#Get a list of names of 75_75 genes
yes75_75<-names(length7575[1,])

write.table(yes75_75,file.path(wd,"1.0_hybpiper",paste0("genelist_7575.txt")),quote=F,row.names=F,col.names=F)


#calculate percentage of genes having >75% reconstructed
r = rowSums(percent.len.75 != 0) / ncol(seq.len)
TAG_length7575=cbind(species,r)

TAG_length7575
