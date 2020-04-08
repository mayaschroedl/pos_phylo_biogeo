#!/usr/bin/env Rscript
###########################################################################
# Project: Orania Phylogeny MT
# Script: oneway_gene_select.R
# --- Action: Compares gene trees to Astral tree (one direction) in order to 
# ------------ select the genes that are the most similar to the Astral tree
# ------------ and the most informative (BS > 70%)
# --- Input: species tree, folder with individual gene trees
# --- Output: 
# Author: Maya Schroedl
# Date: 11/2019
###########################################################################
###########################################################################
###########################################################################

rm(list=ls())

# Functions ------------------------------------------------------------------
source(file.path(getwd(), "scripts", "2_phylo_dating","1.0_oneway_gene_select.R"))


# Libraries ---------------------------------------------------------------
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')


# Working directories -----------------------------------------------------
wd = getwd()

# Easy Example ------------------------------------------------------------

### Input ----
sptree = read.tree(text = "(((A,B),(C,D)),((E,F),G));")
plot.phylo(sptree)

genetree = read.tree(text = "(((A,B,C,F)0.95,D)0.4,((E,F)0.1,G)0.8)0.4;")
plot.phylo(genetree)
nodelabels(frame="none")
gene_name = "example"

#### Output ----
output = file.path(wd, "2_phylo_dating","1_oneway_gene_select","easy_example.stats")


#### Get distances ----
gtree_sptree_distance(genetree, gene_name, sptree, output)

# Actual application ------------------------------------------------------

### Input ---
## Rooted gene trees ----
genetrees_folder = file.path(wd, "1_phylo_reconstruction", "3_gene_trees", "2", "4_collapsed") #where are the gene trees
gene_list = read.table(file.path(wd, "1_phylo_reconstruction", "genelist_7575.txt"))[,1] #list of all the genes
suffix = "_rooted.raxml.support.coll" #how are the gene tree files named
 
## Rooted species tree ----
sptree_file = file.path(wd, "1_phylo_reconstruction","4_coalescent_trees", "2","coalescent_lpp_rooted.tree" ) #where is the species tree

sptree = read.tree(sptree_file) #read species tree

### Output ----
output = file.path(wd, "2_phylo_dating","1_oneway_gene_select","my_genetrees.stats")

#### Get distances ----
file.remove(output) #delete output file
for (gene in gene_list){
  genetree = read.tree(file.path(genetrees_folder,paste0(gene, suffix)))
  gtree_sptree_distance(genetree,gene, sptree, output)
  }

#### Select genes ----

stats=read.table(output, h=T)

stats_sorted = stats %>% 
  arrange(desc(gnd_agree_perc), gnd_disagree_perc)

View(stats_sorted)
