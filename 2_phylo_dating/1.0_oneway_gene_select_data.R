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

# Minimal Example ------------------------------------------------------------

### Input ----
sptree = read.tree(text = "((((A,B),C),(D,E)),((F,G),H));")
plot.phylo(sptree,main="sptree")

genetree = read.tree(text = "(((A,B,C,D)95,E)40,((F,G)10,H)80);")
plot.phylo(genetree,main="genetree")
nodelabels(text = genetree$node.label,frame="none")
gene_name = "example"

#### Output ----
output = file.path(wd, "2_phylo_dating","1_oneway_gene_select","minimal_example.stats")

#### Get distances ----
file.remove(output) #delete output file
#signal_support_stats(genetree, gene_name, sptree, output)

example_stats=read.table(output, h=T)
View(example_stats)

# Actual application ------------------------------------------------------

### Input ---
## Rooted gene trees ----
genetrees_folder = file.path(wd, "1_phylo_reconstruction", "3_gene_trees", "2", "4_collapsed") #where are the gene trees
gene_list = read.table(file.path(wd, "1_phylo_reconstruction", "genelist_7575.txt"))[,1] #list of all the genes
suffix = ".raxml.support.coll_rooted" #how are the gene tree files named
 
## Rooted species tree ----
sptree_file = file.path(wd, "1_phylo_reconstruction","4_coalescent_trees", "2","coalescent_lpp.tree_rooted" ) #where is the species tree

sptree = read.tree(sptree_file) #read species tree

### Output ----
output = file.path(wd, "2_phylo_dating","1_oneway_gene_select","my_genetrees.stats")

#### Get distances ----
if (file.exists(output)){file.remove(output)} #delete output file
for (gene in gene_list){
  genetree = read.tree(file.path(genetrees_folder,paste0(gene, suffix)))
  #genetree$edge.length = NULL #make relations better visible (remove branch lengths for this analysis)
  #plot(genetree, main = gene) #have a look how tree looks
  #nodelabels(text=genetree$node.label,frame = "none") #add bootstrap support
  
  signal_support_stats(genetree,gene, sptree, output)
  
  }

#### Select genes ----

stats=read.table(output, h=T)

View(stats)

# we would like to chose the genes that have:
#   - the most good nodes in general (to not only take "bad trees")
#   - the least good nodes disagreeing with sptree
#   - the most good nodes agreeing with sptree
#   - the most clocklike


stats_sorted = stats %>%
  arrange(gnd_disagree_perc,desc(gnd_agree_perc),dist_root)

stats_sorted = stats %>%
  arrange(gnd_disagree_perc,dist_root)

#####
# 1) Genes with 0 disagreeing good nodes
no_gdis = stats$genetree[which(stats$gnd_disagree_perc == 0 & stats$gnd_agree_perc != 0)]

# 2) no_gdis + genes with most difference between disagreeing good nodes & agreeing good nodes

#get the genes which aren't included in 1)
stats_nogdis = stats[-which(stats$gnd_disagree_perc == 0 & stats$gnd_agree_perc != 0),]

# sort stats by difference between gnd_agree_perc and gnd_disagree_perc
stats_nogdis_sorted = stats_nogdis %>%
 arrange(desc(stats_nogdis$diff_ga_gd))

most_diff = stats_nogdis_sorted$genetree[1:10]

# 3) possible clock-like genes (having 0 disagreeing good nodes). need to be tested with Bayes factor beast

#sort stats by ascending root distance
stats_dist = stats %>%
  arrange(gnd_disagree_perc,dist_root)

one_pclock = stats_dist$genetree[1]
three_pclock = stats_dist$genetree[1:3]
nine_pclock = stats_dist$genetree[1:9]
