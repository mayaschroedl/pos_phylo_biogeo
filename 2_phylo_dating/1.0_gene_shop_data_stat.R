#!/usr/bin/env Rscript
###########################################################################
# Project: Orania Phylogeny MT
# Script: gene_shop_data_stat.R
# --- Action: Compares gene trees to Astral tree (one direction) in order to 
# ------------ select the genes that have the least "good" (BS > 75%) nodes disagreeing with the species tree. Also select among these selected genes the most clock-like genes.
# --- Input: species tree, folder with individual rooted gene trees (with outgroup! + collapsed nodes (BS < 10) and collapsed nodes (branch length < 0.00002))
# --- Output: sorted statistic based on the signal_support_stats function of the script gene_shop_fct.R and selected genes for different scenarii
# Author: Maya Schroedl
###########################################################################

rm(list=ls())

# Functions ------------------------------------------------------------------
source(file.path(getwd(), "scripts", "2_phylo_dating","1.0_gene_shop_fct.R"))


# Libraries ---------------------------------------------------------------
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')


# Working directories -----------------------------------------------------
wd = getwd()

# Application ------------------------------------------------------

### Input ---
## Rooted gene trees ----
genetrees_folder = file.path(wd, "1_phylo_reconstruction", "3_gene_trees", "2", "4_collapsed") #where are the gene trees
gene_list = read.table(file.path(wd, "1_phylo_reconstruction", "genelist_7575.txt"))[,1] #list of all the genes
suffix = ".raxml.support.coll_rooted" #how are the gene tree files named
 
## Rooted species tree ----
sptree_file = file.path(wd, "1_phylo_reconstruction","4_coalescent_trees", "2","coalescent_lpp.tree_rooted" ) #where is the species tree

sptree = read.tree(sptree_file) #read species tree

### Output ----
output = file.path(wd, "2_phylo_dating","1_gene_shop","my_genetrees.stats")

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

stats_sorted = stats %>%
  arrange(gnd_disagree_perc, desc(gnd_agree_perc))

# we would like to chose the genes that have:
#   - the least good nodes disagreeing with sptree
#   - the most clocklike

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

#combine no_gdis + most_diff
no_gdis_most_diff = c(no_gdis, most_diff)


# 3) possible clock-like genes (having 0 disagreeing good nodes). These genes need to be tested with Bayes factor beast, to see whether they are really clock-like and can be run with a strict model or not.

#sort stats by ascending root distance
stats_clock = stats[which(stats$genetree %in% no_gdis),] %>% #only have a look at no_gdis genes
  arrange(dist_root)

one_clock = stats_clock$genetree[1] #the "best" gene
three_clock = stats_clock$genetree[1:3] #the 3 "best" gene
nine_clock = stats_clock$genetree[1:9] #the 9 "best" gene

# Write selected genes to files -------------------------------------------
selected_dir = file.path(wd,"2_phylo_dating","1_gene_shop","selected_genes")
if (!dir.exists(selected_dir)){dir.create(selected_dir)}

write.table(no_gdis, file.path(selected_dir,"no_gdis.txt"),col.names=F,row.names = F,quote=F)
write.table(no_gdis_most_diff, file.path(selected_dir,"most_diff.txt"),col.names=F,row.names = F,quote=F)
write.table(one_clock, file.path(selected_dir,"clock_one.txt"),col.names=F,row.names = F,quote=F)
write.table(three_clock, file.path(selected_dir,"clock_three.txt"),col.names=F,row.names = F,quote=F)
write.table(nine_clock, file.path(selected_dir,"clock_nine.txt"),col.names=F,row.names = F,quote=F)
