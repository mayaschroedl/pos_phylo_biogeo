#!/usr/bin/env Rscript
###########################################################################
# Project: Orania Phylogeny MT
# Script: change_tiplables_plot_trees.R
# --- Action: Changes the tiplabels of a tree from TAGs to individuals (species)
# --- Input: tree with TAGS tiplables (newick) + dataframe indicating which
# ---------- TAG belongs to which individual
# --- Output: tree with individual tiplables (newick)
# Author: Maya Schroedl
# Date: 11/2019
###########################################################################

rm(list=ls())

# Libraries ---------------------------------------------------------------
if (!require('ape')) install.packages('ape'); library('ape')
if (!require('phylotools')) install.packages('phylotools'); library('phylotools')



# Function -----------------------------------------------------------------
change_tip_labs=function(input_tree,tag_indiv_df=read.table(file.path(wd,"renamed_reads","tags_indiv.txt"),h=T) ){
  tree=read.tree(input_tree) #read the newick tree
  lab_tree=sub.taxa.label(tree, tag_indiv_df) #replace tiplables
  #output
  write.tree(lab_tree, paste0(substr(input_tree, 1, nchar(input_tree)-5),"_lab.tree"))
  plot=plot.phylo(lab_tree,main=basename(input_tree),show.node.label=T,use.edge.length=F)
  return(plot) #plot #basename: get only last part of path (here: filename)
}


# WD  --------------------------------------------------------------
gwd=getwd() # global working directory
wd=file.path(getwd(),"1_phylo_reconstruction") #local working directory
t=2 # threshold


# Make directories --------------------------------------------------------
dir_plots_phylo=file.path(gwd,"plots","phylogenies")
if (!dir.exists(dir_plots_phylo)){dir.create(dir_plots_phylo, recursive=T)}

# Species tree -------------------------------------------------------------
sptree_rooted=file.path(wd, "4_coalescent_trees", t, "coalescent_lpp.tree")

dev.off()
pdf(file.path(dir_plots_phylo,"coalescent_lpp_rooted.pdf"))
change_tip_labs(sptree_rooted)
dev.off()

# Gene trees --------------------------------------------------------------
genelist=read.table(file.path(wd,"genelist_7575.txt"))$V1

gene_tree_list=lapply(genelist,function(gene){return(file.path(wd, "3_gene_trees", t, "4_collapsed",paste0(gene,".raxml.support.coll")))})

dev.off()
pdf(file.path(dir_plots_phylo,"gene_trees_rooted.pdf"))
lapply(gene_tree_list,change_tip_labs)
dev.off()

