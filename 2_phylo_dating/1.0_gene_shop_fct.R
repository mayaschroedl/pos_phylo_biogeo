#!/usr/bin/env Rscript
###########################################################################
# Project: Orania Phylogeny MT
# Script: gene_shop_fct.R
# --- Action: Function that compares gene trees to Astral tree (one direction) in order to 
# ------------ select the genes that have the least "good" (BS > 75%) nodes disagreeing with the species tree. Also select among these selected genes the most clock-like genes.
# --- Input: species tree, folder with individual rooted gene trees (with outgroup! + collapsed nodes (BS < 10) and collapsed nodes (branch length < 0.00002))
# --- Output: different statistics on how much the genetree agrees/disagrees with sptree; and how clock-like the genetree is in relation to others
# Author: Maya Schroedl
###########################################################################

rm(list = ls()) #clear environment

addTaskCallback(function(...){set.seed(42);TRUE}) #set seed to 42 for whole doc

# Libraries ---------------------------------------------------------------
if (!require('ape')) install.packages('ape'); library('ape')
if (!require('phylotools')) install.packages('phylotools'); library('phylotools')
if (!require('adephylo')) install.packages('adephylo'); library('adephylo')
if (!require('gtools')) install.packages('gtools'); library('gtools')

# Function ----------------------------------------------------------------

signal_support_stats=function(genetree, gene_name, sptree, output)
{
  #genetree: newick tree file of the genetree, rooted with ape(edgelabel = TRUE,resolve.root=TRUE)
  #gene_name: name of the gene
  #sptree:  newick tree file of the sptree, rooted with ape(edgelabel = TRUE,resolve.root=TRUE)
  #output: directory of the outputfile
  
  ### Prepare data and dataframes ----
  
  genetree_original = genetree # to store original genetree, which we will not modify
  
  #--- Get node numbers
  genetree$node.nums = seq(min(genetree$edge[,1]),max(genetree$edge[,1]))# extract sequence of node numbers
  genetree$node.label # bootstrap values for each node
  
  # dataframe that summarizes node numbers + support values (node.label)
  nodes_df = data.frame(node.nums = genetree$node.nums, node.label = genetree$node.label)
  
  # do not take into consideration the nodes that have no bootstrap value (outgroup)
    no_support_num =   which(nodes_df$node.label == "")
    no_support = length(no_support_num) # how many nodes do not have bootstrap value
    nodes_df = nodes_df[-no_support_num,]
  
    # correct node numbers and node labels in gene tree
    genetree$node.nums = nodes_df$node.nums
    genetree$node.label = nodes_df$node.label
  
  #--- Weight relative percentages with polytomies
  
  # Get "theoretical number of nodes (without polytomies)
  node_num_theo = multi2di(genetree)$Nnode - no_support # resolve randomly polytomies and get the number of theoretical nodes (excluding the nodes without support)
  
  # Node dataframe with weighted polytomies (node with polytomy counts x2 (for trichotomy), x3 etc.) to calculate percentage
  polyt=as.data.frame(table(genetree$edge[,1])-1)[-no_support_num,] # df of which nodes present a polytomy (Freq = 2 means that the node is trichotomous)
  merged = merge(nodes_df,polyt, by.x="node.nums", by.y= "Var1") # merde node df and polytomie df
  nodes_df_weighted = as.data.frame(lapply(merged, rep, merged$Freq)) # df where nodes are as many times represented as "theoretical" nodes would be [if polytomy resolved]
  
  
  ### 1) How informative is the gene tree ----
  
  #--- Which nodes are "good" (BS >= 75)
  gnodes = nodes_df$node.nums[which(nodes_df$node.label >= 75)] # select nodes that have a bootstrap value > 75
  
  gnodes_weighted = nodes_df_weighted$node.nums[which(nodes_df_weighted$node.label >= 75)] # list of good nodes, weighted by the polytomies (node numbers multiplicated)
  
  #--- How many "good" nodes (BS >= 75) [as a percentage relative to theoretical number of nodes] (weighted by polytomies)
  
  gnodes_num = length(which(nodes_df_weighted$node.label >= 75)) # how many nodes are "good" (BS > 75)
  gnodes_perc = gnodes_num/node_num_theo*100 # relative percentage of good nodes for this gene tree (weighted by polytomies)
  
  ### 2) How many nodes agree/disagree with sptree? ----
  #--- Is the clade corresponding to node X monophyletic in the species tree? (agrees)
    #Not on polytomies: we only look at the crown node of the polytomy and see wether its tiplabels are a monophyletic group in the sp tree. the result (T or F) is weighted by the number of nodes a resolved polytomy would have (nodes_df_weighted). !ATTENTION! this analysis is very sensitive to when one species is not in the polytomy clade in the sptree, then the polytomy-node is counted as disagreeing with sptree

      sptree_agrees = function(node_num){
        selected_clade_tree = extract.clade(genetree,node_num)
        return(is.monophyletic(sptree,selected_clade_tree$tip.label))}
  
  nd_agree = sapply(nodes_df_weighted$node.nums, sptree_agrees) # which nodes agree with sptree (weighted)
  nd_agree_perc = length(which(nd_agree==T))/node_num_theo*100 # percentage of nodes that agree with sptree (weighted)
  
  nd_disagree = !(sapply(nodes_df_weighted$node.nums, sptree_agrees)) # which nodes disagree with sptree (weighted)
  nd_disagree_perc = length(which(nd_disagree==T))/node_num_theo*100 # percentage of nodes that disagree with sptree (weighted)
  

  ### 3) How many nodes are "good" (BS >= 75) and agree/disagree with the species tree ---
  gnd_agree = sapply(gnodes_weighted, sptree_agrees) # which good nodes agree with sptree (weighted)
  gnd_agree_perc = length(which(gnd_agree==T))/node_num_theo*100 # percentage of nodes that are good nodes and agree with sptree (weighted)
  
  #--- Is the clade corresponding to node X not monophyletic in the species tree? (disagrees)
  gnd_disagree = !(sapply(gnodes_weighted, sptree_agrees)) # which good nodes disagree with sptree (weighted)
  gnd_disagree_perc = length(which(gnd_disagree==T))/node_num_theo*100 # percentage of nodes that are good nodes and disagree with sptree (weighted)
  gnd_disagree_nbr = length(which(gnd_disagree==T)) # how many good nodes disagree


  ###  clock-likeness -----
  # Root-tip-variance as proxy of "clock-likeness". Check: Smith, S. A., Brown, J. W., & Walker, J. F. (2018). So many genes, so little time: a practical approach to divergence-time estimation in the genomic era. PloS one, 13(5). https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0197433 
  dist_root = var(adephylo::distRoot(genetree_original))
  
  ### Export results ----
  resultdf=data.frame("genetree" = gene_name, "gnodes_perc" = round(gnodes_perc,3),"nd_agree_perc" = nd_agree_perc, "nd_disagree_perc" = nd_disagree_perc,"gnd_agree_perc" = round(gnd_agree_perc,3),"gnd_disagree_perc" = round(gnd_disagree_perc,3), "gnd_disagree_nbr" = round(gnd_disagree_nbr,3), "diff_ga_gd" = round(gnd_agree_perc - gnd_disagree_perc,3), "dist_root" = round(dist_root*100000,3) )
  
  # if output file exists and is not empty:
  if ((file.exists(output)) & (file.info(output)$size != 0)){
    # append resultdf to existing file
    write.table(resultdf,output,append = T, quote = F,  row.names=FALSE, col.names=FALSE)} else {
    # create a new file with resultdf
    write.table(resultdf,output, quote = F, row.names=FALSE)
    }
  
} # end of function

