###########################################################################
# Project: Orania Phylogeny MT
# Script: phyparts.sh
# --- Action:  plot first alternative topology for each node (results from phyparts)
# --- Input: out.hist.alts output from phyparts
# --- Output: first alternative topology plotted for each node
# Author: Maya Schroedl (maya.schroedl@bios.au.dk)

###########################################################################

rm(list=ls())


# WD + file preparation ---------------------------------------------------

wd = getwd()

#out.hist.alts of PhyParts needs to be manually transformed into table + newick format

alt_topos = read.table(file.path(wd,"1_phylo_reconstruction","5_phypartstopiecharts","2","out.hist.alts.sum.trees", h=T))


# First alternative -------------------------------------------------------
#pdf(file.path(wd,"1_phylo_reconstruction","5_phypartstopiecharts","2",phypart_first_alts.pdf"))
for (node in unique(alt_topos$Node_num)){
  
  selected = alt_topos[alt_topos$Node_num == node,] # select the rows which correspond to the node
  
  first_alt_topo = as.character(selected[which.max(selected$Num_genetr_supp),]$Topology) # select the topology that is the most represented in the gene trees (first alternative topology)
  
  plot(read.tree(text=first_alt_topo), main = paste0("First alternative topology for node ", node, ". \n Number of agreeing gene trees: ", max(selected$Num_genetr_supp)))#plot this topology
  
}

#dev.off()


# other alternatives for PS clade crown node ------------------------------
# just have a quick look on alternatives for node Podococcus-Sclerosperma

node=1 #ps clade crown
selected = alt_topos[alt_topos$Node_num == node,]

View(selected)

plot(read.tree(text=as.character(selected$Topology[4]))) 
#4: one of the other alternatives (not the major one) # you can check with the other topologies, Podococcus is always placed as sister with a few Oranias # but always supported only by one gene tree!

