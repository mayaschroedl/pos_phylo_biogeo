#quick script to plot first alternative topology for each node (results from phyparts)
#Maya Schroedl

#out.hist.alts of phypart needs to be manually transformed into table + newick format

rm(list=ls())

alt_topos = read.table("D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/1_phylo_reconstruction/5_phypartstopiecharts/2/out.hist.alts.sum.trees", h=T)

pdf("D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/1_phylo_reconstruction/5_phypartstopiecharts/2/phypart_first_alts.pdf")
for (node in unique(alt_topos$Node_num)){
 selected = alt_topos[alt_topos$Node_num == node,] # select the rows which correspond to the node
  first_alt_topo = as.character(selected[which.max(selected$Num_genetr_supp),]$Topology) # select the topology that is the most represented in the gene trees (first alternative topology)
  plot(read.tree(text=first_alt_topo), main = paste0("First alternative topology for node ", node, ". \n Number of agreeing gene trees: ", max(selected$Num_genetr_supp)))#plot this topology
}

dev.off()
