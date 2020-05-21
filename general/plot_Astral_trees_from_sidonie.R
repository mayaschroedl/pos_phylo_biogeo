#minor changes to script from Sidonie Bellot
#Maya Schrödl, 22/04

library(ape)
library(phytools)
library(gtools)
library(ggtree)
library(cowplot)
library(ggimage)
library(treeio) #read.astral function


# you probably don't need all the above packages but you certainly need ggimage, which requires an up-to-date R

setwd("D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/1_phylo_reconstruction/4_coalescent_trees/2") # set the working directory to the folder containing your tree

### FUNCTION

# import annotated ASTRAL tree (I like to use a tree rooted with pxrr directly after getting it from ASTRAL, before any editing)


astral_plot_qs_parts = function(AS_all){
# plot and ladderize the tree, without using the ASTRAL branch lengths
p <- ggtree(AS_all@phylo, ladderize=T, branch.length = "none") + geom_tiplab(size=2.5, hjust= -0.05) + xlim_tree(70)# + ggtitle("Astral tree of species QS parts")

# check node labels if necessary
#p + geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size = 4)

# format data to make pie charts (or any plot) of the Q scores (could do it for any other data from AS_all@data)
Q1 <- as.numeric(AS_all@data$q1) * 100
Q <- as.data.frame(Q1)
Q$Q2 <- as.numeric(AS_all@data$q2) * 100
Q$Q3 <- as.numeric(AS_all@data$q3) * 100
Q$node <- AS_all@data$node

# make barplots
# bars <- nodebar(Q, cols=1:3, position='dodge', color=c(Q1='red', Q2='cyan', Q3='gray'))
# inset(p, bars, x='node', width=8, height=2)

# make pie charts
pies <- nodepie(Q, cols=1:3, color=c(Q1='blue', Q2='green', Q3='red')) #, alpha=.6) # change style as needed, alpha is for transparency
inset(p, pies, width=0, height=0, hjust=.6)

# pies give the quartet support: percentage of quartets agreeing: 
# with the branch (blue), 
# with the second alternative RS|LO (green), 
# and with the last alternative RO|LS (red).
}


setwd("D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/plots/phylogenies")

#dev.off()
#pdf("ASTRAL_qs_parts.pdf")
############FOR SPECIES tree
setwd("D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/1_phylo_reconstruction/4_coalescent_trees/2")
sp_tree <- read.astral("coalescent_qs_all_idmerg.tree_rooted") #newick tree
astral_plot_qs_parts(sp_tree) + ggtitle("Astral tree of species QS part")

############FOR INDIVIDUALS tree
setwd("D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/1_phylo_reconstruction/4_coalescent_trees/2")
indiv_tree <- read.astral("coalescent_qs_all.tree_lab_rooted") #newick tree
astral_plot_qs_parts(indiv_tree) + ggtitle("Astral tree of individuals QS part")

#dev.off()
