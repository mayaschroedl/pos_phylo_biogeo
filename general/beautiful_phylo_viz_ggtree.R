
#beautiful phylo viz

library("ggplot2")
library("ggtree")


# WD  --------------------------------------------------------------
gwd=getwd() # global working directory
wd=file.path(getwd(),"1_phylo_reconstruction") #local working directory
t=2 # threshold

# Make directories --------------------------------------------------------
dir_plots_phylo=file.path(gwd,"plots","phylogenies")
if (!dir.exists(dir_plots_phylo)){dir.create(dir_plots_phylo, recursive=T)}

# Species tree -------------------------------------------------------------
sptree_rooted=file.path(wd, "4_coalescent_trees", t, "coalescent_lpp_rooted_lab.tree")

tree=read.tree(sptree_rooted)
ggtree(tree)+
  geom_nodelab()+
  geom_tiplab() +
  #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
  xlim(0, 25)+

  geom_cladelabel(node=44, label="Madagascar", align=T, color='red', offset = 5) +
  geom_cladelabel(node=48, label="Africa", align=T, color='blue', offset = 5) +
  geom_cladelabel(node=31, label="SE Asia", align=T, color='green', offset = 5)
