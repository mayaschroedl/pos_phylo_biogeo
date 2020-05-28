
#beautiful phylo viz

library("ggplot2")
library("ggtree")

library("ape")

# WD  --------------------------------------------------------------
gwd=getwd() # global working directory
wd=file.path(getwd(),"1_phylo_reconstruction") #local working directory
t=2 # threshold

# Directories --------------------------------------------------------
dir_plots_phylo=file.path(gwd,"plots","phylogenies")
if (!dir.exists(dir_plots_phylo)){dir.create(dir_plots_phylo, recursive=T)}

# Astral trees -------------------------------------------------------------
#dev.off()
#pdf(file.path(gwd,"plots","phylogenies","all_astral_trees.pdf"))
### ASTRAL IDMERG
dev.off()
svg(file.path(gwd,"plots","phylogenies","astral_sp_final.svg"))
#lpp
astral_idmerg_lpp = read.tree (file.path(gwd,"1_phylo_reconstruction", "4_coalescent_trees", t, "coalescent_lpp_idmerg.tree_rooted" ))

p=ggtree(astral_idmerg_lpp)+geom_tiplab(offset = .9, hjust = .5)#+geom_text2(aes(subset = !isTip, label=label))

p = rotate(p,20)
p = rotate(p,22)
p = rotate(p,31)

p

#plot(astral_idmerg_lpp, main = "Astral tree of species LP", use.edge.length = T)
#nodelabels(astral_idmerg_lpp$node.label, frame= "none")

dev.off()


### ASTRAL INDIV
dev.off()
svg(file.path(gwd,"plots","phylogenies","astral_indiv_final.svg"))
#lpp
astral_idmerg_lpp = read.tree (file.path(gwd,"1_phylo_reconstruction", "4_coalescent_trees", t, "coalescent_lpp.tree_lab_rooted" ))

p=ggtree(astral_idmerg_lpp)+geom_tiplab(offset = .9, hjust = .5)#+geom_text2(aes(subset = !isTip, label=label))

p = rotate(p,29)
p = rotate(p,30)
p = rotate(p,50)
p = rotate(p,37)
p = rotate(p,48)

p

dev.off()

# Gene trees --------------------------------------------------------------

genelist = read.table(file.path(wd, "genelist_7575.txt"))[,1]

pdf(file.path(dir_plots_phylo,"genetrees_rooted.pdf"))
for (gene in genelist){
  
  genetree_rooted_file=file.path(wd, "3_gene_trees", t, "4_collapsed", paste0(gene, ".raxml.support.coll_lab_rooted"))
  tree_rooted=read.tree(genetree_rooted_file)
  
  plot.phylo(tree_rooted, main= gene)
  nodelabels(text=tree_rooted$node.label, frame="none")
  
}

dev.off()


# Genes selected for dating -----------------------------------------------

genelist = read.table(file.path(gwd, "2_phylo_dating","1_gene_shop","selected_genes", "no_gdis.txt"))[,1]

pdf(file.path(dir_plots_phylo,"no_gdis_rooted.pdf"))
for (gene in genelist){
  
  genetree_rooted_file=file.path(wd, "3_gene_trees", t, "4_collapsed", paste0(gene, ".raxml.support.coll_lab_rooted"))
  tree_rooted=read.tree(genetree_rooted_file)
  
  plot.phylo(tree_rooted, main= gene)
  nodelabels(text=tree_rooted$node.label, frame="none")
  
}

dev.off()


# ggtree(tree)+
#   geom_nodelab()+
#   geom_tiplab() +
#   #geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
#   xlim(0, 25)#+
# 
#   #geom_cladelabel(node=44, label="Madagascar", align=T, color='red', offset = 5) +
#   #geom_cladelabel(node=48, label="Africa", align=T, color='blue', offset = 5) +
#   #geom_cladelabel(node=31, label="SE Asia", align=T, color='green', offset = 5)


