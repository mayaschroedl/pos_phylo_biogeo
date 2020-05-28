###########################################################################
# Project: Orania Phylogeny MT
# Script: get_best_indiv.R
# --- Action: for the species where we have two sequences, we will only chose the sequence with the alignments with least gaps (for monophyletic sequences) for dating and two sequences for non-monophyletic sequences
# --- Input: 
# --- Output: 
# Author: Maya Schroedl
###########################################################################

rm(list=ls())

# Packages ----------------------------------------------------------------
if (!require('seqinr')) install.packages('seqinr'); library('seqinr')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('rlist')) install.packages('rlist'); library('rlist')
if (!require('ape')) install.packages('ape'); library('ape')

#not in
'%!in%' <- function(x,y)!('%in%'(x,y))

# Working directory -------------------------------------------------------
gwd = getwd()
wd = file.path(getwd(),"2_phylo_dating")

t = 2 #coverage trimming threshold

# Get best alignments for double names -------------------------------------------------

# Get double names 
tag_sp = read.table(file.path(gwd,"1_phylo_reconstruction","input_reads_and_info","tags_indiv_sp.txt"),h=T, stringsAsFactors = F) #table with which tag corresponds to which species

# Get alignments
all_alignments_concat = read.fasta(file.path(wd, "1_alignment", "concat", "all_genes_concat.fasta")) #get concatenated alignments for all individuals

get_gap_perc = function(indiv_num){ #get per individuum the percentage of gaps in alignment
  indiv_name = attr(all_alignments_concat[indiv_num],"name")
  empty = length(which(all_alignments_concat[indiv_num][[1]]=="-"))+length(which(all_alignments_concat[indiv_num][[1]]=="n")) #number of gaps in alignment
  
  perc = empty/length(all_alignments_concat[indiv_num][[1]])
  
  return(c(indiv_name, perc))
}


gap_perc = data.frame(do.call("rbind", lapply(1:length(all_alignments_concat), get_gap_perc)))

colnames(gap_perc) = c("tags","gap_perc")

gap_perc_merg = merge(gap_perc,tag_sp, by="tags") #merge with corresponding names


### GET THE INDIVIDUAL WITH THE LOWEST NUMBER OF GAPS IN ALIGNMENT
### Except for O. lauterbachiana & O. palindan. Because their individuals are not monophyletic. 

gap_perc_merg_without_ol_op = gap_perc_merg[-which(gap_perc_merg$indiv %in% c("Orania_lauterbachiana","Orania_palindan")),]


filtered_list_of_indiv = gap_perc_merg_without_ol_op %>%
  group_by(indiv) %>%
  summarize(gap_perc_min = min(gap_perc), tag_min = tags[which.min(gap_perc)])

final_vector_indiv = filtered_list_of_indiv$tag_min

## for O. palidan: keep both because only two individuals

final_vector_indiv = c(final_vector_indiv, tag_sp$tags[which(tag_sp$indiv == "Orania_palindan")])

## for O. lauterbachiana: keep TAG-24; and the one with the least gaps of the monophyletic group (TAG-22,TAG-23)

gap_perc_merg_laut_1_2 = gap_perc_merg[which(gap_perc_merg$tags %in% c("TAG-22","TAG-23")),]
laut_1_2_min = gap_perc_merg_laut_1_2$tags[which.min(gap_perc_merg_laut_1_2$gap_perc)]

final_vector_indiv = c(final_vector_indiv, "TAG-24", laut_1_2_min)


# Astral trees -------------------------------------------------------------
astral_gene_dir = file.path(wd, "2_astral_gene_trees")
if (!dir.exists(astral_gene_dir)){dir.create(astral_gene_dir)}

astral_tree = read.tree(file.path(gwd,"1_phylo_reconstruction", "4_coalescent_trees","2", "coalescent_lpp.tree_rooted"))

# drop tips that are not in final_vector_indiv

#which tips shall be droped:
drop_these = tag_sp$tags[which(tag_sp$tags %!in% final_vector_indiv)]

write.table(drop_these, file.path(astral_gene_dir, "astral_dropped_tags.txt"),row.names = F,col.names = F,quote = F)

astral_dropped = drop.tip(astral_tree,drop_these)

#write tree for dating
write.tree(astral_dropped, file.path(astral_gene_dir,"astral_for_dating.tree"))


# Gene trees --------------------------------------------------------------
# 
gene_dir = file.path(wd, "2_astral_gene_trees","gene_trees_drop")
if (!dir.exists(gene_dir)){dir.create(gene_dir)}

# drop tips in genetree (so we can compare genetree - sptree)
genetrees_folder = file.path(gwd, "1_phylo_reconstruction", "3_gene_trees", "2", "4_collapsed") #where are the gene trees
gene_list = read.table(file.path(gwd, "1_phylo_reconstruction", "genelist_7575.txt"))[,1] #list of all the genes
suffix = ".raxml.support.coll_rooted" #how are the gene tree files named


for (gene in gene_list){
  genetree = read.tree(file.path(genetrees_folder,paste0(gene, suffix)))
  genetree_dropped = drop.tip(genetree,drop_these)
  write.tree(genetree_dropped, file.path(gene_dir, paste0(gene, "_drop.tree")))
}





