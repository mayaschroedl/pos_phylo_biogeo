###########################################################################
# Project: Orania Phylogeny MT
# Script: beast_babette.R
# --- Action: for the species where we have two sequences, we will only chose the sequence with the "most-informative" alignments (least "-" and "n") for dating
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



# Get double names -------------------------------------------------------------------
tag_sp_indiv = read.table(file.path(gwd,"1_phylo_reconstruction","input_reads_and_info","tags_indiv.txt"),h=T, stringsAsFactors = F) #table with which tag corresponds to which species

# Astral tree -------------------------------------------------------------

astral_idmerg = read.tree (file.path(gwd,"1_phylo_reconstruction", "4_coalescent_trees", t, "coalescent_lpp_idmerg.tree_rooted" ))

plot(astral_idmerg, main = "Astral tree of species")
nodelabels(astral_idmerg$node.label, frame= "none")

astral_dir = file.path(wd, "2_astral_tree")

if (!dir.exists(astral_dir)){dir.create(astral_dir)}

file.copy(file.path(gwd,"1_phylo_reconstruction", "4_coalescent_trees", t, "coalescent_lpp.tree_rooted" ),astral_dir)


#change tip labels
source(file.path(gwd,"scripts","general","change_tiplabels.R"))

change_tip_labs(file.path(astral_dir,"coalescent_lpp.tree_rooted"), file.path(astral_dir,"coalescent_lpp.tree_lab_rooted"),tag_sp_indiv)

astral_sp = read.tree(file.path(astral_dir,"coalescent_lpp.tree_lab_rooted"))

plot(astral_sp, main = "Astral tree of individuals")
nodelabels(astral_sp$node.label, frame= "none")



# get indiv most similar to astral tree in top ----------------------------

# we exclude those sequences that are placed in a different place in the astral tree without merged individuals (coalescent_lpp.tree_rooted) than in the astral tree where indivs are merged (coalescent_lpp_imerg.tree_rooted) [astral "-a" option]

# in our case, this is Orania_palidan_2 (TAG-26)

# one special case in Orania lauterbachiana. All three individual sequences are in a different place in coalescent_lpp.tree_rooted than in coalescent_lpp_imerg.tree_rooted. I did not find what the astral "-a" option exactly does, so I do not know how it deals with individuals that are in very different positions (as here in Orania lauterbachiana).
# Because Orania lauterbachiana has two individuals that occuras a sister group to O. tabubilensis, and only one [O. lauterbachiana_3] that is sister to O. deflexa/archboldiana/macropetala, we excluded this last individual.



# Get best alignments for double names -------------------------------------------------

all_alignments_concat = read.fasta(file.path(wd, "1_alignment", "concat", "all_genes_concat.fasta"))

perc_empty = function(indiv_num){
  indiv_name = attr(all_alignments_concat[indiv_num],"name")
  empty = length(which(all_alignments_concat[indiv_num][[1]]=="-"))+length(which(all_alignments_concat[indiv_num][[1]]=="n"))
  
  perc = empty/length(all_alignments_concat[indiv_num][[1]])
  
  return(c(indiv_name, perc))
}


all_names_pe = data.frame(do.call("rbind", lapply(1:length(all_alignments_concat), perc_empty)))




# Best alignments for double names ----------------------------------------

get_best = function(elem){
  return(
    all_names_pe[,1][which(all_names_pe[,1] %in% elem)
                     [which.max(all_names_pe[,2]
                                [which(all_names_pe[,1] %in% elem)
                                  ]
                     )
                       ]
                     ]
    
  )
}

list2 = c()
len = length(list)
for (i in 1:len){
  elem = list[[i]][[2]]
  list2 = c(list2,as.character(get_best(elem)))
}

##### astral tree

# drop tips that are not in list2



#which tips shall be droped:
drop_these = tag_sp$tags[which(tag_sp$tags %!in% list2)]

astral_dropped = drop.tip(astral_tree,drop_these)

write.tree(astral_dropped, file.path(astral_dir,"coalescent_lpp_drop.tree_rooted"))


#change tip labels
source(file.path(gwd,"scripts","general","change_tiplabels.R"))

change_tip_labs(file.path(astral_dir,"coalescent_lpp_drop.tree_rooted"), file.path(astral_dir,"coalescent_lpp_drop.tree_lab_rooted"),tag_sp)


plot(read.tree(file.path(astral_dir,"coalescent_lpp_drop.tree_lab_rooted")))

plot(read.tree(file.path(gwd,"1_phylo_reconstruction", "4_coalescent_trees", t, "coalescent_lpp_idmerg.tree_rooted" )))





