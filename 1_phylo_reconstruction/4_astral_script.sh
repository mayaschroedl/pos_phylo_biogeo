#!/bin/bash
###########################################################################
# Project: Orania Phylogeny MT
# Script: alignment_script.sh
# --- Action: This script builts a coalescent tree using Astral.
# --- Input: genetrees from raxml; file with species names
# --- Output: coalescent tree
# Author: Maya Schroedl
# Date: 11/2019
###########################################################################


###################
#----ARGUMENTS----#
###################


#t:if you want trimmed contigs, based on coverage and threshold, enter t:threshold
#o: outgroup (TAG)
#l: labeled outgroup (individual name)
while getopts ":t:o:l:" opt; do
  case ${opt} in
    t)
      echo "-t was triggered, File: $OPTARG"
      t=$OPTARG
      ;;
o)
      echo "-o was triggered, Outgroup: $OPTARG"
      outgroup=$OPTARG
      ;;
	  
	  l)
      echo "-l was triggered, labeled Outgroup: $OPTARG"
      outgroup_lab=$OPTARG
      ;;
	  
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac

done

#################################
#----DIRECTORIES & VARIABLES----#
#################################

#working directories
GWD=$PWD #global working directory, with subprojects and scripts
WD="$GWD"/1_phylo_reconstruction #current working directory

if ! [ -z "$t" ]; #if there is a threshold definded, then work in the threshold directory
then
        dir_value=$t/
fi

mkdir -p $WD/4_coalescent_trees/$dir_value

################
#----ASTRAL----#
################

cd $WD/4_coalescent_trees/$dir_value
rm *

#ASTRAL tree with tips representing individuals (individuals of same species are not merged into one species tip) (extended sptree)

#execute with branch support as local posterior probabilities (lpp) [default]
java -jar $GWD/programs/Astral/astral.5.6.3.jar -i $WD/3_gene_trees/"$dir_value"4_collapsed/all_genes.raxml.support.coll -o $WD/4_coalescent_trees/"$dir_value"coalescent_lpp.tree 2>$WD/4_coalescent_trees/"$dir_value"coalescent_lpp.log 

#execute with branch support as all quartet supports [all alternatives] (qs) [-t 8]
java -jar $GWD/programs/Astral/astral.5.6.3.jar -t 8 -i $WD/3_gene_trees/"$dir_value"4_collapsed/all_genes.raxml.support.coll -o $WD/4_coalescent_trees/"$dir_value"coalescent_qs_all.tree 2>$WD/4_coalescent_trees/"$dir_value"coalescent_qs_all.log

#execute with branch support as quartet supports (qs) [-t 1]
java -jar $GWD/programs/Astral/astral.5.6.3.jar -t 1 -i $WD/3_gene_trees/"$dir_value"4_collapsed/all_genes.raxml.support.coll -o $WD/4_coalescent_trees/"$dir_value"coalescent_qs.tree 2>$WD/4_coalescent_trees/"$dir_value"coalescent_qs.log



#ASTRAL tree with tips representing species (individuals of same species are merged into one species tip) ("-a" otion = "ASTRAL-multi")

#execute with branch support as local posterior probabilities (lpp) [default]
java -jar $GWD/programs/Astral/astral.5.6.3.jar -i $WD/3_gene_trees/"$dir_value"4_collapsed/all_genes.raxml.support.coll -o $WD/4_coalescent_trees/"$dir_value"coalescent_lpp_idmerg.tree -a $WD/input_reads_and_info/astral_name_mapping.txt 2>$WD/4_coalescent_trees/"$dir_value"coalescent_lpp_idmerg.log 

#execute with branch support as all quartet supports [all alternatives] (qs) [-t 8]
java -jar $GWD/programs/Astral/astral.5.6.3.jar -t 8 -i $WD/3_gene_trees/"$dir_value"4_collapsed/all_genes.raxml.support.coll -o $WD/4_coalescent_trees/"$dir_value"coalescent_qs_all_idmerg.tree -a $WD/input_reads_and_info/astral_name_mapping.txt 2>$WD/4_coalescent_trees/"$dir_value"coalescent_qs_all_idmerg.log

#execute with branch support as quartet supports (qs) [-t 1]
java -jar $GWD/programs/Astral/astral.5.6.3.jar -t 1 -i $WD/3_gene_trees/"$dir_value"4_collapsed/all_genes.raxml.support.coll -o $WD/4_coalescent_trees/"$dir_value"coalescent_qs_idmerg.tree -a $WD/input_reads_and_info/astral_name_mapping.txt 2>$WD/4_coalescent_trees/"$dir_value"coalescent_qs_idmerg.log

#############################
#----LABEL SPECIES TREES----#
#############################

#### unlabeled individual trees (extended sp tree) --> add label (species names)
#with branch support as local posterior probabilities (lpp) [default]
Rscript $GWD/scripts/general/change_tiplabels.R $WD/4_coalescent_trees/"$dir_value"coalescent_lpp.tree $WD/4_coalescent_trees/"$dir_value"coalescent_lpp.tree_lab $WD/input_reads_and_info/tags_indiv.txt

#with branch support as all quartet supports [all alternatives] (qs) [-t 8]
Rscript $GWD/scripts/general/change_tiplabels.R $WD/4_coalescent_trees/"$dir_value"coalescent_qs_all.tree $WD/4_coalescent_trees/"$dir_value"coalescent_qs_all.tree_lab $WD/input_reads_and_info/tags_indiv.txt

#with branch support as quartet supports (qs) [-t 1]
Rscript $GWD/scripts/general/change_tiplabels.R $WD/4_coalescent_trees/"$dir_value"coalescent_qs.tree $WD/4_coalescent_trees/"$dir_value"coalescent_qs.tree_lab $WD/input_reads_and_info/tags_indiv.txt


######################
#----REROOT TREES----#
######################

######for unlabeled tree of individuals (extended sptree)
Rscript $GWD/scripts/general/root_tree.R $WD/4_coalescent_trees/"$dir_value"coalescent_lpp.tree $WD/4_coalescent_trees/"$dir_value"coalescent_lpp.tree_rooted $outgroup
 
Rscript $GWD/scripts/general/root_tree.R $WD/4_coalescent_trees/"$dir_value"coalescent_qs_all.tree $WD/4_coalescent_trees/"$dir_value"coalescent_qs_all.tree_rooted  $outgroup

Rscript $GWD/scripts/general/root_tree.R $WD/4_coalescent_trees/"$dir_value"coalescent_qs.tree $WD/4_coalescent_trees/"$dir_value"coalescent_qs.tree_rooted  $outgroup

#####for labeled tree of individuals (extended sptree)
Rscript $GWD/scripts/general/root_tree.R $WD/4_coalescent_trees/"$dir_value"coalescent_lpp.tree_lab $WD/4_coalescent_trees/"$dir_value"coalescent_lpp.tree_lab_rooted $outgroup_lab
 
Rscript $GWD/scripts/general/root_tree.R $WD/4_coalescent_trees/"$dir_value"coalescent_qs_all.tree_lab $WD/4_coalescent_trees/"$dir_value"coalescent_qs_all.tree_lab_rooted  $outgroup_lab

Rscript $GWD/scripts/general/root_tree.R $WD/4_coalescent_trees/"$dir_value"coalescent_qs.tree_lab $WD/4_coalescent_trees/"$dir_value"coalescent_qs.tree_lab_rooted  $outgroup_lab

#####for species tree ("-a" otion = "ASTRAL-multi")
Rscript $GWD/scripts/general/root_tree.R $WD/4_coalescent_trees/"$dir_value"coalescent_lpp_idmerg.tree $WD/4_coalescent_trees/"$dir_value"coalescent_lpp_idmerg.tree_rooted $outgroup_lab
 
Rscript $GWD/scripts/general/root_tree.R $WD/4_coalescent_trees/"$dir_value"coalescent_qs_idmerg.tree $WD/4_coalescent_trees/"$dir_value"coalescent_qs_idmerg.tree_rooted  $outgroup_lab

Rscript $GWD/scripts/general/root_tree.R $WD/4_coalescent_trees/"$dir_value"coalescent_qs_all_idmerg.tree $WD/4_coalescent_trees/"$dir_value"coalescent_qs_all_idmerg.tree_rooted  $outgroup_lab






cd $GWD