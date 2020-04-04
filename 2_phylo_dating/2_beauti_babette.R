###########################################################################
# Project: Orania Phylogeny MT
# Script: beast_babette.R
# --- Action: This script generates a .xml file as input for beast.
# ----------- works like BEAUTi
# --- Input: concatenated alignments of RF selected gene sequences
# --- Output: BEAUTi xml file
# Author: Maya Schroedl
# Date: 11/2019
###########################################################################

install.packages("babette")
library(babette)

install.packages("XML")
library(XML)

install.packages("xml2")
library(xml2)

setwd("D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/2_phylo_dating/babette")
fasta="or,orscldyp,sclpod_RF_selected_concat.fasta"

#BEAUTi
out = create_beast2_input_file(
  input_filename=fasta,
  output_filename = "babette.xml",
  mcmc=create_mcmc(chain_length = 10000, store_every = 1000),
  mrca_prior=create_mrca_prior(
    is_monophyletic=T,
    mrca_distr=create_normal_distr(
      mean=create_mean_param(57),
      sigma=create_sigma_param(value=6.12))),
  site_model=create_hky_site_model(),
  clock_model=create_relaxed_clock_model(),
  tree_prior=create_yule_tree_prior()
)

#### ADD MULTI MONOPHYLETIC CONSTRAINT MANUALLY:####
# after the node: <distribution id="prior" spec="util.CompoundDistribution">
#add: <distribution id="MultiMonophyleticConstraint.t:or,orscldyp,sclpod_RF_selected_concat" spec="beast.math.distributions.MultiMonophyleticConstraint" isBinary="false" newick="(TAG-47,(((TAG-44,(TAG-43,TAG-41)),(TAG-46,(TAG-42,(TAG-39,TAG-40)))),(((TAG-38,(TAG-37,TAG-36)),(TAG-30,(TAG-23,(TAG-28,(TAG-26,(TAG-25,(TAG-22,(TAG-29,((TAG-27,TAG-24),(TAG-21,(TAG-20,TAG-33))))))))))),(TAG-31,(TAG-35,(TAG-34,(TAG-45,TAG-32))))))); " tree="@Tree.t:or,orscldyp,sclpod_RF_selected_concat"/>
