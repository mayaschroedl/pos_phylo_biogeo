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

setwd("D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/3_phylogeny_dating/GENOMEDK_Server/beast.or,orscldyp,sclpod_RF/babette")
fasta="or,orscldyp,sclpod_RF_selected_concat.fasta"

#BEAUTi
out = create_beast2_input_file(
  input_filename=fasta,
  output_filename = "babette.xml",
  mcmc=create_mcmc(chain_length = 10000, store_every = 1000),
  mrca_prior=create_mrca_prior(
    is_monophyletic=T,
    mrca_distr=create_normal_distr(
      mean=create_mean_param(67),
      sigma=create_sigma_param(value=9.1))),
  site_model=create_hky_site_model(),
  clock_model=create_strict_clock_model(),
  tree_prior=create_yule_tree_prior()
)

