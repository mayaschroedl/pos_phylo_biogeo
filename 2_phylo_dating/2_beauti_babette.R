###########################################################################
# Project: Orania Phylogeny MT
# Script: beast_babette.R
# --- Action: This script generates a .xml file as input for beast.
# ----------- works like BEAUTi
# --- Input: concatenated alignments of RF selected gene sequences
# --- Output: BEAUTi xml file
# Author: Maya Schroedl
###########################################################################

rm(list=ls())
# Packages ----------------------------------------------------------------
if (!require('babette')) install.packages('babette'); library('babette')


# Working directory -------------------------------------------------------
wd = file.path(getwd(),"2_phylo_dating")
dir_babette = file.path(wd, "4_babette")

if (!dir.exists(dir_babette)){dir.create(dir_babette)}


# Prepare astral starting tree --------------------------------------------

# Starting tree has to be an ultrametric tree. Starting tree has to conform with calibrations. Therefore, we have to scale the astral tree to the right scale.

# the script "2.1_chronopl_starting_tree.R" does this


# ALL GENES ---------------------------------------------------------------
all_genes_fsta = file.path(wd, "1_alignment","concat", "all_genes_concat.fasta")
all_genes_xml_out = file.path(dir_babette, "all_genes_concat.xml")

#BEAUTi
create_beast2_input_file(
  input_filename = all_genes_fsta,
  output_filename = all_genes_xml_out,
  mcmc=create_mcmc(chain_length = 10000000, store_every = 10000),
  mrca_prior=create_mrca_prior(
    is_monophyletic=T,
    mrca_distr=create_normal_distr(
      mean=create_mean_param(57),
      sigma=create_sigma_param(value=6.12))),
  site_model=create_hky_site_model(),
  clock_model=create_rln_clock_model(),
  tree_prior=create_yule_tree_prior()
)

#### ADD MANUALLY:####

# sometimes "mean=create_mean_param(57)" does not work 
# --> replace: "<parameter id="RealParameter.57" estimate="false" name="mean">0</parameter>" 
# --> with: "<parameter id="RealParameter.57" estimate="false" name="mean">57</parameter>"

### MULTI MONOPHYLETIC CONSTRAINT (MMC)

# after the node: <distribution id="prior" spec="util.CompoundDistribution">
# add: 
# <distribution id="MultiMonophyleticConstraint.t:all_genes_concat" spec="beast.math.distributions.MultiMonophyleticConstraint" isBinary="false" newick="((((TAG-44,TAG-41),(TAG-42,TAG-40)),((TAG-38,((TAG-28,(TAG-25,(TAG-24,(TAG-33,(TAG-20,TAG-21))))),((TAG-29,TAG-30),(TAG-26,(TAG-27,TAG-22))))),(TAG-31,(TAG-35,TAG-34)))),TAG-47);" tree="@Tree.t:all_genes_concat"/>

### STARTING TREE

#replace: <init id="RandomTree.t:all_genes_concat" spec="beast.evolution.tree.RandomTree" estimate="false" initial="@Tree.t:all_genes_concat" taxa="@all_genes_concat">
#<populationModel id="ConstantPopulation0.t:all_genes_concat" spec="ConstantPopulation">
#  <parameter id="randomPopSize.t:all_genes_concat" name="popSize">1.0</parameter>
#  </populationModel>
#  </init>
#(after </state>)
#with:
#<init id="NewickTree.t:all_genes_concat" spec="beast.util.TreeParser" IsLabelledNewick='true' adjustTipHeights='true' estimate="false" initial="@Tree.t:all_genes_concat" newick="((((TAG-44:17.00138032,TAG-41:17.00138032):25.01380443,(TAG-42:17.00475558,TAG-40:17.00475558):25.01042917):14.98481525,((TAG-38:24.62369906,((TAG-28:9.160026002,(TAG-25:7.287618064,(TAG-24:4.169573049,(TAG-33:1.48255441,(TAG-20:0.0008422456155,TAG-21:0.0008422456155):1.481712164):2.687018639):3.118045016):1.872407938):0.3282063006,((TAG-29:3.979652636,TAG-30:3.979652636):2.870878253,(TAG-26:3.287585618,TAG-27:3.287585618,TAG-22:3.287585618):3.562945271):2.637701414):15.13546675):15.23396707,(TAG-31:18.08655712,(TAG-35:9.343375489,TAG-34:9.343375489):8.743181633):21.771109):17.14233388):0,TAG-47:57);"/>


# Genes with no sptree discordant good nodes ------------------------------
no_gdis_fsta = file.path(wd, "1_alignment","concat", "no_gdis_concat.fasta")
no_gdis_xml_out = file.path(wd, "4_babette", "no_gdis_concat_try.xml")

#BEAUTi
create_beast2_input_file(
  input_filename = no_gdis_fsta,
  output_filename = no_gdis_xml_out,
  mcmc=create_mcmc(chain_length = 10000000, store_every = 10000),
  mrca_prior=create_mrca_prior(
    is_monophyletic=T,
    mrca_distr=create_normal_distr(
      mean=create_mean_param(57),
      sigma=create_sigma_param(value=6.12))),
  site_model=create_hky_site_model(),
  clock_model=create_rln_clock_model(),
  tree_prior=create_yule_tree_prior()
)

#### ADD MANUALLY:####

# like above

# Genes with some discordant good nodes -----------------------------------
# addititionally to no_gdis, here are the genes that have the greatest difference between the number of good nodes agreeing & disagreeing (see script gene_shop.R)

most_diff_fsta = file.path(wd, "1_alignment","concat", "most_diff_concat.fasta")
most_diff_xml_out = file.path(wd, "4_babette", "most_diff_concat.xml")

#BEAUTi
create_beast2_input_file(
  input_filename = most_diff_fsta,
  output_filename = most_diff_xml_out,
  mcmc=create_mcmc(chain_length = 10000000, store_every = 10000),
  mrca_prior=create_mrca_prior(
    is_monophyletic=T,
    mrca_distr=create_normal_distr(
      mean=create_mean_param(57),
      sigma=create_sigma_param(value=6.12))),
  site_model=create_hky_site_model(),
  clock_model=create_rln_clock_model(),
  tree_prior=create_yule_tree_prior()
)

#### ADD MANUALLY:####
#like above


# Model Selection ---------------------------------------------------------

# Genes with no sptree discordant good nodes ------------------------------
# 
###### ONE POSSIBLE CLOCK-LIKE GENE #####

clock_fsta = file.path(wd, "1_alignment","concat", "clock_one_concat.fasta")

## Relaxed clock ----
clock_relaxed_xml_out = file.path(wd, "4_babette", "clock_one_relaxed_mods.xml")

#BEAUTi
create_beast2_input_file(
  input_filename = clock_fsta,
  output_filename = clock_relaxed_xml_out,
  mcmc=create_mcmc(chain_length = 10000000, store_every = 10000),
  #mrca_prior=create_mrca_prior(
   # is_monophyletic=T,
   # mrca_distr=create_normal_distr(
   #   mean=create_mean_param(57),
   #   sigma=create_sigma_param(value=6.12))),
  site_model=create_hky_site_model(),
  clock_model=create_rln_clock_model(), #RELAXED MODEL
  tree_prior=create_yule_tree_prior()
)

# sometimes "mean=create_mean_param(57)" does not work 
# --> replace: "<parameter id="RealParameter.57" estimate="false" name="mean">0</parameter>" 
# --> with: "<parameter id="RealParameter.57" estimate="false" name="mean">57</parameter>"

# replace all "clock_one_concat" to "clock_one_relaxed" in doc

# !! ADD MMC + STARTING TREE MANUALLY (see above)

### FOR MODEL SELECTION:
# replace : <run id="mcmc" spec="MCMC" chainLength="1000000">
# with : <run id="mcmc" spec="beast.gss.NS" chainLength="20000" particleCount="1" subChainLength="5000">

## Strict clock ----
clock_strict_xml_out = file.path(wd, "4_babette", "clock_one_strict_mods.xml")

create_beast2_input_file(
  input_filename = clock_fsta,
  output_filename = clock_strict_xml_out,
  mcmc=create_mcmc(chain_length = 10000000, store_every = 10000),
  mrca_prior=create_mrca_prior(
    is_monophyletic=T,
    mrca_distr=create_normal_distr(
      mean=create_mean_param(57),
      sigma=create_sigma_param(value=6.12))),
  site_model=create_hky_site_model(),
  clock_model=create_strict_clock_model(), #STRICT MODEL
  tree_prior=create_yule_tree_prior()
)

# sometimes "mean=create_mean_param(57)" does not work 
# --> replace: "<parameter id="RealParameter.57" estimate="false" name="mean">0</parameter>" 
# --> with: "<parameter id="RealParameter.57" estimate="false" name="mean">57</parameter>"

# replace all "clock_one_concat" to "clock_one_strict" in doc

# !! ADD MMC + STARTING TREE MANUALLY (see above)

### FOR MODEL SELECTION:
# replace : <run id="mcmc" spec="MCMC" chainLength="1000000">
# with : <run id="mcmc" spec="beast.gss.NS" chainLength="20000" particleCount="1" subChainLength="5000">
# 
# 



###### THREE POSSIBLE CLOCK-LIKE GENES #####

clock_fsta = file.path(wd, "1_alignment","concat", "clock_three_concat.fasta")

## Relaxed clock ----
clock_relaxed_xml_out = file.path(wd, "4_babette", "clock_three_relaxed_mods.xml")

#BEAUTi
create_beast2_input_file(
  input_filename = clock_fsta,
  output_filename = clock_relaxed_xml_out,
  mcmc=create_mcmc(chain_length = 10000000, store_every = 10000),
  mrca_prior=create_mrca_prior(
    is_monophyletic=T,
    mrca_distr=create_normal_distr(
      mean=create_mean_param(57),
      sigma=create_sigma_param(value=6.12))),
  site_model=create_hky_site_model(),
  clock_model=create_rln_clock_model(), #RELAXED MODEL
  tree_prior=create_yule_tree_prior()
)

# sometimes "mean=create_mean_param(57)" does not work 
# --> replace: "<parameter id="RealParameter.57" estimate="false" name="mean">0</parameter>" 
# --> with: "<parameter id="RealParameter.57" estimate="false" name="mean">57</parameter>"

# replace all "clock_three_concat" to "clock_three_relaxed" in doc

# !! ADD MMC + STARTING TREE MANUALLY (see above)

### FOR MODEL SELECTION:
# replace : <run id="mcmc" spec="MCMC" chainLength="1000000">
# with : <run id="mcmc" spec="beast.gss.NS" chainLength="20000" particleCount="1" subChainLength="5000">

## Strict clock ----
clock_strict_xml_out = file.path(wd, "4_babette", "clock_three_strict_mods.xml")

create_beast2_input_file(
  input_filename = clock_fsta,
  output_filename = clock_strict_xml_out,
  mcmc=create_mcmc(chain_length = 10000000, store_every = 10000),
  #mrca_prior=create_mrca_prior(
   # is_monophyletic=T,
   # mrca_distr=create_normal_distr(
   #   mean=create_mean_param(57),
   #   sigma=create_sigma_param(value=6.12))),
  site_model=create_hky_site_model(),
  clock_model=create_strict_clock_model(), #STRICT MODEL
  tree_prior=create_yule_tree_prior()
)

# sometimes "mean=create_mean_param(57)" does not work 
# --> replace: "<parameter id="RealParameter.57" estimate="false" name="mean">0</parameter>" 
# --> with: "<parameter id="RealParameter.57" estimate="false" name="mean">57</parameter>"

# replace all "clock_three_concat" to "clock_three_strict" in doc

# !! ADD MMC + STARTING TREE MANUALLY (see above)

### FOR MODEL SELECTION:
# replace : <run id="mcmc" spec="MCMC" chainLength="1000000">
# with : <run id="mcmc" spec="beast.gss.NS" chainLength="20000" particleCount="1" subChainLength="5000">
# 



####### NINE POSSIBLE CLOCK-LIKE GENES #######

clock_fsta = file.path(wd, "1_alignment","concat", "clock_nine_concat.fasta")

## Relaxed clock ----
clock_relaxed_xml_out = file.path(wd, "4_babette", "clock_nine_relaxed_mods.xml")

#BEAUTi
create_beast2_input_file(
  input_filename = clock_fsta,
  output_filename = clock_relaxed_xml_out,
  mcmc=create_mcmc(chain_length = 10000000, store_every = 10000),
  mrca_prior=create_mrca_prior(
    is_monophyletic=T,
    mrca_distr=create_normal_distr(
      mean=create_mean_param(57),
      sigma=create_sigma_param(value=6.12))),
  site_model=create_hky_site_model(),
  clock_model=create_rln_clock_model(), #RELAXED MODEL
  tree_prior=create_yule_tree_prior()
)

# sometimes "mean=create_mean_param(57)" does not work 
# --> replace: "<parameter id="RealParameter.57" estimate="false" name="mean">0</parameter>" 
# --> with: "<parameter id="RealParameter.57" estimate="false" name="mean">57</parameter>"


# replace all "clock_nine_concat" to "clock_nine_relaxed" in doc

# ADD MMC + STARTING TREE MANUALLY

### FOR MODEL SELECTION:
# replace : <run id="mcmc" spec="MCMC" chainLength="1000000">
# with : <run id="mcmc" spec="beast.gss.NS" chainLength="20000" particleCount="1" subChainLength="5000">

## Strict clock ----
clock_strict_xml_out = file.path(wd, "4_babette", "clock_nine_strict_mods.xml")

create_beast2_input_file(
  input_filename = clock_fsta,
  output_filename = clock_strict_xml_out,
  mcmc=create_mcmc(chain_length = 10000000, store_every = 10000),
  mrca_prior=create_mrca_prior(
    is_monophyletic=T,
    mrca_distr=create_normal_distr(
      mean=create_mean_param(57),
      sigma=create_sigma_param(value=6.12))),
  site_model=create_hky_site_model(),
  clock_model=create_strict_clock_model(), #STRICT MODEL
  tree_prior=create_yule_tree_prior()
)

# sometimes "mean=create_mean_param(57)" does not work 
# --> replace: "<parameter id="RealParameter.57" estimate="false" name="mean">0</parameter>" 
# --> with: "<parameter id="RealParameter.57" estimate="false" name="mean">57</parameter>"


# replace all "clock_nine_concat" to "clock_nine_strict" in doc

# !! ADD MMC + STARTING TREE MANUALLY (see above)

### FOR MODEL SELECTION:
# replace : <run id="mcmc" spec="MCMC" chainLength="1000000">
# with : <run id="mcmc" spec="beast.gss.NS" chainLength="20000" particleCount="1" subChainLength="5000">
# 


