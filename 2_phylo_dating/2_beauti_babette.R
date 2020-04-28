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

# Packages ----------------------------------------------------------------
if (!require('babette')) install.packages('babette'); library('babette')


# Working directory -------------------------------------------------------
wd = file.path(getwd(),"2_phylo_dating")
dir_babette = file.path(wd, "3_babette")

if (!dir.exists(dir_babette)){dir.create(dir_babette)}


# Prepare astral starting tree --------------------------------------------

# Starting tree has to be an ultrametric tree. Starting tree has to conform with calibrations. Therefore, we have to scale the astral tree to the right scale.

# the script "2.1_chronopl_starting_tree.R" does this


# ALL GENES ---------------------------------------------------------------
all_genes_fsta = file.path(wd, "2_alignment","concat", "all_genes_concat.fasta")
all_genes_xml_out = file.path(wd, "3_babette", "all_genes_concat.xml")

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
# <distribution id="MultiMonophyleticConstraint.t:all_genes_concat" spec="beast.math.distributions.MultiMonophyleticConstraint" isBinary="false" newick="(((((TAG-41,TAG-43),TAG-44),(TAG-42,(TAG-46,(TAG-39,TAG-40)))),((TAG-38,(TAG-37,((TAG-28,(TAG-25,(TAG-24,(TAG-33,(TAG-20,TAG-21))))),((TAG-29,TAG-30),((TAG-23,TAG-22),(TAG-26,TAG-27)))))),(TAG-31,(TAG-35,(TAG-34,(TAG-45,TAG-32)))))),TAG-47);" tree="@Tree.t:all_genes_concat"/>

### STARTING TREE

#replace: <init id="RandomTree.t:all_genes_concat" spec="beast.evolution.tree.RandomTree" estimate="false" initial="@Tree.t:all_genes_concat" taxa="@all_genes_concat">
#<populationModel id="ConstantPopulation0.t:all_genes_concat" spec="ConstantPopulation">
#  <parameter id="randomPopSize.t:all_genes_concat" name="popSize">1.0</parameter>
#  </populationModel>
#  </init>
#(after </state>)
#with:
#<init id="NewickTree.t:all_genes_concat" spec="beast.util.TreeParser" IsLabelledNewick='true' adjustTipHeights='true' estimate="false" initial="@Tree.t:all_genes_concat" newick="(((((TAG-41:14.09551088,TAG-43:14.09551088)1:10.72297743,TAG-44:24.81848831)1:19.39735748,(TAG-42:18.733924,(TAG-46:11.67897662,(TAG-39:5.056527248,TAG-40:5.056527248)0.97:6.622449377)0.85:7.054947371)1:25.48192179)0.81:10.83781918,((TAG-38:18.52086344,(TAG-37:16.80008428,((TAG-28:5.129279758,(TAG-25:4.942861574,(TAG-24:2.464106376,(TAG-33:0.805392906,(TAG-20:6.108079699e-05,TAG-21:6.108079699e-05)0.91:0.8053318252)0.99:1.65871347)0.99:2.478755198)0.43:0.1864181843)0.43:0.4987030302,((TAG-29:1.06229447,TAG-30:1.06229447)0.47:3.44312657,((TAG-23:6.010985656e-05,TAG-22:6.010985656e-05)1:3.488963273,(TAG-26:1.523530584,TAG-27:1.523530584)0.49:1.965492799)0.76:1.016397657)0.69:1.122561749)1:11.17210149)0.47:1.720779157)1:18.7067502,(TAG-31:18.03998218,(TAG-35:10.9702256,(TAG-34:3.932632282,(TAG-45:1.386045447,TAG-32:1.386045447)0.84:2.546586835)1:7.037593316)1:7.069756581)1:19.18763146)1:17.82605133):0,TAG-47:55.05366497);"/>


# Genes with no sptree discordant good nodes ------------------------------
no_gdis_fsta = file.path(wd, "2_alignment","concat", "no_gdis_concat.fasta")
no_gdis_xml_out = file.path(wd, "3_babette", "no_gdis_concat.xml")

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

# sometimes "mean=create_mean_param(57)" does not work 
# --> replace: "<parameter id="RealParameter.57" estimate="false" name="mean">0</parameter>" 
# --> with: "<parameter id="RealParameter.57" estimate="false" name="mean">57</parameter>"

### MULTI MONOPHYLETIC CONSTRAINT (MMC)

# after the node: <distribution id="prior" spec="util.CompoundDistribution">
# add: 
# <distribution id="MultiMonophyleticConstraint.t:no_gdis_concat" spec="beast.math.distributions.MultiMonophyleticConstraint" isBinary="false" newick="(((((TAG-41,TAG-43),TAG-44),(TAG-42,(TAG-46,(TAG-39,TAG-40)))),((TAG-38,(TAG-37,((TAG-28,(TAG-25,(TAG-24,(TAG-33,(TAG-20,TAG-21))))),((TAG-29,TAG-30),((TAG-23,TAG-22),(TAG-26,TAG-27)))))),(TAG-31,(TAG-35,(TAG-34,(TAG-45,TAG-32)))))),TAG-47);" tree="@Tree.t:no_gdis_concat"/>

### STARTING TREE

#replace: <init id="RandomTree.t:no_gdis_concat" spec="beast.evolution.tree.RandomTree" estimate="false" initial="@Tree.t:no_gdis_concat" taxa="@no_gdis_concat">
#<populationModel id="ConstantPopulation0.t:no_gdis_concat" spec="ConstantPopulation">
#  <parameter id="randomPopSize.t:no_gdis_concat" name="popSize">1.0</parameter>
#  </populationModel>
#  </init>
#(after </state>)
#with:
#<init id="NewickTree.t:no_gdis_concat" spec="beast.util.TreeParser" IsLabelledNewick='true' adjustTipHeights='true' estimate="false" initial="@Tree.t:no_gdis_concat" newick="(((((TAG-41:14.09551088,TAG-43:14.09551088)1:10.72297743,TAG-44:24.81848831)1:19.39735748,(TAG-42:18.733924,(TAG-46:11.67897662,(TAG-39:5.056527248,TAG-40:5.056527248)0.97:6.622449377)0.85:7.054947371)1:25.48192179)0.81:10.83781918,((TAG-38:18.52086344,(TAG-37:16.80008428,((TAG-28:5.129279758,(TAG-25:4.942861574,(TAG-24:2.464106376,(TAG-33:0.805392906,(TAG-20:6.108079699e-05,TAG-21:6.108079699e-05)0.91:0.8053318252)0.99:1.65871347)0.99:2.478755198)0.43:0.1864181843)0.43:0.4987030302,((TAG-29:1.06229447,TAG-30:1.06229447)0.47:3.44312657,((TAG-23:6.010985656e-05,TAG-22:6.010985656e-05)1:3.488963273,(TAG-26:1.523530584,TAG-27:1.523530584)0.49:1.965492799)0.76:1.016397657)0.69:1.122561749)1:11.17210149)0.47:1.720779157)1:18.7067502,(TAG-31:18.03998218,(TAG-35:10.9702256,(TAG-34:3.932632282,(TAG-45:1.386045447,TAG-32:1.386045447)0.84:2.546586835)1:7.037593316)1:7.069756581)1:19.18763146)1:17.82605133):0,TAG-47:55.05366497);"/>


# Genes with some discordant good nodes -----------------------------------
# addititionally to no_gdis, here are the genes that have the greatest difference between the number of good nodes agreeing & disagreeing (see script gene_shop.R)

most_diff_fsta = file.path(wd, "2_alignment","concat", "most_diff_concat.fasta")
most_diff_xml_out = file.path(wd, "3_babette", "most_diff_concat.xml")

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

# sometimes "mean=create_mean_param(57)" does not work 
# --> replace: "<parameter id="RealParameter.57" estimate="false" name="mean">0</parameter>" 
# --> with: "<parameter id="RealParameter.57" estimate="false" name="mean">57</parameter>"


### MULTI MONOPHYLETIC CONSTRAINT (MMC)

# after the node: <distribution id="prior" spec="util.CompoundDistribution">
# add: 
# <distribution id="MultiMonophyleticConstraint.t:most_diff_concat" spec="beast.math.distributions.MultiMonophyleticConstraint" isBinary="false" newick="(((((TAG-41,TAG-43),TAG-44),(TAG-42,(TAG-46,(TAG-39,TAG-40)))),((TAG-38,(TAG-37,((TAG-28,(TAG-25,(TAG-24,(TAG-33,(TAG-20,TAG-21))))),((TAG-29,TAG-30),((TAG-23,TAG-22),(TAG-26,TAG-27)))))),(TAG-31,(TAG-35,(TAG-34,(TAG-45,TAG-32)))))),TAG-47);" tree="@Tree.t:most_diff_concat"/>

### STARTING TREE

#replace: <init id="RandomTree.t:most_diff_concat" spec="beast.evolution.tree.RandomTree" estimate="false" initial="@Tree.t:most_diff_concat" taxa="@most_diff_concat">
#<populationModel id="ConstantPopulation0.t:most_diff_concat" spec="ConstantPopulation">
#  <parameter id="randomPopSize.t:most_diff_concat" name="popSize">1.0</parameter>
#  </populationModel>
#  </init>
#(after </state>)
#with:
#<init id="NewickTree.t:most_diff_concat" spec="beast.util.TreeParser" IsLabelledNewick='true' adjustTipHeights='true' estimate="false" initial="@Tree.t:most_diff_concat" newick="(((((TAG-41:14.09551088,TAG-43:14.09551088)1:10.72297743,TAG-44:24.81848831)1:19.39735748,(TAG-42:18.733924,(TAG-46:11.67897662,(TAG-39:5.056527248,TAG-40:5.056527248)0.97:6.622449377)0.85:7.054947371)1:25.48192179)0.81:10.83781918,((TAG-38:18.52086344,(TAG-37:16.80008428,((TAG-28:5.129279758,(TAG-25:4.942861574,(TAG-24:2.464106376,(TAG-33:0.805392906,(TAG-20:6.108079699e-05,TAG-21:6.108079699e-05)0.91:0.8053318252)0.99:1.65871347)0.99:2.478755198)0.43:0.1864181843)0.43:0.4987030302,((TAG-29:1.06229447,TAG-30:1.06229447)0.47:3.44312657,((TAG-23:6.010985656e-05,TAG-22:6.010985656e-05)1:3.488963273,(TAG-26:1.523530584,TAG-27:1.523530584)0.49:1.965492799)0.76:1.016397657)0.69:1.122561749)1:11.17210149)0.47:1.720779157)1:18.7067502,(TAG-31:18.03998218,(TAG-35:10.9702256,(TAG-34:3.932632282,(TAG-45:1.386045447,TAG-32:1.386045447)0.84:2.546586835)1:7.037593316)1:7.069756581)1:19.18763146)1:17.82605133):0,TAG-47:55.05366497);"/>


# Model Selection ---------------------------------------------------------

# Genes with no sptree discordant good nodes ------------------------------

####### NINE POSSIBLE CLOCK-LIKE GENES #######

clock_fsta = file.path(wd, "2_alignment","concat", "clock_nine_concat.fasta")

## Relaxed clock ----
clock_relaxed_xml_out = file.path(wd, "3_babette", "clock_nine_relaxed_mods.xml")

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
clock_strict_xml_out = file.path(wd, "3_babette", "clock_nine_strict_mods.xml")

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


###### THREE POSSIBLE CLOCK-LIKE GENES #####

clock_fsta = file.path(wd, "2_alignment","concat", "clock_three_concat.fasta")

## Relaxed clock ----
clock_relaxed_xml_out = file.path(wd, "3_babette", "clock_three_relaxed_mods.xml")

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
clock_strict_xml_out = file.path(wd, "3_babette", "clock_three_strict_mods.xml")

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

# replace all "clock_three_concat" to "clock_three_strict" in doc

# !! ADD MMC + STARTING TREE MANUALLY (see above)

### FOR MODEL SELECTION:
# replace : <run id="mcmc" spec="MCMC" chainLength="1000000">
# with : <run id="mcmc" spec="beast.gss.NS" chainLength="20000" particleCount="1" subChainLength="5000">
# 


###### ONE POSSIBLE CLOCK-LIKE GENE #####

clock_fsta = file.path(wd, "2_alignment","concat", "clock_one_concat.fasta")

## Relaxed clock ----
clock_relaxed_xml_out = file.path(wd, "3_babette", "clock_one_relaxed_mods.xml")

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

# replace all "clock_one_concat" to "clock_one_relaxed" in doc

# !! ADD MMC + STARTING TREE MANUALLY (see above)

### FOR MODEL SELECTION:
# replace : <run id="mcmc" spec="MCMC" chainLength="1000000">
# with : <run id="mcmc" spec="beast.gss.NS" chainLength="20000" particleCount="1" subChainLength="5000">

## Strict clock ----
clock_strict_xml_out = file.path(wd, "3_babette", "clock_one_strict_mods.xml")

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



