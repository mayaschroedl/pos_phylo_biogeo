# Supplementary material of the Master's thesis "Colonization and diversification history of Madagascan palms with new phylogenomic evidence from the genus *Orania* (Arecaceae)" by Maya Schrödl


## PART B: Testing the 'Time-for-Speciation hypothesis' in Madagascan palms

#### **`Fig_S1_S-B2.jpg`**

**Figure S1:** Representation of the natural logarithm of Madagascan species richness
as a function of arrival time to Madagascar, and histogram representing
the frequencies of different standardized slopes for the correlations ln(species richness) ~ arrival time. Same as Figure 1 in main text (cf. thesis), but for the dataset S-B2 (selected sources & two Madagascan colonization events for Borassus). For a detailed description see the figure caption.

#### **`Fig_S2_A-B1.jpg`**

**Figure S2:** Representation of the natural logarithm of Madagascan species richness
as a function of arrival time to Madagascar, and histogram representing
the frequencies of different standardized slopes for the correlations ln(species richness) ~ arrival time. Same as Figure 1 in main text (cf. thesis), but for the dataset A-B1 (all sources & only one Madagascan colonization event for Borassus). For a detailed description see the figure caption.

#### **`Fig_S3_A-B2.jpg`**

**Figure S3:** Representation of the natural logarithm of Madagascan species richness
as a function of arrival time to Madagascar, and histogram representing
the frequencies of different standardized slopes for the correlations ln(species richness) ~ arrival time. Same as Figure 1 in main text (cf. thesis), but for the dataset A-B2 (all sources & two Madagascan colonization events for Borassus). For a detailed description see the figure caption.


## PART C: Phylogenomics & Biogeography of the genus *Orania*

### DATA

#### **`Data_S1_Alignments.zip`**

**Data S1: Alignments.** 63 fasta files of alignements generated with MAFFT,
trimmed with GBlocks and manually cleaned. Each fasta file corresponds to the
alignment of one gene. The faster headers correspond to the individual names in
"TAG notation". A text file is given for the conversion of the TAG notation to the
species name for each individual.

#### **`Data_S2_Gene_trees.zip`**

**Data S2: Gene trees.** 63 newick files of gene trees generated with RAxMLng.
Each newick file corresponds to the gene tree of one gene. The faster headers
correspond to the individual names in "TAG notation". A text file is given for the
conversion of the TAG notation to the species name for each individual. A pdf file
is given showing all visualized gene trees.


### TABLES

#### **`Tab_S1_HybPiper_stats.pdf`**

**Table S1: HybPiper statistics.** Summary of target enrichment and gene recovery
efficiency for all 20 samples, output by HybPiper. The columns correspond
respectively to: Individual name, Number of reads, Number of reads on target, Percent
reads on target, Number of genes with reads, Number of genes with contigs,
Number of genes with sequences, Number of genes with sequences > 25% of the
target length, Number of genes with sequences > 50% of the target length, Number
of genes with sequences > 75% of the target length, Number of genes with sequences
$\lt$ 150% of the target length, Number of genes with paralog warnings, as described
in the HybPiper Wiki (https://github.com/mossmatters/HybPiper/wiki).

#### **`Tab_S2_Biogeo_model_select_no_constr.pdf`**

**Table S2: Biogeographic model fit statistics from BioGeoBEARS for the
estimation without node constraint.** The first part of the model name corresponds
to the range evolution model: dec (DEC), dl (DIVALIKE), b (BAYAREALIKE);
and the second part of the model name to the dispersal probability models:
no dispersal probability constraints (M0), dispersal probability constraints as in
Baker & Couvreur (2013a) (B: MB), dispersal probability constraints as in Federman
& al. (2015) (F: MF ). LnL: Log likelihood; numparams: number of parameters
in the model; d: estimated dispersal rate; e: estimated extinction (=extirpation)
rate; AICc: AICc score; AICc_wt: AICc weight.

#### **`Tab_S3_Biogeo_model_select.pdf`**

**Table S3: Biogeographic model fit statistics from BioGeoBEARS for the
estimation with node constraint.** The node between the outgroup (core Arecoids)
and the POS clade was constraint to the range Southeast Asia/Eurasia, based
on Baker & Couvreur (2013a). The first part of the model name corresponds to the
range evolution model: dec (DEC), dl (DIVALIKE), b (BAYAREALIKE); and the
second part of the model name to the dispersal probability models: no dispersal
probability constraints (M0), dispersal probability constraints as in Baker & Couvreur
(2013a) (B: MB), dispersal probability constraints as in Federman & al. (2015)
(F: MF ). LnL: Log likelihood; numparams: number of parameters in the model; d:
estimated dispersal rate; e: estimated extinction (=extirpation) rate; AICc: AICc
score; AICc_wt: AICc weight.

### FIGURES

#### **`Fig_S1_HybPiper_Heatmap.pdf`**

**Figure S1: HybPiper heatmap. Visualization of recovery efficiency, figure created
by HybPiper.** Each row shows a sample, and each column is a gene. The
amount of shading in each box corresponds to the length of the gene recovered for
that sample, relative to the length of the target coding sequence (cf. HybPiper Wiki
https://github.com/mossmatters/HybPiper/wiki).

#### **`Fig_S2_PhyParts_as_pie_charts.pdf`**

**Figure S2: ASTRAL cladogram of individuals with PhyParts as pie charts
with node numbers and gene tree support values.** The pie charts are based
on results from PhyParts and show the proportion of gene trees that are concordant
with the species tree topology at the corresponding node (green), that support the
main alternative topology (orange), that support other alternative topologies (red),
and the proportion of gene trees that are uninformative for that node (black). The
values above and below the branches represent respectively the number of concordant
and conflicting gene trees. The numbers on grey represent the node numbers.


#### **`Fig_S3_Phyparts_top_alternatives.pdf`**

**Figures S3: Each figure corresponds to one node: Top alternative topology
at corresponding node generated with PhyParts.** Node numbers correspond
to the node numbers on the phylogeny in Figure S2. The number of agreeing
gene trees with this alternative topology corresponds to the orange part of the pie
chart in Figure 1b and Figure S2.


#### **`Fig_S4_Clado_monoconstr_lauterb.pdf`**

**Figure S4: Cladogram of individuals when constraining the individuals
of O. lauterbachiana to be monophyletic** with the ASTRAL option "-j".
Support values are given as local posterior probabilities on top of each branch. The
red box indicates the O. lauterbachiana clade.

#### **`Fig_S5_Clado_monoconstr_palindan.pdf`**

**Figure S5: Cladogram of individuals when constraining the individuals of
O. palindan to be monophyletic** with the ASTRAL option "-j". Support
values are given as local posterior probabilities on top of each branch. The red box
indicates the O. palindan clade.

#### **`Fig_S6_chrono_23_genes.pdf`**

**Figure S6: Maximum clade credibility chronogram of the POS clade using
BEAST based on 23 selected genes.** These genes were selected based on the
least topology conflict of the gene trees with the species tree. The blue bars indicate
the 95% highest posterior density intervals. The number above each node indicates
the node age and the number below the node indicates the lowest and highest 95%
highest posterior density.


#### **`Fig_S7_chrono_63_genes.pdf`**

**Figure S7: Maximum clade credibility chronogram of the POS clade using
BEAST based on all 63 genes.** The blue bars indicate the 95% highest
posterior density intervals. The number above each node indicates the node age and
the number below the node indicates the lowest and highest 95% highest posterior
density.

#### **`Fig_S8_Biogeo_no_node_constr.pdf`**

**Figure S8: Biogeographic ancestral range estimates using BioGeoBears
without biogeographic node constraint** on the chronogram of Figure 2 without
the outgroup (i.e. Dypsis mananjarensis). At the top left is represented which colour
corresponds to which geographical area: grey (Africa); orange (Madagascar); blue
(Southeast Asia = Malesia + Papuasia). The colours of the tip labels correspond
to the current range of each species. The parts within the pie charts correspond
to the range probabilities at the corresponding node. The striped parts indicate
the presence in multiple areas. The different sizes of the pie charts are only for
visualization purposes. Pie charts on the stem nodes of the tips are not represented;
the range probability is 1 for all corresponding tip ranges at these nodes. The nodes
of interest are labelled with Roman numerals and grey boxes correspond to the
numbers added to those in the cladogram (Figure 1). The Pliocene is abbreviated
as "Plio." and the Pleistocene as "Pleisto.".

### TEXT

#### **`Text_S1_Biogeo_dispersal_matrices.pdf`**

**Text S1: Dispersal probability matrices for ancestral range estimation.**
Text explaining how dispersal probability matrices were obtained for the two models:
MB (matrices adapted from Baker & Couvreur (2013a)) and MF (matrices adapted
from Federman & al. (2015)); and the resulting dispersal probability matrices for
each model.