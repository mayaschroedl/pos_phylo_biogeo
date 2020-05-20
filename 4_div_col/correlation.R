#!/usr/bin/env Rscript
###########################################################################
# Project: Div_col_palm_Madag
# Script: 
# --- Action: 
# --- Input: 
# --- Output:
# Author: Maya Schroedl
###########################################################################

rm(list=ls())

# Libraries ---------------------------------------------------------------
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2');theme_set(theme_classic())
if (!require("RColorBrewer")) install.packages("RColorBrewer"); library("RColorBrewer")


# Seed --------------------------------------------------------------------

set.seed(42)

# WD + Data ----------------------------------------------------------------------
wd = file.path(getwd(), "4_div_col")

all = read.table(file.path(wd, "data", "all.txt"),h=T,sep="\t") # all sources
selected = read.table(file.path(wd, "data", "selected.txt"),h=T,sep="\t") # selected sources

# Colours
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

custom.col = rev(gg_color_hue(9))
custom.col[5]= "#95F985"
custom.col[7] = "gold"

all$col = custom.col[1:length(all$lineage)]
selected$col = custom.col[1:length(selected$lineage)]

### All; both Borassus in one group:
A_B1 = all[-which(all$lineage == "Borassus_both"),]

### All; each Borassus in one group
A_B2 = all[-which(all$lineage %in% c("Borassus_madagascariensis","Borassus_aethiopium_Madagascar")),]

### Selected; both Borassus in one group:
S_B1 = all[-which(selected$lineage == "Borassus_both"),]


### Selected; each Borassus in one group
S_B2 = all[-which(selected$lineage %in% c("Borassus_madagascariensis","Borassus_aethiopium_Madagascar")),]

# Plot --------------------------------------------------------------------
div_time_plot=function(dataset){
  title=deparse(substitute(dataset)) #get input variable name for plot title
  g = ggplot(dataset, aes(earliest, mg_species_num_log)) +
    geom_errorbarh(aes(xmin=earliest,xmax=latest, colour= col,alpha=0.6),size=2)+
    scale_colour_identity()    + 
    theme_classic() +
    theme(text = element_text(size=20))+
    theme(legend.position = "none") +
    ggtitle(title) +
    xlab("\nArrival time [mya]") + 
    ylab("ln (Malagasy species richness)\n")+
    theme(axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20))+
    scale_x_continuous(breaks=seq(0, 70, 10))+
    scale_y_continuous(breaks=seq(0, 5, 1))
return(g)
}

# Correlation -----------------------------------------------------
cor_sampling = function(dataset,rep_nb=100000){
  # We want to see how "probable" a correlation is according to these data
  # Therefore, we sample randomly an arrival time for each group on its arrival time interval (when group possibly arrived to Madagascar). For each series (one sample per group), we test the correlation between the log(Malagasy species richness) and the arrival time. This is repeated $sample_num times and then a proportion of significant correlations and significant positive (as expected) correlations is calculated over all repetitions.
  rep_num = rep_nb
  # Seed 
  

  samp_unif = function(group, rep_nb = rep_num){
    # Sample randomly for the selected group an arrival time on a uniform distribution over the interval; repeat this for $sample_num times
 
    # We suppose a uniform distribution of arrival time over the interval for each group, because no probability distribution is available. We cannot really do any better. 

    random_vec = runif(rep_nb,dataset$latest[group],dataset$earliest[group]) # random sample on uniform distribution
    
    return(random_vec) # return a vector of samples for this group
  }

  # get a vector of samples for each group of species
  sample_mat = sapply(1:length(dataset[,1]),samp_unif)

  # test correlation between the arrival time and log (Malagasy species richness) for each sampling repetition (row of sample_mat)
  cor_test = function(row){
  est = cor.test(sample_mat[row,], dataset$mg_species_num_log)$estimate
  p = cor.test(sample_mat[row,], dataset$mg_species_num_log)$p.value
  return(c(est,p)) #return correlation estimate and p-value
  }
  
# get correlations for each sampling repetition (row of sample_mat)
all_cor = sapply(1:rep_nb, cor_test)

#transform to dataframe
all_cor_df = as.data.frame(t((all_cor)))
colnames(all_cor_df) = c("rho", "p")

# STATS
sign = length(which(all_cor_df$p<0.05))/rep_nb # percentage of siginficant correlations
sign_pos = length(which(all_cor_df$p<0.05 && all_cor_df$rho>0))/rep_nb # percentage of siginficant positive correlations

return(list(sign,sign_pos,all_cor_df))
}


# Application -------------------------------------------------------------

# dataset A_B1: times calculated over all sources + one colonization event for borassus
div_time_plot(A_B1)
A_B1_cor = cor_sampling(A_B1)
A_B1_cor[c(1,2)] # percentage of significant correlations; percentage of significant positive correlations

# dataset A_B2: times calculated over all sources + two colonization events for borassus
div_time_plot(A_B2)
A_B2_cor = cor_sampling(A_B2)
A_B2_cor[c(1,2)] # percentage of significant correlations; percentage of significant positive correlations

# dataset S_B1: times calculated over selected sources (Bellot in prep. instead of Baker & Couvreur (2013a)) + one colonization event for borassus
div_time_plot(S_B1)
S_B1_cor = cor_sampling(S_B1)
S_B1_cor[c(1,2)] # percentage of significant correlations; percentage of significant positive correlations

# dataset S_B2: times calculated over selected sources (Bellot in prep. instead of Baker & Couvreur (2013a)) + one colonization event for borassus
div_time_plot(S_B2)
S_B2_cor = cor_sampling(S_B2)
S_B2_cor[c(1,2)] # percentage of significant correlations; percentage of significant positive correlations


###### EXPORT PLOT dataset S_B1

plot.dir = file.path(getwd(), "plots", "div_col")
if(!dir.exists(plot.dir)){dir.create(plot.dir)}

pdf(file.path(plot.dir,"S_B1.pdf"))
div_time_plot(S_B1)
dev.off()

##### MAKE HISTOGRAM of slopes dataset S_B1
##### 
pdf(file.path(plot.dir,"S_B1_hist.pdf"))
ggplot(S_B1_cor[3][[1]])+
  geom_histogram(aes(S_B1_cor[3][[1]]$rho))+
  xlim(c(-1,1))+
  xlab("\nPearsons' rho")+
  ylab("Frequency\n")+
  theme(text = element_text(size=20))+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20))
  
dev.off()
