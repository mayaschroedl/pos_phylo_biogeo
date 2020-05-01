#!/usr/bin/env Rscript
###########################################################################
# Project: Div_col_palm_Madag
# Script: interval_regression.R
# --- Action: 
# --- Input: 
# --- Output:
# Author: Maya Schroedl
# Date: 11/2019
###########################################################################

rm(list=ls())

# Libraries ---------------------------------------------------------------
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require("RColorBrewer")) install.packages("RColorBrewer"); library("RColorBrewer")


# WD + Data ----------------------------------------------------------------------
wd = file.path(getwd(), "4_div_col")

all = read.table(file.path(wd, "data", "all.txt"),h=T,sep="\t") # all sources
best = read.table(file.path(wd, "data", "selected.txt"),h=T,sep="\t") # selected sources

# Colours
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

custom.col = rev(gg_color_hue(9))
custom.col[5]= "#95F985"
custom.col[7] = "gold"

all$col = custom.col[1:length(all$lineage)]
best$col = custom.col[1:length(best$lineage)]

### All; both Borassus in one group:
A1 = all[-which(all$lineage %in% c("Borassus_madagascariensis","Borassus_aethiopium_Madagascar")),]

### All; each Borassus in one group
A2 = all[-which(all$lineage == "Borassus_both"),]

### All; both Borassus in one group:
S1 = all[-which(best$lineage %in% c("Borassus_madagascariensis","Borassus_aethiopium_Madagascar")),]

### All; each Borassus in one group
S2 = all[-which(best$lineage == "Borassus_both"),]

# Plot --------------------------------------------------------------------
div_time_plot=function(scenario){
  title=deparse(substitute(scenario)) #get input variable name for plot title
  g = ggplot(scenario, aes(earliest, mg_species_num_log)) +
    geom_errorbarh(aes(xmin=earliest,xmax=latest, colour= col,alpha=0.6),size=2)+
    scale_colour_identity()    + 
    theme_classic() +
    theme(text = element_text(size=20))+
    theme(legend.position = "none") +
    ggtitle(title) +
    xlab("\nMalagasy arrival time [mya]") + 
    ylab("log (Malagasy species richness)\n")+
    theme(axis.text.x = element_text(size=20),
          axis.text.y = element_text(size=20))+
    scale_x_continuous(breaks=seq(0, 70, 10))+
    scale_y_continuous(breaks=seq(0, 5, 1))
return(g)
}

# Correlation -----------------------------------------------------
cor_sampling = function(scenario,sample_num=1000){
  # We want to see how "probable" a correlation is according to these data
  # Therefore, we sample randomly an arrival time for each group on its arrival time interval (when group possibly arrived to Madagascar). For each series (one sample per group), we test the correlation between the log(Malagasy species richness) and the arrival time. This is repeated $sample_num times and then a proportion of significant correlations and significant positive (as expected) correlations is calculated over all repetitions.

  # Seed 
  set.seed(42)

  samp_unif = function(group, sample_num = sample_num){
    # Sample randomly for the selected group an arrival time on a uniform distribution; repeat this for $sample_num times
 
    # We suppose a uniform distribution of arrival time over the interval for each group, because no probability distribution is available. We cannot really do any better. 

    random_vec = runif(sample_num,scenario$latest[group],scenario$earliest[group]) # random sample on uniform distribution
    
    return(random_vec) # return a vector of samples for this group
  }

# get a vector of samples for each group of species
sample_mat = sapply(1:length(scenario[,1]),samp_unif)

# test correlation between the arrival time and log (Malagasy species richness) for each sampling repetition (row of sample_mat)
cor_test = function(row){
  est = cor.test(sample_mat[row,], scenario$mg_species_num_log)$estimate
  p = cor.test(sample_mat[row,], scenario$mg_species_num_log)$p.value
  return(c(est,p)) #return correlation estimate and p-value
  }
  
# get correlations for each sampling repetition (row of sample_mat)
all_cor = sapply(1:sample_num, cor_test)

# STATS
sign = length(which(correlations[2,]<0.05))/sample_num # percentage of siginficant correlations
sign_pos = length(which(correlations[2,]<0.05 && correlations[1,]>0))/sample_num # percentage of siginficant positive correlations

return(c(sign,sign_pos))
}


# Application -------------------------------------------------------------

# Scenario A1: times calculated over all sources + one colonization event for borassus
div_time_plot(A1)
cor_num(A1)

# Scenario A2: times calculated over all sources + two colonization events for borassus
div_time_plot(A2)
cor_num(A2)

# Scenario S1: times calculated over selected sources (Bellot in prep. instead of Baker & Couvreur (2013a)) + one colonization event for borassus
div_time_plot(S1)
cor_num(S1)

# Scenario S2: times calculated over selected sources (Bellot in prep. instead of Baker & Couvreur (2013a)) + one colonization event for borassus
div_time_plot(S2)
cor_num(S2)

###### EXPORT PLOT SCENARIO S1

plot.dir = file.path(getwd(), "plots", "div_col")
if(!dir.exists(plot.dir)){dir.create(plot.dir)}

pdf(file.path(plot.dir,"S1.pdf"))
div_time_plot(S1)
dev.off()

   