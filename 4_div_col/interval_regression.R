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


# WD ----------------------------------------------------------------------
wd = file.path(getwd(), "4_div_col")
all_bb = read.table(file.path(wd, "data", "all_bb.txt"),h=T,sep="\t")
all_b2 = read.table(file.path(wd, "data", "all_b2.txt"),h=T,sep="\t")

# Seed --------------------------------------------------------------------
set.seed(42)

# Plot --------------------------------------------------------------------
ggplot(all_bb, aes(earliest, mg_species_num_log,colour=as.factor(lineage)))+geom_errorbarh(aes(xmin=earliest,xmax=latest),alpha=0.5,size=2)+ theme_classic()


# Correlation -----------------------------------------------------
cor_num = function(date_sp_df,sample_num=100){
  samp_unif = function(group){
    random_vec = runif(sample_num,date_sp_df$latest[group],date_sp_df$earliest[group])
    return(random_vec)
  }

mat = sapply(1:length(date_sp_df[,1]),samp_unif)
  
  cor_test = function(row){
  est = cor.test(mat[row,], date_sp_df$mg_species_num_log)$estimate
  p = cor.test(mat[row,], date_sp_df$mg_species_num_log)$p.value
  return(c(est,p))}
  
length(which(sapply(1:sample_num, cor_test)[2,]<0.05))/sample_num
}

cor_num(all_b2,1000)
  