#!/usr/bin/env Rscript
###########################################################################
# Project: Orania Phylogeny MT
# Script: visualize_covmaps_stats.R
# --- Action: This script computes a "coverage map" for each supercontig. 
# ----------- A coverage map represents how many reads support each
# ----------- position over the supercontig. Maps reads to supercontigs.
# --- Input: supercontigs from script "hybpiper_script.sh"
# --- Output: coverage map for each supercontig (per gene & tagecies)
# Author: Maya Schroedl
# Date: 11/2019
###########################################################################

rm(list=ls())

# Libraries ---------------------------------------------------------------
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
theme_set(theme_classic())
if (!require('rlist')) install.packages('rlist'); library('rlist')
if (!require('plyr')) install.packages('plyr'); library('plyr') #important to load plyr first to not overwrite dplyr
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')


# Basic data --------------------------------------------------------------
wd=file.path(getwd(),"1_phylo_reconstruction")

#read namelist; = list of TAGs
namelist=as.character(read.table(file.path(wd,"namelist.txt"))$V1)

#read genelist; = list of gene names
genelist=as.character(read.table(file.path(wd,"genelist.txt"))$V1)


# Example -----------------------------------------------------------------
tag="TAG-29"
gene="E1007"
t = 2 #threshold for number of reads


# Get coverage maps -------------------------------------------------------
# This function reads the coverage maps which were produced by a part of the script "1.1.0_coverage_maps.sh".

# The function adds also a boolean according to the threshold. 
#over_t==T means that this position has a coverage depth higher than the threshold.

# Each line corresponds to a base pair [bp] of a specific gene for a specific TAG.

coverage=function(tag){
  file=file.path(wd,"1.1_coverage_maps",paste(tag,"_supercontigs_sorted.coverage",sep=""))#the coverage file corresponding to the TAG
  
  if (file.exists(file) & file.info(file)$size != 0) #if file exists and is not empty
    {
  cover=read.table(file, sep="\t", header=F) #read coverage map
  colnames(cover)=c(V1="gene", V2="bp", V3="depth") # renames the header
  cover$tag=rep(tag, length(cover$gene))#add column with the tag name
  cover$over_t=cover$depth>t #add column with whether position have a depth over threshold (TRUE) or not (FALSE)
  return(cover)
    }
  }

cov.ex=coverage(tag)#coverage map for example

#### Add coverage maps of TAGs to one df ####
cover_all=lapply(namelist,coverage) %>% bind_rows() #ignore warnings

write.table(cover_all,file.path(wd,"1.1_coverage_maps","all_covmaps.txt"),row.names = F,quote = F)

# Stats -------------------------------------------------------------------
#---- Simple stat summary ----

cover_all=read.table(file.path(wd,"1.1_coverage_maps","all_covmaps.txt"),h=T)

mean(cover_all$depth)
sd(cover_all$depth)
min(cover_all$depth)
max(cover_all$depth)

#density plot for all basepairs
plot(density(cover_all$depth[which(cover_all$depth<50)]))
abline(v=t,col="red")

#### proportion of genes having a base with a depth <= threshold
length(unique(cover_all$gene[which(cover_all$over_t==FALSE)]))/length(genelist)

#### proportion of contigs having a base with a depth <= threshold
cov.contig=cover_all %>%
  group_by(gene,tag)%>% #for each gene/tag combination:
  summarize(have.tF=(FALSE %in% over_t)) #have.tF: booleean - TRUE = gene/tag combi has at least one base pair whose depth is under the threshold (over_thr==FALSE)

length(which(cover_all$over_t==FALSE))/length(cover_all$gene)

#### percentage of bases having a depth <= threshold (how many bases do we get rid off)
length(which(cov.contig$have.tF==TRUE))/length(cov.contig$gene)

##### percentage [per gene] of bases having a depth <= threshold
cov_per_gene=cover_all %>%
  group_by(gene)%>% #calculate for each gene
  summarize(perc=length(which(over_t==FALSE))/length(over_t)*100) #perc: percentage of base pairs having a depth under threshold (over_t == FALSE)

hist(cov_per_gene$perc) #histogram of "bad base percentage" over genes

#gene with most bases <=t (over all TAGs)
worst_gene=cov_per_gene$gene[which(max(cov_per_gene$perc)==cov_per_gene$perc)] 
worst_gene


# Position stats ----------------------------------------------------------
# This function determines how many positions have a depth under threshold in the flanking regions, and in the middle respectively.

pos.stat=function(gene,tag){
  cover=cover_all[which(cover_all$tag==tag & cover_all$gene==gene),]
  if (length(cover$gene)==0){return(NA)}#if this combination tag/gene does not exist - then function returns NA
  
  len=length(cover$over_t) #total gene length
  over_t_true=grep("TRUE",cover$over_t) #which positions have a depth over threshold
  flan_l=over_t_true[1]-1
  # flan_l: length of the left flanking region
  #     from 0 to position of the first time the depth is > threshold (over_t==TRUE), indicates the end of the flanking region 1
  
  flan_r=length(cover$gene)-over_t_true[length(over_t_true)]
  #flan_r: length of the right flanking region
  #     from position of the last time the depth is > threshold (over_t_true[length(over_t_true)]) to the last position of the gene length(cov.ex$gene)
  
  flan.len=
    flan_l+flan_r #flan.tot: how many positions have a depth<t in the flanking regions
  #sum of left flanking region + right flanking region lengths
  mid.len=
    length(grep("FALSE",cover$over_t))-flan_l-flan_r #mid.len: how many positions have a depth<t in the middle regions
  #length of how many postitions have a depth < threshold (length(grep("FALSE",cov.ex$over_t))), but the lengths of the flanking regions need to be substracted to only have the middle.
  
  perc.flan=round(flan.len/len*100,3) #percentage of positions having a depth < threshold in the flanking regions
  
  perc.mid=round(mid.len/len*100,3) #percentage of positions having a depth < threshold in the middle of the gene
  
  perc.true=round(length(over_t_true)/len*100,3) #percentage of positions having a depth > threshold in the flanking regions
  
  return(c(perc.flan,perc.mid,perc.true))}

pos.stat(gene,tag)

#Summary of the posititon stats
combi_list=expand.grid(genelist,namelist)#dataframe of combinations of gene/tag
colnames(combi_list)=c("gene","tag")

pos.stat_sum=plyr::mdply(combi_list,pos.stat) #for each gene/tag-combination, calculate pos.stat
colnames(pos.stat_sum)=c("gene","tag","perc.flan","perc.mid","perc.good")

#doc with stats for each gene/tag-combination
setwd(file.path(wd,"1.1_coverage_maps"))
write.table(pos.stat_sum,paste("supercontigs_covmap_pos_stats_",t,".txt",sep=""),row.names = F)



# Coverage Maps -----------------------------------------------------------
# This function plots of read coverage (=depth) as a function of position on the gene
# for each gene and tag

covmap_plot=function(gene,tag){
  cover=cover_all[which(cover_all$gene==gene & cover_all$tag==tag),]
  if (length(cover$gene)==0){return(NA)}#if this combination tag/gene does not exist - then function returns NA
    #plot coverage map (depth as a function of position [bp])
    g=ggplot(cover,
             aes(x=bp, y=depth)) +
        geom_ribbon(aes(ymin=pmin(cover$depth,t), ymax=t), 
                  fill="red", col="red", alpha=0.5) + #all values lower than threshold t are coloured in red
        geom_ribbon(aes(ymin=t, ymax=pmax(cover$depth,t)), 
                  fill="green", col="green", alpha=0.5) + #all values higher than threshold t are coloured in green
        geom_line(aes(y=t),col="red")+ggtitle(paste(gene, tag, "t =", t)) + #abline with threshold t
        ylab("depth (reads)")
  return(g)}

#example of coverage map plot with worst gene
covmap_plot(worst_gene,cover_all$tag[which(cover_all$gene==worst_gene)[1]])

#plot all plots for all tags for one gene at the same time
#to have a better overall understanding of the read coverage for this specific gene
covmap_plot_alltags=function(gene){
  taglist=unique(cover_all[which(cover_all$gene==gene),]$tag)#list of tags for this gene
  plots_alltags=sapply(taglist,covmap_plot, gene=gene,simplify=F)#make a list of all the plots for this gene
  return(plots_alltags)
}

covmap_plot_alltags(worst_gene)




