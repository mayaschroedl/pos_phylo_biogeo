###########################################################################
# Project: Orania Phylogeny MT
# Script: alignment_script.sh
# --- Action: This script builts a coalescent tree using Astral.
# --- Input: genetrees from raxml; file with species names
# --- Output: coalescent tree
# Author: Maya Schroedl
# Date:11/2019
###########################################################################


WD=system("$PWD")

## reorganize coverage data to genes
splist=as.character(read.table(file.path(WD,"1_hybpiper",ID,paste0("namelist_",ID,".txt"))
for (tag in splist){
TAG=read.table(file.path(WD,"1.1_coverage_maps",ID,paste0(tag,"_contigs_sorted.coverage")))

for (gene in as.character(unique(TAG$V1))){
select=TAG[which(TAG$V1==gene),]
select$V1=rep(tag,length(select$V1))

write.table(select,file.path("or_bwa_to_contigs_genereorga",paste0(gene,"_contigs_sorted.coverage")),
            row.names = F,
            col.names = F,
            quote=F,
            append=T)  
}}          




