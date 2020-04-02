
wd="D:/Programme/Google_Drive/Aarhus/Orania/Orania_project/2_Assembly_Alignment/Contigmaps_correction"
setwd(wd)

## reorganize data to genes
splist=as.character(read.table(file.path("or_bwa_to_contigs","splist_thomas.txt"))$V1)
for (tag in splist){
TAG=read.table(file.path("or_bwa_to_contigs",paste0(tag,"_contigs_sorted.coverage")))

for (ex in as.character(unique(TAG$V1))){
select=TAG[which(TAG$V1==ex),]
select$V1=rep(tag,length(select$V1))

write.table(select,file.path("or_bwa_to_contigs_genereorga",paste0(ex,"_contigs_sorted.coverage")),
            row.names = F,
            col.names = F,
            quote=F,
            append=T)  
}}          




