#select genes with lowest RF

matrix=read.table("D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/2_hyb_align_reconstruction/AU_server/supercontigs/6_hashrf/or,orscldyp,sclpod_RF_matrix.txt")
genes=read.table("D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/2_hyb_align_reconstruction/AU_server/supercontigs/6_hashrf/or,orscldyp,sclpod_RF_genelist.txt")$V1
colnames(matrix)=genes
rownames(matrix)=genes 
View(matrix)
hist(matrix$astral[2:length(matrix$astral)])

#which ones to take?!

#here:all except last quartile
matrix.ord=matrix[with(matrix,order(matrix$astral)),] #sort by size
sel.genes=rownames(matrix.ord[which(matrix.ord$astral<=quantile(matrix.ord$astral)[4]),])

write.table(sel.genes[2:length(sel.genes)],"D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/2_hyb_align_reconstruction/AU_server/supercontigs/6_hashrf/or,orscldyp,sclpod_RF_selgenes.txt",quote = F,row.names = F,col.names = F)

#maybe only the ones RF = 0 or = 1?

rownames(matrix)[which(matrix$astral==4)] #=4 (lowest)
