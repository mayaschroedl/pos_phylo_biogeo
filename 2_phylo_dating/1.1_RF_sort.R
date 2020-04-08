#select genes with lowest RF

RF_matrix=read.table("D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/2_phylo_dating/1_hashrf/RF_matrix_rooted.txt")
genes=read.table("D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/1_phylo_reconstruction/genelist_7575.txt")$V1

RF_df=data.frame(name=as.character(genes), RF=RF_matrix[2:nrow(RF_matrix),1])

min(RF_df$RF)
hist(RF_df$RF)

#which ones to take?!

RF_df[order(RF_df$RF),]

write.table(sel.genes[2:length(sel.genes)],"D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/2_hyb_align_reconstruction/AU_server/supercontigs/6_hashrf/or,orscldyp,sclpod_RF_selgenes.txt",quote = F,row.names = F,col.names = F)

#maybe only the ones RF = 0 or = 1?

rownames(matrix)[which(matrix$astral==4)] #=4 (lowest)
