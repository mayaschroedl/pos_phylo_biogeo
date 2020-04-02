or=read.table("D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/2_hyb_align_reconstruction/AU_server/1_hybpiper/or/stats_or.txt",h=T)
sclpod=read.table("D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/2_hyb_align_reconstruction/AU_server/1_hybpiper/sclpod/stats_sclpod.txt",h=T)
orscldyp=read.table("D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/2_hyb_align_reconstruction/AU_server/1_hybpiper/orscldyp/stats_orscldyp.txt",h=T)

all=rbind(or,sclpod,orscldyp)
all

plot(all$ReadsMapped,all$GenesMapped)
