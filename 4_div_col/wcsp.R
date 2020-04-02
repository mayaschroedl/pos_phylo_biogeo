#####WCSP#####
rm(list=ls())

wd="D:/ONEDRIVE_AU/OneDrive - Aarhus Universitet/Orania_project/5_div_col/"

dyps=c("Dypsis","Lemurophoenix", "Marojeyja","Masoala")
list="Voanioala"

#list of all palms (WCSP, needs to be updates) 
palms=read.csv(file.path(wd,"WCSP_published_names_19_10_2018_palms2.csv"),sep=";",h=T)
pal.list=as.character(unique(palms$accepted_name[which(palms$genus %in% list & palms$taxon_status_description=="Accepted" & palms$species!="NA")]))


madagascar=read.table(file.path(wd,"Madagascar_palm_WCSP.txt"),h=T,sep=",")
madagascar$name=paste(madagascar$Genus,madagascar$Species)
mad.list=madagascar$name[which(madagascar$Genus %in% list & madagascar$Distribution!="Introduced")]

new_sp=length(which(!(mad.list %in% pal.list))) #new species which are not yet in this version of WCSP, but are online

#how many species for this genus
length(pal.list)+new_sp #overall
length(mad.list) #only on madagascar

#which species do not occur in madagascar
pal.list[which(!(pal.list %in% mad.list))]

