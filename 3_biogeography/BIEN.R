library(BIEN)
library(ape)
library(sp)

wd = file.path(getwd(), "3_biogeography")

if (!dir.exists(file.path(wd,"BIEN"))){dir.create(file.path(wd,"BIEN"))}

or_occ = BIEN_occurrence_genus("Orania")
pod_occ = BIEN_occurrence_genus("Podococcus")
scl_occ = BIEN_occurrence_genus("Sclerosperma")

write.csv(rbind(or_occ,pod_occ,scl_occ), file.path(wd, "BIEN","orsclpod_bien.csv"), quote = F)
