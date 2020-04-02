#plot trees
install.packages("ape")
install.packages("phangorn")
install.packages("phytools")
install.packages("geiger")

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}
BiocManager::install("ggtree")

install.packages("remotes")
remotes::install_github("thackl/thacklr")


library(ape)
library(phangorn)
library(phytools)
library(thacklr)


faurby.tree<-read.nexus(file="faurby2016/Phylogeny_NoCon_Checklist_2.nex")

extr=function(n){
  node<-fastMRCA(faurby.tree[[n]],"Orania_ravaka",
                 "Podococcus_barteri")
  faurby.tree.extr<-extract.clade(faurby.tree[[n]],node)
  return (faurby.tree.extr)
}

result=lapply(1:length(faurby.tree), extr)

#' Convert a List into multiPhylo Object
as.multiPhylo.list <- function(x, ...){
  if(!all(sapply(x, function(y) inherits(y, "phylo"))))
    stop("Need a list of phylo objects")
  thacklr::set_class(x, "multiPhylo")
}

resm=as.multiPhylo.list(result)

plotTree(resm[[42]],ftype="i")

faurby.tree.cons=consensusNet(resm)
