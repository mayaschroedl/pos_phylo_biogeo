install.packages("ape")
install.packages("phangorn")
install.packages("phytools")
install.packages("geiger")

library(ape)
library(phangorn)
library(phytools)
library(geiger)

### tree orsclpoddyp
text.string="(Sclerosperma_mannii:0.18859116980755017,(Sclerosperma_walkeri..mybe:0.22148814541376946,((Dypsis_mananjarensis,(Podoccocus_barteri:0.5519097071358653,Podoccocus_acaulis)1:2.942502398287406)0.48:0.03110072401971957,(((Orania_sylvicola:0.0,Orania_sp.)1:0.6498953102528635,(Orania_timikae,(Orania_deflexa,(Orania_dafonsoroensis,(Orania_palindan:0.0,(Orania_macropetala,(Orania_lauterbachiana:0.0,(Orania_grandiflora,(Orania_tabubilensis,Orania_archboldiana)1:0.25770583617170323)0.98:0.1562448642235)0.49:0.07013277183398067)0.98:0.15304459150265387)1:0.2702654097167512)0.71:0.07282018020462648)1:0.3529123705734845)1:1.1644871784667303)1:1.8859949162056466,((Orania_ravaka:0.2332974836543001,Orania_trispatha)0.81:0.10695886217476851,Orania_longisquama)1:1.7773670585586028)1:2.6130476026770078)1:2.408405801941171)0.99:0.18859116980755017);
"
as.phylo(text.string)
vert.tree<-read.tree(text=text.string)
plot(vert.tree,no.margin=TRUE,edge.width=2)


##really nice:

pdffn = "tree.pdf"
pdf(file=pdffn, width=9, height=12)

tr = read.tree(trfn)
tr
plot(tr)
title("Example Psychotria phylogeny from Ree & Smith (2008)")
axisPhylo() # plots timescale

dev.off()
cmdstr = paste0("open ", pdffn)
system(cmdstr)