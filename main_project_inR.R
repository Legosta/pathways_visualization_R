

library(pathview)
library(gage)
library(GEOquery)
library('KEGGREST')
library(readxl)
library("org.Hs.eg.db")
library(clusterProfiler)

#========================================
#--------------PATHWAY SECTION
#========================================

setwd("/home/legosta/visualization_R")
data_16357 <- read.csv("CancerData.csv")
data_16357_1 <- data_16357[- grep("///", data_16357 $`Gene.ID.title`),]
rownames(data_16357_1) <- data_16357_1$`Gene.ID.title`

data_16357_to_pathview <- as.vector(data_16357_1$logFC)
names(data_16357_to_pathview) <- rownames(data_16357_1)
species <- 'hsa'
chosen_pathway = 'hsa05200'

# - png version
pv.out <- pathview(gene.data = data_16357_to_pathview, gene.idtype = "KEGG", 
                   pathway.id = chosen_pathway, species = species, out.suffix = chosen_pathway, keys.align = "y", 
                   kegg.native = T, sign.pos = "bottomleft", same.layer = T)
#- pdf version (broken, "missing interaction" mistake)
pv.out <- pathview(gene.data = data_16357_to_pathview, gene.idtype = "KEGG", 
                   pathway.id = chosen_pathway, species = species, out.suffix = chosen_pathway, keys.align = "y", 
                   kegg.native = F, sign.pos = "bottomleft", same.layer = T)

#========================================
#--------------GAGE SECTION
#========================================

# - proccesing

setwd("/home/legosta/visualization_R/rubbish")
data_16357_inv <- read.csv("CancerData_remastered.csv")
data_16357_inv <- data_16357_inv[- grep("///", data_16357_inv$`Gene.ID.title`),] 
rownames(data_16357_inv) <- data_16357_inv$`Gene.ID.title`
data_16357_inv_1 <- data_16357_inv[, -c(1:3)]
cn <- colnames(data_16357_inv_1)
lec_con <- grep('LEC.control', cn, ignore.case =TRUE)
lec_vir <- grep('LEC.KSHV', cn, ignore.case =TRUE)

# - gage

data(kegg.gs)
data(go.gs)
inv16357.kegg.p <- gage(data_16357_inv_1, gsets = kegg.gs, ref = lec_con, samp = lec_vir)
inv16357.go.p <- gage(data_16357_inv_1, gsets = go.gs, ref = lec_con, samp = lec_vir)
inv16357.kegg.2d.p <- gage(data_16357_inv_1, gsets = kegg.gs, ref = lec_con, samp = lec_vir, same.dir = F)

# - запись данных, содержащих значимые результаты (q.value < 0.05)

inv16357.kegg.sig<-sigGeneSet(inv16357.kegg.p, outname="can3.kegg")
inv16357.kegg.2d.sig<-sigGeneSet(inv16357.kegg.2d.p, outname="inv16357_2d.kegg")
five_sig_greater_hsa <- rownames(inv16357.kegg.sig$greater)[1:5]
five_sig_less_hsa <- rownames(inv16357.kegg.sig$less)[1:5]
five_sig_2d_greater_hsa <- rownames(inv16357.kegg.2d.sig$greater)[1:5]


write.csv(inv16357.kegg.sig$greater, file = "SIG.inv16357.kegg.p.greater.csv")
Q.SIG.inv16357.kegg.p.greater <- read.csv("SIG.inv16357.kegg.p.greater.csv")
Q.SIG.inv16357.kegg.p.greater <- subset(Q.SIG.inv16357.kegg.p.greater, Q.SIG.inv16357.kegg.p.greater$q.val < 0.05)
write.csv(Q.SIG.inv16357.kegg.p.greater, file = "Q.SIG.inv16357.kegg.p.greater.csv")

#___________________________

write.csv(inv16357.kegg.sig$less, file = "SIG.inv16357.kegg.p.less.csv")
Q.SIG.inv16357.kegg.p.less <- read.csv("SIG.inv16357.kegg.p.less.csv")
Q.SIG.inv16357.kegg.p.less <- subset(Q.SIG.inv16357.kegg.p.less, Q.SIG.inv16357.kegg.p.less$q.val < 0.05)
write.csv(Q.SIG.inv16357.kegg.p.less, file = "Q.SIG.inv16357.kegg.p.less.csv")

#___________________________

write.csv(inv16357.kegg.2d.sig$greater, file = "SIG.inv16357.kegg.2d.greater.csv")
Q.SIG.inv16357.kegg.2d.greater <- read.csv("SIG.inv16357.kegg.2d.greater.csv")
Q.SIG.inv16357.kegg.2d.greater <- subset(Q.SIG.inv16357.kegg.2d.greater, Q.SIG.inv16357.kegg.2d.greater$q.val < 0.05)
write.csv(Q.SIG.inv16357.kegg.2d.greater, file = "Q.SIG.inv16357.kegg.2d.greater.csv")

#___txt__________________________

inv16357.kegg.esg.up <- esset.grp(inv16357.kegg.p$greater, data_16357_inv_1, gsets = kegg.gs, ref = lec_con, samp = lec_vir, test4up = T, output = T, outname = "inv16357.kegg.up", make.plot = F)
inv16357.kegg.esg.dn <- esset.grp(inv16357.kegg.p$less, data_16357_inv_1, gsets = kegg.gs, ref = lec_con, samp = lec_vir, test4up = F, output = T, outname = "inv16357.kegg.dn", make.plot = F)

#============gage heatmaps===============

ref1=1:6
samp1=7:12

# -- $greater_spec

gs=unique(unlist(kegg.gs[rownames(inv16357.kegg.p$greater)[1:5]]))
essData=essGene(gs, data_16357_inv_1, ref =lec_con, samp =lec_vir)

for (gs in rownames(inv16357.kegg.p$greater)[1:5]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  geneData(genes = kegg.gs[[gs]], exprs = essData, ref = ref1,
           samp = samp1, outname = outname, txt = T, heatmap = T, Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}

# -- $greater_all

for (gs in rownames(inv16357.kegg.p$greater)[1:5]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  outname = paste(outname, "all", sep=".")
  geneData(genes = kegg.gs[[gs]], exprs = data_16357_inv_1, ref = lec_con,
           samp = lec_vir, outname = outname, txt = T, heatmap = T,
           Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}

# -- $less_spec (Mistake! The number of genes found in exprs is 0 or 1, no need to proceed)

gs=unique(unlist(kegg.gs[rownames(inv16357.kegg.p$less)[1:5]]))
essData=essGene(gs, data_16357_inv_1, ref =lec_con, samp =lec_vir)

for (gs in rownames(inv16357.kegg.p$less)[1:5]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  geneData(genes = kegg.gs[[gs]], exprs = essData, ref = ref1,
           samp = samp1, outname = outname, txt = T, heatmap = T, Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}

# -- $less_all(It works!)

for (gs in rownames(inv16357.kegg.p$less)[1:5]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  outname = paste(outname, "all", sep=".")
  geneData(genes = kegg.gs[[gs]], exprs = data_16357_inv_1, ref = lec_con,
           samp = lec_vir, outname = outname, txt = T, heatmap = T,
           Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}


# -- $2d.greater_spec

gs=unique(unlist(kegg.gs[rownames(inv16357.kegg.2d.p$greater)[1:5]]))
essData=essGene(gs, data_16357_inv_1, ref =lec_con, samp =lec_vir)

for (gs in rownames(inv16357.kegg.2d.p$greater)[1:5]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  geneData(genes = kegg.gs[[gs]], exprs = essData, ref = ref1,
           samp = samp1, outname = outname, txt = T, heatmap = T, Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}

# -- $2d.greater_all

for (gs in rownames(inv16357.kegg.2d.p$greater)[1:3]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  outname = paste(outname, "all", sep=".")
  geneData(genes = kegg.gs[[gs]], exprs = data_16357_inv_1, ref = lec_con,
           samp = lec_vir, outname = outname, txt = T, heatmap = T,
           Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}



#============gage pathways===============

data_16357_inv_1.d <- data_16357_inv_1[ ,lec_vir] - data_16357_inv_1[ ,lec_con] # distinct (six) pathview sections colouring
data_16357_inv_1.d  <- gagePrep(data_16357_inv_1, ref = c(1, 3, 5, 7, 9, 11), samp = c(2, 4, 6, 8, 10, 12), compare="as.group") # uniform pathview sections colouring

# -- for 5 greater significant pathways 

path.ids=five_sig_greater_hsa 
path.ids2 <- substr(path.ids, 1, 8)
pv.out.list_png <- sapply(path.ids2, function(pid) pathview(gene.data = data_16357_inv_1.d, pathway.id = pid, species = "hsa"))
pv.out.list_pdf <- sapply(path.ids2, function(pid) pathview(gene.data = data_16357_inv_1.d, pathway.id = pid, species = "hsa", kegg.native=F,sign.pos="bottomleft"))

# -- for 5 less significant pathways

path.ids=five_sig_less_hsa 
path.ids2 <- substr(path.ids, 1, 8)
pv.out.list_png <- sapply(path.ids2, function(pid) pathview(gene.data = data_16357_inv_1.d, pathway.id = pid, species = "hsa"))
pv.out.list_pdf <- sapply(path.ids2, function(pid) pathview(gene.data = data_16357_inv_1.d, pathway.id = pid, species = "hsa", kegg.native=F,sign.pos="bottomleft"))

# -- for 5 two-direction greater significant pathways 

path.ids=five_sig_2d_greater_hsa
path.ids2 <- substr(path.ids, 1, 8)
pv.out.list <- sapply(path.ids2, function(pid) pathview(gene.data = data_16357_inv_1.d, pathway.id = pid, species = "hsa"))
pv.out.list_1 <- sapply(path.ids2, function(pid) pathview(gene.data = data_16357_inv_1.d, pathway.id = pid, species = "hsa", kegg.native=F,sign.pos="bottomleft"))
