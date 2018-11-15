library(pathview)
library(gage)
library(GEOquery)
library('KEGGREST')
library(readxl)
library("org.Hs.eg.db")
library(clusterProfiler)

can1 <- read.csv("CancerData.csv")
can1 <- can1[- grep("///", can1 $`Gene.ID.title`),] 
rownames(can1) <- can1$`Gene.ID.title`

can_to_pathview <- as.vector(can1$logFC)
names(can_to_pathview) <- rownames(can1)
species <- 'hsa'
chosen_pathway = 'hsa05200'
#png version
pv.out <- pathview(gene.data = can_to_pathview, gene.idtype = "KEGG", 
                   pathway.id = chosen_pathway, species = species, out.suffix = chosen_pathway, keys.align = "y", 
                   kegg.native = T, sign.pos = "bottomleft", same.layer = T)
#pdf version (broken)
pv.out <- pathview(gene.data = can_to_pathview, gene.idtype = "KEGG", 
                   pathway.id = chosen_pathway, species = species, out.suffix = chosen_pathway, keys.align = "y", 
                   kegg.native = F, sign.pos = "bottomleft", same.layer = T)

#Basic Analysis

can2 <- read.csv("CancerData_remastered.csv")
can2 <- can2[- grep("///", can2$`Gene.ID.title`),] 
rownames(can2) <- can2$`Gene.ID.title`
can3 <- can2[, -c(1:3)]
cn=colnames(can3)
lec_con <- grep('LEC.control',cn, ignore.case =TRUE)
lec_vir <- grep('LEC.KSHV',cn, ignore.case =TRUE)


data(kegg.gs)
data(go.gs)
can3.kegg.p <- gage(can3, gsets = kegg.gs, ref = lec_con, samp = lec_vir)
can3.go.p <- gage(can3, gsets = go.gs, ref = lec_con, samp = lec_vir)
can3.kegg.2d.p <- gage(can3, gsets = kegg.gs, ref = lec_con, samp = lec_vir, same.dir = F)

#Result Presentation and Intepretation (p.s: tables are needed to be corrected by playing with 'sep')

write.table(can3.kegg.2d.p$greater, file = "can3.kegg.2d.p.txt", sep = "\t")
write.table(rbind(can3.kegg.p$greater, can3.kegg.p$less), file = "can3.kegg.p.txt", sep = "\t")
can3.kegg.sig<-sigGeneSet(can3.kegg.p, outname="can3.kegg")
can3.kegg.2d.sig<-sigGeneSet(can3.kegg.2d.p, outname="can3_1.kegg")
write.table(can3.kegg.2d.sig$greater, file = "can3.kegg.2d.sig.txt", sep = "\t")
write.table(rbind(can3.kegg.sig$greater, can3.kegg.sig$less), file = "can3.kegg.sig.txt", sep = "\t")

can3.kegg.esg.up <- esset.grp(can3.kegg.p$greater, can3, gsets = kegg.gs, ref = lec_con, samp = lec_vir, test4up = T, output = T, outname = "can3.kegg.up", make.plot = F)
can3.kegg.esg.dn <- esset.grp(can3.kegg.p$less, can3, gsets = kegg.gs, ref = lec_con, samp = lec_vir, test4up = F, output = T, outname = "can3.kegg.dn", make.plot = F)

gs=unique(unlist(kegg.gs[rownames(can3.kegg.p$greater)[1:3]]))
essData=essGene(gs, can3, ref =lec_con, samp =lec_vir)

ref1=1:6
samp1=7:12
#png
for (gs in rownames(can3.kegg.p$greater)[1:3]) {
outname = gsub(" |:|/", "_", substr(gs, 10, 100))
geneData(genes = kegg.gs[[gs]], exprs = essData, ref = ref1,
samp = samp1, outname = outname, txt = T, heatmap = T, Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}
#pdf
for (gs in rownames(can3.kegg.p$greater)[1:3]) {
outname = gsub(" |:|/", "_", substr(gs, 10, 100))
outname = paste(outname, "all", sep=".")
geneData(genes = kegg.gs[[gs]], exprs = can3, ref = lec_con,
samp = lec_vir, outname = outname, txt = T, heatmap = T,
Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}
# the same with gagePrep
can3.d <- can3[ ,lec_vir] - can3[ ,lec_con]

#to specific pathways
path.ids=c("hsa04142 Lysosome", "hsa04670 Leukocyte transendothelial migration", "hsa04630 Jak-STAT signaling pathway", "hsa04972 Pancreatic secretion")
path.ids2 <- substr(path.ids, 1, 8)
pv.out.list <- sapply(path.ids2, function(pid) pathview(gene.data = can3.d[,1:2], pathway.id = pid, species = "hsa"))
pv.out.list_1 <- sapply(path.ids2, function(pid) pathview(gene.data = can3.d[,1:2], pathway.id = pid, species = "hsa", kegg.native=F,sign.pos="bottomleft"))

#to all selected pathways in a batch (not so good, get stopped after 50 warnings counted)
sel <- can3.kegg.p$greater[, "q.val"] < 0.1 & !is.na(can3.kegg.p$greater[, "q.val"])
path.ids <- rownames(can3.kegg.p$greater)[sel]
path.ids2 <- substr(path.ids, 1, 8)
pv.out.list <- sapply(path.ids2, function(pid) pathview(gene.data = can3.d[,1:2], pathway.id = pid, species = "hsa"))
