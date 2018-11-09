library(pathview)
library(gage)
library(GEOquery)
library('KEGGREST')
library(readxl)
library("org.Hs.eg.db")
library(clusterProfiler)

source("https://bioconductor.org/biocLite.R")
biocLite("clusterProfiler")


cancerdata <- read.csv("CancerData.csv")
cancerdata <- cancerdata[- grep("///", cancerdata $`Gene.ID.title`),] 
rownames(cancerdata) <- cancerdata$`Gene.ID.title`

can_to_pathview <- as.vector(cancerdata$logFC)
names(can_to_pathview) <- rownames(cancerdata)

species <- 'hsa'
chosen_pathway = 'hsa05200'
pv.out <- pathview(gene.data = can_to_pathview, gene.idtype = "KEGG", 
                   pathway.id = chosen_pathway, species = species, out.suffix = chosen_pathway, keys.align = "y", 
                   kegg.native = T, sign.pos = "bottomleft", same.layer = T)

pv.out <- pathview(gene.data = can_to_pathview, gene.idtype = "KEGG", 
                   pathway.id = chosen_pathway, species = species, out.suffix = chosen_pathway, keys.align = "y", 
                   kegg.native = F, sign.pos = "bottomleft", same.layer = T)


cancerdata_invert <- read.csv("CancerData_remastered.csv")
cancerdata_invert <- cancerdata_invert[- grep("///", cancerdata_invert$`Gene.ID.title`),] 
rownames(cancerdata_invert) <- cancerdata_invert$`Gene.ID.title`
cancerdata_invert_2 <- cancerdata_invert[, -c(1:3)]
cn_1=colnames(cancerdata_invert_2)
lec_con <- grep('LEC.control',cn_1, ignore.case =TRUE)
lec_vir <- grep('LEC.KSHV',cn_1, ignore.case =TRUE)
data("kegg.gs")
cancer_1.kegg.p <- gage(cancerdata_invert_2, gsets = kegg.gs,
                        ref = lec_con, samp = lec_vir)
cancer_1.kegg.2d.p <- gage(cancerdata_invert_2, gsets = kegg.gs,
                           ref = lec_con, samp = lec_vir, same.dir = FALSE)
