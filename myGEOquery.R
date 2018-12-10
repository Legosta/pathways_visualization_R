#--------for Linix terminal----------
args <- commandArgs(trailingOnly = TRUE)
gse_name = as.character(args[1])
gpl_name <- as.character(args[2])
norm_colmn <- as.integer(args[3])
dis_colmn <- as.integer(args[4])
#------------------------------------

dir.create('~/Data')
setwd('~/Data')

library(GEOquery)
library(dplyr)
library(pathview)
library(gage)
library("org.Hs.eg.db")
library(clusterProfiler)

gse_name <-  "GSE16357"
platform_name <- 'GPL570'

#======GeoQuery========

gse <- getGEO(gse_name, GSEMatrix=FALSE)

gsmplatforms <- lapply(GSMList(gse),function(x) {Meta(x)$platform_id})

gsmlist = Filter(function(gsm) {Meta(gsm)$platform_id==platform_name},GSMList(gse))
gene_id <- Table(GPLList(gse)[[1]])$ENTREZ_GENE_ID
probesets <- Table(GPLList(gse)[[1]])$ID

data.matrix <- do.call('cbind',lapply(gsmlist,function(x) 
{tab <- Table(x)
mymatch <- match(probesets,tab$ID_REF)
return(tab$VALUE[mymatch])
}))
data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})

#========Data processing==========
rownames(data.matrix) <- gene_id
data.matrix <- data.matrix[- grep("///", rownames(data.matrix)),]

write.csv(data.matrix, "RAW_DATA.csv")
express_data_ <- read.csv("RAW_DATA.csv")

express_data_  <-  filter(express_data_, X != '') 
express_data_ <- express_data_ %>%
  group_by(X) %>%
  summarise_all(mean)
rownames(express_data_) <- express_data_$X
express_data <- express_data_[, -c(1)]
express_data <- as.data.frame(express_data)
rownames(express_data) <- rownames(express_data_)
write.csv(express_data, "SEMI_RAW_DATA.csv")

norm_colmn <- 1:6
dis_colmn <- 7:12

data(kegg.gs)
data(go.gs)

exprss_data.kegg.p <- gage(express_data, gsets = kegg.gs, ref = norm_colmn, samp = dis_colmn)
exprss_data.go.p <- gage(express_data, gsets = go.gs, ref = norm_colmn, samp = dis_colmn)
exprss_data.2d.kegg.p <- gage(express_data, gsets = kegg.gs, ref = norm_colmn, samp = dis_colmn, same.dir = F)

exprss_data.kegg.sig<-sigGeneSet(exprss_data.kegg.p, outname="exprss_data.kegg")
exprss_data.kegg.2d.sig<-sigGeneSet(exprss_data.2d.kegg.p, outname="exprss_data_2d.kegg")

five_sig_greater_hsa <- rownames(exprss_data.kegg.sig$greater)[1:5]
five_sig_less_hsa <- rownames(exprss_data.kegg.sig$less)[1:5]
five_sig_2d_greater_hsa <- rownames(exprss_data.kegg.2d.sig$greater)[1:5]

dir.create('~/Data/greater')
setwd('~/Data/greater')

write.csv(exprss_data.kegg.sig$greater, file = "SIG.exprss_data.kegg.p.greater.csv")
Q.SIG.exprss_data.kegg.p.greater <- read.csv("SIG.exprss_data.kegg.p.greater.csv")
Q.SIG.exprss_data.kegg.p.greater <- subset(Q.SIG.exprss_data.kegg.p.greater, Q.SIG.exprss_data.kegg.p.greater$q.val < 0.05)
write.csv(Q.SIG.exprss_data.kegg.p.greater, file = "Q.SIG.exprss_data.kegg.p.greater.csv")

dir.create('~/Data/less')
setwd('~/Data/less')

write.csv(exprss_data.kegg.sig$less, file = "SIG.exprss_data.kegg.p.less.csv")
Q.SIG.exprss_data.kegg.p.less <- read.csv("SIG.exprss_data.kegg.p.less.csv")
Q.SIG.exprss_data.kegg.p.less <- subset(Q.SIG.exprss_data.kegg.p.less, Q.SIG.exprss_data.kegg.p.less$q.val < 0.05)
write.csv(Q.SIG.exprss_data.kegg.p.less, file = "Q.SIG.exprss_data.kegg.p.less.csv")

dir.create('~/Data/2d.greater')
setwd('~/Data/2d.greater')

write.csv(exprss_data.kegg.2d.sig$greater, file = "SIG.exprss_data.kegg.2d.greater.csv")
Q.SIG.exprss_data.kegg.2d.greater <- read.csv("SIG.exprss_data.kegg.2d.greater.csv")
Q.SIG.exprss_data.kegg.2d.greater <- subset(Q.SIG.exprss_data.kegg.2d.greater, Q.SIG.exprss_data.kegg.2d.greater$q.val < 0.05)
write.csv(Q.SIG.exprss_data.kegg.2d.greater, file = "Q.SIG.exprss_data.kegg.2d.greater.csv")

setwd('~/Data')
exprss_data.kegg.esg.up <- esset.grp(exprss_data.kegg.p$greater, express_data, gsets = kegg.gs, ref = norm_colmn, samp = dis_colmn, test4up = T, output = T, outname = "exprss_data.kegg.up", make.plot = F)
exprss_data.kegg.esg.dn <- esset.grp(exprss_data.kegg.p$less, express_data, gsets = kegg.gs, ref = norm_colmn, samp = dis_colmn, test4up = F, output = T, outname = "exprss_data.kegg.dn", make.plot = F)

gs=unique(unlist(kegg.gs[rownames(exprss_data.kegg.p$greater)[1:5]]))
essData=essGene(gs, express_data, ref = norm_colmn, samp = dis_colmn)

dir.create('~/Data/greater/special')
setwd('~/Data/greater/special')

for (gs in rownames(exprss_data.kegg.p$greater)[1:5]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  geneData(genes = kegg.gs[[gs]], exprs = essData, ref = norm_colmn,
           samp = dis_colmn, outname = outname, txt = T, heatmap = T, Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}

dir.create('~/Data/greater/all')
setwd('~/Data/greater/all')

for (gs in rownames(exprss_data.kegg.p$greater)[1:5]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  outname = paste(outname, "all", sep=".")
  geneData(genes = kegg.gs[[gs]], exprs = express_data, ref = norm_colmn,
           samp = dis_colmn, outname = outname, txt = T, heatmap = T,
           Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}

gs=unique(unlist(kegg.gs[rownames(exprss_data.kegg.p$less)[1:5]]))

dir.create('~/Data/less/special')
setwd('~/Data/less/special')

for (gs in rownames(exprss_data.kegg.p$less)[1:5]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  geneData(genes = kegg.gs[[gs]], exprs = essData, ref = norm_colmn,
           samp = dis_colmn, outname = outname, txt = T, heatmap = T, Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}

dir.create('~/Data/less/all')
setwd('~/Data/less/all')

for (gs in rownames(exprss_data.kegg.p$less)[1:5]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  outname = paste(outname, "all", sep=".")
  geneData(genes = kegg.gs[[gs]], exprs = express_data, ref = norm_colmn,
           samp = dis_colmn, outname = outname, txt = T, heatmap = T,
           Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}

gs=unique(unlist(kegg.gs[rownames(exprss_data.2d.kegg.p$greater)[1:5]]))

dir.create('~/Data/2d.greater/special')
setwd('~/Data/2d.greater/special')

for (gs in rownames(exprss_data.2d.kegg.p$greater)[1:5]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  geneData(genes = kegg.gs[[gs]], exprs = essData, ref = norm_colmn,
           samp = dis_colmn, outname = outname, txt = T, heatmap = T, Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}

dir.create('~/Data/2d.greater/all')
setwd('~/Data/2d.greater/all')

for (gs in rownames(exprss_data.2d.kegg.p$greater)[1:5]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  outname = paste(outname, "all", sep=".")
  geneData(genes = kegg.gs[[gs]], exprs = express_data, ref = norm_colmn,
           samp = dis_colmn, outname = outname, txt = T, heatmap = T,
           Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}


express_data.d  <- gagePrep(express_data, ref = norm_colmn, samp = dis_colmn, compare="as.group") # uniform pathview sections colouring

setwd('~/Data/greater')

path.ids=five_sig_greater_hsa 
path.ids2 <- substr(path.ids, 1, 8)
pv.out.list_png <- sapply(path.ids2, function(pid) pathview(gene.data = express_data.d, pathway.id = pid, species = "hsa"))
pv.out.list_pdf <- sapply(path.ids2, function(pid) pathview(gene.data = express_data.d, pathway.id = pid, species = "hsa", kegg.native=F,sign.pos="bottomleft"))

setwd('~/Data/less')

path.ids=five_sig_less_hsa
path.ids2 <- substr(path.ids, 1, 8)
pv.out.list_png <- sapply(path.ids2, function(pid) pathview(gene.data = express_data.d, pathway.id = pid, species = "hsa"))
pv.out.list_pdf <- sapply(path.ids2, function(pid) pathview(gene.data = express_data.d, pathway.id = pid, species = "hsa", kegg.native=F,sign.pos="bottomleft"))

setwd('~/Data/2d.greater')

path.ids=five_sig_2d_greater_hsa 
path.ids2 <- substr(path.ids, 1, 8)
pv.out.list_png <- sapply(path.ids2, function(pid) pathview(gene.data = express_data.d, pathway.id = pid, species = "hsa"))
pv.out.list_pdf <- sapply(path.ids2, function(pid) pathview(gene.data = express_data.d, pathway.id = pid, species = "hsa", kegg.native=F,sign.pos="bottomleft"))

