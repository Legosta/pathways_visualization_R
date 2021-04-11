library(pathview)
library(gage)
library(data.table)

# data processing 
dat <- data.frame(data.table::fread("~/rnaseqProjects/Fall2018BI/data/init_data.csv", header = T, colClasses = c(V1="character")))
rownames(dat) <- dat$V1
dat <- dat[,-1]

# set column numbers for control (normal) and experimental (virus infected) samples
nrm_clmns <- 1:6
virus_clmns <- 7:12
# upload pathways data
data(kegg.gs)

exprss_data.2d.kegg.p <- gage(dat, gsets = kegg.gs, ref = nrm_clmns, samp = virus_clmns, same.dir = F)
#choose significant gene sets
exprss_data.kegg.2d.sig <- sigGeneSet(exprss_data.2d.kegg.p, outname="exprss_data_2d.kegg")

three.sig.2d_greater_hsa <- rownames(exprss_data.kegg.2d.sig$greater)[1:3]

dir.create('~/rnaseqProjects/Fall2018BI/data/kegg2D/greater', recursive = T)
setwd('~/rnaseqProjects/Fall2018BI/data/kegg2D/greater')

write.csv(exprss_data.kegg.2d.sig$greater, file = "sigExprssData.kegg2D.greater.csv")
sigExprssData.kegg2D.greater <- as.data.frame(exprss_data.kegg.2d.sig$greater)
sigExprssData.kegg2D.greaterQValFltred <- subset(sigExprssData.kegg2D.greater, q.val < 0.05)
write.csv(sigExprssData.kegg2D.greaterQValFltred, file = "qFltredExprssData.kegg.greater.csv")

setwd('~/rnaseqProjects/Fall2018BI/data')
exprss_data.kegg2D.esg.up <- esset.grp(exprss_data.2d.kegg.p$greater, dat, gsets = kegg.gs, ref = nrm_clmns, samp = virus_clmns, test4up = T, output = F, outname = "exprss_data.kegg.up", make.plot = F)

gs=unique(unlist(kegg.gs[rownames(sigExprssData.kegg2D.greaterQValFltred)[1:3]]))
essData=essGene(gs, dat, ref = nrm_clmns, samp = virus_clmns)

dir.create('~/rnaseqProjects/Fall2018BI/data/kegg2D/greater/special')
setwd('~/rnaseqProjects/Fall2018BI/data/kegg2D/greater/special')

for (gs in rownames(sigExprssData.kegg2D.greaterQValFltred)[1:3]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  geneData(genes = kegg.gs[[gs]], exprs = essData, ref = nrm_clmns,
           samp = virus_clmns, outname = outname, txt = T, heatmap = T, Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}

dir.create('~/rnaseqProjects/Fall2018BI/data/kegg2D/greater/all')
setwd('~/rnaseqProjects/Fall2018BI/data/kegg2D/greater/all')

for (gs in rownames(sigExprssData.kegg2D.greaterQValFltred)[1:3]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  outname = paste(outname, "all", sep=".")
  geneData(genes = kegg.gs[[gs]], exprs = dat, ref = nrm_clmns,
           samp = virus_clmns, outname = outname, txt = T, heatmap = T,
           Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}

express_data.d  <- gagePrep(dat, ref = nrm_clmns, samp = virus_clmns, compare="as.group") # uniform pathview sections coloring

setwd('~/rnaseqProjects/Fall2018BI/data/kegg2D/greater/')

path.ids=three.sig.2d_greater_hsa 
path.ids2 <- substr(path.ids, 1, 8)
pv.out.list_png <- sapply(path.ids2, function(pid) pathview(gene.data = express_data.d, pathway.id = pid, species = "hsa"))
pv.out.list_pdf <- sapply(path.ids2, function(pid) pathview(gene.data = express_data.d, pathway.id = pid, species = "hsa", kegg.native=F,sign.pos="bottomleft"))
