library(pathview)
library(gage)
library(data.table)

dat <- data.frame(data.table::fread("~/rnaseqProjects/Fall2018BI/data/init_data.csv", header = T, colClasses = c(V1="character")))
rownames(dat) <- dat$V1
dat <- dat[,-1]

# set column numbers for control (normal) and experimental (virus infected) samples
nrm_clmns <- 1:6
virus_clmns <- 7:12
# upload pathways data
data(kegg.gs)

exprss_data.kegg <- gage(dat, gsets = kegg.gs, ref = nrm_clmns, samp = virus_clmns)
exprss_data.kegg.sig <- sigGeneSet(exprss_data.kegg, outname="exprss_data.kegg")
threeSigGreater <- rownames(exprss_data.kegg.sig$greater)[1:3]
threeSigLess <- rownames(exprss_data.kegg.sig$less)[1:3]

##__activated_genes_in_virus_infected_samples__
dir.create('~/rnaseqProjects/Fall2018BI/data/kegg/greater', recursive = T)
setwd('~/rnaseqProjects/Fall2018BI/data/kegg/greater')
write.csv(exprss_data.kegg.sig$greater, file = "sigExprssData.kegg.greater.csv")

sigExprssData.kegg.greater <- as.data.frame(exprss_data.kegg.sig$greater)
sigExprssData.kegg.greaterQValFltred <- subset(sigExprssData.kegg.greater, q.val < 0.05)
write.csv(sigExprssData.kegg.greaterQValFltred, file = "qFltredExprssData.kegg.greater.csv")

setwd('~/rnaseqProjects/Fall2018BI/data')
exprss_data.kegg.esg.up <- esset.grp(exprss_data.kegg$greater, dat, gsets = kegg.gs, ref = nrm_clmns, samp = virus_clmns, test4up = T, output = F, outname = "exprss_data.kegg.up", make.plot = F)

gs=unique(unlist(kegg.gs[rownames(sigExprssData.kegg.greaterQValFltred)[1:3]]))
essData=essGene(gs, dat, ref = nrm_clmns, samp = virus_clmns)

dir.create('~/rnaseqProjects/Fall2018BI/data/kegg/greater/special')
setwd('~/rnaseqProjects/Fall2018BI/data/kegg/greater/special')

for (gs in rownames(sigExprssData.kegg.greaterQValFltred)[1:3]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  geneData(genes = kegg.gs[[gs]], exprs = essData, ref = nrm_clmns,
           samp = virus_clmns, outname = outname, txt = T, heatmap = T, Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}

dir.create('~/rnaseqProjects/Fall2018BI/data/kegg/greater/all')
setwd('~/rnaseqProjects/Fall2018BI/data/kegg/greater/all')

for (gs in rownames(sigExprssData.kegg.greaterQValFltred)[1:3]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  outname = paste(outname, "all", sep=".")
  geneData(genes = kegg.gs[[gs]], exprs = dat, ref = nrm_clmns,
           samp = virus_clmns, outname = outname, txt = T, heatmap = T,
           Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}

express_data.d  <- gagePrep(dat, ref = nrm_clmns, samp = virus_clmns, compare="as.group") # uniform pathview sections coloring

setwd('~/rnaseqProjects/Fall2018BI/data/kegg/greater/')

path.ids=threeSigGreater 
path.ids2 <- substr(path.ids, 1, 8)
pv.out.list_png <- sapply(path.ids2, function(pid) pathview(gene.data = express_data.d, pathway.id = pid, species = "hsa"))
pv.out.list_pdf <- sapply(path.ids2, function(pid) pathview(gene.data = express_data.d, pathway.id = pid, species = "hsa", kegg.native=F,sign.pos="bottomleft"))

##__inhibited_genes_in_virus_infected_samples__

write.csv(exprss_data.kegg.sig$less, file = "sigExprssData.kegg.less.csv")
sigExprssData.kegg.less <- as.data.frame(exprss_data.kegg.sig$less)
sigExprssData.kegg.lessQValFltred <- subset(sigExprssData.kegg.less, sigExprssData.kegg.less$q.val < 0.05)
write.csv(sigExprssData.kegg.lessQValFltred, file = "qFltredExprssData.kegg.less.csv")

setwd('~/rnaseqProjects/Fall2018BI/data')
exprss_data.kegg.esg.down <- esset.grp(exprss_data.kegg$less, dat, gsets = kegg.gs, ref = nrm_clmns, samp = virus_clmns, test4up = F, output = F, outname = "exprss_data.kegg.down", make.plot = F)

gs=unique(unlist(kegg.gs[rownames(sigExprssData.kegg.lessQValFltred)[1:3]]))

dir.create('~/rnaseqProjects/Fall2018BI/data/kegg/less/special')
setwd('~/rnaseqProjects/Fall2018BI/data/kegg/less/special')

for (gs in rownames(sigExprssData.kegg.lessQValFltred)[1:3]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  geneData(genes = kegg.gs[[gs]], exprs = essData, ref = nrm_clmns,
           samp = virus_clmns, outname = outname, txt = T, heatmap = T, Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}

dir.create('~/rnaseqProjects/Fall2018BI/data/kegg/less/all')
setwd('~/rnaseqProjects/Fall2018BI/data/kegg/less/all')

for (gs in rownames(sigExprssData.kegg.lessQValFltred)[1:3]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  outname = paste(outname, "all", sep=".")
  geneData(genes = kegg.gs[[gs]], exprs = dat, ref = nrm_clmns,
           samp = virus_clmns, outname = outname, txt = T, heatmap = T,
           Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}
setwd('~/rnaseqProjects/Fall2018BI/data/kegg/less/')

path.ids=threeSigLess
path.ids2 <- substr(path.ids, 1, 8)
pv.out.list_png <- sapply(path.ids2, function(pid) pathview(gene.data = express_data.d, pathway.id = pid, species = "hsa"))
pv.out.list_pdf <- sapply(path.ids2, function(pid) pathview(gene.data = express_data.d, pathway.id = pid, species = "hsa", kegg.native=T,sign.pos="bottomleft"))
