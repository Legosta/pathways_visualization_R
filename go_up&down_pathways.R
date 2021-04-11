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
data(go.gs)

exprss_data.go <- gage(as.matrix(dat), gsets = go.gs, ref = nrm_clmns, samp = virus_clmns)
exprss_data.go.sig <- sigGeneSet(exprss_data.go, outname="exprss_data.go")
threeSigGreater <- rownames(exprss_data.go.sig$greater)[1:3]
threeSigLess <- rownames(exprss_data.go.sig$less)[1:3]

##__activated_genes_in_virus_infected_samples__
dir.create('~/rnaseqProjects/Fall2018BI/data/go/greater', recursive = T)
setwd('~/rnaseqProjects/Fall2018BI/data/go/greater')
write.csv(exprss_data.go.sig$greater, file = "sigExprssData.go.greater.csv")

sigExprssData.go.greater <- as.data.frame(exprss_data.go.sig$greater)
sigExprssData.go.greaterQValFltred <- subset(sigExprssData.go.greater, q.val < 0.05)
write.csv(sigExprssData.go.greaterQValFltred, file = "qFltredExprssData.go.greater.csv")

setwd('~/rnaseqProjects/Fall2018BI/data')
exprss_data.go.esg.up <- esset.grp(exprss_data.go$greater, dat, gsets = go.gs, ref = nrm_clmns, samp = virus_clmns, test4up = T, output = F, outname = "exprss_data.go.up", make.plot = F)

gs=unique(unlist(go.gs[rownames(sigExprssData.go.greaterQValFltred)[1:3]]))
essData=essGene(gs, dat, ref = nrm_clmns, samp = virus_clmns)

dir.create('~/rnaseqProjects/Fall2018BI/data/go/greater/special')
setwd('~/rnaseqProjects/Fall2018BI/data/go/greater/special')

for (gs in rownames(sigExprssData.go.greaterQValFltred)[1:3]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  geneData(genes = go.gs[[gs]], exprs = essData, ref = nrm_clmns,
           samp = virus_clmns, outname = outname, txt = T, heatmap = T, Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}

dir.create('~/rnaseqProjects/Fall2018BI/data/go/greater/all')
setwd('~/rnaseqProjects/Fall2018BI/data/go/greater/all')

for (gs in rownames(sigExprssData.go.greaterQValFltred)[1:3]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  outname = paste(outname, "all", sep=".")
  geneData(genes = go.gs[[gs]], exprs = dat, ref = nrm_clmns,
           samp = virus_clmns, outname = outname, txt = T, heatmap = T,
           Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}

##__inhibited_genes_in_virus_infected_samples__
dir.create('~/rnaseqProjects/Fall2018BI/data/go/less', recursive = T)
setwd('~/rnaseqProjects/Fall2018BI/data/go/less')
write.csv(exprss_data.go.sig$less, file = "sigExprssData.go.less.csv")

sigExprssData.go.less <- as.data.frame(exprss_data.go.sig$greater)
sigExprssData.go.lessQValFltred <- subset(sigExprssData.go.less, q.val < 0.05)
write.csv(sigExprssData.go.lessQValFltred, file = "qFltredExprssData.go.less.csv")

setwd('~/rnaseqProjects/Fall2018BI/data')
exprss_data.go.esg.down <- esset.grp(exprss_data.go$less, dat, gsets = go.gs, ref = nrm_clmns, samp = virus_clmns, test4up = F, output = F, outname = "exprss_data.go.less", make.plot = F)

gs=unique(unlist(go.gs[rownames(sigExprssData.go.lessQValFltred)[1:3]]))
essData=essGene(gs, dat, ref = nrm_clmns, samp = virus_clmns)

dir.create('~/rnaseqProjects/Fall2018BI/data/go/less/special', recursive = T)
setwd('~/rnaseqProjects/Fall2018BI/data/go/less/special')

for (gs in rownames(sigExprssData.go.lessQValFltred)[1:3]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  geneData(genes = go.gs[[gs]], exprs = essData, ref = nrm_clmns,
           samp = virus_clmns, outname = outname, txt = T, heatmap = T, Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}

dir.create('~/rnaseqProjects/Fall2018BI/data/go/less/all')
setwd('~/rnaseqProjects/Fall2018BI/data/go/less/all')

for (gs in rownames(sigExprssData.go.lessQValFltred)[1:3]) {
  outname = gsub(" |:|/", "_", substr(gs, 10, 100))
  outname = paste(outname, "all", sep=".")
  geneData(genes = go.gs[[gs]], exprs = dat, ref = nrm_clmns,
           samp = virus_clmns, outname = outname, txt = T, heatmap = T,
           Colv = F, Rowv = F, dendrogram = "none", limit = 3, scatterplot = T)
}


# express_data.d  <- gagePrep(dat, ref = nrm_clmns, samp = virus_clmns, compare="as.group") # uniform pathview sections coloring
# 
# setwd('~/rnaseqProjects/Fall2018BI/data/go/greater/')
# 
# path.ids=threeSigGreater 
# path.ids2 <- substr(path.ids, 1, 8)
# pv.out.list_png <- sapply(path.ids2, function(pid) pathview(gene.data = express_data.d, pathway.id = pid, species = "hsa"))
# pv.out.list_pdf <- sapply(path.ids2, function(pid) pathview(gene.data = express_data.d, pathway.id = pid, species = "hsa", kegg.native=F,sign.pos="bottomleft"))
