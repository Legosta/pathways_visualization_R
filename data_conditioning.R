# #--------for Linix terminal----------
# args <- commandArgs(trailingOnly = TRUE)
# gse_name = as.character(args[1])
# gpl_name <- as.character(args[2])
# norm_colmn <- as.integer(args[3])
# dis_colmn <- as.integer(args[4])
# #------------------------------------

# packages

library(GEOquery)
library(dplyr)
#library(pathview)
library(gage)
library("org.Hs.eg.db")
#library(clusterProfiler)

setwd('~/rnaseqProjects/Fall2018BI/data')

# init data (later will be change)

gse_name <-  "GSE16357"
platform_name <- 'GPL570'

# work with GEO files

gse <- getGEO(gse_name, GSEMatrix=FALSE)

gsmplatforms <- lapply(GSMList(gse),function(x) {Meta(x)$platform_id})

gsmlist = Filter(function(gsm) {Meta(gsm)$platform_id==platform_name},GSMList(gse))
gene_id <- Table(GPLList(gse)[[1]])$ENTREZ_GENE_ID
probesets <- Table(GPLList(gse)[[1]])$ID

data.matrix <- do.call(cbind,lapply(gsmlist, function(x) 
{tab <- Table(x)
mymatch <- match(probesets,tab$ID_REF)
return(tab$VALUE[mymatch])
}))

data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})

# data conditioning

rownames(data.matrix) <- gene_id
data.matrix <- data.matrix[- grep("///", rownames(data.matrix)),] # remove duplicates

# remove rows with empty gene names (just X)
dat <- as.data.frame(data.matrix) %>% 
  filter(row.names(.) != "X") %>% 
  group_by(rownames(.)) %>%
  summarise_all(mean)

colnames(dat)[1] = "gene_id"
dat$gene_id <- gsub("(^X)", "", dat$gene_id)
dat <- as.data.frame(dat)
rownames(dat) <- dat$gene_id
dat <- dat[, -1]
rm(dat_)
data.table::fwrite(dat, file="init_data.csv", row.names = T)
