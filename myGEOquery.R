#--------for Linix terminal----------
args <- commandArgs(trailingOnly = TRUE)
gse_name = as.character(args[1])
gpl_name <- as.character(args[2])
#------------------------------------

library(GEOquery)
library(dplyr)
gse_name <-  "GSE16357"
platform_name <- 'GPL570'

#======GeoQuery========

gse <- getGEO(gse_name, GSEMatrix=FALSE)
gse_as_matrx <- getGEO(gse_name, GSEMatrix=TRUE)

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

#========Обработка данных==========
rownames(data.matrix) <- gene_id
data.matrix <- data.matrix[- grep("///", rownames(data.matrix)),]
colnames(data.matrix) <- gse_as_matrx$`GSE16357-GPL570_series_matrix.txt.gz`$title

write.csv(data.matrix, "RAW_DATA.csv")
express_data_ <- read.csv("RAW_DATA.csv")

express_data_ <- rename(express_data_, Gene.ID.title = X)
express_data_  <-  filter(express_data_, Gene.ID.title != '') 
express_data_ <- express_data_ %>%
  group_by(Gene.ID.title) %>%
  summarise_all(mean)
rownames(express_data_) <- express_data_$Gene.ID.title
express_data <- express_data_[, -c(1)]
rownames(express_data) <- rownames(express_data_)
write.csv(express_data, "SEMI_RAW_DATA.csv")
