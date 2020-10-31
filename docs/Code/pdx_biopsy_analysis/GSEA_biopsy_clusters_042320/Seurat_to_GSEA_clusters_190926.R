#GSEA biopsy CLUSTERS (to estimate positional differences in cell states) 9/26/19
#start with scaled data - change to cluster ident - merge here
library(Seurat)
library(dplyr)
library(Matrix)
library(stats)
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/")
biopsy <- readRDS(file = "filtered_biopsy_SO_final_042220.rds")
Idents(biopsy) <- "seurat_clusters"

# add 'ALL' ident - 'ALL' simply denotes 'all cells ---------
x = biopsy
x = SetIdent(x, cells = NULL, value = rep("ALL",nrow(x@meta.data)))
inf = x@meta.data
nam=sapply(strsplit(rownames(inf), "_"), function(x) paste(x[1],x[2], sep="_"))
nam=paste("ALL", nam,1:nrow(inf), sep="_")
x <- RenameCells(x, add.cell.id = nam)
mergedData <- merge(biopsy, x, project="ALL")
mergedData$new.ident <- Idents(mergedData)

#save(mergedData, file = "GSEA_biopsy_clusters_042320/data/mergedData.rda")
#data(mergedData)

#set identity to clusters
mergedData <- SetIdent(mergedData, value = 'new.ident')

# 3 calculate log fold change diferences (compared to whole sphere) and means for each gene in each cluster -------- 
#log(mean(expm1(x))+1) is used to un-log the log normalized data, take the mean, and then take the natural logorithm. then we get rid of zeros. then we take the diference in mean between each layer and all layers pooled.
calc.logMean = function(x) log(mean(expm1(x))+1)
means = as.data.frame(t(apply(mergedData@assays$RNA@data,1,tapply, mergedData@active.ident, calc.logMean)))
sel <- rowSums(means) > 0
means <- means[sel,]
difs = as.data.frame(means - means[,'ALL'])

# 4 save file as 'preranked' rda file
save(difs, means, file = "GSEA_biopsy_clusters_042320/data/PrerankBiopsy.rda")

# now, go back and remove the same genes from the scaled/normalized data so that the lists match
#data("Data_scaled_clustered")
genes.use <- rownames(difs)
subset.matrix <- biopsy@assays$RNA@data[genes.use, ]
biopsy_GSEA <- CreateSeuratObject(subset.matrix)
old_matrix_idents <- biopsy@meta.data[, "seurat_clusters"]
biopsy_GSEA <- AddMetaData(object = biopsy_GSEA, metadata = old_matrix_idents, col.name = "seurat_clusters")
Idents(object = biopsy_GSEA) <- "seurat_clusters"
#and save it for later import into GSEA
save(biopsy_GSEA, file = "GSEA_biopsy_clusters_042320/data/scaledClusterdDataGSEA.rda")


#functions (from source)
source("GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")

gmtlist(path.gmtFolder="GSEA_biopsy_clusters_042320/MSigDB/MSigDB_subsets",
        cut.version='.v6.1.symbols.gmt',
        rda.file="GSEA_biopsy_clusters_042320/MSigDB/MSigDB_subsets.rda")

setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320")
# load pre-ranked ave_logFC for cluster vs. total spheroid
gsea_0 <- pGSEA(rank.metric=difs$'0', names.metric=rownames(difs),msigdb="MSigDB/MSigDB_subsets.rda")
gsea_1 <- pGSEA(rank.metric=difs$'1', names.metric=rownames(difs),msigdb="MSigDB/MSigDB_subsets.rda")
gsea_2 <- pGSEA(rank.metric=difs$'2', names.metric=rownames(difs),msigdb="MSigDB/MSigDB_subsets.rda")
gsea_3 <- pGSEA(rank.metric=difs$'3', names.metric=rownames(difs),msigdb="MSigDB/MSigDB_subsets.rda")
gsea_4 <- pGSEA(rank.metric=difs$'4', names.metric=rownames(difs),msigdb="MSigDB/MSigDB_subsets.rda")
gsea_5 <- pGSEA(rank.metric=difs$'5', names.metric=rownames(difs),msigdb="MSigDB/MSigDB_subsets.rda")


save(gsea_0, gsea_1, gsea_2, gsea_3, gsea_4, gsea_5,
     file = 'data/pGSEA.rda')

write.gsea(location='0', gsea.list=gsea_0)
write.gsea(location='1', gsea.list=gsea_1)
write.gsea(location='2', gsea.list=gsea_2)
write.gsea(location='3', gsea.list=gsea_3)
write.gsea(location='4', gsea.list=gsea_4)
write.gsea(location='5', gsea.list=gsea_5)


