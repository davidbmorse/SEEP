#GSEA biopsy CLUSTERS (to estimate positional differences in cell states) 9/26/19
#start with scaled data - change to cluster ident - merge here
library(Seurat)
library(dplyr)
library(Matrix)
library(stats)

setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/")
biopsy <- readRDS(file = "filtered_biopsy_SO_final_042220.rds")
Idents(biopsy) <- "orig.ident"

# add 'ALL' ident - 'ALL' simply denotes 'all cells ---------
x = biopsy
x = SetIdent(x, cells = NULL, value = rep("ALL",nrow(x@meta.data)))
inf = x@meta.data
nam=sapply(strsplit(rownames(inf), "_"), function(x) paste(x[1],x[2], sep="_"))
nam=paste("ALL", nam,1:nrow(inf), sep="_")
x <- RenameCells(x, add.cell.id = nam)
mergedData <- merge(biopsy, x, project="ALL")
mergedData$new.ident <- Idents(mergedData)

#save(mergedData, file = "GSEA_biopsy_layers_042320/data/mergedData.rda")
#data(mergedData)

#set identity to clusters
mergedData <- SetIdent(mergedData, value = 'new.ident')

# 3 calculate log fold change diferences (compared to whole sphere) and means for each gene in each cluster
#log(mean(expm1(x))+1) is used to un-log the log normalized data, take the mean, and then take the natural logorithm. then we get rid of zeros. then we take the diference in mean between each layer and all layers pooled.
calc.logMean = function(x) log(mean(expm1(x))+1)
means = as.data.frame(t(apply(mergedData@assays$RNA@data,1,tapply, mergedData@active.ident, calc.logMean)))
sel <- rowSums(means) > 0
means <- means[sel,]
difs = as.data.frame(means - means[,'ALL'])

# 4 save file as 'preranked' rda file
save(difs, means, file = "GSEA_biopsy_layers_042320/data/PrerankBiopsy.rda")

# now, go back and remove the same genes from the scaled/normalized data so that the lists match
#data("Data_scaled_clustered")
genes.use <- rownames(difs)
subset.matrix <- biopsy@assays$RNA@data[genes.use, ]
biopsy_GSEA <- CreateSeuratObject(subset.matrix)
old_matrix_idents <- biopsy@meta.data[, "orig.ident"]
biopsy_GSEA <- AddMetaData(object = biopsy_GSEA, metadata = old_matrix_idents, col.name = "orig.ident")
Idents(object = biopsy_GSEA) <- "orig.ident"
#and save it for later import into GSEA
save(biopsy_GSEA, file = "GSEA_biopsy_layers_042320/data/scaledClusterdDataGSEA.rda")


#functions (from source)
source("GSEA_biopsy_layers_042320/Codes/01GSEA_pGSEAfunction.R")

gmtlist(path.gmtFolder="GSEA_biopsy_layers_042320/MSigDB/MSigDB_subsets",
        cut.version='.v6.1.symbols.gmt',
        rda.file="GSEA_biopsy_layers_042320/MSigDB/MSigDB_subsets.rda")

setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_layers_042320")
# load pre-ranked ave_logFC for cluster vs. total spheroid
gsea_center <- pGSEA(rank.metric=difs$'center', names.metric=rownames(difs),msigdb="MSigDB/MSigDB_subsets.rda")
gsea_middle <- pGSEA(rank.metric=difs$'middle', names.metric=rownames(difs),msigdb="MSigDB/MSigDB_subsets.rda")
gsea_surface <- pGSEA(rank.metric=difs$'surface', names.metric=rownames(difs),msigdb="MSigDB/MSigDB_subsets.rda")

save(gsea_center, gsea_middle, gsea_surface,
     file = 'data/pGSEA.rda')

write.gsea(location='center', gsea.list=gsea_center)
write.gsea(location='middle', gsea.list=gsea_middle)
write.gsea(location='surface', gsea.list=gsea_surface)



