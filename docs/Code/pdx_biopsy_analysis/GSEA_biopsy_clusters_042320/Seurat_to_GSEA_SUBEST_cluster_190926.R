#GSEA biopsy CLUSTERS (to estimate positional differences in cell states) 9/26/19
#start with scaled data - change to cluster ident - merge here

# use seurat v3
#install.packages('Seurat')
library(Seurat)
library(dplyr)
library(Matrix)
library(stats)

setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/190926_GSEA_PDX_biopsy23_clusters")
getwd()

#load data used to make 'figure2_biopsy' saved here as 'Data_scaled_clustered.rda'
PDX_biopsy2_3 <- readRDS(file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_exports/r_data/4th_merge_mouse/PDX_biopsy2_3_MERGE_final.rds")

#set identity to clusters
PDX_biopsy2_3 <- SetIdent(PDX_biopsy2_3, value = "orig.ident")

# add 'BACKGROUND' ident - 'BACKGROUND' denotes a downsampling of 150 cells from each layer ---------
x <- subset(PDX_biopsy2_3, idents = c("c", "m", "s"), downsample = 150)

#change idents back to clusters
PDX_biopsy2_3 <- SetIdent(PDX_biopsy2_3, value = "clusterID")
x <- SetIdent(x, value = "clusterID")

x = SetIdent(x, cells = NULL, value = rep("BACKGROUND",nrow(x@meta.data)))
inf = x@meta.data
nam=sapply(strsplit(rownames(inf), "_"), function(x) paste(x[1],x[2], sep="_"))
nam=paste("BACKGROUND", nam,1:nrow(inf), sep="_")
x <- RenameCells(x, add.cell.id = nam)
mergedData <- merge(PDX_biopsy2_3, x, project="BACKGROUND")
mergedData$new.ident <- Idents(mergedData)

save(mergedData, file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/190926_GSEA_PDX_biopsy23_clusters/data_evenBackground/mergedData.rda")

#data_evenBackground(mergedData)
#set identity to layers
mergedData <- SetIdent(mergedData, value = 'new.ident')

# 3 calculate log fold change diferences (compared to whole sphere) and means for each gene in each cluster -------- 
#log(mean(expm1(x))+1) is used to un-log the log normalized data, take the mean, and then take the natural logorithm. then we get rid of zeros. then we take the diference in mean between each layer and all layers pooled.
calc.logMean = function(x) log(mean(expm1(x))+1)
means = as.data.frame(t(apply(mergedData@assays$RNA@data,1,tapply, mergedData@active.ident, calc.logMean)))
sel <- rowSums(means) > 0
means <- means[sel,]
difs = as.data.frame(means - means[,'BACKGROUND'])

# 4 save file as 'preranked' rda file
# save to rda file
save(difs, means, file = "data_evenBackground/PrerankBiopsy.rda")

# now, go back and remove the same genes from the scaled/normalized data so that the lists match
#data_evenBackground("Data_scaled_clustered")
genes.use <- rownames(difs)
subset.matrix <- PDX_biopsy2_3@assays$RNA@data[genes.use, ]
PDX_biopsy2_3_GSEA <- CreateSeuratObject(subset.matrix)
old_matrix_idents <- PDX_biopsy2_3@meta.data[, "clusterID"]
PDX_biopsy2_3_GSEA <- AddMetaData(object = PDX_biopsy2_3_GSEA, metadata = old_matrix_idents, col.name = "clusterID")
Idents(object = PDX_biopsy2_3_GSEA) <- "clusterID"
#and save it for later import into GSEA
save(PDX_biopsy2_3_GSEA, file = "data_evenBackground/scaledClusterdDataGSEA.rda")


#functions (from source)
source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/190926_GSEA_PDX_biopsy23_clusters/Codes/01GSEA_pGSEAfunction.R")

gmtlist(path.gmtFolder="/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/190926_GSEA_PDX_biopsy23_clusters/MSigDB/MSigDB_subsets",
        cut.version='.v6.1.symbols.gmt',
        rda.file="MSigDB/MSigDB_subsets.rda")


# load pre-ranked ave_logFC for cluster vs. total spheroid
gsea_0 <- pGSEA(rank.metric=difs$'0', names.metric=rownames(difs),msigdb="MSigDB/MSigDB_subsets.rda")
gsea_1 <- pGSEA(rank.metric=difs$'1', names.metric=rownames(difs),msigdb="MSigDB/MSigDB_subsets.rda")
gsea_2 <- pGSEA(rank.metric=difs$'2', names.metric=rownames(difs),msigdb="MSigDB/MSigDB_subsets.rda")
gsea_3 <- pGSEA(rank.metric=difs$'3', names.metric=rownames(difs),msigdb="MSigDB/MSigDB_subsets.rda")
gsea_4 <- pGSEA(rank.metric=difs$'4', names.metric=rownames(difs),msigdb="MSigDB/MSigDB_subsets.rda")
gsea_5 <- pGSEA(rank.metric=difs$'5', names.metric=rownames(difs),msigdb="MSigDB/MSigDB_subsets.rda")


save(gsea_0, gsea_1, gsea_2, gsea_3, gsea_4, gsea_5,
     file = 'data_evenBackground/pGSEA.rda')

write.gsea(location='0_bkgd', gsea.list=gsea_0)
write.gsea(location='1_bkgd', gsea.list=gsea_1)
write.gsea(location='2_bkgd', gsea.list=gsea_2)
write.gsea(location='3_bkgd', gsea.list=gsea_3)
write.gsea(location='4_bkgd', gsea.list=gsea_4)
write.gsea(location='5_bkgd', gsea.list=gsea_5)


