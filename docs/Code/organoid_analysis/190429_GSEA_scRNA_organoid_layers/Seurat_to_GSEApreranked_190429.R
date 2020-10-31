#rm(list=ls())

library(Seurat)
library(dplyr)
library(Matrix)

setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190429_GSEA_scRNA_organoid_layers/")

data(mergedData)
#load('data/mergedData.rda')


#with aleksandra may 3 - getting rid of genes with no counts

# 3 calculate log fold change diferences (compared to whole sphere) and means for each gene in each layer -------- 
#log(mean(expm1(x))+1) is used to un-log the log normalized data, take the mean, and then take the natural logorithm. then we get rid of zeros. then we take the diference in mean between each layer and all layers pooled.
calc.logMean = function(x) log(mean(expm1(x))+1)
means = as.data.frame(t(apply(mergedData@assays$RNA@data,1,tapply, mergedData@active.ident, calc.logMean)))
sel <- rowSums(means) > 0
means <- means[sel,]
difs = as.data.frame(means - means[,'SPHERE'])

# 4 save file as 'preranked' rda file
# save to rda file
save(difs, means, file = "R_scripts/190429_GSEA_scRNA_organoid_layers/data/PrerankOrg.rda")

# now, go back and remove the same genes from the scaled/normalized data so that the lists match
data(scaledData)
genes.use <- rownames(difs)
subset.matrix <- organoid2@assays$RNA@data[genes.use, ]
organoid_GSEA <- CreateSeuratObject(subset.matrix)
old_matrix_idents <- organoid2@meta.data[, "orig.ident"]
organoid_GSEA <- AddMetaData(object = organoid_GSEA, metadata = old_matrix_idents, col.name = "orig.ident")
Idents(object = organoid_GSEA) <- "orig.ident"
#and save it for later import into GSEA
save(organoid_GSEA, file = "data/scaledDataGSEA.rda")

library(Seurat)
library(Matrix)
library(dplyr)
library(stats)

#functions (from source)
source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190429_GSEA_scRNA_organoid_layers/Codes/01GSEA_pGSEAfunction.R")

gmtlist(path.gmtFolder="/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190429_GSEA_scRNA_organoid_layers/MSigDB/MSigDB_subsets",
        cut.version='.v6.1.symbols.gmt',
        rda.file="MSigDB/MSigDB_subsets.rda")


# load pre-ranked ave_logFC for layer vs. total sphere
gsea_CENTRAL <- pGSEA(rank.metric=difs$c, names.metric=rownames(difs),msigdb="MSigDB/MSigDB_subsets.rda")
gsea_MIDDLE <- pGSEA(rank.metric=difs$m, names.metric=rownames(difs), msigdb="MSigDB/MSigDB_subsets.rda")
gsea_SURFACE <- pGSEA(rank.metric=difs$s, names.metric=rownames(difs),msigdb="MSigDB/MSigDB_subsets.rda")


save(gsea_CENTRAL, gsea_MIDDLE, gsea_SURFACE,
     file = 'data/pGSEA.rda')

write.gsea(location='CENTRAL', gsea.list=gsea_CENTRAL)
write.gsea(location='MIDDLE', gsea.list=gsea_MIDDLE)
write.gsea(location='SURFACE', gsea.list=gsea_SURFACE)



