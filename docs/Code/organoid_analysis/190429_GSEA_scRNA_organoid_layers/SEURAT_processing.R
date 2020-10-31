library(Seurat)
library(dplyr)
library(Matrix)

#1 load and scale data in Seurat -----------------
getwd()
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/")

# Load the spheroid dataset
raw_center <- Read10X("counts_matrices/oC_center")
raw_middle <- Read10X("counts_matrices/oC_middle/")
raw_surface <- Read10X("counts_matrices/oC_surface/")
# Initialize the Seurat object with the raw (non-normalized data).
center <- CreateSeuratObject(counts = raw_center, min.cells = 3, min.features = 300, project = "center")
middle <- CreateSeuratObject(counts = raw_middle, min.cells = 3, min.features = 300, project = "middle")
surface <- CreateSeuratObject(counts = raw_surface, min.cells = 3, min.features= 300, project = "surface")
#select mito genes, and filter layers

#seusrat3
center[["percent.mt"]] <- PercentageFeatureSet(center, pattern = "^MT-")
middle[["percent.mt"]] <- PercentageFeatureSet(middle, pattern = "^MT-")
surface[["percent.mt"]] <- PercentageFeatureSet(surface, pattern = "^MT-")


#filter Seurat3
center <- subset(center, subset = nFeature_RNA > 580 & nFeature_RNA < 5000 & percent.mt < 30 & percent.mt > 4.5 & nCount_RNA > 489 & nCount_RNA < 15000)
middle <- subset(middle, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 30 & percent.mt > 5.5 & nCount_RNA > 478 & nCount_RNA < 15000)
surface <- subset(surface, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 30& percent.mt > 5 & nCount_RNA > 610 & nCount_RNA < 15000)

# merge layers
# merge datasets organoid2 = merged seurat object (temp_merge is a placeholder to combine all 3 objects)
temp_merg <- merge(x = center, y = middle, project = "ORGANOID")
organoid2 <- merge(x = temp_merg, y = surface, project = "ORGANOID")


#SEURAT 3
organoid2[["percent.mt"]] <- PercentageFeatureSet(center, pattern = "^MT-")



# normalize
organoid2 <- NormalizeData(object = organoid2, normalization.method = "LogNormalize", scale.factor = 1e4)

# detection of variable genes
organoid2 <- FindVariableFeatures(object = organoid2, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
# SCALE: regress out cell-cell variation by number of detected molecules percentage mitochondrial gene content.
organoid2 <- ScaleData(object = organoid2, vars.to.regress = c("nCounts_RNA", "percent.mt"))

#SEURAT3
all.genes <- rownames(organoid2)
organoid2 <- ScaleData(organoid2, features = all.genes)


# 2 add 'SPHERE' ident - 'SPHERE' simply denotes 'all cells ---------
# load("mergedData_scaled.rda") # actually called 'organoid2' in this case
# a) add SPHERE ident
x = organoid2
x = SetIdent(x, cells.use = NULL, value = rep("SPHERE",nrow(x@meta.data)))
inf = x@meta.data
nam=sapply(strsplit(rownames(inf), "_"), function(x) paste(x[1],x[2], sep="_"))
nam=paste("SPHERE", nam,1:nrow(inf), sep="_")
mergedData <- merge(organoid2, x, project="SPHERE", min.cells=0, min.genes=0, do.normalize=FALSE, do.center=FALSE, do.scale=FALSE, add.cell.id2=nam)

save(organoid2, file = "R_scripts/190429_GSEA_scRNA_organoid_layers/data/scaledData.rda")
save(mergedData, file = "R_scripts/190429_GSEA_scRNA_organoid_layers/data/mergedData.rda")

