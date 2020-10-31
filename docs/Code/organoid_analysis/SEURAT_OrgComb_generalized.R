#Setup----
#organoids still processed under Seurat V2: Please depreciate your seurat to V2 rather than V3 prior to running the flowwing code.
#for depreciating to seurat version 2.3.4
source("https://z.umn.edu/archived-seurat")

library(Seurat)
library(dplyr)
library(Matrix)
library(data.table) #for extra data table wrangling

getwd()
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/")

# Load the spheroid dataset
raw_center <- Read10X("counts_matrices/oC_center")
raw_middle <- Read10X("counts_matrices/oC_middle/")
raw_surface <- Read10X("counts_matrices/oC_surface/")

# Initialize the Seurat object with the raw (non-normalized data).
# Keep all genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at least 200 detected genes
center <- CreateSeuratObject(raw.data = raw_center, min.cells = 3, min.genes = 300, project = "center")
middle <- CreateSeuratObject(raw.data = raw_middle, min.cells = 3, min.genes = 300, project = "middle")
surface <- CreateSeuratObject(raw.data = raw_surface, min.cells = 3, min.genes = 300, project = "surface")

#filter center cells ----------------------------------------------

##center
# calculate the percentage of mitochondrial genes and store it in percent.mito using AddMetaData.
mito.genes <- grep(pattern = "^MT-", x = rownames(x = center@data), value = TRUE)
#percent.mito <- colSums(center@raw.data[mito.genes, ]) / colSums(center@raw.data)
percent.mito <- Matrix::colSums(center@raw.data[mito.genes, ])/Matrix::colSums(center@raw.data)
# AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
center <- AddMetaData(object = center, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = center, point.size.use = .01, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

#plot genes/mito and UMI/Genes
par(mfrow = c(1, 2))
GenePlot(object = center, gene1 = "nUMI", gene2 = "percent.mito", pch.use = ".")
GenePlot(object = center, gene1 = "nUMI", gene2 = "nGene", pch.use = ".")

#check out gene and UMI abundance
#UMI
hist(center@meta.data$nUMI, breaks = 2500, xlim = c(0,2000), main = "UMI's per cell") 
hist(center@meta.data$nGene, breaks = 440, xlim = c(0,1500), main = "genes per cell")
hist(center@meta.data$percent.mito, breaks = 150, xlim = c(0,0.5), main = "fraction of mitochondrial reads per cell")

#find real cells vs background
barcode_rank <- rank(-center@meta.data$nUMI)
plot(barcode_rank, center@meta.data$nUMI, xlim=c(1,8000))

#make logoritmic
log_lib_size <- log10(center@meta.data$nUMI)
plot(barcode_rank, log_lib_size, xlim=c(1,8000))

#find inflection point
o <- order(barcode_rank)
log_lib_size <- log_lib_size[o]
barcode_rank <- barcode_rank[o]

rawdiff <- diff(log_lib_size)/diff(barcode_rank)
inflection <- which(rawdiff == min(rawdiff[100:length(rawdiff)], na.rm=TRUE))

plot(barcode_rank, log_lib_size, main = "cells ranked by number of UMIs - inflection test")
abline(v=inflection, col="red", lwd=2)
#SEEMS too high


#CellRanger (10X) technique
n_cells <- center@data@Dim[2]
totals <- center@meta.data$nUMI
totals <- sort(totals, decreasing = TRUE)
# 99th percentile of top n_cells divided by 10
thresh = totals[round(0.01*n_cells)]/10
plot(totals, main = "cells ranked by number of UMIs - 10-fold test")
abline(h=thresh, col="red", lwd=2)


#fix a mixture model using mixtools
set.seed(-92497)
require("mixtools")
mix <- normalmixEM(log_lib_size)
plot(mix, breaks = 100, which=2, xlab2="log(mol per cell)") + abline(v = split, col="red")
#split defined below
p1 <- dnorm(log_lib_size, mean=mix$mu[1], sd=mix$sigma[1])
p2 <- dnorm(log_lib_size, mean=mix$mu[2], sd=mix$sigma[2])
if (mix$mu[1] < mix$mu[2]) {
  split <- min(log_lib_size[p2 > p1])
} else {
  split <- min(log_lib_size[p1 > p2])
}
#identify cells using this split point and calculate the TPR and Recall.

threshold <- 10^split
umi_per_barcode <- data.frame(center@cell.names, center@meta.data$nUMI)
cells <- umi_per_barcode[umi_per_barcode[,2] > threshold,1]
TPR <- sum(cells %in% center@cell.names)/length(cells)
Recall <- sum(cells %in% center@cell.names)/length(center@cell.names)

#Filter here
par(mfrow = c(1, 2))
GenePlot(object = center, gene1 = "nUMI", gene2 = "percent.mito", pch.use = ".")
GenePlot(object = center, gene1 = "nUMI", gene2 = "nGene", pch.use = ".")
# We filter out cells that have UMI counts over 489
center <- FilterCells(object = center, subset.names = c("nGene", "percent.mito", "nUMI"), low.thresholds = c(5, -Inf, 489), high.thresholds = c(6000, 0.4, 15000))

#after filtering for UMI, the bimodality of the gene expression becomes visible. At this point, filter genes removing drops with an obvious lower number of genes expressed.
hist(center@meta.data$nUMI, breaks = 2500, xlim = c(0,6000), main = "UMI's per cell") 
hist(center@meta.data$nGene, breaks = 440, xlim = c(0,1500), main = "genes per cell")
hist(center@meta.data$percent.mito, breaks = 150, xlim = c(0,0.3), main = "fraction of mitochondrial reads per cell") +  abline(v = 0.04, col="red")

#use a fixed mixture model to do this using mixtools on number of genes in the drops of the layer
log_gene_size <- log10(center@meta.data$nGene)
#gene_rank <- rank(-center@meta.data$nGene)

#fix a mixture model using mixtools
set.seed(-92497)
require("mixtools")
mix_gene <- normalmixEM(log_gene_size)
plot(mix_gene, breaks = 100, which=2, xlab2="log(genes per cell)") + abline(v = split_gene, col="red")
#split defined below
p1_gene <- dnorm(log_gene_size, mean=mix_gene$mu[1], sd=mix_gene$sigma[1])
p2_gene <- dnorm(log_gene_size, mean=mix_gene$mu[2], sd=mix_gene$sigma[2])
if (mix_gene$mu[1] < mix_gene$mu[2]) {
  split_gene <- min(log_gene_size[p2_gene > p1_gene])
} else {
  split_gene <- min(log_gene_size[p1_gene > p2_gene])
}

threshold_gene <- 10^split_gene

#refilter now focusing on gene filter
center <- FilterCells(object = center, subset.names = c("nGene", "percent.mito", "nUMI"), low.thresholds = c(398, -Inf, 489), high.thresholds = c(3000, 0.3, 10000))



#filter middle cells ----------------------------------------------

##MIDDLE
# calculate the percentage of mitochondrial genes and store it in percent.mito using AddMetaData.
mito.genes <- grep(pattern = "^MT-", x = rownames(x = middle@data), value = TRUE)
#percent.mito <- colSums(middle@raw.data[mito.genes, ]) / colSums(middle@raw.data)
percent.mito <- Matrix::colSums(middle@raw.data[mito.genes, ])/Matrix::colSums(middle@raw.data)
# AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
middle <- AddMetaData(object = middle, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = middle, point.size.use = .01, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

#plot genes/mito and UMI/Genes
par(mfrow = c(1, 2))
GenePlot(object = middle, gene1 = "nUMI", gene2 = "percent.mito", pch.use = ".")
GenePlot(object = middle, gene1 = "nUMI", gene2 = "nGene", pch.use = ".")

#check out gene and UMI abundance
#UMI
hist(middle@meta.data$nUMI, breaks = 2500, xlim = c(0,10000), main = "UMI's per cell") 
hist(middle@meta.data$nGene, breaks = 440, xlim = c(0,5000), main = "genes per cell")
hist(middle@meta.data$percent.mito, breaks = 150, xlim = c(0,0.5), main = "fraction of mitochondrial reads per cell")

#find real cells vs background
barcode_rank <- rank(-middle@meta.data$nUMI)
plot(barcode_rank, middle@meta.data$nUMI, xlim=c(1,8000))

#make logoritmic
log_lib_size <- log10(middle@meta.data$nUMI)
plot(barcode_rank, log_lib_size, xlim=c(1,8000))

#find inflection point
o <- order(barcode_rank)
log_lib_size <- log_lib_size[o]
barcode_rank <- barcode_rank[o]

rawdiff <- diff(log_lib_size)/diff(barcode_rank)
inflection <- which(rawdiff == min(rawdiff[100:length(rawdiff)], na.rm=TRUE))

plot(barcode_rank, log_lib_size)
abline(v=inflection, col="red", lwd=2)
#SEEMS too high


#CellRanger (10X) technique
n_cells <- middle@data@Dim[2]
totals <- middle@meta.data$nUMI
totals <- sort(totals, decreasing = TRUE)
# 99th percentile of top n_cells divided by 10
thresh = totals[round(0.01*n_cells)]/10
plot(totals)
abline(h=thresh, col="red", lwd=2)


#fix a mixture model using mixtools
set.seed(-92497)
require("mixtools")
mix <- normalmixEM(log_lib_size)
plot(mix, breaks = 100, which=2, xlab2="log(mol per cell)") + abline(v = split, col="red")
#split defined below
p1 <- dnorm(log_lib_size, mean=mix$mu[1], sd=mix$sigma[1])
p2 <- dnorm(log_lib_size, mean=mix$mu[2], sd=mix$sigma[2])
if (mix$mu[1] < mix$mu[2]) {
  split <- min(log_lib_size[p2 > p1])
} else {
  split <- min(log_lib_size[p1 > p2])
}
#identify cells using this split point and calculate the TPR and Recall.

threshold <- 10^split
umi_per_barcode <- data.frame(middle@cell.names, middle@meta.data$nUMI)
cells <- umi_per_barcode[umi_per_barcode[,2] > threshold,1]
TPR <- sum(cells %in% middle@cell.names)/length(cells)
Recall <- sum(cells %in% middle@cell.names)/length(middle@cell.names)
c(TPR, Recall)

par(mfrow = c(1, 2))
GenePlot(object = middle, gene1 = "nUMI", gene2 = "percent.mito", pch.use = ".")
GenePlot(object = middle, gene1 = "nUMI", gene2 = "nGene", pch.use = ".")
# We filter out cells that have unique gene counts over 4,000 or less than 200
middle <- FilterCells(object = middle, subset.names = c("nGene", "percent.mito", "nUMI"), low.thresholds = c(5, -Inf, 1000), high.thresholds = c(10000, 0.4, Inf))

#after filtering for UMI, the bimodality of the gene expression becomes visible. At this point, filter genes removing drops with an obvious lower number of genes expressed.
hist(middle@meta.data$nUMI, breaks = 2500, xlim = c(0,15000), main = "UMI's per cell") + abline(v = 1000, col="red")
hist(middle@meta.data$nGene, breaks = 440, xlim = c(0,5000), main = "genes per cell")
hist(middle@meta.data$percent.mito, breaks = 150, xlim = c(0,0.5), main = "fraction of mitochondrial reads per cell") + abline(v = 0.055, col="red")

#use a fixed mixture model to do this using mixtools on number of genes in the drops of the layer
log_gene_size <- log10(middle@meta.data$nGene)
#gene_rank <- rank(-middle@meta.data$nGene)

#fix a mixture model using mixtools
set.seed(-92497)
require("mixtools")
mix_gene <- normalmixEM(log_gene_size)
plot(mix_gene, breaks = 100, which=2, xlab2="log(genes per cell)") + abline(v = split_gene, col="red")
#split defined below
p1_gene <- dnorm(log_gene_size, mean=mix_gene$mu[1], sd=mix_gene$sigma[1])
p2_gene <- dnorm(log_gene_size, mean=mix_gene$mu[2], sd=mix_gene$sigma[2])
if (mix_gene$mu[1] < mix_gene$mu[2]) {
  split_gene <- min(log_gene_size[p2_gene > p1_gene])
} else {
  split_gene <- min(log_gene_size[p1_gene > p2_gene])
}

threshold_gene <- 10^split_gene

#refilter now focusing on gene filter
middle <- FilterCells(object = middle, subset.names = c("nGene", "percent.mito", "nUMI"), low.thresholds = c(589, 0.055, 478), high.thresholds = c(5000, 0.3, 25000))




#filter surface cells ---------------------------

##SURFACE
# calculate the percentage of mitochondrial genes and store it in percent.mito using AddMetaData.
mito.genes <- grep(pattern = "^MT-", x = rownames(x = surface@data), value = TRUE)
#percent.mito <- colSums(surface@raw.data[mito.genes, ]) / colSums(surface@raw.data)
percent.mito <- Matrix::colSums(surface@raw.data[mito.genes, ])/Matrix::colSums(surface@raw.data)
# AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
surface <- AddMetaData(object = surface, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = surface, point.size.use = .01, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

#plot genes/mito and UMI/Genes
par(mfrow = c(1, 2))
GenePlot(object = surface, gene1 = "nUMI", gene2 = "percent.mito", pch.use = ".")
GenePlot(object = surface, gene1 = "nUMI", gene2 = "nGene", pch.use = ".")

#check out gene and UMI abundance
#UMI
hist(surface@meta.data$nUMI, breaks = 200, xlim = c(0,15000))
hist(surface@meta.data$nGene, breaks = 250, xlim = c(0,4000)) + abline(v = threshold_gene, col="red")
hist(surface@meta.data$percent.mito, breaks = 100) + abline(v = 0.05, col="red")

#find real cells vs background
barcode_rank <- rank(-surface@meta.data$nUMI)
plot(barcode_rank, surface@meta.data$nUMI, xlim=c(1,8000))

#make logoritmic
log_lib_size <- log10(surface@meta.data$nUMI)
plot(barcode_rank, log_lib_size, xlim=c(1,8000))

#find inflection point
o <- order(barcode_rank)
log_lib_size <- log_lib_size[o]
barcode_rank <- barcode_rank[o]

rawdiff <- diff(log_lib_size)/diff(barcode_rank)
inflection <- which(rawdiff == min(rawdiff[100:length(rawdiff)], na.rm=TRUE))

plot(barcode_rank, log_lib_size)
abline(v=inflection, col="red", lwd=2)
#SEEMS too high

#CellRanger (10X) technique
n_cells <- surface@data@Dim[2]
totals <- surface@meta.data$nUMI
totals <- sort(totals, decreasing = TRUE)
# 99th percentile of top n_cells divided by 10
thresh = totals[round(0.01*n_cells)]/10
plot(totals)
abline(h=thresh, col="red", lwd=2)

#fix a mixture model using mixtools
set.seed(-92497)
require("mixtools")
mix <- normalmixEM(log_lib_size)
plot(mix, breaks = 100, which=2, xlab2="log(mol per cell)") + abline(v = split, col="red")
#split defined below
p1 <- dnorm(log_lib_size, mean=mix$mu[1], sd=mix$sigma[1])
p2 <- dnorm(log_lib_size, mean=mix$mu[2], sd=mix$sigma[2])
if (mix$mu[1] < mix$mu[2]) {
  split <- min(log_lib_size[p2 > p1])
} else {
  split <- min(log_lib_size[p1 > p2])
}
#identify cells using this split point and calculate the TPR and Recall.

threshold <- 10^split
umi_per_barcode <- data.frame(surface@cell.names, surface@meta.data$nUMI)
cells <- umi_per_barcode[umi_per_barcode[,2] > threshold,1]
TPR <- sum(cells %in% surface@cell.names)/length(cells)
Recall <- sum(cells %in% surface@cell.names)/length(surface@cell.names)
c(TPR, Recall)


par(mfrow = c(1, 2))
GenePlot(object = surface, gene1 = "nUMI", gene2 = "percent.mito", pch.use = ".")
GenePlot(object = surface, gene1 = "nUMI", gene2 = "nGene", pch.use = ".")
# We filter out cells that have unique gene counts over 4,000 or less than 200

surface <- FilterCells(object = surface, subset.names = c("nGene", "percent.mito", "nUMI"), low.thresholds = c(5, -Inf, 610), high.thresholds = c(6000, 0., 20000))


#after filtering for UMI, the bimodality of the gene expression becomes visible. At this point, filter genes removing drops with an obvious lower number of genes expressed.
#use a fixed mixture model to do this using mixtools on number of genes in the drops of the layer
log_gene_size <- log10(surface@meta.data$nGene)
#gene_rank <- rank(-surface@meta.data$nGene)

#fix a mixture model using mixtools
set.seed(-92497)
require("mixtools")
mix_gene <- normalmixEM(log_gene_size)
plot(mix_gene, breaks = 100, which=2, xlab2="log(genes per cell)") + abline(v = split_gene, col="red")
#split defined below
p1_gene <- dnorm(log_gene_size, mean=mix_gene$mu[1], sd=mix_gene$sigma[1])
p2_gene <- dnorm(log_gene_size, mean=mix_gene$mu[2], sd=mix_gene$sigma[2])
if (mix_gene$mu[1] < mix_gene$mu[2]) {
  split_gene <- min(log_gene_size[p2_gene > p1_gene])
} else {
  split_gene <- min(log_gene_size[p1_gene > p2_gene])
}

threshold_gene <- 10^split_gene

#refilter now focusing on gene filter
surface <- FilterCells(object = surface, subset.names = c("nGene", "percent.mito", "nUMI"), low.thresholds = c(592, 0.05, 610), high.thresholds = c(5000, 0.3, 12000))



#filtering summary--------

# calculate the percentage of mitochondrial genes and store it in percent.mito using AddMetaData.
#center
mito.genes <- grep(pattern = "^MT-", x = rownames(x = center@data), value = TRUE)
percent.mito <- Matrix::colSums(center@raw.data[mito.genes, ])/Matrix::colSums(center@raw.data)
center <- AddMetaData(object = center, metadata = percent.mito, col.name = "percent.mito")
#middle
mito.genes <- grep(pattern = "^MT-", x = rownames(x = middle@data), value = TRUE)
percent.mito <- Matrix::colSums(middle@raw.data[mito.genes, ])/Matrix::colSums(middle@raw.data)
middle <- AddMetaData(object = middle, metadata = percent.mito, col.name = "percent.mito")
#surface
mito.genes <- grep(pattern = "^MT-", x = rownames(x = surface@data), value = TRUE)
percent.mito <- Matrix::colSums(surface@raw.data[mito.genes, ])/Matrix::colSums(surface@raw.data)
surface <- AddMetaData(object = surface, metadata = percent.mito, col.name = "percent.mito")

#Filter
center <- FilterCells(object = center, subset.names = c("nGene", "percent.mito", "nUMI"), low.thresholds = c(580, 0.045, 489), high.thresholds = c(5000, 0.3, 15000))
middle <- FilterCells(object = middle, subset.names = c("nGene", "percent.mito", "nUMI"), low.thresholds = c(589, 0.055, 478), high.thresholds = c(5000, 0.3, 15000))
surface <- FilterCells(object = surface, subset.names = c("nGene", "percent.mito", "nUMI"), low.thresholds = c(592, 0.05, 610), high.thresholds = c(5000, 0.3, 15000))



#merge layers --------------------------------------------------------
#MERGE
# merge datasets organoid2 = merged seurat object (temp_merge is a placeholder to combine all 3 objects)
temp_merg <- MergeSeurat(object1 = center, object2 = middle, project = "ORGANOID")
organoid2 <- MergeSeurat(object1 = temp_merg, object2 = surface, project = "ORGANOID")

# calculate the percentage of mitochondrial genes and store it in percent.mito using AddMetaData.
mito.genes <- grep(pattern = "^MT-", x = rownames(x = organoid2@data), value = TRUE)
#percent.mito <- colSums(organoid2@raw.data[mito.genes, ]) / colSums(organoid2@raw.data)
percent.mito <- Matrix::colSums(organoid2@raw.data[mito.genes, ])/Matrix::colSums(organoid2@raw.data)
# AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
organoid2 <- AddMetaData(object = organoid2, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = organoid2, point.size.use = .01, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

#scatterplot
par(mfrow = c(1, 2))
GenePlot(object = organoid2, gene1 = "nUMI", gene2 = "percent.mito", pch.use = ".")
GenePlot(object = organoid2, gene1 = "nUMI", gene2 = "nGene", pch.use = ".")

# combined filter if desired
#organoid2 <- FilterCells(object = organoid2, subset.names = c("nGene", "percent.mito", "nUMI"), low.thresholds = c(500, -Inf, 300), high.thresholds = c(4000, 0.375, Inf))

# normalize
organoid2 <- NormalizeData(object = organoid2, normalization.method = "LogNormalize", scale.factor = 1e4)

#Detection of variable genes
organoid2 <- FindVariableGenes(object = organoid2, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

dev.off() #resets graphing window layout so that the variable gene plot takes up the width of the entire window
VariableGenePlot(organoid2, do.text = TRUE,  do.contour = FALSE)
length(x = organoid2@var.genes)

#regress out cell-cell variation based on the number of detected molecules per cell as well as the percentage mitochondrial gene content.  
organoid2 <- ScaleData(object = organoid2, vars.to.regress = c("nUMI", "percent.mito"))

# Perform linear dimensional reduction (PCA)
organoid2 <- RunPCA(object = organoid2, pc.genes = organoid2@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 10)

# visualizing both cells and genes that define the PCA
PrintPCA(object = organoid2, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
VizPCA(object = organoid2, pcs.use = 1:2)
PCAPlot(object = organoid2, dim.1 = 1, dim.2 = 2)

# ProjectPCA scores each gene (including genes not included in the PCA) based on their correlation with the calculated components
organoid2 <- ProjectPCA(object = organoid2, do.print = FALSE)

#`PCHeatmap`can be useful when trying to decide which PCs to include for further downstream analyses.
PCHeatmap(object = organoid2, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, remove.key = TRUE)

# multi-heatmap
PCHeatmap(object = organoid2, pc.use = 1:12, cells.use = 200, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)

#Determine statistically significant principal components
#organoid2 <- JackStraw(object = organoid2, num.replicate = 100, do.print = FALSE)
#JackStrawPlot(object = organoid2, PCs = 1:12)

# PC Elbow plot
PCElbowPlot(object = organoid2)

# Cluster the cells
# save.SNN = T saves the SNN so that the clustering algorithm can be rerun using the same graph but with a different resolution value)
organoid2 <- FindClusters(object = organoid2, reduction.type = "pca", dims.use = 1:5, resolution = 0.6, print.output = 0, save.SNN = TRUE)

#`PrintFindClustersParams` prints a summary of the parameters that were chosen. 
PrintFindClustersParams(object = organoid2)

# Run Non-linear dimensional reduction (tSNE)
organoid2 <- RunTSNE(object = organoid2, dims.use = 1:5, do.fast = TRUE)

# note that you can set do.label=T to help label individual clusters
TSNEPlot(object = organoid2)
TSNEPlot(object = organoid2, group.by = "orig.ident")

#save the object at this point so that it can easily be loaded back without having to rerun the computationally intensive steps performed above
saveRDS(organoid2, file = "R_scripts/organoids_comb.rds")


##################Save and load RDS file
organoid2 <- readRDS("organoids_comb.rds")

#Finding differentially expressed genes (cluster biomarkers)
# find all markers of cluster 1
cluster1.markers <- FindMarkers(object = organoid2, ident.1 = 1, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 5))
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(object = organoid2, ident.1 = 5, ident.2 = c(0,3), min.pct = 0.25)
print(x = head(x = cluster5.markers, n = 5))
# find markers for every cluster compared to all remaining cells, report only the positive ones
organoid2.markers <- FindAllMarkers(object = organoid2, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
organoid2.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)


#Jan 6, 2018: attempting to set identity of cells acording to original identity (layers), and then find marker genes for each layer ---------
organoid2 = SetAllIdent(organoid2, "orig.ident")
#reset identity to cluster identity
organoid2 = SetAllIdent(organoid2, "res.0.6")
layers.markers <- FindAllMarkers(object = organoid2, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
# take a peak
head(layers.markers)
top10logFC <- layers.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(organoid2, genes.use = top10logFC$gene, slim.col.label = TRUE, remove.key = TRUE)


# select markers to use for a heatmap (positive markers with high discriminatory power)
# avg_diff doesnt seem to work here: markers.use=subset(layers.markers, avg_diff >0&power>0.8)$gene 
markers.use=subset(layers.markers, avg_logFC > 0 & power > 0.8)$gene # power function doesn't work and only displays C and O layers - maybe I and S are too weak?

#try new marker test:
layers.markers.poisson <- FindAllMarkers(organoid2, logfc.threshold = 0.1,test.use = "poisson", do.print = TRUE)
# take a peak
head(layers.markers.poisson)

top10logFC <- layers.markers.poisson %>% group_by(cluster) %>% top_n(10, avg_logFC)
top10Pval <- layers.markers.tobit %>% group_by(cluster) %>% top_n(10, p_val)
top10AdjPval <- layers.markers %>% group_by(cluster) %>% top_n(10, p_val_adj)

DoHeatmap(organoid2, genes.use = top10logFC$gene, group.by = organoid2@meta.data$orig.ident, slim.col.label = TRUE, remove.key = TRUE)
DoHeatmap(organoid2, genes.use = top10Pval$gene, slim.col.label = TRUE, remove.key = TRUE)
DoHeatmap(organoid2, genes.use = top10AdjPval$gene, slim.col.label = TRUE, remove.key = TRUE)

#again
layers.markers <- FindAllMarkers(object = organoid2, only.pos = TRUE, min.pct = 0.2, thresh.use = 0.2)
#markers.use=subset(layers.markers, avg_logFC > 0 & power > 0.8)

top10 <- layers.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(organoid2, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)

#manually find genes that differ between differnt layers
clusterC.markers <- FindMarkers(object = organoid2, ident.1 = "C", ident.2 = c("I", "O", "S"), 
                                min.pct = 0.25)
print(x = head(x = clusterC.markers, n = 5))

clusterI.markers <- FindMarkers(object = organoid2, ident.1 = "I", ident.2 = c("C", "O", "S"), 
                                min.pct = 0.25)
print(x = head(x = clusterI.markers, n = 5))

clusterO.markers <- FindMarkers(object = organoid2, ident.1 = "O", ident.2 = c("C", "I", "S"), 
                                min.pct = 0.25)
print(x = head(x = clusterO.markers, n = 5))

clusterS.markers <- FindMarkers(object = organoid2, ident.1 = "S", ident.2 = c("C", "I", "O"), 
                                min.pct = 0.25)
print(x = head(x = clusterS.markers, n = 5))

# TRY TO MOVE THE GENE NAMES FROM THE ROW NAME TO A COLUMN, THEN COMBINE
clusterC.markers <- setDT(clusterC.markers, keep.rownames = TRUE)[]
names(clusterC.markers)[names(clusterC.markers) == 'rn'] <- 'gene'

clusterI.markers <- setDT(clusterI.markers, keep.rownames = TRUE)[]
names(clusterI.markers)[names(clusterI.markers) == 'rn'] <- 'gene'

clusterO.markers <- setDT(clusterO.markers, keep.rownames = TRUE)[]
names(clusterO.markers)[names(clusterO.markers) == 'rn'] <- 'gene'

clusterS.markers <- setDT(clusterS.markers, keep.rownames = TRUE)[]
names(clusterS.markers)[names(clusterS.markers) == 'rn'] <- 'gene'

#try to combine these lists of genes
man.top.genes <- dplyr::bind_rows(list(C = clusterC.markers, I = clusterI.markers, O = clusterO.markers, S = clusterS.markers), .id = 'cluster')

# choose only top10 genes manually from each layer
# sorted by avg logFC here, try pvalue next maybe
man.top10 <- man.top.genes %>% group_by(cluster) %>% top_n(10, avg_logFC) 
DoHeatmap(organoid2, genes.use = man.top10$gene, draw.line = TRUE, slim.col.label = TRUE, remove.key = TRUE)



# Seurat has four tests for differential expression which can be set with the test.use parameter: ROC test ("roc"), t-test ("t"), LRT test based on zero-inflated data ("bimod", default), LRT test based on tobit-censoring models ("tobit") The ROC test returns the 'classification power' for any individual marker (ranging from 0 - random, to 1 - perfect).

# markersroc, fig.height=8, fig.width=15,
cluster1.markers <- FindMarkers(object = organoid2, ident.1 = 0, thresh.use = 0.25, test.use = "roc", only.pos = TRUE)

# We include several tools for visualizing marker expression. `VlnPlot` (shows expression probability distributions across clusters), and `FeaturePlot` (visualizes gene expression on a tSNE or PCA plot) are our most commonly used visualizations. We also suggest exploring `JoyPlot`, `CellPlot`, and `DotPlot` as additional methods to view your dataset.

# markerplots, fig.height=8, fig.width=15,
VlnPlot(object = organoid2, features.plot = c("MS4A1", "CD79A"))

# you can plot raw UMI counts as well
VlnPlot(object = organoid2, features.plot = c("NKG7", "PF4"), use.raw = TRUE, y.log = TRUE)

#surface
FeaturePlot(object = organoid2, features.plot = c("SOD2", "IFITM1", "NFKBIA", "RSAD2", "HLA-B", "TXNIP", "HLA-B", "IFI27", "PSME1"), cols.use = c("grey", "red"), reduction.use = "tsne")

#surface
FeaturePlot(object = organoid2, features.plot = c("MT2A", "BTG1", "NCOA7", "PFKP", "LY6E", "SAMD9L", "EIF4E3", "NAMPT", "SSPN"), cols.use = c("grey", "red"), reduction.use = "tsne")

#CLUSTER1
FeaturePlot(object = organoid2, features.plot = c("HIST1H4C", "S100A6", "RPS3A", "GOLGA4", "RPL39", "EIF3E"), cols.use = c("grey", "red"), reduction.use = "tsne")



# `DoHeatmap` generates an expression heatmap for given cells and genes. In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.

# clusterHeatmap, fig.height=8, fig.width=15, message=FALSE, warning=FALSE
top10 <- organoid2.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
DoHeatmap(object = organoid2, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)








#further SEURAT analysis----------


# labelplot, fig.height=5, fig.width=9, warning = FALSE
current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", "FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Megakaryocytes")
organoid2@ident <- plyr::mapvalues(x = organoid2@ident, from = current.cluster.ids, to = new.cluster.ids)
TSNEPlot(object = organoid2, do.label = TRUE, pt.size = 0.5)

### Further subdivisions within cell types

# First lets stash our identities for later
organoid2 <- StashIdent(object = organoid2, save.name = "ClusterNames_0.6")

# Note that if you set save.snn=T above, you don't need to recalculate the SNN, and can simply put: 
# organoid2 <- FindClusters(organoid2,resolution = 0.8)
organoid2 <- FindClusters(object = organoid2, reduction.type = "pca", dims.use = 1:10, resolution = 0.8, print.output = FALSE)

# Demonstration of how to plot two tSNE plots side by side, and how to color points based on different criteria
plot1 <- TSNEPlot(object = organoid2, do.return = TRUE, no.legend = TRUE, do.label = TRUE)
plot2 <- TSNEPlot(object = organoid2, do.return = TRUE, group.by = "ClusterNames_0.6", no.legend = TRUE, do.label = TRUE)
plot_grid(plot1, plot2)

# Find discriminating markers
tcell.markers <- FindMarkers(object = organoid2, ident.1 = 0, ident.2 = 1)

# Most of the markers tend to be expressed in C1 (i.e. S100A4). However, we can see that CCR7 is upregulated in 
# C0, strongly indicating that we can differentiate memory from naive CD4 cells.
# cols.use demarcates the color palette from low to high expression
FeaturePlot(object = organoid2, features.plot = c("S100A4", "CCR7"), cols.use = c("green", "blue"))


The memory/naive split is bit weak, and we would probably benefit from looking at more cells to see if this becomes more convincing. In the meantime, we can restore our old cluster identities for downstream processing.

# restore
organoid2 <- SetAllIdent(object = organoid2, id = "ClusterNames_0.6")
save(organoid2, file = "~/Projects/datasets/organoid23k_final.Rda")







