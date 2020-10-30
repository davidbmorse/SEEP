#re-analysis on biopsy after subtracting mouse reads
library(Seurat)
library(dplyr)
library(Matrix)
library(Matrix.utils)
library(data.table)
library(ggplot2)
library(DropletUtils)
library(patchwork)
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/")

biopsy <- readRDS(file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/raw_filtered_biopsy_seurat.rds")

#write.table(as.matrix(GetAssayData(object = biopsy, slot = "counts")), 'counts.tsv', sep = '\t', row.names = T, col.names = T, quote = F)

biopsy[["percent.mt"]] <- PercentageFeatureSet(biopsy, pattern = "^MT-")

VlnPlot(biopsy, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(biopsy, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(biopsy, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

#split seurat, filter, join
center <- subset(biopsy, subset = orig.ident == "center")
middle <- subset(biopsy, subset = orig.ident == "middle")
surface <- subset(biopsy, subset = orig.ident == "surface")
                                    #NOTE - CONSIDER RASING THRESHOLDS BELOW
center <- subset(x = center, subset = percent.mt < 40 & nFeature_RNA < 2000 & nCount_RNA > 250 & nCount_RNA < 10000)
middle <- subset(x = middle, subset = percent.mt < 40 & nFeature_RNA < 3000 & nCount_RNA > 350 & nCount_RNA < 15000)
surface <- subset(x = surface, subset = percent.mt < 40 & nFeature_RNA < 4000 & nCount_RNA > 400 & nCount_RNA < 15000)

biopsy <-merge(x = center, y = c(middle, surface))
biopsy$orig.ident <- Idents(biopsy)
#saveRDS(biopsy, file="filtered_biopsy_SO.rds")
#biopsy <- readRDS(file = "filtered_biopsy_SO.rds")

#POST FILTER - normalization, scaling, and analysis
#cell cycle identity from Tirosh et al, 2015
cc.genes <- readLines(con = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

#normalize
biopsy <- NormalizeData(biopsy)
biopsy = FindVariableFeatures(biopsy, verbose = TRUE, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))

#find variable
top10 <- head(VariableFeatures(biopsy), 10)
plot1 <- VariableFeaturePlot(biopsy)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scale
all.genes <- rownames(biopsy)
biopsy <- ScaleData(biopsy, features = all.genes)

#PCA
biopsy <- RunPCA(biopsy, features = VariableFeatures(biopsy))
biopsy <- JackStraw(biopsy, num.replicate = 100, dims = 50)
biopsy <- ScoreJackStraw(biopsy, dims = 1:50)
JackStrawPlot(biopsy, dims = 1:20)
#export as 4x6 pdf
ElbowPlot(biopsy)
#export as 3X3 pdf

#clustering
biopsy <- FindNeighbors(biopsy, dims = 1:10)
biopsy <- FindClusters(biopsy, resolution = 0.5)

#non-linear dimenstional reductions
biopsy <- RunUMAP(biopsy, dims = 1:10)
DimPlot(biopsy, reduction = "umap")
DimPlot(biopsy, reduction = "umap", group.by = "orig.ident")

biopsy <- RunTSNE(biopsy, dims = 1:10, perplexity=80, exaggeration_factor=20)
DimPlot(biopsy, reduction = "tsne")
DimPlot(biopsy, reduction = "tsne", group.by = "orig.ident")

#find differentially expressed genes
biopsy.markers <- FindAllMarkers(biopsy, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.20)
biopsy.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#change idents
Idents(biopsy) <- "orig.ident"
Idents(biopsy) <- "seurat_clusters"

biopsy.markers.layers <- FindAllMarkers(biopsy, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
biopsy.markers.layers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#plot for clusters
top20 <- biopsy.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
DoHeatmap(biopsy, features = top20$gene) + NoLegend()

#plot for layers
top20.layers <- biopsy.markers.layers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
DoHeatmap(biopsy, features = top20.layers$gene) + NoLegend()

pdf("figures/Hmap_full_biop.pdf",  width = 4, height = 5) # width = 2.1, height = 4.3)
Hmap
dev.off()




#save the gene lists as csv files
write.csv(biopsy.markers.layers, file = "layer_markers_025_025.csv")
write.csv(biopsy.markers, file = "cluster_markers_010_020.csv")

#space

#save finalized (clustered and dimensionally reduced)
saveRDS(biopsy, file="/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/filtered_biopsy_SO_final_042220.rds")









