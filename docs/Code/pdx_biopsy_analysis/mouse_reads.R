#merging mouse and human biopsy analysis
library(Seurat)
library(dplyr)
library(Matrix)
library(data.table) #for extra data table wrangling

setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/")
getwd()

# Load MOUSE ALLIGNMENT biopsy2_3 dataset -----
raw_center1 <- Read10X("counts_matrices/4th_MERGE_mouses/mouse_allignment/m1_MERGE/")
raw_middle1 <- Read10X("counts_matrices/4th_MERGE_mouses/mouse_allignment/m2_MERGE/")
raw_surface1 <- Read10X("counts_matrices/4th_MERGE_mouses/mouse_allignment/m3_MERGE/")
raw_center2 <- Read10X("counts_matrices/4th_MERGE_mouses/mouse_allignment/2m1_MERGE/")
raw_middle2 <- Read10X("counts_matrices/4th_MERGE_mouses/mouse_allignment/2m2_MERGE/")
raw_surface2 <- Read10X("counts_matrices/4th_MERGE_mouses/mouse_allignment/2m3_MERGE/")

#inititate seurat objects ------
PDX_C1 <- CreateSeuratObject(counts = raw_center1, min.cells = 1, assay = "RNA",names.field = 1,
                             names.delim = "_",
                             project = "PDX_1")
PDX_M1 <- CreateSeuratObject(counts = raw_middle1, min.cells = 1, assay = "RNA",names.field = 1,
                             names.delim = "_",
                             project = "PDX_1")
PDX_S1 <- CreateSeuratObject(counts = raw_surface1, min.cells = 1, assay = "RNA",names.field = 1,
                             names.delim = "_",
                             project = "PDX_1")

PDX_C2 <- CreateSeuratObject(counts = raw_center2, min.cells = 1, assay = "RNA",names.field = 1,
                             names.delim = "_",
                             project = "PDX_1")
PDX_M2 <- CreateSeuratObject(counts = raw_middle2, min.cells = 1, assay = "RNA",names.field = 1,
                             names.delim = "_",
                             project = "PDX_1")
PDX_S2 <- CreateSeuratObject(counts = raw_surface2, min.cells = 1, assay = "RNA",names.field = 1,
                             names.delim = "_",
                             project = "PDX_1")

#merge by layer before filtering
PDX_C <- merge(x = PDX_C1, y = c(PDX_C2), 
               add.cell.ids = c( "center", "center"))
PDX_M <- merge(x = PDX_M1, y = c(PDX_M2), 
               add.cell.ids = c( "middle", "middle"))
PDX_S <- merge(x = PDX_S1, y = c(PDX_S2), 
               add.cell.ids = c( "surface", "surface"))

#filter1---------
#do i need to do scrublet before this?
PDX_Cs <- SubsetByBarcodeInflections(CalculateBarcodeInflections(PDX_C))
PDX_Ms <- SubsetByBarcodeInflections(CalculateBarcodeInflections(PDX_M))
PDX_Ss <- SubsetByBarcodeInflections(CalculateBarcodeInflections(PDX_S))

PDX_Cs[["which_assay"]] <- "indrop"
PDX_Ms[["which_assay"]] <- "indrop"
PDX_Ss[["which_assay"]] <- "indrop"
PDX_Cs[["percent.mt"]] <- PercentageFeatureSet(PDX_Cs, pattern = "^mt-")
PDX_Ms[["percent.mt"]] <- PercentageFeatureSet(PDX_Ms, pattern = "^mt-")
PDX_Ss[["percent.mt"]] <- PercentageFeatureSet(PDX_Ss, pattern = "^mt-")

#take a peek -------
VlnPlot(PDX_Cs, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(PDX_Cs, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(PDX_Cs, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

VlnPlot(PDX_Ms, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(PDX_Ms, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(PDX_Ms, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

VlnPlot(PDX_Ss, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(PDX_Ss, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(PDX_Ss, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

PDX_Cs <- subset(x = PDX_Cs, subset = percent.mt < 10 & nCount_RNA > 300 & nCount_RNA < 100000)
PDX_Ms <- subset(x = PDX_Ms, subset = percent.mt < 10 & nCount_RNA > 300 & nCount_RNA < 100000)
PDX_Ss <- subset(x = PDX_Ss, subset = percent.mt < 10 & nCount_RNA > 300 & nCount_RNA < 100000)

#merge datasets-------
PDX_mouse <-merge(x = PDX_Cs, y = c(PDX_Ms,PDX_Ss), 
                      add.cell.ids = c( "center", "middle", "surface"))

#now, try to merge mouse and human datasets
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/")
PDX_biopsy2_3 <- readRDS(file = "R_exports/r_data/4th_merge_mouse/PDX_biopsy2_3_MERGE_final.rds")


xeno <- merge(x = PDX_biopsy2_3, y = PDX_mouse, add.cell.ids = NULL)

#this (above) didn't merge the cells - it created a new cell for different mouse/human reads

#SCTRANFORM ---------
library(sctransform)
xeno <- SCTransform(xeno, vars.to.regress = c("percent.mt"), verbose = FALSE)
#additional 'variables to regress': "S.Score","G2M.Score"

xeno <- RunPCA(xeno, verbose = FALSE)
xeno <- RunUMAP(xeno, dims = 1:50, verbose = FALSE)

xeno <- FindNeighbors(xeno, dims = 1:50, verbose = FALSE)
xeno <- FindClusters(xeno, verbose = FALSE)
DimPlot(xeno, label = TRUE) + NoLegend()
DimPlot(xeno, reduction = "umap", group.by = 'orig.ident')

xeno <- RunTSNE(xeno, dims = 1:50, verbose = FALSE)
DimPlot(xeno, reduction = "tsne", label = TRUE) + NoLegend()
DimPlot(xeno, reduction = "tsne", group.by = 'orig.ident')

#exploring gene expression------
xeno.markers <- FindAllMarkers(xeno, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
xeno.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top20 <- xeno.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
DoHeatmap(xeno, features = top20$gene) + NoLegend()



