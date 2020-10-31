library(Seurat)
library(dplyr)
library(Matrix)
library(data.table) #for extra data table wrangling

setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/")
getwd()

#Read in RDS files if necessary ------
PDX_biopsy2_3 <- readRDS(file = "R_exports/r_data/4th_merge_mouse/PDX_biopsy2_3_MERGE_final.rds")

# Load biopsy2_3 dataset -----
raw_center1 <- Read10X("counts_matrices/4th_MERGE_mouses/m1_MERGE/")
raw_middle1 <- Read10X("counts_matrices/4th_MERGE_mouses/m2_MERGE/")
raw_surface1 <- Read10X("counts_matrices/4th_MERGE_mouses/m3_MERGE/")
raw_center2 <- Read10X("counts_matrices/4th_MERGE_mouses/2m1_MERGE/")
raw_middle2 <- Read10X("counts_matrices/4th_MERGE_mouses/2m2_MERGE/")
raw_surface2 <- Read10X("counts_matrices/4th_MERGE_mouses/2m3_MERGE/")

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
PDX_Cs[["percent.mt"]] <- PercentageFeatureSet(PDX_Cs, pattern = "^MT-")
PDX_Ms[["percent.mt"]] <- PercentageFeatureSet(PDX_Ms, pattern = "^MT-")
PDX_Ss[["percent.mt"]] <- PercentageFeatureSet(PDX_Ss, pattern = "^MT-")

#NO FILTERING YET
PDX_biopsy2_3 <-merge(x = PDX_Cs, y = c(PDX_Ms,PDX_Ss), 
                      add.cell.ids = c( "center", "middle", "surface"))

PDX_biopsy2_3 <- subset(x = PDX_biopsy2_3, subset = percent.mt < 40 & nCount_RNA > 200 & nCount_RNA < 10000 & nFeature_RNA < 2000)
#cell cycle identity from Tirosh et al, 2015
cc.genes <- readLines(con = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/regev_lab_cell_cycle_genes.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
#PDX_biopsy2_3 <- CellCycleScoring(object = PDX_biopsy2_3 , s.features = s.genes, g2m.features = g2m.genes,
                                  set.ident = TRUE)

#SCTRANFORM ---------
library(sctransform)
PDX_biopsy2_3 <- PercentageFeatureSet(PDX_biopsy2_3, pattern = "^MT-", col.name = "percent.mt")
PDX_biopsy2_3 <- SCTransform(PDX_biopsy2_3, vars.to.regress = c("percent.mt"), verbose = FALSE)
#additional 'variables to regress': "S.Score","G2M.Score"

PDX_biopsy2_3 <- RunPCA(PDX_biopsy2_3, verbose = FALSE)
PDX_biopsy2_3 <- RunUMAP(PDX_biopsy2_3, dims = 1:50, verbose = FALSE)

PDX_biopsy2_3 <- FindNeighbors(PDX_biopsy2_3, dims = 1:50, verbose = FALSE)
PDX_biopsy2_3 <- FindClusters(PDX_biopsy2_3, verbose = FALSE)
DimPlot(PDX_biopsy2_3, label = TRUE) + NoLegend()
DimPlot(PDX_biopsy2_3, reduction = "umap", group.by = 'orig.ident')



