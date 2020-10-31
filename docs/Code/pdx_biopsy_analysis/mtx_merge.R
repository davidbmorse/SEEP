#merging mouse and human reads
library(Seurat)
library(dplyr)
library(Matrix)
library(Matrix.utils)
library(data.table)
library(ggplot2)
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/")

#mouse genome cells
raw_center1_mouse <- Read10X("counts_matrices/4th_MERGE_mouses/mouse_allignment/m1_MERGE/")
#human genome cells
raw_center1_human <- Read10X("counts_matrices/4th_MERGE_mouses/m1_MERGE/")
#merge these by rBind - this will merge each cell, but will keep ALL genes
# [,1] gene names, [1,] cell names
h_m_join <- rBind.fill(raw_center1_mouse, raw_center1_human)


#load human allignment
raw_center1_h <- Read10X("counts_matrices/4th_MERGE_mouses/m1_MERGE/")
raw_middle1_h <- Read10X("counts_matrices/4th_MERGE_mouses/m2_MERGE/")
raw_surface1_h <- Read10X("counts_matrices/4th_MERGE_mouses/m3_MERGE/")
raw_center2_h <- Read10X("counts_matrices/4th_MERGE_mouses/2m1_MERGE/")
raw_middle2_h <- Read10X("counts_matrices/4th_MERGE_mouses/2m2_MERGE/")
raw_surface2_h <- Read10X("counts_matrices/4th_MERGE_mouses/2m3_MERGE/")

#load mouse allignment
raw_center1_m <- Read10X("counts_matrices/4th_MERGE_mouses/mouse_allignment/m1_MERGE/")
raw_middle1_m <- Read10X("counts_matrices/4th_MERGE_mouses/mouse_allignment/m2_MERGE/")
raw_surface1_m <- Read10X("counts_matrices/4th_MERGE_mouses/mouse_allignment/m3_MERGE/")
raw_center2_m <- Read10X("counts_matrices/4th_MERGE_mouses/mouse_allignment/2m1_MERGE/")
raw_middle2_m <- Read10X("counts_matrices/4th_MERGE_mouses/mouse_allignment/2m2_MERGE/")
raw_surface2_m <- Read10X("counts_matrices/4th_MERGE_mouses/mouse_allignment/2m3_MERGE/")

#join reads
raw_center1 <- rBind.fill(raw_center1_h, raw_center1_m)
raw_middle1 <- rBind.fill(raw_middle1_h, raw_middle1_m)
raw_surface1 <- rBind.fill(raw_surface1_h, raw_surface1_m)
raw_center2 <- rBind.fill(raw_center2_h, raw_center2_m)
raw_middle2 <- rBind.fill(raw_middle2_h, raw_middle2_m)
raw_surface2 <- rBind.fill(raw_surface2_h, raw_surface2_m)

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
PDX_C <- merge(x = PDX_C1, y = c(PDX_C2))
PDX_M <- merge(x = PDX_M1, y = c(PDX_M2))
PDX_S <- merge(x = PDX_S1, y = c(PDX_S2))

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

#take a peek -------
VlnPlot(PDX_biopsy2_3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(PDX_biopsy2_3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(PDX_biopsy2_3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

PDX_biopsy2_3 <- subset(x = PDX_biopsy2_3, subset = percent.mt < 40 & nCount_RNA > 200 & nCount_RNA < 40000 & nFeature_RNA < 8000)

#plot
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

PDX_biopsy2_3 <- RunTSNE(PDX_biopsy2_3, dims = 1:50, verbose = FALSE)
DimPlot(PDX_biopsy2_3, reduction = "tsne", label = TRUE) + NoLegend()
DimPlot(PDX_biopsy2_3, reduction = "tsne", group.by = 'orig.ident')

#exploring gene expression------
PDX_biopsy2_3.markers <- FindAllMarkers(PDX_biopsy2_3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
PDX_biopsy2_3.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)



###try to make a species mixing plot by selecting gene names with lowercase letters
#1 export all counts a sparse matrix using droplet utils
library(DropletUtils)
#write10xCounts(x = PDX_biopsy2_3@assays$RNA@counts, path = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/mouse_human_counts/output")

human_mouse <- Read10X("R_scripts/4th_MERGE_mouses/SubtractMouseReads/mouse_human_counts/output/")
gene_species <- rownames(human_mouse)
human_inds <- gene_species == gene_species[1:23429]
mouse_inds <- !human_inds
#mouse_inds <- gene_species == gene_species[23430:57138]
table(mouse_inds)["TRUE"]

# mark cells as mouse or human
cell_species <- tibble(n_mouse_umi = Matrix::colSums(human_mouse[mouse_inds,]),
                       n_human_umi = Matrix::colSums(human_mouse[human_inds,]),
                       tot_umi = Matrix::colSums(human_mouse),
                       prop_mouse = n_mouse_umi / tot_umi,
                       prop_human = n_human_umi / tot_umi)
# Classify species based on proportion of UMI, with cutoff of 90%
cell_species <- cell_species %>% 
  mutate(species = case_when(
    prop_mouse > 0.5 ~ "mouse",
    prop_human > 0.75 ~ "human",
    TRUE ~ "mixed"
  ))
ggplot(cell_species, aes(n_human_umi, n_mouse_umi, color = species)) +
  geom_point(size = 0.5)

#summary of mouse v human cells
cell_species %>% 
  dplyr::count(species) %>% 
  mutate(proportion = n / ncol(human_mouse))

# species     n     proportion
#   <chr>   <int>      <dbl>
# 1 human    2805    0.940  
# 2 mixed      29    0.00972
# 3 mouse     151    0.0506 

#explore
seu <- CreateSeuratObject(human_mouse, min.cells = 1, assay = "RNA") %>% 
  NormalizeData(verbose = FALSE) %>% 
  ScaleData(verbose = FALSE) %>% 
  FindVariableFeatures(verbose = FALSE)

# Add species to meta data
seu <- AddMetaData(seu, metadata = cell_species$species, col.name = "species")

VlnPlot(seu, c("nCount_RNA", "nFeature_RNA"), group.by = "species",
        pt.size = 0.1)

ggplot(seu@meta.data, aes(nCount_RNA, nFeature_RNA, color = species)) +
  geom_point(alpha = 0.7, size = 0.5) +
  labs(x = "Total UMI counts per cell", y = "Number of genes detected")

seu <- RunPCA(seu, verbose = FALSE, npcs = 50)
DimPlot(seu, reduction = "pca", pt.size = 0.5, group.by = "species")

seu <- RunTSNE(seu, dims = 1:20, check_duplicates = FALSE)
DimPlot(seu, reduction = "tsne", pt.size = 0.5, group.by = "species")


#add cell names (and edit them to match pdx object) to cell_species table
names_species <- tibble(row.names = colnames(human_mouse), cell_species)
names_species$row.names[1:426] <- paste0('center_', names_species$row.names[1:426])
names_species$row.names[427:1550] <- paste0('middle_', names_species$row.names[427:1550])
names_species$row.names[1551:2985] <- paste0('surface_', names_species$row.names[1551:2985])
names_species <- data.frame(row.names = names_species$row.names, species = names_species$cell_species$species)

#return to PDX_biopsy2_3 sample
#import data as RDS file (from figure 3 in paper) -----
PDX_biopsy2_3 <- readRDS(file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_exports/r_data/4th_merge_mouse/PDX_biopsy2_3_MERGE_final.rds")

PDX_biopsy2_3 <- AddMetaData(PDX_biopsy2_3, metadata = cell_species$species, col.name = "species")
PDX_biopsy2_3 <- AddMetaData(PDX_biopsy2_3, metadata = names_species, col.name = "species")

DimPlot(PDX_biopsy2_3, label = TRUE) + NoLegend()
DimPlot(PDX_biopsy2_3, reduction = "umap", group.by = 'orig.ident')
DimPlot(PDX_biopsy2_3, reduction = "umap", group.by = 'species')

DimPlot(PDX_biopsy2_3, reduction = "tsne", label = TRUE) + NoLegend()
DimPlot(PDX_biopsy2_3, reduction = "tsne", group.by = 'seurat_clusters')
DimPlot(PDX_biopsy2_3, reduction = "tsne", group.by = 'orig.ident')
DimPlot(PDX_biopsy2_3, reduction = "tsne", group.by = 'species')

test <- data.table(PDX_biopsy2_3@meta.data$species)
test$V1 <- as.character(test$V1)

test %>%
  dplyr::count(V1) %>% 
  mutate(proportion = n / nrow(test))

#remove mouse and mixed cells from biopsy
PDX_biopsy2_3 <- subset(PDX_biopsy2_3, subset = species == "human")

PDX_biopsy2_3 <- SetIdent(PDX_biopsy2_3, value = "clusterID")
PDX_biopsy2_3.markers <- FindAllMarkers(PDX_biopsy2_3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
PDX_biopsy2_3.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

PDX_biopsy2_3 <- SetIdent(PDX_biopsy2_3, value = "old.ident")
PDX_biopsy2_3.layers.markers <- FindAllMarkers(PDX_biopsy2_3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
PDX_biopsy2_3.layers.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)







                              