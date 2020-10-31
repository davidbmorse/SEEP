#re-analysis on biopsy after subtracting mouse reads
library(Seurat)
library(dplyr)
library(Matrix)
library(Matrix.utils)
library(data.table)
library(ggplot2)
library(DropletUtils)
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/")

#import human mouse reads as sparse matrix
human_mouse <- Read10X("R_scripts/4th_MERGE_mouses/SubtractMouseReads/mouse_human_counts/output/")
#filtering for this dataset was: 
#PDX_biopsy2_3 <- subset(x = PDX_biopsy2_3, subset = percent.mt < 40 & nCount_RNA > 200 & nCount_RNA < 40000 & nFeature_RNA < 8000)

#remove mouse genes
human_genes_only <- human_mouse[1:23429,]
biopsy <- CreateSeuratObject(human_genes_only, min.cells = 1, assay = "RNA")


# differentiate cells
gene_species <- rownames(human_mouse)
human_inds <- gene_species == gene_species[1:23429]
mouse_inds <- !human_inds

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

# summary of mouse v human cells
cell_species %>% 
  dplyr::count(species) %>% 
  mutate(proportion = n / ncol(human_mouse))

# add species identity to seurat object
biopsy <- AddMetaData(biopsy, metadata = cell_species$species, col.name = "species")
VlnPlot(biopsy, c("nCount_RNA", "nFeature_RNA"), group.by = "species",
        pt.size = 0.1)
ggplot(biopsy@meta.data, aes(nCount_RNA, nFeature_RNA, color = species)) +
  geom_point(alpha = 0.7, size = 0.5) +
  labs(x = "Total UMI counts per cell", y = "Number of genes detected")

#remove mouse and mixed cells from biopsy
biopsy <- subset(biopsy, subset = species == "human")

saveRDS(biopsy, file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/raw_filtered_biopsy_seurat.rds")







