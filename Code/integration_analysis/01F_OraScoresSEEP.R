rm(list=ls())
setwd("C://PROJECTS/P2022/SEEP_Manuscript")

library(tidyverse); library(dplyr); library(Seurat); library(UCell)

# data
data("SEEP_ORA")
#data("SEEP_Regions")
data("SS2trans_Regions")

iso = isoSEEP_ref
DefaultAssay(iso) = "RNA"
iso$Regions = factor(iso$DA_regions_layer)
Idents(iso) = iso$Regions
colo2=ggsci::pal_d3()(4)
names(colo2) = levels(Idents(iso))

lpaths = lorasig_main
lgenes = lapply(lpaths, function(x) {
  y = data.frame(x)
  ly = y$overlapGenes
  names(ly) = y$pathway
  return(ly)
})
lgenes = lapply(lgenes, function(x){
  x[sapply(x, length) >= 5]
})
#lgenes_pos = lapply(lgenes, function(x) lapply(x, paste0, "+"))

so = iso
for(i in names(lgenes)){
so = AddModuleScore_UCell(so, features=lgenes[[i]])
colnames(so@meta.data) =
  ifelse(grepl("_UCell$", colnames(so@meta.data)),
         paste(colnames(so@meta.data), i, sep="."),
         colnames(so@meta.data))
}

mat = t(so@meta.data %>% select(contains(c("UCell"))))
scores = gsub("_", "-", rownames(mat))
rownames(mat) = scores

sof <- CreateSeuratObject(counts = mat, meta.data = so@meta.data)
sof@assays$RNA@var.features = rownames(sof)
sof <- ScaleData(object = sof, do.center=TRUE, do.scale=TRUE, scale.max = Inf)
coord = Embeddings(iso, "umap")
sof[["umap"]] <- CreateDimReducObject(embeddings = coord, key = "UMAP_", assay = "RNA")
coord = Embeddings(iso, "pca")
sof[["pca"]] <- CreateDimReducObject(embeddings = coord, key = "PC_", assay = "RNA")
Idents(sof) = sof@meta.data$Regions

iso2_paths = sof
iso2_scores = so
seep_pathwayScores = mat

save(iso2_paths, iso2_scores, seep_pathwayScores, file='data/SEEP_ORApaths.rda')

