rm(list=ls())
setwd("C://PROJECTS/P2022/SEEP_Manuscript")

library(Seurat)

# Read in data sources - Seurat objects v3
sph=readRDS("Input/spheroid_CIOS_AllGenes.rds")
org=readRDS("Input/organoids_comb_v3.rds")
bio= readRDS("Input/PDX_biopsy2_3_MERGE_final.rds")

# Update Seurat objects to Seurat v4
sph = UpdateSeuratObject(sph)
org = UpdateSeuratObject(org)
bio = UpdateSeuratObject(bio)

# Unify/find common metadata columns across data sources
org$percent.mt = org$percent.mito
org$seurat_clusters = org$res.0.6
cols=intersect(intersect(colnames(sph@meta.data),colnames(org@meta.data)),colnames(bio@meta.data))

# Merge data sources into one Seurat object
obj=merge(sph, y = c(org, bio), add.cell.ids = c("sph", "org", "bio"), project = "SEEP")

# Remove middle layers; add/modify data sources metadata
for(i in setdiff(colnames(obj@meta.data), cols)) obj[[i]] = NULL
obj$source = substring(rownames(obj@meta.data),1,3)
obj$orig.ident = toupper(obj$orig.ident)
obj$layer <- ifelse(obj$orig.ident %in% c("I","O"), "M", paste(obj$orig.ident))
obj$layer_source = paste(obj$layer,obj$source,sep="_")
Idents(obj) = obj$layer_source
so <- subset(obj, subset=orig.ident %in% c("C","S"))
so$long_layer <- factor(ifelse(so$layer == "C", "Center", "Surface"))
so$long_source <- factor(so$source)
levels(so$long_source) <- c("Biopsy","Organoid","Spheroid")
so$seurat_clusters <- NULL

# Save data
seep = so
save(seep, file="data/SEEP_Merge.rda")
