rm(list=ls())
setwd("C://PROJECTS/P2022/SEEP_Manuscript")

# data
load("data/SEEP_Integrate.rda")

library(Seurat); library(tidyverse);library(DAseq)

# run reductions ====
iso = isoSEEP

# PCA/UMAP
iso <- RunPCA(iso, npcs=50, verbose = FALSE)
iso <- RunUMAP(iso, reduction="pca", dims=1:20, return.model = TRUE)
iso3D <- RunUMAP(iso, reduction="pca", dims=1:20, return.model = TRUE, n.components=3)
pca = data.frame(Embeddings(iso, "pca"))[,1:20]
umap = data.frame(Embeddings(iso, "umap"))[,1:2]

# prep labels
label.info = data.frame(
label=levels(factor(iso$layer_source)),
condition=rep(c("Center",'Surface'), each=3))
labelsCenter <- label.info[label.info$condition == "Center", "label"]
labelsSurface <- label.info[label.info$condition == "Surface", "label"]

# DA scores
da_cells = getDAcells( X = pca,
cell.labels = iso$layer_source,
labels.1 = labelsCenter,
labels.2 = labelsSurface,
k.vector = seq(50, 500, 50),
plot.embedding = pca,
n.runs = 10,
n.rand = 5)

cutoff = 0.8
condition_1 = "Center"
condition_2 = "Surface"

da_cells <- updateDAcells(
X = da_cells, pred.thres = c(-cutoff,cutoff),
plot.embedding = umap
)

iso$DA_pred = da_cells$da.pred
iso$DA_sign = ifelse(da_cells$da.pred < 0, paste0("DA_", condition_1), paste0("DA_", condition_2))

da_regions <- getDAregion(
X = pca,
da.cells = da_cells,
cell.labels = iso$layer_source,
labels.1 = labelsCenter,
labels.2 = labelsSurface,
resolution = 0.2,
plot.embedding = umap,
min.cell = 50
)

table(da_regions$da.region.label)
iso$DA_regions=da_regions$da.region.label
iso$DA_regions_layer = ifelse(iso$DA_regions == 0, paste(iso$DA_regions, "UNC", sep=":"), paste(iso$DA_regions, iso$DA_sign, sep=":"))

isoreg = subset(iso, subset = DA_regions > 0)

iso2D = iso
Idents(iso2D) = factor(iso$DA_regions_layer)
iso2D_regions = isoreg
Idents(iso2D_regions) = factor(isoreg$DA_regions_layer)

save(iso2D, iso2D_regions, da_regions, da_cells, file="data/SEEP_Regions.rda")
