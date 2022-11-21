rm(list=ls())
setwd("C://PROJECTS/P2022/SEEP_Manuscript")

library(Seurat);library(ggplot2);library(GGally); library(tidyverse)

# data
load("data/SEEP_Regions.rda")
sopat=readRDS("data/SO_GSE146026_SS2.rds")

iso <- iso2D_regions
iso <- RunPCA(iso, ndims=30)
pct <- iso[["pca"]]@stdev / sum(iso[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine last point where change of % of variation is more than 0.1%.
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# Minimum of the two calculation
n <- min(co1, co2)
iso <- RunUMAP(iso, dims = 1:n, reduction = "pca", return.model = TRUE)

# map UMAP and transfer labels
anchors <- FindTransferAnchors(reference = iso, query = sopat, reduction = 'cca', npcs=n, dims=1:n, k.anchor=20)
predictions <- TransferData(anchorset = anchors, refdata = iso$DA_regions_layer, weight.reduction = 'cca')
sopat <- AddMetaData(sopat, metadata = predictions)
sopat <- MapQuery(anchorset = anchors, reference = iso, query = sopat, refdata = list(clusters = "DA_regions_layer"), reduction.model = "umap")

# scale good prediction data (RNA assay) for all genes
prediction = sopat@meta.data %>% select(contains("prediction.score"), -ends_with("max")) %>%
  select(contains("Surface"), contains("Center"))
good = apply(prediction, 1, function(x){
  index = order(x, decreasing=TRUE)
  max1 = x[index[1]]
  max2 = x[index[2]]
  dif = max1-max2
  return(c(max1, max2, dif))
})
good = data.frame(t(good))
colnames(good) = c("prediction.max1", "prediction.max2", "prediction.maxdiff")
prediction = data.frame(prediction, good)
sopat = AddMetaData(sopat, metadata=good)
tsoSS2 = sopat
tsoSS2_scores = prediction

sel.good = sopat$prediction.score.max > 0.5 & sopat$prediction.maxdiff>=0.25
goodSS2 = subset(sopat, cells=sopat$Cell_ID[sel.good])

isoSEEP_ref=iso

save(tsoSS2, tsoSS2_scores, goodSS2, isoSEEP_ref, file="data/SS2trans_Regions.rda")

