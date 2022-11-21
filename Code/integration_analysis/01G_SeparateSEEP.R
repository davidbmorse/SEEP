rm(list=ls())
setwd("C://PROJECTS/P2022/SEEP_Manuscript")

library(Seurat); library(tidyverse); library(vcd); library(SeuratObject)


s.genes = cc.genes$s.genes; s.genes[s.genes=="MLF1IP"] = "CENPU"
g2m.genes = cc.genes$g2m.genes
colo3 = ggsci::pal_startrek()(3)
colo2 = ggsci::pal_npg()(2)

load("data/SEEP_Merge.rda")

# normalize across source
sol = list(seep=seep)
sol = lapply(sol, CellCycleScoring, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
sol <- lapply(sol, NormalizeData)
sol <- lapply(sol, FindVariableFeatures, selection.method = "vst", nfeatures=2000)

# scale data
sso <- ScaleData(sol$seep, vars.to.regress=c("percent.mt","G2M.Score","S.Score"), features=rownames(sol$seep))

# run PCA reduction and identify clusters
sso <- RunPCA(sso, npcs=50, verbose = FALSE)

# run tSNE
sso <- RunTSNE(object=sso, dims=1:20, dim.embed=2)
sso3D <- RunTSNE(object=sso, dims=1:20, dim.embed=3)

# run UMAP
sso2D <- RunUMAP(sso, reduction="pca", dims=1:20, return.model = FALSE)
sso3D <- RunUMAP(sso3D, reduction="pca", dims=1:20, return.model = FALSE, n.components=3)

save(sso2D, sso3D, file="data/SEEP_Separate.rda")
