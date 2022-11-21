rm(list=ls())
setwd("C://PROJECTS/P2022/SEEP_Manuscript")

library(Seurat); library(tidyverse); library(SeuratObject)

load("data/SEEP_Merge.rda")

# normalize each source ====
sol = SplitObject(seep, split.by = 'source')
s.genes = cc.genes$s.genes; s.genes[s.genes=="MLF1IP"] = "CENPU"; g2m.genes = cc.genes$g2m.genes
sol = lapply(sol, CellCycleScoring, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
sol <- lapply(sol, NormalizeData)
sol = lapply(sol, function(x) ScaleData(x, vars.to.regress=c("percent.mt","G2M.Score","S.Score"), features=rownames(x)))
sph=GetAssayData(sol$sph[["RNA"]], slot="scale.data")
org=GetAssayData(sol$org[["RNA"]], slot="scale.data")
bio=GetAssayData(sol$bio[["RNA"]], slot="scale.data")
scaled_RNAassay = cbind(sph, org, bio)
sol <- lapply(sol, FindVariableFeatures, selection.method = "vst", nfeatures=2000)

# find anchors and integrate ====
seep.features <- SelectIntegrationFeatures(object.list = sol, nfeatures=2000)
seep.anchors <- FindIntegrationAnchors(object.list = sol, anchor.features = seep.features, verbose = TRUE)
iso <- IntegrateData(anchorset = seep.anchors, verbose = TRUE)
# scale integrated data
iso <- ScaleData(iso, vars.to.regress=c("percent.mt","G2M.Score","S.Score"))
iso@assays$RNA@scale.data = scaled_RNAassay

# save output ====
isoSEEP = iso
save(isoSEEP, file="data/SEEP_Integrate.rda")

