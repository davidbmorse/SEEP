rm(list=ls())
setwd("C://PROJECTS/P2022/SEEP_Manuscript")

library(Seurat)
library(tidyverse)
library(dplyr)
library(ggsci)
library(patchwork)

load("Input/GSE146026_10x.rda")

tpm1 = GSE146026_10x_logdata %>% column_to_rownames("Gene")
tpm1 = apply(tpm1,2,function(x) as.numeric(x))
rownames(tpm1) = GSE146026_10x_logdata$Gene
rownames(tpm1) = sub("^HLA\\.","HLA-", rownames(tpm1))
tpm = (2^tpm1-1)/10
tpm = log1p(tpm)

add_metadata = read.csv("Input/GSE146026_metadata_10x_clean.csv")
colnames(add_metadata) = sub("ï..", "", colnames(add_metadata))
GSE146026_10x_metadata$patient = as.numeric(GSE146026_10x_metadata$patient)
metadata = inner_join(GSE146026_10x_metadata, add_metadata, by=c("Cell_ID"="Sample_title","patient"="Patient_id")) %>% dplyr::rename("tSNE1_orig"="TSNE_x", "tSNE2_orig"="TSNE_y") %>%
  dplyr::mutate(Patient = ifelse(Time_point==1,paste(patient,Time_point, sep='.'), patient)) %>%
  dplyr::filter(! clst %in% c(-1,19:21) ) %>% dplyr::rename(Tissue=Cell_type)
  rownames(metadata) = metadata$Cell_ID
Cell_type = factor(metadata$clst)
levels(Cell_type) = c("Cancer", rep("Macrophage", 4), rep("Immune_cells", 5),rep("Cancer",4), rep("Fibroblasts",4))
Cell_types = factor(metadata$clst)
levels(Cell_types) = c("Cancer", rep("Macrophage", 4), "DC","DC","B","T", "Ery.", rep("Cancer",4), rep("Fibroblasts",4))
metadata$Cell_type = Cell_type
metadata$Cell_types = Cell_types
metadata$npg_clusters = as.numeric(metadata$clst)

# filter removed cells (unknown cluster assignment compared to the paper figure)
tpm = tpm[,match(metadata$Cell_ID, colnames(tpm))]

so <- CreateSeuratObject(counts = tpm, meta.data = metadata)
so <- FindVariableFeatures(object = so, nfeatures=5000)
s.genes = cc.genes$s.genes; #s.genes[s.genes=="MLF1IP"] = "CENPU"
g2m.genes = cc.genes$g2m.genes
so = CellCycleScoring(so,s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
so = ScaleData(so, features=rownames(so),  vars.to.regress=c("G2M.Score","S.Score"))
so <- RunPCA(object = so)
so <- FindNeighbors(object = so)
so <- FindClusters(object = so)
so <- RunTSNE(object = so, dim.embed=2, perplexity=30)
so <- RunUMAP(so, reduction = "pca", dims = 1:30)

orig = so@meta.data[,c("tSNE1_orig","tSNE2_orig")]
orig = t(apply(orig, 1, as.numeric))
colnames(orig) = c("tSNE_1","tSNE_2")
so[["tsne2d"]] <- CreateDimReducObject(embeddings = orig, key = "tSNE2d_", assay = DefaultAssay(so))


p1a=DimPlot(object = so, label.size=4, pt.size=0.5, reduction = "tsne2d", group.by="npg_clusters", label.color="black",label=TRUE,label.box=F, repel=F, cols=pal_ucscgb()(26)) + theme(aspect.ratio=1)
p1b=DimPlot(object = so, label.size=4, pt.size=0.5, reduction = "tsne2d", group.by="npg_clusters", split.by="Patient", label.color="black",label=F,label.box=F, repel=F, ncol=4,cols=pal_ucscgb()(26)) + theme(aspect.ratio=1)

p2a=DimPlot(object = so, label.size=4, pt.size=0.5, reduction = "tsne2d", group.by="Cell_type", label.color="black",label=F, ncol=1, repel=F, cols=pal_ucscgb()(26)[c(1,4,5,7)]) + theme(aspect.ratio=1)
p2aa=DimPlot(object = so, label.size=4, pt.size=0.5, reduction = "tsne2d", group.by="Cell_types", label.color="black",label=F, ncol=1, repel=F, cols=pal_ucscgb()(26)[-3]) + theme(aspect.ratio=1)
p2b=DimPlot(object = so, label.size=4, pt.size=.5, reduction = "tsne2d", group.by="Cell_type", split.by="Patient", label.color="black",label=F, ncol=4, repel=F, cols=pal_ucscgb()(26)[c(1,4,5,7)]) + theme(aspect.ratio=1)

png("10x_tSNE_npg_clusters.png", res=300, wid=10, hei=6, units='in')
print(p1a)# + guide_area() + plot_layout(widths=c(1), guides="collect"))
dev.off()

png("10x_tSNE_npg_cell_types.png", res=300, wid=10, hei=6, units='in')

print(p2aa) # + guide_area() + plot_layout(widths=c(1), guides="collect"))
dev.off()

saveRDS(so, file="data/SO_GSE146026_10x.rds")
