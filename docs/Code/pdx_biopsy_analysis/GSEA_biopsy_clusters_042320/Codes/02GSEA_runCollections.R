setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190429_GSEA_scRNA_organoid_layers")

library(Seurat)
library(Matrix)
library(dplyr)
library(stats)

#functions (from source)
source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190429_GSEA_scRNA_organoid_layers/Codes/01GSEA_pGSEAfunction.R")

gmtlist(path.gmtFolder="/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190429_GSEA_scRNA_organoid_layers/MSigDB/MSigDB_subsets",
              cut.version='.v6.1.symbols.gmt',
                rda.file="MSigDB/MSigDB_subsets.rda")


# load pre-ranked ave_logFC for layer vs. total sphere
data("PrerankOrg")
gsea_CENTRAL <- pGSEA(rank.metric=difs$c, names.metric=rownames(difs),msigdb="MSigDB/MSigDB_subsets.rda")
gsea_MIDDLE <- pGSEA(rank.metric=difs$m, names.metric=rownames(difs), msigdb="MSigDB/MSigDB_subsets.rda")
gsea_SURFACE <- pGSEA(rank.metric=difs$s, names.metric=rownames(difs),msigdb="MSigDB/MSigDB_subsets.rda")

gsea <- pGSEA(rank.metric=x, names.metric=names(x),msigdb="MSigDB/MSigDB_subsets.rda")
save(gsea_c, gsea_m, gsea_s,
     file = 'data/pGSEA.rda')

