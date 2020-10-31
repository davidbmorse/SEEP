setwd("C://PROJECTS/DM_scRNa")

library(Seurat)
library(Matrix)
library(dplyr)

gmtlist(path.gmtFolder="C:/PROJECTS/DM_scRNA/MSigDB/MSigDB_subsets",
              cut.version='.v6.1.symbols.gmt',
                rda.file="MSigDB/MSigDB_subsets.rda")


# load pre-ranked ave_logFC for layer vs. total sphere
data("Preranked")
gsea_CENTRAL <- pGSEA(rank.metric=difs$CENTRAL, names.metric=rownames(difs),msigdb="MSigDB/MSigDB_subsets.rda")
gsea_INNER <- pGSEA(rank.metric=difs$INNER, names.metric=rownames(difs), msigdb="MSigDB/MSigDB_subsets.rda")
gsea_OUTER <- pGSEA(rank.metric=difs$OUTER, names.metric=rownames(difs),msigdb="MSigDB/MSigDB_subsets.rda")
gsea_SURFACE <- pGSEA(rank.metric=difs$SURFACE, names.metric=rownames(difs),msigdb="MSigDB/MSigDB_subsets.rda")

save(gsea_CENTRAL, gsea_INNER, gsea_OUTER, gsea_SURFACE,
     file = 'data/pGSEA.rda')



