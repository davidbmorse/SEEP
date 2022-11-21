rm(list=ls())
setwd("C://PROJECTS/P2022/SEEP_Manuscript")

# data
load("data/SEEP_Regions.rda")

library(Seurat); library(tidyverse)

isoreg = iso2D_regions

# Set default Seurat assay to RNA (individual sources, all genes)
DefaultAssay(isoreg) <- "RNA"
Idents(isoreg)=factor(isoreg$DA_regions_layer)

# Find cluster vs. all markers
clusters = levels(Idents(isoreg))
n=length(clusters)
list_markers = vector('list',n)
names(list_markers)=clusters

for(i in 1:n){
list_markers[[i]] <- FindConservedMarkers(isoreg, ident.1 = clusters[i], grouping.var = "source", verbose = FALSE, only.pos=TRUE, logfc=0, min.pct=0.1)
}

# Add variables
list_markers_all = lapply(list_markers, function(x) {

  colnames(x)[colnames(x) == "minimump_p_val" ] = "min_pval"
  x$min_padj = apply(x[,grep("adj$", colnames(x))], 1, min, na.rm=T)
  x$max_padj = apply(x[,grep("adj$", colnames(x))], 1, max, na.rm=T)
  x$comb_pval = apply(x[,grep("_p_val$", colnames(x))], 1, function(f) metap::sumlog(f)$p)
  x$comb_padj = apply(x[,grep("_val_adj$$", colnames(x))], 1, function(f) metap::sumlog(f)$p)
  x$min_avg_log2FC = apply(x[,grep("avg_log2FC$", colnames(x))], 1, min, na.rm=T)
  x$max_avg_log2FC = apply(x[,grep("avg_log2FC$", colnames(x))], 1, max, na.rm=T)
  x$min_pct.1 = apply(x[,grep("pct.1$", colnames(x))], 1, min, na.rm=T)
  x
}
)
list_markers_significant = lapply(list_markers_all, function(x){
  x = x[x$max_pval <= 0.05,]
  x$rank_comb_pval = rank(x[,"comb_pval"], ties.method = 'min')
  x$rank_comb_padj = rank(x[,"comb_padj"], ties.method = 'min')
  x$rank_max_pval = rank(x[,"max_pval"], ties.method = 'min')
  x$rank_max_padj = rank(x[,"max_padj"], ties.method = 'min')
  x$rank_min_avg_log2FC = rank(x[,"min_avg_log2FC"], ties.method = 'min')
  x$rank_max_avg_log2FC = rank(x[,"max_avg_log2FC"], ties.method = 'min')
  x$rank_min_pct.1 = rank(x[,"min_pct.1"], ties.method = 'min')
  x
})
nsig = sapply(list_markers_significant,nrow)

allcol = colnames(list_markers_significant$`1:DA_Surface`)
miscol = setdiff(allcol,colnames(list_markers_significant$`3:DA_Surface`))
misdat = data.frame(matrix(NA,nrow(list_markers_significant$`3:DA_Surface`), length(miscol)))
colnames(misdat) = miscol
dat = data.frame(list_markers_significant$`3:DA_Surface`, misdat)
list_markers_significant$`3:DA_Surface` = dat[,match(allcol,colnames(dat))]

# final marker data set
markers = do.call(rbind, list_markers_significant)

markers <- markers %>%  tibble::rownames_to_column("Comparison_id") %>%
  dplyr::mutate(
    Comparison_cluster=sapply(strsplit(Comparison_id, "\\."), function(x) x[1]),
    Gene = sapply(strsplit(Comparison_id, "\\."), function(x) x[2])) %>%

  dplyr::select(Comparison_cluster,
                starts_with("adj_pval"),
                Gene,
                ends_with("log2FC"),
                ends_with("p_val"), ends_with("p_val_adj"),
                starts_with(c("min","max")),
                everything())

duplicates = markers$Gene[duplicated(markers$Gene)]
markers$unique_marker = ifelse(markers$Gene %in% duplicates, FALSE, TRUE)
ltmp = split(markers, markers$Gene)
ltmp = lapply(ltmp, function(x) x[x$max_pval == min(x$max_pval),])
markers = data.frame(do.call(rbind,ltmp)) %>% arrange(Comparison_cluster, max_pval)

write.table(markers, file="SEEPregions_markers.txt", sep="\t", row=FALSE, quote=FALSE)

seep_regionMarkers = markers

# Save marker lists
save(seep_regionMarkers, file='data/SEEP_Markers.rda')

