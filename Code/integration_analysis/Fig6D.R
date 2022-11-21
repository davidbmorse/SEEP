rm(list=ls())
setwd("C://PROJECTS/P2022/SEEP_Manuscript")

library(tidyverse); library(dplyr); library(Seurat); library(UCell); library(patchwork)

# data
data("SEEP_ORApaths")
data("SEEP_Markers")
data("SEEP_ORAgenenet")

# markers ====
iso = iso2_scores
print(table(iso$DA_regions_layer))
liso = SplitObject(iso, "DA_regions_layer")
nams = names(liso)
liso = sapply(nams, function(y, n.min=200, n.max=2000, perc=95){
  x=liso[[y]]
  if(ncol(x)<=n.min){
    nc = ncol(x)
  } else {
    nc=c(perc*ncol(x)/100,n.max)
    nc = nc[which.min(nc)]
  }
  return(TopCells(x, dim=1, red='umap', ncells=nc))
})
top_cells = unlist(liso)
print(sapply(liso,length))
iso=subset(iso, cells=top_cells)
print(table(iso$DA_regions_layer))

markers = seep_regionMarkers
levels(Idents(iso))=c("Surface\nCluster 1", "Surface\nCluster 2", "Surface\nCluster 3", "Center\nCluster 1")
markers = group_by(markers, Comparison_cluster) %>%
  mutate(rpval = rank(max_pval),
         rfc = rank(-min_avg_log2FC)) %>%
  ungroup()

nodes = unique(unlist(sapply(lnodes, function(x) x$name)))
mark = markers %>% filter(Gene %in% nodes)
markGeneHeatmap = mark

# suppl table ====
stab <- mark[,1:27] %>% select(-Comparison_id) %>%
  mutate(Comparison_cluster = factor(Comparison_cluster))
stab_marker <- stab
  levels(stab$Comparison_cluster) = c("Suface Cluster 1", "Surface Cluster 2", "Surface Cluster 3", "Center Cluster 1")
  colnames(stab)[1] = "SEEP_DA_Cluster"
write.table(stab, "Figures/Supplement/Conserved_Markers.tsv", sep="\t", row=F, quote=F)


# fig6D ====
ran = c(-3,3)
colo2=ggsci::pal_d3()(4)
names(colo2) = levels(Idents(iso))

tiff("Figures/Fig6D_11.20.22.tiff", compress="lzw", res=300, units='in', width=14, height=3)
heatmap=DoHeatmap(iso, label=F, cells=top_cells, features=mark$Gene, group.bar.height = 0.04, lines.width=1, raster=FALSE, disp.min=ran[1], disp.max=ran[2], slot='scale.data',group.colors=colo2)+theme(axis.text.y  = element_blank()) + scale_fill_gradient2(high='yellow',low='magenta', mid='black', midpoint=0, limits=ran, , na.value = "white", name='sd')
#high='#FFFF00',low='magenta'
heatmap = heatmap + scale_color_manual(values=colo2) +
  theme(legend.key.height = unit(0.6, 'cm'),
        legend.key.width = unit(0.45, 'cm'),
        legend.text=element_text(size=7)) + guides(colour = guide_legend(override.aes = list(size=3, shape=15, alpha = 1), title=NULL)) + geom_hline(yintercept = cumsum(rev(table(mark$Comparison_cluster)))+0.5, color='white', size=0.3)
print(heatmap)
dev.off()

pdf("Figures/Fig6D_11.20.22.pdf", width=14, height=3)
heatmap=DoHeatmap(iso, label=F, cells=top_cells, features=mark$Gene, group.bar.height = 0.04, lines.width=1, raster=FALSE, disp.min=ran[1], disp.max=ran[2], slot='scale.data',group.colors=colo2)+theme(axis.text.y  = element_blank()) + scale_fill_gradient2(high='yellow',low='magenta', mid='black', midpoint=0, limits=ran, , na.value = "white", name='sd')
#high='#FFFF00',low='magenta'
heatmap = heatmap + scale_color_manual(values=colo2) +
  theme(legend.key.height = unit(0.6, 'cm'),
        legend.key.width = unit(0.45, 'cm'),
        legend.text=element_text(size=7)) + guides(colour = guide_legend(override.aes = list(size=3, shape=15, alpha = 1), title=NULL)) + geom_hline(yintercept = cumsum(rev(table(mark$Comparison_cluster)))+0.5, color='white', size=0.3)
print(heatmap)
dev.off()

tiff("Figures/Additional//Fig6DnoLines_11.20.22.tiff", compress="lzw", res=300, units='in', width=14, height=3)
heatmap=DoHeatmap(iso, label=F, cells=top_cells, features=mark$Gene, group.bar.height = 0.04, lines.width=1, raster=FALSE, disp.min=ran[1], disp.max=ran[2], slot='scale.data',group.colors=colo2)+theme(axis.text.y  = element_blank()) + scale_fill_gradient2(high='yellow',low='magenta', mid='black', midpoint=0, limits=ran, , na.value = "white", name="sd")
#high='#FFFF00',low='magenta'
heatmap = heatmap + scale_color_manual(values=colo2) +
  theme(legend.key.height = unit(0.6, 'cm'),
        legend.key.width = unit(0.45, 'cm'),
        legend.text=element_text(size=7)) + guides(colour = guide_legend(override.aes = list(size=3, shape=15, alpha = 1), title=NULL)) #+ geom_hline(yintercept = cumsum(rev(table(mark$Comparison_cluster)))+0.5, color='white', size=0.5)
print(heatmap)
dev.off()


# fc ====
nodes = unique(unlist(sapply(lnodes, function(x) x$name)))
mark = markers %>% filter(Gene %in% nodes) %>% filter(rfc<=25)

tiff("Figures/Additional/SEEP_HeatmapPathwayGenes_FC25.tiff", compress="lzw", res=300, units='in', width=16, height=5)
heatmap=DoHeatmap(iso, label=F, cells=top_cells, features=mark$Gene, group.bar.height = 0.04, lines.width=1, raster=FALSE, disp.min=ran[1], disp.max=ran[2], slot='scale.data',group.colors=colo2)+theme(axis.text.y  = element_text(size=5)) + scale_fill_gradient2(high='yellow',low='magenta', mid='black', midpoint=0, limits=ran, , na.value = "white", name='sd')
#high='#FFFF00',low='magenta'
heatmap + theme(legend.key.height = unit(0.5, 'cm'),
                legend.key.width = unit(0.4, 'cm'),
                legend.text=element_text(size=8)) +
  guides(colour = guide_legend(override.aes = list(size=3, shape=15, alpha = 1), title=NULL))+ geom_hline(yintercept = cumsum(rev(table(mark$Comparison_cluster)))+0.5, color='white', size=0.5) + theme(axis.text.y = element_text(colour = rep(rev(colo2), rev(table(mark$Comparison_cluster))))) + scale_color_manual(values = colo2)
dev.off()

# pval ====
nodes = unique(unlist(sapply(lnodes, function(x) x$name)))
mark = markers %>% filter(Gene %in% nodes) %>% filter(rpval<=25)

tiff("Figures/Additional/SEEP_HeatmapPathwayGenes_Pval25.tiff", compress="lzw", res=300, units='in', width=16, height=5)
heatmap=DoHeatmap(iso, label=F, cells=top_cells, features=mark$Gene, group.bar.height = 0.04, draw.lines=TRUE, lines.width=1, raster=FALSE, disp.min=ran[1], disp.max=ran[2], slot='scale.data',group.colors=colo2)+theme(axis.text.y  = element_text(size=5)) + scale_fill_gradient2(high='yellow',low='magenta', mid='black', midpoint=0, limits=ran, na.value = "white", name='sd')
#high='#FFFF00',low='magenta'
heatmap=heatmap + theme(legend.key.height = unit(0.5, 'cm'),
                legend.key.width = unit(0.4, 'cm'),
                legend.text=element_text(size=8)) +
  guides(colour = guide_legend(override.aes = list(size=3, shape=15, alpha = 1), title=NULL)) + geom_hline(yintercept = cumsum(rev(table(mark$Comparison_cluster)))+0.5, color='white', size=0.5) + theme(axis.text.y = element_text(colour = rep(rev(colo2), rev(table(mark$Comparison_cluster))))) + scale_color_manual(values = colo2)
print(heatmap)
dev.off()

isoHeatmap = iso
save(markGeneHeatmap, top_cells, isoHeatmap , stab_marker, file='data/Fig6D_heatmap.rda')