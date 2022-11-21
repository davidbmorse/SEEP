rm(list=ls())
# SEEP ====

library(Seurat);library(dplyr);library(ggplot2); library(patchwork)

data("Fig6D_heatmap")
data("SEEP_ORApaths")
data("SEEP_ORAgenenet")
data("SEEP_Markers")
data("SEEP_ORA")

perc=95
n.min=200
n.max=2000
ran = c(-3,3)

# select ====
iso = iso2_paths
print(table(iso$DA_regions_layer))
liso = SplitObject(iso, "DA_regions_layer")
nams = names(liso)
liso = sapply(nams, function(y){
  x=liso[[y]]
  if(ncol(x)<=n.min){
    nc = ncol(x)
  } else {
    nc=c(perc*ncol(x)/100,n.max)
    nc = nc[which.min(nc)]
  }
  return(TopCells(x, dim=1, red='umap', ncells=nc))
})
top1 = unlist(liso)
print(sapply(liso,length))
iso=subset(iso, cells=top1)
markers = seep_pathwayScores
levels(Idents(iso))=c("Surface\nCluster 1", "Surface\nCluster 2", "Surface\nCluster 3", "Center\nCluster 1")
print(table(iso$DA_regions_layer))

colo2=ggsci::pal_d3()(4)
names(colo2) = levels(Idents(iso))

lora = sapply(names(lorasig_main), function(x) {
  y = lorasig_main[[x]] %>% slice_min(pval, n=3)
  y = paste0(y$pathway,"-UCell.", x)
  gsub("_","-",y)
})
oras = unlist(lora)

paths = rownames(iso)
paths_region = sapply(strsplit(paths, "\\."), function(x) x[2])
paths_names = gsub("HALLMARK-|-UCell", "", sapply(strsplit(paths, "\\."), function(x) x[1]))
tab=table(paths_names)
unipath=names(tab)[tab>=1]
selpath = paths_names %in% unipath


tiff("Figures/Additional/SEEP_HeatmapORA.tiff", compress="lzw", res=300, units='in', width=16, height=8)
heatmap=DoHeatmap(iso, raster=FALSE,label=F, disp.min=ran[1], disp.max=ran[2], group.bar.height = 0.07, lines.width=2,features=paths[selpath],slot='scale.data',group.colors=colo2)+theme(axis.text.y  = element_text(size=3)) + scale_fill_gradient2(high='yellow',low='darkmagenta', mid='black', midpoint=0, limits=ran, na.value = 'white', name='sd')

heatmap = heatmap + theme(axis.text.y = element_text(size=8, colour = rep(rev(colo2), rev(table(paths_region[selpath]))))) + scale_color_manual(values = colo2) +
  theme(legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10)) +
  guides(colour = guide_legend(override.aes = list(size=3, shape=15, alpha = 1), title=NULL)) +
  scale_y_discrete(labels = rev(paths_names[selpath])) + geom_hline(yintercept = cumsum(rev(table(paths_region[selpath])))+0.5, color='white', size=0.5) + ggtitle("SEEP HGSOC Cells") +  theme(plot.title = element_text(size = 12, face = "plain"))
print(heatmap)
seep = heatmap
dev.off()


# SS2 ====

data("SS2_ORApaths")

# select ====
iso = SS2_paths
print(table(iso$predicted.clusters))
liso = SplitObject(iso, "predicted.clusters")
nams = names(liso)
liso = sapply(nams, function(y){
  x=liso[[y]]
  if(ncol(x)<=n.min){
    nc = ncol(x)
  } else {
    nc=c(perc*ncol(x)/100,n.max)
    nc = nc[which.min(nc)]
  }
  return(TopCells(x, dim=1, red='ref.umap', ncells=nc))
})
top1 = unlist(liso)
print(sapply(liso,length))
iso=subset(iso, cells=intersect(top1, goodSS2_cells))
markers = SS2_pathwayScores
levels(Idents(iso))=c("Surface\nCluster 1", "Surface\nCluster 2", "Surface\nCluster 3", "Center\nCluster 1")
print(table(iso$predicted.clusters))

colo2=ggsci::pal_d3()(4)
names(colo2) = levels(Idents(iso))

lora = sapply(names(lorasig_main), function(x) {
  y = lorasig_main[[x]] %>% slice_min(pval, n=3)
  y = paste0(y$pathway,"-UCell.", x)
  gsub("_","-",y)
})
oras = unlist(lora)

paths = rownames(iso)
paths_region = sapply(strsplit(paths, "\\."), function(x) x[2])
paths_names = gsub("HALLMARK-|-UCell", "", sapply(strsplit(paths, "\\."), function(x) x[1]))
tab=table(paths_names)
unipath=names(tab)[tab>=1]
selpath = paths_names %in% unipath


tiff("Figures/Additional/SS2_HeatmapORA.tiff", compress="lzw", res=300, units='in', width=16, height=8)
heatmap=DoHeatmap(iso, raster=FALSE,label=F, disp.min=ran[1], disp.max=ran[2], group.bar.height = 0.07, lines.width=1,features=paths[selpath],slot='scale.data',group.colors=colo2)+theme(axis.text.y  = element_text(size=3)) + scale_fill_gradient2(high='yellow',low='darkmagenta', mid='black', midpoint=0, limits=ran, na.value = 'white', name='sd')

heatmap = heatmap + theme(axis.text.y = element_text(size=8, colour = rep(rev(colo2), rev(table(paths_region[selpath]))))) + scale_color_manual(values = colo2) +
  theme(legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10)) +
  guides(colour = guide_legend(override.aes = list(size=3, shape=15, alpha = 1), title=NULL)) +
  scale_y_discrete(labels = rev(paths_names[selpath])) + geom_hline(yintercept = cumsum(rev(table(paths_region[selpath])))+0.5, color='white', size=0.5) + ggtitle("HGSOC Ascites Tumor Cells Isolated") +  theme(plot.title = element_text(size = 12, face = "plain"))
print(heatmap)
ss2 = heatmap
dev.off()

# 10X

data("10X_ORApaths")

# 10X ====

# select ====
iso = X10_paths
print(table(iso$predicted.clusters))
liso = SplitObject(iso, "predicted.clusters")
nams = names(liso)
liso = sapply(nams, function(y){
  x=liso[[y]]
  if(ncol(x)<=n.min){
    nc = ncol(x)
  } else {
    nc=c(perc*ncol(x)/100,n.max)
    nc = nc[which.min(nc)]
  }
  return(TopCells(x, dim=1, red='ref.umap', ncells=nc))
})
top1 = unlist(liso)
print(sapply(liso,length))
iso=subset(iso, cells=intersect(top1, good10X_cells))
markers = X10_pathwayScores
levels(Idents(iso))=c("Surface\nCluster 1", "Surface\nCluster 2", "Surface\nCluster 3", "Center\nCluster 1")
print(table(iso$predicted.clusters))

colo2=ggsci::pal_d3()(4)
names(colo2) = levels(Idents(iso))

lora = sapply(names(lorasig_main), function(x) {
  y = lorasig_main[[x]] %>% slice_min(pval, n=3)
  y = paste0(y$pathway,"-UCell.", x)
  gsub("_","-",y)
})
oras = unlist(lora)

paths = rownames(iso)
paths_region = sapply(strsplit(paths, "\\."), function(x) x[2])
paths_names = gsub("HALLMARK-|-UCell", "", sapply(strsplit(paths, "\\."), function(x) x[1]))
tab=table(paths_names)
unipath=names(tab)[tab>=1]
selpath = paths_names %in% unipath
keep_selpath = selpath

tiff(sprintf("Figures/Additional/10X_HeatmapORA.tiff"), compress="lzw", res=300, units='in', width=16, height=8)
heatmap=DoHeatmap(iso, raster=FALSE,label=F, disp.min=ran[1], disp.max=ran[2], group.bar.height = 0.07, lines.width=3,features=paths[selpath],slot='scale.data',group.colors=colo2)+theme(axis.text.y  = element_text(size=3)) + scale_fill_gradient2(high='yellow',low='darkmagenta', mid='black', midpoint=0, limits=ran, na.value = 'white', name='sd')

heatmap = heatmap + theme(axis.text.y = element_text(size=8, colour = rep(rev(colo2), rev(table(paths_region[selpath]))))) + scale_color_manual(values = colo2) +
  theme(legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(0.5, 'cm'),
        legend.text=element_text(size=10),
        legend.title=element_text(size=10)) +
        guides(colour = guide_legend(override.aes = list(size=3, shape=15, alpha = 1), title=NULL)) +
        scale_y_discrete(labels = rev(paths_names[selpath])) + geom_hline(yintercept = cumsum(rev(table(paths_region[selpath])))+0.5, color='white', size=0.5) +
  ggtitle("HGSOC Ascites Cells") + theme(plot.title = element_text(size = 12, face = "plain"))
print(heatmap)
x10 = heatmap
dev.off()

# Fig7 ====

#tiff(sprintf("Figures/Additional/HeatmapORA.tiff"), compress="lzw", res=300, units='in', width=14, height=20)
Fig7DEF = (seep/x10/ss2) +
  plot_annotation(tag_levels = list(c("D", "E", "F"))) + plot_layout(guides = "collect") &
  theme(plot.tag = element_text(size = 20))
#print(Fig7DEF)
#dev.off()

# cell type ====
iso = ScaleData(iso, scale.max = Inf, split.by="Cell_type")
lcells = split(iso$Cell_ID, iso$Cell_type)

for(i in names(lcells)){
  tab=table(subset(iso, cells=lcells[[i]])$predicted.clusters)
  if(! "3:DA_Surface" %in% names(tab)){
    selpath = selpath & paths_region != "3:DA-Surface"
    colo = colo2[-3]
  } else {
    colo = colo2
    selpath = keep_selpath
  }
  tiff(sprintf("Figures/Additional/10X_HeatmapORA_%s.tiff", i), compress="lzw", res=300, units='in', width=16, height=8)
  heatmap=DoHeatmap(iso, cells=lcells[[i]], raster=FALSE,label=F, disp.min=ran[1], disp.max=ran[2], group.bar.height = 0.07, lines.width=2,features=paths[selpath],slot='scale.data',group.colors=colo)+theme(axis.text.y  = element_text(size=3)) + scale_fill_gradient2(high='yellow',low='darkmagenta', mid='black', midpoint=0, limits=ran, na.value = 'white', name='sd')

  print(heatmap + theme(axis.text.y = element_text(size=8, colour = rep(rev(colo), rev(table(paths_region[selpath]))))) + scale_color_manual(values = colo) +
          theme(legend.key.height = unit(1, 'cm'),
                legend.key.width = unit(0.5, 'cm'),
                legend.text=element_text(size=8),
                legend.title=element_text(size=8)) +
          guides(colour = guide_legend(override.aes = list(size=3, shape=15, alpha = 1), title=NULL)) +
          scale_y_discrete(labels = rev(paths_names[selpath])) + geom_hline(yintercept = cumsum(rev(table(paths_region[selpath])))+0.5, color='white', size=0.5))
  dev.off()
}
