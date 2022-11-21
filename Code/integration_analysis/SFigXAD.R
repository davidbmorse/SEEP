rm(list=ls())

library(Seurat);library(dplyr);library(ggplot2); library(patchwork)

data("SEEP_ORApaths")
data("SS2_ORApaths")
data("10X_ORApaths")
data("top_pathways")
#data("SEEP_Markers")

isoreg = subset(iso2_paths, cells=TopCells(iso2_paths, ncells = ncol(iso2_paths)*99/100, dim=1, red='umap'))
ss2 = subset(SS2_paths,  cells=TopCells(SS2_paths, ncells = ncol(SS2_paths)*99/100, dim=1, red='ref.umap'))
x10 = subset(X10_paths, cells=TopCells(X10_paths, ncells = ncol(X10_paths)*99/100, dim=1, red='ref.umap'))

levels(isoreg$Regions) = c("Surface\nCluster 1", "Surface\nCluster 2", "Surface\nCluster 3", "Center\nCluster 1")
levels(ss2$Regions) = c("Surface\nCluster 1", "Surface\nCluster 2", "Surface\nCluster 3", "Center\nCluster 1")
levels(x10$Regions) = c("Surface\nCluster 1", "Surface\nCluster 2", "Surface\nCluster 3", "Center\nCluster 1")
Idents(isoreg) = isoreg$Regions
Idents(ss2) = ss2$Regions
Idents(x10) = x10$Regions


#FeaturePlot(x10, feature=rownames(isoreg)[2],pt.size = 0.3, order=T, min.cutoff = 'q5',max.cutoff = 'q95') + scale_colour_gradientn(colours = kovesi.linear_kry_5_98_c75(50))

oras = unlist(lora)
oras = c("HALLMARK-TNFA-SIGNALING-VIA-NFKB-UCell.1:DA-Surface", "HALLMARK-INTERFERON-GAMMA-RESPONSE-UCell.1:DA-Surface", "HALLMARK-E2F-TARGETS-UCell.2:DA-Surface",  "HALLMARK-G2M-CHECKPOINT-UCell.2:DA-Surface")
names(oras) = NULL
bcolo = '#00000055'
colopal = "magma"
nams = sub("HALLMARK-", "", sapply(strsplit(oras, "-UCell"), function(x) x[1]))
fnams = gsub("-", "_", oras)
fnams = gsub(":|\\.", "_", fnams)


for(i in 1:length(oras)){

  tiff(sprintf("Figures/SFigXAD_%s_11.20.22.tiff", nams[i]), width=5, height=8, compression = 'lzw', units='in', res=300)

seep <- FeaturePlot(isoreg, feature=oras[i],pt.size = 0.1, order=T, min.cutoff = 'q5', max.cutoff='q95') +
  scale_colour_viridis_c(option=colopal, begin=0, end=1, direction=1) + theme_classic() +
  theme(panel.background = element_rect(fill = bcolo)) +
  ggtitle("HGSOC Model Cells")

p10 <- FeaturePlot(x10, feature=oras[i],pt.size = 0.1, order=T, min.cutoff = 'q5', max.cutoff='q95') +
  scale_colour_viridis_c(option=colopal, begin=0, end=1, direction=1) + theme_classic() +
  theme(panel.background = element_rect(fill = bcolo)) +
  ggtitle("HGSOC Ascites Cells")

p2 <- FeaturePlot(ss2, feature=oras[i],pt.size = 0.1, order=T, min.cutoff = 'q5', max.cutoff='q95') +
  scale_colour_viridis_c(option=colopal, begin=0, end=1, direction=1) + theme_classic() +
  theme(panel.background = element_rect(fill = bcolo)) +
  ggtitle("HGSOC Ascites Tumor Cells Isolated")

print( (seep + p10 + p2) + plot_annotation(title = oras[i], theme = theme(plot.title = element_text(size = 10)))  &
  labs(colour = 'Score', x = "UMAP1", y = "UMAP2") &
  ylim(-8,6) &
  theme(legend.key.height = unit(0.5, 'cm'),
        legend.key.width = unit(0.3, 'cm'),
        legend.text=element_text(size=8),
        legend.title=element_text(size=8))
)

dev.off()

}

pdf(sprintf("Figures/SFigX_11.20.22.pdf"), width=5, height=8)

for(i in 1:length(oras)){

  seep <- FeaturePlot(isoreg, feature=oras[i],pt.size = 0.001, order=T, min.cutoff = 'q5', max.cutoff='q95') +
    scale_colour_viridis_c(option=colopal, begin=0, end=1, direction=1) + theme_classic() +
    theme(panel.background = element_rect(fill = bcolo)) +
    ggtitle("HGSOC Model Cells")

  p10 <- FeaturePlot(x10, feature=oras[i],pt.size = 0.001, order=T, min.cutoff = 'q5', max.cutoff='q95') +
    scale_colour_viridis_c(option=colopal, begin=0, end=1, direction=1) + theme_classic() +
    theme(panel.background = element_rect(fill = bcolo)) +
    ggtitle("HGSOC Ascites Cells")

  p2 <- FeaturePlot(ss2, feature=oras[i],pt.size = 0.001, order=T, min.cutoff = 'q5', max.cutoff='q95') +
    scale_colour_viridis_c(option=colopal, begin=0, end=1, direction=1) + theme_classic() +
    theme(panel.background = element_rect(fill = bcolo)) +
    ggtitle("HGSOC Ascites Tumor Cells Isolated")

  print( (seep + p10 + p2) + plot_annotation(title = oras[i], theme = theme(plot.title = element_text(size = 10))) &
           labs(colour = 'Score', x = "UMAP1", y = "UMAP2") &
           ylim(-8,6) &
           theme(legend.key.height = unit(0.5, 'cm'),
                 legend.key.width = unit(0.3, 'cm'),
                 legend.text=element_text(size=8),
                 legend.title=element_text(size=8))
  )
}
dev.off()

oras = unlist(lora)
names(oras) = NULL
bcolo = '#00000055'
colopal = "magma"
nams = sub("HALLMARK-", "", sapply(strsplit(oras, "-UCell"), function(x) x[1]))
fnams = gsub("-", "_", oras)
fnams = gsub(":|\\.", "_", fnams)

for(i in 1:length(oras)){

  tiff(sprintf("Figures/Additional/%s.tiff", fnams[i]), width=5, height=8, compression = 'lzw', units='in', res=300)

  seep <- FeaturePlot(isoreg, feature=oras[i],pt.size = 0.1, order=T, min.cutoff = 'q5', max.cutoff='q95') +
    scale_colour_viridis_c(option=colopal, begin=0, end=1, direction=1) + theme_classic() +
    theme(panel.background = element_rect(fill = bcolo)) +
    ggtitle("HGSOC Model Cells")

  p10 <- FeaturePlot(x10, feature=oras[i],pt.size = 0.1, order=T, min.cutoff = 'q5', max.cutoff='q95') +
    scale_colour_viridis_c(option=colopal, begin=0, end=1, direction=1) + theme_classic() +
    theme(panel.background = element_rect(fill = bcolo)) +
    ggtitle("HGSOC Ascites Cells")

  p2 <- FeaturePlot(ss2, feature=oras[i],pt.size = 0.1, order=T, min.cutoff = 'q5', max.cutoff='q95') +
    scale_colour_viridis_c(option=colopal, begin=0, end=1, direction=1) + theme_classic() +
    theme(panel.background = element_rect(fill = bcolo)) +
    ggtitle("HGSOC Ascites Tumor Cells Isolated")

  print( (seep + p10 + p2) + plot_annotation(title = oras[i], theme = theme(plot.title = element_text(size = 10)))  &
           labs(colour = 'Score', x = "UMAP1", y = "UMAP2") &
           ylim(-8,6) &
           theme(legend.key.height = unit(0.5, 'cm'),
                 legend.key.width = unit(0.3, 'cm'),
                 legend.text=element_text(size=8),
                 legend.title=element_text(size=8))
  )

  dev.off()

}
