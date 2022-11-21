rm(list=ls())

library(Seurat);library(dplyr);library(ggplot2); library(patchwork)

data("SEEP_ORApaths")
data("SS2_ORApaths")
data("10X_ORApaths")
data("top_cells")
data("SEEP_Markers")

isoreg = subset(iso2_scores, cells = seep_top_cells)
ss2 = subset(SS2_scores,  cells = ss2_top_cells)
x10 = subset(X10_scores, cells = x10_top_cells)

levels(isoreg$Regions) = c("Surface\nCluster 1", "Surface\nCluster 2", "Surface\nCluster 3", "Center\nCluster 1")
levels(ss2$Regions) = c("Surface\nCluster 1", "Surface\nCluster 2", "Surface\nCluster 3", "Center\nCluster 1")
levels(x10$Regions) = c("Surface\nCluster 1", "Surface\nCluster 2", "Surface\nCluster 3", "Center\nCluster 1")
Idents(isoreg) = isoreg$Regions
Idents(ss2) = ss2$Regions
Idents(x10) = x10$Regions

Marker = "ZWINT"
sopat = x10
markers = seep_regionMarkers
x = markers %>% filter(Gene %in% Marker)
colo2=ggsci::pal_d3()(4)
colo2 = sub("FF$","55", colo2)
names(colo2) = levels(isoreg$Regions)
qua <- function(x){
  out <- quantile(x, probs = c(0.95))
  names(out) <- c("ymax")
  return(out)
}

seep_violin = VlnPlot(object=isoreg, pt.size=0.1, log=T, slot='data', sort=T, features=x$Gene, split.by="Regions") +
  ggtitle("HGSOC Model Cells") +
  xlab("") +
  geom_boxplot(alpha=0.5, outlier.shape = NA, width=0.3, col='black', show.legend = FALSE) +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.5, colour = "black", show.legend = FALSE) +
  ylab("ZWINT Log Expression")  +
  theme(plot.title = element_text(size = 12, face = "plain")) +
  scale_fill_manual(values=colo2) + #  stat_summary(fun=qua,geom='point', shape=21, stroke=1, show.legend = F, size=3) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(angle=60, size=12), axis.ticks = element_blank())

x10_violin = VlnPlot(object=x10, log=T,pt.size=0.1,  slot='data', sort=T, features=x$Gene, split.by="Regions") +
  ggtitle("HGSOC Ascites Cells") +
  xlab("") + #geom_violin(alpha=0.8) +
  geom_boxplot(alpha=0.5, outlier.shape = NA, width=0.3, col='black', show.legend = FALSE) +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.5, colour = "black", show.legend = FALSE) +
  ylab("ZWINT Log Expression")  +
  theme(plot.title = element_text(size = 12, face = 'plain')) +
  scale_fill_manual(values=colo2) +
  #stat_summary(fun=qua,geom='point', shape=21, stroke=1, show.legend = F, size=3) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(angle=60, size=12), axis.ticks = element_blank())

ss2_violin = VlnPlot(object=ss2, pt.size=0.1, log=T, slot='data', sort=T, features=x$Gene, split.by="Regions") +
  ggtitle("HGSOC Ascites Tumor Cells Isolated") +
  xlab("") +
  geom_boxplot(alpha=0.5, outlier.shape = NA, width=0.3, col='black', show.legend = FALSE) +
  stat_summary(fun = "mean", geom = "crossbar", width = 0.5, colour = "black", show.legend = FALSE) +   ylab("ZWINT Log Expression")  +
  theme(plot.title = element_text(size = 12, face = 'plain')) +
  scale_fill_manual(values=colo2)+
  #stat_summary(fun=qua,geom='point', shape=21, stroke=1, show.legend = F, size=3) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(angle=60, size=12), axis.ticks = element_blank())

tiff("Figures/SFigXE_11.20.22.tiff", width=12, height=5, compression = "lzw", res=300, units='in')
  print(
    (seep_violin | x10_violin | ss2_violin)
     # + plot_annotation(tag_levels = list(c("A","B","C"))) & theme(plot.tag = element_text(size = 16, face='plain'))
  )
dev.off()

pdf("Figures/SFigXE_11_20.22.pdf", width=12, height=5)
print(
    (seep_violin | x10_violin | ss2_violin)
)
dev.off()

