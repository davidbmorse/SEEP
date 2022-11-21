rm(list=ls())

setwd("C://PROJECTS/P2022/SEEP_Manuscript")

library(Seurat);library(ggplot2); library(tidyverse); library(patchwork)

data("10Xtrans_Regions")
data("SS2trans_Regions")

colo = ggsci::pal_d3()(4)
perc = 99

iso = isoSEEP_ref
nc=ncol(iso)*perc/100
top = TopCells(iso, dim=1, red='umap', n=nc)
iso = subset(iso, cells=top)
iso$color_by = factor(iso$DA_regions_layer)
levels(iso$color_by) = c("Surface\nCluster 1", "Surface\nCluster 2", "Surface\nCluster 3", "Center\nCluster 1")
dat <- data.frame(iso@reductions$umap@cell.embeddings) %>%
  rownames_to_column('Cells') %>%
  left_join(iso@meta.data %>% rownames_to_column("Cells"))

seep <-
  ggplot(dat, aes(x = UMAP_1, y = UMAP_2, col=color_by)) +
  geom_point(size=0.2) + scale_color_manual(name="", values=colo) +
  theme_classic() +
  labs(x="UMAP1", y="UMAP2", title="HGSOC Model Cells") +
  theme(plot.title = element_text(size = 9, face = "plain")) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme(legend.text=element_text(size = 7), legend.key.size = unit(0.7, 'cm')) +
  coord_cartesian(xlim = c(-12, 12), ylim=c(-7,6)) +
  theme(axis.title = element_text(size = 9))
dat_seep = dat

X10 <- tso10X
nc=ncol(X10)*perc/100
top = TopCells(X10, dim=1, red='ref.umap', n=nc)
X10 = subset(good10X, cells=top)
X10$color_by <- factor(X10$predicted.clusters)
levels(X10$color_by) = c("Surface\nCluster 1", "Surface\nCluster 2", "Surface\nCluster 3", "Center\nCluster 1")
dat <- data.frame(X10@reductions$ref.umap@cell.embeddings) %>%
  rownames_to_column('Cells') %>%
  left_join(X10@meta.data %>% rownames_to_column("Cells"))
dat_x10= dat

x10 <-
  ggplot(dat %>% arrange(desc(color_by)), aes(x = refUMAP_1, y = refUMAP_2, col=color_by)) +
  geom_point(size=0.2) + scale_color_manual(name="", values=colo) +
  theme_classic() +
  labs(x="ref.UMAP1", y="ref.UMAP2", title="HGSOC Ascites Cells")  +
  theme(plot.title = element_text(size = 9, face = "plain")) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme(legend.text=element_text(size = 7), legend.key.size = unit(0.7, 'cm')) +
  coord_cartesian(xlim = c(-12, 12), ylim=c(-7,6)) +
  theme(axis.title = element_text(size = 9))

SS2 <- tsoSS2
nc=ncol(SS2)*perc/100
top = TopCells(SS2, dim=1, red='ref.umap', n=nc)
SS2 = subset(goodSS2, cells=top)
SS2$color_by <- factor(SS2$predicted.clusters)
levels(SS2$color_by) = c("Surface\nCluster 1", "Surface\nCluster 2", "Surface\nCluster 3", "Center\nCluster 1")
dat <- data.frame(SS2@reductions$ref.umap@cell.embeddings) %>%
  rownames_to_column('Cells') %>%
  left_join(SS2@meta.data %>% rownames_to_column("Cells"))
dat_ss2 = dat

ss2 <-
  ggplot(dat %>% arrange(desc(color_by)), aes(x = refUMAP_1, y = refUMAP_2, col=color_by)) +
  geom_point(size=0.2) + scale_color_manual(name="", values=colo) +
  theme_classic() +
  labs(x="ref.UMAP1", y="ref.UMAP2", title="HGSOC Ascites Tumor Cells Isolated")  +
  theme(plot.title = element_text(size = 9, face = "plain")) +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme(legend.text=element_text(size = 7), legend.key.size = unit(0.7, 'cm')) +
  coord_cartesian(xlim = c(-12, 12), ylim=c(-7,6)) +
  theme(axis.title = element_text(size = 9))


tiff("Figures/Fig7ABC_11.20.22.tiff", width=10, height=2.5, compression = "lzw", res=300, units='in')
print(seep + x10 + ss2 +
  plot_layout(nrow=1) +
  #plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect"))
dev.off()

pdf("Figures/Fig7ABC_11.20.22.pdf", width=10, height=2.5)
print(seep + x10 + ss2 +
        plot_layout(nrow=1) +
        #plot_annotation(tag_levels = "A") +
        plot_layout(guides = "collect"))
dev.off()

tiff("Figures/Additional/UMAPtransferred.tiff", width=10, height=8, compression = "lzw", res=300, units='in')
print(((seep + geom_point(size=1) + facet_wrap(~long_source, nrow=1) + theme(strip.background = element_blank())) /
(x10 + geom_point(size=1) + facet_wrap(~Cell_type, nrow=1) + theme(strip.background = element_blank())) / (ss2 + geom_point(size=1) + facet_wrap(~Cell_type, nrow=1)+ theme(strip.background = element_blank())) )+
  #plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect"))
dev.off()


