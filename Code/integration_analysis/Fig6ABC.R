rm(list=ls())
setwd("C://PROJECTS/P2022/SEEP_Manuscript")

library(Seurat); library(tidyverse); library(SeuratObject)
library(ggplot2); library(patchwork)

data("SEEP_Separate")
data("SEEP_Regions")
colo = ggsci::pal_d3()(4)

so = sso2D

so$color_by = factor(so$long_source)
dat <- data.frame(so@reductions$umap@cell.embeddings) %>%
  rownames_to_column('Cells') %>%
  left_join(so@meta.data %>% rownames_to_column("Cells"))

model_source <-
ggplot(dat, aes(x = UMAP_1, y = UMAP_2, col=color_by)) +
  geom_point(size=0.05) + #scale_color_manual(name="", values=colo) +
  theme_classic() +
  labs(x="UMAP1", y="UMAP2", title="Model Source") +
  theme(plot.title = element_text(size = 10, face = "plain")) +
  guides(color = guide_legend(override.aes = list(size = 2, fill=NA))) +
  theme(legend.text=element_text(size = 7), legend.key.height = unit(0.5, 'cm'), legend.title = element_blank()) +
  #coord_cartesian(xlim = c(-12, 12), ylim=c(-7,6)) +
  theme(axis.title = element_text(size = 9), plot.title=element_text(hjust=0.5))

so$color_by = factor(so$long_layer)
dat <- data.frame(so@reductions$umap@cell.embeddings) %>%
  rownames_to_column('Cells') %>%
  left_join(so@meta.data %>% rownames_to_column("Cells"))

layer_source <-
  ggplot(dat, aes(x = UMAP_1, y = UMAP_2, col=color_by)) +
  geom_point(size=0.05) + #scale_color_manual(name="", values=colo) +
  theme_classic() +
  labs(x="UMAP1", y="UMAP2", title="Layer Source") +
  theme(plot.title = element_text(size = 10, face = "plain")) +
  guides(color = guide_legend(override.aes = list(size = 2, fill=NA))) +
  theme(legend.text=element_text(size = 7), legend.key.height = unit(0.5, 'cm'), legend.title = element_blank()) +
  #coord_cartesian(xlim = c(-12, 12), ylim=c(-7,6)) +
  theme(axis.title = element_text(size = 9), plot.title=element_text(hjust=0.5))

# integ
iso = iso2D

iso$color_by = factor(iso$long_source)
dat <- data.frame(iso@reductions$umap@cell.embeddings) %>%
  rownames_to_column('Cells') %>%
  left_join(iso@meta.data %>% rownames_to_column("Cells"))

integ_model <-
  ggplot(dat, aes(x = UMAP_1, y = UMAP_2, col=color_by)) +
  geom_point(size=0.05) + #scale_color_manual(name="", values=colo) +
  theme_classic() +
  labs(x="UMAP1", y="UMAP2", title="Integrated Model") +
  theme(plot.title = element_text(size = 10, face = "plain")) +
  guides(color = guide_legend(override.aes = list(size = 2, fill=NA))) +
  theme(legend.text=element_text(size = 7), legend.key.height = unit(0.5, 'cm'), legend.title = element_blank()) +
  #coord_cartesian(xlim = c(-12, 12), ylim=c(-7,6)) +
  theme(axis.title = element_text(size = 9), plot.title=element_text(hjust=0.5))

iso$color_by = factor(iso$long_layer)
dat <- data.frame(iso@reductions$umap@cell.embeddings) %>%
  rownames_to_column('Cells') %>%
  left_join(iso@meta.data %>% rownames_to_column("Cells"))

integ_layer <-
  ggplot(dat, aes(x = UMAP_1, y = UMAP_2, col=color_by)) +
  geom_point(size=0.05) + #scale_color_manual(name="", values=colo) +
  theme_classic() +
  labs(x="UMAP1", y="UMAP2", title="Integrated Model") +
  theme(plot.title = element_text(size = 10, face = "plain")) +
  guides(color = guide_legend(override.aes = list(size = 2, fill=NA))) +
  theme(legend.text=element_text(size = 7), legend.key.height = unit(0.5, 'cm'), legend.title = element_blank()) +
  #coord_cartesian(xlim = c(-12, 12), ylim=c(-7,6)) +
  theme(axis.title = element_text(size = 9), plot.title=element_text(hjust=0.5))

# arrows
df2 <- expand.grid(
  lineend = c('square'), linejoin = c('mitre'), stringsAsFactors = FALSE)
df2 <- data.frame(df2, y = 1)

arr1=ggplot(df2, aes(x = 0.7, y = 1, xend = 0.7, yend = 0, label = "")) +
  geom_segment(
    lineend = df2$lineend, linejoin = df2$linejoin,
    size = 0.8, arrow = arrow(type="closed", length = unit(0.2, "cm"))) +
  geom_text() + annotate(
    "text", label = "Seurat\nIntegration", x = 1, y = 0.5, size = 3, colour = "black",
  hjust=0.5, fontface='plain') +  xlim(0, 2) + ylim(-0.2, 1.1) + theme_void()

arr2=ggplot(df2, aes(x = 0.7+0.38, y = y, xend = 0.7+0.38, yend = 0, label = "")) +
  geom_segment(
    lineend = df2$lineend, linejoin = df2$linejoin,
    size = 0.8, arrow = arrow(type="closed", length = unit(0.2, "cm"))) +
  geom_text() + annotate(
    "text", label = "Seurat\nIntegration", x = 1+0.38, y = 0.5, size = 3, colour = "black",
    hjust=0.5, fontface='plain') +  xlim(0, 2) + ylim(-0.2, 1.1) + theme_void()

arr3=ggplot(df2, aes(x = 0.7+0.74, y = y, xend = 0.7+0.74, yend = 0, label = "")) +
  geom_segment(
    lineend = df2$lineend, linejoin = df2$linejoin,
    size = 0.8, arrow = arrow(type="closed", length = unit(0.2, "cm"))) +
  geom_text() + annotate(
    "text", label = "DA-seq\nClustering", x = 1+0.74, y = 0.5, size = 3, colour = "black",
    hjust=0.5, fontface='plain') +  xlim(0, 2) + ylim(-0.2, 1.1) + theme_void()


##################

iso = iso2D

iso$color_by = factor(iso$DA_regions_layer)
levels(iso$color_by) = c("UNC", "Surface\nCluster 1", "Surface\nCluster 2", "Surface\nCluster 3", "Center\nCluster 1")

dat <- data.frame(iso@reductions$umap@cell.embeddings) %>%
  rownames_to_column('Cells') %>%
  left_join(iso@meta.data %>% rownames_to_column("Cells"))
colo0 = "antiquewhite2" #'#ffffff00'
colos = c(colo0, ggsci::pal_d3()(max(iso$DA_regions)))
colo = colos[-1]

da_clust <-
  DimPlot(iso, repel=T,group.by = "color_by", cols=colos, reduction = 'umap', pt.size = 0.1, label=F, label.size = 3) + theme_classic() +
  labs(x="UMAP1", y="UMAP2", title="Differentially Abundant Subpopulations") +
  theme(plot.title = element_text(size = 10, face = "plain")) +
  guides(color = guide_legend(override.aes = list(size = 2, fill=NA))) +
  theme(legend.text=element_text(size = 7), legend.key.height = unit(0.7, 'cm'), legend.title = element_blank()) +
  #coord_cartesian(xlim = c(-12, 12), ylim=c(-7,6)) +
  theme(axis.title = element_text(size = 9), plot.title=element_text(hjust=0.5))

da_pred <-
  da_cells$pred.plot +
  theme_classic() +
  labs(x="UMAP1", y="UMAP2", title="Integrated Model") +
  theme(plot.title = element_text(size = 10, face = "plain")) +
  guides(fill = guide_legend(override.aes = list(size = 2, fill=NA))) +
  theme(legend.text=element_text(size = 7), legend.key.height = unit(0.4, 'cm'), legend.key.width = unit(0.3, 'cm')) +
  theme(axis.title = element_text(size = 9), plot.title=element_text(hjust=0.5))
da_pred$layers[[1]]$aes_params$size = 0.1


R1 = model_source + layer_source + da_pred
R2 = arr1 + arr2 + arr3 #+ #plot_layout(tag_level = "new")
R3 = integ_model + integ_layer + da_clust #+ plot_layout(tag_level = 'new')
tiff("Figures/Fig6ABC_11.20.22.tiff", width=11, height=5.5, compression = "lzw", res=300, units='in')
print((R1/R2/R3) + plot_layout(heights = unit(c(1, 1.5, 1), c('null', 'cm','null'))))  #+ plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size=24)))
  #plot_layout(guides = "collect")
dev.off()
pdf("Figures/Fig6ABC_11.20.22.pdf", width=11, height=5.5)
print((R1/R2/R3) + plot_layout(heights = unit(c(1, 1.5, 1), c('null', 'cm','null'))))  #+ plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size=24)))
#plot_layout(guides = "collect")
dev.off()

