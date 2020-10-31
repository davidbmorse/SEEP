#re-analysis on biopsy after subtracting mouse reads
library(Seurat)
library(dplyr)
library(Matrix)
library(Matrix.utils)
library(data.table)
library(ggplot2)
library(DropletUtils)
library(patchwork)
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/")
biopsy <- readRDS(file = "filtered_biopsy_SO_final_042220.rds")
biopsy.markers <- read.csv(file = "cluster_markers_025_025.csv", row.names = 1)
biopsy.markers.layers <- read.csv(file = "layer_markers_025_025.csv", row.names = 1)

#1 plot percent layers barplot -------
counts_per_cluster <- data.frame (cells = biopsy@assays$RNA@scale.data[2,], cluster = biopsy@meta.data$seurat_clusters, layer = biopsy@meta.data$orig.ident)

counts_per_cluster <- counts_per_cluster %>% 
  group_by(cluster, layer) %>%
  summarise(count = n())

show_percentages <- counts_per_cluster %>%
  group_by(cluster) %>%
  mutate(countT = sum(count)) %>%
  group_by(layer, add = TRUE) %>%
  mutate(per = paste0(round(100*count/countT,2),'%'))

counts_per_cluster$cluster <- factor(counts_per_cluster$cluster, levels = c(3,0,2,5,1,4))

cluster_by_layer <- ggplot(counts_per_cluster, aes(x = cluster, y = count, fill = layer)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_light() +
  labs(y = "percent layer", x = "clusters") +
  theme(axis.text=element_text(size=6), axis.title=element_blank(), legend.text=element_text(size=8)) +
  scale_fill_manual(values = c("coral3","green4", "skyblue3"))
cluster_by_layer

pdf("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/figures/cluster_by_layer.pdf", width = 2.5, height = 2)
cluster_by_layer
dev.off()


#2 plot chi-square analysis
library(vcd)
meta = data.frame(Layer = biopsy$orig.ident, Cluster=biopsy$seurat_clusters)
tab = table(meta$Layer, meta$Cluster)
names(dimnames(tab))=c("Layer","Cluster")
rownames(tab) = c('CENTER','MIDDLE','SURFACE')

tiff("figures/Biopsy_chi_square.tiff", res=300, units='in', wid=6, hei=8, compress='lzw')
assoc(tab, gp=shading_hcl, compress=T, xscale=0.7)
dev.off()

pdf("figures/Biopsy_chi_square.pdf", width = 5, height = 7, useDingbats=FALSE)
assoc(tab, gp=shading_hcl, compress=T, xscale=0.7)
dev.off()

stats = assocstats(tab)
V= stats$cramer

#3 tsne plots
tSNE_cluster <- DimPlot(biopsy, reduction = "tsne", group.by = "seurat_clusters", label = FALSE,
                 cols = "Set1", pt.size = .001)
tSNE_cluster <- tSNE_cluster + 
  theme(axis.text=element_text(size=6), axis.title=element_blank(), legend.position = "none") + 
  labs(x = "tSNE 1", y = "tSNE 2")
tSNE_cluster

tSNE_layer <- DimPlot(biopsy, reduction = "tsne", group.by = "orig.ident", 
                 label = FALSE, pt.size = .001)
tSNE_layer <- tSNE_layer + 
  theme(axis.text=element_text(size=6), axis.title=element_blank(), legend.position = "none") + 
  labs(x = "tSNE 1", y = "tSNE 2")
tSNE_layer

pdf("figures/Biopsy_tSNE_clusters.pdf", width = 1.875, height = 1.5, useDingbats = FALSE)
tSNE_cluster
dev.off()

pdf("figures/Biopsy_tSNE_layers.pdf", width = 1.875, height = 1.5, useDingbats = FALSE)
tSNE_layer
dev.off()

#4 UMAP plots (alternative)
umap_cluster <- DimPlot(biopsy, reduction = "umap", group.by = "seurat_clusters", label = FALSE,
                        cols = "Set1", pt.size = .001)
umap_cluster <- umap_cluster + 
  theme(axis.text=element_text(size=6), axis.title=element_blank(), legend.position = "none") + 
  labs(x = "umap 1", y = "umap 2")
umap_cluster

umap_layer <- DimPlot(biopsy, reduction = "umap", group.by = "orig.ident", 
                      label = FALSE, pt.size = .001)
umap_layer <- umap_layer + 
  theme(axis.text=element_text(size=6), axis.title=element_blank(), legend.position = "none") + 
  labs(x = "umap 1", y = "umap 2")
umap_layer

pdf("figures/Biopsy_umap_clusters.pdf", width = 1.875, height = 1.5)
umap_cluster
dev.off()

pdf("figures/Biopsy_umap_layers.pdf", width = 1.875, height = 1.5)
umap_layer
dev.off()


#5 heatmap
Idents(biopsy) <- "seurat_clusters"
heatmap_subset <- subset(biopsy, idents = c("3","2","5","1","4"))
my_levels <- c(3,2,5,1,4)

levels(heatmap_subset) <- my_levels
#identify markers
heatmap.markers <- FindAllMarkers(object = heatmap_subset, only.pos = TRUE, min.pct = 0.2, thresh.use = 0.2)
heatmap.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top20_heatmap <- heatmap.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
Hmap302514 <- DoHeatmap(heatmap_subset, features = top20_heatmap$gene, raster = FALSE, size = 6, disp.max = 2, disp.min = -2, slot = "scale.data") + scale_fill_gradientn(colors = c("dodgerblue", "black", "yellow"), na.value = "white")
Hmap302514

pdf("figures/Hmap32514_025_025.pdf",  width = 4, height = 5)
Hmap302514
dev.off()

Hmap302514_NO_ANNO <- DoHeatmap(heatmap_subset, features = top20_heatmap$gene, raster = FALSE, size = 6, disp.max = 2, disp.min = -2, slot = "scale.data", group.bar = FALSE) + scale_fill_gradientn(colors = c("dodgerblue", "grey10", "yellow"), na.value = "white") + NoLegend() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())
Hmap302514_NO_ANNO

pdf("figures/Hmap32514_025_025_NO_ANNO.pdf",  width = 2.1, height = 4.3)
Hmap302514_NO_ANNO
dev.off()



#6 layer vs layer scatter plots: plot a regression line and find the residuals
library(ggrepel)
require(stats)
require(stats)
biopsy <- SetIdent(biopsy, value = "orig.ident")
avg.biopsy.cells <- log1p(AverageExpression(biopsy, verbose = FALSE)$RNA)
avg.biopsy.cells$gene <- rownames(avg.biopsy.cells)
reg <- lm(center ~ surface, data = avg.biopsy.cells)
avg.biopsy.cells$res <- residuals(reg)
write.csv(avg.biopsy.cells, file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/biop_avgPDXcells.csv")


#surace genes - inflamation and EMT, center genes MTOR and hypoxia
inflamation <- c("NFKBIA", "BTG1")
EMT <- c("SAT1","VIM", "CXCL8")
ANTI_APOP <- c("MTRNR2L12", "MTRNR2L1", "MTRNR2L8", "MTRNR2L10")
hypoxia <- c("JUN", "DDIT4", "TPI1", "PCK1")
P53 <- c("NINJ1", "IRAK1", "CTSF")

#surface vs center
surf_cent <- ggplot(avg.biopsy.cells, aes(center, surface, label = gene)) + 
  geom_abline(intercept = 0.00, slope = 1.00, linetype = "dashed") + #intercept = 0.0089, slope = 1.0035 +
  geom_point(color = ifelse(avg.biopsy.cells$res < -0.3 | avg.biopsy.cells$res > 0.3, "red", "black"), size = .01, alpha = .7) +
  #SURFACE SIDE
  geom_label_repel(
    data          = subset(avg.biopsy.cells, rownames(avg.biopsy.cells) %in% inflamation),
    segment.size  = 0.2,
    segment.color = "grey50",
    fill          = "lightcoral",
    size          = 2,
    label.padding = 0.13,
    segment.alpha = 0.7,
    nudge_y       = 1.5,
    box.padding   = 0.1,
    point.padding = 0.1,
    direction     = "both",
    force         = 10
  ) +
  geom_label_repel(
    data          = subset(avg.biopsy.cells, rownames(avg.biopsy.cells) %in% EMT),
    segment.size  = 0.2,
    segment.color = "grey50",
    fill          = "palegreen",
    size          = 2,
    label.padding = 0.13,
    segment.alpha = 0.7,
    nudge_y       = 1.5,
    box.padding   = 0.1,
    point.padding = 0.1,
    direction     = "both",
    force         = 10
  ) +
  geom_label_repel(
    data          = subset(avg.biopsy.cells, rownames(avg.biopsy.cells) %in% ANTI_APOP),
    segment.size  = 0.2,
    segment.color = "grey50",
    fill          = "lightskyblue",
    size          = 2,
    label.padding = 0.13,
    segment.alpha = 0.7,
    nudge_y       = 1.5,
    box.padding   = 0.1,
    point.padding = 0.1,
    direction     = "both",
    force         = 10
  ) +
  #CENTER SIDE
  geom_label_repel(
    data          = subset(avg.biopsy.cells, rownames(avg.biopsy.cells) %in% hypoxia),
    segment.size  = 0.2,
    segment.color = "grey50",
    fill          = "navajowhite",
    size          = 2,
    label.padding = 0.13,
    segment.alpha = 0.7,
    nudge_x       = 1.5,
    #    xlim          = c(3,5),
    box.padding   = 0.1,
    point.padding = 0.1,
    direction     = "both",
    force         = 10
  ) +
  geom_label_repel(
    data          = subset(avg.biopsy.cells, rownames(avg.biopsy.cells) %in% P53),
    segment.size  = 0.2,
    segment.color = "grey50",
    fill          = "plum2",
    size          = 2,
    label.padding = 0.13,
    segment.alpha = 0.7,
    nudge_x       = 1.5,
    #    xlim          = c(1.5,3),
    box.padding   = 0.1,
    point.padding = 0.1,
    direction     = "both",
    force         = 10
  ) +
  scale_x_continuous(limits = c(NA, 5)) +
  scale_y_continuous(limits = c(NA, 5)) + 
  theme_classic(base_size = 12) +
  labs(y = "surface cells (log(1+x))", x = "center cells (log(1+x))") +
  theme(axis.text=element_text(size=6), axis.title=element_blank())

surf_cent

pdf("figures/Biopsy_surface_center.pdf", width = 1.7, height = 1.7)
surf_cent
dev.off()


#SI figures ----
#full heatmap

#plot for layers
biopsy.markers <- FindAllMarkers(biopsy, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.20)
biopsy.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
top20_heatmap <- biopsy.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)

#plot_heatmap
Hmap_full <- DoHeatmap(biopsy, features = top20_heatmap$gene, raster = TRUE, size = 6, disp.max = 2, disp.min = -2, slot = "scale.data") + scale_fill_gradientn(colors = c("dodgerblue", "grey10", "yellow"), na.value = "white")

pdf("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Manuscripts/paint_sort_sequencing/supplement/biopsy_filtered/Hmap_all_genes.pdf",  width = 6, height = 8)
Hmap_full
dev.off()

Hmap_no_annot <- DoHeatmap(biopsy, features = top20_heatmap$gene, raster = TRUE, size = 6, disp.max = 2, disp.min = -2, slot = "scale.data", group.bar = FALSE) + scale_fill_gradientn(colors = c("dodgerblue", "grey10", "yellow"), na.value = "white") + NoLegend() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank())
Hmap_no_annot

pdf("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Manuscripts/paint_sort_sequencing/supplement/biopsy_filtered/Hmap_all_genes_No_Anno.pdf",  width = 6, height = 8)
Hmap_no_annot
dev.off()







