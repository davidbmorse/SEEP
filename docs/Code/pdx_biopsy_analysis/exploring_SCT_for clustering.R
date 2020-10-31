library(sctransform)
biopsy <- PercentageFeatureSet(biopsy, pattern = "^MT-", col.name = "percent.mt")
biopsy <- SCTransform(biopsy, vars.to.regress = c("percent.mt"), verbose = TRUE, variable.features.n = NULL, variable.features.rv.th = 1.2, return.only.var.genes = FALSE)


 #additional 'variables to regress': "S.Score","G2M.Score"

biopsy <- RunPCA(biopsy, verbose = FALSE)
biopsy <- RunUMAP(biopsy, dims = 1:30, verbose = FALSE)

biopsy <- FindNeighbors(biopsy, dims = 1:30, verbose = FALSE)
biopsy <- FindClusters(biopsy, verbose = FALSE)
DimPlot(biopsy, label = TRUE) + NoLegend()
DimPlot(biopsy, reduction = "umap", group.by = 'orig.ident')

biopsy <- RunTSNE(biopsy, dims = 1:30, verbose = FALSE)
DimPlot(biopsy, reduction = "tsne", group.by = 'seurat_clusters')
DimPlot(biopsy, reduction = "tsne", group.by = 'orig.ident')


Idents(biopsy) <- "seurat_clusters"
heatmap_subset <- subset(biopsy, idents = c("4","0","3","1","2","5", "6", "7"))
my_levels <- c(0,5,7,2,3,1,4,6)
levels(heatmap_subset) <- my_levels

biopsy.markers <- FindAllMarkers(heatmap_subset, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.2)
biopsy.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

top20 <- biopsy.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
DoHeatmap(heatmap_subset, features = top20$gene, raster = FALSE, size = 6, disp.max = 2, disp.min = -2, slot = "scale.data") + scale_fill_gradientn(colors = c("dodgerblue", "black", "yellow"), na.value = "white")

write.csv(biopsy.markers, file = "cluster_markers_SCT_010_020.csv")






counts_per_cluster <- data.frame (cells = biopsy@assays$SCT@scale.data[2,], cluster = biopsy@meta.data$seurat_clusters, layer = biopsy@meta.data$orig.ident)

counts_per_cluster <- counts_per_cluster %>% 
  group_by(cluster, layer) %>%
  summarise(count = n())

show_percentages <- counts_per_cluster %>%
  group_by(cluster) %>%
  mutate(countT = sum(count)) %>%
  group_by(layer, add = TRUE) %>%
  mutate(per = paste0(round(100*count/countT,2),'%'))

counts_per_cluster$cluster <- factor(counts_per_cluster$cluster, levels = c(0,5,7,2,3,1,4,6))

cluster_by_layer <- ggplot(counts_per_cluster, aes(x = cluster, y = count, fill = layer)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_light() +
  labs(y = "percent layer", x = "clusters") +
  theme(axis.text=element_text(size=6), axis.title=element_blank(), legend.text=element_text(size=8)) +
  scale_fill_manual(values = c("coral3","green4", "skyblue3"))
cluster_by_layer
