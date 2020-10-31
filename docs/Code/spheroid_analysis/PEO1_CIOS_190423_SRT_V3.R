# intro, clustering, and vizualization --------------------------------------
#loading data and workspace
library(Seurat)
library(dplyr)
library(tidyr)
library(wesanderson)
library(Matrix)
library(data.table)
library(ggplot2)

CIOS <- Read10X("/Users/morsedb/Documents/Experiments/2017/4July_August/20170711_inDrop_spheroid_layers/data_analysis/counts_matrices/spheroid_CIOS_new")
spheroid_CIOS <- CreateSeuratObject(counts = CIOS, min.cells = 3, min.features = 200, project = "indrop_spheroid_CIOS")

#pre-processing
spheroid_CIOS[["percent.mt"]] <- PercentageFeatureSet(spheroid_CIOS, pattern = "^MT-")

VlnPlot(spheroid_CIOS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(spheroid_CIOS, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(spheroid_CIOS, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

#Filter
spheroid_CIOS <- subset(spheroid_CIOS, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 18)
spheroid_CIOS <- NormalizeData(spheroid_CIOS, normalization.method = "LogNormalize", scale.factor = 10000)
spheroid_CIOS <- FindVariableFeatures(spheroid_CIOS, verbose = TRUE, selection.method = "mean.var.plot", mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.1, Inf))
#removed from find variable gene command (5/9/19)
    #mean.function = ExpMean
    #dispersion.function = LogVMR
top10 <- head(VariableFeatures(spheroid_CIOS), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(spheroid_CIOS)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

spheroid_CIOS <- ScaleData(spheroid_CIOS, vars.to.regress = c("nCount_RNA", "percent.mt"))

#############

#processing
spheroid_CIOS <- RunPCA(spheroid_CIOS, features = VariableFeatures(object = spheroid_CIOS), do.print = TRUE, pcs.print = 1:5, genes.print = 5)
DimHeatmap(spheroid_CIOS, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(object = spheroid_CIOS)

#clustering & visualization
spheroid_CIOS <- FindNeighbors(spheroid_CIOS, dims = 1:12)
spheroid_CIOS <- FindClusters(spheroid_CIOS, resolution = 0.45)
head(Idents(spheroid_CIOS), 5)

#old clustering
  #spheroid_CIOS <- FindClusters(object = spheroid_CIOS, reduction.type = "pca", dims.use = 1:12, resolution = 0.4, print.output = 0, save.SNN = TRUE)
  #spheroid_CIOS <- RunTSNE(object = spheroid_CIOS, dims.use = 1:12, resolution = 0.4, do.fast = TRUE, perplexity = 100, exaggeration_factor=30)

#UMAP
spheroid_CIOS <- RunUMAP(spheroid_CIOS, dims = 1:20, n_neighbors = 30, min_dist = 0.0, metric = 'cosine')
DimPlot(spheroid_CIOS, reduction = "umap")
DimPlot(spheroid_CIOS, reduction = "umap", group.by = "orig.ident")

#tSNE
spheroid_CIOS <- RunTSNE(spheroid_CIOS, dims = 1:12, resolution = 0.45, do.fast = TRUE, perplexity = 100, exaggeration_factor=30)
plot1 <- DimPlot(spheroid_CIOS, reduction = "tsne", label = TRUE, cols = "Set1", pt.size = .001)
plot1 <- plot1 + theme(axis.text=element_text(size=6), axis.title=element_text(size=10)) + labs(x = "tSNE 1", y = "tSNE 2")
plot2 <- DimPlot(spheroid_CIOS, reduction = "tsne", group.by = "orig.ident", label = TRUE, pt.size = .001)
plot2 <- plot2 + theme(axis.text=element_text(size=6), axis.title=element_text(size=10)) + labs(x = "tSNE 1", y = "tSNE 2")

tsne_plot <- CombinePlots(plots = list(plot1, plot2), legend = 'none')
tsne_plot

pdf("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/reAnalysis_190412/plots/tSNE/tsne_spheroid2.pdf", width = 5.5, height = 2.5)
tsne_plot
dev.off()


#Save RDS file and provide option for loading
saveRDS(spheroid_CIOS, file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/reAnalysis_190412/data/spheroid_CIOS.rds")
#read back the RDS file
spheroid_CIOS <- readRDS(file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/reAnalysis_190412/data/spheroid_CIOS.rds")


#Find marker genes of each cluster
spheroid_CIOS.markers <- FindAllMarkers(object = spheroid_CIOS, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
spheroid_CIOS.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top10 <- spheroid_CIOS.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(spheroid_CIOS, features = top10$gene, raster = FALSE) + NoLegend()

#save full heatmap (clusters 0 - 6) gene lists
write.csv(spheroid_CIOS.markers, "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/reAnalysis_190412/plots/heatmaps/all_clusters_Markers.csv")
write.csv(top10, "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/reAnalysis_190412/plots/heatmaps/all_clusters_20_Markers.csv")

#try to subset to only clusters with heavy weighting towards a layer
#subset for clusters 4,5,6
heatmap_subset <- subset(spheroid_CIOS, idents = c("1", "3", "4", "5", "6"))
#set group ordering
my_levels <- c(3, 1, 4, 5, 6)
levels(heatmap_subset) <- my_levels
#identify markers
heatmap.markers <- FindAllMarkers(object = heatmap_subset, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
heatmap.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top20_heatmap <- heatmap.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
#plot heatmap
DoHeatmap(heatmap_subset, features = top20_heatmap$gene, raster = FALSE) + NoLegend()
#try plotting highly variable cells - how to find most variable cells in a cluster?
DoHeatmap(heatmap_subset, features = heatmap.markers$gene, raster = FALSE) + NoLegend()

#save heatmap markers and top 10
write.csv(heatmap.markers, "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/reAnalysis_190412/plots/heatmaps/Markers.csv")
write.csv(top20_heatmap, "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/reAnalysis_190412/plots/heatmaps/20_Markers.csv")

#analyse in STRING db and in reactome


#save spheroid markers and top 10
write.csv(spheroid_CIOS.markers, "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/reAnalysis_190412/plots/SpheroidMarkers.csv")
write.csv(top10, "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/reAnalysis_190412/plots/top10.csv")



#Find Markers between layers SET NEW IDENT-----------------------
layers.markers <- FindAllMarkers(spheroid_CIOS)
top10.layers <- layers.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
#plot heatmop of identifying genes
DoHeatmap(object = spheroid_CIOS, features = top10.layers$gene)
PCHeatmap(object = spheroid_CIOS, pc.use = 1, cells.use = 400, do.balanced = TRUE, label.columns = FALSE, num.genes = 30, remove.key = TRUE)
#compare principle components coloring by layer identity
PCAPlot(object = spheroid_CIOS, dim.1 = 1, dim.2 = 2, group.by = "orig.ident", pt.size = 0.1)


#Find Markers between clusters -----------------------------------------
clusters.markers <- FindAllMarkers(spheroid_CIOS)
top10.clusters <- clusters.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
#plot heatmop of identifying genes
DoHeatmap(object = spheroid_CIOS, features = top10.clusters$gene)
PCHeatmap(object = spheroid_CIOS, pc.use = 1, cells.use = 400, do.balanced = TRUE, label.columns = FALSE, num.genes = 30, remove.key = TRUE)
#compare principle components coloring by layer identity
PCAPlot(object = spheroid_CIOS, dim.1 = 1, dim.2 = 2, group.by = "orig.ident", pt.size = 0.1)


# UMAP visualization  --------------------------------------
library(reticulate)
spheroid_CIOS <- FindClusters(object = spheroid_CIOS, reduction.type = "pca", dims.use = 1:20, resolution = 0.6, print.output = 0, save.SNN = TRUE)

spheroid_CIOS <- RunUMAP(spheroid_CIOS, reduction.use = "pca", dims = 1:20, n_neighbors = 10, min_dist = 0.0, metric = 'cosine')

DimPlot(spheroid_CIOS, reduction.use = "umap", plot.title = "dims 20, metric cosine, min_dist 0.0, neigh 10", pt.size = .3)
DimPlot(spheroid_CIOS, reduction.use = "umap", group.by = "orig.ident", plot.title = "dims 20, metric cosine, min_dist 0.0, neigh 10", pt.size = .3)
#best so far is 20 dims, 15/10 neighboors, 0.0 min dist, and cosine or correlation metric, also euclidean distance metric



# comparison scatterplots  --------------------------------------
#reason for using log1p [this is log(1+x)]. http://onbiostatistics.blogspot.com/2012/05/logx1-data-transformation.html; also, https://scrnaseq-course.cog.sanger.ac.uk/website/cleaning-the-expression-matrix.html
library(ggrepel)
#1 stash cluster identities
spheroid_CIOS[["clusterID"]] <- Idents(object = spheroid_CIOS)
# or this to stash ident
spheroid_CIOS[["clusterID"]] <- spheroid_CIOS@meta.data$seurat_clusters
#set identity to layers (orig.ident)
spheroid_CIOS <- SetIdent(spheroid_CIOS, value = "orig.ident")
#set identity to clusters (clusterID)
spheroid_CIOS <- SetIdent(spheroid_CIOS, value = "clusterID")


#avg.CIOS.cells <- log1p(AverageExpression(spheroid_CIOS))
avg.CIOS.cells <- (AverageExpression(spheroid_CIOS))
avg.CIOS.cells <- log1p(avg.CIOS.cells$RNA)
avg.CIOS.cells$gene <- rownames(avg.CIOS.cells)
#plot a regression line and find the residuals (distance between the abline and the datapoints)
require(stats)
#plot regression
reg <- lm(C ~ S, data = avg.CIOS.cells)
#and residuals to the natural log expression data
avg.CIOS.cells$res <- residuals(reg)

setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/reAnalysis_190412/plots")

#surface/center genes to label
surf_AVGgenes <- c("IL1R2", "SOD2", "ECM1", "SLPI", "SERF2", "TNFAIP2")
cent_AVGgenes <- c("AKAP12", "KRT5", "PODXL", "CD24", "SP1", "INSIG1", "CCDC80")
#surface vs center
surf_cent <- ggplot(avg.CIOS.cells, aes(C, S, label = gene)) + 
  geom_abline(intercept = 0.00, slope = 1.00, linetype = "dashed") + #intercept = 0.0089, slope = 1.0035 +
  geom_point(color = ifelse(avg.CIOS.cells$res < -0.37 | avg.CIOS.cells$res > 0.37, "red", "black"), size = 1, alpha = .7) +
  geom_label_repel(data = subset(avg.CIOS.cells, rownames(avg.CIOS.cells) %in% surf_AVGgenes),
  #geom_label_repel(data = subset(avg.CIOS.cells, res < -0.47),
                   nudge_x = -.5,
                   nudge_y = 0.8,
                   size = 2,
                   box.padding = 0.25,
                   point.padding = 0.1,
                   force = 3,
                   segment.size  = 0.4,
                   direction = "y",
                   segment.alpha = 0.7,
                   segment.color = "grey50") +
  #geom_label_repel(data = subset(avg.CIOS.cells, res > 0.47),
  geom_label_repel(data = subset(avg.CIOS.cells, rownames(avg.CIOS.cells) %in% cent_AVGgenes),
                   nudge_x = 1.2,
                   size = 2,
                   box.padding = .5,
                   point.padding = 0.1,
                   force = 3,
                   segment.size  = 0.4,
                   segment.alpha = 0.7,
                   segment.color = "grey50") +
  #scale_x_continuous(expand = c(0.05, 0.05)) +
  scale_x_continuous(limits = c(NA, 5)) +
  scale_y_continuous(limits = c(NA, 5)) + 
  theme_classic(base_size = 12) +
  labs(y = "surface cells (log(1+x))", x = "center cells (log(1+x))") +
  theme(axis.title = element_text(size = 12, face = "bold"))

pdf("Layer_ScatterPlot/surface_center2.pdf", width = 2, height = 2)
surf_cent
dev.off()


#surface vs outer
surf_AVGgenesVout <- c("TUBA1B", "ACTN1", "ILF2", "RAP2B", "KRT13", "SFRP1", "PTMA", "AHNAK")
out_AVGgenesVsurf <- c("SLPI", "SOD2", "IFI27", "KLF10", "ECM1", "SAT1", "FTH1", "MIF", "HLA-C")
#calculate abline
require(stats)
reg_os <- lm(O ~ S, data = avg.CIOS.cells)
avg.CIOS.cells$res_os <- residuals(reg_os)
#plot
surf_out <- ggplot(avg.CIOS.cells, aes(O, S, label = gene)) + 
  geom_abline(intercept = 0.00, slope = 1.00, linetype = "dashed") + #intercept = 0.0089, slope = 1.0035
  geom_point(color = ifelse(avg.CIOS.cells$res_os < -0.37 | avg.CIOS.cells$res_os > 0.37, "red", "black"),  size = 1, alpha = .7) +
  #geom_label_repel(data = subset(avg.CIOS.cells, res_os < -0.37),
  geom_label_repel(data = subset(avg.CIOS.cells, rownames(avg.CIOS.cells) %in% out_AVGgenesVsurf),
                   nudge_x = -.5,
                   nudge_y = 0.8,
                   size = 2,
                   box.padding = 0.25,
                   point.padding = 0.1,
                   force = 3,
                   segment.size  = 0.4,
                   direction = "y",
                   segment.alpha = 0.7,
                   segment.color = "grey50") +
  geom_label_repel(data = subset(avg.CIOS.cells, rownames(avg.CIOS.cells) %in% surf_AVGgenesVout),
                   nudge_x = 1.2,
                   size = 2,
                   box.padding = .5,
                   point.padding = 0.1,
                   force = 3,
                   segment.size  = 0.4,
                   segment.alpha = 0.7,
                   segment.color = "grey50") +
  #scale_x_continuous(expand = c(0.05, 0.05)) +
  scale_x_continuous(limits = c(NA, 5)) +
  scale_y_continuous(limits = c(NA, 5)) + 
  theme_classic(base_size = 12) +
  labs(y = "surface cells (log(1+x))", x = "outer cells (log(1+x))") +
  theme(axis.title = element_text(size = 12, face = "bold"))
#save
pdf("surface_outer.pdf", width = 3, height = 3)
surf_out
dev.off()




#surface vs inner
#calculate abline
require(stats)
reg_si <- lm(I ~ S, data = avg.CIOS.cells)
avg.CIOS.cells$res_si <- residuals(reg_si)
#plot
ggplot(avg.CIOS.cells, aes(I, S, label = gene)) + 
  geom_point(color = ifelse(avg.CIOS.cells$res_si < -0.47 | avg.CIOS.cells$res_si > 0.47, "red", "black")) +
  geom_label_repel(data = subset(avg.CIOS.cells, res_si < -0.47),
                   nudge_x = -.5,
                   nudge_y = 0.5,
                   box.padding = 0.5,
                   point.padding = 0.35,
                   force = 3,
                   segment.size  = 0.4,
                   direction = "y",
                   segment.alpha = 0.3,
                   segment.color = "grey50") +
  geom_label_repel(data = subset(avg.CIOS.cells, res_si > 0.47),
                   nudge_x = 0.8,
                   box.padding = 0.5,
                   point.padding = 0.35,
                   force = 3,
                   segment.size  = 0.4,
                   segment.alpha = 0.3,
                   segment.color = "grey50") +
  scale_x_continuous(expand = c(0.05, 0.05)) +
  scale_y_continuous(limits = c(NA, 6)) + 
  theme_classic(base_size = 12) +
  geom_abline(intercept = 0.0089, slope = 1.0035, linetype = "dashed") + 
  labs(y = "surface cells (ln)", x = "inner cells (ln)") +
  theme(axis.title = element_text(size = 14, face = "bold"))



#center vs inner
#calculate abline
require(stats)
reg_ci <- lm(C ~ I, data = avg.CIOS.cells)
avg.CIOS.cells$res_ci <- residuals(reg_ci)
#plot
ggplot(avg.CIOS.cells, aes(C, I, label = gene)) + 
  geom_point(color = ifelse(avg.CIOS.cells$res_ci < -0.47 | avg.CIOS.cells$res_ci > 0.47, "red", "black")) +
  geom_label_repel(data = subset(avg.CIOS.cells, res_ci < -0.47),
                   nudge_x = -.5,
                   nudge_y = 0.5,
                   box.padding = 0.5,
                   point.padding = 0.35,
                   force = 3,
                   segment.size  = 0.4,
                   direction = "y",
                   segment.alpha = 0.3,
                   segment.color = "grey50") +
  geom_label_repel(data = subset(avg.CIOS.cells, res_ci > 0.47),
                   nudge_x = 0.8,
                   box.padding = 0.5,
                   point.padding = 0.35,
                   force = 3,
                   segment.size  = 0.4,
                   segment.alpha = 0.3,
                   segment.color = "grey50") +
  scale_x_continuous(expand = c(0.05, 0.05)) +
  scale_y_continuous(limits = c(NA, 6)) + 
  theme_classic(base_size = 12) +
  geom_abline(intercept = 0.0089, slope = 1.0035, linetype = "dashed") + 
  labs(y = "inner cells (ln)", x = "center cells (ln)") +
  theme(axis.title = element_text(size = 14, face = "bold"))









# Feature plots (tSNE and UMAP) --------------------------------------
#tsne
#surface
FeaturePlot(object = spheroid_CIOS, features.plot = c("SOD2", "IFITM1", "NFKBIA", "RSAD2", "HLA-B", "TXNIP", "HLA-B", "IFI27", "PSME1"), cols.use = c("grey", "red"), reduction.use = "tsne")
#center
FeaturePlot(object = spheroid_CIOS, features.plot = c("MT2A", "BTG1", "NCOA7", "PFKP", "LY6E", "SAMD9L", "EIF4E3", "NAMPT", "SSPN"), cols.use = c("grey", "red"), reduction.use = "tsne")

#UMAP
#enriched in center
FeaturePlot(object = spheroid_CIOS, features.plot = c("PARP1", "STAT1"), reduction.use = "umap", cols.use = c("blue", "green"), dark.theme = TRUE, pt.size = 0.5)
#enriched in surface
FeaturePlot(object = spheroid_CIOS, features.plot = c("HLA-A",	"SLPI",	"ADIRF",	"SERF2"), reduction.use = "umap", cols.use = c("blue", "green"), dark.theme = TRUE, pt.size = 0.5)











# clusters by layer composition  --------------------------------------
# old organizing -------
# test <- data.frame (cells = spheroid_CIOS@scale.data[2,], cluster = spheroid_CIOS@meta.data$clusterID, layer = spheroid_CIOS@meta.data$orig.ident)
# 
# counts_per_cluster <- 
#   test %>% 
#   group_by(cluster, layer) %>%
#   summarise(count = n())

#new cluster seurat V3 -------
test <- data.frame (cells = spheroid_CIOS@assays$RNA@scale.data[2,], cluster = spheroid_CIOS@meta.data$clusterID, layer = spheroid_CIOS@meta.data$orig.ident)
counts_per_cluster <- 
  test %>% 
  group_by(cluster, layer) %>%
  summarise(count = n())
counts_per_cluster$cluster <- factor(counts_per_cluster$cluster, levels = c(3,4,2,0,1,5,6))

cluster_by_layer <- ggplot(counts_per_cluster, aes(x = cluster, y = count, fill = layer)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "percent layer", x = "clusters") +
  theme(axis.title = element_text(size = 12, face = "bold")) +
  scale_fill_manual(values = c("blue2","green4", "yellow2", "red2"))

pdf("cluster_by_layer.pdf", width = 3.5, height = 2)
cluster_by_layer
dev.off()



# cell cycle attributes per layer --------------------------------------
#make a gene list containing cell cycle genes:
cell_cycle_genes <- read.table("/Users/morsedb/Documents/Experiments/2017/4July_August/20170711_inDrop_spheroid_layers/data_analysis/NCATS_analysis_171213/regev_lab_cell_cycle_genes.txt")
cell_cycle_gene_list <- cell_cycle_genes$V1
s.genes <- cell_cycle_gene_list[1:43]
g2m.genes <- cell_cycle_gene_list[44:97]

#Assign a cell cycle score 
spheroid_CIOS <- CellCycleScoring(object = spheroid_CIOS, s.features = s.genes, g2m.features =  g2m.genes, set.ident = TRUE)
# view cell cycle scores and phase assignments
head(x = spheroid_CIOS@meta.data)
RidgePlot(object = spheroid_CIOS, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), nCol = 2)

#reset cell identity
spheroid_CIOS <- SetIdent(spheroid_CIOS, value = "orig.ident")
#spheroid_CIOS = SetAllIdent(spheroid_CIOS, "clusterID")

#create a table describing cell cycle position of each layer
layers_cc <- tibble(spheroid_CIOS@active.ident, spheroid_CIOS@meta.data$Phase, .name_repair = c("check_unique", "unique", "universal", "minimal"))
layer_list <- spheroid_CIOS@active.ident
phase_list <- spheroid_CIOS@meta.data$Phase
layers_cc <- tibble(layer_list, phase_list)

#summarise counts for each layers
phase_intensity <- 
  layers_cc %>% 
  group_by(layer_list, phase_list) %>% 
  summarise(n())

#spread this data
cc_phase <-
  spread(phase_intensity, phase_list, 'n()')
#plot
cc_phase %>% 
  gather(key, value, -layer_list) %>% 
  ggplot(aes(layer_list, value, fill = key)) +
  geom_col(position="fill") +
  theme_bw() +
  labs(fill = "Cell Cycle Phases") +
  xlab("layer") +
  ylab("fraction of cells")


#rename G1 phase to include it with S phase
cc <- phase_intensity
cc$phase_list[cc$phase_list == "G1"] <- "S"
cc <- cc %>% group_by(layer_list, phase_list) %>% summarise_each(list(sum))
cc <- spread(cc, phase_list, 'n()')
#plot cycling vs non cycling cells
categories <- c("G2/M phases", "G1 & S phases")
cols <- c( "#56B4E9", "#E69F00")
# more colors: "#5121E0" #105415, "#E61515", "#5121E0"

cc %>% 
  gather(key, value, -layer_list) %>% 
  ggplot(aes(layer_list, value, fill = key)) +
  geom_col(position="fill") +
  #  coord_flip() +
  theme_bw() +
  theme(legend.position = "right") +
  scale_fill_manual(labels = categories, values = cols) +
  labs(fill = "Cell Cycle Phases") +
  xlab("layer") +
  ylab("fraction of cells")







#UMAP
spheroid_CIOS <- RunUMAP(spheroid_CIOS, dims = 1:20, n_neighbors = 10, min_dist = 0.0, metric = 'cosine')
DimPlot(spheroid_CIOS, reduction = "umap", group.by = "clusterID")
DimPlot(spheroid_CIOS, reduction = "umap", group.by = "orig.ident")

#tSNE
spheroid_CIOS <- RunTSNE(spheroid_CIOS, dims = 1:12, resolution = 0.4, do.fast = TRUE, perplexity = 100, exaggeration_factor=30)
plot1 <- DimPlot(spheroid_CIOS, reduction = "tsne", label = TRUE)
plot2 <- DimPlot(spheroid_CIOS, reduction = "tsne", group.by = "orig.ident", label = TRUE)
tsne_plot <- CombinePlots(plots = list(plot1, plot2), legend = 'none')

pdf("tsne_spheroid.pdf", width = 7, height = 3)
tsne_plot
dev.off()









