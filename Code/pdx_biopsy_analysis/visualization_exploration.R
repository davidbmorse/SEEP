library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(wesanderson)

#dowload and load Seurat v3
#install.packages('Seurat')
#dowload and load Seurat v2
#install.packages('devtools')
#devtools::install_version(package = 'Seurat', version = package_version('2.3.4'))

library(Seurat)

#tSNE plotting variation

#1 load RDS file
PDX_biopsy2_3 <- readRDS("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_exports/r_data/4th_merge_mouse/PDX_biopsy2_3_MERGE_final.rds")

PrintFindClustersParams(object = PDX_biopsy2_3)

#2 confirm number of PC's used to cluster
PDX_biopsy2_3 <- FindClusters(object = PDX_biopsy2_3, reduction.type = "pca", dims.use = 1:5, resolution = c(0.6), print.output = 0, save.SNN = TRUE)

#reassigning resolution - run if using same dimensions as previously run
PDX_biopsy2_3 <- FindClusters(PDX_biopsy2_3, resolution = 0.8, print.output = FALSE)


#run tSNE
#extra parameters for tSNE: https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html  learning_rate=# doesn't seem to do anything
# can we use the early_exaggeration parameter?
PDX_biopsy2_3 <- RunTSNE(object = PDX_biopsy2_3, dims.use = 1:5, perplexity=200, exaggeration_factor=30)
#exaggeration factor of 24, 20, was good

#11/11/18 perplexity looks good at 200/300

#plot tSNE
TSNEPlot(object = PDX_biopsy2_3, plot.title = "perplexity 200, dims 5, res 0.8, e_f 30")
TSNEPlot(object = PDX_biopsy2_3, group.by = "orig.ident", plot.title = "perplexity 200, dims 5, res 0.8, e_f 30")
PrintTSNEParams(object = PDX_biopsy2_3, raw = TRUE)

plot1 <- DimPlot(PDX_biopsy2_3, reduction = "tsne", label = TRUE)
plot2 <- DimPlot(PDX_biopsy2_3, reduction = "tsne", group.by = "orig.ident", label = TRUE)
tsne_plot <- CombinePlots(plots = list(plot1, plot2), legend = 'none')

pdf("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PDX_biopsy2_3sComb_190226/R_exports/quality_control/layers_merged/tSNE/tsne_PDX_biopsy2_3.pdf", width = 7, height = 3.5)
tsne_plot
dev.off()





#UMAP exploration: 

library(reticulate)

PDX_biopsy2_3 <- FindClusters(object = PDX_biopsy2_3, reduction.type = "pca", dims.use = 1:20, resolution = 0.6, print.output = 0, save.SNN = TRUE)

#run UMAP
PDX_biopsy2_3 <- RunUMAP(PDX_biopsy2_3, reduction.use = "pca", dims.use = 1:20, n_neighbors = 10, min_dist = 0.0, metric = 'cosine')

#plot
DimPlot(PDX_biopsy2_3, reduction.use = "umap", plot.title = "dims 20, metric cosine, min_dist 0.0, neigh 10", pt.size = .1)
DimPlot(PDX_biopsy2_3, reduction.use = "umap", group.by = "orig.ident", plot.title = "dims 20, metric cosine, min_dist 0.0, neigh 10", pt.size = .1)
#best so far is 20 dims, 15/10 neighboors, 0.0 min dist, and cosine or correlation metric, also euclidian distance metric


#Finding differentially expressed genes (cluster biomarkers)
# find all markers of cluster 1
cluster0.markers <- FindMarkers(object = PDX_biopsy2_3, ident.1 = 0, min.pct = 0.25)
print(x = head(x = cluster0.markers, n = 20))
# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster3.markers <- FindMarkers(object = PDX_biopsy2_3, ident.1 = 3, ident.2 = c(0,5), min.pct = 0.25)
print(x = head(x = cluster3.markers, n = 10))
# find markers for every cluster compared to all remaining cells, report only the positive ones
PDX_biopsy2_3.markers <- FindAllMarkers(object = PDX_biopsy2_3, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
PDX_biopsy2_3.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
top10 <- PDX_biopsy2_3.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)






#Take average expression of layers and generate the scatter plots, highlighting genes that exhibit dramatic differences
#center
#1 stash cluster identities
PDX_biopsy2_3[["clusterID"]] <- Idents(object = PDX_biopsy2_3)
#set identity to layers (orig.ident)
Idents(PDX_biopsy2_3) <- "orig.ident"

#full PDX_biopsy2_3
avg.PDX_biopsy2_3.cells <- log1p(AverageExpression(PDX_biopsy2_3, verbose = FALSE)$RNA)
avg.PDX_biopsy2_3.cells$gene <- rownames(avg.PDX_biopsy2_3.cells)


#plot average expression differences ----------
library(ggplot2)
library(ggrepel)

#plot a regression line and find the residuals (distance between the abline and the datapoints)
require(stats)
#plot regression
reg <- lm(c ~ s, data = avg.PDX_biopsy2_3.cells)
#and residuals to the natural log expression data
avg.PDX_biopsy2_3.cells$res <- residuals(reg)

### geom_label_repel
#natural label position
# ggplot(avg.PDX_biopsy2_3.cells, aes(c, s, label = gene)) + 
#   geom_point(color = ifelse(avg.PDX_biopsy2_3.cells$res < -0.38 | avg.PDX_biopsy2_3.cells$res > 0.38, "red", "black")) +
#   geom_label_repel(data = subset(avg.PDX_biopsy2_3.cells, res < -0.45 | res > 0.5),
#                   box.padding = 0.5,
#                   point.padding = 0.35,
#                   force = 1,
#                   segment.size  = 0.5,
#                   segment.color = "grey50") +
#   scale_x_continuous(expand = c(0.05, 0.05)) +
#   scale_y_continuous(limits = c(NA, 6)) + 
#   theme_classic(base_size = 12) +
#   geom_abline(intercept = 0.0089, slope = 1.0035, linetype = "dashed")


#surface vs center
ggplot(avg.PDX_biopsy2_3.cells, aes(c, s, label = gene)) + 
  geom_point(color = ifelse(avg.PDX_biopsy2_3.cells$res < -0.4 | avg.PDX_biopsy2_3.cells$res > 0.4, "red", "black")) +
  geom_label_repel(data = subset(avg.PDX_biopsy2_3.cells, res < -0.5),
                   nudge_x = -.5,
                   nudge_y = 0.5,
                   box.padding = 0.5,
                   point.padding = 0.35,
                   force = 3,
                   segment.size  = 0.4,
                   direction = "y",
                   segment.alpha = 0.3,
                   segment.color = "grey50") +
  geom_label_repel(data = subset(avg.PDX_biopsy2_3.cells, res > 0.48),
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
  labs(y = "surface cells (ln)", x = "center cells (ln)") +
  theme(axis.title = element_text(size = 14, face = "bold"))


#surface vs middle
  #calculate abline
require(stats)
reg_sm <- lm(m ~ s, data = avg.PDX_biopsy2_3.cells)
avg.PDX_biopsy2_3.cells$res_sm <- residuals(reg_sm)
  #plot
ggplot(avg.PDX_biopsy2_3.cells, aes(m, s, label = gene)) + 
  geom_point(color = ifelse(avg.PDX_biopsy2_3.cells$res_sm < -0.35 | avg.PDX_biopsy2_3.cells$res_sm > 0.35, "red", "black")) +
  geom_label_repel(data = subset(avg.PDX_biopsy2_3.cells, res_sm < -0.4),
                   nudge_x = -.5,
                   nudge_y = 0.5,
                   box.padding = 0.5,
                   point.padding = 0.35,
                   force = 3,
                   segment.size  = 0.4,
                   direction = "y",
                   segment.alpha = 0.3,
                   segment.color = "grey50") +
  geom_label_repel(data = subset(avg.PDX_biopsy2_3.cells, res_sm > 0.4),
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
  labs(y = "surface cells (ln)", x = "middle cells (ln)") +
  theme(axis.title = element_text(size = 14, face = "bold"))



#center vs middle
#calculate abline
require(stats)
reg_cm <- lm(c ~ m, data = avg.PDX_biopsy2_3.cells)
avg.PDX_biopsy2_3.cells$res_cm <- residuals(reg_cm)
#plot
ggplot(avg.PDX_biopsy2_3.cells, aes(c, m, label = gene)) + 
  geom_point(color = ifelse(avg.PDX_biopsy2_3.cells$res_cm < -0.35 | avg.PDX_biopsy2_3.cells$res_cm > 0.35, "red", "black")) +
  geom_label_repel(data = subset(avg.PDX_biopsy2_3.cells, res_cm < -0.35),
                   nudge_x = -.5,
                   nudge_y = 0.5,
                   box.padding = 0.5,
                   point.padding = 0.35,
                   force = 3,
                   segment.size  = 0.4,
                   direction = "y",
                   segment.alpha = 0.3,
                   segment.color = "grey50") +
  geom_label_repel(data = subset(avg.PDX_biopsy2_3.cells, res_cm > 0.35),
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
  labs(y = "middle cells (ln)", x = "center cells (ln)") +
  theme(axis.title = element_text(size = 14, face = "bold"))





#Find Markers between layers ------
layers.markers <- FindAllMarkers(PDX_biopsy2_3)
top10.layers <- layers.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)

#plot heatmop of identifying genes
#DoHeatmap(object = PDX_biopsy2_3, genes.use = top10.layers$gene, slim.col.label = TRUE, remove.key = TRUE, num.genes = 30)
DoHeatmap(object = PDX_biopsy2_3, features = top10.layers$gene)

PCHeatmap(object = PDX_biopsy2_3, dims = 1, cells = 400, balanced = TRUE, nfeatures = 30)
#compare principle components coloring by layer identity
PCAPlot(object = PDX_biopsy2_3, dim.1 = 1, dim.2 = 2, group.by = "orig.ident", pt.size = 0.5)

#Find % of layer in each luster --------------------------------------
test <- data.frame (cells = PDX_biopsy2_3@assays$RNA@scale.data[2,], cluster = PDX_biopsy2_3@meta.data$clusterID, layer = PDX_biopsy2_3@meta.data$orig.ident)

counts_per_cluster <- 
  test %>% 
  group_by(cluster, layer) %>%
  summarise(count = n())

#SHOW PERCENTAGES SEPT 13, 2019
show_percentages <- counts_per_cluster %>%
  group_by(cluster) %>%
  mutate(countT= sum(count)) %>%
  group_by(layer, add=TRUE) %>%
  mutate(per=paste0(round(100*count/countT,2),'%'))

counts_per_cluster$cluster <- factor(counts_per_cluster$cluster, levels = c(1,2,4,3,0,5))

ggplot(counts_per_cluster, aes(x = cluster, y = count, fill = layer)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(y = "percent layer per cluster", x = "clusters ordered by composition of surface layer") +
  theme(axis.title = element_text(size = 12, face = "bold")) +
  scale_fill_manual(values = c("blue2","green4", "red2"))





#find cell cycle differences between each layer ------------------------
#make a gene list containing cell cycle genes:
cell_cycle_genes <- read.table("/Users/morsedb/Documents/Experiments/2017/4July_August/20170711_inDrop_spheroid_layers/data_analysis/NCATS_analysis_171213/regev_lab_cell_cycle_genes.txt")
cell_cycle_gene_list <- cell_cycle_genes$V1
s.genes <- cell_cycle_gene_list[1:43]
g2m.genes <- cell_cycle_gene_list[44:97]

#Assign a cell cycle score 
PDX_biopsy2_3 <- CellCycleScoring(object = PDX_biopsy2_3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# view cell cycle scores and phase assignments
head(x = PDX_biopsy2_3@meta.data)
RidgePlot(object = PDX_biopsy2_3, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), nCol = 2)

#reset cell identity
Idents(PDX_biopsy2_3) <- "orig.ident"
#Idents(PDX_biopsy2_3) <- "clusterID"

#create a table describing cell cycle position of each layer
layers_cc <- tibble(PDX_biopsy2_3@active.ident, PDX_biopsy2_3@meta.data$Phase, .name_repair = c("check_unique", "unique", "universal", "minimal"))
layer_list <- PDX_biopsy2_3@active.ident
phase_list <- PDX_biopsy2_3@meta.data$Phase
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
  labs(fill = "Cell Cycle Phases")


#rename G1 phase to include it with S phase
cc <- phase_intensity
cc$phase_list[cc$phase_list == "G1"] <- "S"
cc <- cc %>% group_by(layer_list, phase_list) %>% summarise_each(list(sum))
cc <- spread(cc, phase_list, 'n()')
#plot cycling vs non cycling cells
categories <- c("G1 & S phases", "G2/M phases")
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


#regress out cell cycle differences ---------------------
#scale data agaoin
PDX_biopsy2_3 <- PDX_biopsy2_32

PDX_biopsy2_3 <- RunPCA(object = PDX_biopsy2_3, pc.genes = PDX_biopsy2_3@var.genes, pcs.print = 1:4, genes.print = 10)
PDX_biopsy2_3 <- CellCycleScoring(object = PDX_biopsy2_3, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = TRUE)
PDX_biopsy2_3 <- ScaleData(object = PDX_biopsy2_3, display.progress = FALSE)

PDX_biopsy2_3 <- RunPCA(object = PDX_biopsy2_3, pc.genes = c(s.genes, g2m.genes), do.print = FALSE, pcs.compute = 5)
PCAPlot(object = PDX_biopsy2_3)

all(g2m.genes %in% rownames(PDX_biopsy2_3@data))
g2m.genes[!(g2m.genes %in% rownames(PDX_biopsy2_3@scale.data))]
#I am missing MLF1IP from the s.genes, so lets remove this from the s.genes list
s.genes <- s.genes[-15]



#regressing out the difference between the G2M and S phase scores. Signals separating non-cycling cells and cycling cells will be maintained, but differences in cell cycle phase amongst proliferating cells, will be regressed out of the data
PDX_biopsy2_3@meta.data$CC.Difference <- PDX_biopsy2_3@meta.data$S.Score - PDX_biopsy2_3@meta.data$G2M.Score
marrow <- ScaleData(object = marrow, vars.to.regress = "CC.Difference", display.progress = FALSE)

# cell cycle effects strongly mitigated in PCA
PDX_biopsy2_3 <- RunPCA(object = PDX_biopsy2_3, pc.genes = PDX_biopsy2_3@var.genes, genes.print = 10)

#check out plotting differences:
#cluster
PDX_biopsy2_3 <- FindClusters(object = PDX_biopsy2_3, reduction.type = "pca", dims.use = 1:20, resolution = 0.6, print.output = 0, save.SNN = TRUE)
#run UMAP
PDX_biopsy2_3 <- RunUMAP(PDX_biopsy2_3, reduction.use = "pca", dims.use = 1:20, n_neighbors = 10, min_dist = 0.0, metric = 'cosine')
#plot
DimPlot(PDX_biopsy2_3, reduction.use = "umap", plot.title = "dims 20, metric cosine, min_dist 0.0, neigh 10", pt.size = .1)
DimPlot(PDX_biopsy2_3, reduction.use = "umap", group.by = "orig.ident", plot.title = "dims 20, metric cosine, min_dist 0.0, neigh 10", pt.size = .1)

#removing cell cycle variation through subtraction of cell cycle phase doesnt seem to affect clustering 














#old notes -------------------


#parameters that work well :) Nov 12th, 2018
# PDX_biopsy2_3 <- FindClusters(object = PDX_biopsy2_3, reduction.type = "pca", dims.use = 1:8, resolution = 0.6, print.output = 0, save.SNN = TRUE)
# 
# PDX_biopsy2_3 <- RunUMAP(PDX_biopsy2_3, reduction.use = "pca", dims.use = 1:8, n_neighbors = 100, min_dist = 0.0, metric = 'hamming')
# 
# DimPlot(PDX_biopsy2_3, reduction.use = "umap")
# DimPlot(PDX_biopsy2_3, reduction.use = "umap", group.by = "orig.ident")



# try UMAP on unclustered data:
PDX_biopsy2_3 <- RunUMAP(PDX_biopsy2_3, genes.use = PDX_biopsy2_3@var.genes, n_neighbors = 200, min_dist = 0.5, metric = 'hamming')
DimPlot(PDX_biopsy2_3, reduction.use = "umap")
DimPlot(PDX_biopsy2_3, reduction.use = "umap", group.by = "orig.ident")



# metrics for filtered cells:
table(PDX_biopsy2_3@ident)
# C    M    S          Total:
# 5850 3369 3472      12,691 cells


  
  
ggplot(data=nUMI, aes(x = nUMI)) + geom_histogram()

  







#old code:--------------------
CDH1), vimentin (VIM), N-cadherin (CDH2), and/or fibronectin 1 (FN1)
#0
FeaturePlot(object = PDX_biopsy2_3, features.plot = c("IFI27", "MTRNR2L1"), reduction.use = "umap", cols.use = c("blue", "green"), dark.theme = TRUE, pt.size = 0.5)
#1
FeaturePlot(object = PDX_biopsy2_3, features.plot = c("S100A2", "MT-RNR2"), reduction.use = "umap", cols.use = c("blue", "green"), dark.theme = TRUE, pt.size = 0.5)
#2
FeaturePlot(object = PDX_biopsy2_3, features.plot = c("RPL41", "MT-ND2"), reduction.use = "umap", cols.use = c("blue", "green"), dark.theme = TRUE, pt.size = 0.5)
FeaturePlot(object = PDX_biopsy2_3, features.plot = c("ANXA4", "PI3"), reduction.use = "umap", cols.use = c("blue", "green"), dark.theme = TRUE, pt.size = 0.5)
#3
FeaturePlot(object = PDX_biopsy2_3, features.plot = c("TUBA1B", "UBE2C"), reduction.use = "umap", cols.use = c("blue", "green"), dark.theme = TRUE, pt.size = 0.5)
#4 MUCI forms a mucus barrier on epitelia proteins! seems cool to be on the surface and middle cells, but not the center ones!

FeaturePlot(object = PDX_biopsy2_3, features.plot = c("MUC1", "MUC16"), reduction.use = "umap", cols.use = c("blue", "green"), dark.theme = TRUE, pt.size = 0.5)
#5
FeaturePlot(object = PDX_biopsy2_3, features.plot = c("CRIP1", "TMSB4X"), reduction.use = "umap", cols.use = c("blue", "green"), dark.theme = TRUE, pt.size = 0.5)
FeaturePlot(object = PDX_biopsy2_3, features.plot = c("CDH1", "CDC42"), reduction.use = "umap", cols.use = c("blue", "green"), dark.theme = TRUE, pt.size = 0.5)
#6
FeaturePlot(object = PDX_biopsy2_3, features.plot = c("MT-ATP6", "MT-CYB"), reduction.use = "umap", cols.use = c("green", "blue"), dark.theme = TRUE, max.cutoff = 'q8', pt.size = 0.5)
#7
FeaturePlot(object = PDX_biopsy2_3, features.plot = c("KIFC1"), reduction.use = "umap", cols.use = c("blue", "green"), dark.theme = TRUE, pt.size = 0.5)

#enriched in center
FeaturePlot(object = PDX_biopsy2_3, features.plot = c("PARP1", "STAT1"), reduction.use = "umap", cols.use = c("blue", "green"), dark.theme = TRUE, pt.size = 0.5)

#enriched in surface
FeaturePlot(object = PDX_biopsy2_3, features.plot = c("CTGF",	"APEX1",	"HSPA5",	"WNT7A"), reduction.use = "umap", cols.use = c("blue", "green"), dark.theme = TRUE, pt.size = 0.5)

"BBC3",	"CCND1",	"CCNE1",	"CDK1",	"CDKN2A",	"CYCS",	"RRM2",	"SESN2",	"SFN",

FeaturePlot(object = PDX_biopsy2_3, features.plot = c("PPP1R15A", "ANXA2"), reduction.use = "umap", cols.use = c("grey", "red", "blue", "green"), overlay = TRUE, no.legend = FALSE, dark.theme = TRUE, pt.size = 0.5)

FeaturePlot(object = PDX_biopsy2_3, features.plot = c("IFI27", "CXCL8"), reduction.use = "umap", cols.use = c("grey", "red", "blue", "green"), overlay = TRUE, no.legend = FALSE, dark.theme = TRUE, pt.size = 0.5)

#feature plot cluster 4, inflamatory response and surface cells
FeaturePlot(object = PDX_biopsy2_3, features.plot = c('CXCL1',	'CXCL2',	'CXCL3',	'CXCL8',	'TNFAIP3',	'PI3',	'CXCL5',	'TNF',	'TNFAIP2',	'IFI27'), reduction.use = "umap", cols.use = c("lightgrey", "blue"), dark.theme = TRUE, pt.size = 0.5)

#violin plots:
#HYPOXIA GENES
VlnPlot(object = PDX_biopsy2_3, features.plot = c("ANXA2", "PPP1R15A"), point.size.use = 0.2, group.by = "orig.ident")
VlnPlot(object = PDX_biopsy2_3, features.plot = c("MTRNR2L1", "MT-RNR2"), point.size.use = 0.2)

VlnPlot(object = PDX_biopsy2_3, features.plot = c("ANXA2", "PPP1R15A"), point.size.use = 0.2, group.by = "orig.ident")


#heatmap
DoHeatmap(object = PDX_biopsy2_3, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
('CXCL1',	'CXCL2',	'CXCL3',	'CXCL8',	'TNFAIP3',	'PI3',	'CXCL5',	'TNF',	'TNFAIP2',	'IFI27')


hypoxia_genes <- read.table("/Users/morsedb/Downloads/geneset (2).txt", header = TRUE)
#remove genes from hypoxia genes that are not in the highly variable genes - didn't seem to remove anything
test = hypoxia_genes[na.omit(match(PDX_biopsy2_3@var.genes, hypoxia_genes)),]
#as list
hypox_vector <- hypoxia_genes[['HALLMARK_HYPOXIA']]

#plot hypoxia gene expression in each layer
DotPlot(object = PDX_biopsy2_3, genes.plot = c("ADM",	"ADORA2B",	"AK4",	"AKAP12",	"ALDOA",	"ALDOB",	"ALDOC",	"AMPD3",	"ANGPTL4",	"ANKZF1",	"ANXA2",	"ATF3",	"ATP7A",	"B3GALT6",	"B4GALNT2",	"BCAN",	"BCL2",	"BGN",	"BHLHE40",	"BNIP3L",	"BRS3",	"BTG1",	"CA12",	"CASP6",	"CAV1",	"CCNG2",	"CDKN1A",	"CDKN1B",	"CDKN1C",	"CHST2",	"CHST3",	"CITED2",	"COL5A1",	"CP",	"CSRP2",	"CTGF",	"CYR61",	"DCN",	"DDIT3",	"DDIT4",	"DPYSL4",	"DTNA",	"DUSP1",	"EDN2",	"EFNA1",	"EFNA3",	"EGFR",	"ENO1",	"ENO2",	"ENO3",	"ERRFI1",	"ETS1",	"EXT1",	"F3",	"FAM162A",	"FBP1",	"FOS",	"FOSL2",	"FOXO3",	"GAA",	"GALK1",	"GAPDH",	"GAPDHS",	"GBE1",	"GCK",	"GCNT2",	"GLRX",	"GPC1",	"GPC3",	"GPC4",	"GPI",	"GRHPR",	"GYS1",	"HAS1",	"HDLBP",	"HEXA",	"HK1",	"HK2",	"HMOX1",	"HOXB9",	"HS3ST1",	"HSPA5",	"IDS",	"IER3",	"IGFBP1",	"IGFBP3",	"IL6",	"ILVBL",	"INHA",	"IRS2",	"ISG20",	"JMJD6",	"JUN",	"KDELR3",	"KDM3A",	"KIF5A",	"KLF6",	"KLF7",	"KLHL24",	"LALBA",	"LDHA",	"LDHC",	"LOX",	"MAFF",	"MAP3K1",	"MIF",	"MT1E",	"MT2A",	"MXI1",	"MYH9",	"NAGK",	"NCAN",	"NDRG1",	"NDST1",	"NDST2",	"NEDD4L",	"NFIL3",	"NR3C1",	"P4HA1",	"P4HA2",	"PAM",	"PCK1",	"PDGFB",	"PDK1",	"PDK3",	"PFKFB3",	"PFKL",	"PFKP",	"PGAM2",	"PGF",	"PGK1",	"PGM1",	"PGM2",	"PHKG1",	"PIM1",	"PKLR",	"PKP1",	"PLAC8",	"PLAUR",	"PLIN2",	"PNRC1",	"PPARGC1A",	"PPFIA4",	"PPP1R15A",	"PPP1R3C",	"PRDX5",	"PRKCA",	"PRKCDBP",	"PTRF",	"PYGM",	"RBPJ",	"RORA",	"RRAGD",	"S100A4",	"SAP30",	"SCARB1",	"SDC2",	"SDC3",	"SDC4",	"SELENBP1",	"SERPINE1",	"SIAH2",	"SLC25A1",	"SLC2A1",	"SLC2A3",	"SLC2A5",	"SLC37A4",	"SLC6A6",	"SRPX",	"STBD1",	"STC1",	"STC2",	"SULT2B1",	"TES",	"TGFB3",	"TGFBI",	"TGM2",	"TIPARP",	"TKTL1",	"TMEM45A",	"TNFAIP3",	"TPBG",	"TPD52",	"TPI1",	"TPST2",	"UGP2",	"VEGFA",	"VHL",	"VLDLR",	"WISP2",	"WSB1",	"XPNPEP1",	"ZFP36",	"ZNF292"), plot.legend = TRUE, group.by = "orig.ident", x.lab.rot = TRUE)


#upregulated genes in cisplatin resistant hypoxia ovarian tumours
#Gap Junction
test <- c("GNAI1",	"GUCY1A3",	"GUCY1B3",	"ITPR3",	"PDGFC",	"PDGFA",	"PRKCA",	"PRKCB",	"TUBB4A")
# Pathways in Cancer
test <- c("FAS",	"JAK1",	"KITLG",	"AR",	"ARNT2",	"CTNNA3",	"FGF1",	"FGF10",	"FGFR2",	"ITGA6",	"JUN",	"PPARG",	"PLD1",	"VEGFC")
# Calcium Signalling
test <- c("ATP2B4",	"CHRNA7",	"CACNA1H",	"CAMK4",	"CYSLTR2",	"GNAL",	"PTGER3",	"P2RX5",	"ERBB3")
#PPAR Signalling
test <- c("CD36",	"ACSL1",	"CPT1A",	"FABP5",	"MMP1",	"SLC27A2")
# Long-term depression
test <- c("PLA2G3")


DotPlot(object = PDX_biopsy2_3, genes.plot = c("MUC16"), plot.legend = TRUE, group.by = "orig.ident", x.lab.rot = TRUE)

#selected hypoxia genes
DotPlot(object = PDX_biopsy2_3, genes.plot = c("ANXA2",	"SDC4",	"CYR61",	"ATF3",	"PPP1R15A",	"NDRG1",	"TNFAIP3",	"SLC2A1",	"MAFF",	"DDIT3",	"S100A4",	"SLC25A1",	"AKAP12",	"TGFB3",	"GAPDH",	"MIF",	"TPI1",	"MYH9",	"FAM162A",	"KLHL24",	"RBPJ",	"PNRC1",	"MXI1",	"LDHA",	"AK4",	"KDM3A",	"UGP2",	"FOXO3",	"KDELR3"), plot.legend = TRUE, group.by = "orig.ident", x.lab.rot = TRUE)






