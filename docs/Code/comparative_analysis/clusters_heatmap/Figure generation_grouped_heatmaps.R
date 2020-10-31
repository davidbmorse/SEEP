#heatmap figure generation
setwd('/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/comparative_analysis_filt_biop/clusters_heatmap')
library(data.table)
library(ComplexHeatmap)
library(scales)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
col.biology <- list("biology" = c("c" = "blue", "m" = "green", "s" = "red"))
cols.use <- colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100) # reversed RdBu, creates blue/white/red



#1 average expression across hallmark gene sets ---------
load(file = "comparison_data/gene_sets_avg_exp/avg_gene_set_matrixes.rda")
load(file = "comparison_data/gene_sets_avg_exp/annotations/gene_set_annotations.rda")

# set up matrix for avg_gene_set plotting----
scaled_avg_s644 <- avg_s644[,-c('category','subcategory')]
scaled_avg_s644 <- t(scale(t(scaled_avg_s644)))
scaled_avg_s644 <- rescale(scaled_avg_s644, to=c(-2,2))
rownames(scaled_avg_s644) <- avg_s644$subcategory

scaled_avg_s535 <- avg_s535[,-c('category','subcategory')]
scaled_avg_s535 <- t(scale(t(scaled_avg_s535)))
scaled_avg_s535 <- rescale(scaled_avg_s535, to=c(-2,2))
rownames(scaled_avg_s535) <- avg_s535$subcategory

scaled_avg_m422 <- avg_m422[,-c('category','subcategory')]
scaled_avg_m422 <- t(scale(t(scaled_avg_m422)))
scaled_avg_m422 <- rescale(scaled_avg_m422, to=c(-2,2))
rownames(scaled_avg_m422) <- avg_m422$subcategory

scaled_avg_c363 <- avg_c363[,-c('category','subcategory')]
scaled_avg_c363 <- t(scale(t(scaled_avg_c363)))
scaled_avg_c363 <- rescale(scaled_avg_c363, to=c(-2,2))
rownames(scaled_avg_c363) <- avg_c363$subcategory

#now group on average gene sets as a whole and drop some clusters
drops <- c("0_s", "2_s", "1_s", "0_o", "1_o", "5_o", "0_b", "1_b")

#subset into hallmark sets ONLY and remove certain clusters ----
ss_scaled_avg_s644 <- as.data.frame(scaled_avg_s644)
ss_scaled_avg_s644 <- ss_scaled_avg_s644[ , !(names(ss_scaled_avg_s644) %in% drops)]
ss_scaled_avg_s644_h.all <- as.data.frame(t(ss_scaled_avg_s644)) %>% select(starts_with("HALL"))
ss_scaled_avg_s644_h.all <- t(ss_scaled_avg_s644_h.all)
rownames(ss_scaled_avg_s644_h.all) = sapply(strsplit(rownames(ss_scaled_avg_s644_h.all),"_"), function(x) paste(x[-1], collapse='_'))

ss_scaled_avg_s535 <- as.data.frame(scaled_avg_s535)
ss_scaled_avg_s535 <- ss_scaled_avg_s535[ , !(names(ss_scaled_avg_s535) %in% drops)]
ss_scaled_avg_s535_h.all <- as.data.frame(t(ss_scaled_avg_s535)) %>% select(starts_with("HALL"))
ss_scaled_avg_s535_h.all <- t(ss_scaled_avg_s535_h.all)
rownames(ss_scaled_avg_s535_h.all) = sapply(strsplit(rownames(ss_scaled_avg_s535_h.all),"_"), function(x) paste(x[-1], collapse='_'))

ss_scaled_avg_m422 <- as.data.frame(scaled_avg_m422)
ss_scaled_avg_m422 <- ss_scaled_avg_m422[ , !(names(ss_scaled_avg_m422) %in% drops)]
ss_scaled_avg_m422_h.all <- as.data.frame(t(ss_scaled_avg_m422)) %>% select(starts_with("HALL"))
ss_scaled_avg_m422_h.all <- t(ss_scaled_avg_m422_h.all)
rownames(ss_scaled_avg_m422_h.all) = sapply(strsplit(rownames(ss_scaled_avg_m422_h.all),"_"), function(x) paste(x[-1], collapse='_'))

ss_scaled_avg_c363 <- as.data.frame(scaled_avg_c363)
ss_scaled_avg_c363 <- ss_scaled_avg_c363[ , !(names(ss_scaled_avg_c363) %in% drops)]
ss_scaled_avg_c363_h.all <- as.data.frame(t(ss_scaled_avg_c363)) %>% select(starts_with("HALL"))
ss_scaled_avg_c363_h.all <- t(ss_scaled_avg_c363_h.all)
rownames(ss_scaled_avg_c363_h.all) = sapply(strsplit(rownames(ss_scaled_avg_c363_h.all),"_"), function(x) paste(x[-1], collapse='_'))

subset2_cluster_positions <- data.frame(row.names = colnames(ss_scaled_avg_s644))
subset2_cluster_positions$biology <- c("c", "m", "s", "s", "s", "c", "m", "s", "m", "c", "s", "s")


#plot heatmap ------
h1 <- Heatmap(as.matrix(ss_scaled_avg_s644_h.all), row_names_side = "right",
              column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
              column_names_side = "top",col = cols.use, show_column_names = T,
              show_row_names = T,
              cluster_rows = T, cluster_columns = T, cluster_column_slices = FALSE, row_names_max_width = unit(30, 'cm'),
              row_names_gp = gpar(fontsize = c(9)),
              heatmap_legend_param = list(legend_width = unit(4, "cm"), title='avg. scaled expression', direction = "horizontal"))

h2 <- Heatmap(as.matrix(ss_scaled_avg_s535_h.all), row_names_side = "right",
              column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
              column_names_side = "top",col = cols.use, show_column_names = T,
              show_row_names = T,
              cluster_rows = T, cluster_columns = T, cluster_column_slices = FALSE, row_names_max_width = unit(30, 'cm'),
              row_names_gp = gpar(fontsize = c(9)), show_heatmap_legend = FALSE)
#heatmap_legend_param = list(legend_height = unit(9, "cm"), title='Row Z-Score\nAvg. expression'))

h3 <- Heatmap(as.matrix(ss_scaled_avg_m422_h.all), row_names_side = "right",
              column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
              column_names_side = "top",col = cols.use, show_column_names = T,
              show_row_names = T,
              cluster_rows = T, cluster_columns = T, cluster_column_slices = FALSE, row_names_max_width = unit(30, 'cm'),
              row_names_gp = gpar(fontsize = c(9)), show_heatmap_legend = FALSE)

h5 <- Heatmap(as.matrix(ss_scaled_avg_c363_h.all), row_names_side = "right",
              column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
              column_names_side = "top",col = cols.use, show_column_names = T,
              show_row_names = T,
              cluster_rows = T, cluster_columns = T, cluster_column_slices = FALSE, row_names_max_width = unit(30, 'cm'),
              row_names_gp = gpar(fontsize = c(9)), show_heatmap_legend = FALSE)

draw(h1 %v% h2 %v% h3 %v% h5, column_title = "average expression per hallmark gene set", heatmap_legend_side = "bottom")





#drop middle clusters (and gene sets) 4/26/20 #########
drops_middle <- c("0_s", "1_s", "2_s", "4_s", "0_o", "1_o", "2_o", "5_o", "0_b", "1_b", "2_b")

#subset into hallmark sets ONLY and remove certain clusters ----
ss_scaled_avg_s644 <- as.data.frame(scaled_avg_s644)
ss_scaled_avg_s644 <- ss_scaled_avg_s644[ , !(names(ss_scaled_avg_s644) %in% drops_middle)]
ss_scaled_avg_s644_h.all <- as.data.frame(t(ss_scaled_avg_s644)) %>% select(starts_with("HALL"))
ss_scaled_avg_s644_h.all <- t(ss_scaled_avg_s644_h.all)
rownames(ss_scaled_avg_s644_h.all) = sapply(strsplit(rownames(ss_scaled_avg_s644_h.all),"_"), function(x) paste(x[-1], collapse='_'))

ss_scaled_avg_s535 <- as.data.frame(scaled_avg_s535)
ss_scaled_avg_s535 <- ss_scaled_avg_s535[ , !(names(ss_scaled_avg_s535) %in% drops_middle)]
ss_scaled_avg_s535_h.all <- as.data.frame(t(ss_scaled_avg_s535)) %>% select(starts_with("HALL"))
ss_scaled_avg_s535_h.all <- t(ss_scaled_avg_s535_h.all)
rownames(ss_scaled_avg_s535_h.all) = sapply(strsplit(rownames(ss_scaled_avg_s535_h.all),"_"), function(x) paste(x[-1], collapse='_'))

ss_scaled_avg_c363 <- as.data.frame(scaled_avg_c363)
ss_scaled_avg_c363 <- ss_scaled_avg_c363[ , !(names(ss_scaled_avg_c363) %in% drops_middle)]
ss_scaled_avg_c363_h.all <- as.data.frame(t(ss_scaled_avg_c363)) %>% select(starts_with("HALL"))
ss_scaled_avg_c363_h.all <- t(ss_scaled_avg_c363_h.all)
rownames(ss_scaled_avg_c363_h.all) = sapply(strsplit(rownames(ss_scaled_avg_c363_h.all),"_"), function(x) paste(x[-1], collapse='_'))

subset2_cluster_positions <- data.frame(row.names = colnames(ss_scaled_avg_s644))
subset2_cluster_positions$biology <- c("c", "s", "s", "s", "c", "s", "c", "s", "s")

#plot heatmap ------
h1 <- Heatmap(as.matrix(ss_scaled_avg_s644_h.all), row_names_side = "right",
              column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
              column_names_side = "top",col = cols.use, show_column_names = T,
              show_row_names = T,
              cluster_rows = T, cluster_columns = T, cluster_column_slices = FALSE, row_names_max_width = unit(30, 'cm'),
              row_names_gp = gpar(fontsize = c(9)),
              heatmap_legend_param = list(legend_width = unit(4, "cm"), title='avg. scaled expression', direction = "horizontal"))

h2 <- Heatmap(as.matrix(ss_scaled_avg_s535_h.all), row_names_side = "right",
              column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
              column_names_side = "top",col = cols.use, show_column_names = T,
              show_row_names = T,
              cluster_rows = T, cluster_columns = T, cluster_column_slices = FALSE, row_names_max_width = unit(30, 'cm'),
              row_names_gp = gpar(fontsize = c(9)), show_heatmap_legend = FALSE)

h3 <- Heatmap(as.matrix(ss_scaled_avg_c363_h.all), row_names_side = "right",
              column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
              column_names_side = "top",col = cols.use, show_column_names = T,
              show_row_names = T,
              cluster_rows = T, cluster_columns = T, cluster_column_slices = FALSE, row_names_max_width = unit(30, 'cm'),
              row_names_gp = gpar(fontsize = c(9)), show_heatmap_legend = FALSE)

draw(h1 %v% h2 %v% h3, column_title = "average expression per hallmark gene set", heatmap_legend_side = "bottom")

pdf("exports/avg_hallmrk_per_set_SvC.pdf", width = 4.2, height = 5)
draw(h1 %v% h2 %v% h3, column_title = "average expression per hallmark gene set", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()








#try to subset to remove certain row (repeated gene sets)------
subset2_avg_s644_h.all <- ss_scaled_avg_s644_h.all
subset2_avg_s535_h.all <- ss_scaled_avg_s535_h.all
subset2_avg_m422_h.all <- ss_scaled_avg_m422_h.all[-c(4), ]
subset2_avg_c363_h.all <- ss_scaled_avg_c363_h.all

subset2_avg_s644_h.all <- rescale(subset2_avg_s644_h.all, to=c(-2,2))
subset2_avg_s535_h.all <- rescale(subset2_avg_s535_h.all, to=c(-2,2))
subset2_avg_m422_h.all <- rescale(subset2_avg_m422_h.all, to=c(-2,2))
subset2_avg_c363_h.all <- rescale(subset2_avg_c363_h.all, to=c(-2,2))


#plot subset heatmap
h1.s <- Heatmap(as.matrix(subset2_avg_s644_h.all), row_names_side = "right",
              column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
              column_names_side = "top",col = cols.use, show_column_names = T,
              show_row_names = T, row_dend_width = unit(4, "mm"), column_dend_height = unit(5, "mm"),
              cluster_rows = T, cluster_columns = T, cluster_column_slices = FALSE, row_names_max_width = unit(30, 'cm'),
              row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)),
              heatmap_legend_param = list(legend_width = unit(3, "cm"), grid_height = unit(3, "mm"), title_position = "topcenter",
                                          labels_gp = gpar(fontsize = 8), title='avg. scaled expression', direction = "horizontal"))

h2.s <- Heatmap(as.matrix(subset2_avg_s535_h.all), row_names_side = "right",
              column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
              column_names_side = "top",col = cols.use, show_column_names = T,
              show_row_names = T, row_dend_width = unit(4, "mm"), column_dend_height = unit(5, "mm"),
              cluster_rows = T, cluster_columns = T, cluster_column_slices = FALSE, row_names_max_width = unit(30, 'cm'),
              row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)), show_heatmap_legend = FALSE)

h3.s <- Heatmap(as.matrix(subset2_avg_m422_h.all), row_names_side = "right",
              column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
              column_names_side = "top",col = cols.use, show_column_names = T,
              show_row_names = T, row_dend_width = unit(4, "mm"), column_dend_height = unit(5, "mm"),
              cluster_rows = T, cluster_columns = T, cluster_column_slices = FALSE, row_names_max_width = unit(30, 'cm'),
              row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)), show_heatmap_legend = FALSE)

h5.s <- Heatmap(as.matrix(subset2_avg_c363_h.all), row_names_side = "right",
              column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
              column_names_side = "top",col = cols.use, show_column_names = T,
              show_row_names = T, row_dend_width = unit(4, "mm"), column_dend_height = unit(5, "mm"),
              cluster_rows = T, cluster_columns = T, cluster_column_slices = FALSE, row_names_max_width = unit(30, 'cm'),
              row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)), show_heatmap_legend = FALSE)

draw(h1.s %v% h2.s %v% h3.s %v% h5.s, column_title = "average expression per hallmark gene set", heatmap_legend_side = "bottom")

pdf("exports/avg_hallmrk_per_set2.pdf", width = 4, height = 5.3)
draw(h1.s %v% h2.s %v% h3.s %v% h5.s, column_title = "average expression per hallmark gene set", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()








#1 Individual gene expression for specific hallmark gene sets ---------
#set column ordering
cluster_order <- c("6_s", "4_b", "4_o", "3_o", "5_s", "5_b", "2_b", "2_o", "4_s", "3_s", "6_o", "3_b")
#new
cluster_order <- c("6_s", "4_b", "4_o", "3_o", "5_s", "5_b", "3_s", "6_o", "3_b")
#load
s_644_matrices <- readRDS("comparison_data/surface644.rds")
s_535_matrices <- readRDS("comparison_data/surface535.rds")
m_422_matrices <- readRDS("comparison_data/middle422.rds")
c_363_matrices <- readRDS("comparison_data/center363.rds")
cols.use <- colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100) # reversed RdBu, creates Blue-white-red
cols.use2 <- colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdGy")))(100) # reversed RdBu, creates Blue-white-red
drops <- c("0_s", "1_s", "2_s", "0_o", "1_o", "5_o", "0_b", "1_b")
#new
drops <- c("0_s", "1_s", "2_s", "4_s", "0_o", "1_o", "2_o", "5_o", "0_b", "1_b", "2_b")



#s644 #####
#INFy
s644_INFy <- s_644_matrices$h.all$HALLMARK_INTERFERON_GAMMA_RESPONSE
scld_s644_INFy <- t(scale(t(s644_INFy)))
scld_s644_INFy <- rescale(scld_s644_INFy, to=c(-2,2))

#subset into hallmark sets ONLY and remove certain clusters ----
scld_s644_INFy <- as.data.frame(scld_s644_INFy)
scld_s644_INFy <- scld_s644_INFy[ , !(names(scld_s644_INFy) %in% drops)]
#sum columns and plot as lower annotation to display gene set (LE) density
scld_s644_INFy_ColSums <- as.data.frame(colSums(scld_s644_INFy))
colnames(scld_s644_INFy_ColSums) <- c("sums")

scld_s644_INFy_ColSums <- as.vector(colSums(scld_s644_INFy))

subset2_cluster_positions <- data.frame(row.names = colnames(scld_s644_INFy))
subset2_cluster_positions$biology <- c("c", "s", "s",  "s", "c", "s", "c", "s", "s")

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_s644_INFy), baseline = 0), height = unit(1, "cm"))
HM_s644_INFy <- Heatmap(as.matrix(scld_s644_INFy), row_names_side = "right",
                column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                column_names_side = "top",col = cols.use2, show_column_names = F,
                cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                #row_dend_gp = gpar(unit(10, "cm")),
                row_dend_width = unit(2, "mm"),
                heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                column_gap = unit(0, "mm"), border = TRUE,
                bottom_annotation = hbplot, column_order = cluster_order,
                row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_s644_INFy, column_title = "INFy gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

#####
pdf("plot_HM_s644_INFy.pdf", width = 2.2, height = 4.8)
draw(HM_s644_INFy, column_title = "INFy gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()


#INFa
s644_INFa <- s_644_matrices$h.all$HALLMARK_INTERFERON_ALPHA_RESPONSE
scld_s644_INFa <- t(scale(t(s644_INFa)))
scld_s644_INFa <- rescale(scld_s644_INFa, to=c(-2,2))

scld_s644_INFa <- as.data.frame(scld_s644_INFa)
scld_s644_INFa <- scld_s644_INFa[ , !(names(scld_s644_INFa) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_s644_INFa), baseline = 0), height = unit(1, "cm"))
HM_s644_INFa <- Heatmap(as.matrix(scld_s644_INFa), row_names_side = "right",
                        column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                        column_names_side = "top",col = cols.use2, show_column_names = F,
                        cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                        #row_dend_gp = gpar(unit(10, "cm")),
                        row_dend_width = unit(2, "mm"),
                        heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                    title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                        column_gap = unit(0, "mm"), border = TRUE,
                        bottom_annotation = hbplot, column_order = cluster_order,
                        row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_s644_INFa, column_title = "INFa gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("plot_HM_s644_INFa.pdf", width = 2.2, height = 3.6)
draw(HM_s644_INFa, column_title = "INFa gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()

#TNFa_via_NFKB
s644_TNFa_NFkb <- s_644_matrices$h.all$HALLMARK_TNFA_SIGNALING_VIA_NFKB
scld_s644_TNFa_NFkb <- t(scale(t(s644_TNFa_NFkb)))
scld_s644_TNFa_NFkb <- rescale(scld_s644_TNFa_NFkb, to=c(-2,2))

scld_s644_TNFa_NFkb <- as.data.frame(scld_s644_TNFa_NFkb)
scld_s644_TNFa_NFkb <- scld_s644_TNFa_NFkb[ , !(names(scld_s644_TNFa_NFkb) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_s644_TNFa_NFkb), baseline = 0), height = unit(1, "cm"))
HM_s644_TNFa_NFkb <- Heatmap(as.matrix(scld_s644_TNFa_NFkb), row_names_side = "right",
                        column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                        column_names_side = "top",col = cols.use2, show_column_names = F,
                        cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                        #row_dend_gp = gpar(unit(10, "cm")),
                        row_dend_width = unit(2, "mm"),
                        heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                    title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                        column_gap = unit(0, "mm"), border = TRUE,
                        bottom_annotation = hbplot, column_order = cluster_order,
                        row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_s644_TNFa_NFkb, column_title = "TNFa_NFkb gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("plot_HM_s644_TNFa_NFkb.pdf", width = 2.2, height = 5.5)
draw(HM_s644_TNFa_NFkb, column_title = "TNFa_NFkb gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()


#IL6_JAK_STAT
s644_IL6_JAK <- s_644_matrices$h.all$HALLMARK_IL6_JAK_STAT3_SIGNALING
scld_s644_IL6_JAK <- t(scale(t(s644_IL6_JAK)))
scld_s644_IL6_JAK <- rescale(scld_s644_IL6_JAK, to=c(-2,2))

scld_s644_IL6_JAK <- as.data.frame(scld_s644_IL6_JAK)
scld_s644_IL6_JAK <- scld_s644_IL6_JAK[ , !(names(scld_s644_IL6_JAK) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_s644_IL6_JAK), baseline = 0), height = unit(1, "cm"))
HM_s644_IL6_JAK <- Heatmap(as.matrix(scld_s644_IL6_JAK), row_names_side = "right",
                             column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                             column_names_side = "top",col = cols.use2, show_column_names = F,
                           cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                           #row_dend_gp = gpar(unit(10, "cm")),
                           row_dend_width = unit(2, "mm"),
                           heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                       title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                           column_gap = unit(0, "mm"), border = TRUE,
                           bottom_annotation = hbplot, column_order = cluster_order,
                             row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_s644_IL6_JAK, column_title = "IL6_JAK gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("plot_HM_s644_IL6_JAK.pdf", width = 2.2, height = 2.6)
draw(HM_s644_IL6_JAK, column_title = "IL6_JAK gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()


#INFLAM
s644_INFLAM <- s_644_matrices$h.all$HALLMARK_INFLAMMATORY_RESPONSE
scld_s644_INFLAM <- t(scale(t(s644_INFLAM)))
scld_s644_INFLAM <- rescale(scld_s644_INFLAM, to=c(-2,2))

scld_s644_INFLAM <- as.data.frame(scld_s644_INFLAM)
scld_s644_INFLAM <- scld_s644_INFLAM[ , !(names(scld_s644_INFLAM) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_s644_INFLAM), baseline = 0), height = unit(1, "cm"))
HM_s644_INFLAM <- Heatmap(as.matrix(scld_s644_INFLAM), row_names_side = "right",
                           column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                           column_names_side = "top",col = cols.use2, show_column_names = F,
                          cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                          #row_dend_gp = gpar(unit(10, "cm")),
                          row_dend_width = unit(2, "mm"),
                          heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                      title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                          column_gap = unit(0, "mm"), border = TRUE,
                          bottom_annotation = hbplot, column_order = cluster_order,
                           row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_s644_INFLAM, column_title = "INFLAM gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("plot_HM_s644_INFLAM.pdf", width = 2.2, height = 4)
draw(HM_s644_INFLAM, column_title = "INFLAM gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()

#APOPTOSIS
s644_APOPTOSIS <- s_644_matrices$h.all$HALLMARK_APOPTOSIS
scld_s644_APOPTOSIS <- t(scale(t(s644_APOPTOSIS)))
scld_s644_APOPTOSIS <- rescale(scld_s644_APOPTOSIS, to=c(-2,2))

scld_s644_APOPTOSIS <- as.data.frame(scld_s644_APOPTOSIS)
scld_s644_APOPTOSIS <- scld_s644_APOPTOSIS[ , !(names(scld_s644_APOPTOSIS) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_s644_APOPTOSIS), baseline = 0), height = unit(1, "cm"))
HM_s644_APOPTOSIS <- Heatmap(as.matrix(scld_s644_APOPTOSIS), row_names_side = "right",
                          column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                          column_names_side = "top",col = cols.use2, show_column_names = F,
                          cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                          #row_dend_gp = gpar(unit(10, "cm")),
                          row_dend_width = unit(2, "mm"),
                          heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                      title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                          column_gap = unit(0, "mm"), border = TRUE,
                          bottom_annotation = hbplot, column_order = cluster_order,
                          row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_s644_APOPTOSIS, column_title = "APOPTOSIS gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("plot_HM_s644_APOPTOSIS.pdf", width = 2.2, height = 3.4)
draw(HM_s644_APOPTOSIS, column_title = "APOPTOSIS gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()


#ALLOGRAFT
s644_ALLOGRAFT <- s_644_matrices$h.all$HALLMARK_ALLOGRAFT_REJECTION
scld_s644_ALLOGRAFT <- t(scale(t(s644_ALLOGRAFT)))
scld_s644_ALLOGRAFT <- rescale(scld_s644_ALLOGRAFT, to=c(-2,2))

scld_s644_ALLOGRAFT <- as.data.frame(scld_s644_ALLOGRAFT)
scld_s644_ALLOGRAFT <- scld_s644_ALLOGRAFT[ , !(names(scld_s644_ALLOGRAFT) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_s644_ALLOGRAFT), baseline = 0), height = unit(1, "cm"))
HM_s644_ALLOGRAFT <- Heatmap(as.matrix(scld_s644_ALLOGRAFT), row_names_side = "right",
                             column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                             column_names_side = "top",col = cols.use2, show_column_names = F,
                             cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                             #row_dend_gp = gpar(unit(10, "cm")),
                             row_dend_width = unit(2, "mm"),
                             heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                         title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                             column_gap = unit(0, "mm"), border = TRUE,
                             bottom_annotation = hbplot, column_order = cluster_order,
                             row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_s644_ALLOGRAFT, column_title = "ALLOGRAFT gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("plot_HM_s644_ALLOGRAFT.pdf", width = 2.2, height = 4)
draw(HM_s644_ALLOGRAFT, column_title = "ALLOGRAFT gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()


#KRAS
s644_KRAS <- s_644_matrices$h.all$HALLMARK_KRAS_SIGNALING_UP
scld_s644_KRAS <- t(scale(t(s644_KRAS)))
scld_s644_KRAS <- rescale(scld_s644_KRAS, to=c(-2,2))

scld_s644_KRAS <- as.data.frame(scld_s644_KRAS)
scld_s644_KRAS <- scld_s644_KRAS[ , !(names(scld_s644_KRAS) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_s644_KRAS), baseline = 0), height = unit(1, "cm"))
HM_s644_KRAS <- Heatmap(as.matrix(scld_s644_KRAS), row_names_side = "right",
                        column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                        column_names_side = "top",col = cols.use2, show_column_names = F,
                        cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                        #row_dend_gp = gpar(unit(10, "cm")),
                        row_dend_width = unit(2, "mm"),
                        heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                    title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                        column_gap = unit(0, "mm"), border = TRUE,
                        bottom_annotation = hbplot, column_order = cluster_order,
                        row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_s644_KRAS, column_title = "KRAS gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("plot_HM_s644_KRAS.pdf", width = 2.2, height = 3.15)
draw(HM_s644_KRAS, column_title = "KRAS gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()


#P53
s644_P53 <- s_644_matrices$h.all$HALLMARK_P53_PATHWAY
scld_s644_P53 <- t(scale(t(s644_P53)))
scld_s644_P53 <- rescale(scld_s644_P53, to=c(-2,2))

scld_s644_P53 <- as.data.frame(scld_s644_P53)
scld_s644_P53 <- scld_s644_P53[ , !(names(scld_s644_P53) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_s644_P53), baseline = 0), height = unit(1, "cm"))
HM_s644_P53 <- Heatmap(as.matrix(scld_s644_P53), row_names_side = "right",
                        column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                        column_names_side = "top",col = cols.use2, show_column_names = F,
                       cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                       #row_dend_gp = gpar(unit(10, "cm")),
                       row_dend_width = unit(2, "mm"),
                       heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                   title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                       column_gap = unit(0, "mm"), border = TRUE,
                       bottom_annotation = hbplot, column_order = cluster_order,
                        row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_s644_P53, column_title = "P53 gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("plot_HM_s644_P53.pdf", width = 2.2, height = 3.2)
draw(HM_s644_P53, column_title = "P53 gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()


#COMPLEMENT
s644_COMPLEMENT <- s_644_matrices$h.all$HALLMARK_COMPLEMENT
scld_s644_COMPLEMENT <- t(scale(t(s644_COMPLEMENT)))
scld_s644_COMPLEMENT <- rescale(scld_s644_COMPLEMENT, to=c(-2,2))

scld_s644_COMPLEMENT <- as.data.frame(scld_s644_COMPLEMENT)
scld_s644_COMPLEMENT <- scld_s644_COMPLEMENT[ , !(names(scld_s644_COMPLEMENT) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_s644_COMPLEMENT), baseline = 0), height = unit(1, "cm"))
HM_s644_COMPLEMENT <- Heatmap(as.matrix(scld_s644_COMPLEMENT), row_names_side = "right",
                       column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                       column_names_side = "top",col = cols.use2, show_column_names = F,
                       cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                       #row_dend_gp = gpar(unit(10, "cm")),
                       row_dend_width = unit(2, "mm"),
                       heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                   title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                       column_gap = unit(0, "mm"), border = TRUE,
                       bottom_annotation = hbplot, column_order = cluster_order,
                       row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_s644_COMPLEMENT, column_title = "COMPLEMENT gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("plot_HM_s644_COMPLEMENT.pdf", width = 2.2, height = 3.4)
draw(HM_s644_COMPLEMENT, column_title = "COMPLEMENT gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()

#COAGULATION
s644_COAGULATION <- s_644_matrices$h.all$HALLMARK_COAGULATION
scld_s644_COAGULATION <- t(scale(t(s644_COAGULATION)))
scld_s644_COAGULATION <- rescale(scld_s644_COAGULATION, to=c(-2,2))

scld_s644_COAGULATION <- as.data.frame(scld_s644_COAGULATION)
scld_s644_COAGULATION <- scld_s644_COAGULATION[ , !(names(scld_s644_COAGULATION) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_s644_COAGULATION), baseline = 0), height = unit(1, "cm"))
HM_s644_COAGULATION <- Heatmap(as.matrix(scld_s644_COAGULATION), row_names_side = "right",
                       column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                       column_names_side = "top",col = cols.use2, show_column_names = F,
                       cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                       #row_dend_gp = gpar(unit(10, "cm")),
                       row_dend_width = unit(2, "mm"),
                       heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                   title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                       column_gap = unit(0, "mm"), border = TRUE,
                       bottom_annotation = hbplot, column_order = cluster_order,
                       row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_s644_COAGULATION, column_title = "COAGULATION gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("plot_HM_s644_COAGULATION.pdf", width = 2.2, height = 3.2)
draw(HM_s644_COAGULATION, column_title = "COAGULATION gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()


#EPITHELIAL_MESENCHYMAL_TRANSITION
s644_EPITHELIAL_MESENCHYMAL_TRANSITION <- s_644_matrices$h.all$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION
scld_s644_EPITHELIAL_MESENCHYMAL_TRANSITION <- t(scale(t(s644_EPITHELIAL_MESENCHYMAL_TRANSITION)))
scld_s644_EPITHELIAL_MESENCHYMAL_TRANSITION <- rescale(scld_s644_EPITHELIAL_MESENCHYMAL_TRANSITION, to=c(-2,2))

scld_s644_EPITHELIAL_MESENCHYMAL_TRANSITION <- as.data.frame(scld_s644_EPITHELIAL_MESENCHYMAL_TRANSITION)
scld_s644_EPITHELIAL_MESENCHYMAL_TRANSITION <- scld_s644_EPITHELIAL_MESENCHYMAL_TRANSITION[ , !(names(scld_s644_EPITHELIAL_MESENCHYMAL_TRANSITION) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_s644_EPITHELIAL_MESENCHYMAL_TRANSITION), baseline = 0), height = unit(1, "cm"))
HM_s644_EPITHELIAL_MESENCHYMAL_TRANSITION <- Heatmap(as.matrix(scld_s644_EPITHELIAL_MESENCHYMAL_TRANSITION), row_names_side = "right",
                       column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                       column_names_side = "top",col = cols.use2, show_column_names = F,
                       cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                       #row_dend_gp = gpar(unit(10, "cm")),
                       row_dend_width = unit(2, "mm"),
                       heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                   title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                       column_gap = unit(0, "mm"), border = TRUE,
                       bottom_annotation = hbplot, column_order = cluster_order,
                       row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_s644_EPITHELIAL_MESENCHYMAL_TRANSITION, column_title = "EPITHELIAL_MESENCHYMAL_TRANSITION gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("plot_HM_s644_EPITHELIAL_MESENCHYMAL_TRANSITION.pdf", width = 2.2, height = 2.7)
draw(HM_s644_EPITHELIAL_MESENCHYMAL_TRANSITION, column_title = "EPITHELIAL_MESENCHYMAL_TRANSITION gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()



#s535 #####
#MYC1
s535_MYC1 <- s_535_matrices$h.all$HALLMARK_MYC_TARGETS_V1
scld_s535_MYC1 <- t(scale(t(s535_MYC1)))
scld_s535_MYC1 <- rescale(scld_s535_MYC1, to=c(-2,2))

scld_s535_MYC1 <- as.data.frame(scld_s535_MYC1)
scld_s535_MYC1 <- scld_s535_MYC1[ , !(names(scld_s535_MYC1) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_s535_MYC1), baseline = 0), height = unit(1, "cm"))
HM_s535_MYC1 <- Heatmap(as.matrix(scld_s535_MYC1), row_names_side = "right",
                       column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                       column_names_side = "top",col = cols.use2, show_column_names = F,
                       cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                       #row_dend_gp = gpar(unit(10, "cm")),
                       row_dend_width = unit(2, "mm"),
                       heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                   title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                       column_gap = unit(0, "mm"), border = TRUE,
                       bottom_annotation = hbplot, column_order = cluster_order,
                       row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_s535_MYC1, column_title = "MYC1 gene expression\nfor surface535", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("exports/s535_sets/plot_HM_s535_MYC1.pdf", width = 2.2, height = 4.8)
draw(HM_s535_MYC1, column_title = "MYC1 gene expression\nfor surface535", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()


#MYC2
s535_MYC2 <- s_535_matrices$h.all$HALLMARK_MYC_TARGETS_V2
scld_s535_MYC2 <- t(scale(t(s535_MYC2)))
scld_s535_MYC2 <- rescale(scld_s535_MYC2, to=c(-2,2))

scld_s535_MYC2 <- as.data.frame(scld_s535_MYC2)
scld_s535_MYC2 <- scld_s535_MYC2[ , !(names(scld_s535_MYC2) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_s535_MYC2), baseline = 0), height = unit(1, "cm"))
HM_s535_MYC2 <- Heatmap(as.matrix(scld_s535_MYC2), row_names_side = "right",
                        column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                        column_names_side = "top",col = cols.use2, show_column_names = F,
                        cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                        #row_dend_gp = gpar(unit(10, "cm")),
                        row_dend_width = unit(2, "mm"),
                        heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                    title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                        column_gap = unit(0, "mm"), border = TRUE,
                        bottom_annotation = hbplot, column_order = cluster_order,
                        row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_s535_MYC2, column_title = "MYC2 gene expression\nfor surface535", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("exports/s535_sets/plot_HM_s535_MYC2.pdf", width = 2.2, height = 2.35)
draw(HM_s535_MYC2, column_title = "MYC2 gene expression\nfor surface535", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()

#DNA_REPAIR
s535_DNA_REPAIR <- s_535_matrices$h.all$HALLMARK_DNA_REPAIR
scld_s535_DNA_REPAIR <- t(scale(t(s535_DNA_REPAIR)))
scld_s535_DNA_REPAIR <- rescale(scld_s535_DNA_REPAIR, to=c(-2,2))

scld_s535_DNA_REPAIR <- as.data.frame(scld_s535_DNA_REPAIR)
scld_s535_DNA_REPAIR <- scld_s535_DNA_REPAIR[ , !(names(scld_s535_DNA_REPAIR) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_s535_DNA_REPAIR), baseline = 0), height = unit(1, "cm"))
HM_s535_DNA_REPAIR <- Heatmap(as.matrix(scld_s535_DNA_REPAIR), row_names_side = "right",
                        column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                        column_names_side = "top",col = cols.use2, show_column_names = F,
                        cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                        #row_dend_gp = gpar(unit(10, "cm")),
                        row_dend_width = unit(2, "mm"),
                        heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                    title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                        column_gap = unit(0, "mm"), border = TRUE,
                        bottom_annotation = hbplot, column_order = cluster_order,
                        row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_s535_DNA_REPAIR, column_title = "DNA_REPAIR gene expression\nfor surface535", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("exports/s535_sets/plot_HM_s535_DNA_REPAIR.pdf", width = 2.2, height = 3.5)
draw(HM_s535_DNA_REPAIR, column_title = "DNA_REPAIR gene expression\nfor surface535", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()

#E2F
s535_E2F_TARGETS <- s_535_matrices$h.all$HALLMARK_E2F_TARGETS
scld_s535_E2F_TARGETS <- t(scale(t(s535_E2F_TARGETS)))
scld_s535_E2F_TARGETS <- rescale(scld_s535_E2F_TARGETS, to=c(-2,2))

scld_s535_E2F_TARGETS <- as.data.frame(scld_s535_E2F_TARGETS)
scld_s535_E2F_TARGETS <- scld_s535_E2F_TARGETS[ , !(names(scld_s535_E2F_TARGETS) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_s535_E2F_TARGETS), baseline = 0), height = unit(1, "cm"))
HM_s535_E2F_TARGETS <- Heatmap(as.matrix(scld_s535_E2F_TARGETS), row_names_side = "right",
                              column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                              column_names_side = "top",col = cols.use2, show_column_names = F,
                              cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                              #row_dend_gp = gpar(unit(10, "cm")),
                              row_dend_width = unit(2, "mm"),
                              heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                          title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                              column_gap = unit(0, "mm"), border = TRUE,
                              bottom_annotation = hbplot, column_order = cluster_order,
                              row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_s535_E2F_TARGETS, column_title = "E2F_TARGETS gene expression\nfor surface535", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("exports/s535_sets/plot_HM_s535_E2F_TARGETS.pdf", width = 2.2, height = 6.5)
draw(HM_s535_E2F_TARGETS, column_title = "E2F_TARGETS gene expression\nfor surface535", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()


#G2M
s535_G2M_CHECKPOINT <- s_535_matrices$h.all$HALLMARK_G2M_CHECKPOINT
scld_s535_G2M_CHECKPOINT <- t(scale(t(s535_G2M_CHECKPOINT)))
scld_s535_G2M_CHECKPOINT <- rescale(scld_s535_G2M_CHECKPOINT, to=c(-2,2))

scld_s535_G2M_CHECKPOINT <- as.data.frame(scld_s535_G2M_CHECKPOINT)
scld_s535_G2M_CHECKPOINT <- scld_s535_G2M_CHECKPOINT[ , !(names(scld_s535_G2M_CHECKPOINT) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_s535_G2M_CHECKPOINT), baseline = 0), height = unit(1, "cm"))
HM_s535_G2M_CHECKPOINT <- Heatmap(as.matrix(scld_s535_G2M_CHECKPOINT), row_names_side = "right",
                               column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                               column_names_side = "top",col = cols.use2, show_column_names = F,
                               cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                               #row_dend_gp = gpar(unit(10, "cm")),
                               row_dend_width = unit(2, "mm"),
                               heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                           title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                               column_gap = unit(0, "mm"), border = TRUE,
                               bottom_annotation = hbplot, column_order = cluster_order,
                               row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_s535_G2M_CHECKPOINT, column_title = "G2M_CHECKPOINT gene expression\nfor surface535", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("exports/s535_sets/plot_HM_s535_G2M_CHECKPOINT.pdf", width = 2.2, height = 7)
draw(HM_s535_G2M_CHECKPOINT, column_title = "G2M_CHECKPOINT gene expression\nfor surface535", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()



#m422 #####
#INFy
m422_INFy <- m_422_matrices$h.all$HALLMARK_INTERFERON_ALPHA_RESPONSE
scld_m422_INFy <- t(scale(t(m422_INFy)))
scld_m422_INFy <- rescale(scld_m422_INFy, to=c(-2,2))

scld_m422_INFy <- as.data.frame(scld_m422_INFy)
scld_m422_INFy <- scld_m422_INFy[ , !(names(scld_m422_INFy) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_m422_INFy), baseline = 0), height = unit(1, "cm"))
HM_m422_INFy <- Heatmap(as.matrix(scld_m422_INFy), row_names_side = "right",
                        column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                        column_names_side = "top",col = cols.use2, show_column_names = F,
                        cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                        #row_dend_gp = gpar(unit(10, "cm")),
                        row_dend_width = unit(2, "mm"),
                        heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                    title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                        column_gap = unit(0, "mm"), border = TRUE,
                        bottom_annotation = hbplot, column_order = cluster_order,
                        row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_m422_INFy, column_title = "INFy gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("exports/m422_sets/plot_HM_m422_INFy.pdf", width = 2.2, height = 2.3)
draw(HM_m422_INFy, column_title = "INFy gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()

#PROT_SEC
m422_PROT_SEC <- m_422_matrices$h.all$HALLMARK_PROTEIN_SECRETION
scld_m422_PROT_SEC <- t(scale(t(m422_PROT_SEC)))
scld_m422_PROT_SEC <- rescale(scld_m422_PROT_SEC, to=c(-2,2))

scld_m422_PROT_SEC <- as.data.frame(scld_m422_PROT_SEC)
scld_m422_PROT_SEC <- scld_m422_PROT_SEC[ , !(names(scld_m422_PROT_SEC) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_m422_PROT_SEC), baseline = 0), height = unit(1, "cm"))
HM_m422_PROT_SEC <- Heatmap(as.matrix(scld_m422_PROT_SEC), row_names_side = "right",
                        column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                        column_names_side = "top",col = cols.use2, show_column_names = F,
                        cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                        #row_dend_gp = gpar(unit(10, "cm")),
                        row_dend_width = unit(2, "mm"),
                        heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                    title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                        column_gap = unit(0, "mm"), border = TRUE,
                        bottom_annotation = hbplot, column_order = cluster_order,
                        row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_m422_PROT_SEC, column_title = "PROT_SEC gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("exports/m422_sets/plot_HM_m422_PROT_SEC.pdf", width = 2.2, height = 3.3)
draw(HM_m422_PROT_SEC, column_title = "PROT_SEC gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()

#COMPLEMENT
m422_COMPLEMENT <- m_422_matrices$h.all$HALLMARK_COMPLEMENT
scld_m422_COMPLEMENT <- t(scale(t(m422_COMPLEMENT)))
scld_m422_COMPLEMENT <- rescale(scld_m422_COMPLEMENT, to=c(-2,2))

scld_m422_COMPLEMENT <- as.data.frame(scld_m422_COMPLEMENT)
scld_m422_COMPLEMENT <- scld_m422_COMPLEMENT[ , !(names(scld_m422_COMPLEMENT) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_m422_COMPLEMENT), baseline = 0), height = unit(1, "cm"))
HM_m422_COMPLEMENT <- Heatmap(as.matrix(scld_m422_COMPLEMENT), row_names_side = "right",
                            column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                            column_names_side = "top",col = cols.use2, show_column_names = F,
                            cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                            heatmap_legend_param = list(legend_height = unit(3, "cm"), title='scaled expression'),
                            column_gap = unit(0, "mm"), border = TRUE,
                            bottom_annotation = hbplot, column_order = cluster_order,
                            row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_m422_COMPLEMENT, column_title = "COMPLEMENT gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("exports/m422_sets/plot_HM_m422_COMPLEMENT.pdf", width = 2.2, height = 2.8)
draw(HM_m422_COMPLEMENT, column_title = "COMPLEMENT gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()

#HYPOXIA
m422_HYPOXIA <- m_422_matrices$h.all$HALLMARK_HYPOXIA
scld_m422_HYPOXIA <- t(scale(t(m422_HYPOXIA)))
scld_m422_HYPOXIA <- rescale(scld_m422_HYPOXIA, to=c(-2,2))

scld_m422_HYPOXIA <- as.data.frame(scld_m422_HYPOXIA)
scld_m422_HYPOXIA <- scld_m422_HYPOXIA[ , !(names(scld_m422_HYPOXIA) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_m422_HYPOXIA), baseline = 0), height = unit(1, "cm"))
HM_m422_HYPOXIA <- Heatmap(as.matrix(scld_m422_HYPOXIA), row_names_side = "right",
                              column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                              column_names_side = "top",col = cols.use2, show_column_names = F,
                           cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                           #row_dend_gp = gpar(unit(10, "cm")),
                           row_dend_width = unit(2, "mm"),
                           heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                       title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                           column_gap = unit(0, "mm"), border = TRUE,
                           bottom_annotation = hbplot, column_order = cluster_order,
                              row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_m422_HYPOXIA, column_title = "HYPOXIA gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("exports/m422_sets/plot_HM_m422_HYPOXIA.pdf", width = 2.2, height = 2.85)
draw(HM_m422_HYPOXIA, column_title = "HYPOXIA gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()


#HEME_MET
m422_HEME_MET <- m_422_matrices$h.all$HALLMARK_HEME_METABOLISM
scld_m422_HEME_MET <- t(scale(t(m422_HEME_MET)))
scld_m422_HEME_MET <- rescale(scld_m422_HEME_MET, to=c(-2,2))

scld_m422_HEME_MET <- as.data.frame(scld_m422_HEME_MET)
scld_m422_HEME_MET <- scld_m422_HEME_MET[ , !(names(scld_m422_HEME_MET) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_m422_HEME_MET), baseline = 0), height = unit(1, "cm"))
HM_m422_HEME_MET <- Heatmap(as.matrix(scld_m422_HEME_MET), row_names_side = "right",
                           column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                           column_names_side = "top",col = cols.use2, show_column_names = F,
                           cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                           #row_dend_gp = gpar(unit(10, "cm")),
                           row_dend_width = unit(2, "mm"),
                           heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                       title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                           column_gap = unit(0, "mm"), border = TRUE,
                           bottom_annotation = hbplot, column_order = cluster_order,
                           row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_m422_HEME_MET, column_title = "HEME_MET gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("exports/m422_sets/plot_HM_m422_HEME_MET.pdf", width = 2.2, height = 3.2)
draw(HM_m422_HEME_MET, column_title = "HEME_MET gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()


#TNFa_via_NFKB
m422_TNFa_via_NFKB <- m_422_matrices$h.all$HALLMARK_TNFA_SIGNALING_VIA_NFKB
scld_m422_TNFa_via_NFKB <- t(scale(t(m422_TNFa_via_NFKB)))
scld_m422_TNFa_via_NFKB <- rescale(scld_m422_TNFa_via_NFKB, to=c(-2,2))

scld_m422_TNFa_via_NFKB <- as.data.frame(scld_m422_TNFa_via_NFKB)
scld_m422_TNFa_via_NFKB <- scld_m422_TNFa_via_NFKB[ , !(names(scld_m422_TNFa_via_NFKB) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_m422_TNFa_via_NFKB), baseline = 0), height = unit(1, "cm"))
HM_m422_TNFa_via_NFKB <- Heatmap(as.matrix(scld_m422_TNFa_via_NFKB), row_names_side = "right",
                            column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                            column_names_side = "top",col = cols.use2, show_column_names = F,
                            cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                            heatmap_legend_param = list(legend_height = unit(3, "cm"), title='scaled expression'),
                            column_gap = unit(0, "mm"), border = TRUE,
                            bottom_annotation = hbplot, column_order = cluster_order,
                            row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_m422_TNFa_via_NFKB, column_title = "TNFa_via_NFKB gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("exports/m422_sets/plot_HM_m422_TNFa_via_NFKB.pdf", width = 2.2, height = 4.8)
draw(HM_m422_TNFa_via_NFKB, column_title = "TNFa_via_NFKB gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()



#c363 #####
#UV_DN
c363_UV_DN <- c_363_matrices$h.all$HALLMARK_UV_RESPONSE_DN
scld_c363_UV_DN <- t(scale(t(c363_UV_DN)))
scld_c363_UV_DN <- rescale(scld_c363_UV_DN, to=c(-2,2))

scld_c363_UV_DN <- as.data.frame(scld_c363_UV_DN)
scld_c363_UV_DN <- scld_c363_UV_DN[ , !(names(scld_c363_UV_DN) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_c363_UV_DN), baseline = 0), height = unit(1, "cm"))
HM_c363_UV_DN <- Heatmap(as.matrix(scld_c363_UV_DN), row_names_side = "right",
                                 column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                                 column_names_side = "top",col = cols.use2, show_column_names = F,
                         cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                         #row_dend_gp = gpar(unit(10, "cm")),
                         row_dend_width = unit(2, "mm"),
                         heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                     title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                         column_gap = unit(0, "mm"), border = TRUE,
                         bottom_annotation = hbplot, column_order = cluster_order,
                                 row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_c363_UV_DN, column_title = "UV_DN gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("exports/c363_sets/plot_HM_c363_UV_DN.pdf", width = 2.2, height = 3.5)
draw(HM_c363_UV_DN, column_title = "UV_DN gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()


#HEME_MET
c363_HEME_MET <- c_363_matrices$h.all$HALLMARK_HEME_METABOLISM
scld_c363_HEME_MET <- t(scale(t(c363_HEME_MET)))
scld_c363_HEME_MET <- rescale(scld_c363_HEME_MET, to=c(-2,2))

scld_c363_HEME_MET <- as.data.frame(scld_c363_HEME_MET)
scld_c363_HEME_MET <- scld_c363_HEME_MET[ , !(names(scld_c363_HEME_MET) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_c363_HEME_MET), baseline = 0), height = unit(1, "cm"))
HM_c363_HEME_MET <- Heatmap(as.matrix(scld_c363_HEME_MET), row_names_side = "right",
                         column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                         column_names_side = "top",col = cols.use2, show_column_names = F,
                         cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                         #row_dend_gp = gpar(unit(10, "cm")),
                         row_dend_width = unit(2, "mm"),
                         heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                     title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                         column_gap = unit(0, "mm"), border = TRUE,
                         bottom_annotation = hbplot, column_order = cluster_order,
                         row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_c363_HEME_MET, column_title = "HEME_MET gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("exports/c363_sets/plot_HM_c363_HEME_MET.pdf", width = 2.2, height = 3.3)
draw(HM_c363_HEME_MET, column_title = "HEME_MET gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()


#MYOGENESIS
c363_MYOGENESIS <- c_363_matrices$h.all$HALLMARK_MYOGENESIS
scld_c363_MYOGENESIS <- t(scale(t(c363_MYOGENESIS)))
scld_c363_MYOGENESIS <- rescale(scld_c363_MYOGENESIS, to=c(-2,2))

scld_c363_MYOGENESIS <- as.data.frame(scld_c363_MYOGENESIS)
scld_c363_MYOGENESIS <- scld_c363_MYOGENESIS[ , !(names(scld_c363_MYOGENESIS) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_c363_MYOGENESIS), baseline = 0), height = unit(1, "cm"))
HM_c363_MYOGENESIS <- Heatmap(as.matrix(scld_c363_MYOGENESIS), row_names_side = "right",
                            column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                            column_names_side = "top",col = cols.use2, show_column_names = F,
                            cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                            #row_dend_gp = gpar(unit(10, "cm")),
                            row_dend_width = unit(2, "mm"),
                            heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                        title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                            column_gap = unit(0, "mm"), border = TRUE,
                            bottom_annotation = hbplot, column_order = cluster_order,
                            row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_c363_MYOGENESIS, column_title = "MYOGENESIS gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("exports/c363_sets/plot_HM_c363_MYOGENESIS.pdf", width = 2.2, height = 2.35)
draw(HM_c363_MYOGENESIS, column_title = "MYOGENESIS gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()


#PROTEIN_SECRETION
c363_PROTEIN_SECRETION <- c_363_matrices$h.all$HALLMARK_PROTEIN_SECRETION
scld_c363_PROTEIN_SECRETION <- t(scale(t(c363_PROTEIN_SECRETION)))
scld_c363_PROTEIN_SECRETION <- rescale(scld_c363_PROTEIN_SECRETION, to=c(-2,2))

scld_c363_PROTEIN_SECRETION <- as.data.frame(scld_c363_PROTEIN_SECRETION)
scld_c363_PROTEIN_SECRETION <- scld_c363_PROTEIN_SECRETION[ , !(names(scld_c363_PROTEIN_SECRETION) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_c363_PROTEIN_SECRETION), baseline = 0), height = unit(1, "cm"))
HM_c363_PROTEIN_SECRETION <- Heatmap(as.matrix(scld_c363_PROTEIN_SECRETION), row_names_side = "right",
                              column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                              column_names_side = "top",col = cols.use2, show_column_names = F,
                              cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                              #row_dend_gp = gpar(unit(10, "cm")),
                              row_dend_width = unit(2, "mm"),
                              heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                          title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                              column_gap = unit(0, "mm"), border = TRUE,
                              bottom_annotation = hbplot, column_order = cluster_order,
                              row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_c363_PROTEIN_SECRETION, column_title = "PROTEIN_SECRETION gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("exports/c363_sets/plot_HM_c363_PROTEIN_SECRETION.pdf", width = 2.2, height = 3.2)
draw(HM_c363_PROTEIN_SECRETION, column_title = "PROTEIN_SECRETION gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()


#PI3K_AKT_MTOR_SIGNALING
c363_PI3K_AKT_MTOR_SIGNALING <- c_363_matrices$h.all$HALLMARK_PI3K_AKT_MTOR_SIGNALING
scld_c363_PI3K_AKT_MTOR_SIGNALING <- t(scale(t(c363_PI3K_AKT_MTOR_SIGNALING)))
scld_c363_PI3K_AKT_MTOR_SIGNALING <- rescale(scld_c363_PI3K_AKT_MTOR_SIGNALING, to=c(-2,2))

scld_c363_PI3K_AKT_MTOR_SIGNALING <- as.data.frame(scld_c363_PI3K_AKT_MTOR_SIGNALING)
scld_c363_PI3K_AKT_MTOR_SIGNALING <- scld_c363_PI3K_AKT_MTOR_SIGNALING[ , !(names(scld_c363_PI3K_AKT_MTOR_SIGNALING) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_c363_PI3K_AKT_MTOR_SIGNALING), baseline = 0), height = unit(1, "cm"))
HM_c363_PI3K_AKT_MTOR_SIGNALING <- Heatmap(as.matrix(scld_c363_PI3K_AKT_MTOR_SIGNALING), row_names_side = "right",
                              column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                              column_names_side = "top",col = cols.use2, show_column_names = F,
                              cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                              #row_dend_gp = gpar(unit(10, "cm")),
                              row_dend_width = unit(2, "mm"),
                              heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                          title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                              column_gap = unit(0, "mm"), border = TRUE,
                              bottom_annotation = hbplot, column_order = cluster_order,
                              row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_c363_PI3K_AKT_MTOR_SIGNALING, column_title = "PI3K_AKT_MTOR_SIGNALING gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("exports/c363_sets/plot_HM_c363_PI3K_AKT_MTOR_SIGNALING.pdf", width = 2.2, height = 3.8)
draw(HM_c363_PI3K_AKT_MTOR_SIGNALING, column_title = "PI3K_AKT_MTOR_SIGNALING gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()


#UNFOLDED_PROTEIN_RESPONSE
c363_UNFOLDED_PROTEIN_RESPONSE <- c_363_matrices$h.all$HALLMARK_UNFOLDED_PROTEIN_RESPONSE
scld_c363_UNFOLDED_PROTEIN_RESPONSE <- t(scale(t(c363_UNFOLDED_PROTEIN_RESPONSE)))
scld_c363_UNFOLDED_PROTEIN_RESPONSE <- rescale(scld_c363_UNFOLDED_PROTEIN_RESPONSE, to=c(-2,2))

scld_c363_UNFOLDED_PROTEIN_RESPONSE <- as.data.frame(scld_c363_UNFOLDED_PROTEIN_RESPONSE)
scld_c363_UNFOLDED_PROTEIN_RESPONSE <- scld_c363_UNFOLDED_PROTEIN_RESPONSE[ , !(names(scld_c363_UNFOLDED_PROTEIN_RESPONSE) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_c363_UNFOLDED_PROTEIN_RESPONSE), baseline = 0), height = unit(1, "cm"))
HM_c363_UNFOLDED_PROTEIN_RESPONSE <- Heatmap(as.matrix(scld_c363_UNFOLDED_PROTEIN_RESPONSE), row_names_side = "right",
                              column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                              column_names_side = "top",col = cols.use2, show_column_names = F,
                              cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                              #row_dend_gp = gpar(unit(10, "cm")),
                              row_dend_width = unit(2, "mm"),
                              heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                          title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                              column_gap = unit(0, "mm"), border = TRUE,
                              bottom_annotation = hbplot, column_order = cluster_order,
                              row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_c363_UNFOLDED_PROTEIN_RESPONSE, column_title = "UNFOLDED_PROTEIN_RESPONSE gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("exports/c363_sets/plot_HM_c363_UNFOLDED_PROTEIN_RESPONSE.pdf", width = 2.2, height = 2.6)
draw(HM_c363_UNFOLDED_PROTEIN_RESPONSE, column_title = "UNFOLDED_PROTEIN_RESPONSE gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()

#MTORC1_SIGNALING
c363_MTORC1_SIGNALING <- c_363_matrices$h.all$HALLMARK_MTORC1_SIGNALING
scld_c363_MTORC1_SIGNALING <- t(scale(t(c363_MTORC1_SIGNALING)))
scld_c363_MTORC1_SIGNALING <- rescale(scld_c363_MTORC1_SIGNALING, to=c(-2,2))

scld_c363_MTORC1_SIGNALING <- as.data.frame(scld_c363_MTORC1_SIGNALING)
scld_c363_MTORC1_SIGNALING <- scld_c363_MTORC1_SIGNALING[ , !(names(scld_c363_MTORC1_SIGNALING) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_c363_MTORC1_SIGNALING), baseline = 0), height = unit(1, "cm"))
HM_c363_MTORC1_SIGNALING <- Heatmap(as.matrix(scld_c363_MTORC1_SIGNALING), row_names_side = "right",
                              column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                              column_names_side = "top",col = cols.use2, show_column_names = F,
                              cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                              #row_dend_gp = gpar(unit(10, "cm")),
                              row_dend_width = unit(2, "mm"),
                              heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                          title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                              column_gap = unit(0, "mm"), border = TRUE,
                              bottom_annotation = hbplot, column_order = cluster_order,
                              row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_c363_MTORC1_SIGNALING, column_title = "MTORC1_SIGNALING gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("exports/c363_sets/plot_HM_c363_MTORC1_SIGNALING.pdf", width = 2.2, height = 4)
draw(HM_c363_MTORC1_SIGNALING, column_title = "MTORC1_SIGNALING gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()


#ANDROGEN_RESPONSE
c363_ANDROGEN_RESPONSE <- c_363_matrices$h.all$HALLMARK_ANDROGEN_RESPONSE
scld_c363_ANDROGEN_RESPONSE <- t(scale(t(c363_ANDROGEN_RESPONSE)))
scld_c363_ANDROGEN_RESPONSE <- rescale(scld_c363_ANDROGEN_RESPONSE, to=c(-2,2))

scld_c363_ANDROGEN_RESPONSE <- as.data.frame(scld_c363_ANDROGEN_RESPONSE)
scld_c363_ANDROGEN_RESPONSE <- scld_c363_ANDROGEN_RESPONSE[ , !(names(scld_c363_ANDROGEN_RESPONSE) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_c363_ANDROGEN_RESPONSE), baseline = 0), height = unit(1, "cm"))
HM_c363_ANDROGEN_RESPONSE <- Heatmap(as.matrix(scld_c363_ANDROGEN_RESPONSE), row_names_side = "right",
                              column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                              column_names_side = "top",col = cols.use2, show_column_names = F,
                              cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                              #row_dend_gp = gpar(unit(10, "cm")),
                              row_dend_width = unit(2, "mm"),
                              heatmap_legend_param = list(legend_height = unit(3, "cm"), grid_height = unit(3, "mm"), 
                                                          title='scaled expression', direction = "horizontal", title_position = "topcenter"),
                              column_gap = unit(0, "mm"), border = TRUE,
                              bottom_annotation = hbplot, column_order = cluster_order,
                              row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_c363_ANDROGEN_RESPONSE, column_title = "ANDROGEN_RESPONSE gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("exports/c363_sets/plot_HM_c363_ANDROGEN_RESPONSE.pdf", width = 2.2, height = 2.35)
draw(HM_c363_ANDROGEN_RESPONSE, column_title = "ANDROGEN_RESPONSE gene expression\nfor middle422", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()





GO_INFLAMMATORY_RESPONSE
GO_CYTOKINE_MEDIATED_SIGNALING_PATHWAY
GO_NEGATIVE_REGULATION_OF_IMMUNE_RESPONSE
GO_REGULATION_OF_CYTOKINE_BIOSYNTHETIC_PROCESS


ZHANG_RESPONSE_TO_IKK_INHIBITOR_AND_TNF_UP
SANA_TNF_SIGNALING_UP
SEKI_INFLAMMATORY_RESPONSE_LPS_UP
HINATA_NFKB_TARGETS_KERATINOCYTE_UP
PHONG_TNF_TARGETS_UP
ZWANG_CLASS_3_TRANSIENTLY_INDUCED_BY_EGF
BROWNE_INTERFERON_RESPONSIVE_GENES
BASSO_CD40_SIGNALING_UP










#testing grounds
s644_IL6_JAK <- s_644_matrices$h.all$HALLMARK_COMPLEMENT
scld_s644_IL6_JAK <- t(scale(t(s644_IL6_JAK)))
scld_s644_IL6_JAK <- rescale(scld_s644_IL6_JAK, to=c(-2,2))

scld_s644_IL6_JAK <- as.data.frame(scld_s644_IL6_JAK)
scld_s644_IL6_JAK <- scld_s644_IL6_JAK[ , !(names(scld_s644_IL6_JAK) %in% drops)]

hbplot = HeatmapAnnotation(sum = anno_barplot(colSums(scld_s644_IL6_JAK), baseline = 0), height = unit(1, "cm"))
HM_s644_IL6_JAK <- Heatmap(as.matrix(scld_s644_IL6_JAK), row_names_side = "right",
                           column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
                           column_names_side = "top",col = cols.use2, show_column_names = F,
                           cluster_rows = T, cluster_columns = F,  cluster_column_slices = FALSE,
                           heatmap_legend_param = list(legend_height = unit(3, "cm"), title='scaled expression'),
                           column_gap = unit(0, "mm"), border = TRUE,
                           bottom_annotation = hbplot, column_order = cluster_order,
                           row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)))
draw(HM_s644_IL6_JAK, column_title = "IL6_JAK gene expression\nfor surface644", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")





















