#testing ground for different heatmaps



#plot heatmap - hallmark genes in clusters: some clusters romoved-----
drops <- c("5_o","2_b", "2_s", "0_o", "4_s", "0_s")

ss_hallmark_s_640 <- as.data.frame(hallmark_s_640)
ss_hallmark_s_640 <- ss_hallmark_s_640[ , !(names(ss_hallmark_s_640) %in% drops)]

ss_hallmark_s_533 <- as.data.frame(hallmark_s_533)
ss_hallmark_s_533 <- ss_hallmark_s_533[ , !(names(ss_hallmark_s_533) %in% drops)]

ss_hallmark_m_124 <- as.data.frame(hallmark_m_124)
ss_hallmark_m_124 <- ss_hallmark_m_124[ , !(names(ss_hallmark_m_124) %in% drops)]

ss_hallmark_m_204 <- as.data.frame(hallmark_m_204)
ss_hallmark_m_204 <- ss_hallmark_m_204[ , !(names(ss_hallmark_m_204) %in% drops)]

ss_hallmark_c_361 <- as.data.frame(hallmark_c_361)
ss_hallmark_c_361 <- ss_hallmark_c_361[ , !(names(ss_hallmark_c_361) %in% drops)]


hm.s640 = Heatmap(as.matrix(ss_hallmark_s_640), name = '640', 
                  column_split = factor(subset_cluster_positions$biology, levels = c("s", "m", "c")),
                  row_km = 1, cluster_column_slices = FALSE, show_row_names = FALSE)
hm.s533 = Heatmap(as.matrix(ss_hallmark_s_533), name = '533', 
                  column_split = factor(subset_cluster_positions$biology, levels = c("s", "m", "c")),
                  row_km = 1, cluster_column_slices = FALSE, show_row_names = FALSE)
hm.m124 = Heatmap(as.matrix(ss_hallmark_m_124), name = '124', 
                  column_split = factor(subset_cluster_positions$biology, levels = c("s", "m", "c")),
                  row_km = 1, cluster_column_slices = FALSE, show_row_names = FALSE)
hm.m204 = Heatmap(as.matrix(ss_hallmark_m_204), name = '204', 
                  column_split = factor(subset_cluster_positions$biology, levels = c("s", "m", "c")),
                  row_km = 1, cluster_column_slices = FALSE, show_row_names = FALSE)
hm.c361 = Heatmap(as.matrix(ss_hallmark_c_361), name = '361', 
                  column_split = factor(subset_cluster_positions$biology, levels = c("s", "m", "c")),
                  row_km = 1, cluster_column_slices = FALSE, show_row_names = FALSE)
hallmark_layers = hm.s640 %v% hm.s533 %v% hm.m124 %v% hm.c361

draw(hallmark_layers)



#now group on average gene sets as a whole and drop some clusters
drops <- c("5_o", "1_o", "5_b", "2_b", "2_s", "0_o", "4_s", "0_s")

ss_scaled_avg_s640 <- as.data.frame(scaled_avg_s640)
ss_scaled_avg_s640 <- ss_scaled_avg_s640[ , !(names(ss_scaled_avg_s640) %in% drops)]
ss_scaled_avg_s640_h.all <- as.data.frame(t(ss_scaled_avg_s640)) %>% select(starts_with("HALL"))
ss_scaled_avg_s640_h.all <- t(ss_scaled_avg_s640_h.all)
rownames(ss_scaled_avg_s640_h.all) = sapply(strsplit(rownames(ss_scaled_avg_s640_h.all),"_"), function(x) paste(x[-1], collapse='_'))

ss_scaled_avg_s533 <- as.data.frame(scaled_avg_s533)
ss_scaled_avg_s533 <- ss_scaled_avg_s533[ , !(names(ss_scaled_avg_s533) %in% drops)]
ss_scaled_avg_s533_h.all <- as.data.frame(t(ss_scaled_avg_s533)) %>% select(starts_with("HALL"))
ss_scaled_avg_s533_h.all <- t(ss_scaled_avg_s533_h.all)
rownames(ss_scaled_avg_s533_h.all) = sapply(strsplit(rownames(ss_scaled_avg_s533_h.all),"_"), function(x) paste(x[-1], collapse='_'))

ss_scaled_avg_m124 <- as.data.frame(scaled_avg_m124)
ss_scaled_avg_m124 <- ss_scaled_avg_m124[ , !(names(ss_scaled_avg_m124) %in% drops)]
ss_scaled_avg_m124_h.all <- as.data.frame(t(ss_scaled_avg_m124)) %>% select(starts_with("HALL"))
ss_scaled_avg_m124_h.all <- t(ss_scaled_avg_m124_h.all)
rownames(ss_scaled_avg_m124_h.all) = sapply(strsplit(rownames(ss_scaled_avg_m124_h.all),"_"), function(x) paste(x[-1], collapse='_'))

ss_scaled_avg_m204 <- as.data.frame(scaled_avg_m204)
ss_scaled_avg_m204 <- ss_scaled_avg_m204[ , !(names(ss_scaled_avg_m204) %in% drops)]
ss_scaled_avg_m204_h.all <- as.data.frame(t(ss_scaled_avg_m204)) %>% select(starts_with("HALL"))
ss_scaled_avg_m204_h.all <- t(ss_scaled_avg_m204_h.all)
rownames(ss_scaled_avg_m204_h.all) = sapply(strsplit(rownames(ss_scaled_avg_m204_h.all),"_"), function(x) paste(x[-1], collapse='_'))

ss_scaled_avg_c361 <- as.data.frame(scaled_avg_c361)
ss_scaled_avg_c361 <- ss_scaled_avg_c361[ , !(names(ss_scaled_avg_c361) %in% drops)]
ss_scaled_avg_c361_h.all <- as.data.frame(t(ss_scaled_avg_c361)) %>% select(starts_with("HALL"))
ss_scaled_avg_c361_h.all <- t(ss_scaled_avg_c361_h.all)
rownames(ss_scaled_avg_c361_h.all) = sapply(strsplit(rownames(ss_scaled_avg_c361_h.all),"_"), function(x) paste(x[-1], collapse='_'))

subset2_cluster_positions <- data.frame(row.names = colnames(ss_scaled_avg_s640))
subset2_cluster_positions$biology <- c("m", "c", "s", "s", "s", "c", "m", "s", "s", "c", "s", "m")




h1 <- Heatmap(as.matrix(ss_scaled_avg_s640_h.all), row_names_side = "right",
              column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
              column_names_side = "top",col = cols.use, show_column_names = T,
              show_row_names = T,
              cluster_rows = T, cluster_columns = T, cluster_column_slices = FALSE, row_names_max_width = unit(30, 'cm'),
              row_names_gp = gpar(fontsize = c(9)),
              heatmap_legend_param = list(legend_width = unit(4, "cm"), title='avg. scaled expression', direction = "horizontal"))

h2 <- Heatmap(as.matrix(ss_scaled_avg_s533_h.all), row_names_side = "right",
              column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
              column_names_side = "top",col = cols.use, show_column_names = T,
              show_row_names = T,
              cluster_rows = T, cluster_columns = T, cluster_column_slices = FALSE, row_names_max_width = unit(30, 'cm'),
              row_names_gp = gpar(fontsize = c(9)), show_heatmap_legend = FALSE)
              #heatmap_legend_param = list(legend_height = unit(9, "cm"), title='Row Z-Score\nAvg. expression'))

h3 <- Heatmap(as.matrix(ss_scaled_avg_m124_h.all), row_names_side = "right",
              column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
              column_names_side = "top",col = cols.use, show_column_names = T,
              show_row_names = T,
              cluster_rows = T, cluster_columns = T, cluster_column_slices = FALSE, row_names_max_width = unit(30, 'cm'),
              row_names_gp = gpar(fontsize = c(9)), show_heatmap_legend = FALSE)

h4 <- Heatmap(as.matrix(ss_scaled_avg_m204_h.all), row_names_side = "right",
              column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
              column_names_side = "top",col = cols.use, show_column_names = T,
              show_row_names = T,
              cluster_rows = T, cluster_columns = T, cluster_column_slices = FALSE, row_names_max_width = unit(30, 'cm'),
              row_names_gp = gpar(fontsize = c(9)), show_heatmap_legend = FALSE)

h5 <- Heatmap(as.matrix(ss_scaled_avg_c361_h.all), row_names_side = "right",
              column_split = factor(subset2_cluster_positions$biology, levels = c("s", "m", "c")),
              column_names_side = "top",col = cols.use, show_column_names = T,
              show_row_names = T,
              cluster_rows = T, cluster_columns = T, cluster_column_slices = FALSE, row_names_max_width = unit(30, 'cm'),
              row_names_gp = gpar(fontsize = c(9)), show_heatmap_legend = FALSE)

draw(h1 %v% h2 %v% h3 %v% h5, column_title = "average expression per hallmark gene set", heatmap_legend_side = "bottom")















#aleksandra's code:
heat = Heatmap(unigen, cluster_columns = TRUE)

tab = annot[,grep('HALLMARK',names(annot))]#[,1:2]
colo = vector('list',ncol(tab))
names(colo) = colnames(tab)
colo = lapply(colo, function(x) c("FALSE"='white', "TRUE"='grey') )
tab = HeatmapAnnotation(tab, which='row', show_legend = FALSE, show_annotation_name = TRUE, col = colo, annotation_name_side='top',annotation_name_rot=90 )

hlist=heat+tab
draw(hlist, padding = unit(c(0, 0, 3, 0), "cm"))






annot = readRDS("/Volumes/ncats/NCATS Matrix Screening/single-cell data/layer_comparison/GSEA_across_layers/surface_640_annotation.rds")
unigen = readRDS("/Volumes/ncats/NCATS Matrix Screening/single-cell data/layer_comparison/GSEA_across_layers/all_genes_unique_surface640.rds")
surf = readRDS("/Volumes/ncats/NCATS Matrix Screening/single-cell data/layer_comparison/GSEA_across_layers/surface640.rds")

library(ComplexHeatmap)

heat = Heatmap(unigen, cluster_columns = TRUE, row_names_gp = gpar(fontsize = 6), column_title = "HALLMARK")

tab = annot[,grep('HALLMARK',names(annot))]#[,1:2]
colnames(tab) = sapply(strsplit(colnames(tab),"_"), function(x) paste(x[-1], collapse='_'))
colo = vector('list',ncol(tab))
names(colo) = colnames(tab)
colo = lapply(colo, function(x) c("FALSE"='white', "TRUE"='grey') )
colo = list(tab = c("FALSE"="white", "TRUE"="grey"))
tab_annot = HeatmapAnnotation(test = as.matrix(tab), which='row', show_legend = FALSE, show_annotation_name = TRUE, col = colo, annotation_name_side='top',annotation_name_rot=90, annotation_name_gp=gpar(fontsize=5), width=unit(3,'cm'))

hlist=heat+tab_annot
draw(hlist, padding = unit(c(0, 0, 3, 0), "cm"))






