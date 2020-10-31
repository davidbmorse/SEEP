#grouped biology/layer heatmaps

setwd('/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/comparative_analysis_filt_biop/clusters_heatmap')
library(data.table)
library(ComplexHeatmap)
library(scales)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
cols.use <- colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100) # reversed RdBu, creates blue/white/red
cols.use3 <- colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdGy")))(100) # reversed RdBu, creates blue/white/red

#load data/annotations

s_644_matrices <- readRDS("comparison_data/surface644.rds")
s_535_matrices <- readRDS("comparison_data/surface535.rds")
m_422_matrices <- readRDS("comparison_data/middle422.rds")
c_363_matrices <- readRDS("comparison_data/center363.rds")

surface644 <- readRDS('comparison_data/all_genes_unique_surface644.rds')
surface535 <- readRDS('comparison_data/all_genes_unique_surface535.rds')
middle422 <- readRDS('comparison_data/all_genes_unique_middle422.rds')
center363 <- readRDS('comparison_data/all_genes_unique_center363.rds')

surface644_annotation <- readRDS('comparison_data/surface_644_annotation.rds')
surface535_annotation <- readRDS('comparison_data/surface_535_annotation.rds')
middle422_annotation <- readRDS('comparison_data/middle_422_annotation.rds')
center363_annotation <- readRDS('comparison_data/center_363_annotation.rds')

#extract only hallmark genes
hallmark_s_644 <- unique(do.call(rbind, s_644_matrices[[1]]))
hallmark_s_535 <- unique(do.call(rbind, s_535_matrices[[1]]))
hallmark_m_422 <- unique(do.call(rbind, m_422_matrices[[1]]))
hallmark_c_363 <- unique(do.call(rbind, c_363_matrices[[1]]))

#define cluster positions
cluster_positions <- data.frame(row.names = colnames(surface644))
cluster_positions$biology <- c("m", "m", "m", "c", "m", "s", "s", "c", "s", "m", "s", "c", "m", "s", "c", "s", "m", "c", "s", "s")
col.biology <- list("biology" = c("c" = "blue", "m" = "green", "s" = "red"))
cols.use <- colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100) # reversed RdBu, creates blue/white/red


#plot heatmap - hallmark genes in clusters-----
hm.s644 = Heatmap(hallmark_s_644, name = '644', 
                   column_split = factor(cluster_positions$biology, levels = c("s", "m", "c")),
                   row_km = 1, cluster_column_slices = FALSE)
hm.s535 = Heatmap(hallmark_s_535, name = '535', 
                  column_split = factor(cluster_positions$biology, levels = c("s", "m", "c")),
                  row_km = 1, cluster_column_slices = FALSE)
hm.m422 = Heatmap(hallmark_m_422, name = '422', 
                  column_split = factor(cluster_positions$biology, levels = c("s", "m", "c")),
                  row_km = 1, cluster_column_slices = FALSE)
hm.c363 = Heatmap(hallmark_c_363, name = '363', 
                  column_split = factor(cluster_positions$biology, levels = c("s", "m", "c")),
                  row_km = 1, cluster_column_slices = FALSE)
hallmark_layers = hm.s644 %v% hm.s535 %v% hm.m422 %v% hm.c363

draw(hallmark_layers)


#plot heatmap - hallmark genes in groups with cluster annotation-----
hm.s644 = Heatmap(hallmark_s_644, name = '644', 
                  #column_split = factor(cluster_positions$biology, levels = c("s", "m", "c")),
                  row_km = 1, cluster_column_slices = FALSE)
hm.s535 = Heatmap(hallmark_s_535, name = '535', 
                  #column_split = factor(cluster_positions$biology, levels = c("s", "m", "c")),
                  row_km = 1, cluster_column_slices = FALSE)
hm.m422 = Heatmap(hallmark_m_422, name = '422', 
                  #column_split = factor(cluster_positions$biology, levels = c("s", "m", "c")),
                  row_km = 1, cluster_column_slices = FALSE)
hm.c363 = Heatmap(hallmark_c_363, name = '363', 
                  #column_split = factor(cluster_positions$biology, levels = c("s", "m", "c")),
                  row_km = 1, cluster_column_slices = FALSE)
cluster_labels = bottom_annotation = HeatmapAnnotation(df = cluster_positions, col = col.biology)
hallmark_layers2 = hm.s644 %v% hm.s535 %v% hm.m422 %v% hm.c363 %v% cluster_labels

draw(hallmark_layers2, column_km = 1)

#full heatmap - genes in clusters-----
drops <- NULL
drops <- c("0_s", "1_s", "2_s", "0_o", "1_o", "5_o", "0_b", "1_b")
drops_no_middle <- c("0_s", "1_s", "2_s", "4_s", "0_o", "1_o", "2_o", "5_o", "0_b", "1_b", "2_b")

surface644 <- as.data.frame(surface644)
surface644 <- surface644[ , !(names(surface644) %in% drops_no_middle)]
surface535 <- as.data.frame(surface535)
surface535 <- surface535[ , !(names(surface535) %in% drops_no_middle)]
middle422 <- as.data.frame(middle422)
middle422 <- middle422[ , !(names(middle422) %in% drops_no_middle)]
center363 <- as.data.frame(center363)
center363 <- center363[ , !(names(center363) %in% drops_no_middle)]
#rescale if you want
surface644 <- t(scale(t(surface644)))
surface644 <- rescale(surface644, to=c(-2,2))
surface535 <- t(scale(t(surface535)))
surface535 <- rescale(surface535, to=c(-2,2))
middle422 <- t(scale(t(middle422)))
middle422 <- rescale(middle422, to=c(-2,2))
center363 <- t(scale(t(center363)))
center363 <- rescale(center363, to=c(-2,2))

column_subset <- data.frame(row.names = colnames(surface644))
#column_subset$biology <- c("c", "m", "s", "s", "s", "c", "m", "s", "m", "c", "s", "s")
#for no middle
column_subset$biology <- c("c", "s", "s", "s", "c", "s", "c", "s", "s")

#label my favorite genes
to_label_s644 <- which(rownames(surface644) %in% c("IFI27", "TNF", "HLA-E", "NFKBIA", "SAT1", "IL1B"))
to_label_s535 <- which(rownames(surface535) %in% c("RFC4", "DUT", "COX17", "SNRPG", "POLE4", "PTEN", "E2F1", "KIF2C"))
to_label_m422 <- which(rownames(middle422) %in% c("VEGFA", "JUNB", "KLF4", "ID1"))
to_label_c363 <- which(rownames(center363) %in% c("PTEN", "FOXJ1", "TOP1", "AK4", "MAPK1", "YWHAB", "NCK1", "HSPA4", "VAMP4", "XPOT", "LIFR", "GNAS", "BET1"))

r_644 = rowAnnotation(test = anno_mark(at = to_label_s644, labels = rownames(surface644)[to_label_s644]))
r_535 = rowAnnotation(test = anno_mark(at = to_label_s535, labels = rownames(surface535)[to_label_s535]))
r_422 = rowAnnotation(test = anno_mark(at = to_label_m422, labels = rownames(middle422)[to_label_m422]))
r_363 = rowAnnotation(test = anno_mark(at = to_label_c363, labels = rownames(center363)[to_label_c363]))

s644 = Heatmap(as.matrix(surface644), name = '644', 
               column_split = factor(column_subset$biology, levels = c("s", "m", "c")), column_names_side = "top",
               show_heatmap_legend = T, show_row_names = F, cluster_column_slices = F,
               right_annotation = r_644, row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)), col = cols.use3,
               heatmap_legend_param = list(legend_width = unit(2, "cm"), title='avg. scaled expression', direction = "horizontal"))
s535 = Heatmap(as.matrix(surface535), name = '535', 
               column_split = factor(column_subset$biology, levels = c("s", "m", "c")), column_names_side = "top",
               show_heatmap_legend = F, show_row_names = F, cluster_column_slices = F,
               right_annotation = r_535, row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)), col = cols.use3)
m422 = Heatmap(as.matrix(middle422), name = '422', 
               column_split = factor(column_subset$biology, levels = c("s", "m", "c")), column_names_side = "top",
               show_heatmap_legend = F, show_row_names = F, cluster_column_slices = F,
               right_annotation = r_422, row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)), col = cols.use3)
c363 = Heatmap(as.matrix(center363), name = '363', 
               column_split = factor(column_subset$biology, levels = c("s", "m", "c")), column_names_side = "top",
               show_heatmap_legend = F, show_row_names = F, cluster_column_slices = F,
               right_annotation = r_363, row_names_gp = gpar(fontsize = c(8)), column_names_gp = gpar(fontsize = c(8)), col = cols.use3)

layers = s644 %v% s535 %v% m422 %v% c363
layers = s644 %v% s535 %v% c363

draw(layers, column_title = "expression per cell cluster", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")

pdf("exports/gene_exp_S_C.pdf", width = 3.5, height = 6.5)
draw(layers, column_title = "expression per cell cluster", column_title_gp = gpar(fontsize = 9), heatmap_legend_side = "bottom")
dev.off()






#rescale data
scaled_mat <- t(scale(t(hallmark_s_644)))
rescaled_mat <- rescale(scaled_mat, to=c(-2,2))

scaled.H.s644 <- Heatmap(rescaled_mat, row_names_side = "right",
                       column_names_side = "top",col = cols.use, show_column_names = T,
                       #david add 2 lines below
                       #column_split=split,
                       #bottom_annotation = HeatmapAnnotation(df = cluster_positions, col = col.biology),
                       cluster_rows = T, cluster_columns = T,
                       heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-Score'))
scaled.H.s644


#RESCALE and plot heatmap (also drop some columns)- hallmark genes in clusters-----
drops <- c("0_s", "2_s", "4_s", "0_o", "1_o", "5_o", "0_b", "1_b")

scaled_1 <- t(scale(t(hallmark_s_644)))
scaled_1 <- rescale(scaled_1, to=c(-2,2))
scaled_1 <- as.data.frame(scaled_1)
scaled_1 <- scaled_1[ , !(names(scaled_1) %in% drops)]

scaled_2 <- t(scale(t(hallmark_s_535)))
scaled_2 <- rescale(scaled_2, to=c(-2,2))
scaled_2 <- as.data.frame(scaled_2)
scaled_2 <- scaled_2[ , !(names(scaled_2) %in% drops)]

scaled_3 <- t(scale(t(hallmark_m_422)))
scaled_3 <- rescale(scaled_3, to=c(-2,2))
scaled_3 <- as.data.frame(scaled_3)
scaled_3 <- scaled_3[ , !(names(scaled_3) %in% drops)]

scaled_5 <- t(scale(t(hallmark_c_363)))
scaled_5 <- rescale(scaled_5, to=c(-2,2))
scaled_5 <- as.data.frame(scaled_5)
scaled_5 <- scaled_5[ , !(names(scaled_5) %in% drops)]

subset_cluster_positions <- data.frame(row.names = colnames(scaled_1))
subset_cluster_positions$biology <- c("m", "c", "s", "s", "s", "c", "m", "s", "m", "c", "s", "s")

scaled_hm.s644 = Heatmap(as.matrix(scaled_1), name = '644', 
                  column_split = factor(subset_cluster_positions$biology, levels = c("s", "m", "c")),
                  row_km = 1, cluster_column_slices = FALSE)

scaled_hm.s535 = Heatmap(as.matrix(scaled_2), name = '535', 
                  column_split = factor(subset_cluster_positions$biology, levels = c("s", "m", "c")),
                  row_km = 1, cluster_column_slices = FALSE)

scaled_hm.m422 = Heatmap(as.matrix(scaled_3), name = '422', 
                  column_split = factor(subset_cluster_positions$biology, levels = c("s", "m", "c")),
                  row_km = 1, cluster_column_slices = FALSE)

scaled_hm.c363 = Heatmap(as.matrix(scaled_5), name = '363', 
                  column_split = factor(subset_cluster_positions$biology, levels = c("s", "m", "c")),
                  row_km = 1, cluster_column_slices = FALSE)
scaled_hallmark_layers = scaled_hm.s644 %v% scaled_hm.s535 %v% scaled_hm.m422 %v% scaled_hm.c363
scaled_hallmark_layers = scaled_hm.s644 %v% scaled_hm.s535

draw(scaled_hallmark_layers)


combined_hallmark = hallmark_layers + hallmark_layers





#all categories into one heatmap with average gene expression per sample per term
# reorganize the list of lists by getting average per sample per term----

surface <- s_644_matrices
set.merged <- data.table()
samples <- colnames(surface644)

for (category in 1:length(surface)){
  
  category.name <- names(surface[category])
  
  for (sub in 1:length(surface[[category]])){
    
    subcat.name <- names(surface[[category]][sub])
    subcat.data <- as.data.frame(surface[[category]][sub])
    names(subcat.data) <- samples
    subcat.colmeans <- as.data.table(t(colMeans(subcat.data)))
    subcat.colmeans[,category:= paste0(category.name)]
    subcat.colmeans[,subcategory:= paste0(subcat.name)]
    
    if(nrow(set.merged) > 0){
      set.merged <- rbind(set.merged, subcat.colmeans)
    }
    if(nrow(set.merged) == 0){set.merged <- subcat.colmeans}
    
  }
  
}

set.merged
#fwrite(set.merged, 'comparison_data/surface644_averages.csv')


# set up row annotation df-----
myrowanno <- as.data.frame(set.merged$category)
rownames(myrowanno) <- set.merged$subcategory
names(myrowanno) <- 'Category'
myrowanno

# set up matrix for plotting----
avg_s644 <- set.merged[,-c('category','subcategory')]
scaled_avg_s644 <- t(scale(t(avg_s644)))
scaled_avg_s644 <- rescale(scaled_avg_s644, to=c(-2,2))
rownames(scaled_avg_s644) <- set.merged$subcategory

#plot
ha <- HeatmapAnnotation(df = myrowanno, which='row', width=unit(1, 'cm'))

h3 <- Heatmap(as.matrix(scaled_avg_s644), row_names_side = "right",
              column_names_side = "top",col = cols.use, show_column_names = T,
              #david add column split
              #column_split = factor(cluster_positions$biology, levels = c("s", "m", "c")),
              cluster_rows = F, cluster_columns = T, row_names_max_width = unit(30, 'cm'),
              row_names_gp = gpar(fontsize = c(9)),
              heatmap_legend_param = list(legend_height = unit(9, "cm"), title='Row Z-Score\nAvg. expression'))
h3

draw(h3 + ha, column_title = "surface644")


#now do this through all matrices
#all categories into one heatmap with average gene expression per sample per term
# reorganize the list of lists by getting average per sample per term----

#s644 matrix averaging -----
through <- s_644_matrices
set.merged <- data.table()
samples <- colnames(surface644)

for (category in 1:length(through)){
  
  category.name <- names(through[category])
  
  for (sub in 1:length(through[[category]])){
    
    subcat.name <- names(through[[category]][sub])
    subcat.data <- as.data.frame(through[[category]][sub])
    names(subcat.data) <- samples
    subcat.colmeans <- as.data.table(t(colMeans(subcat.data)))
    subcat.colmeans[,category:= paste0(category.name)]
    subcat.colmeans[,subcategory:= paste0(subcat.name)]
    
    if(nrow(set.merged) > 0){
      set.merged <- rbind(set.merged, subcat.colmeans)
    }
    if(nrow(set.merged) == 0){set.merged <- subcat.colmeans}
    
  }
  
}

avg_s644 <- set.merged
#fwrite(set.merged, 'comparison_data/surface604_averages.csv')

# set annotation df
myrowanno <- as.data.frame(set.merged$category)
rownames(myrowanno) <- set.merged$subcategory
names(myrowanno) <- 'Category'
s644_anno <- myrowanno


#535 matrix averaging -----
through <- s_535_matrices
set.merged <- data.table()
samples <- colnames(surface644)

for (category in 1:length(through)){
  
  category.name <- names(through[category])
  
  for (sub in 1:length(through[[category]])){
    
    subcat.name <- names(through[[category]][sub])
    subcat.data <- as.data.frame(through[[category]][sub])
    names(subcat.data) <- samples
    subcat.colmeans <- as.data.table(t(colMeans(subcat.data)))
    subcat.colmeans[,category:= paste0(category.name)]
    subcat.colmeans[,subcategory:= paste0(subcat.name)]
    
    if(nrow(set.merged) > 0){
      set.merged <- rbind(set.merged, subcat.colmeans)
    }
    if(nrow(set.merged) == 0){set.merged <- subcat.colmeans}
    
  }
  
}

avg_s535 <- set.merged
#fwrite(set.merged, 'comparison_data/surface535_averages.csv')

# set annotation df
myrowanno <- as.data.frame(set.merged$category)
rownames(myrowanno) <- set.merged$subcategory
names(myrowanno) <- 'Category'
s535_anno <- myrowanno



#m422 matrix averaging -----
through <- m_422_matrices

set.merged <- data.table()
samples <- colnames(surface644)

for (category in 1:length(through)){
  
  category.name <- names(through[category])
  
  for (sub in 1:length(through[[category]])){
    
    subcat.name <- names(through[[category]][sub])
    subcat.data <- as.data.frame(through[[category]][sub])
    names(subcat.data) <- samples
    subcat.colmeans <- as.data.table(t(colMeans(subcat.data)))
    subcat.colmeans[,category:= paste0(category.name)]
    subcat.colmeans[,subcategory:= paste0(subcat.name)]
    
    if(nrow(set.merged) > 0){
      set.merged <- rbind(set.merged, subcat.colmeans)
    }
    if(nrow(set.merged) == 0){set.merged <- subcat.colmeans}
    
  }
  
}

avg_m422 <- set.merged
#fwrite(set.merged, 'comparison_data/middle422_averages.csv')

# set annotation df
myrowanno <- as.data.frame(set.merged$category)
rownames(myrowanno) <- set.merged$subcategory
names(myrowanno) <- 'Category'
m422_anno <- myrowanno


#c363 matrix averaging -----
through <- c_363_matrices

set.merged <- data.table()
samples <- colnames(surface644)

for (category in 1:length(through)){
  
  category.name <- names(through[category])
  
  for (sub in 1:length(through[[category]])){
    
    subcat.name <- names(through[[category]][sub])
    subcat.data <- as.data.frame(through[[category]][sub])
    names(subcat.data) <- samples
    subcat.colmeans <- as.data.table(t(colMeans(subcat.data)))
    subcat.colmeans[,category:= paste0(category.name)]
    subcat.colmeans[,subcategory:= paste0(subcat.name)]
    
    if(nrow(set.merged) > 0){
      set.merged <- rbind(set.merged, subcat.colmeans)
    }
    if(nrow(set.merged) == 0){set.merged <- subcat.colmeans}
    
  }
  
}

avg_c363 <- set.merged
#fwrite(set.merged, 'comparison_data/center363_averages.csv')

# set annotation df
myrowanno <- as.data.frame(set.merged$category)
rownames(myrowanno) <- set.merged$subcategory
names(myrowanno) <- 'Category'
c363_anno <- myrowanno

#save matrices and annotations -----
#save(avg_s644, avg_s535, avg_m422, avg_c363, file = "comparison_data/gene_sets_avg_exp/avg_gene_set_matrixes.rda")
#save(s644_anno, s535_anno, m422_anno, c363_anno, file = "comparison_data/gene_sets_avg_exp/annotations/gene_set_annotations.rda")

# set up gene set matrices for plotting----
#load averages:
load(file = "comparison_data/gene_sets_avg_exp/avg_gene_set_matrixes.rda")
#load annotations
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


#plot avg_gene_set

ha_s644 <- HeatmapAnnotation(df = s644_anno, which='row', width=unit(1, 'cm'))
ha_s535 <- HeatmapAnnotation(df = s535_anno, which='row', width=unit(1, 'cm'))
ha_m422 <- HeatmapAnnotation(df = m422_anno, which='row', width=unit(1, 'cm'))
ha_c363 <- HeatmapAnnotation(df = c363_anno, which='row', width=unit(1, 'cm'))

h1 <- Heatmap(as.matrix(scaled_avg_s644), row_names_side = "right",
              column_names_side = "top",col = cols.use, show_column_names = T,
              show_row_names = F,
              cluster_rows = T, cluster_columns = T, row_names_max_width = unit(30, 'cm'),
              row_names_gp = gpar(fontsize = c(9)))
             #heatmap_legend_param = list(legend_height = unit(9, "cm"), title='Row Z-Score\nAvg. expression'))
h2 <- Heatmap(as.matrix(scaled_avg_s535), row_names_side = "right",
              column_names_side = "top",col = cols.use, show_column_names = T,
              show_row_names = F,
              cluster_rows = T, cluster_columns = T, row_names_max_width = unit(30, 'cm'),
              row_names_gp = gpar(fontsize = c(9)))
              #heatmap_legend_param = list(legend_height = unit(9, "cm"), title='Row Z-Score\nAvg. expression'))
h3 <- Heatmap(as.matrix(scaled_avg_m422), row_names_side = "right",
              column_names_side = "top",col = cols.use, show_column_names = T,
              show_row_names = F,
              cluster_rows = T, cluster_columns = T, row_names_max_width = unit(30, 'cm'),
              row_names_gp = gpar(fontsize = c(9)))
              #heatmap_legend_param = list(legend_height = unit(9, "cm"), title='Row Z-Score\nAvg. expression'))
h5 <- Heatmap(as.matrix(scaled_avg_c363), row_names_side = "right",
              column_names_side = "top",col = cols.use, show_column_names = T,
              show_row_names = F,
              cluster_rows = T, cluster_columns = T, row_names_max_width = unit(30, 'cm'),
              row_names_gp = gpar(fontsize = c(9)))
              #heatmap_legend_param = list(legend_height = unit(9, "cm"), title='Row Z-Score\nAvg. expression'))


draw(h1 %v% h2 %v% h3 %v% h5 , column_title = "surface644")








#end------