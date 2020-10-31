#group modules across biology

setwd('/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/comparative_analysis_filt_biop/clusters_heatmap')
library(data.table)
library(ComplexHeatmap)
library(scales)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
cols.use <- colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100) # reversed RdBu, creates blue/white/red
cols.use3 <- colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdGy")))(100) # reversed RdBu, creates blue/white/red

#load data and annotations
s_644_matrices <- readRDS("comparison_data/surface644.rds")
s_535_matrices <- readRDS("comparison_data/surface535.rds")
m_422_matrices <- readRDS("comparison_data/middle422.rds")
c_363_matrices <- readRDS("comparison_data/center363.rds")
surface644_annotation <- readRDS('comparison_data/surface_644_annotation.rds')
surface535_annotation <- readRDS('comparison_data/surface_535_annotation.rds')
middle422_annotation <- readRDS('comparison_data/middle_422_annotation.rds')
center363_annotation <- readRDS('comparison_data/center_363_annotation.rds')

surface644 <- readRDS('comparison_data/all_genes_unique_surface644.rds')
surface535 <- readRDS('comparison_data/all_genes_unique_surface535.rds')
middle422 <- readRDS('comparison_data/all_genes_unique_middle422.rds')
center363 <- readRDS('comparison_data/all_genes_unique_center363.rds')

#extract only hallmark genes
hallmark_s_644 <- unique(do.call(rbind, s_644_matrices[[1]]))
hallmark_s_535 <- unique(do.call(rbind, s_535_matrices[[1]]))
hallmark_m_422 <- unique(do.call(rbind, m_422_matrices[[1]]))
hallmark_c_363 <- unique(do.call(rbind, c_363_matrices[[1]]))

#all categories into one heatmap with average gene expression per sample per term. Reorganize the list of lists by getting average per sample per term

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
fwrite(set.merged, 'comparison_data/surface640_averages.csv')

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
fwrite(set.merged, 'comparison_data/surface535_averages.csv')

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
fwrite(set.merged, 'comparison_data/middle422_averages.csv')

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
fwrite(set.merged, 'comparison_data/center363_averages.csv')

# set annotation df
myrowanno <- as.data.frame(set.merged$category)
rownames(myrowanno) <- set.merged$subcategory
names(myrowanno) <- 'Category'
c363_anno <- myrowanno

#save matrices and annotations -----
save(avg_s644, avg_s535, avg_m422, avg_c363, file = "comparison_data/gene_sets_avg_exp/avg_gene_set_matrixes.rda")
save(s644_anno, s535_anno, m422_anno, c363_anno, file = "comparison_data/gene_sets_avg_exp/annotations/gene_set_annotations.rda")















