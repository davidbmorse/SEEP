
library(Seurat)
library(dplyr)
library(Matrix)
library(data.table) #for extra data table wrangling

#1 load and scale data in Seurat -----------------
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/")
getwd()

#load preprocessed data set and work with scaled data for the GSEA input -------
PDX_biopsy2_3 <- readRDS(file = "R_exports/r_data/4th_merge_mouse/PDX_biopsy2_3_MERGE_final.rds")
PDX_biopsy2_3 <- SetIdent(PDX_biopsy2_3, value = 'orig.ident')


# 2 add 'ALL' ident - 'ALL' simply denotes 'all cells ---------
# load("mergedData_scaled.rda") # actually called 'organoid2' in this case
# a) add ALL ident
x = PDX_biopsy2_3
x = SetIdent(x, value = rep("ALL",nrow(x@meta.data)))
inf = x@meta.data
nam=sapply(strsplit(rownames(inf), "_"), function(x) paste(x[1],x[2], sep="_"))
nam=paste("ALL", nam,1:nrow(inf), sep="_")
x <- RenameCells(x, add.cell.id = nam)
mergedData <- merge(PDX_biopsy2_3, x, project="ALL")
mergedData$new.ident <- Idents(mergedData)

save(PDX_biopsy2_3, file = "R_scripts/4th_MERGE_mouses/190916_GSEA_PDX_biopsy2_3/data/scaledData.rda")
save(mergedData, file = "R_scripts/4th_MERGE_mouses/190916_GSEA_PDX_biopsy2_3/data/mergedData.rda")
