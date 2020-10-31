#David Morse 8.13.19

library(Seurat)
library(dplyr)
library(Matrix)
library(data.table)

setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Data/sequencing_runs/biopsy2_3_CRUKJuly2019_mouse_human/sequence data/run_190911_SEEP_m1m2m3_MERGE")
getwd()

# step 1: read data matrix from inDrop pipeline which are in a Gene(X) by Cell(Y) matrix, S1.counts.tsv
#fread is a function from library(data.table) that reads in a table quickly

rawDataCenter <- fread("m1/m1.counts.tsv.gz",header=TRUE)
rawDataMiddle <- fread("m2/m2.counts.tsv.gz",header=TRUE)
rawDataSurface <- fread("m3/m3.counts.tsv.gz",header=TRUE)

#change the barcode names to make layer specific
rawDataCenter$barcode = paste0('c_', rawDataCenter$barcode)
rawDataMiddle$barcode = paste0('m_', rawDataMiddle$barcode)
rawDataSurface$barcode = paste0('s_', rawDataSurface$barcode)

#converting input matrix (inDrop form) into Seurat compatable files
#create barcodes.tsv by taking the first column from inDrop data and placing it in 'sparse_S1' folder
write.table(rawDataCenter[,1], file= "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/counts_matrices/4th_MERGE_mouses/m1_MERGE/barcodes.tsv",
            quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)

write.table(rawDataMiddle[,1], file= "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/counts_matrices/4th_MERGE_mouses/m2_MERGE/barcodes.tsv",
            quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)

write.table(rawDataSurface[,1], file= "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/counts_matrices/4th_MERGE_mouses/m3_MERGE/barcodes.tsv",
            quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)

#extract the firt row (gene name) from the unformatted data
genesO2center <- as.character(colnames(rawDataCenter))
genesO2center <- data.frame(geneName = genesO2center, geneName2 = genesO2center)

#save a as a file
write.table(genesO2center[-1,], file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/counts_matrices/4th_MERGE_mouses/m1_MERGE/genes.tsv",
            row.names=FALSE, col.names = FALSE, sep="\t", quote = FALSE)

# convert rawData into a matrix which lacks the column 1 labels.  Make this into a sparse matrix by using the 'Matrix' 
# function (capital M). 
dataO2center <- Matrix(data.matrix(data.frame(rawDataCenter[,-1])))
dataO2middle <- Matrix(data.matrix(data.frame(rawDataMiddle[,-1])))
dataO2surface <- Matrix(data.matrix(data.frame(rawDataSurface[,-1])))

#save sparse_matrix as mtx file and TRASPOSE (t) it so that it is now Cell(X) by Genes(Y)
writeMM(t(dataO2center), "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/counts_matrices/4th_MERGE_mouses/m1_MERGE/matrix.mtx")
writeMM(t(dataO2middle), "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/counts_matrices/4th_MERGE_mouses/m2_MERGE/matrix.mtx")
writeMM(t(dataO2surface), "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/counts_matrices/4th_MERGE_mouses/m3_MERGE/matrix.mtx")


#ABOVE STEP TOO BIG FOR MY MATRIXES: FEB 24, 2019
#tried to do this on Denali - sweet! worked :)

#Test if this file works by trying out the 'Read10X' function and seeing if it will load the file
# Load the 'Seurat prepared' dataset
spheroid.data <- Read10X("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/counts_matricies/m3_surface")

#Examine the memory savings between regular and sparse matrices
dense.size <- object.size(as.matrix(spheroid.data))
dense.size
sparse.size <- object.size(spheroid.data)
sparse.size
dense.size/sparse.size


###code for perfoming transformation on Denali cluster:

#Data location:
#/home/morsedb/R/data/indrop/organoidComb/dm1.counts.tsv.gz
#Output:
#/home/morsedb/R/data/indrop/output/center

rawDataCenter <- fread("organoidComb/dm1.counts.tsv.gz",header=TRUE)
rawDataMiddle <- fread("organoidComb/dm2.counts.tsv.gz",header=TRUE)
rawDataSurface <- fread("organoidComb/dm3.counts.tsv.gz",header=TRUE)

rawDataCenter$barcode = paste0('c_', rawDataCenter$barcode)
rawDataMiddle$barcode = paste0('m_', rawDataMiddle$barcode)
rawDataSurface$barcode = paste0('s_', rawDataSurface$barcode)

write.table(rawDataCenter[,1], file= "output/center/barcodes.tsv",
            quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
write.table(rawDataMiddle[,1], file= "output/middle/barcodes.tsv",
            quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
write.table(rawDataSurface[,1], file= "output/surface/barcodes.tsv",
            quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)

writeMM(t(dataO2center), "output/center/matrix.mtx")
writeMM(t(dataO2middle), "output/middle/matrix.mtx")
writeMM(t(dataO2surface), "output/surface/matrix.mtx")




#now read in data to work on
#/home/morsedb/R/data/indrop/organoidComb/working_files/oC_center
#Output:
#/home/morsedb/R/data/indrop/output/center

raw_center <- Read10X("working_files/oC_center")











