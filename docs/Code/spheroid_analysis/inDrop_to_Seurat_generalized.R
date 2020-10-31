#David Morse 4.12.19

library(Seurat)
library(dplyr)
library(Matrix)
library(data.table)

setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/reAnalysis_190412")
getwd()

# step 1: re-format data (already in sparse matric form) so that the barcode names are unique per layer

/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/reAnalysis_190412/counts_matrices/spheroid_central/CENTRAL/barcodes.tsv

rawDataCenter <- fread("counts_matrices/spheroid_central/CENTRAL/barcodes.tsv", header = FALSE)
rawDataInner <- fread("counts_matrices/spheroid_inner/INNER/barcodes.tsv",header= FALSE)
rawDataOuter <- fread("counts_matrices/spheroid_outer/OUTER/barcodes.tsv",header= FALSE)
rawDataSurface <- fread("counts_matrices/spheroid_surface/SURFACE/barcodes.tsv",header= FALSE)

#change the barcode names to make layer specific
rawDataCenter$V1 = paste0('c_', rawDataCenter$V1)
rawDataInner$V1 = paste0('i_', rawDataInner$V1)
rawDataOuter$V1 = paste0('o_', rawDataOuter$V1)
rawDataSurface$V1 = paste0('s_', rawDataSurface$V1)

#create barcodes.tsv by writing the single column bacrode table to a file
write.table(rawDataCenter[,1], file= "counts_matrices/spheroid_central/CENTRAL/barcodes.tsv",
            quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
write.table(rawDataInner[,1], file= "counts_matrices/spheroid_inner/INNER/barcodes.tsv",
            quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
write.table(rawDataOuter[,1], file= "counts_matrices/spheroid_outer/OUTER/barcodes.tsv",
            quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
write.table(rawDataSurface[,1], file= "counts_matrices/spheroid_surface/SURFACE/barcodes.tsv",
            quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)








#EXTRA CODE HERE
#extract the firt row (gene name) from the unformatted data
genesO2center <- as.character(colnames(rawDataO2center))
genesO2center <- data.frame(geneName = genesO2center, geneName2 = genesO2center)

#save a as a file
write.table(genesO2center[-1,], file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoids/counts_matrices/organoids2/o2_center/genes.tsv",
            row.names=FALSE, col.names = FALSE, sep="\t", quote = FALSE)

# convert rawData into a matrix which lacks the column 1 labels.  Make this into a sparse matrix by using the 'Matrix' 
# function (capital M). 
dataO2center <- Matrix(data.matrix(data.frame(rawDataO2center[,-1])))
dataO2middle <- Matrix(data.matrix(data.frame(rawDataO2middle[,-1])))
dataO2surface <- Matrix(data.matrix(data.frame(rawDataO2surface[,-1])))

#save sparse_matrix as mtx file and TRASPOSE (t) it so that it is now Cell(X) by Genes(Y)
writeMM(t(dataO2center), "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoids/counts_matrices/organoids2/o2_center/matrix.mtx")
writeMM(t(dataO2middle), "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoids/counts_matrices/organoids2/o2_center/matrix.mtx")
writeMM(t(dataO2surface), "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoids/counts_matrices/organoids2/o2_center/matrix.mtx")


#ABOVE STEP TOO BIG FOR MY MATRIXES: FEB 24, 2019
#tried to do this on Denali - sweet! worked :)

#Test if this file works by trying out the 'Read10X' function and seeing if it will load the file
# Load the 'Seurat prepared' dataset
spheroid.data <- Read10X("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/counts_matrices/oC_middle")

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

rawDataO2center <- fread("organoidComb/dm1.counts.tsv.gz",header=TRUE)
rawDataO2middle <- fread("organoidComb/dm2.counts.tsv.gz",header=TRUE)
rawDataO2surface <- fread("organoidComb/dm3.counts.tsv.gz",header=TRUE)

rawDataO2center$barcode = paste0('c_', rawDataO2center$barcode)
rawDataO2middle$barcode = paste0('m_', rawDataO2middle$barcode)
rawDataO2surface$barcode = paste0('s_', rawDataO2surface$barcode)

write.table(rawDataO2center[,1], file= "output/center/barcodes.tsv",
            quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
write.table(rawDataO2middle[,1], file= "output/middle/barcodes.tsv",
            quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
write.table(rawDataO2surface[,1], file= "output/surface/barcodes.tsv",
            quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)

writeMM(t(dataO2center), "output/center/matrix.mtx")
writeMM(t(dataO2middle), "output/middle/matrix.mtx")
writeMM(t(dataO2surface), "output/surface/matrix.mtx")




#now read in data to work on
#/home/morsedb/R/data/indrop/organoidComb/working_files/oC_center
#Output:
#/home/morsedb/R/data/indrop/output/center

raw_center <- Read10X("working_files/oC_center")











