#importing supplimental data files from GEO
#install and load GEOquery from bioconductor (and BiocManager if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("GEOquery")

library(GEOquery)
library(filesstrings)

#download counts matricies from GEO accession GSE157299 and place in directories each containing one gene list, cell barcode list, and .mtx compressed matrix
spheroid_filepaths = getGEOSuppFiles("GSE157299", filter_regex = 'spheroid')
file.rename("GSE157299", "spheroid_10X")

biopsy_filepaths = getGEOSuppFiles("GSE157299", filter_regex = 'biopsy')
file.rename("GSE157299", "biopsy_10X")

#for organoid, we have to divide up matricies based on layer
organoid_filepaths = getGEOSuppFiles("GSE157299", filter_regex = 'organoid')
dir.create("organoid_10X")
dir.create("organoid_10X/center")
dir.create("organoid_10X/middle")
dir.create("organoid_10X/surface")

#select and move files based on positions
file.move(paste("GSE157299/", list.files(path = "./GSE157299", pattern = "center"), sep=""), "organoid_10X/center")
file.move(paste("GSE157299/", list.files(path = "./GSE157299", pattern = "middle"), sep=""), "organoid_10X/middle")
file.move(paste("GSE157299/", list.files(path = "./GSE157299", pattern = "surface"), sep=""), "organoid_10X/surface")

#delete now-empty GSE157299 directory
unlink("GSE157299", recursive = TRUE)

