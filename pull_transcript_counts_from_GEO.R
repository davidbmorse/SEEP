#importing supplimental data files from GEO
#install and load GEOquery from bioconductor (and BiocManager if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("GEOquery")

library(GEOquery)

#download counts matricies from GEO accession GSE157299

spheroid_filepaths = getGEOSuppFiles("GSE157299", filter_regex = 'spheroid')
organoid_filepaths = getGEOSuppFiles("GSE157299", filter_regex = 'organoid')
biopsy_filepaths = getGEOSuppFiles("GSE157299", filter_regex = 'biopsy')

#Aleksandra, can you code this so that each 'set of files' (should be in groups of three) is placed in a seperate directory? This way we can more easily load them into Seurat using the 'load10X' function