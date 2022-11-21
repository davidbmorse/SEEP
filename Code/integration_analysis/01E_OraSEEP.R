rm(list=ls())
setwd("C://PROJECTS/P2022/SEEP_Manuscript")

# data
load("data/SEEP_Regions.rda")
load('data/SEEP_Markers.rda')

markers = seep_regionMarkers
isoreg = iso2D_regions

library(fgsea); library(Seurat); library(RColorBrewer); library(tidyverse);
library(igraph)

DefaultAssay(isoreg) <- "RNA"

lsig <- split(markers, markers$Comparison_cluster)
nsig = sapply(lsig,nrow)
lfeat = lapply(lsig, function(x) x$Gene)

load("Input/MSigDB_subsets.rda")
univer = rownames(isoreg)
db = L_subsets$h.all
lora = lapply(lfeat, function(x) fora(genes=x, pathways=db, universe=univer, minSize=5))
lorasig = lapply(lora, function(x) x[x$padj <= 0.05])
psig = sapply(lorasig, nrow)

lorasig = lapply(lorasig, function(x) {
  y = x[x$overlap >= 1,]
  y$percent = y$overlap/y$size*100
  return(y)
})

lcoll = lorasig
for(i in 1:length(lorasig)){

  lcoll[[i]] =
    collapsePathwaysORA(foraRes=lorasig[[i]], pathways=db, universe=univer,
                        genes=lfeat[[i]])
}

lparent = lapply(names(lcoll), function(x){
  parent = lcoll[[x]]$parentPathways
  parent = ifelse(is.na(parent), names(parent), parent)
  if(length(parent)>0){
  y = data.frame(from=parent, to=names(parent), cluster=x)
  rownames(y) = NULL
  } else {
    y = data.frame(from=NA,to=NA,cluster=x)
    rownames(y) = NULL
  }
  return(y)
})
names(lparent) = names(lcoll)

lorasig_main = lapply(names(lorasig), function(x) {
  y = lorasig[[x]][lorasig[[x]]$pathway %in% lcoll[[x]]$mainPathways,]
  return(y[y$overlap>=5,])
})
names(lorasig_main) = names(lorasig)

save(lorasig_main, lcoll, lparent, lsig, lorasig, file='data/SEEP_ORA.rda')

