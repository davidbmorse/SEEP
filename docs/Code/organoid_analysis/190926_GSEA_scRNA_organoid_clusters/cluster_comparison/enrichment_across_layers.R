#Find LE of gene sets across biology
#load libraries
library(RColorBrewer)
library(lattice)
library(latticeExtra)
library(grid)
library(gridExtra)
library(Seurat)
library(DT)
#Summary infromation for common gene sets enriched in surface layers of spehroids, organoids, and biopsies

#HALLMARK-------
#SURFACE
    #ALL
HALLMARK_APOPTOSIS
HALLMARK_TNFA_SIGNALING_VIA_NFKB
    #SPHERE/ORG
HALLMARK_INTERFERON_ALPHA_RESPONSE
HALLMARK_INTERFERON_GAMMA_RESPONSE
HALLMARK_ALLOGRAFT_REJECTION
HALLMARK_KRAS_SIGNALING_DN
HALLMARK_KRAS_SIGNALING_UP
    #ORG/BIOPSY
HALLMARK_ANGIOGENESIS

#CENTER
    #ALL
HALLMARK_SPERMATOGENESIS
HALLMARK_MTORC1_SIGNALING
HALLMARK_PROTEIN_SECRETION
HALLMARK_G2M_CHECKPOINT
HALLMARK_APICAL_JUNCTION
    #ORG/BIOPSY
HALLMARK_MYC_TARGETS_V1

#C2.CGP ---------
#SURFACE
    #ALL
PECE_MAMMARY_STEM_CELL_UP
FARMER_BREAST_CANCER_CLUSTER_1
SANA_TNF_SIGNALING_UP
ZHANG_INTERFERON_RESPONSE
TIAN_TNF_SIGNALING_VIA_NFKB
    #SPHERE/ORG
BROWNE_INTERFERON_RESPONSIVE_GENES
    #ORG/BIOPSY
BOSCO_EPITHELIAL_DIFFERENTIATION_MODULE

#C6.ALL -------
#SURFACE
    #ALL
HINATA_NFKB_IMMU_INF
KRAS.LUNG.BREAST_UP.V1_UP
EGFR_UP.V1_UP

#C5.BP --------
#SURFACE
    #ALL
GO_PROTEIN_LOCALIZATION_TO_ENDOPLASMIC_RETICULUM
GO_PROTEIN_TARGETING_TO_MEMBRANE












#set global parameters - "gene set" and "location" ---------
collection="h.all"
pathway="HALLMARK_APOPTOSIS"
loc="SURFACE"
#for plotting
ds=5
asp='fill'

#plotting parameters --------
threshold = 0.35 #max(ceiling(range(dat1)))-2
scale.range= c(-threshold, threshold)
zak= seq(scale.range[1],scale.range[2],0.01)
nzak=length(zak)
pal = rev(brewer.pal(11,'Spectral'))
pal[6] = 'white'
pal=c('navyblue',pal[1:6],'white',pal[7:11],'#3f003f')
pal = c(rev(c(brewer.pal(9,"Blues"),'navyblue')), c(brewer.pal(9,'Reds'),'#3f003f'))
colo = colorRampPalette(pal, space="Lab")(nzak) 
#breaks = c(seq(-35, 10, length.out=NBR.COLORS/2), 10, seq(10, 35, length.out=NBR.COLORS/2))

#Load Spheroid data ---------
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_layers")
source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_layers/Codes/01GSEA_pGSEAfunction.R")
data(pGSEA)
padj=c(0.01,0.05, 0.1)
tab = rbind(
  CENTRAL=count.sig(gsea_CENTRAL, fdr.level=padj, group.name='CENTRAL'),
  INNER=count.sig(gsea_INNER, fdr.level=padj, group.name='INNER'),
  OUTER=count.sig(gsea_OUTER, fdr.level=padj, group.name='OUTER'),
  SURFACE=count.sig(gsea_SURFACE, fdr.level=padj, group.name='SURFACE'),
  ALL.Tests=count.sig(gsea_SURFACE, fdr.level=1, group.name='N pathways tested')
)
rownames(tab) = NULL
datatable(tab, options = list(pageLength = 13, dom = 'tip'), rownames = FALSE)

data(Preranked)
data(Data_scaled)
data(pGSEA)
source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_layers/Codes/01GSEA_pGSEAfunction.R")
dx=Data@data
mm=t(apply(dx,1,tapply, Data@ident, calc.logMean))
mms = apply(dx,1,calc.logMean)
mm1= data.frame(mm,TOTAL=mms)
dd_spheroid = mm - mms
dd1=mm1-mms

setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_layers/MSigDB")
load("MSigDB_subsets.rda")
fid=which(gsea_SURFACE$h.all$pathway==pathway)
gseaRes=gsea_SURFACE$h.all[fid, ]
nameGset = gseaRes$pathway
gset = L_subsets[[collection]][[nameGset]]
ranks = difs$SURFACE
names(ranks) = rownames(difs)
LE_sphere = sort(unlist(gseaRes$leadingEdge))
plot.ES(rank.statistic = ranks, gsea=gseaRes, gene.set=gset, name=nameGset, pathwaySource=collection, location=loc)

#create sub-matrix for plotting
mat_sphere=t(dd_spheroid[LE_sphere,])

#plot average difference heatmap 
rowv = as.dendrogram(hclust(dist(mat_sphere),method='ward.D'))
colv = as.dendrogram(hclust(dist(t(mat_sphere)),method='ward.D'))
row.ord = order.dendrogram(rowv)
col.ord = order.dendrogram(colv)
sdat = apply(mat_sphere,2,function(x) ifelse(x<=-threshold, -threshold, x))
sdat = t(apply(sdat,2,function(x) ifelse(x>=threshold, threshold, x)))
#rownames(sdat) = substring(rownames(sdat),1,1) #Dave jan 7th - remove to allow full gene names in the rows
p=levelplot(sdat[col.ord, c("SURFACE", "OUTER", "INNER", "CENTRAL")], col.regions=colo, at=zak,  # THIS NEEDS TO ADJUST THE DENDROGRAM TOO - AS OF NOW, I DON'T THINK IT DOES
            aspect = asp, xlab='', ylab='',
            scales = list(x = list(rot = 90, cex=.23, lwd=0), tck=c(0,0), alternating=2, y = list(cex=0.6)), 
            colorkey = list(space='left', labels=list(cex=0.6)),
            legend = list(
              right = list(fun = dendrogramGrob, args = list(x = rowv, ord = row.ord, side = "right", size = 1)),
              top = list(fun = dendrogramGrob, args = list(x = colv, ord=col.ord, side = "top", size = 2, colour='white'))))
p.diff = p

# plot comparison heatmap only
grid.arrange(p.diff)



#Load organoid data ---------
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190429_GSEA_scRNA_organoid_layers")
source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190429_GSEA_scRNA_organoid_layers/Codes/01GSEA_pGSEAfunction.R")
data(pGSEA)
padj=c(0.01,0.05, 0.1)
tab = rbind(
  CENTRAL=count.sig(gsea_CENTRAL, fdr.level=padj, group.name='CENTRAL'),
  MIDDLE=count.sig(gsea_MIDDLE, fdr.level=padj, group.name='MIDDLE'),
  SURFACE=count.sig(gsea_SURFACE, fdr.level=padj, group.name='SURFACE'),
  ALL.Tests=count.sig(gsea_SURFACE, fdr.level=1, group.name='N pathways tested')
)
rownames(tab) = NULL

data(PrerankOrg)
data(scaledDataGSEA)
data(pGSEA)
source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190429_GSEA_scRNA_organoid_layers/Codes/01GSEA_pGSEAfunction.R")
dx=organoid_GSEA@assays$RNA@data
mm=t(apply(dx,1,tapply, organoid_GSEA@active.ident, calc.logMean))
mms = apply(dx,1,calc.logMean)
mm1= data.frame(mm,TOTAL=mms)
dd_organoid = mm - mms
dd1=mm1-mms

setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190429_GSEA_scRNA_organoid_layers/MSigDB")
load("MSigDB_subsets.rda")
fid=which(gsea_SURFACE$h.all$pathway==pathway)
gseaRes=gsea_SURFACE$h.all[fid, ]
nameGset = gseaRes$pathway
gset = L_subsets[[collection]][[nameGset]]
ranks = difs$s
names(ranks) = rownames(difs)
LE_organoid = sort(unlist(gseaRes$leadingEdge))
plot.ES(rank.statistic = ranks, gsea=gseaRes, gene.set=gset, name=nameGset, pathwaySource=collection, location=loc)

#create sub-matrix for plotting
mat_organoid=t(dd_organoid[LE_organoid,])

#plot average difference heatmap 
rowv = as.dendrogram(hclust(dist(mat_organoid),method='ward.D'))
colv = as.dendrogram(hclust(dist(t(mat_organoid)),method='ward.D'))
row.ord = order.dendrogram(rowv)
col.ord = order.dendrogram(colv)
sdat = apply(mat_organoid,2,function(x) ifelse(x<=-threshold, -threshold, x))
sdat = t(apply(sdat,2,function(x) ifelse(x>=threshold, threshold, x)))
#rownames(sdat) = substring(rownames(sdat),1,1) #Dave jan 7th - remove to allow full gene names in the rows
p=levelplot(sdat[col.ord, c("s", "m", "c")], col.regions=colo, at=zak,  # THIS NEEDS TO ADJUST THE DENDROGRAM TOO - AS OF NOW, I DON'T THINK IT DOES
            aspect = asp, xlab='', ylab='',
            scales = list(x = list(rot = 90, cex=.23, lwd=0), tck=c(0,0), alternating=2, y = list(cex=0.6)), 
            colorkey = list(space='left', labels=list(cex=0.6)),
            legend = list(
              right = list(fun = dendrogramGrob, args = list(x = rowv, ord = row.ord, side = "right", size = 1)),
              top = list(fun = dendrogramGrob, args = list(x = colv, ord=col.ord, side = "top", size = 2, colour='white'))))
p.diff = p

# plot comparison heatmap only
grid.arrange(p.diff)


#Load biopsy data ---------
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/190916_GSEA_PDX_biopsy2_3")
source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/190916_GSEA_PDX_biopsy2_3/Codes/01GSEA_pGSEAfunction.R")
data(pGSEA)
padj=c(0.01,0.05, 0.1)
tab = rbind(
  CENTRAL=count.sig(gsea_CENTRAL, fdr.level=padj, group.name='CENTRAL'),
  MIDDLE=count.sig(gsea_MIDDLE, fdr.level=padj, group.name='MIDDLE'),
  SURFACE=count.sig(gsea_SURFACE, fdr.level=padj, group.name='SURFACE'),
  ALL.Tests=count.sig(gsea_SURFACE, fdr.level=1, group.name='N pathways tested')
)
rownames(tab) = NULL

data(PrerankBiopsy)
data(scaledDataGSEA)
data(pGSEA)
source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/190916_GSEA_PDX_biopsy2_3/Codes/01GSEA_pGSEAfunction.R")
dx=PDX_biopsy2_3_GSEA@assays$RNA@data
mm=t(apply(dx,1,tapply, PDX_biopsy2_3_GSEA@active.ident, calc.logMean))
mms = apply(dx,1,calc.logMean)
mm1= data.frame(mm,TOTAL=mms)
dd_biop = mm - mms
dd1=mm1-mms

setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/190916_GSEA_PDX_biopsy2_3/MSigDB")
load("MSigDB_subsets.rda")
fid=which(gsea_SURFACE$h.all$pathway==pathway)
gseaRes=gsea_SURFACE$h.all[fid, ]
nameGset = gseaRes$pathway
gset = L_subsets[[collection]][[nameGset]]
ranks = difs$s
names(ranks) = rownames(difs)
LE_biopsy = sort(unlist(gseaRes$leadingEdge))
#plot.ES(rank.statistic = ranks, gsea=gseaRes, gene.set=gset, name=nameGset, pathwaySource=collection, location=loc)

#create sub-matrix for plotting
mat_biopsy=t(dd_biop[LE_biopsy,])

#plot average difference heatmap 
rowv = as.dendrogram(hclust(dist(mat_biopsy),method='ward.D'))
colv = as.dendrogram(hclust(dist(t(mat_biopsy)),method='ward.D'))
row.ord = order.dendrogram(rowv)
col.ord = order.dendrogram(colv)
sdat = apply(mat_biopsy,2,function(x) ifelse(x<=-threshold, -threshold, x))
sdat = t(apply(sdat,2,function(x) ifelse(x>=threshold, threshold, x)))
#rownames(sdat) = substring(rownames(sdat),1,1) #Dave jan 7th - remove to allow full gene names in the rows
p=levelplot(sdat[col.ord, c("s", "m", "c")], col.regions=colo, at=zak,  # THIS NEEDS TO ADJUST THE DENDROGRAM TOO - AS OF NOW, I DON'T THINK IT DOES
            aspect = asp, xlab='', ylab='',
            scales = list(x = list(rot = 90, cex=.23, lwd=0), tck=c(0,0), alternating=2, y = list(cex=0.6)), 
            colorkey = list(space='left', labels=list(cex=0.6)),
            legend = list(
              right = list(fun = dendrogramGrob, args = list(x = rowv, ord = row.ord, side = "right", size = 1)),
              top = list(fun = dendrogramGrob, args = list(x = colv, ord=col.ord, side = "top", size = 2, colour='white'))))
p.diff = p

# plot comparison heatmap only
grid.arrange(p.diff)








#find LE overlap between all biology
Common_S_H_Apoptosis <- intersect(LE_sphere, LE_organoid)
Common_S_H_Apoptosis <- intersect(Common_S_H_Apoptosis, LE_biopsy)


#create sub-matricies for plotting
mat_sphere=t(dd_spheroid[Common_S_H_Apoptosis,])
mat_organoid=t(dd_organoid[Common_S_H_Apoptosis,])
mat_biopsy=t(dd_biop[Common_S_H_Apoptosis,])

#for now - dendogram based on BIOPSY gene ordering
#plot average difference heatmap 
rowv = as.dendrogram(hclust(dist(mat_biopsy),method='ward.D'))
colv = as.dendrogram(hclust(dist(t(mat_biopsy)),method='ward.D'))
row.ord = order.dendrogram(rowv)
col.ord = order.dendrogram(colv)
sdat = apply(mat_biopsy,2,function(x) ifelse(x<=-threshold, -threshold, x))
sdat = t(apply(sdat,2,function(x) ifelse(x>=threshold, threshold, x)))
#rownames(sdat) = substring(rownames(sdat),1,1) #Dave jan 7th - remove to allow full gene names in the rows
p=levelplot(sdat[col.ord, c("s", "m", "c")], col.regions=colo, at=zak,  # THIS NEEDS TO ADJUST THE DENDROGRAM TOO - AS OF NOW, I DON'T THINK IT DOES
            aspect = asp, xlab='', ylab='',
            scales = list(x = list(rot = 90, cex=.23, lwd=0), tck=c(0,0), alternating=2, y = list(cex=0.6)), 
            colorkey = list(space='left', labels=list(cex=0.6)),
            legend = list(
              right = list(fun = dendrogramGrob, args = list(x = rowv, ord = row.ord, side = "right", size = 1)),
              top = list(fun = dendrogramGrob, args = list(x = colv, ord=col.ord, side = "top", size = 2, colour='white'))))
p.biop = p
grid.arrange(p.biop)

#organoid
#plot average difference heatmap 
rowv = as.dendrogram(hclust(dist(mat_organoid),method='ward.D'))
colv = as.dendrogram(hclust(dist(t(mat_organoid)),method='ward.D'))
#row.ord = order.dendrogram(rowv)
#col.ord = order.dendrogram(colv)
sdat = apply(mat_organoid,2,function(x) ifelse(x<=-threshold, -threshold, x))
sdat = t(apply(sdat,2,function(x) ifelse(x>=threshold, threshold, x)))
#rownames(sdat) = substring(rownames(sdat),1,1) #Dave jan 7th - remove to allow full gene names in the rows
p=levelplot(sdat[col.ord, c("s", "m", "c")], col.regions=colo, at=zak,  # THIS NEEDS TO ADJUST THE DENDROGRAM TOO - AS OF NOW, I DON'T THINK IT DOES
            aspect = asp, xlab='', ylab='',
            scales = list(x = list(rot = 90, cex=.23, lwd=0), tck=c(0,0), alternating=2, y = list(cex=0.6)), 
            colorkey = list(space='left', labels=list(cex=0.6)),
            legend = list(
              right = list(fun = dendrogramGrob, args = list(x = rowv, ord = row.ord, side = "right", size = 1)),
              top = list(fun = dendrogramGrob, args = list(x = colv, ord=col.ord, side = "top", size = 0, colour='white'))))
p.org = p
grid.arrange(p.org)

#spheroid
#plot average difference heatmap 
rowv = as.dendrogram(hclust(dist(mat_sphere),method='ward.D'))
colv = as.dendrogram(hclust(dist(t(mat_sphere)),method='ward.D'))
#row.ord = order.dendrogram(rowv)
#col.ord = order.dendrogram(colv)
sdat = apply(mat_sphere,2,function(x) ifelse(x<=-threshold, -threshold, x))
sdat = t(apply(sdat,2,function(x) ifelse(x>=threshold, threshold, x)))
#rownames(sdat) = substring(rownames(sdat),1,1) #Dave jan 7th - remove to allow full gene names in the rows
colnames(sdat) = substring(colnames(sdat),1,1)
p=levelplot(sdat[col.ord, c("S", "O", "I", "C")], col.regions=colo, at=zak,  # THIS NEEDS TO ADJUST THE DENDROGRAM TOO - AS OF NOW, I DON'T THINK IT DOES
            aspect = asp, xlab='', ylab='',
            scales = list(x = list(rot = 90, cex=.23, lwd=0), tck=c(0,0), alternating=2, y = list(cex=0.6)), 
            colorkey = list(space='left', labels=list(cex=0.6)),
            legend = list(
              right = list(fun = dendrogramGrob, args = list(x = rowv, ord = row.ord, side = "right", size = 1)),
              top = list(fun = dendrogramGrob, args = list(x = colv, ord=col.ord, side = "top", size = 0, colour='white'))))
p.sphere = p
grid.arrange(p.sphere)

#plot all 3

grid.arrange(p.biop, p.org, p.sphere, ncol=1, heights=c(1,0.84,1))


