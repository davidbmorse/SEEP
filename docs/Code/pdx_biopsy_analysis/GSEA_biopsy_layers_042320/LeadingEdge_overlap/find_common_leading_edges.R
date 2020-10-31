#david morse 9/17/19
#1 load data
#2 find LE and LE overlap
#plot LE overlap for sperhoids, organoids, and biopsy and force similar gene ordering


#environment
library(RColorBrewer)
library(lattice)
library(latticeExtra)
library(grid)
library(gridExtra)


#1 load data
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/190916_GSEA_PDX_biopsy2_3")
getwd()
#biopsy
load("LeadingEdge_overlap/LE_biop_data.rda")
#organoid
load("LeadingEdge_overlap/LE_org_data.rda")
#spheroid
load("LeadingEdge_overlap/LE_sphere_data.rda")

#find LE overlap
Common_S_genes <- intersect(S_inter_sphere, S_inter_org)
Common_S_genes <- intersect(Common_S_genes, M_inter_biop)


#plotting heatmaps

# SET HEATMAP PARAMETERS -----------
asp='fill'
threshold = 0.3 #max(ceiling(range(dat1)))-2
scale.range= c(-threshold, threshold)
zak= seq(scale.range[1],scale.range[2],0.01)
nzak=length(zak)
pal = rev(brewer.pal(11,'Spectral'))
pal[6] = 'white'
pal=c('navyblue',pal[1:6],'white',pal[7:11],'#3f003f')
pal = c(rev(c(brewer.pal(9,"Blues"),'navyblue')), c(brewer.pal(9,'Reds'),'#3f003f'))
colo = colorRampPalette(pal, space="Lab")(nzak) 
#breaks = c(seq(-35, 10, length.out=NBR.COLORS/2), 10, seq(10, 35, length.out=NBR.COLORS/2))

#make spheroid heatmap ----------
mat_sphere=t(dd_sphere[Common_S_genes,])
rowv = as.dendrogram(hclust(dist(mat_sphere),method='ward.D'))
colv = as.dendrogram(hclust(dist(t(mat_sphere)),method='ward.D'))
row.ord = order.dendrogram(rowv)
col.ord = order.dendrogram(colv)
sdat = apply(mat_sphere,2,function(x) ifelse(x<=-threshold, -threshold, x))
sdat = t(apply(sdat,2,function(x) ifelse(x>=threshold, threshold, x)))
#rownames(sdat) = substring(rownames(sdat),1,1) #Dave jan 7th - remove to allow full gene names in the rows
colnames(sdat) = c("CENTRAL", "INNER", "OUTER", "SURFACE") #change column names
p=levelplot(sdat[col.ord, row.ord], col.regions=colo, at=zak,  
            #p=levelplot(sdat[col.ord, c("SURFACE", "MIDDLE", "CENTRAL")], col.regions=colo, at=zak, #change order of rows here with c('x','x','x')
            aspect = asp, xlab='', ylab='',
            scales = list(x = list(rot = 90, cex=.5, lwd=0), tck=c(0,0), alternating=2, y = list(cex=0.8)), 
            colorkey = list(space='left', labels=list(cex=0.6)),
            legend = list(
              #right = list(fun = dendrogramGrob, args = list(x = rowv, ord = row.ord, side = "right", size = 1)),
              top = list(fun = dendrogramGrob, args = list(x = colv, ord=col.ord, side = "top", size = 2, colour='white'))))
p.diff = p
# plot comparison heatmap only
grid.arrange(p.diff)

#make organoid heatmap ----------
mat_org=t(dd_org[Common_S_genes,])
#rowv = as.dendrogram(hclust(dist(mat_org),method='ward.D'))
#colv = as.dendrogram(hclust(dist(t(mat_org)),method='ward.D'))
row.ord = order.dendrogram(rowv)
col.ord = order.dendrogram(colv)
sdat = apply(mat_org,2,function(x) ifelse(x<=-threshold, -threshold, x))
sdat = t(apply(sdat,2,function(x) ifelse(x>=threshold, threshold, x)))
#rownames(sdat) = substring(rownames(sdat),1,1) #Dave jan 7th - remove to allow full gene names in the rows
colnames(sdat) = c("CENTRAL", "MIDDLE", "SURFACE") #change column names
#p=levelplot(sdat[col.ord, row.ord], col.regions=colo, at=zak,  
            p=levelplot(sdat[col.ord, c("SURFACE", "MIDDLE", "CENTRAL")], col.regions=colo, at=zak, #change order of rows here with c('x','x','x')
            aspect = asp, xlab='', ylab='',
            scales = list(x = list(rot = 90, cex=.5, lwd=0), tck=c(0,0), alternating=2, y = list(cex=0.8)), 
            colorkey = list(space='left', labels=list(cex=0.6)),
            legend = list(
              #right = list(fun = dendrogramGrob, args = list(x = rowv, ord = row.ord, side = "right", size = 1)),
              top = list(fun = dendrogramGrob, args = list(x = colv, ord=col.ord, side = "top", size = 2, colour='white'))))
p.diff = p
# plot comparison heatmap only
grid.arrange(p.diff)

#make biopsy heatmap ----------
mat_biop=t(dd_biop[Common_S_genes,])
#rowv = as.dendrogram(hclust(dist(mat_biop),method='ward.D'))
#colv = as.dendrogram(hclust(dist(t(mat_biop)),method='ward.D'))
row.ord = order.dendrogram(rowv)
col.ord = order.dendrogram(colv)
sdat = apply(mat_biop,2,function(x) ifelse(x<=-threshold, -threshold, x))
sdat = t(apply(sdat,2,function(x) ifelse(x>=threshold, threshold, x)))
#rownames(sdat) = substring(rownames(sdat),1,1) #Dave jan 7th - remove to allow full gene names in the rows
colnames(sdat) = c("CENTRAL", "MIDDLE", "SURFACE") #change column names
#p=levelplot(sdat[col.ord, row.ord], col.regions=colo, at=zak,  
p=levelplot(sdat[col.ord, c("SURFACE", "MIDDLE", "CENTRAL")], col.regions=colo, at=zak, #change order of rows here with c('x','x','x')
            aspect = asp, xlab='', ylab='',
            scales = list(x = list(rot = 90, cex=.5, lwd=0), tck=c(0,0), alternating=2, y = list(cex=0.8)), 
            colorkey = list(space='left', labels=list(cex=0.6)),
            legend = list(
              #right = list(fun = dendrogramGrob, args = list(x = rowv, ord = row.ord, side = "right", size = 1)),
              top = list(fun = dendrogramGrob, args = list(x = colv, ord=col.ord, side = "top", size = 2, colour='white'))))
p.diff = p
# plot comparison heatmap only
grid.arrange(p.diff)



#vertical HEATMAPS:

#biop ----------
mat_biop=dd_biop[Common_S_genes,]
rowv = as.dendrogram(hclust(dist(mat_biop),method='ward.D'))
colv = as.dendrogram(hclust(dist(t(mat_biop)),method='ward.D'))
row.ord = order.dendrogram(rowv)
col.ord = order.dendrogram(colv)
sdat = apply(mat_biop,2,function(x) ifelse(x<=-threshold, -threshold, x))
sdat = t(apply(sdat,2,function(x) ifelse(x>=threshold, threshold, x)))
#rownames(sdat) = substring(rownames(sdat),1,1)
rownames(sdat) = c("CENTRAL", "MIDDLE", "SURFACE") #change column names
#p=levelplot(sdat[col.ord, row.ord], col.regions=colo, at=zak,  
p=levelplot(sdat[c("CENTRAL", "MIDDLE", "SURFACE"), row.ord], col.regions=colo, at=zak, #change order of rows here with c('x','x','x')
            aspect = asp, xlab='', ylab='',
            scales = list(x = list(rot = 90, cex=3, lwd=0), tck=c(0,0), alternating=2, y = list(cex=0.9)), 
            colorkey = list(space='left'),
            legend = list(
              right = list(fun = dendrogramGrob, args = list(x = rowv, ord = row.ord, side = "right",size = 6)),
              top = list(fun = dendrogramGrob, args = list(x = colv, ord=col.ord, side = "top", size = 2, colour='white'))))
p.diff = p
grid.arrange(p.diff)

pdf("LeadingEdge_overlap/heatmaps/biopsy_heatmap.pdf", width = 4, height = 9)
grid.arrange(p.diff)
dev.off() 


#spheroid ----------
#keep gene order same as biopsy
mat_sphere=dd_sphere[Common_S_genes,]
rowv = as.dendrogram(hclust(dist(mat_sphere),method='ward.D'))
colv = as.dendrogram(hclust(dist(t(mat_sphere)),method='ward.D'))
#row.ord = order.dendrogram(rowv)
#col.ord = order.dendrogram(colv)
sdat = apply(mat_sphere,2,function(x) ifelse(x<=-threshold, -threshold, x))
sdat = t(apply(sdat,2,function(x) ifelse(x>=threshold, threshold, x)))
#rownames(sdat) = substring(rownames(sdat),1,1)
rownames(sdat) = c("CENTRAL", "INNER", "OUTER", "SURFACE") #change column names
#p=levelplot(sdat[col.ord, row.ord], col.regions=colo, at=zak,  
p=levelplot(sdat[c("CENTRAL", "INNER", "OUTER", "SURFACE"), row.ord], col.regions=colo, at=zak, #change order of rows here with c('x','x','x')
            aspect = asp, xlab='', ylab='',
            scales = list(x = list(rot = 90, cex=3, lwd=0), tck=c(0,0), alternating=2, y = list(cex=0.9)), 
            colorkey = list(space='left'),
            legend = list(
              right = list(fun = dendrogramGrob, args = list(x = rowv, ord = row.ord, side = "right",size = 6)),
              top = list(fun = dendrogramGrob, args = list(x = colv, ord=col.ord, side = "top", size = 2, colour='white'))))
p.diff = p
grid.arrange(p.diff)

pdf("LeadingEdge_overlap/heatmaps/spheroid_heatmap.pdf", width = 4.35, height = 9)
grid.arrange(p.diff)
dev.off() 



#organoid ----------
#keep gene order same as biopsy
mat_org=dd_org[Common_S_genes,]
rowv = as.dendrogram(hclust(dist(mat_org),method='ward.D'))
colv = as.dendrogram(hclust(dist(t(mat_org)),method='ward.D'))
#row.ord = order.dendrogram(rowv)
#col.ord = order.dendrogram(colv)
sdat = apply(mat_org,2,function(x) ifelse(x<=-threshold, -threshold, x))
sdat = t(apply(sdat,2,function(x) ifelse(x>=threshold, threshold, x)))
#rownames(sdat) = substring(rownames(sdat),1,1)
rownames(sdat) = c("CENTRAL", "MIDDLE", "SURFACE") #change column names
#p=levelplot(sdat[col.ord, row.ord], col.regions=colo, at=zak,  
p=levelplot(sdat[c("CENTRAL", "MIDDLE", "SURFACE"), row.ord], col.regions=colo, at=zak, #change order of rows here with c('x','x','x')
            aspect = asp, xlab='', ylab='',
            scales = list(x = list(rot = 90, cex=3, lwd=0), tck=c(0,0), alternating=2, y = list(cex=0.9)), 
            colorkey = list(space='left'),
            legend = list(
              right = list(fun = dendrogramGrob, args = list(x = rowv, ord = row.ord, side = "right",size = 6)),
              top = list(fun = dendrogramGrob, args = list(x = colv, ord=col.ord, side = "top", size = 2, colour='white'))))
p.diff = p
grid.arrange(p.diff)

pdf("LeadingEdge_overlap/heatmaps/organoid_heatmap.pdf", width = 4, height = 9)
grid.arrange(p.diff)
dev.off() 



