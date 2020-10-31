## input params
ds=2
asp='fill'


library(RColorBrewer)
mat=dd[inter,]

threshold = 0.35 #max(ceiling(range(dat1)))-2
scale.range= c(-threshold, threshold)
zak= seq(scale.range[1],scale.range[2],0.01)
nzak=length(zak)
pal = rev(brewer.pal(11,'RdYlBu'))
pal=c('cyan','skyblue', 'darkblue','black','red', 'yellow', 'white')
#c("#DAF7A6","#FFC300","#FF5733","#C70039","#900C3F","#581845")
#pal[6]='whitesmoke'
colo = colorRampPalette(pal, space="Lab")(nzak) 
#breaks = c(seq(-35, 10, length.out=NBR.COLORS/2), 10, seq(10, 35, length.out=NBR.COLORS/2))


rowv = as.dendrogram(hclust(dist(mat),method='ward.D'))
colv = as.dendrogram(hclust(dist(t(mat)),method='ward.D'))
row.ord = order.dendrogram(rowv)
col.ord = order.dendrogram(colv)
sdat = apply(mat,2,function(x) ifelse(x<=-threshold, -threshold, x))
sdat = t(apply(sdat,2,function(x) ifelse(x>=threshold, threshold, x)))
#sdat = sdat[-1,]
#rownames(sdat) = rep("",nrow(sdat))
rownames(sdat) = substring(rownames(sdat),1,1)

library(lattice)
library(latticeExtra)
# prepare plot
p=levelplot(sdat[col.ord, row.ord], col.regions=colo, at=zak,  
            aspect = asp, xlab='', ylab='',
            scales = list(x = list(rot = 0, cex=1, lwd=0), tck=c(0,0), alternating=2, y = list(cex=1)), 
            colorkey = list(space='left'),
            legend = list(
              right = list(fun = dendrogramGrob, args = list(x = rowv, ord = row.ord, side = "right",size = 6)),
              top = list(fun = dendrogramGrob, args = list(x = colv, ord=col.ord, side = "top", size = ds, colour='white'))
            ))
# assign a new variable name
p.diff = p

# display plot


mat=difmat[,Data@ident=='CENTRAL']

threshold = 4#max(ceiling(range(dat1)))-2
scale.range= c(-threshold, threshold)
zak= seq(scale.range[1],scale.range[2],0.01)
nzak=length(zak)
#pal = rev(brewer.pal(11,'RdYlBu'))
#pal[6]='whitesmoke'
colo = colorRampPalette(pal, space="Lab")(nzak) 
#breaks = c(seq(-35, 10, length.out=NBR.COLORS/2), 10, seq(10, 35, length.out=NBR.COLORS/2))


#rowv = as.dendrogram(hclust(dist(mat),method='ward.D'))
colv = as.dendrogram(hclust(dist(t(mat)),method='ward.D'))
#row.ord = order.dendrogram(rowv)
col.ord = order.dendrogram(colv)
sdat = apply(mat,2,function(x) ifelse(x<=-threshold, -threshold, x))
sdat = t(apply(sdat,2,function(x) ifelse(x>=threshold, threshold, x)))
#sdat = sdat[-1,]
rownames(sdat) = c('CENTRAL',rep("",nrow(sdat)-1))

library(lattice)
library(latticeExtra)
# prepare plot
p=levelplot(sdat[col.ord, row.ord], col.regions=colo, at=zak,  
            aspect = asp, xlab='', ylab='',
            scales = list(x = list(rot = 0, cex=1, lwd=0), tck=c(0,0), alternating=2, y = list(cex=1)), 
            colorkey = list(space='left'),
            legend = list(
              top = list(fun = dendrogramGrob, args = list(x = colv, ord=col.ord, side = "top", size = ds, colour='white'))
            ))
# assign a new variable name
p.c = p



mat=difmat[,Data@ident=='INNER']

threshold = 4#max(ceiling(range(dat1)))-2
scale.range= c(-threshold, threshold)
zak= seq(scale.range[1],scale.range[2],0.01)
nzak=length(zak)
#pal = rev(brewer.pal(11,'RdYlBu'))
#pal[6]='whitesmoke'
colo = colorRampPalette(pal, space="Lab")(nzak) 
#breaks = c(seq(-35, 10, length.out=NBR.COLORS/2), 10, seq(10, 35, length.out=NBR.COLORS/2))


#rowv = as.dendrogram(hclust(dist(mat),method='ward.D'))
colv = as.dendrogram(hclust(dist(t(mat)),method='ward.D'))
#row.ord = order.dendrogram(rowv)
col.ord = order.dendrogram(colv)
sdat = apply(mat,2,function(x) ifelse(x<=-threshold, -threshold, x))
sdat = t(apply(sdat,2,function(x) ifelse(x>=threshold, threshold, x)))
#sdat = sdat[-1,]
rownames(sdat) = c('INNER',rep("",nrow(sdat)-1))


library(lattice)
library(latticeExtra)
# prepare plot
p=levelplot(sdat[col.ord, row.ord], col.regions=colo, at=zak,  
            aspect = asp, xlab='', ylab='',
            scales = list(x = list(rot = 0, cex=1, lwd=0), tck=c(0,0), alternating=2, y = list(cex=1)), 
            colorkey = FALSE,
            legend = list(
            top = list(fun = dendrogramGrob, args = list(x = colv, ord=col.ord, side = "top", size = ds, colour='white'))
            ))
# assign a new variable name
p.i = p




mat=difmat[,Data@ident=='OUTER']

threshold = 4#max(ceiling(range(dat1)))-2
scale.range= c(-threshold, threshold)
zak= seq(scale.range[1],scale.range[2],0.01)
nzak=length(zak)
#pal = rev(brewer.pal(11,'RdYlBu'))
#pal[6]='whitesmoke'
colo = colorRampPalette(pal, space="Lab")(nzak) 
#breaks = c(seq(-35, 10, length.out=NBR.COLORS/2), 10, seq(10, 35, length.out=NBR.COLORS/2))


#rowv = as.dendrogram(hclust(dist(mat),method='ward.D'))
colv = as.dendrogram(hclust(dist(t(mat)),method='ward.D'))
#row.ord = order.dendrogram(rowv)
col.ord = order.dendrogram(colv)
sdat = apply(mat,2,function(x) ifelse(x<=-threshold, -threshold, x))
sdat = t(apply(sdat,2,function(x) ifelse(x>=threshold, threshold, x)))
#sdat = sdat[-1,]
rownames(sdat) = c('OUTER',rep("",nrow(sdat)-1))

library(lattice)
library(latticeExtra)
# prepare plot
p=levelplot(sdat[col.ord, row.ord], col.regions=colo, at=zak,  
            aspect = asp, xlab='', ylab='',
            scales = list(x = list(rot = 0, cex=1, lwd=0), tck=c(0,0), alternating=2, y = list(cex=1)), 
            colorkey = FALSE,
            legend = list(
            top = list(fun = dendrogramGrob, args = list(x = colv, ord=col.ord, side = "top", size = ds, colour='white'))
            ))
# assign a new variable name
p.o = p


mat=difmat[,Data@ident=='SURFACE']

threshold = 4#max(ceiling(range(dat1)))-2
scale.range= c(-threshold, threshold)
zak= seq(scale.range[1],scale.range[2],0.01)
nzak=length(zak)
#pal = rev(brewer.pal(11,'RdYlBu'))
#pal[6]='whitesmoke'
colo = colorRampPalette(pal, space="Lab")(nzak) 
#breaks = c(seq(-35, 10, length.out=NBR.COLORS/2), 10, seq(10, 35, length.out=NBR.COLORS/2))


#rowv = as.dendrogram(hclust(dist(mat),method='ward.D'))
colv = as.dendrogram(hclust(dist(t(mat)),method='ward.D'))
#row.ord = order.dendrogram(rowv)
col.ord = order.dendrogram(colv)
sdat = apply(mat,2,function(x) ifelse(x<=-threshold, -threshold, x))
sdat = t(apply(sdat,2,function(x) ifelse(x>=threshold, threshold, x)))
#sdat = sdat[-1,]
rownames(sdat) = c('SURFACE',rep("",nrow(sdat)-1))

library(lattice)
library(latticeExtra)
# prepare plot
p=levelplot(sdat[col.ord, row.ord], col.regions=colo, at=zak,  
            aspect = asp, xlab='', ylab='',
            colorkey=FALSE,
            scales = list(x = list(rot = 0, cex=1, lwd=0), tck=c(0,0), alternating=2, y = list(cex=1)), 
              legend = list(
                        top = list(fun = dendrogramGrob, args = list(x = colv, ord=col.ord, side = "top", size = ds, colour='white'))
            ))
# assign a new variable name
p.s = p

library(grid)
library(gridExtra)
bmp('Plot_cells.bmp', units='in', res=600, wid=40,hei=20)
p=grid.arrange(p.diff, p.c,p.i,p.o,p.s, ncol=5, widths=c(0.2,1,1,1,1))
p
dev.off()
