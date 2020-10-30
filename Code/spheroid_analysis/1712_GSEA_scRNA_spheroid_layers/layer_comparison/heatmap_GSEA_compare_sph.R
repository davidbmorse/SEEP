ds=5
asp='fill'
library(RColorBrewer)
library(lattice)
library(latticeExtra)
library(grid)
library(gridExtra)

# Average difference heatmap 
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

#mat=dd[inter,] david jan 7th transpose heat map
mat=t(dd[inter,])

rowv = as.dendrogram(hclust(dist(mat),method='ward.D'))
colv = as.dendrogram(hclust(dist(t(mat)),method='ward.D'))
row.ord = order.dendrogram(rowv)
col.ord = order.dendrogram(colv)
sdat = apply(mat,2,function(x) ifelse(x<=-threshold, -threshold, x))
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

