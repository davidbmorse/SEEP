# load gmt files to R list object

gmtlist <- function(path.gmtFolder, cut.version='.v6.1.symbols.gmt', rda.file="MSigDB/MSigDB_subsets.rda")
  {
  library(fgsea)
  library(cogena)
  collections = list.files(path.gmtFolder, pattern='gmt', full=TRUE)
  nams = sapply(strsplit(collections,path.gmtFolder), function(x) x[2])
  nams = sub("/", "", nams)
  nams = gsub(cut.version, "", nams)
  lgmt = lapply(collections, gmt2list)
  names(lgmt) = nams
  L_subsets = lgmt
  save(L_subsets, file=rda.file)
  print(paste("rda file saved:", rda.file))
}

# estimate rank metric
calc.logMean = function(x) log(mean(expm1(x))+1)

calc.rankMetric = function(dx, groups, out=c('means','difs'))
  {
  calc.logMean = function(x) log(mean(expm1(x))+1)
  mean.1 <- t(apply(dx, 1, tapply, groups, calc.logMean))
  mean.2 <- apply(dx, 1, calc.logMean)
  dif <- mean.1 - mean.2
  if(out=='means'){
    return(data.frame(mean.1, Total=mean.2))
  } else {return(dif)}
}

# compute gsea test

pGSEA <- function(
  rank.metric, names.metric, msigdb="MSigDB/MSigDB_subsets.rda", collections=NULL,
    min.gset=15, max.gset=500, np=10000, random.seed=1234)
{
  library(fgsea)

  load(msigdb)
  if(is.null(collections)){
  L_subsets = L_subsets
  } else {L_subsets = L_subsets[collections]}

  names(rank.metric) = names.metric
  rank.metric = rank.metric[order(rank.metric)]
  
  gsea <- L_subsets
  for(i in 1:length(L_subsets)){
    set.seed(random.seed)
    gsea[[i]] <- fgsea(pathways = L_subsets[[i]],
                     stats = rank.metric,
                     minSize = min.gset,
                     maxSize = max.gset,
                     nperm = np)
  }
 return(gsea)
}

write.gsea <- function(location, gsea.list){
  dir.create(location)
  ltmp=gsea.list
  for(i in 1:length(ltmp)){
    tmp = data.frame(ltmp[[i]])
    tmp$leadingEdge = sapply(tmp$leadingEdge, function(x) paste(x, collapse=', '))
    tmp=tmp[order(tmp$NES, tmp$pval, decreasing=c(TRUE,FALSE)),]
    rownames(tmp)=NULL
    file.name = paste(location,'/', location,"_", names(ltmp)[i], '.csv', sep='')
    write.csv(tmp, file.name)
  } 
}

count.sig <- function(gsea.list, fdr.level, group.name)
{
  gsea <- gsea.list
  tab = NULL
  for(j in 1:length(fdr.level))
    {
    selGSEA <- lapply(gsea, function(x) x[x$padj<=fdr.level[j],])
    selGSEA <- lapply(selGSEA, function(x) x[order(x$pval),])
    for(i in 1:length(selGSEA)){
      if( nrow(selGSEA[[i]])>0 ){
        selGSEA[[i]] <- data.frame(source=names(selGSEA)[i], selGSEA[[i]])
        }
      }
      tabi <- sapply(selGSEA,nrow)
    tab = rbind(tab,tabi)
    rownames(tab)=NULL
  }
  tab = data.frame(Test=group.name, padj=fdr.level,tab)
  return(tab)
}

color.transparent <- function(color, a = 0.5){ #a[0,1]
  color_new <- col2rgb(color)
  apply(color_new, 2,
        function(x){
          rgb(red = x[1], green = x[2], blue = x[3], alpha = a*255, maxColorValue = 255)
        }
  )
}

plot.ES <- function(rank.statistic, gsea, gene.set, name, pathwaySource,location,
                    file.name=NULL, plot.metric=TRUE){
  
  if( !is.null(file.name) ) {
    tiff(paste(file.name,'tiff', sep='.'), units='in', res=300, wid=11, hei=8.5, compress='lzw')
  }
  library(colorspace)
  
  old.par <- par(no.readonly = TRUE)
  if (plot.metric == TRUE){
  nf <- layout(matrix(c(1,1,
                        2,2),
                        2, 2, byrow = TRUE))
  } else { par(mfrow=c(1,1))
                                }
  
  f.yrange <- function(a){
    if(abs(min(a))> abs(max(a))){
      c(floor(min(a)),max(a),-1)
    } else { c(min(a),ceiling(max(a)),1)}
  }
  
  # Enrichment stat plot
  library(fgsea)
  es.data <- plotEnrichment(pathway=gene.set, stats=rank.statistic, gseaParam = 1)
  x <- es.data$data$x
  y <- es.data$data$y
  
  #par()$mar+c(-5,2,0,0) 
  par(mar=c(0.1, 6.1, 4.1, 2.1), oma=c(0,0,1.5,0))
  plot(x, y, axes=F, ylab='', xlab='', type='n', cex.lab=1.5, ylim=f.yrange(y)[1:2],
       xlim = c(0, length(rank.statistic)))
  axis(2, las=2, cex.axis=1.5, xpd=TRUE)
  abline(h=0,lwd=1)
  mtext(paste(name, sep='/'), adj=0.5, cex=1, out=TRUE, line=-0.3, font=2)
  title(ylab="Enrichment Score (ES)", line=3.5, cex.lab=1.5)
  #title(main=location, cex.main=1.5, line=2.3, adj=0.02)
  txt <- paste(location, ', NES = ', round(gsea$NES,2), ', \tQ-value = ', signif(gsea$padj,3), sep='')
  segments(x, -0.02, x, 0.02)
  h <- ifelse(f.yrange(y)[3]==-1, min(y), max(y))
  abline(h=h, col='black',lty='dashed')
  mtext(out=T, line=-3.5, txt, adj=0.27, cex=1.5)
  es.colo <- ifelse(sign(h) == 1, 'red2', 'blue')
  lines(x, y, col=es.colo, lwd=3)
  box()
  
  if(plot.metric==FALSE){
    if(!is.null(file.name)){
      dev.off()
  }} else {
  
  # rank.statistic plot
  par(mar=c(5.1, 6.1, 0.1, 2.1))  
  #par(mar=par()$mar+c(5,0,-4,0))
  plot(rev(sort(rank.statistic)), col='black', type='h', xlab='Rank in Ordered Dataset',
       axes=F, cex.lab=1.5, ylab='')
  axis(1, seq(0,round(max(es.data$data$x)), 1000), cex.axis=1.5)
  axis(2, las=2, cex.axis=1.5)
  box()
  title(ylab="rank statistic", line=3.5, cex.lab=1.5)
  zero = sum(rank.statistic>0)
  abline(v=zero+0.5)
  text(zero, max(rank.statistic)-0.3, paste('zero crossed at', zero+1), adj=-0.05, cex=1.5)
  
  if(!is.null(file.name)){
    dev.off()
  }}
}

getGSEA <- function(pathway, collection = 'h.all', gsea, layer, fdr=0.1, out=c('LE','stat')) #change fdr to 0.2 rather than 0.1
{
  fid = which(gsea[[collection]]$pathway == pathway)
  gseaRes = gsea[[collection]][fid,]
  
  if(out == 'LE'){
    if(gseaRes$padj <= fdr){
      le = sort(unlist(gseaRes$leadingEdge))
  } else { le = NA }
  return(le)} else if(out == 'stat'){
    gseaRes = data.frame(layer, gseaRes[,1:7])
    return(gseaRes)}
}
  
getHeatOrder = function(dx, genes, threshold=0.35){
  
  library(RColorBrewer)
  library(lattice)
  library(latticeExtra)
  
  mat = calc.rankMetric(dx@assays$RNA@data[genes,], dx@active.ident, out='difs')
  scale.range= c(-threshold, threshold)
  zak= seq(scale.range[1],scale.range[2],0.01)
  nzak=length(zak)
  pal = rev(brewer.pal(11,'Spectral'))
  pal[6] = 'white'
  pal=c('navyblue',pal[1:6],'white',pal[7:11],'#3f003f')
  pal = c(rev(c(brewer.pal(9,"Blues"),'navyblue')), c(brewer.pal(9,'Reds'),'#3f003f'))
  colo = colorRampPalette(pal, space="Lab")(nzak) 
  #breaks = c(seq(-35, 10, length.out=NBR.COLORS/2), 10, seq(10, 35, length.out=NBR.COLORS/2))
  rowv = as.dendrogram(hclust(dist(mat),method='ward.D'))
  colv = as.dendrogram(hclust(dist(t(mat)),method='ward.D'))
  row.ord = order.dendrogram(rowv)
  col.ord = order.dendrogram(colv)
  sdat = apply(mat,2,function(x) ifelse(x<=-threshold, -threshold, x))
  sdat = t(apply(sdat,2,function(x) ifelse(x>=threshold, threshold, x)))
  rownames(sdat) = substring(rownames(sdat),1,1)
  p=levelplot(sdat[col.ord, row.ord], col.regions=colo, at=zak,  
              aspect = 'fill', xlab='', ylab='',
              scales = list(x = list(rot = 0, cex=3, lwd=0), tck=c(0,0), alternating=2, y = list(cex=0.9)), 
              colorkey = list(space='left'),
              legend = list(
              right = list(fun = dendrogramGrob, args = list(x = rowv, ord = row.ord, side = "right",size = 6)),
                top = list(fun = dendrogramGrob, args = list(x = colv, ord=col.ord, side = "top", size = 5, colour='white'))))
  
  return(list(gene.order=rowv, plot=p))
}

getHeat= function(dx, genes, ord, layer, threshold=4){
  
  library(RColorBrewer)
  library(lattice)
  library(latticeExtra)
  
  mat = t(apply(dx@assays$RNA@data[genes, dx@active.ident==layer],1, function(x) x - calc.logMean(x)))
  scale.range= c(-threshold, threshold)
  zak= seq(scale.range[1],scale.range[2],0.01)
  nzak=length(zak)
  pal = rev(brewer.pal(11,'Spectral'))
  pal[6] = 'white'
  pal=c('navyblue',pal[1:6],'white',pal[7:11],'#3f003f')
  pal = c(rev(c(brewer.pal(9,"Blues"),'navyblue')), c(brewer.pal(9,'Reds'),'#3f003f'))
  colo = colorRampPalette(pal, space="Lab")(nzak) 
  #breaks = c(seq(-35, 10, length.out=NBR.COLORS/2), 10, seq(10, 35, length.out=NBR.COLORS/2))
  rowv = ord
  colv = as.dendrogram(hclust(dist(t(mat)),method='ward.D'))
  col.ord = order.dendrogram(colv)
  row.ord = order.dendrogram(rowv)
  sdat = apply(mat,2,function(x) ifelse(x<=-threshold, -threshold, x))
  sdat = t(apply(sdat,2,function(x) ifelse(x>=threshold, threshold, x)))
  rownames(sdat) = rep("",nrow(sdat))
  rownames(sdat)[col.ord][ncol(sdat)/2] = layer
  p=levelplot(sdat[col.ord, row.ord], col.regions=colo, at=zak,  
              aspect = 'fill', xlab='', ylab='',
              scales = list(x = list(rot = 0, cex=3, lwd=0), tck=c(0,0), alternating=2, y = list(cex=0.9)), 
              colorkey = FALSE,
              #colorkey = list(space='left'),
              legend = list(
                right = list(fun = dendrogramGrob, args = list(x = rowv, ord = row.ord, side = "right",size = 6)),
                top = list(fun = dendrogramGrob, args = list(x = colv, ord=col.ord, side = "top", size = 5, colour='white'))))
  
  return(p)
}

