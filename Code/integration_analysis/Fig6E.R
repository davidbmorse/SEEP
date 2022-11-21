# settings ====
rm(list=ls())
setwd("C://PROJECTS/P2022/SEEP_Manuscript")

# data
load('data/SEEP_ORA.rda')
load('data/SEEP_Markers.rda')
load("data/SEEP_Regions.rda")

library(fgsea); library(Seurat); library(RColorBrewer); library(tidyverse); library(igraph)

# SEEP integrated cluster colors
colo2=ggsci::pal_d3()(4)
names(colo2) = names(lparent)
colo2 = sub("FF$","BC", colo2)

lpaths = lorasig_main
lgenes = lapply(lpaths, function(x) {
  y = data.frame(x)
  ly = y$overlapGenes
  names(ly) = y$pathway
  return(ly)
})

lparents = lparent
lparent = lapply(names(lpaths), function(x) {
  keep = lpaths[[x]]$pathway
  lparent[[x]] %>% filter(from %in% keep) %>% filter(to %in% keep)
  })
names(lparent) = names(lparents)

# separate nets ====

lnet = sapply(names(lparent), function(x) NULL)
lnodes = ledges = lnet

for(i in 1:length(lparent)){

  edges = lparent[[i]]
  rownames(edges)=NULL
  edges$from = sub("HALLMARK_","",edges$from)
  edges$to = sub("HALLMARK_","",edges$to)
  edges$from = gsub("_","\n",edges$from)
  edges$to = gsub("_","\n",edges$to)
  edges$color = factor(edges$cluster)
  levels(edges$color) = colo2[i]
  edges$color = paste(edges$color)
  list_edges = split(edges, edges$cluster)
  list_nodes = lapply(list_edges, function(x){
    main = unique(x$from[x$from == x$to])
    data.frame(node=unique(c(x$from,x$to)),cluster=unique(x$cluster)) %>% mutate(collapse=ifelse(node %in% main, "parent","child"))
  })
  edges$interaction='pp'

  gedge = lgenes[[i]]
  gedge=lapply(names(gedge), function(x) {
    pp = sub("HALLMARK_", "", x)
    pp = gsub("_", "\n", pp)
    data.frame(from = pp, to = gedge[[x]], cluster = names(lparent)[i], color = unique(edges$color), interaction="pg")
  })
  gedges = data.frame(do.call(rbind, gedge))


  nodes = list_nodes[[1]]
  rownames(nodes) = NULL
  nodes$color = factor(nodes$cluster)
  levels(nodes$color) = colo2[i]
  nodes$color = paste(nodes$color)
  #nodes <- nodes %>% group_by(node) %>% summarize(
    #nclust = length(cluster),
    #cluster=paste(cluster,collapse = ":"),
    #collapse=paste(collapse,collapse=":"),
    #colors = list(ifelse(colo2 %in% color, 1,0)),
    #color = ifelse(length(color)>1, "lightgrey", color)
  #)
  nodes$pathway = paste("HALLMARK", sub("\n","_", nodes$node), sep="_")
  #nodes$size = lorasig_main[[i]]$percent[match(nodes$pathway, lorasig_main[[i]]$pathway)]
  nodes = nodes %>% rename('name'='node', "node_type"="collapse", "origin"="pathway")

  gnodes = gedges %>% select(-from,-interaction) %>% distinct() %>% rename("name"="to" ) %>%
    mutate(node_type='gene', origin=name) %>% select(colnames(nodes))

  nodes = bind_rows(nodes,gnodes)
  edges = bind_rows(edges, gedges)
  lnodes[[i]] = nodes
  ledges[[i]] = edges

  # Generate/format gene set network
  nets = graph_from_data_frame(as.matrix(edges), directed=FALSE, vertices=nodes) %>%  simplify(remove.multiple = FALSE, remove.loops = TRUE)
  V(nets)$degree=degree(nets)
  V(nets)$size = V(nets)$degree+1
  E(nets)$weight=ifelse(E(nets)$cluster %in% c("1:DA_Surface", "2:DA_Surface", "3:DA_Surface",  "4:DA_Center"),0.5,2)
  V(nets)$label = ifelse(V(nets)$node_type=='gene',NA, paste(V(nets)$name))
   lnet[[i]] = nets

}
names(lnet) = c("Surface Cluster 1", "Surface Cluster 2", "Surface Cluster 3", "Center Cluster 1")

# Plot network  with geneset node labels
## tiffs
lgg = lnet
seeds = c(58567, 58567, 3, 1)
tiff("Figures/Fig6E_%01d__11.20.22.tiff", res=300, compress='lzw', units='in', width=6, height=5)
  par(mar=par()$mar-par()$mar)
  for(i in 1:length(lnet)){
  net=lnet[[i]]
  set.seed(seeds[i])
  lay=layout_with_fr(net, dim=2)
  gg = plot(net, vertex.label=V(net)$label, vertex.label.color='black', vertex.frame.width=0, layout=lay, edge.curved=0.2, vertex.label.font=4, asp=1, vertex.size=V(net)$size, vertex.label.cex = V(net)$size/30, vertex.shape=ifelse(V(net)$node_type=="gene","square", "circle"))
  print(gg)
      }
dev.off()

# pdf
lgg = lnet
seeds = c(58567, 58567, 3, 1)
pdf("Figures/Fig6E_11.20.22.pdf", width=6, height=5)
par(mar=par()$mar-par()$mar)
for(i in 1:length(lnet)){
  net=lnet[[i]]
  set.seed(seeds[i])
  lay=layout_with_fr(net, dim=2)
  gg = plot(net, vertex.label=V(net)$label, vertex.label.color='black', vertex.frame.width=0, layout=lay, edge.curved=0.2, vertex.label.font=4, asp=1, vertex.size=V(net)$size, vertex.label.cex = V(net)$size/30, vertex.shape=ifelse(V(net)$node_type=="gene","square", "circle"))
  print(gg)
}
dev.off()


# Plot network  without geneset node labels
## tiffs
lgg = lnet
seeds = c(58567, 58567, 3, 1)
tiff("Figures/Fig6E_NOLABEL_%01d_11.20.22.tiff", res=300, compress='lzw', units='in', width=6, height=5)
par(mar=par()$mar-par()$mar)
for(i in 1:length(lnet)){
  net=lnet[[i]]
  set.seed(seeds[i])
  lay=layout_with_fr(net, dim=2)
  gg = plot(net, vertex.label=NA, vertex.label.color='black', vertex.frame.width=0, layout=lay, edge.curved=0.2, vertex.label.font=4, asp=1, vertex.size=V(net)$size, vertex.label.cex = V(net)$size/30, vertex.shape=ifelse(V(net)$node_type=="gene","square", "circle"))
  print(gg)
}
dev.off()

# pdf
lgg = lnet
seeds = c(58567, 58567, 3, 1)
pdf("Figures/Fig6E_NOLABEL_11.20.22.pdf", width=6, height=5)
par(mar=par()$mar-par()$mar)
for(i in 1:length(lnet)){
  net=lnet[[i]]
  set.seed(seeds[i])
  lay=layout_with_fr(net, dim=2)
  gg = plot(net, vertex.label=NA, vertex.label.color='black', vertex.frame.width=0, layout=lay, edge.curved=0.2, vertex.label.font=4, asp=1, vertex.size=V(net)$size, vertex.label.cex = V(net)$size/30, vertex.shape=ifelse(V(net)$node_type=="gene","square", "circle"))
  print(gg)
}
dev.off()

# suppl ====
lstab = lorasig_main
names(lstab) = c("Suface Cluster 1", "Surface Cluster 2", "Surface Cluster 3", "Center Cluster 1")
lstab = lapply(names(lstab), function(x) data.frame(SEEP_DA_cluster = x, lstab[[x]]))
stab = data.frame(do.call(rbind,lstab)) %>% select(-percent)
stab_ora = stab
stab <- stab %>% mutate(overlapGenes = sapply(overlapGenes, function(x) paste(x, collapse=",")))
write.table(stab, "Figures/Supplement/Model_Cluster_Enrichment.tsv", sep="\t", row=F, quote=F)

save(ledges, lnodes, lnet, lgg, stab_ora, file="data/SEEP_ORAgenenet.rda")
