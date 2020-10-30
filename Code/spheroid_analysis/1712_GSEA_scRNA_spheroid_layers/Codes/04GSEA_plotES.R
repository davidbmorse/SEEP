

gseaRes=central[18,]
nameGset = gseaRes$pathway
collection="h.all"
gset = L_subsets[[collection]][[nameGset]]
ranks = difs$CENTRAL
names(ranks) = rownames(difs)
plot.ES(rank.statistic = ranks, gsea=gseaRes, gene.set=gset, name=nameGset,
        pathwaySource=collection, location='CENTRAL')

