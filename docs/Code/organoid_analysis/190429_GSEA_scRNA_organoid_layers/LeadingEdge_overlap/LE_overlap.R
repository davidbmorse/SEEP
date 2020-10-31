setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190429_GSEA_scRNA_organoid_layers/MSigDB")
load("MSigDB_subsets.rda")
data(pGSEA)
pathway="HALLMARK_INTERFERON_ALPHA_RESPONSE"
loc="CENTRAL"
fid=which(gsea_CENTRAL$h.all$pathway==pathway)
gseaRes=gsea_CENTRAL$h.all[fid, ]
nameGset = gseaRes$pathway
collection="h.all"
gset = L_subsets[[collection]][[nameGset]]
ranks = difs$c
names(ranks) = rownames(difs)
plot.ES(rank.statistic = ranks, gsea=gseaRes, gene.set=gset, name=nameGset, pathwaySource=collection, location=loc)

loc="MIDDLE"
fid=which(gsea_MIDDLE$h.all$pathway==pathway)
gseaRes=gsea_MIDDLE$h.all[fid, ]
nameGset = gseaRes$pathway
collection="h.all"
gset = L_subsets[[collection]][[nameGset]]
ranks = difs$m
names(ranks) = rownames(difs)
M_alpha = sort(unlist(gseaRes$leadingEdge))
plot.ES(rank.statistic = ranks, gsea=gseaRes, gene.set=gset, name=nameGset, pathwaySource=collection, location=loc)

loc="SURFACE"
fid=which(gsea_SURFACE$h.all$pathway==pathway)
gseaRes=gsea_SURFACE$h.all[fid, ]
nameGset = gseaRes$pathway
collection="h.all"
gset = L_subsets[[collection]][[nameGset]]
ranks = difs$s
names(ranks) = rownames(difs)
S_alpha = sort(unlist(gseaRes$leadingEdge))
plot.ES(rank.statistic = ranks, gsea=gseaRes, gene.set=gset, name=nameGset, pathwaySource=collection, location=loc)



## find gene overlab for inf alpha and beta
#define gene sets to examine
pathway_a="HALLMARK_INTERFERON_ALPHA_RESPONSE"
pathway_y="HALLMARK_INTERFERON_GAMMA_RESPONSE"

#center organoid
find_a = which(gsea_CENTRAL$h.all$pathway==pathway_a)
gseaC_a = gsea_CENTRAL$h.all[find_a, ]
find_y = which(gsea_CENTRAL$h.all$pathway==pathway_y)
gseaC_y = gsea_CENTRAL$h.all[find_y, ] 

LE_OrgC_a <- gseaC_a$leadingEdge
LE_OrgC_y <- gseaC_y$leadingEdge
save(LE_OrgC_a, LE_OrgC_y, file = 'LeadingEdge_overlap/LE_OrgC.rda')

#surface organoid
find_a = which(gsea_SURFACE$h.all$pathway==pathway_a)
gseaS_a = gsea_SURFACE$h.all[find_a, ]
find_y = which(gsea_SURFACE$h.all$pathway==pathway_y)
gseaS_y = gsea_SURFACE$h.all[find_y, ] 

LE_OrgS_a <- gseaS_a$leadingEdge
LE_OrgS_y <- gseaS_y$leadingEdge
save(LE_OrgS_a, LE_OrgS_y, file = 'LeadingEdge_overlap/LE_OrgS.rda')

#Spheroid Data (loaded from different directory)---------------
load("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_layers/data/pGSEA.rda")

#center spheroid
find_a = which(gsea_CENTRAL$h.all$pathway==pathway_a)
gseaC_a = gsea_CENTRAL$h.all[find_a, ]
find_y = which(gsea_CENTRAL$h.all$pathway==pathway_y)
gseaC_y = gsea_CENTRAL$h.all[find_y, ] 

LE_SphC_a <- gseaC_a$leadingEdge
LE_SphC_y <- gseaC_y$leadingEdge
save(LE_SphC_a, LE_SphC_y, file = 'LeadingEdge_overlap/LE_SphC.rda')

#surface spheroid
find_a = which(gsea_SURFACE$h.all$pathway==pathway_a)
gseaS_a = gsea_SURFACE$h.all[find_a, ]
find_y = which(gsea_SURFACE$h.all$pathway==pathway_y)
gseaS_y = gsea_SURFACE$h.all[find_y, ] 

LE_SphS_a <- gseaS_a$leadingEdge
LE_SphS_y <- gseaS_y$leadingEdge
save(LE_SphS_a, LE_SphS_y, file = 'LeadingEdge_overlap/LE_SphS.rda')



#find overlap between interferon alpha and gamma genes in spheroid and organoid center

#try again
#load data sets
load("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190429_GSEA_scRNA_organoid_layers/LeadingEdge_overlap/LE_OrgC.rda")
load("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190429_GSEA_scRNA_organoid_layers/LeadingEdge_overlap/LE_OrgS.rda")
load("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190429_GSEA_scRNA_organoid_layers/LeadingEdge_overlap/LE_SphC.rda")
load("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190429_GSEA_scRNA_organoid_layers/LeadingEdge_overlap/LE_SphS.rda")



#make venn diagram
library(VennDiagram)
#combine gene lists for interferone genes from the same layer
Org_C_INFg <- unlist(
  c(split(LE_OrgC_a, ","), split(LE_OrgC_y, ",")))
Org_S_INFg <- unlist(
  c(split(LE_OrgS_a, ","), split(LE_OrgS_y, ",")))
Sph_C_INFg <- unlist(
  c(split(LE_SphC_a, ","), split(LE_SphC_y, ",")))
Sph_S_INFg <- unlist(
  c(split(LE_SphS_a, ","), split(LE_SphS_y, ",")))

Full_C_INFg <- unique(c(Org_C_INFg, Sph_C_INFg))
Full_S_INFg <- unique(c(Org_S_INFg, Sph_S_INFg))


#find common genes
Common_S_genes <- intersect(Org_S_INFg, Sph_S_INFg)
Common_C_genes <- intersect(Org_C_INFg, Sph_C_INFg)

Common_C_alpha_genes <- intersect(unlist(LE_OrgC_a), unlist(LE_SphC_a))
Common_C_gamma_genes <- intersect(unlist(LE_OrgC_y), unlist(LE_SphC_y))
Common_S_alpha_genes <- intersect(unlist(LE_OrgS_a), unlist(LE_SphS_a))
Common_S_gamma_genes <- intersect(unlist(LE_OrgS_y), unlist(LE_SphS_y))

save(Common_S_genes, Common_C_genes, file = 'LeadingEdge_overlap/LE_overlap.rda')


#make some Venn diagrams
venn.plot <- venn.diagram(
  list(Organoid = Org_S_INFg, Spheroid = Sph_S_INFg), 
  "LeadingEdge_overlap/VennPlots/Venn_surface.tiff", main = "Surface",
  main.cex = 3,
  fill = c("cornflowerblue", "darkorchid1"),
  alpha = 0.4,
  label.col = "white",
  cex = 3,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("cornflowerblue", "darkorchid1"),
  cat.cex = 2,
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  cat.dist = c(0.03, 0.03),
  cat.pos = c(-20, 14))

venn.plot <- venn.diagram(
  list(Organoid = Org_C_INFg, Spheroid = Sph_C_INFg), 
  "LeadingEdge_overlap/VennPlots/Venn_center.tiff", main = "Center",
  main.cex = 3,
  rotation.degree = 180,
  fill = c("cornflowerblue", "darkorchid1"),
  alpha = 0.4,
  label.col = "white",
  cex = 3,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("cornflowerblue", "darkorchid1"),
  cat.cex = 2,
  cat.fontfamily = "serif",
  cat.fontface = "bold",
  cat.dist = c(0.03, 0.03),
  cat.pos = c(160, -166))





