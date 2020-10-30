#code to find common gene sets between clusters from same layers of different biology (spheroid/organoid/biopsy here)
#try to find the intersection of the common gene sets and sort by the sum of the normalized enrichment scores across the layers - this will put things enriched in everything first.

#clusters in centers
#load GSEA scores

#subset data table to only NES and gene set names

#find intersection between mutliple data tables like this... https://stackoverflow.com/questions/32917934/how-to-find-common-rows-between-two-dataframe-in-r

library(tidyverse)
library(dplyr)

###centers sph3, org1&6, biop0&3
# cancer gene neighborhoods
#1 load files-----
sph_3 <- read.csv(file="/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/3/3_c4.cgn.csv", header=TRUE, sep=",")
#remove unused columns and rename NES column
sph_3 <- sph_3 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s3 = NES, padj_s3 = padj)
#IF RENAME PRODUCES AN ERROR, WE NEED TO RELOAD PACKAGES - OTHERWISE USE rename(c("old_name" = "new_name"))
org_1 <- read.csv(file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/1/1_c4.cgn.csv", header=TRUE, sep=",")
org_1 <- org_1 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o1 = NES, padj_o1 = padj)

org_6 <- read.csv(file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/6/6_c4.cgn.csv", header=TRUE, sep=",")
org_6 <- org_6 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o6 = NES, padj_o6 = padj)

biop_0 <- read.csv(file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/0/0_c4.cgn.csv", header=TRUE, sep=",")
biop_0 <- biop_0 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b0 = NES, padj_b0 = padj)

biop_3 <- read.csv(file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/3/3_c4.cgn.csv", header=TRUE, sep=",")
biop_3 <- biop_3 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b3 = NES, padj_b3 = padj)

#2 find intersecions ------

#spheroid and organoid1 clusters
C_s3_o1 <- merge(sph_3, org_1, by = "pathway")
C_s3_o1 <- C_s3_o1 %>%
  mutate(NESsum = NES_o1 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())
#now merge with biopsy0 cluster
C_s3_o1_b0 <- merge(C_s3_o1, biop_0, by = "pathway")
C_s3_o1_b0 <- C_s3_o1_b0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o1 + NES_s3 + NES_b0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
#and biopsy3
C_s3_o1_b3 <- merge(C_s3_o1, biop_3, by = "pathway")
C_s3_o1_b3 <- C_s3_o1_b3 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o1 + NES_s3 + NES_b3) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#spheroid and organoid_6
C_s3_o6 <- merge(sph_3, org_6, by = "pathway")
C_s3_o6 <- C_s3_o6 %>%
  mutate(NESsum = NES_o6 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())
#now merge with biopsy0 cluster
C_s3_o6_b0 <- merge(C_s3_o6, biop_0, by = "pathway")
C_s3_o6_b0 <- C_s3_o6_b0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o6 + NES_s3 + NES_b0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
#and biopsy3
C_s3_o6_b3 <- merge(C_s3_o6, biop_3, by = "pathway")
C_s3_o6_b3 <- C_s3_o6_b3 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o6 + NES_s3 + NES_b3) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#finally (for thouroughness)
C_s3_b0 <- merge(sph_3, biop_0, by = "pathway")
C_s3_b0 <- C_s3_b0 %>%
  mutate(NESsum = NES_b0 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_s3_b3 <- merge(sph_3, biop_3, by = "pathway")
C_s3_b3 <- C_s3_b3 %>%
  mutate(NESsum = NES_b3 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

#added later
C_o1_b0 <- merge(org_1, biop_0, by = "pathway")
C_o1_b0 <- C_o1_b0 %>%
  mutate(NESsum = NES_b0 + NES_o1) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_o1_b3 <- merge(org_1, biop_3, by = "pathway")
C_o1_b3 <- C_o1_b3 %>%
  mutate(NESsum = NES_b3 + NES_o1) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_o2_b0 <- merge(org_2, biop_0, by = "pathway")
C_o2_b0 <- C_o2_b0 %>%
  mutate(NESsum = NES_b0 + NES_o2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_o2_b3 <- merge(org_2, biop_3, by = "pathway")
C_o2_b3 <- C_o2_b3 %>%
  mutate(NESsum = NES_b3 + NES_o2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())





###middles  sph 0124, org 02, biop2
#load --------
sph_0 <- read.csv(file="/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/0/0_c4.cgn.csv", header=TRUE, sep=",")
sph_0 <- sph_0 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s0 = NES, padj_s0 = padj)
sph_1 <- read.csv(file="/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/1/1_c4.cgn.csv", header=TRUE, sep=",")
sph_1 <- sph_1 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s1 = NES, padj_s1 = padj)
sph_2 <- read.csv(file="/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/2/2_c4.cgn.csv", header=TRUE, sep=",")
sph_2 <- sph_2 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s2 = NES, padj_s2 = padj)
sph_4 <- read.csv(file="/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/4/4_c4.cgn.csv", header=TRUE, sep=",")
sph_4 <- sph_4 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s4 = NES, padj_s4 = padj)

org_0 <- read.csv(file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/0/0_c4.cgn.csv", header=TRUE, sep=",")
org_0 <- org_0 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o0 = NES, padj_o0 = padj)
org_2 <- read.csv(file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/2/2_c4.cgn.csv", header=TRUE, sep=",")
org_2 <- org_2 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o2 = NES, padj_o2 = padj)

biop_2 <- read.csv(file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/2/2_c4.cgn.csv", header=TRUE, sep=",")
biop_2 <- biop_2 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b2 = NES, padj_b2 = padj)

#2 find intersections ---------

#biopsy2 and organoid0 clusters
M_b2_o0 <- merge(biop_2, org_0, by = "pathway")
M_b2_o0 <- M_b2_o0 %>%
  mutate(NESsum = NES_o0 + NES_b2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())
#biopsy2 and organoid2 clusters
M_b2_o2 <- merge(biop_2, org_2, by = "pathway")
M_b2_o2 <- M_b2_o2 %>%
  mutate(NESsum = NES_o2 + NES_b2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

#merge organoid 0 with spheroid 0,1,2,4 clusters
M_b2_o0_s0 <- merge(M_b2_o0, sph_0, by = "pathway")
M_b2_o0_s0 <- M_b2_o0_s0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o0_s1 <- merge(M_b2_o0, sph_1, by = "pathway")
M_b2_o0_s1 <- M_b2_o0_s1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o0_s2 <- merge(M_b2_o0, sph_2, by = "pathway")
M_b2_o0_s2 <- M_b2_o0_s2 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s2) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o0_s4 <- merge(M_b2_o0, sph_4, by = "pathway")
M_b2_o0_s4 <- M_b2_o0_s4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#merge organoid 2 with spheroid 0,1,2,4 clusters
M_b2_o2_s0 <- merge(M_b2_o2, sph_0, by = "pathway")
M_b2_o2_s0 <- M_b2_o2_s0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o2_s1 <- merge(M_b2_o2, sph_1, by = "pathway")
M_b2_o2_s1 <- M_b2_o2_s1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o2_s2 <- merge(M_b2_o2, sph_2, by = "pathway")
M_b2_o2_s2 <- M_b2_o2_s2 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s2) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o2_s4 <- merge(M_b2_o2, sph_4, by = "pathway")
M_b2_o2_s4 <- M_b2_o2_s4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())


###surfaces  sph 56, org 345, biop 145
#load --------
sph_5 <- read.csv(file="/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/5/5_c4.cgn.csv", header=TRUE, sep=",")
sph_5 <- sph_5 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s5 = NES, padj_s5 = padj)
sph_6 <- read.csv(file="/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/6/6_c4.cgn.csv", header=TRUE, sep=",")
sph_6 <- sph_6 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s6 = NES, padj_s6 = padj)

org_3 <- read.csv(file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/3/3_c4.cgn.csv", header=TRUE, sep=",")
org_3 <- org_3 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o3 = NES, padj_o3 = padj)
org_4 <- read.csv(file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/4/4_c4.cgn.csv", header=TRUE, sep=",")
org_4 <- org_4 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o4 = NES, padj_o4 = padj)
org_5 <- read.csv(file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/5/5_c4.cgn.csv", header=TRUE, sep=",")
org_5 <- org_5 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o5 = NES, padj_o5 = padj)

biop_1 <- read.csv(file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/1/1_c4.cgn.csv", header=TRUE, sep=",")
biop_1 <- biop_1 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b1 = NES, padj_b1 = padj)
biop_4 <- read.csv(file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/4/4_c4.cgn.csv", header=TRUE, sep=",")
biop_4 <- biop_4 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b4 = NES, padj_b4 = padj)
biop_5 <- read.csv(file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/5/5_c4.cgn.csv", header=TRUE, sep=",")
biop_5 <- biop_5 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b5 = NES, padj_b5 = padj)


#2 find intersecions [sph 56, org 345, biop 145] ------
#spheroid5 and organoid clusters
S_s5_o3 <- merge(sph_5, org_3, by = "pathway")
S_s5_o3 <- S_s5_o3 %>%
  mutate(NESsum = NES_o3 + NES_s5) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s5_o4 <- merge(sph_5, org_4, by = "pathway")
S_s5_o4 <- S_s5_o4 %>%
  mutate(NESsum = NES_o4 + NES_s5) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s5_o5 <- merge(sph_5, org_5, by = "pathway")
S_s5_o5 <- S_s5_o5 %>%
  mutate(NESsum = NES_o5 + NES_s5) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

#spheroid6 and organoid clusters
S_s6_o3 <- merge(sph_6, org_3, by = "pathway")
S_s6_o3 <- S_s6_o3 %>%
  mutate(NESsum = NES_o3 + NES_s6) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s6_o4 <- merge(sph_6, org_4, by = "pathway")
S_s6_o4 <- S_s6_o4 %>%
  mutate(NESsum = NES_o4 + NES_s6) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s6_o5 <- merge(sph_6, org_5, by = "pathway")
S_s6_o5 <- S_s6_o5 %>%
  mutate(NESsum = NES_o5 + NES_s6) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())



#now add biopsy 145 layers:
#spheroid 5, org3
S_s5_o3_b1 <- merge(S_s5_o3, biop_1, by = "pathway")
S_s5_o3_b1 <- S_s5_o3_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s5 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o3_b4 <- merge(S_s5_o3, biop_4, by = "pathway")
S_s5_o3_b4 <- S_s5_o3_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s5 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o3_b5 <- merge(S_s5_o3, biop_5, by = "pathway")
S_s5_o3_b5 <- S_s5_o3_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s5 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

##spheroid 5, org4
S_s5_o4_b1 <- merge(S_s5_o4, biop_1, by = "pathway")
S_s5_o4_b1 <- S_s5_o4_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s5 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o4_b4 <- merge(S_s5_o4, biop_4, by = "pathway")
S_s5_o4_b4 <- S_s5_o4_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s5 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o4_b5 <- merge(S_s5_o4, biop_5, by = "pathway")
S_s5_o4_b5 <- S_s5_o4_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s5 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#spheroid 5, org5
S_s5_o5_b1 <- merge(S_s5_o5, biop_1, by = "pathway")
S_s5_o5_b1 <- S_s5_o5_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s5 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o5_b4 <- merge(S_s5_o5, biop_4, by = "pathway")
S_s5_o5_b4 <- S_s5_o5_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s5 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o5_b5 <- merge(S_s5_o5, biop_5, by = "pathway")
S_s5_o5_b5 <- S_s5_o5_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s5 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#now add biopsy 145 layers:
#spheroid 6, org3
S_s6_o3_b1 <- merge(S_s6_o3, biop_1, by = "pathway")
S_s6_o3_b1 <- S_s6_o3_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s6 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o3_b4 <- merge(S_s6_o3, biop_4, by = "pathway")
S_s6_o3_b4 <- S_s6_o3_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s6 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o3_b5 <- merge(S_s6_o3, biop_5, by = "pathway")
S_s6_o3_b5 <- S_s6_o3_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s6 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

##spheroid 6, org4
S_s6_o4_b1 <- merge(S_s6_o4, biop_1, by = "pathway")
S_s6_o4_b1 <- S_s6_o4_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s6 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o4_b4 <- merge(S_s6_o4, biop_4, by = "pathway")
S_s6_o4_b4 <- S_s6_o4_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s6 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o4_b5 <- merge(S_s6_o4, biop_5, by = "pathway")
S_s6_o4_b5 <- S_s6_o4_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s6 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#spheroid 6, org5
S_s6_o5_b1 <- merge(S_s6_o5, biop_1, by = "pathway")
S_s6_o5_b1 <- S_s6_o5_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s6 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o5_b4 <- merge(S_s6_o5, biop_4, by = "pathway")
S_s6_o5_b4 <- S_s6_o5_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s6 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o5_b5 <- merge(S_s6_o5, biop_5, by = "pathway")
S_s6_o5_b5 <- S_s6_o5_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s6 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())




#-------
#spheroid 6, org3
S_s6_o3_b1 <- merge(S_s6_o3, biop_1, by = "pathway")
S_s6_o3_b1 <- S_s6_o3_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s6 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o3_b4 <- merge(S_s6_o3, biop_4, by = "pathway")
S_s6_o3_b4 <- S_s6_o3_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s6 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o3_b5 <- merge(S_s6_o3, biop_5, by = "pathway")
S_s6_o3_b5 <- S_s6_o3_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s6 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

##spheroid 6, org4
S_s6_o4_b1 <- merge(S_s6_o4, biop_1, by = "pathway")
S_s6_o4_b1 <- S_s6_o4_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s6 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o4_b4 <- merge(S_s6_o4, biop_4, by = "pathway")
S_s6_o4_b4 <- S_s6_o4_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s6 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o4_b5 <- merge(S_s6_o4, biop_5, by = "pathway")
S_s6_o4_b5 <- S_s6_o4_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s6 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#spheroid 6, org5
S_s6_o5_b1 <- merge(S_s6_o5, biop_1, by = "pathway")
S_s6_o5_b1 <- S_s6_o5_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s6 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o5_b4 <- merge(S_s6_o5, biop_4, by = "pathway")
S_s6_o5_b4 <- S_s6_o5_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s6 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o5_b5 <- merge(S_s6_o5, biop_5, by = "pathway")
S_s6_o5_b5 <- S_s6_o5_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s6 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())


#Now, save everything ------
#make sure you create some directories too
#central
central <- grep("C_",names(.GlobalEnv),value=TRUE)
central_list <- do.call("list",mget(central))

middle <- grep("M_", names(.GlobalEnv), value = TRUE)
middle_list <- do.call("list", mget(middle))

surface <- grep("S_", names(.GlobalEnv), value = TRUE)
surface_list <- do.call("list", mget(surface))

setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/comparative_analysis_filt_biop/GSEA_enrch_comparison")
getwd()

#save lists to created directories
for(i in names(central_list)){
  write.csv(central_list[[i]], paste0("cancer_gene_neighborhoods/central/",i,".csv"))
}

for(i in names(middle_list)){
  write.csv(middle_list[[i]], paste0("cancer_gene_neighborhoods/middle/",i,".csv"))
}

for(i in names(surface_list)){
  write.csv(surface_list[[i]], paste0("cancer_gene_neighborhoods/surface/",i,".csv"))
}





  










#############################################################################
#Repeat above analysis on different intersteing gene sets (i.e., hallmark, cancer modules, go_biological process c5.bp)

#HALLMARK - filtered ------------
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis")
getwd()
#
sph_3 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/3/3_h.all.csv", header=TRUE, sep=",")
#remove unused columns and rename NES column
sph_3 <- sph_3 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s3 = NES, padj_s3 = padj)

org_1 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/1/1_h.all.csv", header=TRUE, sep=",")
org_1 <- org_1 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o1 = NES, padj_o1 = padj)

org_6 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/6/6_h.all.csv", header=TRUE, sep=",")
org_6 <- org_6 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o6 = NES, padj_o6 = padj)

biop_1 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/1/1_h.all.csv", header=TRUE, sep=",")
biop_1 <- biop_1 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b1 = NES, padj_b1 = padj)

biop_2 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/2/2_h.all.csv", header=TRUE, sep=",")
biop_2 <- biop_2 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b2 = NES, padj_b2 = padj)

sph_0 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/0/0_h.all.csv", header=TRUE, sep=",")
sph_0 <- sph_0 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s0 = NES, padj_s0 = padj)
sph_1 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/1/1_h.all.csv", header=TRUE, sep=",")
sph_1 <- sph_1 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s1 = NES, padj_s1 = padj)
sph_2 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/2/2_h.all.csv", header=TRUE, sep=",")
sph_2 <- sph_2 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s2 = NES, padj_s2 = padj)
sph_4 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/4/4_h.all.csv", header=TRUE, sep=",")
sph_4 <- sph_4 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s4 = NES, padj_s4 = padj)

org_0 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/0/0_h.all.csv", header=TRUE, sep=",")
org_0 <- org_0 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o0 = NES, padj_o0 = padj)
org_2 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/2/2_h.all.csv", header=TRUE, sep=",")
org_2 <- org_2 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o2 = NES, padj_o2 = padj)

biop_4 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/4/4_h.all.csv", header=TRUE, sep=",")
biop_4 <- biop_4 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b4 = NES, padj_b4 = padj)

sph_5 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/5/5_h.all.csv", header=TRUE, sep=",")
sph_5 <- sph_5 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s5 = NES, padj_s5 = padj)
sph_6 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/6/6_h.all.csv", header=TRUE, sep=",")
sph_6 <- sph_6 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s6 = NES, padj_s6 = padj)

org_3 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/3/3_h.all.csv", header=TRUE, sep=",")
org_3 <- org_3 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o3 = NES, padj_o3 = padj)
org_4 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/4/4_h.all.csv", header=TRUE, sep=",")
org_4 <- org_4 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o4 = NES, padj_o4 = padj)
org_5 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/5/5_h.all.csv", header=TRUE, sep=",")
org_5 <- org_5 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o5 = NES, padj_o5 = padj)

biop_0 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/0/0_h.all.csv", header=TRUE, sep=",")
biop_0 <- biop_0 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b0 = NES, padj_b0 = padj)
biop_3 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/3/3_h.all.csv", header=TRUE, sep=",")
biop_3 <- biop_3 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b3 = NES, padj_b3 = padj)
biop_5 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/5/5_h.all.csv", header=TRUE, sep=",")
biop_5 <- biop_5 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b5 = NES, padj_b5 = padj)

#  intersections -------

#spheroid and organoid1 clusters
C_s3_o1 <- merge(sph_3, org_1, by = "pathway")
C_s3_o1 <- C_s3_o1 %>%
  mutate(NESsum = NES_o1 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())
#now merge with biopsy0 cluster
C_s3_o1_b0 <- merge(C_s3_o1, biop_0, by = "pathway")
C_s3_o1_b0 <- C_s3_o1_b0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o1 + NES_s3 + NES_b0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
#and biopsy3
C_s3_o1_b3 <- merge(C_s3_o1, biop_3, by = "pathway")
C_s3_o1_b3 <- C_s3_o1_b3 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o1 + NES_s3 + NES_b3) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#spheroid and organoid_6
C_s3_o6 <- merge(sph_3, org_6, by = "pathway")
C_s3_o6 <- C_s3_o6 %>%
  mutate(NESsum = NES_o6 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())
#now merge with biopsy0 cluster
C_s3_o6_b0 <- merge(C_s3_o6, biop_0, by = "pathway")
C_s3_o6_b0 <- C_s3_o6_b0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o6 + NES_s3 + NES_b0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
#and biopsy3
C_s3_o6_b3 <- merge(C_s3_o6, biop_3, by = "pathway")
C_s3_o6_b3 <- C_s3_o6_b3 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o6 + NES_s3 + NES_b3) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#finally (for thouroughness)
C_s3_b0 <- merge(sph_3, biop_0, by = "pathway")
C_s3_b0 <- C_s3_b0 %>%
  mutate(NESsum = NES_b0 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_s3_b3 <- merge(sph_3, biop_3, by = "pathway")
C_s3_b3 <- C_s3_b3 %>%
  mutate(NESsum = NES_b3 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_o1_b0 <- merge(org_1, biop_0, by = "pathway")
C_o1_b0 <- C_o1_b0 %>%
  mutate(NESsum = NES_b0 + NES_o1) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_o1_b3 <- merge(org_1, biop_3, by = "pathway")
C_o1_b3 <- C_o1_b3 %>%
  mutate(NESsum = NES_b3 + NES_o1) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_o2_b0 <- merge(org_2, biop_0, by = "pathway")
C_o2_b0 <- C_o2_b0 %>%
  mutate(NESsum = NES_b0 + NES_o2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_o2_b3 <- merge(org_2, biop_3, by = "pathway")
C_o2_b3 <- C_o2_b3 %>%
  mutate(NESsum = NES_b3 + NES_o2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())













#MIDDLE
#biopsy2 and organoid0 clusters
M_b2_o0 <- merge(biop_2, org_0, by = "pathway")
M_b2_o0 <- M_b2_o0 %>%
  mutate(NESsum = NES_o0 + NES_b2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())
#biopsy2 and organoid2 clusters
M_b2_o2 <- merge(biop_2, org_2, by = "pathway")
M_b2_o2 <- M_b2_o2 %>%
  mutate(NESsum = NES_o2 + NES_b2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

#merge organoid 0 with spheroid 0,1,2,4 clusters
M_b2_o0_s0 <- merge(M_b2_o0, sph_0, by = "pathway")
M_b2_o0_s0 <- M_b2_o0_s0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o0_s1 <- merge(M_b2_o0, sph_1, by = "pathway")
M_b2_o0_s1 <- M_b2_o0_s1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o0_s2 <- merge(M_b2_o0, sph_2, by = "pathway")
M_b2_o0_s2 <- M_b2_o0_s2 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s2) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o0_s4 <- merge(M_b2_o0, sph_4, by = "pathway")
M_b2_o0_s4 <- M_b2_o0_s4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#merge organoid 2 with spheroid 0,1,2,4 clusters
M_b2_o2_s0 <- merge(M_b2_o2, sph_0, by = "pathway")
M_b2_o2_s0 <- M_b2_o2_s0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o2_s1 <- merge(M_b2_o2, sph_1, by = "pathway")
M_b2_o2_s1 <- M_b2_o2_s1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o2_s2 <- merge(M_b2_o2, sph_2, by = "pathway")
M_b2_o2_s2 <- M_b2_o2_s2 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s2) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o2_s4 <- merge(M_b2_o2, sph_4, by = "pathway")
M_b2_o2_s4 <- M_b2_o2_s4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#biopsy3 and organoid0 clusters
M_b3_o0 <- merge(biop_3, org_0, by = "pathway")
M_b3_o0 <- M_b3_o0 %>%
  mutate(NESsum = NES_o0 + NES_b3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())
#biopsy2 and organoid2 clusters
M_b3_o2 <- merge(biop_3, org_2, by = "pathway")
M_b3_o2 <- M_b3_o2 %>%
  mutate(NESsum = NES_o2 + NES_b3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

#merge organoid 0 with spheroid 0,1,2,4 clusters
M_b3_o0_s0 <- merge(M_b3_o0, sph_0, by = "pathway")
M_b3_o0_s0 <- M_b3_o0_s0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o0 + NES_s0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o0_s1 <- merge(M_b3_o0, sph_1, by = "pathway")
M_b3_o0_s1 <- M_b3_o0_s1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o0 + NES_s1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o0_s2 <- merge(M_b3_o0, sph_2, by = "pathway")
M_b3_o0_s2 <- M_b3_o0_s2 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o0 + NES_s2) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o0_s4 <- merge(M_b3_o0, sph_4, by = "pathway")
M_b3_o0_s4 <- M_b3_o0_s4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o0 + NES_s4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#merge organoid 2 with spheroid 0,1,2,4 clusters
M_b3_o2_s0 <- merge(M_b3_o2, sph_0, by = "pathway")
M_b3_o2_s0 <- M_b3_o2_s0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o2 + NES_s0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o2_s1 <- merge(M_b3_o2, sph_1, by = "pathway")
M_b3_o2_s1 <- M_b3_o2_s1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o2 + NES_s1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o2_s2 <- merge(M_b3_o2, sph_2, by = "pathway")
M_b3_o2_s2 <- M_b3_o2_s2 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o2 + NES_s2) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o2_s4 <- merge(M_b3_o2, sph_4, by = "pathway")
M_b3_o2_s4 <- M_b3_o2_s4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o2 + NES_s4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())





#SURFACE
#spheroid5 and organoid clusters
S_s5_o3 <- merge(sph_5, org_3, by = "pathway")
S_s5_o3 <- S_s5_o3 %>%
  mutate(NESsum = NES_o3 + NES_s5) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s5_o4 <- merge(sph_5, org_4, by = "pathway")
S_s5_o4 <- S_s5_o4 %>%
  mutate(NESsum = NES_o4 + NES_s5) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s5_o5 <- merge(sph_5, org_5, by = "pathway")
S_s5_o5 <- S_s5_o5 %>%
  mutate(NESsum = NES_o5 + NES_s5) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

#spheroid6 and organoid clusters
S_s6_o3 <- merge(sph_6, org_3, by = "pathway")
S_s6_o3 <- S_s6_o3 %>%
  mutate(NESsum = NES_o3 + NES_s6) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s6_o4 <- merge(sph_6, org_4, by = "pathway")
S_s6_o4 <- S_s6_o4 %>%
  mutate(NESsum = NES_o4 + NES_s6) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s6_o5 <- merge(sph_6, org_5, by = "pathway")
S_s6_o5 <- S_s6_o5 %>%
  mutate(NESsum = NES_o5 + NES_s6) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())



#now add biopsy (035) 1,4,5 layers:
#spheroid 5, org3
S_s5_o3_b1 <- merge(S_s5_o3, biop_1, by = "pathway")
S_s5_o3_b1 <- S_s5_o3_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s5 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o3_b4 <- merge(S_s5_o3, biop_4, by = "pathway")
S_s5_o3_b4 <- S_s5_o3_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s5 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o3_b5 <- merge(S_s5_o3, biop_5, by = "pathway")
S_s5_o3_b5 <- S_s5_o3_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s5 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

##spheroid 5, org4
S_s5_o4_b1 <- merge(S_s5_o4, biop_1, by = "pathway")
S_s5_o4_b1 <- S_s5_o4_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s5 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o4_b4 <- merge(S_s5_o4, biop_4, by = "pathway")
S_s5_o4_b4 <- S_s5_o4_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s5 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o4_b5 <- merge(S_s5_o4, biop_5, by = "pathway")
S_s5_o4_b5 <- S_s5_o4_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s5 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#spheroid 5, org5
S_s5_o5_b1 <- merge(S_s5_o5, biop_1, by = "pathway")
S_s5_o5_b1 <- S_s5_o5_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s5 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o5_b4 <- merge(S_s5_o5, biop_4, by = "pathway")
S_s5_o5_b4 <- S_s5_o5_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s5 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o5_b5 <- merge(S_s5_o5, biop_5, by = "pathway")
S_s5_o5_b5 <- S_s5_o5_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s5 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#now add biopsy 035 layers:
#spheroid 6, org3
S_s6_o3_b1 <- merge(S_s6_o3, biop_1, by = "pathway")
S_s6_o3_b1 <- S_s6_o3_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s6 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o3_b4 <- merge(S_s6_o3, biop_4, by = "pathway")
S_s6_o3_b4 <- S_s6_o3_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s6 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o3_b5 <- merge(S_s6_o3, biop_5, by = "pathway")
S_s6_o3_b5 <- S_s6_o3_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s6 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

##spheroid 6, org4
S_s6_o4_b1 <- merge(S_s6_o4, biop_1, by = "pathway")
S_s6_o4_b1 <- S_s6_o4_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s6 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o4_b4 <- merge(S_s6_o4, biop_4, by = "pathway")
S_s6_o4_b4 <- S_s6_o4_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s6 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o4_b5 <- merge(S_s6_o4, biop_5, by = "pathway")
S_s6_o4_b5 <- S_s6_o4_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s6 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#spheroid 6, org5
S_s6_o5_b1 <- merge(S_s6_o5, biop_1, by = "pathway")
S_s6_o5_b1 <- S_s6_o5_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s6 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o5_b4 <- merge(S_s6_o5, biop_4, by = "pathway")
S_s6_o5_b4 <- S_s6_o5_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s6 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o5_b5 <- merge(S_s6_o5, biop_5, by = "pathway")
S_s6_o5_b5 <- S_s6_o5_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s6 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())




#  save ------
#central
central <- grep("C_",names(.GlobalEnv),value=TRUE)
central_list <- do.call("list",mget(central))

middle <- grep("M_", names(.GlobalEnv), value = TRUE)
middle_list <- do.call("list", mget(middle))

surface <- grep("S_", names(.GlobalEnv), value = TRUE)
surface_list <- do.call("list", mget(surface))

setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/comparative_analysis_filt_biop/GSEA_enrch_comparison")
getwd()

#save lists to created directories
for(i in names(central_list)){
  write.csv(central_list[[i]], paste0("hallmark/central/",i,".csv"))
}

for(i in names(middle_list)){
  write.csv(middle_list[[i]], paste0("hallmark/middle/",i,".csv"))
}

for(i in names(surface_list)){
  write.csv(surface_list[[i]], paste0("hallmark/surface/",i,".csv"))
}





#CANCER MODULES - filtered ----------------
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis")
getwd()
#
sph_3 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/3/3_c4.cm.csv", header=TRUE, sep=",")
#remove unused columns and rename NES column
sph_3 <- sph_3 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s3 = NES, padj_s3 = padj)

org_1 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/1/1_c4.cm.csv", header=TRUE, sep=",")
org_1 <- org_1 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o1 = NES, padj_o1 = padj)

org_6 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/6/6_c4.cm.csv", header=TRUE, sep=",")
org_6 <- org_6 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o6 = NES, padj_o6 = padj)

biop_1 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/1/1_c4.cm.csv", header=TRUE, sep=",")
biop_1 <- biop_1 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b1 = NES, padj_b1 = padj)

biop_2 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/2/2_c4.cm.csv", header=TRUE, sep=",")
biop_2 <- biop_2 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b2 = NES, padj_b2 = padj)

sph_0 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/0/0_c4.cm.csv", header=TRUE, sep=",")
sph_0 <- sph_0 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s0 = NES, padj_s0 = padj)
sph_1 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/1/1_c4.cm.csv", header=TRUE, sep=",")
sph_1 <- sph_1 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s1 = NES, padj_s1 = padj)
sph_2 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/2/2_c4.cm.csv", header=TRUE, sep=",")
sph_2 <- sph_2 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s2 = NES, padj_s2 = padj)
sph_4 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/4/4_c4.cm.csv", header=TRUE, sep=",")
sph_4 <- sph_4 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s4 = NES, padj_s4 = padj)

org_0 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/0/0_c4.cm.csv", header=TRUE, sep=",")
org_0 <- org_0 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o0 = NES, padj_o0 = padj)
org_2 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/2/2_c4.cm.csv", header=TRUE, sep=",")
org_2 <- org_2 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o2 = NES, padj_o2 = padj)

biop_4 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/4/4_c4.cm.csv", header=TRUE, sep=",")
biop_4 <- biop_4 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b4 = NES, padj_b4 = padj)

sph_5 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/5/5_c4.cm.csv", header=TRUE, sep=",")
sph_5 <- sph_5 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s5 = NES, padj_s5 = padj)
sph_6 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/6/6_c4.cm.csv", header=TRUE, sep=",")
sph_6 <- sph_6 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s6 = NES, padj_s6 = padj)

org_3 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/3/3_c4.cm.csv", header=TRUE, sep=",")
org_3 <- org_3 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o3 = NES, padj_o3 = padj)
org_4 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/4/4_c4.cm.csv", header=TRUE, sep=",")
org_4 <- org_4 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o4 = NES, padj_o4 = padj)
org_5 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/5/5_c4.cm.csv", header=TRUE, sep=",")
org_5 <- org_5 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o5 = NES, padj_o5 = padj)

biop_0 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/0/0_c4.cm.csv", header=TRUE, sep=",")
biop_0 <- biop_0 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b0 = NES, padj_b0 = padj)
biop_3 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/3/3_c4.cm.csv", header=TRUE, sep=",")
biop_3 <- biop_3 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b3 = NES, padj_b3 = padj)
biop_5 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/5/5_c4.cm.csv", header=TRUE, sep=",")
biop_5 <- biop_5 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b5 = NES, padj_b5 = padj)

#  interactions-------
#spheroid and organoid1 clusters
C_s3_o1 <- merge(sph_3, org_1, by = "pathway")
C_s3_o1 <- C_s3_o1 %>%
  mutate(NESsum = NES_o1 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())
#now merge with biopsy0 cluster
C_s3_o1_b0 <- merge(C_s3_o1, biop_0, by = "pathway")
C_s3_o1_b0 <- C_s3_o1_b0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o1 + NES_s3 + NES_b0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
#and biopsy3
C_s3_o1_b3 <- merge(C_s3_o1, biop_3, by = "pathway")
C_s3_o1_b3 <- C_s3_o1_b3 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o1 + NES_s3 + NES_b3) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#spheroid and organoid_6
C_s3_o6 <- merge(sph_3, org_6, by = "pathway")
C_s3_o6 <- C_s3_o6 %>%
  mutate(NESsum = NES_o6 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())
#now merge with biopsy0 cluster
C_s3_o6_b0 <- merge(C_s3_o6, biop_0, by = "pathway")
C_s3_o6_b0 <- C_s3_o6_b0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o6 + NES_s3 + NES_b0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
#and biopsy3
C_s3_o6_b3 <- merge(C_s3_o6, biop_3, by = "pathway")
C_s3_o6_b3 <- C_s3_o6_b3 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o6 + NES_s3 + NES_b3) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#finally (for thouroughness)
C_s3_b0 <- merge(sph_3, biop_0, by = "pathway")
C_s3_b0 <- C_s3_b0 %>%
  mutate(NESsum = NES_b0 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_s3_b3 <- merge(sph_3, biop_3, by = "pathway")
C_s3_b3 <- C_s3_b3 %>%
  mutate(NESsum = NES_b3 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_o1_b0 <- merge(org_1, biop_0, by = "pathway")
C_o1_b0 <- C_o1_b0 %>%
  mutate(NESsum = NES_b0 + NES_o1) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_o1_b3 <- merge(org_1, biop_3, by = "pathway")
C_o1_b3 <- C_o1_b3 %>%
  mutate(NESsum = NES_b3 + NES_o1) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_o2_b0 <- merge(org_2, biop_0, by = "pathway")
C_o2_b0 <- C_o2_b0 %>%
  mutate(NESsum = NES_b0 + NES_o2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_o2_b3 <- merge(org_2, biop_3, by = "pathway")
C_o2_b3 <- C_o2_b3 %>%
  mutate(NESsum = NES_b3 + NES_o2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())


#MIDDLE
#biopsy2 and organoid0 clusters
M_b2_o0 <- merge(biop_2, org_0, by = "pathway")
M_b2_o0 <- M_b2_o0 %>%
  mutate(NESsum = NES_o0 + NES_b2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())
#biopsy2 and organoid2 clusters
M_b2_o2 <- merge(biop_2, org_2, by = "pathway")
M_b2_o2 <- M_b2_o2 %>%
  mutate(NESsum = NES_o2 + NES_b2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

#merge organoid 0 with spheroid 0,1,2,4 clusters
M_b2_o0_s0 <- merge(M_b2_o0, sph_0, by = "pathway")
M_b2_o0_s0 <- M_b2_o0_s0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o0_s1 <- merge(M_b2_o0, sph_1, by = "pathway")
M_b2_o0_s1 <- M_b2_o0_s1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o0_s2 <- merge(M_b2_o0, sph_2, by = "pathway")
M_b2_o0_s2 <- M_b2_o0_s2 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s2) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o0_s4 <- merge(M_b2_o0, sph_4, by = "pathway")
M_b2_o0_s4 <- M_b2_o0_s4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#merge organoid 2 with spheroid 0,1,2,4 clusters
M_b2_o2_s0 <- merge(M_b2_o2, sph_0, by = "pathway")
M_b2_o2_s0 <- M_b2_o2_s0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o2_s1 <- merge(M_b2_o2, sph_1, by = "pathway")
M_b2_o2_s1 <- M_b2_o2_s1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o2_s2 <- merge(M_b2_o2, sph_2, by = "pathway")
M_b2_o2_s2 <- M_b2_o2_s2 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s2) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o2_s4 <- merge(M_b2_o2, sph_4, by = "pathway")
M_b2_o2_s4 <- M_b2_o2_s4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#biopsy3 and organoid0 clusters
M_b3_o0 <- merge(biop_3, org_0, by = "pathway")
M_b3_o0 <- M_b3_o0 %>%
  mutate(NESsum = NES_o0 + NES_b3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())
#biopsy2 and organoid2 clusters
M_b3_o2 <- merge(biop_3, org_2, by = "pathway")
M_b3_o2 <- M_b3_o2 %>%
  mutate(NESsum = NES_o2 + NES_b3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

#merge organoid 0 with spheroid 0,1,2,4 clusters
M_b3_o0_s0 <- merge(M_b3_o0, sph_0, by = "pathway")
M_b3_o0_s0 <- M_b3_o0_s0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o0 + NES_s0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o0_s1 <- merge(M_b3_o0, sph_1, by = "pathway")
M_b3_o0_s1 <- M_b3_o0_s1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o0 + NES_s1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o0_s2 <- merge(M_b3_o0, sph_2, by = "pathway")
M_b3_o0_s2 <- M_b3_o0_s2 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o0 + NES_s2) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o0_s4 <- merge(M_b3_o0, sph_4, by = "pathway")
M_b3_o0_s4 <- M_b3_o0_s4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o0 + NES_s4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#merge organoid 2 with spheroid 0,1,2,4 clusters
M_b3_o2_s0 <- merge(M_b3_o2, sph_0, by = "pathway")
M_b3_o2_s0 <- M_b3_o2_s0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o2 + NES_s0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o2_s1 <- merge(M_b3_o2, sph_1, by = "pathway")
M_b3_o2_s1 <- M_b3_o2_s1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o2 + NES_s1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o2_s2 <- merge(M_b3_o2, sph_2, by = "pathway")
M_b3_o2_s2 <- M_b3_o2_s2 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o2 + NES_s2) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o2_s4 <- merge(M_b3_o2, sph_4, by = "pathway")
M_b3_o2_s4 <- M_b3_o2_s4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o2 + NES_s4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())





#SURFACE
#spheroid5 and organoid clusters
S_s5_o3 <- merge(sph_5, org_3, by = "pathway")
S_s5_o3 <- S_s5_o3 %>%
  mutate(NESsum = NES_o3 + NES_s5) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s5_o4 <- merge(sph_5, org_4, by = "pathway")
S_s5_o4 <- S_s5_o4 %>%
  mutate(NESsum = NES_o4 + NES_s5) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s5_o5 <- merge(sph_5, org_5, by = "pathway")
S_s5_o5 <- S_s5_o5 %>%
  mutate(NESsum = NES_o5 + NES_s5) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

#spheroid6 and organoid clusters
S_s6_o3 <- merge(sph_6, org_3, by = "pathway")
S_s6_o3 <- S_s6_o3 %>%
  mutate(NESsum = NES_o3 + NES_s6) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s6_o4 <- merge(sph_6, org_4, by = "pathway")
S_s6_o4 <- S_s6_o4 %>%
  mutate(NESsum = NES_o4 + NES_s6) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s6_o5 <- merge(sph_6, org_5, by = "pathway")
S_s6_o5 <- S_s6_o5 %>%
  mutate(NESsum = NES_o5 + NES_s6) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())



#now add biopsy (035) 1,4,5 layers:
#spheroid 5, org3
S_s5_o3_b1 <- merge(S_s5_o3, biop_1, by = "pathway")
S_s5_o3_b1 <- S_s5_o3_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s5 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o3_b4 <- merge(S_s5_o3, biop_4, by = "pathway")
S_s5_o3_b4 <- S_s5_o3_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s5 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o3_b5 <- merge(S_s5_o3, biop_5, by = "pathway")
S_s5_o3_b5 <- S_s5_o3_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s5 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

##spheroid 5, org4
S_s5_o4_b1 <- merge(S_s5_o4, biop_1, by = "pathway")
S_s5_o4_b1 <- S_s5_o4_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s5 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o4_b4 <- merge(S_s5_o4, biop_4, by = "pathway")
S_s5_o4_b4 <- S_s5_o4_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s5 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o4_b5 <- merge(S_s5_o4, biop_5, by = "pathway")
S_s5_o4_b5 <- S_s5_o4_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s5 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#spheroid 5, org5
S_s5_o5_b1 <- merge(S_s5_o5, biop_1, by = "pathway")
S_s5_o5_b1 <- S_s5_o5_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s5 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o5_b4 <- merge(S_s5_o5, biop_4, by = "pathway")
S_s5_o5_b4 <- S_s5_o5_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s5 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o5_b5 <- merge(S_s5_o5, biop_5, by = "pathway")
S_s5_o5_b5 <- S_s5_o5_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s5 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#now add biopsy 035 layers:
#spheroid 6, org3
S_s6_o3_b1 <- merge(S_s6_o3, biop_1, by = "pathway")
S_s6_o3_b1 <- S_s6_o3_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s6 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o3_b4 <- merge(S_s6_o3, biop_4, by = "pathway")
S_s6_o3_b4 <- S_s6_o3_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s6 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o3_b5 <- merge(S_s6_o3, biop_5, by = "pathway")
S_s6_o3_b5 <- S_s6_o3_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s6 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

##spheroid 6, org4
S_s6_o4_b1 <- merge(S_s6_o4, biop_1, by = "pathway")
S_s6_o4_b1 <- S_s6_o4_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s6 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o4_b4 <- merge(S_s6_o4, biop_4, by = "pathway")
S_s6_o4_b4 <- S_s6_o4_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s6 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o4_b5 <- merge(S_s6_o4, biop_5, by = "pathway")
S_s6_o4_b5 <- S_s6_o4_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s6 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#spheroid 6, org5
S_s6_o5_b1 <- merge(S_s6_o5, biop_1, by = "pathway")
S_s6_o5_b1 <- S_s6_o5_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s6 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o5_b4 <- merge(S_s6_o5, biop_4, by = "pathway")
S_s6_o5_b4 <- S_s6_o5_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s6 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o5_b5 <- merge(S_s6_o5, biop_5, by = "pathway")
S_s6_o5_b5 <- S_s6_o5_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s6 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())




#  save ------
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/comparative_analysis_filt_biop/GSEA_enrch_comparison")

central <- grep("C_",names(.GlobalEnv),value=TRUE)
central_list <- do.call("list",mget(central))
middle <- grep("M_", names(.GlobalEnv), value = TRUE)
middle_list <- do.call("list", mget(middle))
surface <- grep("S_", names(.GlobalEnv), value = TRUE)
surface_list <- do.call("list", mget(surface))

#save lists to created directories
for(i in names(central_list)){
  write.csv(central_list[[i]], paste0("cancer_modules/central/",i,".csv"))
}

for(i in names(middle_list)){
  write.csv(middle_list[[i]], paste0("cancer_modules/middle/",i,".csv"))
}

for(i in names(surface_list)){
  write.csv(surface_list[[i]], paste0("cancer_modules/surface/",i,".csv"))
}







#GO BIOLOGICAL PROCESSES - filtered ----------
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis")
getwd()
#
sph_3 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/3/3_c5.bp.csv", header=TRUE, sep=",")
#remove unused columns and rename NES column
sph_3 <- sph_3 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s3 = NES, padj_s3 = padj)

org_1 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/1/1_c5.bp.csv", header=TRUE, sep=",")
org_1 <- org_1 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o1 = NES, padj_o1 = padj)

org_6 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/6/6_c5.bp.csv", header=TRUE, sep=",")
org_6 <- org_6 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o6 = NES, padj_o6 = padj)

biop_1 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/1/1_c5.bp.csv", header=TRUE, sep=",")
biop_1 <- biop_1 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b1 = NES, padj_b1 = padj)

biop_2 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/2/2_c5.bp.csv", header=TRUE, sep=",")
biop_2 <- biop_2 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b2 = NES, padj_b2 = padj)

sph_0 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/0/0_c5.bp.csv", header=TRUE, sep=",")
sph_0 <- sph_0 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s0 = NES, padj_s0 = padj)
sph_1 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/1/1_c5.bp.csv", header=TRUE, sep=",")
sph_1 <- sph_1 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s1 = NES, padj_s1 = padj)
sph_2 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/2/2_c5.bp.csv", header=TRUE, sep=",")
sph_2 <- sph_2 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s2 = NES, padj_s2 = padj)
sph_4 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/4/4_c5.bp.csv", header=TRUE, sep=",")
sph_4 <- sph_4 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s4 = NES, padj_s4 = padj)

org_0 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/0/0_c5.bp.csv", header=TRUE, sep=",")
org_0 <- org_0 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o0 = NES, padj_o0 = padj)
org_2 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/2/2_c5.bp.csv", header=TRUE, sep=",")
org_2 <- org_2 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o2 = NES, padj_o2 = padj)

biop_4 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/4/4_c5.bp.csv", header=TRUE, sep=",")
biop_4 <- biop_4 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b4 = NES, padj_b4 = padj)

sph_5 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/5/5_c5.bp.csv", header=TRUE, sep=",")
sph_5 <- sph_5 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s5 = NES, padj_s5 = padj)
sph_6 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/6/6_c5.bp.csv", header=TRUE, sep=",")
sph_6 <- sph_6 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s6 = NES, padj_s6 = padj)

org_3 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/3/3_c5.bp.csv", header=TRUE, sep=",")
org_3 <- org_3 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o3 = NES, padj_o3 = padj)
org_4 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/4/4_c5.bp.csv", header=TRUE, sep=",")
org_4 <- org_4 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o4 = NES, padj_o4 = padj)
org_5 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/5/5_c5.bp.csv", header=TRUE, sep=",")
org_5 <- org_5 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o5 = NES, padj_o5 = padj)

biop_0 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/0/0_c5.bp.csv", header=TRUE, sep=",")
biop_0 <- biop_0 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b0 = NES, padj_b0 = padj)
biop_3 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/3/3_c5.bp.csv", header=TRUE, sep=",")
biop_3 <- biop_3 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b3 = NES, padj_b3 = padj)
biop_5 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/5/5_c5.bp.csv", header=TRUE, sep=",")
biop_5 <- biop_5 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b5 = NES, padj_b5 = padj)

#  interactions-------
#spheroid and organoid1 clusters
C_s3_o1 <- merge(sph_3, org_1, by = "pathway")
C_s3_o1 <- C_s3_o1 %>%
  mutate(NESsum = NES_o1 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())
#now merge with biopsy0 cluster
C_s3_o1_b0 <- merge(C_s3_o1, biop_0, by = "pathway")
C_s3_o1_b0 <- C_s3_o1_b0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o1 + NES_s3 + NES_b0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
#and biopsy3
C_s3_o1_b3 <- merge(C_s3_o1, biop_3, by = "pathway")
C_s3_o1_b3 <- C_s3_o1_b3 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o1 + NES_s3 + NES_b3) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#spheroid and organoid_6
C_s3_o6 <- merge(sph_3, org_6, by = "pathway")
C_s3_o6 <- C_s3_o6 %>%
  mutate(NESsum = NES_o6 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())
#now merge with biopsy0 cluster
C_s3_o6_b0 <- merge(C_s3_o6, biop_0, by = "pathway")
C_s3_o6_b0 <- C_s3_o6_b0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o6 + NES_s3 + NES_b0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
#and biopsy3
C_s3_o6_b3 <- merge(C_s3_o6, biop_3, by = "pathway")
C_s3_o6_b3 <- C_s3_o6_b3 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o6 + NES_s3 + NES_b3) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#finally (for thouroughness)
C_s3_b0 <- merge(sph_3, biop_0, by = "pathway")
C_s3_b0 <- C_s3_b0 %>%
  mutate(NESsum = NES_b0 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_s3_b3 <- merge(sph_3, biop_3, by = "pathway")
C_s3_b3 <- C_s3_b3 %>%
  mutate(NESsum = NES_b3 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_o1_b0 <- merge(org_1, biop_0, by = "pathway")
C_o1_b0 <- C_o1_b0 %>%
  mutate(NESsum = NES_b0 + NES_o1) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_o1_b3 <- merge(org_1, biop_3, by = "pathway")
C_o1_b3 <- C_o1_b3 %>%
  mutate(NESsum = NES_b3 + NES_o1) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_o2_b0 <- merge(org_2, biop_0, by = "pathway")
C_o2_b0 <- C_o2_b0 %>%
  mutate(NESsum = NES_b0 + NES_o2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_o2_b3 <- merge(org_2, biop_3, by = "pathway")
C_o2_b3 <- C_o2_b3 %>%
  mutate(NESsum = NES_b3 + NES_o2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())


#MIDDLE
#biopsy2 and organoid0 clusters
M_b2_o0 <- merge(biop_2, org_0, by = "pathway")
M_b2_o0 <- M_b2_o0 %>%
  mutate(NESsum = NES_o0 + NES_b2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())
#biopsy2 and organoid2 clusters
M_b2_o2 <- merge(biop_2, org_2, by = "pathway")
M_b2_o2 <- M_b2_o2 %>%
  mutate(NESsum = NES_o2 + NES_b2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

#merge organoid 0 with spheroid 0,1,2,4 clusters
M_b2_o0_s0 <- merge(M_b2_o0, sph_0, by = "pathway")
M_b2_o0_s0 <- M_b2_o0_s0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o0_s1 <- merge(M_b2_o0, sph_1, by = "pathway")
M_b2_o0_s1 <- M_b2_o0_s1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o0_s2 <- merge(M_b2_o0, sph_2, by = "pathway")
M_b2_o0_s2 <- M_b2_o0_s2 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s2) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o0_s4 <- merge(M_b2_o0, sph_4, by = "pathway")
M_b2_o0_s4 <- M_b2_o0_s4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#merge organoid 2 with spheroid 0,1,2,4 clusters
M_b2_o2_s0 <- merge(M_b2_o2, sph_0, by = "pathway")
M_b2_o2_s0 <- M_b2_o2_s0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o2_s1 <- merge(M_b2_o2, sph_1, by = "pathway")
M_b2_o2_s1 <- M_b2_o2_s1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o2_s2 <- merge(M_b2_o2, sph_2, by = "pathway")
M_b2_o2_s2 <- M_b2_o2_s2 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s2) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o2_s4 <- merge(M_b2_o2, sph_4, by = "pathway")
M_b2_o2_s4 <- M_b2_o2_s4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#biopsy3 and organoid0 clusters
M_b3_o0 <- merge(biop_3, org_0, by = "pathway")
M_b3_o0 <- M_b3_o0 %>%
  mutate(NESsum = NES_o0 + NES_b3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())
#biopsy2 and organoid2 clusters
M_b3_o2 <- merge(biop_3, org_2, by = "pathway")
M_b3_o2 <- M_b3_o2 %>%
  mutate(NESsum = NES_o2 + NES_b3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

#merge organoid 0 with spheroid 0,1,2,4 clusters
M_b3_o0_s0 <- merge(M_b3_o0, sph_0, by = "pathway")
M_b3_o0_s0 <- M_b3_o0_s0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o0 + NES_s0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o0_s1 <- merge(M_b3_o0, sph_1, by = "pathway")
M_b3_o0_s1 <- M_b3_o0_s1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o0 + NES_s1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o0_s2 <- merge(M_b3_o0, sph_2, by = "pathway")
M_b3_o0_s2 <- M_b3_o0_s2 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o0 + NES_s2) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o0_s4 <- merge(M_b3_o0, sph_4, by = "pathway")
M_b3_o0_s4 <- M_b3_o0_s4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o0 + NES_s4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#merge organoid 2 with spheroid 0,1,2,4 clusters
M_b3_o2_s0 <- merge(M_b3_o2, sph_0, by = "pathway")
M_b3_o2_s0 <- M_b3_o2_s0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o2 + NES_s0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o2_s1 <- merge(M_b3_o2, sph_1, by = "pathway")
M_b3_o2_s1 <- M_b3_o2_s1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o2 + NES_s1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o2_s2 <- merge(M_b3_o2, sph_2, by = "pathway")
M_b3_o2_s2 <- M_b3_o2_s2 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o2 + NES_s2) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o2_s4 <- merge(M_b3_o2, sph_4, by = "pathway")
M_b3_o2_s4 <- M_b3_o2_s4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o2 + NES_s4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())





#SURFACE
#spheroid5 and organoid clusters
S_s5_o3 <- merge(sph_5, org_3, by = "pathway")
S_s5_o3 <- S_s5_o3 %>%
  mutate(NESsum = NES_o3 + NES_s5) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s5_o4 <- merge(sph_5, org_4, by = "pathway")
S_s5_o4 <- S_s5_o4 %>%
  mutate(NESsum = NES_o4 + NES_s5) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s5_o5 <- merge(sph_5, org_5, by = "pathway")
S_s5_o5 <- S_s5_o5 %>%
  mutate(NESsum = NES_o5 + NES_s5) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

#spheroid6 and organoid clusters
S_s6_o3 <- merge(sph_6, org_3, by = "pathway")
S_s6_o3 <- S_s6_o3 %>%
  mutate(NESsum = NES_o3 + NES_s6) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s6_o4 <- merge(sph_6, org_4, by = "pathway")
S_s6_o4 <- S_s6_o4 %>%
  mutate(NESsum = NES_o4 + NES_s6) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s6_o5 <- merge(sph_6, org_5, by = "pathway")
S_s6_o5 <- S_s6_o5 %>%
  mutate(NESsum = NES_o5 + NES_s6) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())



#now add biopsy (035) 1,4,5 layers:
#spheroid 5, org3
S_s5_o3_b1 <- merge(S_s5_o3, biop_1, by = "pathway")
S_s5_o3_b1 <- S_s5_o3_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s5 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o3_b4 <- merge(S_s5_o3, biop_4, by = "pathway")
S_s5_o3_b4 <- S_s5_o3_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s5 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o3_b5 <- merge(S_s5_o3, biop_5, by = "pathway")
S_s5_o3_b5 <- S_s5_o3_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s5 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

##spheroid 5, org4
S_s5_o4_b1 <- merge(S_s5_o4, biop_1, by = "pathway")
S_s5_o4_b1 <- S_s5_o4_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s5 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o4_b4 <- merge(S_s5_o4, biop_4, by = "pathway")
S_s5_o4_b4 <- S_s5_o4_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s5 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o4_b5 <- merge(S_s5_o4, biop_5, by = "pathway")
S_s5_o4_b5 <- S_s5_o4_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s5 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#spheroid 5, org5
S_s5_o5_b1 <- merge(S_s5_o5, biop_1, by = "pathway")
S_s5_o5_b1 <- S_s5_o5_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s5 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o5_b4 <- merge(S_s5_o5, biop_4, by = "pathway")
S_s5_o5_b4 <- S_s5_o5_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s5 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o5_b5 <- merge(S_s5_o5, biop_5, by = "pathway")
S_s5_o5_b5 <- S_s5_o5_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s5 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#now add biopsy 035 layers:
#spheroid 6, org3
S_s6_o3_b1 <- merge(S_s6_o3, biop_1, by = "pathway")
S_s6_o3_b1 <- S_s6_o3_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s6 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o3_b4 <- merge(S_s6_o3, biop_4, by = "pathway")
S_s6_o3_b4 <- S_s6_o3_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s6 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o3_b5 <- merge(S_s6_o3, biop_5, by = "pathway")
S_s6_o3_b5 <- S_s6_o3_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s6 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

##spheroid 6, org4
S_s6_o4_b1 <- merge(S_s6_o4, biop_1, by = "pathway")
S_s6_o4_b1 <- S_s6_o4_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s6 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o4_b4 <- merge(S_s6_o4, biop_4, by = "pathway")
S_s6_o4_b4 <- S_s6_o4_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s6 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o4_b5 <- merge(S_s6_o4, biop_5, by = "pathway")
S_s6_o4_b5 <- S_s6_o4_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s6 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#spheroid 6, org5
S_s6_o5_b1 <- merge(S_s6_o5, biop_1, by = "pathway")
S_s6_o5_b1 <- S_s6_o5_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s6 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o5_b4 <- merge(S_s6_o5, biop_4, by = "pathway")
S_s6_o5_b4 <- S_s6_o5_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s6 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o5_b5 <- merge(S_s6_o5, biop_5, by = "pathway")
S_s6_o5_b5 <- S_s6_o5_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s6 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())




#  save ------
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/comparative_analysis_filt_biop/GSEA_enrch_comparison")

central <- grep("C_",names(.GlobalEnv),value=TRUE)
central_list <- do.call("list",mget(central))
middle <- grep("M_", names(.GlobalEnv), value = TRUE)
middle_list <- do.call("list", mget(middle))
surface <- grep("S_", names(.GlobalEnv), value = TRUE)
surface_list <- do.call("list", mget(surface))

#save lists to created directories
for(i in names(central_list)){
  write.csv(central_list[[i]], paste0("go_biological_processes/central/",i,".csv"))
}

for(i in names(middle_list)){
  write.csv(middle_list[[i]], paste0("go_biological_processes/middle/",i,".csv"))
}

for(i in names(surface_list)){
  write.csv(surface_list[[i]], paste0("go_biological_processes/surface/",i,".csv"))
}







#CANONICAL PATHWAYS - filtered -----
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis")
getwd()
#
sph_3 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/3/3_c2.cp.csv", header=TRUE, sep=",")
#remove unused columns and rename NES column
sph_3 <- sph_3 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s3 = NES, padj_s3 = padj)

org_1 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/1/1_c2.cp.csv", header=TRUE, sep=",")
org_1 <- org_1 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o1 = NES, padj_o1 = padj)

org_6 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/6/6_c2.cp.csv", header=TRUE, sep=",")
org_6 <- org_6 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o6 = NES, padj_o6 = padj)

biop_1 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/1/1_c2.cp.csv", header=TRUE, sep=",")
biop_1 <- biop_1 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b1 = NES, padj_b1 = padj)

biop_2 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/2/2_c2.cp.csv", header=TRUE, sep=",")
biop_2 <- biop_2 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b2 = NES, padj_b2 = padj)

sph_0 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/0/0_c2.cp.csv", header=TRUE, sep=",")
sph_0 <- sph_0 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s0 = NES, padj_s0 = padj)
sph_1 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/1/1_c2.cp.csv", header=TRUE, sep=",")
sph_1 <- sph_1 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s1 = NES, padj_s1 = padj)
sph_2 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/2/2_c2.cp.csv", header=TRUE, sep=",")
sph_2 <- sph_2 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s2 = NES, padj_s2 = padj)
sph_4 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/4/4_c2.cp.csv", header=TRUE, sep=",")
sph_4 <- sph_4 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s4 = NES, padj_s4 = padj)

org_0 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/0/0_c2.cp.csv", header=TRUE, sep=",")
org_0 <- org_0 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o0 = NES, padj_o0 = padj)
org_2 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/2/2_c2.cp.csv", header=TRUE, sep=",")
org_2 <- org_2 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o2 = NES, padj_o2 = padj)

biop_4 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/4/4_c2.cp.csv", header=TRUE, sep=",")
biop_4 <- biop_4 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b4 = NES, padj_b4 = padj)

sph_5 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/5/5_c2.cp.csv", header=TRUE, sep=",")
sph_5 <- sph_5 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s5 = NES, padj_s5 = padj)
sph_6 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/6/6_c2.cp.csv", header=TRUE, sep=",")
sph_6 <- sph_6 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s6 = NES, padj_s6 = padj)

org_3 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/3/3_c2.cp.csv", header=TRUE, sep=",")
org_3 <- org_3 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o3 = NES, padj_o3 = padj)
org_4 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/4/4_c2.cp.csv", header=TRUE, sep=",")
org_4 <- org_4 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o4 = NES, padj_o4 = padj)
org_5 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/5/5_c2.cp.csv", header=TRUE, sep=",")
org_5 <- org_5 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o5 = NES, padj_o5 = padj)

biop_0 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/0/0_c2.cp.csv", header=TRUE, sep=",")
biop_0 <- biop_0 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b0 = NES, padj_b0 = padj)
biop_3 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/3/3_c2.cp.csv", header=TRUE, sep=",")
biop_3 <- biop_3 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b3 = NES, padj_b3 = padj)
biop_5 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/5/5_c2.cp.csv", header=TRUE, sep=",")
biop_5 <- biop_5 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b5 = NES, padj_b5 = padj)

#  interactions-------
#spheroid and organoid1 clusters
C_s3_o1 <- merge(sph_3, org_1, by = "pathway")
C_s3_o1 <- C_s3_o1 %>%
  mutate(NESsum = NES_o1 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())
#now merge with biopsy0 cluster
C_s3_o1_b0 <- merge(C_s3_o1, biop_0, by = "pathway")
C_s3_o1_b0 <- C_s3_o1_b0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o1 + NES_s3 + NES_b0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
#and biopsy3
C_s3_o1_b3 <- merge(C_s3_o1, biop_3, by = "pathway")
C_s3_o1_b3 <- C_s3_o1_b3 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o1 + NES_s3 + NES_b3) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#spheroid and organoid_6
C_s3_o6 <- merge(sph_3, org_6, by = "pathway")
C_s3_o6 <- C_s3_o6 %>%
  mutate(NESsum = NES_o6 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())
#now merge with biopsy0 cluster
C_s3_o6_b0 <- merge(C_s3_o6, biop_0, by = "pathway")
C_s3_o6_b0 <- C_s3_o6_b0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o6 + NES_s3 + NES_b0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
#and biopsy3
C_s3_o6_b3 <- merge(C_s3_o6, biop_3, by = "pathway")
C_s3_o6_b3 <- C_s3_o6_b3 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o6 + NES_s3 + NES_b3) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#finally (for thouroughness)
C_s3_b0 <- merge(sph_3, biop_0, by = "pathway")
C_s3_b0 <- C_s3_b0 %>%
  mutate(NESsum = NES_b0 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_s3_b3 <- merge(sph_3, biop_3, by = "pathway")
C_s3_b3 <- C_s3_b3 %>%
  mutate(NESsum = NES_b3 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_o1_b0 <- merge(org_1, biop_0, by = "pathway")
C_o1_b0 <- C_o1_b0 %>%
  mutate(NESsum = NES_b0 + NES_o1) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_o1_b3 <- merge(org_1, biop_3, by = "pathway")
C_o1_b3 <- C_o1_b3 %>%
  mutate(NESsum = NES_b3 + NES_o1) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_o2_b0 <- merge(org_2, biop_0, by = "pathway")
C_o2_b0 <- C_o2_b0 %>%
  mutate(NESsum = NES_b0 + NES_o2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_o2_b3 <- merge(org_2, biop_3, by = "pathway")
C_o2_b3 <- C_o2_b3 %>%
  mutate(NESsum = NES_b3 + NES_o2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())


#MIDDLE
#biopsy2 and organoid0 clusters
M_b2_o0 <- merge(biop_2, org_0, by = "pathway")
M_b2_o0 <- M_b2_o0 %>%
  mutate(NESsum = NES_o0 + NES_b2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())
#biopsy2 and organoid2 clusters
M_b2_o2 <- merge(biop_2, org_2, by = "pathway")
M_b2_o2 <- M_b2_o2 %>%
  mutate(NESsum = NES_o2 + NES_b2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

#merge organoid 0 with spheroid 0,1,2,4 clusters
M_b2_o0_s0 <- merge(M_b2_o0, sph_0, by = "pathway")
M_b2_o0_s0 <- M_b2_o0_s0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o0_s1 <- merge(M_b2_o0, sph_1, by = "pathway")
M_b2_o0_s1 <- M_b2_o0_s1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o0_s2 <- merge(M_b2_o0, sph_2, by = "pathway")
M_b2_o0_s2 <- M_b2_o0_s2 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s2) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o0_s4 <- merge(M_b2_o0, sph_4, by = "pathway")
M_b2_o0_s4 <- M_b2_o0_s4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#merge organoid 2 with spheroid 0,1,2,4 clusters
M_b2_o2_s0 <- merge(M_b2_o2, sph_0, by = "pathway")
M_b2_o2_s0 <- M_b2_o2_s0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o2_s1 <- merge(M_b2_o2, sph_1, by = "pathway")
M_b2_o2_s1 <- M_b2_o2_s1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o2_s2 <- merge(M_b2_o2, sph_2, by = "pathway")
M_b2_o2_s2 <- M_b2_o2_s2 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s2) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o2_s4 <- merge(M_b2_o2, sph_4, by = "pathway")
M_b2_o2_s4 <- M_b2_o2_s4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#biopsy3 and organoid0 clusters
M_b3_o0 <- merge(biop_3, org_0, by = "pathway")
M_b3_o0 <- M_b3_o0 %>%
  mutate(NESsum = NES_o0 + NES_b3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())
#biopsy2 and organoid2 clusters
M_b3_o2 <- merge(biop_3, org_2, by = "pathway")
M_b3_o2 <- M_b3_o2 %>%
  mutate(NESsum = NES_o2 + NES_b3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

#merge organoid 0 with spheroid 0,1,2,4 clusters
M_b3_o0_s0 <- merge(M_b3_o0, sph_0, by = "pathway")
M_b3_o0_s0 <- M_b3_o0_s0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o0 + NES_s0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o0_s1 <- merge(M_b3_o0, sph_1, by = "pathway")
M_b3_o0_s1 <- M_b3_o0_s1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o0 + NES_s1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o0_s2 <- merge(M_b3_o0, sph_2, by = "pathway")
M_b3_o0_s2 <- M_b3_o0_s2 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o0 + NES_s2) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o0_s4 <- merge(M_b3_o0, sph_4, by = "pathway")
M_b3_o0_s4 <- M_b3_o0_s4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o0 + NES_s4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#merge organoid 2 with spheroid 0,1,2,4 clusters
M_b3_o2_s0 <- merge(M_b3_o2, sph_0, by = "pathway")
M_b3_o2_s0 <- M_b3_o2_s0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o2 + NES_s0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o2_s1 <- merge(M_b3_o2, sph_1, by = "pathway")
M_b3_o2_s1 <- M_b3_o2_s1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o2 + NES_s1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o2_s2 <- merge(M_b3_o2, sph_2, by = "pathway")
M_b3_o2_s2 <- M_b3_o2_s2 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o2 + NES_s2) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o2_s4 <- merge(M_b3_o2, sph_4, by = "pathway")
M_b3_o2_s4 <- M_b3_o2_s4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o2 + NES_s4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())





#SURFACE
#spheroid5 and organoid clusters
S_s5_o3 <- merge(sph_5, org_3, by = "pathway")
S_s5_o3 <- S_s5_o3 %>%
  mutate(NESsum = NES_o3 + NES_s5) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s5_o4 <- merge(sph_5, org_4, by = "pathway")
S_s5_o4 <- S_s5_o4 %>%
  mutate(NESsum = NES_o4 + NES_s5) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s5_o5 <- merge(sph_5, org_5, by = "pathway")
S_s5_o5 <- S_s5_o5 %>%
  mutate(NESsum = NES_o5 + NES_s5) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

#spheroid6 and organoid clusters
S_s6_o3 <- merge(sph_6, org_3, by = "pathway")
S_s6_o3 <- S_s6_o3 %>%
  mutate(NESsum = NES_o3 + NES_s6) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s6_o4 <- merge(sph_6, org_4, by = "pathway")
S_s6_o4 <- S_s6_o4 %>%
  mutate(NESsum = NES_o4 + NES_s6) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s6_o5 <- merge(sph_6, org_5, by = "pathway")
S_s6_o5 <- S_s6_o5 %>%
  mutate(NESsum = NES_o5 + NES_s6) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())



#now add biopsy (035) 1,4,5 layers:
#spheroid 5, org3
S_s5_o3_b1 <- merge(S_s5_o3, biop_1, by = "pathway")
S_s5_o3_b1 <- S_s5_o3_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s5 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o3_b4 <- merge(S_s5_o3, biop_4, by = "pathway")
S_s5_o3_b4 <- S_s5_o3_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s5 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o3_b5 <- merge(S_s5_o3, biop_5, by = "pathway")
S_s5_o3_b5 <- S_s5_o3_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s5 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

##spheroid 5, org4
S_s5_o4_b1 <- merge(S_s5_o4, biop_1, by = "pathway")
S_s5_o4_b1 <- S_s5_o4_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s5 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o4_b4 <- merge(S_s5_o4, biop_4, by = "pathway")
S_s5_o4_b4 <- S_s5_o4_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s5 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o4_b5 <- merge(S_s5_o4, biop_5, by = "pathway")
S_s5_o4_b5 <- S_s5_o4_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s5 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#spheroid 5, org5
S_s5_o5_b1 <- merge(S_s5_o5, biop_1, by = "pathway")
S_s5_o5_b1 <- S_s5_o5_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s5 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o5_b4 <- merge(S_s5_o5, biop_4, by = "pathway")
S_s5_o5_b4 <- S_s5_o5_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s5 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o5_b5 <- merge(S_s5_o5, biop_5, by = "pathway")
S_s5_o5_b5 <- S_s5_o5_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s5 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#now add biopsy 035 layers:
#spheroid 6, org3
S_s6_o3_b1 <- merge(S_s6_o3, biop_1, by = "pathway")
S_s6_o3_b1 <- S_s6_o3_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s6 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o3_b4 <- merge(S_s6_o3, biop_4, by = "pathway")
S_s6_o3_b4 <- S_s6_o3_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s6 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o3_b5 <- merge(S_s6_o3, biop_5, by = "pathway")
S_s6_o3_b5 <- S_s6_o3_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s6 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

##spheroid 6, org4
S_s6_o4_b1 <- merge(S_s6_o4, biop_1, by = "pathway")
S_s6_o4_b1 <- S_s6_o4_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s6 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o4_b4 <- merge(S_s6_o4, biop_4, by = "pathway")
S_s6_o4_b4 <- S_s6_o4_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s6 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o4_b5 <- merge(S_s6_o4, biop_5, by = "pathway")
S_s6_o4_b5 <- S_s6_o4_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s6 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#spheroid 6, org5
S_s6_o5_b1 <- merge(S_s6_o5, biop_1, by = "pathway")
S_s6_o5_b1 <- S_s6_o5_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s6 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o5_b4 <- merge(S_s6_o5, biop_4, by = "pathway")
S_s6_o5_b4 <- S_s6_o5_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s6 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o5_b5 <- merge(S_s6_o5, biop_5, by = "pathway")
S_s6_o5_b5 <- S_s6_o5_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s6 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())




#  save ------
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/comparative_analysis_filt_biop/GSEA_enrch_comparison")

central <- grep("C_",names(.GlobalEnv),value=TRUE)
central_list <- do.call("list",mget(central))
middle <- grep("M_", names(.GlobalEnv), value = TRUE)
middle_list <- do.call("list", mget(middle))
surface <- grep("S_", names(.GlobalEnv), value = TRUE)
surface_list <- do.call("list", mget(surface))

#save lists to created directories
for(i in names(central_list)){
  write.csv(central_list[[i]], paste0("canonical_pathways/central/",i,".csv"))
}

for(i in names(middle_list)){
  write.csv(middle_list[[i]], paste0("canonical_pathways/middle/",i,".csv"))
}

for(i in names(surface_list)){
  write.csv(surface_list[[i]], paste0("canonical_pathways/surface/",i,".csv"))
}

#CHEMICAL & GENETIC PERTURBATIONS - filtered -----
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis")
getwd()
#
sph_3 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/3/3_c2.cgp.csv", header=TRUE, sep=",")
#remove unused columns and rename NES column
sph_3 <- sph_3 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s3 = NES, padj_s3 = padj)

org_1 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/1/1_c2.cgp.csv", header=TRUE, sep=",")
org_1 <- org_1 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o1 = NES, padj_o1 = padj)

org_6 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/6/6_c2.cgp.csv", header=TRUE, sep=",")
org_6 <- org_6 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o6 = NES, padj_o6 = padj)

biop_1 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/1/1_c2.cgp.csv", header=TRUE, sep=",")
biop_1 <- biop_1 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b1 = NES, padj_b1 = padj)

biop_2 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/2/2_c2.cgp.csv", header=TRUE, sep=",")
biop_2 <- biop_2 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b2 = NES, padj_b2 = padj)

sph_0 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/0/0_c2.cgp.csv", header=TRUE, sep=",")
sph_0 <- sph_0 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s0 = NES, padj_s0 = padj)
sph_1 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/1/1_c2.cgp.csv", header=TRUE, sep=",")
sph_1 <- sph_1 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s1 = NES, padj_s1 = padj)
sph_2 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/2/2_c2.cgp.csv", header=TRUE, sep=",")
sph_2 <- sph_2 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s2 = NES, padj_s2 = padj)
sph_4 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/4/4_c2.cgp.csv", header=TRUE, sep=",")
sph_4 <- sph_4 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s4 = NES, padj_s4 = padj)

org_0 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/0/0_c2.cgp.csv", header=TRUE, sep=",")
org_0 <- org_0 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o0 = NES, padj_o0 = padj)
org_2 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/2/2_c2.cgp.csv", header=TRUE, sep=",")
org_2 <- org_2 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o2 = NES, padj_o2 = padj)

biop_4 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/4/4_c2.cgp.csv", header=TRUE, sep=",")
biop_4 <- biop_4 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b4 = NES, padj_b4 = padj)

sph_5 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/5/5_c2.cgp.csv", header=TRUE, sep=",")
sph_5 <- sph_5 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s5 = NES, padj_s5 = padj)
sph_6 <- read.csv(file="PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/6/6_c2.cgp.csv", header=TRUE, sep=",")
sph_6 <- sph_6 %>%
  select(pathway, padj, NES) %>%
  rename(NES_s6 = NES, padj_s6 = padj)

org_3 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/3/3_c2.cgp.csv", header=TRUE, sep=",")
org_3 <- org_3 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o3 = NES, padj_o3 = padj)
org_4 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/4/4_c2.cgp.csv", header=TRUE, sep=",")
org_4 <- org_4 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o4 = NES, padj_o4 = padj)
org_5 <- read.csv(file = "organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/5/5_c2.cgp.csv", header=TRUE, sep=",")
org_5 <- org_5 %>%
  select(pathway, padj, NES) %>%
  rename(NES_o5 = NES, padj_o5 = padj)

biop_0 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/0/0_c2.cgp.csv", header=TRUE, sep=",")
biop_0 <- biop_0 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b0 = NES, padj_b0 = padj)
biop_3 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/3/3_c2.cgp.csv", header=TRUE, sep=",")
biop_3 <- biop_3 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b3 = NES, padj_b3 = padj)
biop_5 <- read.csv(file = "biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/5/5_c2.cgp.csv", header=TRUE, sep=",")
biop_5 <- biop_5 %>%
  select(pathway, padj, NES) %>%
  rename(NES_b5 = NES, padj_b5 = padj)

#  interactions-------
#spheroid and organoid1 clusters
C_s3_o1 <- merge(sph_3, org_1, by = "pathway")
C_s3_o1 <- C_s3_o1 %>%
  mutate(NESsum = NES_o1 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())
#now merge with biopsy0 cluster
C_s3_o1_b0 <- merge(C_s3_o1, biop_0, by = "pathway")
C_s3_o1_b0 <- C_s3_o1_b0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o1 + NES_s3 + NES_b0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
#and biopsy3
C_s3_o1_b3 <- merge(C_s3_o1, biop_3, by = "pathway")
C_s3_o1_b3 <- C_s3_o1_b3 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o1 + NES_s3 + NES_b3) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#spheroid and organoid_6
C_s3_o6 <- merge(sph_3, org_6, by = "pathway")
C_s3_o6 <- C_s3_o6 %>%
  mutate(NESsum = NES_o6 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())
#now merge with biopsy0 cluster
C_s3_o6_b0 <- merge(C_s3_o6, biop_0, by = "pathway")
C_s3_o6_b0 <- C_s3_o6_b0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o6 + NES_s3 + NES_b0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
#and biopsy3
C_s3_o6_b3 <- merge(C_s3_o6, biop_3, by = "pathway")
C_s3_o6_b3 <- C_s3_o6_b3 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o6 + NES_s3 + NES_b3) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#finally (for thouroughness)
C_s3_b0 <- merge(sph_3, biop_0, by = "pathway")
C_s3_b0 <- C_s3_b0 %>%
  mutate(NESsum = NES_b0 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_s3_b3 <- merge(sph_3, biop_3, by = "pathway")
C_s3_b3 <- C_s3_b3 %>%
  mutate(NESsum = NES_b3 + NES_s3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_o1_b0 <- merge(org_1, biop_0, by = "pathway")
C_o1_b0 <- C_o1_b0 %>%
  mutate(NESsum = NES_b0 + NES_o1) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_o1_b3 <- merge(org_1, biop_3, by = "pathway")
C_o1_b3 <- C_o1_b3 %>%
  mutate(NESsum = NES_b3 + NES_o1) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_o2_b0 <- merge(org_2, biop_0, by = "pathway")
C_o2_b0 <- C_o2_b0 %>%
  mutate(NESsum = NES_b0 + NES_o2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

C_o2_b3 <- merge(org_2, biop_3, by = "pathway")
C_o2_b3 <- C_o2_b3 %>%
  mutate(NESsum = NES_b3 + NES_o2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())


#MIDDLE
#biopsy2 and organoid0 clusters
M_b2_o0 <- merge(biop_2, org_0, by = "pathway")
M_b2_o0 <- M_b2_o0 %>%
  mutate(NESsum = NES_o0 + NES_b2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())
#biopsy2 and organoid2 clusters
M_b2_o2 <- merge(biop_2, org_2, by = "pathway")
M_b2_o2 <- M_b2_o2 %>%
  mutate(NESsum = NES_o2 + NES_b2) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

#merge organoid 0 with spheroid 0,1,2,4 clusters
M_b2_o0_s0 <- merge(M_b2_o0, sph_0, by = "pathway")
M_b2_o0_s0 <- M_b2_o0_s0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o0_s1 <- merge(M_b2_o0, sph_1, by = "pathway")
M_b2_o0_s1 <- M_b2_o0_s1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o0_s2 <- merge(M_b2_o0, sph_2, by = "pathway")
M_b2_o0_s2 <- M_b2_o0_s2 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s2) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o0_s4 <- merge(M_b2_o0, sph_4, by = "pathway")
M_b2_o0_s4 <- M_b2_o0_s4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o0 + NES_s4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#merge organoid 2 with spheroid 0,1,2,4 clusters
M_b2_o2_s0 <- merge(M_b2_o2, sph_0, by = "pathway")
M_b2_o2_s0 <- M_b2_o2_s0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o2_s1 <- merge(M_b2_o2, sph_1, by = "pathway")
M_b2_o2_s1 <- M_b2_o2_s1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o2_s2 <- merge(M_b2_o2, sph_2, by = "pathway")
M_b2_o2_s2 <- M_b2_o2_s2 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s2) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b2_o2_s4 <- merge(M_b2_o2, sph_4, by = "pathway")
M_b2_o2_s4 <- M_b2_o2_s4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b2 + NES_o2 + NES_s4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#biopsy3 and organoid0 clusters
M_b3_o0 <- merge(biop_3, org_0, by = "pathway")
M_b3_o0 <- M_b3_o0 %>%
  mutate(NESsum = NES_o0 + NES_b3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())
#biopsy2 and organoid2 clusters
M_b3_o2 <- merge(biop_3, org_2, by = "pathway")
M_b3_o2 <- M_b3_o2 %>%
  mutate(NESsum = NES_o2 + NES_b3) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

#merge organoid 0 with spheroid 0,1,2,4 clusters
M_b3_o0_s0 <- merge(M_b3_o0, sph_0, by = "pathway")
M_b3_o0_s0 <- M_b3_o0_s0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o0 + NES_s0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o0_s1 <- merge(M_b3_o0, sph_1, by = "pathway")
M_b3_o0_s1 <- M_b3_o0_s1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o0 + NES_s1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o0_s2 <- merge(M_b3_o0, sph_2, by = "pathway")
M_b3_o0_s2 <- M_b3_o0_s2 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o0 + NES_s2) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o0_s4 <- merge(M_b3_o0, sph_4, by = "pathway")
M_b3_o0_s4 <- M_b3_o0_s4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o0 + NES_s4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#merge organoid 2 with spheroid 0,1,2,4 clusters
M_b3_o2_s0 <- merge(M_b3_o2, sph_0, by = "pathway")
M_b3_o2_s0 <- M_b3_o2_s0 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o2 + NES_s0) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o2_s1 <- merge(M_b3_o2, sph_1, by = "pathway")
M_b3_o2_s1 <- M_b3_o2_s1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o2 + NES_s1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o2_s2 <- merge(M_b3_o2, sph_2, by = "pathway")
M_b3_o2_s2 <- M_b3_o2_s2 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o2 + NES_s2) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())
M_b3_o2_s4 <- merge(M_b3_o2, sph_4, by = "pathway")
M_b3_o2_s4 <- M_b3_o2_s4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_b3 + NES_o2 + NES_s4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())





#SURFACE
#spheroid5 and organoid clusters
S_s5_o3 <- merge(sph_5, org_3, by = "pathway")
S_s5_o3 <- S_s5_o3 %>%
  mutate(NESsum = NES_o3 + NES_s5) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s5_o4 <- merge(sph_5, org_4, by = "pathway")
S_s5_o4 <- S_s5_o4 %>%
  mutate(NESsum = NES_o4 + NES_s5) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s5_o5 <- merge(sph_5, org_5, by = "pathway")
S_s5_o5 <- S_s5_o5 %>%
  mutate(NESsum = NES_o5 + NES_s5) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

#spheroid6 and organoid clusters
S_s6_o3 <- merge(sph_6, org_3, by = "pathway")
S_s6_o3 <- S_s6_o3 %>%
  mutate(NESsum = NES_o3 + NES_s6) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s6_o4 <- merge(sph_6, org_4, by = "pathway")
S_s6_o4 <- S_s6_o4 %>%
  mutate(NESsum = NES_o4 + NES_s6) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())

S_s6_o5 <- merge(sph_6, org_5, by = "pathway")
S_s6_o5 <- S_s6_o5 %>%
  mutate(NESsum = NES_o5 + NES_s6) %>%
  arrange(-NESsum) %>%
  select(1, 2, 4, everything())



#now add biopsy (035) 1,4,5 layers:
#spheroid 5, org3
S_s5_o3_b1 <- merge(S_s5_o3, biop_1, by = "pathway")
S_s5_o3_b1 <- S_s5_o3_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s5 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o3_b4 <- merge(S_s5_o3, biop_4, by = "pathway")
S_s5_o3_b4 <- S_s5_o3_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s5 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o3_b5 <- merge(S_s5_o3, biop_5, by = "pathway")
S_s5_o3_b5 <- S_s5_o3_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s5 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

##spheroid 5, org4
S_s5_o4_b1 <- merge(S_s5_o4, biop_1, by = "pathway")
S_s5_o4_b1 <- S_s5_o4_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s5 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o4_b4 <- merge(S_s5_o4, biop_4, by = "pathway")
S_s5_o4_b4 <- S_s5_o4_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s5 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o4_b5 <- merge(S_s5_o4, biop_5, by = "pathway")
S_s5_o4_b5 <- S_s5_o4_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s5 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#spheroid 5, org5
S_s5_o5_b1 <- merge(S_s5_o5, biop_1, by = "pathway")
S_s5_o5_b1 <- S_s5_o5_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s5 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o5_b4 <- merge(S_s5_o5, biop_4, by = "pathway")
S_s5_o5_b4 <- S_s5_o5_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s5 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s5_o5_b5 <- merge(S_s5_o5, biop_5, by = "pathway")
S_s5_o5_b5 <- S_s5_o5_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s5 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#now add biopsy 035 layers:
#spheroid 6, org3
S_s6_o3_b1 <- merge(S_s6_o3, biop_1, by = "pathway")
S_s6_o3_b1 <- S_s6_o3_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s6 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o3_b4 <- merge(S_s6_o3, biop_4, by = "pathway")
S_s6_o3_b4 <- S_s6_o3_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s6 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o3_b5 <- merge(S_s6_o3, biop_5, by = "pathway")
S_s6_o3_b5 <- S_s6_o3_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o3 + NES_s6 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

##spheroid 6, org4
S_s6_o4_b1 <- merge(S_s6_o4, biop_1, by = "pathway")
S_s6_o4_b1 <- S_s6_o4_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s6 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o4_b4 <- merge(S_s6_o4, biop_4, by = "pathway")
S_s6_o4_b4 <- S_s6_o4_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s6 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o4_b5 <- merge(S_s6_o4, biop_5, by = "pathway")
S_s6_o4_b5 <- S_s6_o4_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o4 + NES_s6 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

#spheroid 6, org5
S_s6_o5_b1 <- merge(S_s6_o5, biop_1, by = "pathway")
S_s6_o5_b1 <- S_s6_o5_b1 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s6 + NES_b1) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o5_b4 <- merge(S_s6_o5, biop_4, by = "pathway")
S_s6_o5_b4 <- S_s6_o5_b4 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s6 + NES_b4) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())

S_s6_o5_b5 <- merge(S_s6_o5, biop_5, by = "pathway")
S_s6_o5_b5 <- S_s6_o5_b5 %>%
  select(-NESsum) %>%
  mutate(NESsum = NES_o5 + NES_s6 + NES_b5) %>%
  arrange(-NESsum) %>%
  select(1:3, 6, everything())




#  save ------
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/comparative_analysis_filt_biop/GSEA_enrch_comparison")

central <- grep("C_",names(.GlobalEnv),value=TRUE)
central_list <- do.call("list",mget(central))
middle <- grep("M_", names(.GlobalEnv), value = TRUE)
middle_list <- do.call("list", mget(middle))
surface <- grep("S_", names(.GlobalEnv), value = TRUE)
surface_list <- do.call("list", mget(surface))

#save lists to created directories
for(i in names(central_list)){
  write.csv(central_list[[i]], paste0("CG_perturbations/central/",i,".csv"))
}

for(i in names(middle_list)){
  write.csv(middle_list[[i]], paste0("CG_perturbations/middle/",i,".csv"))
}

for(i in names(surface_list)){
  write.csv(surface_list[[i]], paste0("CG_perturbations/surface/",i,".csv"))
}









