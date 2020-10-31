#generate a matrix and data frame that can be used to generate a heatmap showing correlation between clusters

library(RColorBrewer)
library(lattice)
library(latticeExtra)
library(grid)
library(gridExtra)
library(Seurat)
library(DT)
library(ComplexHeatmap)
library(tidyverse)
#define a function to change the name of an RDA file object loaded
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#load data ------

#spheroid
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
data("pGSEA")
data("PrerankSphere")
data("scaledClusterdDataGSEA")
dx=spheroid_CIOS_GSEA@assays$RNA@data
mm=t(apply(dx,1,tapply, spheroid_CIOS_GSEA@active.ident, calc.logMean))
mms = apply(dx,1,calc.logMean)
dd_spheroid = mm - mms
S_subsets <- loadRData("MSigDB/MSigDB_subsets.rda")
#give unique names
difs_s <- difs
means_s <- means
gsea_0s <- gsea_0
gsea_1s <- gsea_1
gsea_2s <- gsea_2
gsea_3s <- gsea_3
gsea_4s <- gsea_4
gsea_5s <- gsea_5
gsea_6s <- gsea_6


#organoid
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
data("pGSEA")
data("PrerankOrg")
data("scaledClusterdDataGSEA")
dx=Organoid_GSEA@assays$RNA@data
mm=t(apply(dx,1,tapply, Organoid_GSEA@active.ident, calc.logMean))
mms = apply(dx,1,calc.logMean)
dd_organoid = mm - mms
O_subsets <- loadRData("MSigDB/MSigDB_subsets.rda")
#give unique names
difs_o <- difs
means_o <- means
gsea_0o <- gsea_0
gsea_1o <- gsea_1
gsea_2o <- gsea_2
gsea_3o <- gsea_3
gsea_4o <- gsea_4
gsea_5o <- gsea_5
gsea_6o <- gsea_6



#biopsy
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320")
source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
data("pGSEA")
data("PrerankBiopsy")
data("scaledClusterdDataGSEA")
dx=biopsy_GSEA@assays$RNA@data
mm=t(apply(dx,1,tapply, biopsy_GSEA@active.ident, calc.logMean))
mms = apply(dx,1,calc.logMean)
dd_biop = mm - mms
B_subsets <- loadRData("MSigDB/MSigDB_subsets.rda")
#give unique names
difs_b <- difs
means_b <- means
gsea_0b <- gsea_0
gsea_1b <- gsea_1
gsea_2b <- gsea_2
gsea_3b <- gsea_3
gsea_4b <- gsea_4
gsea_5b <- gsea_5
#gsea_6b <- gsea_6

#Run loaded data through loops over all 'interesteing' gene sets from all gene collections -----
surface_clusters535_matrices <- list()
for (collection in c("h.all", "c5.bp", "c2.cgp", "c2.cp", "c4.cm", "c4.cgn")) {
  if (collection == "h.all"){
    hallmark_matrices <- list()
    for (pathway in c("HALLMARK_MYC_TARGETS_V1",
                      "HALLMARK_DNA_REPAIR",
                      "HALLMARK_MYC_TARGETS_V2",
                      "HALLMARK_E2F_TARGETS",
                      "HALLMARK_G2M_CHECKPOINT")) {
      #Spheroid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 5
      # point to correct 'gsea_x' object!!
      fid=which(gsea_5s$h.all$pathway==pathway)
      gseaRes=gsea_5s$h.all[fid, ]
      nameGset = gseaRes$pathway
      gset = S_subsets[[collection]][[nameGset]]
      ranks = difs_s$'5'
      names(ranks) = rownames(difs_s)
      LE_sphere = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_sphere = dd_spheroid[LE_sphere,]
      #Heatmap(mat_sphere)
      
      
      #Organoid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 3
      fid=which(gsea_3o$h.all$pathway==pathway)
      gseaRes=gsea_3o$h.all[fid, ]
      nameGset = gseaRes$pathway
      gset = O_subsets[[collection]][[nameGset]]
      ranks = difs_o$'3'
      names(ranks) = rownames(difs_o)
      LE_organoid = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_organoid=dd_organoid[LE_organoid,]
      #Heatmap(mat_organoid)
      
      #Biopsy ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 5
      fid=which(gsea_5b$h.all$pathway==pathway)
      gseaRes=gsea_5b$h.all[fid, ]
      nameGset = gseaRes$pathway
      gset = B_subsets[[collection]][[nameGset]]
      ranks = difs_b$'5'
      names(ranks) = rownames(difs_b)
      LE_biopsy = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_biopsy=dd_biop[LE_biopsy,]
      #Heatmap(mat_biopsy)
      
      
      #find LE overlap and combine in a matrix --------
      Common_LE <- intersect(LE_sphere, LE_organoid)
      Common_LE <- intersect(Common_LE, LE_biopsy)
      
      #next if Common_LE is too short
      if (length(Common_LE) < 3) {
        next
      }
      
      #create sub-matricies for plotting
      mat_sphere=t(dd_spheroid[Common_LE,])
      mat_organoid=t(dd_organoid[Common_LE,])
      mat_biopsy=t(dd_biop[Common_LE,])
      
      #rename layers by biology
      rownames(mat_sphere) <- paste(rownames(mat_sphere), "s", sep = "_")
      rownames(mat_organoid) <- paste(rownames(mat_organoid), "o", sep = "_")
      rownames(mat_biopsy) <- paste(rownames(mat_biopsy), "b", sep = "_")
      
      #bind new matrix
      BIND.MTX <- rbind(mat_sphere, mat_organoid, mat_biopsy)
      BIND.MTX <- t(BIND.MTX)
      hallmark_matrices[[pathway]] <- BIND.MTX
    }
    surface_clusters535_matrices[[collection]] <- hallmark_matrices
  }
  else if (collection == "c5.bp"){
    GObp_matrices <- list()
    for (pathway in c("GO_ANAPHASE_PROMOTING_COMPLEX_DEPENDENT_CATABOLIC_PROCESS",
                      "GO_CENTROMERE_COMPLEX_ASSEMBLY",
                      "GO_HISTONE_EXCHANGE",
                      "GO_DNA_REPLICATION_INDEPENDENT_NUCLEOSOME_ORGANIZATION",
                      "GO_POSITIVE_REGULATION_OF_LIGASE_ACTIVITY",
                      "GO_REGULATION_OF_PROTEIN_UBIQUITINATION_INVOLVED_IN_UBIQUITIN_DEPENDENT_PROTEIN_CATABOLIC_PROCESS",
                      "GO_REGULATION_OF_LIGASE_ACTIVITY",
                      "GO_PROTEIN_HETEROTETRAMERIZATION",
                      "GO_ATP_DEPENDENT_CHROMATIN_REMODELING",
                      "GO_DNA_REPLICATION_DEPENDENT_NUCLEOSOME_ORGANIZATION",
                      "GO_REGULATION_OF_MEGAKARYOCYTE_DIFFERENTIATION",
                      "GO_DNA_BIOSYNTHETIC_PROCESS",
                      "GO_DNA_STRAND_ELONGATION_INVOLVED_IN_DNA_REPLICATION",
                      "GO_MITOTIC_RECOMBINATION",
                      "GO_INTERSTRAND_CROSS_LINK_REPAIR",
                      "GO_NEGATIVE_REGULATION_OF_PROTEIN_MODIFICATION_BY_SMALL_PROTEIN_CONJUGATION_OR_REMOVAL",
                      "GO_TELOMERE_ORGANIZATION",
                      "GO_SPLICEOSOMAL_SNRNP_ASSEMBLY",
                      "GO_TELOMERE_MAINTENANCE_VIA_RECOMBINATION",
                      "GO_DNA_STRAND_ELONGATION")) {
      #Spheroid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 5
      # point to correct 'gsea_x' object!!
      fid=which(gsea_5s$c5.bp$pathway==pathway)
      gseaRes=gsea_5s$c5.bp[fid, ]
      nameGset = gseaRes$pathway
      gset = S_subsets[[collection]][[nameGset]]
      ranks = difs_s$'5'
      names(ranks) = rownames(difs_s)
      LE_sphere = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_sphere = dd_spheroid[LE_sphere,]
      #Heatmap(mat_sphere)
      
      
      #Organoid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 3
      fid=which(gsea_3o$c5.bp$pathway==pathway)
      gseaRes=gsea_3o$c5.bp[fid, ]
      nameGset = gseaRes$pathway
      gset = O_subsets[[collection]][[nameGset]]
      ranks = difs_o$'3'
      names(ranks) = rownames(difs_o)
      LE_organoid = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_organoid=dd_organoid[LE_organoid,]
      #Heatmap(mat_organoid)
      
      #Biopsy ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 5
      fid=which(gsea_5b$c5.bp$pathway==pathway)
      gseaRes=gsea_5b$c5.bp[fid, ]
      nameGset = gseaRes$pathway
      gset = B_subsets[[collection]][[nameGset]]
      ranks = difs_b$'5'
      names(ranks) = rownames(difs_b)
      LE_biopsy = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_biopsy=dd_biop[LE_biopsy,]
      #Heatmap(mat_biopsy)
      
      
      #find LE overlap and plot as combined matrices-------
      Common_LE <- intersect(LE_sphere, LE_organoid)
      Common_LE <- intersect(Common_LE, LE_biopsy)
      
      #next if Common_LE is too short
      if (length(Common_LE) < 3) {
        next
      }
      
      #create sub-matricies for plotting
      mat_sphere=t(dd_spheroid[Common_LE,])
      mat_organoid=t(dd_organoid[Common_LE,])
      mat_biopsy=t(dd_biop[Common_LE,])
      
      #rename layers by biology
      rownames(mat_sphere) <- paste(rownames(mat_sphere), "s", sep = "_")
      rownames(mat_organoid) <- paste(rownames(mat_organoid), "o", sep = "_")
      rownames(mat_biopsy) <- paste(rownames(mat_biopsy), "b", sep = "_")
      
      #bind new matrix
      BIND.MTX <- rbind(mat_sphere, mat_organoid, mat_biopsy)
      BIND.MTX <- t(BIND.MTX)
      GObp_matrices[[pathway]] <- BIND.MTX
    }
    surface_clusters535_matrices[[collection]] <- GObp_matrices
  }
  else if (collection == "c2.cgp"){
    CGP_matrices <- list()
    for (pathway in c("SOTIRIOU_BREAST_CANCER_GRADE_1_VS_3_UP",
                      "GRAHAM_NORMAL_QUIESCENT_VS_NORMAL_DIVIDING_DN",
                      "GRAHAM_CML_DIVIDING_VS_NORMAL_QUIESCENT_UP",
                      "ZHOU_CELL_CYCLE_GENES_IN_IR_RESPONSE_24HR",
                      "RHODES_UNDIFFERENTIATED_CANCER",
                      "CROONQUIST_IL6_DEPRIVATION_DN",
                      "MORI_LARGE_PRE_BII_LYMPHOCYTE_UP",
                      "CROONQUIST_NRAS_SIGNALING_DN",
                      "WINNEPENNINCKX_MELANOMA_METASTASIS_UP",
                      "BENPORATH_PROLIFERATION",
                      "WONG_EMBRYONIC_STEM_CELL_CORE",
                      "GAVIN_FOXP3_TARGETS_CLUSTER_P6",
                      "TARTE_PLASMA_CELL_VS_PLASMABLAST_DN",
                      "FERREIRA_EWINGS_SARCOMA_UNSTABLE_VS_STABLE_UP",
                      "MORI_PRE_BI_LYMPHOCYTE_UP",
                      "LY_AGING_OLD_DN",
                      "YU_MYC_TARGETS_UP",
                      "FERRANDO_T_ALL_WITH_MLL_ENL_FUSION_DN",
                      "CONCANNON_APOPTOSIS_BY_EPOXOMICIN_DN",
                      "BURTON_ADIPOGENESIS_PEAK_AT_24HR",
                      "CHANG_CORE_SERUM_RESPONSE_UP",
                      "MORI_MATURE_B_LYMPHOCYTE_DN",
                      "GARCIA_TARGETS_OF_FLI1_AND_DAX1_DN",
                      "CHICAS_RB1_TARGETS_LOW_SERUM",
                      "PAL_PRMT5_TARGETS_UP")) {
      #Spheroid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 5
      # point to correct 'gsea_x' object!!
      fid=which(gsea_5s$c2.cgp$pathway==pathway)
      gseaRes=gsea_5s$c2.cgp[fid, ]
      nameGset = gseaRes$pathway
      gset = S_subsets[[collection]][[nameGset]]
      ranks = difs_s$'5'
      names(ranks) = rownames(difs_s)
      LE_sphere = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_sphere = dd_spheroid[LE_sphere,]
      #Heatmap(mat_sphere)
      
      
      #Organoid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 3
      fid=which(gsea_3o$c2.cgp$pathway==pathway)
      gseaRes=gsea_3o$c2.cgp[fid, ]
      nameGset = gseaRes$pathway
      gset = O_subsets[[collection]][[nameGset]]
      ranks = difs_o$'3'
      names(ranks) = rownames(difs_o)
      LE_organoid = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_organoid=dd_organoid[LE_organoid,]
      #Heatmap(mat_organoid)
      
      #Biopsy ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 5
      fid=which(gsea_5b$c2.cgp$pathway==pathway)
      gseaRes=gsea_5b$c2.cgp[fid, ]
      nameGset = gseaRes$pathway
      gset = B_subsets[[collection]][[nameGset]]
      ranks = difs_b$'5'
      names(ranks) = rownames(difs_b)
      LE_biopsy = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_biopsy=dd_biop[LE_biopsy,]
      #Heatmap(mat_biopsy)
      
      
      #find LE overlap and plot as combined matrices-------
      Common_LE <- intersect(LE_sphere, LE_organoid)
      Common_LE <- intersect(Common_LE, LE_biopsy)
      
      #next if Common_LE is too short
      if (length(Common_LE) < 3) {
        next
      }
      
      #create sub-matricies for plotting
      mat_sphere=t(dd_spheroid[Common_LE,])
      mat_organoid=t(dd_organoid[Common_LE,])
      mat_biopsy=t(dd_biop[Common_LE,])
      
      #rename layers by biology
      rownames(mat_sphere) <- paste(rownames(mat_sphere), "s", sep = "_")
      rownames(mat_organoid) <- paste(rownames(mat_organoid), "o", sep = "_")
      rownames(mat_biopsy) <- paste(rownames(mat_biopsy), "b", sep = "_")
      
      #bind new matrix
      BIND.MTX <- rbind(mat_sphere, mat_organoid, mat_biopsy)
      BIND.MTX <- t(BIND.MTX)
      CGP_matrices[[pathway]] <- BIND.MTX
    }
    surface_clusters535_matrices[[collection]] <- CGP_matrices
  }
  else if (collection == "c2.cp"){
    CanonPW_matrices <- list()
    for (pathway in c("REACTOME_REGULATION_OF_MITOTIC_CELL_CYCLE",
                      "REACTOME_DNA_REPLICATION",
                      "REACTOME_CELL_CYCLE_CHECKPOINTS",
                      "REACTOME_CELL_CYCLE",
                      "REACTOME_MITOTIC_M_M_G1_PHASES",
                      "REACTOME_APC_C_CDC20_MEDIATED_DEGRADATION_OF_MITOTIC_PROTEINS",
                      "REACTOME_CELL_CYCLE_MITOTIC",
                      "REACTOME_APC_C_CDH1_MEDIATED_DEGRADATION_OF_CDC20_AND_OTHER_APC_C_CDH1_TARGETED_PROTEINS_IN_LATE_MITOSIS_EARLY_G1",
                      "REACTOME_SYNTHESIS_OF_DNA",
                      "REACTOME_G1_S_TRANSITION",
                      "REACTOME_M_G1_TRANSITION",
                      "REACTOME_S_PHASE",
                      "REACTOME_FORMATION_OF_TUBULIN_FOLDING_INTERMEDIATES_BY_CCT_TRIC",
                      "REACTOME_MITOTIC_G1_G1_S_PHASES",
                      "REACTOME_ASSEMBLY_OF_THE_PRE_REPLICATIVE_COMPLEX",
                      "REACTOME_CHROMOSOME_MAINTENANCE",
                      "REACTOME_DEPOSITION_OF_NEW_CENPA_CONTAINING_NUCLEOSOMES_AT_THE_CENTROMERE",
                      "REACTOME_ORC1_REMOVAL_FROM_CHROMATIN",
                      "REACTOME_SCFSKP2_MEDIATED_DEGRADATION_OF_P27_P21",
                      "REACTOME_AUTODEGRADATION_OF_CDH1_BY_CDH1_APC_C")) {
      #Spheroid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 5
      # point to correct 'gsea_x' object!!
      fid=which(gsea_5s$c2.cp$pathway==pathway)
      gseaRes=gsea_5s$c2.cp[fid, ]
      nameGset = gseaRes$pathway
      gset = S_subsets[[collection]][[nameGset]]
      ranks = difs_s$'5'
      names(ranks) = rownames(difs_s)
      LE_sphere = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_sphere = dd_spheroid[LE_sphere,]
      #Heatmap(mat_sphere)
      
      
      #Organoid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 3
      fid=which(gsea_3o$c2.cp$pathway==pathway)
      gseaRes=gsea_3o$c2.cp[fid, ]
      nameGset = gseaRes$pathway
      gset = O_subsets[[collection]][[nameGset]]
      ranks = difs_o$'3'
      names(ranks) = rownames(difs_o)
      LE_organoid = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_organoid=dd_organoid[LE_organoid,]
      #Heatmap(mat_organoid)
      
      #Biopsy ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 5
      fid=which(gsea_5b$c2.cp$pathway==pathway)
      gseaRes=gsea_5b$c2.cp[fid, ]
      nameGset = gseaRes$pathway
      gset = B_subsets[[collection]][[nameGset]]
      ranks = difs_b$'5'
      names(ranks) = rownames(difs_b)
      LE_biopsy = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_biopsy=dd_biop[LE_biopsy,]
      #Heatmap(mat_biopsy)
      
      
      #find LE overlap and plot as combined matrices-------
      Common_LE <- intersect(LE_sphere, LE_organoid)
      Common_LE <- intersect(Common_LE, LE_biopsy)
      
      #next if Common_LE is too short
      if (length(Common_LE) < 3) {
        next
      }
      
      #create sub-matricies for plotting
      mat_sphere=t(dd_spheroid[Common_LE,])
      mat_organoid=t(dd_organoid[Common_LE,])
      mat_biopsy=t(dd_biop[Common_LE,])
      
      #rename layers by biology
      rownames(mat_sphere) <- paste(rownames(mat_sphere), "s", sep = "_")
      rownames(mat_organoid) <- paste(rownames(mat_organoid), "o", sep = "_")
      rownames(mat_biopsy) <- paste(rownames(mat_biopsy), "b", sep = "_")
      
      #bind new matrix
      BIND.MTX <- rbind(mat_sphere, mat_organoid, mat_biopsy)
      BIND.MTX <- t(BIND.MTX)
      CanonPW_matrices[[pathway]] <- BIND.MTX
    }
    surface_clusters535_matrices[[collection]] <- CanonPW_matrices
  }
  else if (collection == "c4.cm"){
    CancerMod_matrices <- list()
    for (pathway in c("MODULE_54",
                      "MODULE_219",
                      "MODULE_91",
                      "MODULE_158",
                      "MODULE_28",
                      "MODULE_125",
                      "MODULE_388",
                      "MODULE_102",
                      "MODULE_61",
                      "MODULE_299")) {
      #Spheroid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 5
      # point to correct 'gsea_x' object!!
      fid=which(gsea_5s$c4.cm$pathway==pathway)
      gseaRes=gsea_5s$c4.cm[fid, ]
      nameGset = gseaRes$pathway
      gset = S_subsets[[collection]][[nameGset]]
      ranks = difs_s$'5'
      names(ranks) = rownames(difs_s)
      LE_sphere = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_sphere = dd_spheroid[LE_sphere,]
      #Heatmap(mat_sphere)
      
      
      #Organoid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 3
      fid=which(gsea_3o$c4.cm$pathway==pathway)
      gseaRes=gsea_3o$c4.cm[fid, ]
      nameGset = gseaRes$pathway
      gset = O_subsets[[collection]][[nameGset]]
      ranks = difs_o$'3'
      names(ranks) = rownames(difs_o)
      LE_organoid = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_organoid=dd_organoid[LE_organoid,]
      #Heatmap(mat_organoid)
      
      #Biopsy ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 5
      fid=which(gsea_5b$c4.cm$pathway==pathway)
      gseaRes=gsea_5b$c4.cm[fid, ]
      nameGset = gseaRes$pathway
      gset = B_subsets[[collection]][[nameGset]]
      ranks = difs_b$'5'
      names(ranks) = rownames(difs_b)
      LE_biopsy = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_biopsy=dd_biop[LE_biopsy,]
      #Heatmap(mat_biopsy)
      
      
      #find LE overlap and plot as combined matrices-------
      Common_LE <- intersect(LE_sphere, LE_organoid)
      Common_LE <- intersect(Common_LE, LE_biopsy)
      
      #next if Common_LE is too short
      if (length(Common_LE) < 3) {
        next
      }
      
      #create sub-matricies for plotting
      mat_sphere=t(dd_spheroid[Common_LE,])
      mat_organoid=t(dd_organoid[Common_LE,])
      mat_biopsy=t(dd_biop[Common_LE,])
      
      #rename layers by biology
      rownames(mat_sphere) <- paste(rownames(mat_sphere), "s", sep = "_")
      rownames(mat_organoid) <- paste(rownames(mat_organoid), "o", sep = "_")
      rownames(mat_biopsy) <- paste(rownames(mat_biopsy), "b", sep = "_")
      
      #bind new matrix
      BIND.MTX <- rbind(mat_sphere, mat_organoid, mat_biopsy)
      BIND.MTX <- t(BIND.MTX)
      CancerMod_matrices[[pathway]] <- BIND.MTX
    }
    surface_clusters535_matrices[[collection]] <- CancerMod_matrices
  }
  else if (collection == "c4.cgn"){
    CancerGeneNet_matrices <- list()
    for (pathway in c("GNF2_RRM1",
                      "GNF2_RAN",
                      "GNF2_RFC4",
                      "GNF2_PA2G4",
                      "MORF_PCNA",
                      "GNF2_CKS1B",
                      "GNF2_MCM4",
                      "GNF2_CKS2",
                      "GNF2_RFC3",
                      "MORF_BUB3",
                      "MORF_FEN1",
                      "GNF2_ESPL1",
                      "MORF_FBL",
                      "GNF2_BUB1",
                      "MORF_CSNK2B",
                      "MORF_UNG",
                      "MORF_RAD23A",
                      "MORF_HAT1",
                      "MORF_MAP2K2",
                      "MORF_ANP32B")) {
      #Spheroid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 5
      # point to correct 'gsea_x' object!!
      fid=which(gsea_5s$c4.cgn$pathway==pathway)
      gseaRes=gsea_5s$c4.cgn[fid, ]
      nameGset = gseaRes$pathway
      gset = S_subsets[[collection]][[nameGset]]
      ranks = difs_s$'5'
      names(ranks) = rownames(difs_s)
      LE_sphere = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_sphere = dd_spheroid[LE_sphere,]
      #Heatmap(mat_sphere)
      
      
      #Organoid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 3
      fid=which(gsea_3o$c4.cgn$pathway==pathway)
      gseaRes=gsea_3o$c4.cgn[fid, ]
      nameGset = gseaRes$pathway
      gset = O_subsets[[collection]][[nameGset]]
      ranks = difs_o$'3'
      names(ranks) = rownames(difs_o)
      LE_organoid = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_organoid=dd_organoid[LE_organoid,]
      #Heatmap(mat_organoid)
      
      #Biopsy ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 5
      fid=which(gsea_5b$c4.cgn$pathway==pathway)
      gseaRes=gsea_5b$c4.cgn[fid, ]
      nameGset = gseaRes$pathway
      gset = B_subsets[[collection]][[nameGset]]
      ranks = difs_b$'5'
      names(ranks) = rownames(difs_b)
      LE_biopsy = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_biopsy=dd_biop[LE_biopsy,]
      #Heatmap(mat_biopsy)
      
      
      #find LE overlap and plot as combined matrices-------
      Common_LE <- intersect(LE_sphere, LE_organoid)
      Common_LE <- intersect(Common_LE, LE_biopsy)
      
      #next if Common_LE is too short
      if (length(Common_LE) < 3) {
        next
      }
      
      #create sub-matricies for plotting
      mat_sphere=t(dd_spheroid[Common_LE,])
      mat_organoid=t(dd_organoid[Common_LE,])
      mat_biopsy=t(dd_biop[Common_LE,])
      
      #rename layers by biology
      rownames(mat_sphere) <- paste(rownames(mat_sphere), "s", sep = "_")
      rownames(mat_organoid) <- paste(rownames(mat_organoid), "o", sep = "_")
      rownames(mat_biopsy) <- paste(rownames(mat_biopsy), "b", sep = "_")
      
      #bind new matrix
      BIND.MTX <- rbind(mat_sphere, mat_organoid, mat_biopsy)
      BIND.MTX <- t(BIND.MTX)
      CancerGeneNet_matrices[[pathway]] <- BIND.MTX
    }
    surface_clusters535_matrices[[collection]] <- CancerGeneNet_matrices
  } 
  else{
    print("no loop found for this pathway")
  }
}

#save list of lists as RDS -------
saveRDS(surface_clusters535_matrices, file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/comparative_analysis_filt_biop/clusters_heatmap/comparison_data/surface535.rds")

#unlist and bind into a large matrix called surface_535
grouped_families <- list()
for (i in seq_along(surface_clusters535_matrices)) {
  grouped_families[[i]] <- do.call(rbind, surface_clusters535_matrices[[i]])
}
surface_535 <- do.call(rbind, grouped_families)
surface_535 <- unique(surface_535)
Heatmap(surface_535)
Heatmap(unique(grouped_families[[1]]))
hallmark_genes <- grouped_families[[1]]
hallmark_genes2 <- unique(grouped_families[[1]])

Heatmap(hallmark_genes2)

all_genes_unique_surface535 <- surface_535
saveRDS(all_genes_unique_surface535, file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/comparative_analysis_filt_biop/clusters_heatmap/comparison_data/all_genes_unique_surface535.rds")


#make a data frame that annotates the overlapping LE from important gene sets
surface_535_annotation <- data.frame(row.names = rownames(surface_535))
n=0
for (i in seq_along(surface_clusters535_matrices)) {
  for (j in seq_along(surface_clusters535_matrices[[i]])) {
    n=n+1
    surface_535_annotation[[n]] <- rownames(surface_535_annotation) %in% rownames(surface_clusters535_matrices[[i]][[j]])
    colnames(surface_535_annotation)[n] <- names(surface_clusters535_matrices[[i]][j])
  }
}
n=0

saveRDS(surface_535_annotation, file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/comparative_analysis_filt_biop/clusters_heatmap/comparison_data/surface_535_annotation.rds")


#try to plot heatmap with layers clustered
cluster_positions <- data.frame(row.names = colnames(surface_535))
cluster_positions$biology <- c("m", "m", "m", "c", "m", "s", "s", "c", "s", "m", "s", "c", "m", "s", "c", "s", "m", "c", "s", "s")
col.biology <- list("biology" = c("c" = "blue", "m" = "green3", "s" = "red"))

Heatmap(hallmark_genes2, name = '535', 
        column_split = factor(cluster_positions$biology, levels = c("s", "m", "c")),
        row_km = 4,
        cluster_column_slices = FALSE)

#reorder columns for kmeans clustering
kclus <- kmeans(t(surface_535), 2)
kclus$cluster

#custom ordering
split <- factor(paste0("Cluster\n", kclus$cluster), levels=c("Cluster\n2","Cluster\n1"))
#reorder.hmap <- Heatmap(surface_535, column_split=split)
Heatmap(surface_535, name = '535',
        column_split=split, 
        row_km = 4,
        bottom_annotation = HeatmapAnnotation(df = cluster_positions, col = col.biology))

Heatmap(surface_535, name = '535',
        #column_split=split, 
        column_km = 3,
        row_km = 1,
        #row_split = surface_535_annotation[1:2],
        #right_annotation = HeatmapAnnotation(df = surface_535_annotation[5], which = "row"),
        bottom_annotation = HeatmapAnnotation(df = cluster_positions, col = col.biology))



#end-----