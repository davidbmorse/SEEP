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

#load data

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
  #give unique names----
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
  #give unique names----
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
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
data("pGSEA")
data("PrerankBiopsy")
data("scaledClusterdDataGSEA")
dx=biopsy_GSEA@assays$RNA@data
mm=t(apply(dx,1,tapply, biopsy_GSEA@active.ident, calc.logMean))
mms = apply(dx,1,calc.logMean)
dd_biop = mm - mms
B_subsets <- loadRData("MSigDB/MSigDB_subsets.rda")
  #give unique names----
  difs_b <- difs
  means_b <- means
  gsea_0b <- gsea_0
  gsea_1b <- gsea_1
  gsea_2b <- gsea_2
  gsea_3b <- gsea_3
  gsea_4b <- gsea_4
  gsea_5b <- gsea_5
#  gsea_6b <- gsea_6

#Run loaded data through loops over all 'interesteing' gene sets from all gene collections

  surface_clusters644_matrices <- list()
  for (collection in c("h.all", "c5.bp", "c2.cgp", "c2.cp", "c4.cm", "c4.cgn")) {
    if (collection == "h.all"){
      hallmark_matrices <- list()
      for (pathway in c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                        "HALLMARK_IL6_JAK_STAT3_SIGNALING",
                        "HALLMARK_INFLAMMATORY_RESPONSE",
                        "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                        "HALLMARK_APOPTOSIS",
                        "HALLMARK_COMPLEMENT",
                        "HALLMARK_COAGULATION",
                        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                        "HALLMARK_KRAS_SIGNALING_UP",
                        "HALLMARK_ALLOGRAFT_REJECTION",
                        "HALLMARK_P53_PATHWAY")) {
        #Spheroid ---------
        setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
        source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
        
        #cluster 6
        # point to correct 'gsea_x' object!!
        fid=which(gsea_6s$h.all$pathway==pathway)
        gseaRes=gsea_6s$h.all[fid, ]
        nameGset = gseaRes$pathway
        gset = S_subsets[[collection]][[nameGset]]
        ranks = difs_s$'6'
        names(ranks) = rownames(difs_s)
        LE_sphere = sort(unlist(gseaRes$leadingEdge))
        
        #create sub-matrix for plotting
        mat_sphere = dd_spheroid[LE_sphere,]
        #Heatmap(mat_sphere)
        
        
        #Organoid ---------
        setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
        source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
        
        #cluster 4
        fid=which(gsea_4o$h.all$pathway==pathway)
        gseaRes=gsea_4o$h.all[fid, ]
        nameGset = gseaRes$pathway
        gset = O_subsets[[collection]][[nameGset]]
        ranks = difs_o$'4'
        names(ranks) = rownames(difs_o)
        LE_organoid = sort(unlist(gseaRes$leadingEdge))
        
        #create sub-matrix for plotting
        mat_organoid=dd_organoid[LE_organoid,]
        #Heatmap(mat_organoid)
     
        #Biopsy ---------
        setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
        source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
        
        #cluster 4
        fid=which(gsea_4b$h.all$pathway==pathway)
        gseaRes=gsea_4b$h.all[fid, ]
        nameGset = gseaRes$pathway
        gset = B_subsets[[collection]][[nameGset]]
        ranks = difs_b$'4'
        names(ranks) = rownames(difs_b)
        LE_biopsy = sort(unlist(gseaRes$leadingEdge))
        
        #create sub-matrix for plotting
        mat_biopsy=dd_biop[LE_biopsy,]
        #Heatmap(mat_biopsy)
        
        
        #find LE overlap and combine in a matrix --------
        Common_LE <- intersect(LE_sphere, LE_organoid)
        Common_LE <- intersect(Common_LE, LE_biopsy)
        
        #create sub-matricies for plotting
        mat_sphere=t(dd_spheroid[Common_LE,])
        mat_organoid=t(dd_organoid[Common_LE,])
        mat_biopsy=t(dd_biop[Common_LE,])
        
        #rename layers by biology and select only surface and center layers
        rownames(mat_sphere) <- paste(rownames(mat_sphere), "s", sep = "_")
        rownames(mat_organoid) <- paste(rownames(mat_organoid), "o", sep = "_")
        rownames(mat_biopsy) <- paste(rownames(mat_biopsy), "b", sep = "_")
        
        #bind new matrix
        BIND.MTX <- rbind(mat_sphere, mat_organoid, mat_biopsy)
        BIND.MTX <- t(BIND.MTX)
        hallmark_matrices[[pathway]] <- BIND.MTX
      }
      surface_clusters644_matrices[[collection]] <- hallmark_matrices
    }
    else if (collection == "c5.bp"){
      GObp_matrices <- list()
      for (pathway in c("GO_INTERFERON_GAMMA_MEDIATED_SIGNALING_PATHWAY",
                        "GO_RESPONSE_TO_TYPE_I_INTERFERON",
                        "GO_RESPONSE_TO_INTERFERON_GAMMA",
                        "GO_CELLULAR_RESPONSE_TO_INTERFERON_GAMMA",
                        "GO_CYTOKINE_MEDIATED_SIGNALING_PATHWAY",
                        "GO_REGULATION_OF_LEUKOCYTE_MIGRATION",
                        "GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_ENDOGENOUS_ANTIGEN",
                        "GO_MONOCYTE_CHEMOTAXIS",
                        "GO_POSITIVE_REGULATION_OF_INFLAMMATORY_RESPONSE",
                        "GO_INFLAMMATORY_RESPONSE",
                        "GO_POSITIVE_REGULATION_OF_CHEMOTAXIS",
                        "GO_POSITIVE_REGULATION_OF_ACUTE_INFLAMMATORY_RESPONSE",
                        "GO_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION",
                        "GO_NEGATIVE_REGULATION_OF_IMMUNE_RESPONSE",
                        "GO_REGULATION_OF_CHEMOTAXIS",
                        "GO_POSITIVE_REGULATION_OF_CHEMOKINE_PRODUCTION",
                        "GO_POSITIVE_REGULATION_OF_CYTOKINE_BIOSYNTHETIC_PROCESS",
                        "GO_POSITIVE_REGULATION_OF_RESPONSE_TO_EXTERNAL_STIMULUS",
                        "GO_POSITIVE_REGULATION_OF_PHAGOCYTOSIS",
                        "GO_REGULATION_OF_TYPE_2_IMMUNE_RESPONSE",
                        "GO_REGULATION_OF_INFLAMMATORY_RESPONSE",
                        "GO_MYELOID_LEUKOCYTE_MIGRATION",
                        "GO_REGULATION_OF_ALPHA_BETA_T_CELL_PROLIFERATION",
                        "GO_DEFENSE_RESPONSE_TO_BACTERIUM",
                        "GO_ACUTE_INFLAMMATORY_RESPONSE",
                        "GO_POSITIVE_REGULATION_OF_ALPHA_BETA_T_CELL_PROLIFERATION",
                        "GO_ANTIMICROBIAL_HUMORAL_RESPONSE")) {
        #Spheroid ---------
        setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
        source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
        
        #cluster 6
        # point to correct 'gsea_x' object!!
        fid=which(gsea_6s$c5.bp$pathway==pathway)
        gseaRes=gsea_6s$c5.bp[fid, ]
        nameGset = gseaRes$pathway
        gset = S_subsets[[collection]][[nameGset]]
        ranks = difs_s$'6'
        names(ranks) = rownames(difs_s)
        LE_sphere = sort(unlist(gseaRes$leadingEdge))
        
        #create sub-matrix for plotting
        mat_sphere = dd_spheroid[LE_sphere,]
        #Heatmap(mat_sphere)
        
        
        #Organoid ---------
        setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
        source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
        
        #cluster 4
        fid=which(gsea_4o$c5.bp$pathway==pathway)
        gseaRes=gsea_4o$c5.bp[fid, ]
        nameGset = gseaRes$pathway
        gset = O_subsets[[collection]][[nameGset]]
        ranks = difs_o$'4'
        names(ranks) = rownames(difs_o)
        LE_organoid = sort(unlist(gseaRes$leadingEdge))
        
        #create sub-matrix for plotting
        mat_organoid=dd_organoid[LE_organoid,]
        #Heatmap(mat_organoid)
        
        #Biopsy ---------
        setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
        source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
        
        #cluster 4
        fid=which(gsea_4b$c5.bp$pathway==pathway)
        gseaRes=gsea_4b$c5.bp[fid, ]
        nameGset = gseaRes$pathway
        gset = B_subsets[[collection]][[nameGset]]
        ranks = difs_b$'4'
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
        
        #rename layers by biology and select only surface and center layers
        rownames(mat_sphere) <- paste(rownames(mat_sphere), "s", sep = "_")
        rownames(mat_organoid) <- paste(rownames(mat_organoid), "o", sep = "_")
        rownames(mat_biopsy) <- paste(rownames(mat_biopsy), "b", sep = "_")
        
        #bind new matrix
        BIND.MTX <- rbind(mat_sphere, mat_organoid, mat_biopsy)
        BIND.MTX <- t(BIND.MTX)
        GObp_matrices[[pathway]] <- BIND.MTX
      }
      surface_clusters644_matrices[[collection]] <- GObp_matrices
    }
    else if (collection == "c2.cgp"){
      CGP_matrices <- list()
      for (pathway in c("ZHANG_RESPONSE_TO_IKK_INHIBITOR_AND_TNF_UP",
                        "SANA_TNF_SIGNALING_UP",
                        "SEKI_INFLAMMATORY_RESPONSE_LPS_UP",
                        "HINATA_NFKB_TARGETS_KERATINOCYTE_UP",
                        "ALTEMEIER_RESPONSE_TO_LPS_WITH_MECHANICAL_VENTILATION",
                        "MISSIAGLIA_REGULATED_BY_METHYLATION_UP",
                        "PHONG_TNF_TARGETS_UP")) {
        #Spheroid ---------
        setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
        source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
        
        #cluster 6
        # point to correct 'gsea_x' object!!
        fid=which(gsea_6s$c2.cgp$pathway==pathway)
        gseaRes=gsea_6s$c2.cgp[fid, ]
        nameGset = gseaRes$pathway
        gset = S_subsets[[collection]][[nameGset]]
        ranks = difs_s$'6'
        names(ranks) = rownames(difs_s)
        LE_sphere = sort(unlist(gseaRes$leadingEdge))
        
        #create sub-matrix for plotting
        mat_sphere = dd_spheroid[LE_sphere,]
        #Heatmap(mat_sphere)
        
        
        #Organoid ---------
        setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
        source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
        
        #cluster 4
        fid=which(gsea_4o$c2.cgp$pathway==pathway)
        gseaRes=gsea_4o$c2.cgp[fid, ]
        nameGset = gseaRes$pathway
        gset = O_subsets[[collection]][[nameGset]]
        ranks = difs_o$'4'
        names(ranks) = rownames(difs_o)
        LE_organoid = sort(unlist(gseaRes$leadingEdge))
        
        #create sub-matrix for plotting
        mat_organoid=dd_organoid[LE_organoid,]
        #Heatmap(mat_organoid)
        
        #Biopsy ---------
        setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
        source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
        
        #cluster 4
        fid=which(gsea_4b$c2.cgp$pathway==pathway)
        gseaRes=gsea_4b$c2.cgp[fid, ]
        nameGset = gseaRes$pathway
        gset = B_subsets[[collection]][[nameGset]]
        ranks = difs_b$'4'
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
        
        #rename layers by biology and select only surface and center layers
        rownames(mat_sphere) <- paste(rownames(mat_sphere), "s", sep = "_")
        rownames(mat_organoid) <- paste(rownames(mat_organoid), "o", sep = "_")
        rownames(mat_biopsy) <- paste(rownames(mat_biopsy), "b", sep = "_")
        
        #bind new matrix
        BIND.MTX <- rbind(mat_sphere, mat_organoid, mat_biopsy)
        BIND.MTX <- t(BIND.MTX)
        CGP_matrices[[pathway]] <- BIND.MTX
      }
      surface_clusters644_matrices[[collection]] <- CGP_matrices
    }
    else if (collection == "c2.cp"){
      CanonPW_matrices <- list()
      for (pathway in c("KEGG_GRAFT_VERSUS_HOST_DISEASE",
                        "KEGG_ALLOGRAFT_REJECTION",
                        "KEGG_AUTOIMMUNE_THYROID_DISEASE",
                        "PID_IL23_PATHWAY",
                        "KEGG_HEMATOPOIETIC_CELL_LINEAGE",
                        "REACTOME_INTERFERON_GAMMA_SIGNALING",
                        "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
                        "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
                        "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
                        "KEGG_TYPE_I_DIABETES_MELLITUS",
                        "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
                        "PID_INTEGRIN2_PATHWAY",
                        "REACTOME_PEPTIDE_LIGAND_BINDING_RECEPTORS",
                        "PID_CD40_PATHWAY",
                        "REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES",
                        "ST_TUMOR_NECROSIS_FACTOR_PATHWAY",
                        "BIOCARTA_TH1TH2_PATHWAY",
                        "BIOCARTA_IL1R_PATHWAY",
                        "PID_IL27_PATHWAY",
                        "BIOCARTA_CD40_PATHWAY",
                        "BIOCARTA_TNFR2_PATHWAY",
                        "BIOCARTA_INFLAM_PATHWAY",
                        "BIOCARTA_NFKB_PATHWAY",
                        "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES",
                        "PID_TNF_PATHWAY",
                        "BIOCARTA_NTHI_PATHWAY",
                        "KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY",
                        "KEGG_ASTHMA",
                        "NABA_ECM_REGULATORS",
                        "REACTOME_CLASS_A1_RHODOPSIN_LIKE_RECEPTORS")) {
        #Spheroid ---------
        setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
        source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
        
        #cluster 6
        # point to correct 'gsea_x' object!!
        fid=which(gsea_6s$c2.cp$pathway==pathway)
        gseaRes=gsea_6s$c2.cp[fid, ]
        nameGset = gseaRes$pathway
        gset = S_subsets[[collection]][[nameGset]]
        ranks = difs_s$'6'
        names(ranks) = rownames(difs_s)
        LE_sphere = sort(unlist(gseaRes$leadingEdge))
        
        #create sub-matrix for plotting
        mat_sphere = dd_spheroid[LE_sphere,]
        #Heatmap(mat_sphere)
        
        
        #Organoid ---------
        setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
        source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
        
        #cluster 4
        fid=which(gsea_4o$c2.cp$pathway==pathway)
        gseaRes=gsea_4o$c2.cp[fid, ]
        nameGset = gseaRes$pathway
        gset = O_subsets[[collection]][[nameGset]]
        ranks = difs_o$'4'
        names(ranks) = rownames(difs_o)
        LE_organoid = sort(unlist(gseaRes$leadingEdge))
        
        #create sub-matrix for plotting
        mat_organoid=dd_organoid[LE_organoid,]
        #Heatmap(mat_organoid)
        
        #Biopsy ---------
        setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
        source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
        
        #cluster 4
        fid=which(gsea_4b$c2.cp$pathway==pathway)
        gseaRes=gsea_4b$c2.cp[fid, ]
        nameGset = gseaRes$pathway
        gset = B_subsets[[collection]][[nameGset]]
        ranks = difs_b$'4'
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
        
        #rename layers by biology and select only surface and center layers
        rownames(mat_sphere) <- paste(rownames(mat_sphere), "s", sep = "_")
        rownames(mat_organoid) <- paste(rownames(mat_organoid), "o", sep = "_")
        rownames(mat_biopsy) <- paste(rownames(mat_biopsy), "b", sep = "_")
        
        #bind new matrix
        BIND.MTX <- rbind(mat_sphere, mat_organoid, mat_biopsy)
        BIND.MTX <- t(BIND.MTX)
        CanonPW_matrices[[pathway]] <- BIND.MTX
      }
      surface_clusters644_matrices[[collection]] <- CanonPW_matrices
    }
    else if (collection == "c4.cm"){
      CancerMod_matrices <- list()
      for (pathway in c("MODULE_46",
                        "MODULE_75",
                        "MODULE_76",
                        "MODULE_345",
                        "MODULE_5",
                        "MODULE_436",
                        "MODULE_44",
                        "MODULE_79",
                        "MODULE_170",
                        "MODULE_128",
                        "MODULE_223",
                        "MODULE_263",
                        "MODULE_300",
                        "MODULE_6",
                        "MODULE_92",
                        "MODULE_73",
                        "MODULE_27",
                        "MODULE_108",
                        "MODULE_121",
                        "MODULE_521",
                        "MODULE_488",
                        "MODULE_265",
                        "MODULE_241",
                        "MODULE_291",
                        "MODULE_122")) {
        #Spheroid ---------
        setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
        source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
        
        #cluster 6
        # point to correct 'gsea_x' object!!
        fid=which(gsea_6s$c4.cm$pathway==pathway)
        gseaRes=gsea_6s$c4.cm[fid, ]
        nameGset = gseaRes$pathway
        gset = S_subsets[[collection]][[nameGset]]
        ranks = difs_s$'6'
        names(ranks) = rownames(difs_s)
        LE_sphere = sort(unlist(gseaRes$leadingEdge))
        
        #create sub-matrix for plotting
        mat_sphere = dd_spheroid[LE_sphere,]
        #Heatmap(mat_sphere)
        
        
        #Organoid ---------
        setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
        source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
        
        #cluster 4
        fid=which(gsea_4o$c4.cm$pathway==pathway)
        gseaRes=gsea_4o$c4.cm[fid, ]
        nameGset = gseaRes$pathway
        gset = O_subsets[[collection]][[nameGset]]
        ranks = difs_o$'4'
        names(ranks) = rownames(difs_o)
        LE_organoid = sort(unlist(gseaRes$leadingEdge))
        
        #create sub-matrix for plotting
        mat_organoid=dd_organoid[LE_organoid,]
        #Heatmap(mat_organoid)
        
        #Biopsy ---------
        setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
        source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
        
        #cluster 4
        fid=which(gsea_4b$c4.cm$pathway==pathway)
        gseaRes=gsea_4b$c4.cm[fid, ]
        nameGset = gseaRes$pathway
        gset = B_subsets[[collection]][[nameGset]]
        ranks = difs_b$'4'
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
        
        #rename layers by biology and select only surface and center layers
        rownames(mat_sphere) <- paste(rownames(mat_sphere), "s", sep = "_")
        rownames(mat_organoid) <- paste(rownames(mat_organoid), "o", sep = "_")
        rownames(mat_biopsy) <- paste(rownames(mat_biopsy), "b", sep = "_")
        
        #bind new matrix
        BIND.MTX <- rbind(mat_sphere, mat_organoid, mat_biopsy)
        BIND.MTX <- t(BIND.MTX)
        CancerMod_matrices[[pathway]] <- BIND.MTX
      }
      surface_clusters644_matrices[[collection]] <- CancerMod_matrices
    }
    else if (collection == "c4.cgn"){
      CancerGeneNet_matrices <- list()
      for (pathway in c("GNF2_HLA-C",
                        "GNF2_CD14",
                        "GNF2_CD48",
                        "GNF2_PECAM1",
                        "GNF2_CD53",
                        "GNF2_PTPN6",
                        "GNF2_INPP5D",
                        "GNF2_CARD15",
                        "GNF2_GLTSCR2",
                        "GNF2_MMP11",
                        "GNF2_FOS",
                        "GNF2_VAV1",
                        "GNF2_HCK",
                        "GNF2_CDKN1C",
                        "GNF2_CD1D",
                        "GNF2_EGFR",
                        "GNF2_CASP1",
                        "GNF2_TPT1",
                        "GNF2_KISS1",
                        "GNF2_CD33")) {
        #Spheroid ---------
        setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
        source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
        
        #cluster 6
        # point to correct 'gsea_x' object!!
        fid=which(gsea_6s$c4.cgn$pathway==pathway)
        gseaRes=gsea_6s$c4.cgn[fid, ]
        nameGset = gseaRes$pathway
        gset = S_subsets[[collection]][[nameGset]]
        ranks = difs_s$'6'
        names(ranks) = rownames(difs_s)
        LE_sphere = sort(unlist(gseaRes$leadingEdge))
        
        #create sub-matrix for plotting
        mat_sphere = dd_spheroid[LE_sphere,]
        #Heatmap(mat_sphere)
        
        
        #Organoid ---------
        setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
        source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
        
        #cluster 4
        fid=which(gsea_4o$c4.cgn$pathway==pathway)
        gseaRes=gsea_4o$c4.cgn[fid, ]
        nameGset = gseaRes$pathway
        gset = O_subsets[[collection]][[nameGset]]
        ranks = difs_o$'4'
        names(ranks) = rownames(difs_o)
        LE_organoid = sort(unlist(gseaRes$leadingEdge))
        
        #create sub-matrix for plotting
        mat_organoid=dd_organoid[LE_organoid,]
        #Heatmap(mat_organoid)
        
        #Biopsy ---------
        setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
        source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
        
        #cluster 4
        fid=which(gsea_4b$c4.cgn$pathway==pathway)
        gseaRes=gsea_4b$c4.cgn[fid, ]
        nameGset = gseaRes$pathway
        gset = B_subsets[[collection]][[nameGset]]
        ranks = difs_b$'4'
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
        
        #rename layers by biology and select only surface and center layers
        rownames(mat_sphere) <- paste(rownames(mat_sphere), "s", sep = "_")
        rownames(mat_organoid) <- paste(rownames(mat_organoid), "o", sep = "_")
        rownames(mat_biopsy) <- paste(rownames(mat_biopsy), "b", sep = "_")
        
        #bind new matrix
        BIND.MTX <- rbind(mat_sphere, mat_organoid, mat_biopsy)
        BIND.MTX <- t(BIND.MTX)
        CancerGeneNet_matrices[[pathway]] <- BIND.MTX
      }
      surface_clusters644_matrices[[collection]] <- CancerGeneNet_matrices
    } 
    else{
      print("no loop found for this pathway")
    }
  }
  
#save list of lists as RDS
saveRDS(surface_clusters644_matrices, file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/comparative_analysis_filt_biop/clusters_heatmap/comparison_data/surface644.rds")
  
#unlist and bind into a large matrix called surface_644
grouped_families <- list()
for (i in seq_along(surface_clusters644_matrices)) {
  grouped_families[[i]] <- do.call(rbind, surface_clusters644_matrices[[i]])
}
surface_644 <- do.call(rbind, grouped_families)
surface_644 <- unique(surface_644)
Heatmap(surface_644)
Heatmap(unique(grouped_families[[1]]))
hallmark_genes <- grouped_families[[1]]
hallmark_genes2 <- unique(grouped_families[[1]])

Heatmap(hallmark_genes2)

all_genes_unique_surface644 <- surface_644
saveRDS(all_genes_unique_surface644, file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/comparative_analysis_filt_biop/clusters_heatmap/comparison_data/all_genes_unique_surface644.rds")


#make a data frame that annotates the overlapping LE from important gene sets
surface_644_annotation <- data.frame(row.names = rownames(surface_644))
n=0
for (i in seq_along(surface_clusters644_matrices)) {
  for (j in seq_along(surface_clusters644_matrices[[i]])) {
      n=n+1
      surface_644_annotation[[n]] <- rownames(surface_644_annotation) %in% rownames(surface_clusters644_matrices[[i]][[j]])
      colnames(surface_644_annotation)[n] <- names(surface_clusters644_matrices[[i]][j])
  }
}
n=0
  
saveRDS(surface_644_annotation, file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/comparative_analysis_filt_biop/clusters_heatmap/comparison_data/surface_644_annotation.rds")


#try to plot heatmap with layers clustered
cluster_positions <- data.frame(row.names = colnames(surface_644))
cluster_positions$biology <- c("m", "m", "m", "c", "m", "s", "s", "c", "s", "m", "s", "c", "m", "s", "c", "s", "m", "c", "s", "s")
col.biology <- list("biology" = c("c" = "blue", "m" = "green3", "s" = "red"))

Heatmap(surface_644, name = '644', 
        column_split = factor(cluster_positions$biology, levels = c("s", "m", "c")),
        row_km = 2,
        cluster_column_slices = FALSE)

#reorder columns for kmeans clustering
kclus <- kmeans(t(surface_644), 2)
kclus$cluster

#custom ordering
split <- factor(paste0("Cluster\n", kclus$cluster), levels=c("Cluster\n2","Cluster\n1"))
#reorder.hmap <- Heatmap(surface_644, column_split=split)
Heatmap(surface_644, name = '644',
        column_split=split, 
        row_km = 4,
        bottom_annotation = HeatmapAnnotation(df = cluster_positions, col = col.biology))

Heatmap(surface_644, name = '644',
        column_split=split, 
        row_km = 1,
        #row_split = surface_644_annotation[1:2],
        #right_annotation = HeatmapAnnotation(df = surface_644_annotation[5], which = "row"),
        bottom_annotation = HeatmapAnnotation(df = cluster_positions, col = col.biology))





  
#set global parameters - "gene set" and "location" ---------
  #collections h.all, c5.bp, c2.cgp, c4.cm
collection="c4.cm"
pathway="MODULE_436"

#for (pathway in c("HALLMARK_INFLAMMATORY_RESPONSE")) {

#Spheroid ---------
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")

  #cluster 6
  # point to correct 'gsea_x' object!!
  fid=which(gsea_6s$c4.cm$pathway==pathway)
  gseaRes=gsea_6s$c4.cm[fid, ]
  nameGset = gseaRes$pathway
  gset = S_subsets[[collection]][[nameGset]]
  ranks = difs_s$'6'
  names(ranks) = rownames(difs_s)
  LE_sphere = sort(unlist(gseaRes$leadingEdge))
  
  #create sub-matrix for plotting
  mat_sphere = dd_spheroid[LE_sphere,]
  #Heatmap(mat_sphere)
  
  
#Organoid ---------
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")

  #cluster 4
  fid=which(gsea_4o$c4.cm$pathway==pathway)
  gseaRes=gsea_4o$c4.cm[fid, ]
  nameGset = gseaRes$pathway
  gset = O_subsets[[collection]][[nameGset]]
  ranks = difs_o$'4'
  names(ranks) = rownames(difs_o)
  LE_organoid = sort(unlist(gseaRes$leadingEdge))
  
  #create sub-matrix for plotting
  mat_organoid=dd_organoid[LE_organoid,]
  #Heatmap(mat_organoid)
  
#Biopsy ---------
setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
  
  #cluster 4
  fid=which(gsea_4b$c4.cm$pathway==pathway)
  gseaRes=gsea_4b$c4.cm[fid, ]
  nameGset = gseaRes$pathway
  gset = B_subsets[[collection]][[nameGset]]
  ranks = difs_b$'4'
  names(ranks) = rownames(difs_b)
  LE_biopsy = sort(unlist(gseaRes$leadingEdge))
  
  #create sub-matrix for plotting
  mat_biopsy=dd_biop[LE_biopsy,]
  #Heatmap(mat_biopsy)
  

#find LE overlap between all biology--------
Common_LE <- intersect(LE_sphere, LE_organoid)
Common_LE <- intersect(Common_LE, LE_biopsy)
  
#create sub-matricies for plotting
mat_sphere=t(dd_spheroid[Common_LE,])
mat_organoid=t(dd_organoid[Common_LE,])
mat_biopsy=t(dd_biop[Common_LE,])
#Heatmap(mat_organoid)
  

#plot as combined biology heatmaps---------
#rename layers by biology and select only surface and center layers
rownames(mat_sphere) <- paste(rownames(mat_sphere), "s", sep = "_")
rownames(mat_organoid) <- paste(rownames(mat_organoid), "o", sep = "_")
rownames(mat_biopsy) <- paste(rownames(mat_biopsy), "b", sep = "_")

#bind new matrix
BIND.MTX <- rbind(mat_sphere, mat_organoid, mat_biopsy)
BIND.MTX <- t(BIND.MTX)
Heatmap(BIND.MTX)
  
#}


# now, i can load all gene sets and there corresponding enrichements :)
# from here, attempt to pull out ALL the genes that i am intersted in for all the layers and all the gene sets - do this by iterating though all interesting gene sets from all layers and growing this into a list of 'all genes'. now, pull out the values from all these genes for each biology and bind these together in a single matrix. Also, make a corresponding data.frame that has all the annotation for each gene set.  Example, GeneSet1: gene_a, gene_b, ... gene_n. Also, force the heatmap to group cells by layers and then group the clusters within each layer by heirarchial clustering.
# make sure to save the matrix with an intutitive name.













#test for merging the union of gene sets rather than the intersection -----------
#create sub-matricies for plotting
union_LE <- union(LE_sphere, LE_organoid)
union_LE <- union(union_LE, LE_biopsy)
#make sure union_LE genes are present in ALL of the biology - if not remove those genes
union_in_all <- Reduce(intersect, list(union_LE, rownames(dd_organoid), rownames(dd_biop), rownames(dd_spheroid)))

#create sub-matricies for plotting
mat_sphere=t(dd_spheroid[union_in_all,])
mat_organoid=t(dd_organoid[union_in_all,])
mat_biopsy=t(dd_biop[union_in_all,])

#rename layers by biology and select only surface and center layers
rownames(mat_sphere) <- paste(rownames(mat_sphere), "s", sep = "_")
rownames(mat_organoid) <- paste(rownames(mat_organoid), "o", sep = "_")
rownames(mat_biopsy) <- paste(rownames(mat_biopsy), "b", sep = "_")

#testbind new union matrix
unionBIND <- rbind(mat_sphere, mat_organoid, mat_biopsy)
unionBIND <- t(unionBIND)
Heatmap(unionBIND)


#end-----


