#load data from one of the other cluster specific Rscripts (the surface clusters have this 'load' code)
#Run loaded data through loops over all 'interesteing' gene sets from all gene collections -----

middle_clusters422_matrices <- list()
for (collection in c("h.all", "c5.bp", "c2.cgp", "c2.cp", "c4.cm", "c4.cgn")) {
  if (collection == "h.all"){
    hallmark_matrices <- list()
    for (pathway in c("HALLMARK_TGF_BETA_SIGNALING",
                      "HALLMARK_WNT_BETA_CATENIN_SIGNALING",
                      "HALLMARK_HEDGEHOG_SIGNALING",
                      "HALLMARK_ANGIOGENESIS",
                      "HALLMARK_MYOGENESIS",
                      "HALLMARK_NOTCH_SIGNALING")) {
      #Spheroid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 4
      # point to correct 'gsea_x' object!!
      fid=which(gsea_4s$h.all$pathway==pathway)
      gseaRes=gsea_4s$h.all[fid, ]
      nameGset = gseaRes$pathway
      gset = S_subsets[[collection]][[nameGset]]
      ranks = difs_s$'4'
      names(ranks) = rownames(difs_s)
      LE_sphere = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_sphere = dd_spheroid[LE_sphere,]
      #Heatmap(mat_sphere)
      
      
      #Organoid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 2
      fid=which(gsea_2o$h.all$pathway==pathway)
      gseaRes=gsea_2o$h.all[fid, ]
      nameGset = gseaRes$pathway
      gset = O_subsets[[collection]][[nameGset]]
      ranks = difs_o$'2'
      names(ranks) = rownames(difs_o)
      LE_organoid = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_organoid=dd_organoid[LE_organoid,]
      #Heatmap(mat_organoid)
      
      #Biopsy ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 2
      fid=which(gsea_2b$h.all$pathway==pathway)
      gseaRes=gsea_2b$h.all[fid, ]
      nameGset = gseaRes$pathway
      gset = B_subsets[[collection]][[nameGset]]
      ranks = difs_b$'2'
      names(ranks) = rownames(difs_b)
      LE_biopsy = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_biopsy=dd_biop[LE_biopsy,]
      #Heatmap(mat_biopsy)
      
      
      #find LE overlap and combine in a matrix --------
      Common_LE <- intersect(LE_sphere, LE_organoid)
      Common_LE <- intersect(Common_LE, LE_biopsy)
      
      #next if Common_LE is too short
      if (length(Common_LE) < 2) {
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
    middle_clusters422_matrices[[collection]] <- hallmark_matrices
  }
  else if (collection == "c5.bp"){
    GObp_matrices <- list()
    for (pathway in c("GO_POSITIVE_REGULATION_OF_SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR",
                      "GO_PATTERNING_OF_BLOOD_VESSELS",
                      "GO_POSITIVE_REGULATION_OF_KIDNEY_DEVELOPMENT",
                      "GO_POSITIVE_REGULATION_OF_MRNA_METABOLIC_PROCESS",
                      "GO_POSITIVE_T_CELL_SELECTION",
                      "GO_CHONDROCYTE_DEVELOPMENT",
                      "GO_NEGATIVE_REGULATION_OF_CYTOKINE_BIOSYNTHETIC_PROCESS",
                      "GO_REGULATION_OF_KIDNEY_DEVELOPMENT",
                      "GO_NEGATIVE_REGULATION_OF_INTERLEUKIN_6_PRODUCTION",
                      "GO_NEGATIVE_REGULATION_OF_PHOSPHORYLATION",
                      "GO_REGULATION_OF_NOTCH_SIGNALING_PATHWAY",
                      "GO_POSITIVE_REGULATION_OF_MRNA_PROCESSING",
                      "GO_RESPONSE_TO_HEAT",
                      "GO_CANONICAL_WNT_SIGNALING_PATHWAY",
                      "GO_NEGATIVE_REGULATION_OF_OSTEOCLAST_DIFFERENTIATION",
                      "GO_REGULATION_OF_DNA_DAMAGE_RESPONSE_SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR",
                      "GO_MAMMARY_GLAND_LOBULE_DEVELOPMENT",
                      "GO_REGULATION_OF_OSTEOCLAST_DIFFERENTIATION",
                      "GO_CELLULAR_RESPONSE_TO_HEAT",
                      "GO_NEGATIVE_REGULATION_OF_EPITHELIAL_CELL_DIFFERENTIATION",
                      "GO_NEGATIVE_REGULATION_OF_KINASE_ACTIVITY",
                      "GO_POSITIVE_REGULATION_OF_MESENCHYMAL_CELL_PROLIFERATION",
                      "GO_BRANCHING_MORPHOGENESIS_OF_AN_EPITHELIAL_TUBE",
                      "GO_BLOOD_VESSEL_ENDOTHELIAL_CELL_MIGRATION",
                      "GO_NEGATIVE_REGULATION_OF_MAPK_CASCADE",
                      "GO_CELL_FATE_DETERMINATION",
                      "GO_NEGATIVE_REGULATION_OF_PROTEIN_SERINE_THREONINE_KINASE_ACTIVITY",
                      "GO_NEGATIVE_REGULATION_OF_SMOOTH_MUSCLE_CELL_PROLIFERATION",
                      "GO_REGULATION_OF_NEUROBLAST_PROLIFERATION",
                      "GO_MORPHOGENESIS_OF_A_BRANCHING_STRUCTURE")) {
      #Spheroid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 4
      # point to correct 'gsea_x' object!!
      fid=which(gsea_4s$c5.bp$pathway==pathway)
      gseaRes=gsea_4s$c5.bp[fid, ]
      nameGset = gseaRes$pathway
      gset = S_subsets[[collection]][[nameGset]]
      ranks = difs_s$'4'
      names(ranks) = rownames(difs_s)
      LE_sphere = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_sphere = dd_spheroid[LE_sphere,]
      #Heatmap(mat_sphere)
      
      
      #Organoid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 2
      fid=which(gsea_2o$c5.bp$pathway==pathway)
      gseaRes=gsea_2o$c5.bp[fid, ]
      nameGset = gseaRes$pathway
      gset = O_subsets[[collection]][[nameGset]]
      ranks = difs_o$'2'
      names(ranks) = rownames(difs_o)
      LE_organoid = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_organoid=dd_organoid[LE_organoid,]
      #Heatmap(mat_organoid)
      
      #Biopsy ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 2
      fid=which(gsea_2b$c5.bp$pathway==pathway)
      gseaRes=gsea_2b$c5.bp[fid, ]
      nameGset = gseaRes$pathway
      gset = B_subsets[[collection]][[nameGset]]
      ranks = difs_b$'2'
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
    middle_clusters422_matrices[[collection]] <- GObp_matrices
  }
  else if (collection == "c2.cgp"){
    CGP_matrices <- list()
    for (pathway in c("NAGASHIMA_NRG1_SIGNALING_UP",
                      "SMIRNOV_RESPONSE_TO_IR_6HR_DN",
                      "NAGASHIMA_EGF_SIGNALING_UP",
                      "BROCKE_APOPTOSIS_REVERSED_BY_IL6",
                      "CROONQUIST_NRAS_VS_STROMAL_STIMULATION_DN",
                      "BURTON_ADIPOGENESIS_1",
                      "TENEDINI_MEGAKARYOCYTE_MARKERS",
                      "BAKER_HEMATOPOIESIS_STAT3_TARGETS",
                      "GENTILE_UV_HIGH_DOSE_DN",
                      "AMIT_DELAYED_EARLY_GENES",
                      "BEGUM_TARGETS_OF_PAX3_FOXO1_FUSION_DN",
                      "GENTILE_UV_RESPONSE_CLUSTER_D2",
                      "PLASARI_TGFB1_TARGETS_1HR_UP",
                      "WATTEL_AUTONOMOUS_THYROID_ADENOMA_DN",
                      "SARTIPY_NORMAL_AT_INSULIN_RESISTANCE_UP",
                      "BILBAN_B_CLL_LPL_DN",
                      "SIMBULAN_UV_RESPONSE_NORMAL_DN",
                      "GARGALOVIC_RESPONSE_TO_OXIDIZED_PHOSPHOLIPIDS_BLUE_UP",
                      "MARSON_FOXP3_TARGETS_STIMULATED_UP",
                      "WANG_LSD1_TARGETS_DN",
                      "YORDY_RECIPROCAL_REGULATION_BY_ETS1_AND_SP100_DN",
                      "WANG_METHYLATED_IN_BREAST_CANCER",
                      "YIH_RESPONSE_TO_ARSENITE_C1",
                      "PEART_HDAC_PROLIFERATION_CLUSTER_UP",
                      "UEDA_CENTRAL_CLOCK",
                      "DAZARD_RESPONSE_TO_UV_SCC_DN",
                      "KORKOLA_TERATOMA",
                      "PARK_HSC_MARKERS",
                      "SHIN_B_CELL_LYMPHOMA_CLUSTER_2",
                      "BHATI_G2M_ARREST_BY_2METHOXYESTRADIOL_UP")) {
      #Spheroid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 4
      # point to correct 'gsea_x' object!!
      fid=which(gsea_4s$c2.cgp$pathway==pathway)
      gseaRes=gsea_4s$c2.cgp[fid, ]
      nameGset = gseaRes$pathway
      gset = S_subsets[[collection]][[nameGset]]
      ranks = difs_s$'4'
      names(ranks) = rownames(difs_s)
      LE_sphere = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_sphere = dd_spheroid[LE_sphere,]
      #Heatmap(mat_sphere)
      
      
      #Organoid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 2
      fid=which(gsea_2o$c2.cgp$pathway==pathway)
      gseaRes=gsea_2o$c2.cgp[fid, ]
      nameGset = gseaRes$pathway
      gset = O_subsets[[collection]][[nameGset]]
      ranks = difs_o$'2'
      names(ranks) = rownames(difs_o)
      LE_organoid = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_organoid=dd_organoid[LE_organoid,]
      #Heatmap(mat_organoid)
      
      #Biopsy ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 2
      fid=which(gsea_2b$c2.cgp$pathway==pathway)
      gseaRes=gsea_2b$c2.cgp[fid, ]
      nameGset = gseaRes$pathway
      gset = B_subsets[[collection]][[nameGset]]
      ranks = difs_b$'2'
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
    middle_clusters422_matrices[[collection]] <- CGP_matrices
  }
  else if (collection == "c2.cp"){
    CanonPW_matrices <- list()
    for (pathway in c("PID_ATF2_PATHWAY",
                      "KEGG_HEMATOPOIETIC_CELL_LINEAGE",
                      "KEGG_MAPK_SIGNALING_PATHWAY",
                      "KEGG_TGF_BETA_SIGNALING_PATHWAY",
                      "PID_FRA_PATHWAY",
                      "PID_BETA_CATENIN_NUC_PATHWAY",
                      "KEGG_P53_SIGNALING_PATHWAY",
                      "KEGG_HYPERTROPHIC_CARDIOMYOPATHY_HCM",
                      "PID_AMB2_NEUTROPHILS_PATHWAY",
                      "PID_INTEGRIN5_PATHWAY",
                      "PID_INTEGRIN3_PATHWAY",
                      "PID_IL4_2PATHWAY",
                      "PID_GLYPICAN_1PATHWAY",
                      "KEGG_JAK_STAT_SIGNALING_PATHWAY",
                      "PID_SMAD2_3NUCLEAR_PATHWAY",
                      "PID_CMYB_PATHWAY",
                      "PID_TGFBR_PATHWAY",
                      "PID_PDGFRA_PATHWAY",
                      "KEGG_ADIPOCYTOKINE_SIGNALING_PATHWAY",
                      "REACTOME_TRANSCRIPTIONAL_REGULATION_OF_WHITE_ADIPOCYTE_DIFFERENTIATION",
                      "PID_P53_DOWNSTREAM_PATHWAY",
                      "BIOCARTA_IL6_PATHWAY",
                      "PID_IL3_PATHWAY",
                      "ST_GA12_PATHWAY",
                      "PID_MAPK_TRK_PATHWAY",
                      "PID_RHOA_PATHWAY",
                      "PID_HNF3B_PATHWAY",
                      "PID_FGF_PATHWAY",
                      "PID_ERBB1_DOWNSTREAM_PATHWAY",
                      "PID_IL6_7_PATHWAY")) {
      #Spheroid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 4
      # point to correct 'gsea_x' object!!
      fid=which(gsea_4s$c2.cp$pathway==pathway)
      gseaRes=gsea_4s$c2.cp[fid, ]
      nameGset = gseaRes$pathway
      gset = S_subsets[[collection]][[nameGset]]
      ranks = difs_s$'4'
      names(ranks) = rownames(difs_s)
      LE_sphere = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_sphere = dd_spheroid[LE_sphere,]
      #Heatmap(mat_sphere)
      
      
      #Organoid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 2
      fid=which(gsea_2o$c2.cp$pathway==pathway)
      gseaRes=gsea_2o$c2.cp[fid, ]
      nameGset = gseaRes$pathway
      gset = O_subsets[[collection]][[nameGset]]
      ranks = difs_o$'2'
      names(ranks) = rownames(difs_o)
      LE_organoid = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_organoid=dd_organoid[LE_organoid,]
      #Heatmap(mat_organoid)
      
      #Biopsy ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 2
      fid=which(gsea_2b$c2.cp$pathway==pathway)
      gseaRes=gsea_2b$c2.cp[fid, ]
      nameGset = gseaRes$pathway
      gset = B_subsets[[collection]][[nameGset]]
      ranks = difs_b$'2'
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
    middle_clusters422_matrices[[collection]] <- CanonPW_matrices
  }
  else if (collection == "c4.cm"){
    CancerMod_matrices <- list()
    for (pathway in c("MODULE_123",
                      "MODULE_358",
                      "MODULE_182",
                      "MODULE_97",
                      "MODULE_525",
                      "MODULE_254",
                      "MODULE_261",
                      "MODULE_444",
                      "MODULE_287",
                      "MODULE_362",
                      "MODULE_47",
                      "MODULE_416",
                      "MODULE_108",
                      "MODULE_213",
                      "MODULE_38",
                      "MODULE_229",
                      "MODULE_375",
                      "MODULE_411",
                      "MODULE_341",
                      "MODULE_573",
                      "MODULE_145",
                      "MODULE_24",
                      "MODULE_157",
                      "MODULE_280",
                      "MODULE_521",
                      "MODULE_188",
                      "MODULE_174",
                      "MODULE_94",
                      "MODULE_448",
                      "MODULE_489")) {
      #Spheroid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 4
      # point to correct 'gsea_x' object!!
      fid=which(gsea_4s$c4.cm$pathway==pathway)
      gseaRes=gsea_4s$c4.cm[fid, ]
      nameGset = gseaRes$pathway
      gset = S_subsets[[collection]][[nameGset]]
      ranks = difs_s$'4'
      names(ranks) = rownames(difs_s)
      LE_sphere = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_sphere = dd_spheroid[LE_sphere,]
      #Heatmap(mat_sphere)
      
      
      #Organoid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 2
      fid=which(gsea_2o$c4.cm$pathway==pathway)
      gseaRes=gsea_2o$c4.cm[fid, ]
      nameGset = gseaRes$pathway
      gset = O_subsets[[collection]][[nameGset]]
      ranks = difs_o$'2'
      names(ranks) = rownames(difs_o)
      LE_organoid = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_organoid=dd_organoid[LE_organoid,]
      #Heatmap(mat_organoid)
      
      #Biopsy ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 2
      fid=which(gsea_2b$c4.cm$pathway==pathway)
      gseaRes=gsea_2b$c4.cm[fid, ]
      nameGset = gseaRes$pathway
      gset = B_subsets[[collection]][[nameGset]]
      ranks = difs_b$'2'
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
    middle_clusters422_matrices[[collection]] <- CancerMod_matrices
  }
  else if (collection == "c4.cgn"){
    CancerGeneNet_matrices <- list()
    for (pathway in c("GNF2_PTX3",
                      "MORF_CAMK4",
                      "MORF_IL4",
                      "GCM_ERBB2IP",
                      "MORF_SUPT3H",
                      "MORF_MAGEA8",
                      "MORF_CTSB",
                      "MORF_MAP2K7",
                      "MORF_RAD51L3",
                      "MORF_CDH4",
                      "MORF_ATF2",
                      "MORF_ERCC4",
                      "GNF2_MATK",
                      "GCM_PTPRD",
                      "MORF_MDM2",
                      "MORF_MLLT10",
                      "MORF_PAX7",
                      "MORF_FSHR",
                      "MORF_MAGEA9",
                      "MORF_FOSL1",
                      "MORF_TNFRSF6",
                      "MORF_LTK",
                      "GNF2_ITGAL",
                      "MORF_BCL2",
                      "MORF_PTPRB",
                      "GNF2_PTPN4",
                      "GNF2_RAB7L1",
                      "MORF_NOS2A",
                      "MORF_IL16",
                      "MORF_BMPR2")) {
      #Spheroid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 4
      # point to correct 'gsea_x' object!!
      fid=which(gsea_4s$c4.cgn$pathway==pathway)
      gseaRes=gsea_4s$c4.cgn[fid, ]
      nameGset = gseaRes$pathway
      gset = S_subsets[[collection]][[nameGset]]
      ranks = difs_s$'4'
      names(ranks) = rownames(difs_s)
      LE_sphere = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_sphere = dd_spheroid[LE_sphere,]
      #Heatmap(mat_sphere)
      
      
      #Organoid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 2
      fid=which(gsea_2o$c4.cgn$pathway==pathway)
      gseaRes=gsea_2o$c4.cgn[fid, ]
      nameGset = gseaRes$pathway
      gset = O_subsets[[collection]][[nameGset]]
      ranks = difs_o$'2'
      names(ranks) = rownames(difs_o)
      LE_organoid = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_organoid=dd_organoid[LE_organoid,]
      #Heatmap(mat_organoid)
      
      #Biopsy ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 2
      fid=which(gsea_2b$c4.cgn$pathway==pathway)
      gseaRes=gsea_2b$c4.cgn[fid, ]
      nameGset = gseaRes$pathway
      gset = B_subsets[[collection]][[nameGset]]
      ranks = difs_b$'2'
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
    middle_clusters422_matrices[[collection]] <- CancerGeneNet_matrices
  } 
  else{
    print("no loop found for this pathway")
  }
}

#save list of lists as RDS -------
saveRDS(middle_clusters422_matrices, file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/comparative_analysis_filt_biop/clusters_heatmap/comparison_data/middle422.rds")

#unlist and bind into a large matrix called middle_422
grouped_families <- list()
for (i in seq_along(middle_clusters422_matrices)) {
  grouped_families[[i]] <- do.call(rbind, middle_clusters422_matrices[[i]])
}
middle_422 <- do.call(rbind, grouped_families)
middle_422 <- unique(middle_422)
Heatmap(middle_422)
Heatmap(unique(grouped_families[[1]]))
hallmark_genes <- grouped_families[[1]]
hallmark_genes2 <- unique(grouped_families[[1]])

Heatmap(hallmark_genes2)

all_genes_unique_middle422 <- middle_422
saveRDS(all_genes_unique_middle422, file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/comparative_analysis_filt_biop/clusters_heatmap/comparison_data/all_genes_unique_middle422.rds")


#make a data frame that annotates the overlapping LE from important gene sets
middle_422_annotation <- data.frame(row.names = rownames(middle_422))
n=0
for (i in seq_along(middle_clusters422_matrices)) {
  for (j in seq_along(middle_clusters422_matrices[[i]])) {
    n=n+1
    middle_422_annotation[[n]] <- rownames(middle_422_annotation) %in% rownames(middle_clusters422_matrices[[i]][[j]])
    colnames(middle_422_annotation)[n] <- names(middle_clusters422_matrices[[i]][j])
  }
}
n=0

saveRDS(middle_422_annotation, file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/comparative_analysis_filt_biop/clusters_heatmap/comparison_data/middle_422_annotation.rds")


#try to plot heatmap with layers clustered
cluster_positions <- data.frame(row.names = colnames(middle_422))
cluster_positions$biology <- c("m", "m", "m", "c", "m", "s", "s", "c", "s", "m", "s", "c", "m", "s", "c", "s", "m", "c", "s", "s")
col.biology <- list("biology" = c("c" = "blue", "m" = "green", "s" = "red"))

Heatmap(hallmark_genes2, name = '422', 
        column_split = factor(cluster_positions$biology, levels = c("s", "m", "c")),
        row_km = 4,
        cluster_column_slices = FALSE)

#reorder columns for kmeans clustering
kclus <- kmeans(t(middle_422), 2)
kclus$cluster

#custom ordering
split <- factor(paste0("Cluster\n", kclus$cluster), levels=c("Cluster\n2","Cluster\n1"))
#reorder.hmap <- Heatmap(middle_422, column_split=split)
Heatmap(middle_422, name = '422',
        column_split=split, 
        row_km = 4,
        bottom_annotation = HeatmapAnnotation(df = cluster_positions, col = col.biology))

Heatmap(middle_422, name = '422',
        #column_split=split, 
        column_km = 4,
        row_km = 1,
        #row_split = middle_422_annotation[1:2],
        #right_annotation = HeatmapAnnotation(df = middle_422_annotation[5], which = "row"),
        bottom_annotation = HeatmapAnnotation(df = cluster_positions, col = col.biology))



#end-----

