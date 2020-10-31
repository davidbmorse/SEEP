#load data from one of the other cluster specific Rscripts (the surface clusters have this 'load' code)
#Run loaded data through loops over all 'interesteing' gene sets from all gene collections -----

center_clusters363_matrices <- list()
for (collection in c("h.all", "c5.bp", "c2.cgp", "c2.cp", "c4.cm", "c4.cgn")) {
  if (collection == "h.all"){
    hallmark_matrices <- list()
    for (pathway in c("HALLMARK_PROTEIN_SECRETION",
                      "HALLMARK_UV_RESPONSE_DN",
                      "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
                      "HALLMARK_ANDROGEN_RESPONSE",
                      "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
                      "HALLMARK_HEME_METABOLISM",
                      "HALLMARK_MTORC1_SIGNALING",
                      "HALLMARK_MYOGENESIS",
                      "HALLMARK_KRAS_SIGNALING_DN")) {
      #Spheroid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 1
      # point to correct 'gsea_x' object!!
      fid=which(gsea_3s$h.all$pathway==pathway)
      gseaRes=gsea_3s$h.all[fid, ]
      nameGset = gseaRes$pathway
      gset = S_subsets[[collection]][[nameGset]]
      ranks = difs_s$'3'
      names(ranks) = rownames(difs_s)
      LE_sphere = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_sphere = dd_spheroid[LE_sphere,]
      #Heatmap(mat_sphere)
      
      
      #Organoid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 2
      fid=which(gsea_6o$h.all$pathway==pathway)
      gseaRes=gsea_6o$h.all[fid, ]
      nameGset = gseaRes$pathway
      gset = O_subsets[[collection]][[nameGset]]
      ranks = difs_o$'6'
      names(ranks) = rownames(difs_o)
      LE_organoid = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_organoid=dd_organoid[LE_organoid,]
      #Heatmap(mat_organoid)
      
      #Biopsy ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 4
      fid=which(gsea_3b$h.all$pathway==pathway)
      gseaRes=gsea_3b$h.all[fid, ]
      nameGset = gseaRes$pathway
      gset = B_subsets[[collection]][[nameGset]]
      ranks = difs_b$'3'
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
      
      #rename layers by biology and select only surface and center layers
      rownames(mat_sphere) <- paste(rownames(mat_sphere), "s", sep = "_")
      rownames(mat_organoid) <- paste(rownames(mat_organoid), "o", sep = "_")
      rownames(mat_biopsy) <- paste(rownames(mat_biopsy), "b", sep = "_")
      
      #bind new matrix
      BIND.MTX <- rbind(mat_sphere, mat_organoid, mat_biopsy)
      BIND.MTX <- t(BIND.MTX)
      hallmark_matrices[[pathway]] <- BIND.MTX
    }
    center_clusters363_matrices[[collection]] <- hallmark_matrices
  }
  else if (collection == "c5.bp"){
    GObp_matrices <- list()
    for (pathway in c("GO_GLOMERULAR_EPITHELIUM_DEVELOPMENT",
                      "GO_EPITHELIAL_CELL_DIFFERENTIATION_INVOLVED_IN_KIDNEY_DEVELOPMENT",
                      "GO_REGULATION_OF_CIRCADIAN_RHYTHM",
                      "GO_CIRCADIAN_REGULATION_OF_GENE_EXPRESSION",
                      "GO_CENTRAL_NERVOUS_SYSTEM_NEURON_AXONOGENESIS",
                      "GO_PROXIMAL_DISTAL_PATTERN_FORMATION",
                      "GO_NLS_BEARING_PROTEIN_IMPORT_INTO_NUCLEUS",
                      "GO_DEFINITIVE_HEMOPOIESIS",
                      "GO_RESPONSE_TO_ESTRADIOL",
                      "GO_REGULATION_OF_CELL_FATE_COMMITMENT",
                      "GO_DNA_DOUBLE_STRAND_BREAK_PROCESSING",
                      "GO_APPENDAGE_DEVELOPMENT",
                      "GO_REPRODUCTIVE_BEHAVIOR",
                      "GO_ORGANELLE_TRANSPORT_ALONG_MICROTUBULE",
                      "GO_CELLULAR_RESPONSE_TO_ALCOHOL",
                      "GO_ORGAN_GROWTH",
                      "GO_RNA_POLYADENYLATION",
                      "GO_MRNA_3_END_PROCESSING",
                      "GO_PLATELET_MORPHOGENESIS",
                      "GO_MRNA_PROCESSING",
                      "GO_GOLGI_VESICLE_TRANSPORT",
                      "GO_MULTICELLULAR_ORGANISM_GROWTH",
                      "GO_CENTRAL_NERVOUS_SYSTEM_PROJECTION_NEURON_AXONOGENESIS",
                      "GO_ER_TO_GOLGI_VESICLE_MEDIATED_TRANSPORT",
                      "GO_VESICLE_COATING")) {
      #Spheroid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 1
      # point to correct 'gsea_x' object!!
      fid=which(gsea_3s$c5.bp$pathway==pathway)
      gseaRes=gsea_3s$c5.bp[fid, ]
      nameGset = gseaRes$pathway
      gset = S_subsets[[collection]][[nameGset]]
      ranks = difs_s$'3'
      names(ranks) = rownames(difs_s)
      LE_sphere = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_sphere = dd_spheroid[LE_sphere,]
      #Heatmap(mat_sphere)
      
      
      #Organoid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 2
      fid=which(gsea_6o$c5.bp$pathway==pathway)
      gseaRes=gsea_6o$c5.bp[fid, ]
      nameGset = gseaRes$pathway
      gset = O_subsets[[collection]][[nameGset]]
      ranks = difs_o$'6'
      names(ranks) = rownames(difs_o)
      LE_organoid = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_organoid=dd_organoid[LE_organoid,]
      #Heatmap(mat_organoid)
      
      #Biopsy ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 4
      fid=which(gsea_3b$c5.bp$pathway==pathway)
      gseaRes=gsea_3b$c5.bp[fid, ]
      nameGset = gseaRes$pathway
      gset = B_subsets[[collection]][[nameGset]]
      ranks = difs_b$'3'
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
    center_clusters363_matrices[[collection]] <- GObp_matrices
  }
  else if (collection == "c2.cgp"){
    CGP_matrices <- list()
    for (pathway in c("SHEN_SMARCA2_TARGETS_UP",
                      "IKEDA_MIR133_TARGETS_UP",
                      "GABRIELY_MIR21_TARGETS",
                      "JOHNSTONE_PARVB_TARGETS_1_DN",
                      "DACOSTA_UV_RESPONSE_VIA_ERCC3_COMMON_DN",
                      "SENGUPTA_NASOPHARYNGEAL_CARCINOMA_WITH_LMP1_UP",
                      "GENTILE_UV_LOW_DOSE_DN",
                      "RAMALHO_STEMNESS_UP",
                      "GINESTIER_BREAST_CANCER_20Q13_AMPLIFICATION_UP",
                      "CAIRO_HEPATOBLASTOMA_UP",
                      "MARTORIATI_MDM4_TARGETS_FETAL_LIVER_DN",
                      "MILI_PSEUDOPODIA_HAPTOTAXIS_UP",
                      "DAZARD_RESPONSE_TO_UV_NHEK_DN",
                      "ENK_UV_RESPONSE_KERATINOCYTE_DN",
                      "STARK_PREFRONTAL_CORTEX_22Q11_DELETION_UP",
                      "KANG_DOXORUBICIN_RESISTANCE_DN",
                      "YANAGIHARA_ESX1_TARGETS",
                      "GINESTIER_BREAST_CANCER_ZNF217_AMPLIFIED_UP",
                      "DAZARD_UV_RESPONSE_CLUSTER_G6",
                      "JOHNSTONE_PARVB_TARGETS_2_DN")) {
      #Spheroid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 1
      # point to correct 'gsea_x' object!!
      fid=which(gsea_3s$c2.cgp$pathway==pathway)
      gseaRes=gsea_3s$c2.cgp[fid, ]
      nameGset = gseaRes$pathway
      gset = S_subsets[[collection]][[nameGset]]
      ranks = difs_s$'3'
      names(ranks) = rownames(difs_s)
      LE_sphere = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_sphere = dd_spheroid[LE_sphere,]
      #Heatmap(mat_sphere)
      
      
      #Organoid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 2
      fid=which(gsea_6o$c2.cgp$pathway==pathway)
      gseaRes=gsea_6o$c2.cgp[fid, ]
      nameGset = gseaRes$pathway
      gset = O_subsets[[collection]][[nameGset]]
      ranks = difs_o$'6'
      names(ranks) = rownames(difs_o)
      LE_organoid = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_organoid=dd_organoid[LE_organoid,]
      #Heatmap(mat_organoid)
      
      #Biopsy ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 4
      fid=which(gsea_3b$c2.cgp$pathway==pathway)
      gseaRes=gsea_3b$c2.cgp[fid, ]
      nameGset = gseaRes$pathway
      gset = B_subsets[[collection]][[nameGset]]
      ranks = difs_b$'3'
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
    center_clusters363_matrices[[collection]] <- CGP_matrices
  }
  else if (collection == "c2.cp"){
    CanonPW_matrices <- list()
    for (pathway in c("REACTOME_YAP1_AND_WWTR1_TAZ_STIMULATED_GENE_EXPRESSION",
                      "BIOCARTA_CTCF_PATHWAY",
                      "PID_CDC42_PATHWAY",
                      "SIG_INSULIN_RECEPTOR_PATHWAY_IN_CARDIAC_MYOCYTES",
                      "KEGG_LONG_TERM_POTENTIATION",
                      "SIG_CHEMOTAXIS",
                      "PID_ATM_PATHWAY",
                      "PID_NECTIN_PATHWAY",
                      "PID_NCADHERIN_PATHWAY",
                      "KEGG_ADHERENS_JUNCTION",
                      "KEGG_ASCORBATE_AND_ALDARATE_METABOLISM",
                      "KEGG_INSULIN_SIGNALING_PATHWAY",
                      "PID_ERBB4_PATHWAY",
                      "PID_TRKR_PATHWAY",
                      "REACTOME_INSULIN_RECEPTOR_SIGNALLING_CASCADE",
                      "PID_ECADHERIN_STABILIZATION_PATHWAY",
                      "PID_HNF3A_PATHWAY",
                      "BIOCARTA_PTEN_PATHWAY",
                      "PID_ERBB1_RECEPTOR_PROXIMAL_PATHWAY",
                      "REACTOME_DOWNSTREAM_SIGNALING_OF_ACTIVATED_FGFR")) {
      #Spheroid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 1
      # point to correct 'gsea_x' object!!
      fid=which(gsea_3s$c2.cp$pathway==pathway)
      gseaRes=gsea_3s$c2.cp[fid, ]
      nameGset = gseaRes$pathway
      gset = S_subsets[[collection]][[nameGset]]
      ranks = difs_s$'3'
      names(ranks) = rownames(difs_s)
      LE_sphere = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_sphere = dd_spheroid[LE_sphere,]
      #Heatmap(mat_sphere)
      
      
      #Organoid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 2
      fid=which(gsea_6o$c2.cp$pathway==pathway)
      gseaRes=gsea_6o$c2.cp[fid, ]
      nameGset = gseaRes$pathway
      gset = O_subsets[[collection]][[nameGset]]
      ranks = difs_o$'6'
      names(ranks) = rownames(difs_o)
      LE_organoid = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_organoid=dd_organoid[LE_organoid,]
      #Heatmap(mat_organoid)
      
      #Biopsy ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 4
      fid=which(gsea_3b$c2.cp$pathway==pathway)
      gseaRes=gsea_3b$c2.cp[fid, ]
      nameGset = gseaRes$pathway
      gset = B_subsets[[collection]][[nameGset]]
      ranks = difs_b$'3'
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
    center_clusters363_matrices[[collection]] <- CanonPW_matrices
  }
  else if (collection == "c4.cm"){
    CancerMod_matrices <- list()
    for (pathway in c("MODULE_239",
                      "MODULE_323",
                      "MODULE_36",
                      "MODULE_97",
                      "MODULE_160",
                      "MODULE_277",
                      "MODULE_182",
                      "MODULE_567",
                      "MODULE_261",
                      "MODULE_133",
                      "MODULE_352",
                      "MODULE_110",
                      "MODULE_35",
                      "MODULE_361",
                      "MODULE_331",
                      "MODULE_372",
                      "MODULE_453",
                      "MODULE_72",
                      "MODULE_159",
                      "MODULE_432")) {
      #Spheroid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 1
      # point to correct 'gsea_x' object!!
      fid=which(gsea_3s$c4.cm$pathway==pathway)
      gseaRes=gsea_3s$c4.cm[fid, ]
      nameGset = gseaRes$pathway
      gset = S_subsets[[collection]][[nameGset]]
      ranks = difs_s$'3'
      names(ranks) = rownames(difs_s)
      LE_sphere = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_sphere = dd_spheroid[LE_sphere,]
      #Heatmap(mat_sphere)
      
      
      #Organoid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 2
      fid=which(gsea_6o$c4.cm$pathway==pathway)
      gseaRes=gsea_6o$c4.cm[fid, ]
      nameGset = gseaRes$pathway
      gset = O_subsets[[collection]][[nameGset]]
      ranks = difs_o$'6'
      names(ranks) = rownames(difs_o)
      LE_organoid = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_organoid=dd_organoid[LE_organoid,]
      #Heatmap(mat_organoid)
      
      #Biopsy ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 4
      fid=which(gsea_3b$c4.cm$pathway==pathway)
      gseaRes=gsea_3b$c4.cm[fid, ]
      nameGset = gseaRes$pathway
      gset = B_subsets[[collection]][[nameGset]]
      ranks = difs_b$'3'
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
    center_clusters363_matrices[[collection]] <- CancerMod_matrices
  }
  else if (collection == "c4.cgn"){
    CancerGeneNet_matrices <- list()
    for (pathway in c("GCM_RAB10",
                      "GNF2_KPNB1",
                      "GNF2_ELAC2",
                      "GCM_DDX5",
                      "GCM_BAG5",
                      "GCM_HBP1",
                      "GCM_BMPR2",
                      "MORF_SP3",
                      "GNF2_PAK2",
                      "GNF2_DDX5",
                      "GCM_AQP4",
                      "GCM_MAP4K4",
                      "GNF2_BNIP2",
                      "GCM_RAD21",
                      "GCM_RAN",
                      "GCM_CALM1",
                      "GCM_SUFU",
                      "GCM_MYST2",
                      "GCM_CRKL",
                      "MORF_DNMT1")) {
      #Spheroid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/PEO1_spheroids/1712_GSEA_scRNA_spheroid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 1
      # point to correct 'gsea_x' object!!
      fid=which(gsea_3s$c4.cgn$pathway==pathway)
      gseaRes=gsea_3s$c4.cgn[fid, ]
      nameGset = gseaRes$pathway
      gset = S_subsets[[collection]][[nameGset]]
      ranks = difs_s$'3'
      names(ranks) = rownames(difs_s)
      LE_sphere = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_sphere = dd_spheroid[LE_sphere,]
      #Heatmap(mat_sphere)
      
      
      #Organoid ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/organoidsComb_190226/R_scripts/190926_GSEA_scRNA_organoid_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 2
      fid=which(gsea_6o$c4.cgn$pathway==pathway)
      gseaRes=gsea_6o$c4.cgn[fid, ]
      nameGset = gseaRes$pathway
      gset = O_subsets[[collection]][[nameGset]]
      ranks = difs_o$'6'
      names(ranks) = rownames(difs_o)
      LE_organoid = sort(unlist(gseaRes$leadingEdge))
      
      #create sub-matrix for plotting
      mat_organoid=dd_organoid[LE_organoid,]
      #Heatmap(mat_organoid)
      
      #Biopsy ---------
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/SubtractMouseReads/GSEA_biopsy_clusters_042320/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 4
      fid=which(gsea_3b$c4.cgn$pathway==pathway)
      gseaRes=gsea_3b$c4.cgn[fid, ]
      nameGset = gseaRes$pathway
      gset = B_subsets[[collection]][[nameGset]]
      ranks = difs_b$'3'
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
    center_clusters363_matrices[[collection]] <- CancerGeneNet_matrices
  } 
  else{
    print("no loop found for this pathway")
  }
}

#save list of lists as RDS -------
saveRDS(center_clusters363_matrices, file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/comparative_analysis_filt_biop/clusters_heatmap/comparison_data/center363.rds")

#unlist and bind into a large matrix called center_363
grouped_families <- list()
for (i in seq_along(center_clusters363_matrices)) {
  grouped_families[[i]] <- do.call(rbind, center_clusters363_matrices[[i]])
}
center_363 <- do.call(rbind, grouped_families)
center_363 <- unique(center_363)
Heatmap(center_363)
Heatmap(unique(grouped_families[[1]]))
hallmark_genes <- grouped_families[[1]]
hallmark_genes2 <- unique(grouped_families[[1]])

Heatmap(hallmark_genes2)

all_genes_unique_center363 <- center_363
saveRDS(all_genes_unique_center363, file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/comparative_analysis_filt_biop/clusters_heatmap/comparison_data/all_genes_unique_center363.rds")


#make a data frame that annotates the overlapping LE from important gene sets
center_363_annotation <- data.frame(row.names = rownames(center_363))
n=0
for (i in seq_along(center_clusters363_matrices)) {
  for (j in seq_along(center_clusters363_matrices[[i]])) {
    n=n+1
    center_363_annotation[[n]] <- rownames(center_363_annotation) %in% rownames(center_clusters363_matrices[[i]][[j]])
    colnames(center_363_annotation)[n] <- names(center_clusters363_matrices[[i]][j])
  }
}
n=0

saveRDS(center_363_annotation, file = "/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/comparative_analysis_filt_biop/clusters_heatmap/comparison_data/center_363_annotation.rds")


#try to plot heatmap with layers clustered
cluster_positions <- data.frame(row.names = colnames(center_363))
cluster_positions$biology <- c("m", "m", "m", "c", "m", "s", "s", "c", "s", "m", "s", "c", "m", "s", "c", "s", "m", "c", "s", "s")
col.biology <- list("biology" = c("c" = "blue", "m" = "green", "s" = "red"))

Heatmap(hallmark_genes2, name = '363', 
        column_split = factor(cluster_positions$biology, levels = c("s", "m", "c")),
        row_km = 4,
        cluster_column_slices = FALSE)

#reorder columns for kmeans clustering
kclus <- kmeans(t(center_363), 2)
kclus$cluster

#custom ordering
split <- factor(paste0("Cluster\n", kclus$cluster), levels=c("Cluster\n2","Cluster\n1"))
#reorder.hmap <- Heatmap(center_363, column_split=split)
Heatmap(center_363, name = '363',
        column_split=split, 
        row_km = 4,
        bottom_annotation = HeatmapAnnotation(df = cluster_positions, col = col.biology))

Heatmap(center_363, name = '363',
        #column_split=split, 
        column_km = 4,
        row_km = 1,
        #row_split = center_363_annotation[1:2],
        #right_annotation = HeatmapAnnotation(df = center_363_annotation[5], which = "row"),
        bottom_annotation = HeatmapAnnotation(df = cluster_positions, col = col.biology))



#end-----