


surface_clusters640_matrices <- list()
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
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/190926_GSEA_PDX_biopsy23_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/190926_GSEA_PDX_biopsy23_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 0
      fid=which(gsea_0b$h.all$pathway==pathway)
      gseaRes=gsea_0b$h.all[fid, ]
      nameGset = gseaRes$pathway
      gset = B_subsets[[collection]][[nameGset]]
      ranks = difs_b$'0'
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
    surface_clusters640_matrices[[collection]] <- hallmark_matrices
  }
  else if (collection == "c5.bp"){
    GObp_matrices <- list()
    for (pathway in c("GO_INTERFERON_GAMMA_MEDIATED_SIGNALING_PATHWAY",
                      "GO_MONOCYTE_CHEMOTAXIS",
                      "GO_RESPONSE_TO_TYPE_I_INTERFERON",
                      "GO_CYTOKINE_MEDIATED_SIGNALING_PATHWAY",
                      "GO_CELLULAR_RESPONSE_TO_INTERFERON_GAMMA",
                      "GO_INFLAMMATORY_RESPONSE",
                      "GO_REGULATION_OF_ALPHA_BETA_T_CELL_PROLIFERATION",
                      "GO_RESPONSE_TO_INTERFERON_GAMMA",
                      "GO_REGULATION_OF_LEUKOCYTE_MIGRATION",
                      "GO_NEGATIVE_REGULATION_OF_IMMUNE_RESPONSE",
                      "GO_POSITIVE_REGULATION_OF_CYTOKINE_BIOSYNTHETIC_PROCESS",
                      "GO_REGULATION_OF_INFLAMMATORY_RESPONSE",
                      "GO_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_ENDOGENOUS_ANTIGEN",
                      "GO_NEGATIVE_REGULATION_OF_CYTOKINE_PRODUCTION_INVOLVED_IN_IMMUNE_RESPONSE",
                      "GO_POSITIVE_REGULATION_OF_ACUTE_INFLAMMATORY_RESPONSE",
                      "GO_POSITIVE_REGULATION_OF_ALPHA_BETA_T_CELL_PROLIFERATION",
                      "GO_POSITIVE_REGULATION_OF_CHEMOKINE_PRODUCTION",
                      "GO_NEGATIVE_REGULATION_OF_INTERLEUKIN_2_PRODUCTION",
                      "GO_REGULATION_OF_TYPE_2_IMMUNE_RESPONSE",
                      "GO_I_KAPPAB_KINASE_NF_KAPPAB_SIGNALING",
                      "GO_ACUTE_INFLAMMATORY_RESPONSE",
                      "GO_POSITIVE_REGULATION_OF_INFLAMMATORY_RESPONSE",
                      "GO_NEGATIVE_REGULATION_OF_TYPE_I_INTERFERON_PRODUCTION",
                      "GO_REGULATION_OF_CYTOKINE_BIOSYNTHETIC_PROCESS",
                      "GO_NEGATIVE_REGULATION_OF_CYTOKINE_SECRETION")) {
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
        setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/190926_GSEA_PDX_biopsy23_clusters/")
        source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/190926_GSEA_PDX_biopsy23_clusters/Codes/01GSEA_pGSEAfunction.R")

        #cluster 0
        fid=which(gsea_0b$c5.bp$pathway==pathway)
        gseaRes=gsea_0b$c5.bp[fid, ]
        nameGset = gseaRes$pathway
        gset = B_subsets[[collection]][[nameGset]]
        ranks = difs_b$'0'
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
    surface_clusters640_matrices[[collection]] <- GObp_matrices
  }
  else if (collection == "c2.cgp"){
    CGP_matrices <- list()
    for (pathway in c("ZHANG_RESPONSE_TO_IKK_INHIBITOR_AND_TNF_UP",
                      "SANA_TNF_SIGNALING_UP",
                      "SEKI_INFLAMMATORY_RESPONSE_LPS_UP",
                      "HINATA_NFKB_TARGETS_KERATINOCYTE_UP",
                      "ALTEMEIER_RESPONSE_TO_LPS_WITH_MECHANICAL_VENTILATION",
                      "MISSIAGLIA_REGULATED_BY_METHYLATION_UP",
                      "PHONG_TNF_TARGETS_UP",
                      "GAURNIER_PSMD4_TARGETS",
                      "ZWANG_CLASS_3_TRANSIENTLY_INDUCED_BY_EGF",
                      "LIANG_SILENCED_BY_METHYLATION_2",
                      "DIRMEIER_LMP1_RESPONSE_EARLY",
                      "LINDSTEDT_DENDRITIC_CELL_MATURATION_A",
                      "GHANDHI_BYSTANDER_IRRADIATION_UP",
                      "NEMETH_INFLAMMATORY_RESPONSE_LPS_UP",
                      "ICHIBA_GRAFT_VERSUS_HOST_DISEASE_D7_UP",
                      "PHONG_TNF_RESPONSE_VIA_P38_PARTIAL",
                      "HAHTOLA_MYCOSIS_FUNGOIDES_CD4_UP",
                      "WORSCHECH_TUMOR_REJECTION_UP",
                      "BASSO_CD40_SIGNALING_UP",
                      "BROWNE_INTERFERON_RESPONSIVE_GENES")) {
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
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/190926_GSEA_PDX_biopsy23_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/190926_GSEA_PDX_biopsy23_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 0
      fid=which(gsea_0b$c2.cgp$pathway==pathway)
      gseaRes=gsea_0b$c2.cgp[fid, ]
      nameGset = gseaRes$pathway
      gset = B_subsets[[collection]][[nameGset]]
      ranks = difs_b$'0'
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
    surface_clusters640_matrices[[collection]] <- CGP_matrices
  }
  else if (collection == "c2.cp"){
    CanonPW_matrices <- list()
    for (pathway in c("KEGG_GRAFT_VERSUS_HOST_DISEASE",
                      "PID_IL23_PATHWAY",
                      "KEGG_ALLOGRAFT_REJECTION",
                      "KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY",
                      "KEGG_HEMATOPOIETIC_CELL_LINEAGE",
                      "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
                      "KEGG_TYPE_I_DIABETES_MELLITUS",
                      "REACTOME_INTERFERON_GAMMA_SIGNALING",
                      "KEGG_AUTOIMMUNE_THYROID_DISEASE",
                      "ST_TUMOR_NECROSIS_FACTOR_PATHWAY",
                      "BIOCARTA_TNFR2_PATHWAY",
                      "PID_CD40_PATHWAY",
                      "KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION",
                      "BIOCARTA_IL1R_PATHWAY",
                      "PID_IL12_2PATHWAY",
                      "PID_TNF_PATHWAY",
                      "PID_IL27_PATHWAY",
                      "KEGG_CYTOSOLIC_DNA_SENSING_PATHWAY",
                      "PID_INTEGRIN2_PATHWAY",
                      "BIOCARTA_NFKB_PATHWAY",
                      "BIOCARTA_STRESS_PATHWAY",
                      "BIOCARTA_NTHI_PATHWAY",
                      "REACTOME_CHEMOKINE_RECEPTORS_BIND_CHEMOKINES",
                      "BIOCARTA_41BB_PATHWAY",
                      "BIOCARTA_TID_PATHWAY")) {
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
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/190926_GSEA_PDX_biopsy23_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/190926_GSEA_PDX_biopsy23_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 0
      fid=which(gsea_0b$c2.cp$pathway==pathway)
      gseaRes=gsea_0b$c2.cp[fid, ]
      nameGset = gseaRes$pathway
      gset = B_subsets[[collection]][[nameGset]]
      ranks = difs_b$'0'
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
    surface_clusters640_matrices[[collection]] <- CanonPW_matrices
  }
  else if (collection == "c4.cm"){
    CancerMod_matrices <- list()
    for (pathway in c("MODULE_75",
                      "MODULE_46",
                      "MODULE_76",
                      "MODULE_436",
                      "MODULE_345",
                      "MODULE_5",
                      "MODULE_263",
                      "MODULE_300",
                      "MODULE_6",
                      "MODULE_223",
                      "MODULE_108",
                      "MODULE_44",
                      "MODULE_121",
                      "MODULE_488",
                      "MODULE_128",
                      "MODULE_170",
                      "MODULE_79",
                      "MODULE_92",
                      "MODULE_521",
                      "MODULE_340")) {
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
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/190926_GSEA_PDX_biopsy23_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/190926_GSEA_PDX_biopsy23_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 0
      fid=which(gsea_0b$c4.cm$pathway==pathway)
      gseaRes=gsea_0b$c4.cm[fid, ]
      nameGset = gseaRes$pathway
      gset = B_subsets[[collection]][[nameGset]]
      ranks = difs_b$'0'
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
    surface_clusters640_matrices[[collection]] <- CancerMod_matrices
  }
  else if (collection == "c4.cgn"){
    CancerGeneNet_matrices <- list()
    for (pathway in c("GNF2_CD14",
                      "GNF2_PECAM1",
                      "GNF2_FOS",
                      "GNF2_HLA-C",
                      "GNF2_CARD15",
                      "GNF2_CD1D",
                      "GNF2_CD33",
                      "GNF2_MMP11",
                      "GNF2_CDKN1C",
                      "GNF2_CD48",
                      "GNF2_EGFR",
                      "GNF2_CASP1",
                      "GNF2_KISS1",
                      "GNF2_HCK",
                      "GNF2_INPP5D",
                      "GNF2_TIMP2",
                      "GNF2_VAV1",
                      "GNF2_MYL3",
                      "GNF2_MCL1",
                      "CAR_TNFRSF25")) {
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
      setwd("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/190926_GSEA_PDX_biopsy23_clusters/")
      source("/Users/morsedb/Documents/Projects/PSS paint sort sequence/Analysis/biopsy2_3/R_scripts/4th_MERGE_mouses/190926_GSEA_PDX_biopsy23_clusters/Codes/01GSEA_pGSEAfunction.R")
      
      #cluster 0
      fid=which(gsea_0b$c4.cgn$pathway==pathway)
      gseaRes=gsea_0b$c4.cgn[fid, ]
      nameGset = gseaRes$pathway
      gset = B_subsets[[collection]][[nameGset]]
      ranks = difs_b$'0'
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
    surface_clusters640_matrices[[collection]] <- CancerGeneNet_matrices
  } 
  else{
  print("no loop found for this pathway")
  }
}



