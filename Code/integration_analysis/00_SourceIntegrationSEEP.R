rm(list = ls())
setwd("C://PROJECTS/P2022/SEEP_Manuscript")

#cat(sprintf("source('%s', echo=TRUE)\n", list.files("Codes", full=T, rec=T)[-1]))

# Analysis ====

# Combine separate SEEP (no integration)

# Integrate SEEP
source('Codes/01A_SelectSEEP.R', echo = TRUE)    #OK
source('Codes/01B_IntegrateSEEP.R', echo = TRUE) #OK
source('Codes/01C_DaSEEP.R', echo = TRUE)        #OK
source('Codes/01D_MarkersSEEP.R', echo = TRUE)   #OK
source('Codes/01E_OraSEEP.R', echo = TRUE)       #OK

source('Codes/01G_SeparateSEEP.R', echo = TRUE)  #OK

# Transfer SEEP to SS2
source('Codes/02A_ProcessSS2.R', echo = TRUE)    #OK
source('Codes/02B_TransferSS2.R', echo = TRUE)   #OK

# Transfer SEEP to 10X
source('Codes/03A_process10X.R', echo = TRUE)    #OK
source('Codes/03B_Transfer10X.R', echo = TRUE)   #OK

# Calculate pathway scores
source('Codes/01F_OraScoresSEEP.R', echo = TRUE) #OK
source('Codes/02C_OraScoresSS2.R', echo = TRUE)  #OK
source('Codes/03C_OraScores10X.R', echo = TRUE)  #OK

# Manuscript Figures ====

## fig 6
source('Codes/Fig6ABC.R', echo = TRUE)  # OK
source('Codes/Fig6E.R', echo=TRUE)      # OK
source('Codes/Fig6D.R', echo=TRUE)      # OK

## fig 7
source('Codes/Fig7ABC.R', echo=TRUE) # OK
source('Codes/Fig7DEF.R', echo=TRUE) # OK
source('Codes/Fig7GHI.R', echo=TRUE) # OK

# Supplement Figures ====
source('Codes/SFigXAD.R', echo=TRUE) # OK
source('Codes/SFigXE.R', echo=TRUE)  # OK










