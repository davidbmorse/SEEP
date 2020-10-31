#cell cycle study on Spheroid layers
#david Jan 15, 2019

#load a seurat object containing genes of interest
#for examples, load "/Users/morsedb/Documents/Experiments/2017/4July_August/20170711_inDrop_spheroid_layers/data_analysis/seurat scripts/SEURAT_CIOS_Oct2018dm.R"

library(tidyr)
library(dplyr)

#make a gene list containing cell cycle genes:

cell_cycle_genes <- read.table("/Users/morsedb/Documents/Experiments/2017/4July_August/20170711_inDrop_spheroid_layers/data_analysis/NCATS_analysis_171213/regev_lab_cell_cycle_genes.txt")

cell_cycle_gene_list <- cell_cycle_genes$V1
s.genes <- cell_cycle_gene_list[1:43]
g2m.genes <- cell_cycle_gene_list[44:97]


#Assign a cell cycle score 
spheroid_CIOS <- CellCycleScoring(object = spheroid_CIOS, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(x = spheroid_CIOS@meta.data)

RidgePlot(object = spheroid_CIOS, features.plot = c("PCNA", "TOP2A", "MCM6", "MKI67"), nCol = 2)

#set primary identity of cells to layers
spheroid_CIOS = SetAllIdent(spheroid_CIOS, "orig.ident")
#reset primary identity of cells to clusters
spheroid_CIOS = SetAllIdent(spheroid_CIOS, "res.0.4")


layers_cc <- tibble(spheroid_CIOS@ident, spheroid_CIOS@meta.data$Phase, .name_repair = c("check_unique", "unique", "universal", "minimal"))

layer_list <- spheroid_CIOS@ident
phase_list <- spheroid_CIOS@meta.data$Phase
layers_cc <- tibble(layer_list, phase_list)

phase_intensity <- 
  layers_cc %>% 
  group_by(layer_list, phase_list) %>% 
  summarise(n())

write.csv(phase_intensity, file = "/Users/morsedb/Desktop/phase_intensity")

cc_intensity <- read.csv(file = "/Users/morsedb/Documents/Experiments/2017/4July_August/20170711_inDrop_spheroid_layers/data_analysis/NCATS_analysis_171213/cc_intensity_CIOS.csv")

#spread this data

cc_phase <-
  spread(phase_intensity, phase_list, 'n()')

cc <-
  spread(cc_intensity, phase_list, 'n..')

# Stacked Percent phase
ggplot(cc_phase, aes(fill=c(G1, G2M), y='n()', x = layer_list)) + 
  geom_bar( stat="identity")

DoHeatmap(spheroid_CIOS, genes.use = s.genes, draw.line = TRUE, slim.col.label = TRUE, use.scaled = FALSE, remove.key = TRUE)



# create a dataset
specie=c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
condition=rep(c("normal" , "stress" , "Nitrogen") , 4)
value=abs(rnorm(12 , 0 , 15))
data=data.frame(specie,condition,value)


cc_phase %>% 
  gather(key, value, -layer_list) %>% 
  ggplot(aes(layer_list, value, fill = key)) +
  geom_col(position="fill") +
  theme_bw() +
  labs(fill = "Cell Cycle Phases")

categories <- c("G1 & S phases", "G2/M phases")
cols <- c("#E61515", "#5121E0")
# more colors: "#5121E0" #105415

cc %>% 
  gather(key, value, -layer_list) %>% 
  ggplot(aes(layer_list, value, fill = key)) +
  geom_col(position="fill") +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "top") +
  scale_fill_manual(labels = categories, values = cols) +
  labs(fill = "Cell Cycle Phases") +
  xlab("layer") +
  ylab("fraction of cells")

