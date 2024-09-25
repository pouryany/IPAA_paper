rm(list = ls())

library(stringr)
library(ComplexHeatmap)
library(dplyr)


## Shared pathways between transcriptomcis and phosphoproteomics

library(tidyverse)
all.kat.files <- list.files("tables/phospho_path_clean/", full.names = T, recursive = F, pattern = "\\.csv|\\.txt")
all.kat.files <- all.kat.files[grep("Losma",all.kat.files,invert = F)]
all.P38.files <- read.csv(all.kat.files)

gene_pathways <- read.csv("tables/DE_pathways_clean/I45FvsG2B2_diffPathways.csv")

Core_paths0 <-  merge(all.P38.files, gene_pathways,
                      by = c("MSigDB.name","Pathway.Description"), suffixes = c("_Losma","_NoTreatment"),all = F)


Core_paths1 <- Core_paths0[Core_paths0$Direction..AD.vs.Control._Losma != "NONE" &
                               Core_paths0$Direction..AD.vs.Control._NoTreatment != "NONE" ,  ]   

 
table(Core_paths1$Direction..AD.vs.Control._phospho,Core_paths1$Direction..AD.vs.Control._gene)

Core_paths1 <- Core_paths1[Core_paths1$Direction..AD.vs.Control._Losma !=
                               Core_paths1$Direction..AD.vs.Control._NoTreatment,  ]  

write.csv(Core_paths1,"tables/opposite_Losma_phospho_gene_pathways_Abeta42.csv")


Core_paths0 <-  merge(all.P38.files, gene_pathways,
                      by = c("MSigDB.name","Pathway.Description"), suffixes = c("_Lomsa","_NoTreament"),all = F)


Core_paths0$phospho_pathDir <- paste0(Core_paths0$MSigDB.name,
                                      Core_paths0$Direction..AD.vs.Control._phospho)

Core_paths0$gene_pathDir <- paste0(Core_paths0$MSigDB.name,
                                      Core_paths0$Direction..AD.vs.Control._gene)

Core_paths0 <- Core_paths0[Core_paths0$Direction..AD.vs.Control._gene != "NONE" &
                               Core_paths0$Direction..AD.vs.Control._phospho != "NONE" ,  ]   
