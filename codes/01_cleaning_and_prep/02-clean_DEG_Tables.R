rm(list = ls())
library(stringr)
library(dplyr)
library(ggpubr)
temp <- list.files("preprocessed/cell_models/",recursive = T,full.names = T)
temp <- grep(".txt", temp,value = T)

meta.labels <- list()
big.exp     <- list()

fdr.thresh   <- 0.05
logFC.thresh <- 1

plot.list    <- list()
big.tab.list <- list()
for(i in temp) {
    base.ind <- which(temp == i)
    name.tag <- unlist(str_split(i,"/"))
    name.tag <- tail(name.tag,1)
    name.tag <- unlist(str_split(name.tag,pattern = "_"))
    name.tag <- gsub(".txt","",name.tag)
    t.tab    <- read.csv(i,sep = "\t")
    t.tab2   <- t.tab %>% dplyr::select(.,Name,ID,logFC,logCPM,PValue,FDR) %>%
                    arrange(.,PValue) %>%
                    mutate(Direction = logFC/abs(logFC),
                           Direction = factor(Direction, c(-1,1),
                                              c('-1' = 'DOWN', '1' = 'UP')),
                           Direction = as.character(Direction))
                    
    t.tab2$Direction[t.tab2$FDR > fdr.thresh | abs(t.tab2$logFC) < logFC.thresh] = 'NONE'
    
    t.readout       <- as.data.frame.matrix(t(table(t.tab2$Direction)))
    t.readout       <- data.frame(t.readout)
    
    t.assay   <- unlist(strsplit(i, "/"))
    t.assay   <- tail(t.assay,1) 
    t.assay   <- gsub(".txt|_CPM_log2FC_FDR_complete_table","",t.assay)

    
    t.readout$assay <- t.assay
    big.tab.list    <- rbind(big.tab.list,t.readout)
    
    t.address   <- unlist(strsplit(i, "/"))
    t.address   <- tail(t.address,1) 
    t.address   <- paste0("tables/DE_genes_filtered_Mjd/",t.address)
    t.address   <- gsub(".txt",".csv",t.address)
    
    write.csv(t.tab2,t.address,row.names = F)

    }







