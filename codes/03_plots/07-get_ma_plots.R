rm(list = ls())
library(stringr)
library(dplyr)
library(ggpubr)
temp <- list.files("tables/DE_genes_filtered_Mjd/",recursive = T,full.names = T)
temp <- grep("vs_G|F_vs", temp,value = T)

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
    t.tab    <- read.csv(i)
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
    

    colnames(t.tab2)
    ma.data <- t.tab2[,c(4,3,6)]
    
    colnames(ma.data) <- c("baseMeanLog2","log2FoldChange","padj")
    rownames(ma.data) <- t.tab2$ID
    
    name.main <- paste(name.tag[1])
    
    p1 <- ggmaplot(ma.data, 
             main = name.main,
             fdr = 0.05, fc = 2, size = 2,
             alpha = 0.5,
             palette = c("#B31B21", "#1465AC", "darkgray"),
             genenames = rownames(ma.data),
             legend = "top", top = 0,
             font.label = c("bold", 25),
             font.legend = c(15,"bold"),
             font.main = "bold",
             ggtheme = ggplot2::theme_bw(),
             font.x = 24,
             font.y = 24,  
             font.submain = 24,
             font.caption = 24,
             font.title = 24,
             font.subtitle = 40,
             font.tickslab = 20,
             xlim = c(-4,16),
             ylim = c(-13,13))
    
    plot.list[[base.ind]] <- p1
    
    
    
    }


library(gridExtra)

this.plot <- do.call("grid.arrange",c(plot.list,ncol = 3))
ggsave(paste0("figures/S02_b_2023_MA_genes.jpg"),
       plot = this.plot, width = 20, height = 16)

