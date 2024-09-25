rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)


mean_kea_up    <- readxl::read_xlsx("data/Supplementary Table_KEA_up.xlsx",sheet = 2)

mean_kea_up$Protein <- gsub("\\*","",mean_kea_up$Protein) 

p1 <- ggplot(mean_kea_up[30:1,], aes( fct_inorder(Protein), `Mean rank`)) +
    geom_bar(stat = "identity",position = "dodge", width = 0.8, fill = "#67000d") +
    labs(y = "Mean Rank (KEA3)",x = "") + ylim(c(0,80))+
    theme_bw() +
    # scale_fill_gradient(low = "#fb6a4a" ,high = "#67000d")+
    coord_flip()+ 
    theme(panel.grid.major.x =  element_blank(),
          panel.grid.minor.x =  element_blank(),
          panel.border = element_rect(colour = "grey80", fill=NA, size=0),
          axis.title.x = element_text(size = 16),# element_text(size = 20, angle = 0,vjust = 0.5, hjust = 0.5, face = "bold"),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 16,angle = 0,vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 16,angle = 0,vjust = 0.5, hjust = 1),
          legend.position = "none")
    

Cheng_kea_up    <- readxl::read_xlsx("data/Supplementary Table_KEA_up.xlsx",sheet = 3)

Cheng_kea_up$Protein <- gsub("\\*","",Cheng_kea_up$Protein) 

p2 <- ggplot(Cheng_kea_up[30:1,], aes( fct_inorder(Protein), -log10(FDR), fill = "#67000d" )) +
    geom_bar(stat = "identity",position = "dodge", width = 0.8,fill = "#67000d") +
    labs(y = "-log10(FDR) - ChengKSIN",x = "") + ylim(c(0,10))+
    theme_bw() +
    # scale_fill_gradient(low = "#fb6a4a" ,high = "#67000d")+
    coord_flip()+ 
    theme(panel.grid.major.x =  element_blank(),
          panel.grid.minor.x =  element_blank(),
          panel.border = element_rect(colour = "grey80", fill=NA, size=0),
          axis.title.x = element_text(size = 16),# element_text(size = 20, angle = 0,vjust = 0.5, hjust = 0.5, face = "bold"),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 16,angle = 0,vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 16,angle = 0,vjust = 0.5, hjust = 1),
          legend.position = "none")




PhospD_kea_up    <- readxl::read_xlsx("data/Supplementary Table_KEA_up.xlsx",sheet = 4)

PhospD_kea_up$Protein <- gsub("\\*","",PhospD_kea_up$Protein) 

p3 <- ggplot(PhospD_kea_up[30:1,], aes( fct_inorder(Protein), -log10(FDR), fill = "#67000d" )) +
    geom_bar(stat = "identity",position = "dodge", width = 0.8 , fill = "#67000d") +
    labs(y = "-log10(FDR) - PhosDAll ",x = "") + ylim(c(0,10))+
    theme_bw() +
    # scale_fill_gradient(low = "#fb6a4a" ,high = "#67000d")+
    coord_flip()+ 
    theme(panel.grid.major.x =  element_blank(),
          panel.grid.minor.x =  element_blank(),
          panel.border = element_rect(colour = "grey80", fill=NA, size=0),
          axis.title.x = element_text(size = 16),# element_text(size = 20, angle = 0,vjust = 0.5, hjust = 0.5, face = "bold"),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 16,angle = 0,vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 16,angle = 0,vjust = 0.5, hjust = 1),
          legend.position = "none")




PTMsigDB_kea_up    <- readxl::read_xlsx("data/Supplementary Table_KEA_up.xlsx",sheet = 5)

PTMsigDB_kea_up$Protein <- gsub("\\*","",PTMsigDB_kea_up$Protein) 

p4 <- ggplot(PTMsigDB_kea_up[30:1,], aes( fct_inorder(Protein), -log10(FDR), fill = "#67000d" )) +
    geom_bar(stat = "identity",position = "dodge", width = 0.8, fill = "#67000d") +
    labs(y = "-log10(FDR) - PTMSigDB",x = "") + ylim(c(0,10))+
    theme_bw() +
    # scale_fill_gradient(low = "#fb6a4a" ,high = "#67000d")+
    coord_flip()+ 
    theme(panel.grid.major.x =  element_blank(),
          panel.grid.minor.x =  element_blank(),
          panel.border = element_rect(colour = "grey80", fill=NA, size=0),
          axis.title.x = element_text(size = 16),# element_text(size = 20, angle = 0,vjust = 0.5, hjust = 0.5, face = "bold"),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 16,angle = 0,vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 16,angle = 0,vjust = 0.5, hjust = 1),
          legend.position = "none")

library(patchwork)

cairo_pdf("figures/Phospho_KEA3_barplot_up_I45vsCtr.pdf",width = 24,height = 7)
 p1 + p2 + p3 + p4 + plot_layout(nrow = 1)
dev.off()




## Downregulated phosphoproteins


mean_kea_down    <- readxl::read_xlsx("data/Supplementary Table_KEA_down.xlsx",sheet = 2)

mean_kea_down$Protein <- gsub("\\*","",mean_kea_down$Protein) 

p1 <- ggplot(mean_kea_down[30:1,], aes( fct_inorder(Protein), `Mean rank`)) +
    geom_bar(stat = "identity",position = "dodge", width = 0.8, fill = "#08306b" ) +
    labs(y = "Mean Rank (KEA3)",x = "") + ylim(c(0,80))+
    theme_bw() +
    # scale_fill_gradient(low = "#fb6a4a" ,high = "#67000d")+
    coord_flip()+ 
    theme(panel.grid.major.x =  element_blank(),
          panel.grid.minor.x =  element_blank(),
          panel.border = element_rect(colour = "grey80", fill=NA, size=0),
          axis.title.x = element_text(size = 16),# element_text(size = 20, angle = 0,vjust = 0.5, hjust = 0.5, face = "bold"),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 16,angle = 0,vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 16,angle = 0,vjust = 0.5, hjust = 1),
          legend.position = "none")


Cheng_kea_down    <- readxl::read_xlsx("data/Supplementary Table_KEA_down.xlsx",sheet = 3)

Cheng_kea_down$Protein <- gsub("\\*","",Cheng_kea_down$Protein) 

p2 <- ggplot(Cheng_kea_down[30:1,], aes( fct_inorder(Protein), -log10(FDR), fill = "#67000d" )) +
    geom_bar(stat = "identity",position = "dodge", width = 0.8, fill = "#08306b" ) +
    labs(y = "-log10(FDR) - ChengKSIN",x = "") + ylim(c(0,10))+
    theme_bw() +
    # scale_fill_gradient(low = "#fb6a4a" ,high = "#67000d")+
    coord_flip()+ 
    theme(panel.grid.major.x =  element_blank(),
          panel.grid.minor.x =  element_blank(),
          panel.border = element_rect(colour = "grey80", fill=NA, size=0),
          axis.title.x = element_text(size = 16),# element_text(size = 20, angle = 0,vjust = 0.5, hjust = 0.5, face = "bold"),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 16,angle = 0,vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 16,angle = 0,vjust = 0.5, hjust = 1),
          legend.position = "none")




PhospD_kea_down    <- readxl::read_xlsx("data/Supplementary Table_KEA_down.xlsx",sheet = 4)

PhospD_kea_down$Protein <- gsub("\\*","",PhospD_kea_down$Protein) 

p3 <- ggplot(PhospD_kea_down[30:1,], aes( fct_inorder(Protein), -log10(FDR), fill = "#67000d" )) +
    geom_bar(stat = "identity",position = "dodge", width = 0.8 , fill = "#08306b" ) +
    labs(y = "-log10(FDR) - PhosDAll ",x = "") + ylim(c(0,14))+
    theme_bw() +
    # scale_fill_gradient(low = "#fb6a4a" ,high = "#67000d")+
    coord_flip()+ 
    theme(panel.grid.major.x =  element_blank(),
          panel.grid.minor.x =  element_blank(),
          panel.border = element_rect(colour = "grey80", fill=NA, size=0),
          axis.title.x = element_text(size = 16),# element_text(size = 20, angle = 0,vjust = 0.5, hjust = 0.5, face = "bold"),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 16,angle = 0,vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 16,angle = 0,vjust = 0.5, hjust = 1),
          legend.position = "none")




PTMsigDB_kea_down    <- readxl::read_xlsx("data/Supplementary Table_KEA_down.xlsx",sheet = 5)

PTMsigDB_kea_down$Protein <- gsub("\\*","",PTMsigDB_kea_down$Protein) 

p4 <- ggplot(PTMsigDB_kea_down[30:1,], aes( fct_inorder(Protein), -log10(FDR), fill = "#67000d" )) +
    geom_bar(stat = "identity",position = "dodge", width = 0.8, fill = "#08306b" ) +
    labs(y = "-log10(FDR) - PTMSigDB",x = "") + ylim(c(0,10))+
    theme_bw() +
    # scale_fill_gradient(low = "#fb6a4a" ,high = "#67000d")+
    coord_flip()+ 
    theme(panel.grid.major.x =  element_blank(),
          panel.grid.minor.x =  element_blank(),
          panel.border = element_rect(colour = "grey80", fill=NA, size=0),
          axis.title.x = element_text(size = 16),# element_text(size = 20, angle = 0,vjust = 0.5, hjust = 0.5, face = "bold"),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 16,angle = 0,vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 16,angle = 0,vjust = 0.5, hjust = 1),
          legend.position = "none")

library(patchwork)

cairo_pdf("figures/Phospho_KEA3_barplot_down_I45vsCtr.pdf",width = 24,height = 7)
p1 + p2 + p3 + p4 + plot_layout(nrow = 1)
dev.off()



