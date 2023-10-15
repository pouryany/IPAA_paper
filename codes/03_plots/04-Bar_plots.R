
rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)


Core_paths0 <- read.csv("tables/shared_pathways/01_shared_pathways_I45FvsI47_MSBB_Mayo.csv", row.names = NULL)

Core_paths0$I45FvsI47F <- sign(Core_paths0$logFC.x) * -log(Core_paths0$P.Value.x)
Core_paths0$Mayo <- sign(Core_paths0$logFC.y) * -log(Core_paths0$P.Value.y)
Core_paths0$MSBB <- sign(Core_paths0$logFC) * -log(Core_paths0$P.Value)

Core_paths0 <- Core_paths0[order(Core_paths0$Mayo),]
Core_paths0_up <- Core_paths0 %>% dplyr::filter(.,Direction == "UP") %>% dplyr::select(.,X, I45FvsI47F, Mayo, MSBB)

colnames(Core_paths0_up)[2] <- "A\u03B242-H vs A\u03B240-H"

Core_paths0_up <- pivot_longer(Core_paths0_up,c("Mayo","MSBB", "A\u03B242-H vs A\u03B240-H"))
library(viridis)


p1 <- ggplot(Core_paths0_up, aes( fct_inorder(X), value, fill = fct_inorder(name))) +
    geom_bar(stat = "identity",position = "dodge", width = 0.8) +
    labs(y = "Dysregulation score") + ylim(c(0,31))+
    theme_bw() +
    coord_flip() + 
    scale_fill_manual( values=c("#FFB000", 
                                 "#80B1D3",
                                 "#beaed4"))+
    theme(panel.grid.major.x =  element_blank(),
          panel.grid.minor.x =  element_blank(),
          panel.border = element_rect(colour = "grey80", fill=NA, size=0),
          axis.title.x = element_blank(),# element_text(size = 20, angle = 0,vjust = 0.5, hjust = 0.5, face = "bold"),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 16,angle = 0,vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 10,angle = 0,vjust = 0.5, hjust = 1),
          legend.text = element_text(face = "plain"),legend.position = "none",
          legend.title = element_blank())
  


Core_paths0 <- read.csv("tables/shared_pathways/01_shared_pathways_I45FvsI47_MSBB_Mayo.csv", row.names = NULL)

Core_paths0$I45FvsI47F <- sign(Core_paths0$logFC.x) * -log(Core_paths0$P.Value.x)
Core_paths0$Mayo <- sign(Core_paths0$logFC.y) * -log(Core_paths0$P.Value.y)
Core_paths0$MSBB <- sign(Core_paths0$logFC) * -log(Core_paths0$P.Value)

Core_paths0 <- Core_paths0[order(Core_paths0$Mayo),]
Core_paths0_down <- Core_paths0 %>% dplyr::filter(.,Direction == "DOWN") %>% dplyr::select(.,X, I45FvsI47F, Mayo, MSBB)

colnames(Core_paths0_down)[2] <- "A\u03B242-H vs A\u03B240-H"


Core_paths0_down <- Core_paths0_down

Core_paths0_down <- pivot_longer(Core_paths0_down,c("Mayo","MSBB", "A\u03B242-H vs A\u03B240-H"))
library(viridis)




if(any(str_length(Core_paths0_down$X) > 70)){
  longInds <- str_length(Core_paths0_down$X) > 70
  Core_paths0_down[longInds,]$X <- stringr::str_sub(Core_paths0_down[longInds,]$X,1,70)
  Core_paths0_down[longInds,]$X <- paste0(Core_paths0_down[longInds,]$X,"*") 
}







p2 <- ggplot(Core_paths0_down, aes( fct_inorder(X), value, fill = fct_inorder(name))) +
  geom_bar(stat = "identity",position = "dodge", width = 0.8) +
  labs(y = "Dysregulation score") + ylim(c(-31,0))+
  theme_bw() +
  coord_flip() + 
  scale_fill_manual( values=c("#FFB000", 
                              "#80B1D3",
                              "#beaed4"))+
  theme(panel.grid.major.x =  element_blank(),
        panel.grid.minor.x =  element_blank(),
        panel.border = element_rect(colour = "grey80", fill=NA, size=0),
        axis.title.x = element_blank(),#element_text(size = 20, angle = 0,vjust = 0.5, hjust = 0.5, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.x =  element_text(size = 16,angle = 0,vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 10,
                                   angle = 0,
                                   vjust = 0.5, hjust = 1),
        legend.text = element_text(face = "bold"), legend.position = "left",
        legend.title = element_blank()) + 
  scale_x_discrete(position = "top") 

library(patchwork)



cairo_pdf("figures/04_b_shared_barplot.pdf",width = 14,height = 12)   
p1/p2 + patchwork::plot_layout(heights = c(1,2.3))
dev.off()


