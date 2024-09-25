rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)


display_names    <- read.csv("preprocessed/MSigDB_display_names.csv")
rownames(display_names) <- display_names$STANDARD_NAME


Core_paths0 <- read.csv("tables/phospho_path_clean/PhosphoLosmapimod_vs_I45F_diffPathways.csv", row.names = NULL)



Core_paths0$I45FvsControl <-  -log10(Core_paths0$P.Value)

Core_paths0 <- Core_paths0[(order(Core_paths0$logFC,decreasing = T)),]

if(any(str_length(Core_paths0$Pathway.Description) > 80)){
    longInds <- str_length(Core_paths0$Pathway.Description) > 80
    Core_paths0[longInds,]$Pathway.Description <- stringr::str_sub(Core_paths0[longInds,]$Pathway.Description,1,70)
    Core_paths0[longInds,]$Pathway.Description <- paste0(Core_paths0[longInds,]$Pathway.Description,"*") 
}




Core_paths0_up <- Core_paths0 %>% 
    dplyr::filter(.,Direction..AD.vs.Control. == "UP") %>% dplyr::select(.,Pathway.Description, I45FvsControl,logFC)







colnames(Core_paths0_up)[2:3] <- c("signedP","EffectSize")



p1 <- ggplot(Core_paths0_up[30:1,], aes( fct_inorder(Pathway.Description), EffectSize, fill = signedP)) +
    geom_bar(stat = "identity",position = "dodge", width = 0.8) +
    labs(y = "Effect Size (Losmapimod vs A\u03B242-H)", fill = "-log10(adj. p)") + ylim(c(0,7))+
    theme_bw() +
    scale_fill_gradient(low = "#fb6a4a" ,high = "#67000d")+
    coord_flip() + 
    theme(panel.grid.major.x =  element_blank(),
          panel.grid.minor.x =  element_blank(),
          panel.border = element_rect(colour = "grey80", fill=NA, size=0),
          axis.title.x = element_text(size = 16),# element_text(size = 20, angle = 0,vjust = 0.5, hjust = 0.5, face = "bold"),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 16,angle = 0,vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 16,angle = 0,vjust = 0.5, hjust = 1),
          legend.text = element_text(size = 12),
          #legend.title = element_blank()
          )

# cairo_pdf("figures/Phospho_barplot_up_I45vsCtr.pdf",width = 12,height = 12)
# p1
# dev.off()



Core_paths0_down <- Core_paths0 %>% 
    dplyr::filter(.,Direction..AD.vs.Control. == "DOWN") %>% dplyr::select(.,Pathway.Description, I45FvsControl,logFC)



colnames(Core_paths0_down)[2:3] <- c("signedP","EffectSize")

nItems <- nrow(Core_paths0_down)


p2 <- ggplot(Core_paths0_down[(nItems-30):nItems,], aes( fct_inorder(Pathway.Description), EffectSize, fill = signedP)) +
    geom_bar(stat = "identity",position = "dodge", width = 0.8) +
    labs(y = "Effect Size (Losmapimod vs A\u03B242-H)", fill = "-log10(adj. p)") +
    ylim(c(-7,0))+
    theme_bw() +
    scale_fill_gradient(low = "#6baed6" ,high = "#08306b")+
    coord_flip() + 
    theme(panel.grid.major.x =  element_blank(),
          panel.grid.minor.x =  element_blank(),
          panel.border = element_rect(colour = "grey80", fill=NA, size=0),
          axis.title.x = element_text(size = 16),# element_text(size = 20, angle = 0,vjust = 0.5, hjust = 0.5, face = "bold"),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 16,angle = 0,vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 16,angle = 0,vjust = 0.5, hjust = 1),
          legend.text = element_text(size = 12),
          #legend.title = element_blank()
    )

library(patchwork)
cairo_pdf("figures/Phospho_losma_barplot_down_I45vsCtr.pdf",width = 24,height = 8)
p1+p2
dev.off()


