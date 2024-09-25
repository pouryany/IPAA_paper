rm(list = ls())

library(stringr)
library(ComplexHeatmap)
library(dplyr)


## Shared pathways between transcriptomcis and phosphoproteomics

library(tidyverse)
all.kat.files <- list.files("tables/phospho_path_clean", full.names = T, recursive = F, pattern = "\\.csv|\\.txt")
all.kat.files <- all.kat.files[grep("Losma",all.kat.files,invert = T)]
all.P38.files <- read.csv(all.kat.files)

gene_pathways <- read.csv("tables/DE_pathways_clean/I45FvsG2B2_diffPathways.csv")

Core_paths0 <-  merge(all.P38.files, gene_pathways,
                      by = c("MSigDB.name","Pathway.Description"), suffixes = c("_phospho","_gene"),all = F)


Core_paths1 <- Core_paths0[Core_paths0$Direction..AD.vs.Control._gene != "NONE" &
                               Core_paths0$Direction..AD.vs.Control._phospho != "NONE" ,  ]   

 
table(Core_paths1$Direction..AD.vs.Control._phospho,Core_paths1$Direction..AD.vs.Control._gene)

Core_paths1 <- Core_paths1[Core_paths1$Direction..AD.vs.Control._gene ==
                               Core_paths1$Direction..AD.vs.Control._phospho,  ]  

write.csv(Core_paths1,"tables/shared_phospho_gene_pathways_Abeta42.csv")

Core_paths0$phospho <- sign(Core_paths0$logFC_phospho) * -log(Core_paths0$P.Value_phospho)
Core_paths0$genes <- sign(Core_paths0$logFC_gene) * -log(Core_paths0$P.Value_gene)

Core_paths0  <- Core_paths0[,c("Pathway.Description","phospho", "genes")]
rownames(Core_paths0) <- Core_paths0$MSigDB.name


library(ComplexHeatmap)
ht1 <- ComplexHeatmap::Heatmap(as.matrix(Core_paths0[,2:3]), 
                               border_gp = gpar(col = "black", lty = 1),
                               width =  unit(3,"cm"),
                               col = circlize::colorRamp2(c(20,10,2,0,-2,-10, -20), 
                                                          colors = c("#b2182b",
                                                                     "#ef8a62",
                                                                     "#fddbc7",
                                                                     "#f7f7f7",
                                                                     "#d1e5f0",
                                                                     "#67a9cf",
                                                                     "#2166ac")),
                               # name = "-log10(pval)",
                               clustering_distance_rows = "canberra",
                               name = "Dysregulation score",
                               cluster_columns = F,
                               show_row_names = F,
                               show_column_dend = F, 
                               show_row_dend = F,
                               heatmap_legend_param =
                                   list(direction = "horizontal",
                                        title_gp = gpar(fontsize = 10,
                                                        fontface = "bold"))
)


Core_paths0 <-  merge(all.P38.files, gene_pathways,
                      by = c("MSigDB.name","Pathway.Description"), suffixes = c("_phospho","_gene"),all = F)


Core_paths0$phospho_pathDir <- paste0(Core_paths0$MSigDB.name,
                                      Core_paths0$Direction..AD.vs.Control._phospho)

Core_paths0$gene_pathDir <- paste0(Core_paths0$MSigDB.name,
                                      Core_paths0$Direction..AD.vs.Control._gene)

Core_paths0 <- Core_paths0[Core_paths0$Direction..AD.vs.Control._gene != "NONE" &
                               Core_paths0$Direction..AD.vs.Control._phospho != "NONE" ,  ]   

# 
# library(VennDiagram)
# library(RColorBrewer)
# myCol <- brewer.pal(4, "Pastel2")
# 
# futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
# 
# a <- venn.diagram(
#     x = list( Core_paths0$phospho_pathDir,
#               Core_paths0$gene_pathDir),
#     category.names = c("Phospho-proteome" , "Transcriptome"),
#     filename = NULL,
#     output=F,
#     # Numbers,
#     col  = c("#FFB000", "#80B1D3"),#myCol[1:3],
#     fill = c("#FFB000", "#80B1D3"),#myCol[1:3],
#     cex = 3,
#     fontface = "plain",
#     fontfamily = "sans",
#     cat.cex = 3,
#     # Set names
#     cat.fontface = "plain",
#     cat.default.pos = "outer",
#     cat.fontfamily = "sans",
#     main.cex = 2.0,
#     main.fontfamily = "sans",
#     cat.dist = c(0.035, 0.035) ,
#     cat.pos = c(-30, 30) # Modified
#     
# )
# 
# 
# ggsave("figures/phospho_gene_assays_abeta42vsCtrl_venn.PDF",plot =a, width = 12, height = 12)
# 

library(stringr)
library(dplyr)
library(ggpubr)





fix_header_bigTab <- function(big.tab.list){
    
    big_names <- names(big.tab.list)
    
    big_names <- gsub(pattern = ".*/|_QC.*|_diffPathways.*|2021.*","", big_names)
    
    
    big_names<- gsub("vsG2B2","",big_names)
    big_names <- gsub("I45FvsI47F","A\u03B242-H / A\u03B240-H",big_names)
    big_names <- gsub("I45F","A\u03B242-H",big_names)
    big_names <- gsub("I47F","A\u03B240-H",big_names)
    big_names <- gsub("v.*","",big_names)
    big_names <- gsub("_"," ",big_names)
    big_names <- gsub(" PHG","",big_names)
    big_names <- gsub(" DLPFC","",big_names)
    
    big_names <- paste0(big_names, " vs Ctrl")
    big_ind   <- grep("A\u03B242-H / A\u03B240-H",big_names)
    big_names[big_ind] <- "A\u03B242-H vs A\u03B240-H"
    return(big_names)
}


plot_agreement <- function(temp, name_tag, legend_pos = "none"){
    
    nums <- table(temp$Direction.x , temp$Direction.y)
    chi_test <- chisq.test(nums)
    chi_test <- signif(chi_test$p.value,3)
    agreement    <- diag(nums)
    disagreement <- sum(nums) - sum(diag(nums))
    agreement    <- as.data.frame(t(agreement))
    agreement$disagreement <- -1 * disagreement
    agreement <- t(agreement)
    agreement <- tibble::rownames_to_column(as.data.frame(agreement),"Direction")
    agreement$Direction <- c("Down in Aβ42-H (Agree)",
                             "Up in Aβ42-H (Agree)",
                             "Disagree")
    
    p1 <- ggplot(agreement, aes(Direction,V1, fill = Direction)) +
        geom_bar(position= "stack", stat="identity",width = 0.6) +
        theme_minimal() + geom_hline(yintercept = 0, color = "grey30") +
        scale_colour_manual(values = c("#8dd3c7", "#80b1d3", "#fb8072")) +
        scale_fill_manual(values = c("#8073ac", "#4393c3", "#d6604d")) +
        labs(x = name_tag) +
        ylim(c(-40,60)) + 
        geom_text(aes(label = V1, y = V1 + (10 * sign(V1-0.1))), size= 4) + 
        theme(panel.grid.major.x =  element_blank(),
              panel.grid.minor.x =  element_blank(),
              panel.border = element_rect(colour = "grey80", fill=NA, size=2),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size = 20,
                                          angle = 0,
                                          vjust = 0.5, hjust = 1,
                                          face = "plain"),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              legend.text = element_text(face = "plain"),
              legend.position = legend_pos) +
        annotate("text",x=  1, y = 20, hjust = 0,
                 label = paste0("p-value = ", chi_test), size = 6) + 
        coord_flip()
    
    return(p1)
}




Core_paths0 <-  merge(all.P38.files, gene_pathways,
                      by = c("MSigDB.name","Pathway.Description"),all = F)

Core_paths0$Direction.x <- Core_paths0$Direction..AD.vs.Control..x
Core_paths0$Direction.y <- Core_paths0$Direction..AD.vs.Control..y

Core_paths0 <- Core_paths0[Core_paths0$Direction..AD.vs.Control..x != "NONE" &
                               Core_paths0$Direction..AD.vs.Control..y != "NONE" ,  ]   




p1 <-  plot_agreement(Core_paths0, name_tag = NULL, legend_pos = "left")




cairo_pdf("figures/phospho_agreement.pdf", width = 7, height = 1.5, onefile = T)
p1

dev.off()
# 
