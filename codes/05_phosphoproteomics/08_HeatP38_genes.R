rm(list = ls())

library(stringr)
library(ComplexHeatmap)
library(dplyr)


## Generate Source data figure 2



library(tidyverse)
all.kat.files <- list.files("tables/phospho_path_clean//", full.names = T, recursive = F, pattern = "\\.csv|\\.txt")
all.kat.files <- all.kat.files[1]
all.P38.files <- read.csv(all.kat.files)

all.P38.files0 <- grep("P38",all.P38.files[all.P38.files$Direction..AD.vs.Control. != "NONE",]$MSigDB.name,value = T)
all.P38.files0 <- paste0("Pathway.",all.P38.files0)
# The following codes are not provided

# Directory to the pathway gene membership
pathways2 <- readRDS('preprocessed/MSigDBPathGeneTab.RDS')


p38genes <- pathways2[pathways2$Pathway %in%all.P38.files0,]
p38genes$ENTREZID <- NULL
p38genes <- as.data.frame(unique(p38genes))
p38genes$Pathway <- as.character(p38genes$Pathway)

phosphoGenes <- read.csv("tables/DE_phospho//phospho_differential_Abeta42_control.csv")
phosphoGenes <- phosphoGenes[phosphoGenes$ENSEMBL %in% p38genes$ENSEMBL,]
rownames(phosphoGenes) <- phosphoGenes$ENSEMBL

p38genes <- p38genes[p38genes$ENSEMBL %in%rownames(phosphoGenes), ]



n_paths <- length(unique((p38genes$Pathway)))
n_genes <- length(unique((p38genes$ENSEMBL)))


heat_mat <- matrix(0,n_paths,n_genes)
rownames(heat_mat) <- unique((p38genes$Pathway))
colnames(heat_mat) <- unique((p38genes$ENSEMBL))

for(i in 1:nrow(p38genes)){
    heat_mat[p38genes[i,1],p38genes[i,2]] <- 1
}

scaled_mat <- heat_mat
rownames(all.P38.files) <-  paste0("Pathway.",all.P38.files$MSigDB.name)
all.P38.files <- all.P38.files[rownames(scaled_mat),]
all.P38.files <- all.P38.files[order(all.P38.files$P.Value),]
scaled_mat    <- scaled_mat[rownames(all.P38.files),]
rownames(scaled_mat) <- all.P38.files$Pathway.Description

# mat <- as.matrix(enrich_de0)
col_fun = circlize::colorRamp2(c(0,10), c( "#f7f7f7",
                                                  "darkviolet"))


col_fun_2 <- circlize::colorRamp2(c(-15,0,15),
                                  c("dodgerblue4","grey90","#67001f"))

phosphoGenes$signed_p <- -log(phosphoGenes$P.Value) * sign(phosphoGenes$logFC)
phosphoGenes$display_name   <- paste0(phosphoGenes$Gene.Symbol," (",
                                      phosphoGenes$Accession,")")

scaled_mat <- scaled_mat[,rownames(phosphoGenes)]
colnames(scaled_mat) <- phosphoGenes$display_name
DE0 <- phosphoGenes
DE0 <- as.matrix(DE0[,c("signed_p","logFC")])
rownames(DE0) <- (phosphoGenes$display_name)
DE0 <- as.data.frame(DE0)
DE0$pch <- NA
DE0[phosphoGenes[phosphoGenes$Direction !="NONE",]$display_name,]$pch <- "*"

set.seed(1)
ht1 <- Heatmap(scaled_mat,
           width =  unit(6.5,"inch"), 
           height = unit(1.2,"inch"),
           row_title = "",
           column_title = "Pathway membership of genes/proteins",
           clustering_distance_columns =  "binary",
           clustering_method_columns = "ward.D2",
           col = circlize::colorRamp2(c(0,1),c("#e0e0e0","grey20") ),
           show_row_names = T, 
           row_names_side = "left",
           show_column_names = T,
           cluster_rows = F,
           show_row_dend = F,
           top_annotation = HeatmapAnnotation("Dysregulation score" = 
                                     anno_simple(DE0$signed_p,
                                                 col = (col_fun_2),
                                                 pch = DE0$pch,
                                                 pt_gp = gpar(size = 0))),
           show_column_dend = F,
           show_heatmap_legend = F,
           border_gp = gpar(col = "black", lty = 2),
           rect_gp = gpar(col = "white", lwd = 2),
           column_title_gp = gpar(fill = "white",col = "gray35",
                                  border = "white", cex = 1.5),
           row_names_gp = gpar(cex =1.2),
           column_names_gp = gpar(cex =1.2),
           column_title_side = "top",
           column_names_rot =  90,
           heatmap_legend_param = list(direction = "horizontal"),na_col = "grey60")

lgd_pvalue = Legend(title = "Dysregulation score", col_fun = col_fun_2)
# and one for the significant p-values
lgd_sig = Legend(pch = "*", type = "points", labels = "Significant Dysregulation")

draw(ht1, annotation_legend_list = list(lgd_pvalue, lgd_sig))


pdf("figures/phosphoprotein_pathway_membership.pdf", width = 22, height = 10)
draw(ht1, annotation_legend_list = list(lgd_pvalue, lgd_sig))
dev.off()
write.csv(all.P38.files,"tables/p38pathways_phospho.csv")
