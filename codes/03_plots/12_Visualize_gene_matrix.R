rm(list = ls())
library(tidyverse)
library(ComplexHeatmap)
library(dplyr)

file_names <- list.files("Output/gene_similarity/", pattern = "top|spearman", full.names = T)
file_names <- grep("spearman",file_names,value = T)

for (j in file_names){
    
res_df <- read.csv(j)

# remove the below two lines if needed


analysis_name <- gsub(".*similarity//|2021-03-15.*", "", j)  

    res_df_sq <- res_df 
    res_df_sq[, 1] <- gsub(pattern = "2021.*|2022.*|2023.*|.*/|_QC.*|diffPathways_DE_|_diffPathways.*|\\.csv|_CPM_log2FC_.*|\\.txt|deseq2_gene_analysis|_genes|_gene|GSE_dataset|GSEdataset","", res_df_sq[,1])
    res_df_sq[, 1] <- gsub(pattern = "Normal","Control", res_df_sq[,1])
    res_df_sq <- distinct(res_df_sq)
    rownames(res_df_sq) <- res_df_sq[,1]
    res_df_sq <- res_df_sq[,-1] %>% distinct()
    res_df_sq <- t(res_df_sq) %>% as.data.frame() %>% distinct()
    # colnames(res_df_sq) <- gsub(pattern = "\\.txt|\\.csv","", colnames(res_df_sq))
    # colnames(res_df_sq) <- gsub(pattern = ".*\\.|_QC.*|_diffPathways.*|\\.csv|_CPM_log2FC_.*|\\.txt|deseq2_gene_analysis|_genes|_gene|GSE_dataset","", colnames(res_df_sq))
    # rownames(res_df_sq) <- gsub(pattern = ".*\\/|_QC.*|_diffPathways.*|\\.csv|_CPM_log2FC_.*|\\.txt|deseq2_gene_analysis|_genes|_gene|GSE_dataset","", rownames(res_df_sq))
    rownames(res_df_sq) <- colnames(res_df_sq)
    
    colnames(res_df_sq) <- gsub("_vs_G2B2"," vs Ctrl",colnames(res_df_sq))
    colnames(res_df_sq) <- gsub("I45FvsI47F","A\u03B242-H vs A\u03B240-H",colnames(res_df_sq))
    colnames(res_df_sq) <- gsub("I45F","A\u03B242-H",colnames(res_df_sq))
    colnames(res_df_sq) <- gsub("I47F","A\u03B240-H",colnames(res_df_sq))
    #colnames(res_df_sq) <- gsub("v.*","",colnames(res_df_sq))
    colnames(res_df_sq) <- gsub("_"," ",colnames(res_df_sq))
    # colnames(res_df_sq) <- gsub("ROSMAP","DLPFC",colnames(res_df_sq))
    # colnames(res_df_sq) <- gsub("Mayo","TCX",colnames(res_df_sq))
    # colnames(res_df_sq) <- gsub("MSBB","PHG",colnames(res_df_sq))
    # 
    colnames(res_df_sq) <- gsub("vs Control","vs Ctrl",colnames(res_df_sq))
    
    #colnames(res_df_sq) <- paste0(colnames(res_df_sq), "vs Ctrl")
    #colnames(res_df_sq)[5] <- gsub("v.*","vs A\u03B240-H",colnames(res_df_sq)[5])
    
    rownames(res_df_sq) <- colnames(res_df_sq)
    
    
  
 
    DE0 <- rownames(res_df_sq)
    
    DE0 <- data.frame(assay = DE0, source = "Cellular model")
    DE0[grep("ROSMAP|Mayo|MSBB",DE0$assay),]$source <- "Human brain"
    
    
    
    row_ha = rowAnnotation("Source" = DE0$source, 
                           col = list(Source = c("Cellular model" = "peachpuff",
                                                 "Human brain" = "darkgoldenrod")))
    
    
    
    
    # set.seed(i)
    ht5 <- ComplexHeatmap::Heatmap(res_df_sq, 
                                   col = circlize::colorRamp2(c(-1,-0.5,0,0.5, 1), 
                                                              colors = c("darkblue", 
                                                                         "cornflowerblue", 
                                                                         # "grey95",
                                                                         "grey95",
                                                                         # "grey95",
                                                                         "coral",
                                                                         "tomato2")),
                                   # name = "-log10(pval)",
                                   right_annotation = row_ha,
                                   clustering_method_columns = "ward.D2",
                                   clustering_method_rows = "ward.D2",
                                   name = "correlation",
                                   width = unit(1, "cm")*ncol(res_df_sq), 
                                   height = unit(1, "cm")*nrow(res_df_sq),
                                   row_split = 2,
                                   column_split = 2,
                                   row_title = c("",  ""),
                                   column_title = c("",  ""),
                                   row_gap = unit(c(3, 3), "mm"),
                                   column_gap = unit(c(3, 3), "mm"),
                                   show_column_dend = F, 
                                   show_row_dend = F,
                                   cell_fun = function(j, i, x, y, width, height, fill) {
                                       grid.text(sprintf("%.2f", res_df_sq[i, j]), x, y, gp = gpar(fontsize = 10))}
    )
    draw(ht5)
    
    
    
    
    
    cairo_pdf(file = paste0("figures/03_e_Gene-",
                            analysis_name, ".pdf"),
              width = 10,
              height = 10,
              onefile = T)

    draw(ht5)
    dev.off()


}

