rm(list = ls())
library(tidyverse)
library(ComplexHeatmap)
library(dplyr)
library(stringr)
file_names <- list.files("Output/similarity_external/", 
                         pattern = "top|spearman",
                         full.names = T)

file_names <- grep("2023",file_names,value = T)
for (j in file_names){
    
res_df <- read.csv(j)
analysis_name <- gsub(".*similarity//|2021-02-08.*|2021-[0-9]*-[0-p].*", "", j)  
analysis_name <- tail(unlist(str_split(analysis_name,"/")),1)

    res_df_sq <- res_df 
    res_df_sq[, 1] <- gsub(pattern = "2021-03-15|2021-04-28|.*/|_QC.*|diffPathways_DE_|_diffPathways.*|\\.csv|_CPM_log2FC_.*|\\.txt|deseq2_gene_analysis|_genes|_gene|GSE_dataset|GSEdataset","", res_df_sq[,1])
    res_df_sq[, 1] <- gsub(pattern = "Normal","Control", res_df_sq[,1])
    res_df_sq <- distinct(res_df_sq)
    rownames(res_df_sq) <- res_df_sq[,1]
    res_df_sq <- res_df_sq[,-1] %>% distinct()
    res_df_sq <- t(res_df_sq) %>% as.data.frame() %>% distinct()
    # colnames(res_df_sq) <- gsub(pattern = "\\.txt|\\.csv","", colnames(res_df_sq))
    # colnames(res_df_sq) <- gsub(pattern = ".*\\.|_QC.*|_diffPathways.*|\\.csv|_CPM_log2FC_.*|\\.txt|deseq2_gene_analysis|_genes|_gene|GSE_dataset","", colnames(res_df_sq))
    # rownames(res_df_sq) <- gsub(pattern = ".*\\/|_QC.*|_diffPathways.*|\\.csv|_CPM_log2FC_.*|\\.txt|deseq2_gene_analysis|_genes|_gene|GSE_dataset","", rownames(res_df_sq))
    
    colnames(res_df_sq) <- gsub("_.*","",colnames(res_df_sq))
    colnames(res_df_sq) <- gsub("vswt"," vs WT",colnames(res_df_sq))
    colnames(res_df_sq) <- gsub("ROSMAP","ROSMAP AD vs Ctrl",colnames(res_df_sq))
    colnames(res_df_sq) <- gsub("Mayo","Mayo  AD vs Ctrl",colnames(res_df_sq))
    colnames(res_df_sq) <- gsub("MSBB","MSBB  AD vs Ctrl",colnames(res_df_sq))
    colnames(res_df_sq) <- gsub("APP","APP vs WT",colnames(res_df_sq))
    colnames(res_df_sq) <- gsub("PSEN","PSEN1 vs WT",colnames(res_df_sq))
    rownames(res_df_sq) <- colnames(res_df_sq) 

    # set.seed(i)
    set.seed(10)
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
                                   name = "correlation",
                                   width = unit(1, "cm")*ncol(res_df_sq), 
                                   height = unit(1, "cm")*nrow(res_df_sq),
                                   row_split = 2,
                                   column_split = 2,
                                   row_title = c("", ""),
                                   column_title = c("", ""),
                                   show_column_dend = F, 
                                   show_row_dend = F,
                                   cell_fun = function(j, i, x, y, width, height, fill) {
                                       grid.text(sprintf("%.2f", res_df_sq[i, j]), x, y, gp = gpar(fontsize = 8))}
    )
    draw(ht5)
   
    cairo_pdf(file = paste0("figures/S03_heatmap_",
                            analysis_name, ".pdf"),
              width = 10, height = 10,
              onefile = T)

    draw(ht5)
    dev.off()


}

