rm(list = ls())



library(tidyverse)
all.kat.files <- list.files("reduced_pathways_output/", full.names = T,recursive = T,
                            pattern = "diff")

all.kat.files <- setdiff(all.kat.files,grep("Tg|5xFAD|hAb",all.kat.files,value = T))
base_file <- tail(all.kat.files,1)


# limit Rosmap to the common pathways
keep.paths <- read_csv(head(all.kat.files,1))
colnames(keep.paths)[1] <- "pathway"




# rank correlation sorted by t statistic 
res_df <- data.frame(pathway = NA, t = NA, adj.P.Val = NA, file = NA)
res_df <- res_df[-1,]
for (j in all.kat.files){
    pour.paths <- read_csv(j)
    colnames(pour.paths)[1] <- "pathway"
    pour.paths <- pour.paths %>% filter(pathway %in% keep.paths$pathway)
    pour.ups      <- pour.paths[,c("pathway","t","adj.P.Val")]
    pour.ups$file <- j
    head(pour.ups)
    res_df <- rbind(res_df, pour.ups)
}

res_df_sq <- res_df %>% dplyr::select(-adj.P.Val) %>%
    pivot_wider(names_from = file, values_from = t)


cor_sq <- cor(res_df_sq[,-1], method = "spearman", 
              use = "pairwise.complete.obs")
cor_sg <- psych::corr.test(res_df_sq[,-1], method = "spearman", 
              use = "pairwise.complete.obs", adjust = "none")


# Adjusting p-values
# There is a bug for symmetric matrices in psych::corr.test
adjusted_pvals <- cor_sg$p[upper.tri(cor_sg$p)]
adjusted_pvals <- p.adjust(adjusted_pvals, method = "fdr")

temp <- cor_sg$p.adj
temp[upper.tri(temp)] <- adjusted_pvals
temp <- t(temp)
temp[upper.tri(temp)] <- adjusted_pvals
cor_sg$p.adj <- temp

# 
# 
# write.csv(cor_sq, file = paste0("Output/pathway_similarity/spearman_cor_AD_model_sq",
#                                 Sys.Date(), ".csv"))
# 




    res_df_sq <- cor_sq 
    rownames(res_df_sq) <- gsub(pattern = "2021-03-15|2021-02-08|.*/|_QC.*|diffPathways_DE_|_diffPathways.*|\\.csv|_CPM_log2FC_.*|\\.txt|deseq2_gene_analysis|_genes|_gene|GSE_dataset|GSEdataset","", 
                           rownames(res_df_sq))
    rownames(res_df_sq) <- gsub(pattern = "Normal","Control", rownames(res_df_sq))
    colnames(res_df_sq) <- rownames(res_df_sq)
    
    colnames(res_df_sq) <- gsub("vsG2B2","",colnames(res_df_sq))
    colnames(res_df_sq) <- gsub("I45FvsI47F","A\u03B242-H vs A\u03B240-H",colnames(res_df_sq))
    colnames(res_df_sq) <- gsub("I45F","A\u03B242-H",colnames(res_df_sq))
    colnames(res_df_sq) <- gsub("I47F","A\u03B240-H",colnames(res_df_sq))
    #colnames(res_df_sq) <- gsub("v.*","",colnames(res_df_sq))
    colnames(res_df_sq) <- gsub("_"," ",colnames(res_df_sq))
    colnames(res_df_sq) <- gsub("PHG ","",colnames(res_df_sq))
    colnames(res_df_sq) <- gsub("DLPFC ","",colnames(res_df_sq))
    colnames(res_df_sq) <- gsub(" ADvsControl"," AD",colnames(res_df_sq))
    

    colnames(res_df_sq) <- paste0(colnames(res_df_sq), " vs Ctrl")
    colnames(res_df_sq) <- gsub("A\u03B242-H vs A\u03B240-H vs Ctrl","A\u03B242-H vs A\u03B240-H",colnames(res_df_sq))
    
    rownames(res_df_sq) <- colnames(res_df_sq)
    

    
    
    col_fun_2 <- circlize::colorRamp2(c(0,7),c("#f7f7f7", "#67001f") )
    col_fun2  <- c("coral","tomato")
    
    

    library(ComplexHeatmap)
    ha = HeatmapAnnotation(bar = sample(letters[1:3], 10, replace = TRUE),
                           col = list(bar = c("a" = "red", "b" = "green", "c" = "blue")))
    
    DE0 <- rownames(res_df_sq)
    
    DE0 <- data.frame(assay = DE0, source = "Cellular model")
    DE0[grep("ROSMAP|Mayo|MSBB",DE0$assay),]$source <- "Human brain"
    
    
    
    row_ha = rowAnnotation("Source" = DE0$source, 
                           col = list(Source = c("Cellular model" = "peachpuff",
                                                 "Human brain" = "darkgoldenrod")))
    
    

    set.seed(9)
    ht5 <- ComplexHeatmap::Heatmap(res_df_sq,
                                   cluster_rows = T,
                                   cluster_columns = T,
                                   right_annotation = row_ha,
                                   col = circlize::colorRamp2(c(-1,-0.5,0,0.5, 1), 
                                                              colors = c("darkblue", 
                                                                         "cornflowerblue", 
                                                                         # "grey95",
                                                                         "grey95",
                                                                         # "grey95",
                                                                         "coral",
                                                                         "tomato2")),
                                   # name = "-log10(pval)",
                                   clustering_method_columns = "ward.D2",
                                   clustering_method_rows = "ward.D2",
                                   name = "correlation",
                                   width = unit(1, "cm")*ncol(res_df_sq), 
                                   height = unit(1, "cm")*nrow(res_df_sq),
                                   show_column_dend = F, 
                                   show_row_dend = F,
                                   cell_fun = function(j, i, x, y, width, height, fill) {
                                       grid.text(sprintf("%.2f", res_df_sq[i, j]), x, y, gp = gpar(fontsize = 10))}
    )
    
    
   
    
# The next plot will show p-values in the same exact layout as the last heatmap
    # preserving column order for the next plot
    col_order1 <- c(7,6,8,4,5,2,3,1)
    #col_order1 <- c(2,3,1,4,5,8,6,7)
    res_df_sq_new_ord <- res_df_sq[col_order1,col_order1]
    
# p and adjust p are the same and already applied. 
# For raw p-values get raw.pvalue from summary
    res_df_sq22 <- cor_sg$p.adj[col_order1,col_order1]
    rownames(res_df_sq22) <- colnames(res_df_sq22) <- rownames(res_df_sq22)
    res_df_sq3 <- res_df_sq22
       
    res_df_sq22 <- format(res_df_sq3,digits = 3)
    res_df_sq22[res_df_sq3 > 0.05] <- ""
    res_df_sq22[res_df_sq3 < 2.2e-16] <- "2.2e-16"
    
    diag(res_df_sq22) <- ""
    #res_df_sq22[lower.tri(res_df_sq2)] <- ""


    res_df_sq4 <- res_df_sq_new_ord
    diag(res_df_sq4) <- NA
    #res_df_sq4[lower.tri(res_df_sq4)] <- NA
    
    #res_df_sq4 <- res_df_sq4[rownames(res_df_sq3),colnames(res_df_sq3)]
    
    DE0 <- rownames(res_df_sq4)
    
    DE0 <- data.frame(assay = DE0, source = "Cellular model")
    DE0[grep("ROSMAP|Mayo|MSBB",DE0$assay),]$source <- "Human brain"
    
    
    
    row_ha2 = rowAnnotation("Source" = DE0$source, 
                           col = list(Source = c("Cellular model" = "peachpuff",
                                                 "Human brain" = "darkgoldenrod")))
    
    

    
    ht7 <- ComplexHeatmap::Heatmap(res_df_sq4, cluster_rows = F, cluster_columns = F,
                                   right_annotation = row_ha2,
                                   border = "black", na_col = "grey35",
                                   rect_gp = gpar(col = "white", lwd = 2),
                                   col = circlize::colorRamp2(c(-1,-0.5,0,0.5, 1), 
                                                              colors = c("grey95", 
                                                                         "grey95", 
                                                                         # "grey95",
                                                                         "grey95",
                                                                         # "grey95",
                                                                         "grey95",
                                                                         "grey95")),
                                   # name = "-log10(pval)",
                                   name = "correlation",
                                   width = unit(1, "cm")*ncol(res_df_sq), 
                                   height = unit(1, "cm")*nrow(res_df_sq),
                                   show_column_dend = F, 
                                   show_row_dend = F,
                                   cell_fun = function(j, i, x, y, width,
                                                       height, fill) {
                                       grid.text(sprintf((res_df_sq22[i, j])),
                                                 x, y, gp = gpar(fontsize = 6))}
    )
    
    
    draw(ht7)
    
    
    analysis_name <- "pathway_correlation_"
    cairo_pdf(file = paste0("figures/S03_d_reduced_", analysis_name, ".pdf"),
              width = 10, height = 10, onefile = T)
    draw(ht5)
    draw(ht7)
    dev.off()
    
    





