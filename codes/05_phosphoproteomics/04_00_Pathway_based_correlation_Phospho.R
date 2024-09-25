rm(list = ls())



library(tidyverse)
all.kat.files0 <- list.files("output", full.names = T,recursive = T,
                             pattern = "diff")
all.kat.files0 <- grep("I45",all.kat.files0,invert = F,value = T)

all.kat.files <- list.files("output_phospho/", full.names = T,recursive = T,
                            pattern = "diff")
all.kat.files <- grep("Losma",all.kat.files,invert = T,value = T)

all.kat.files <- c(all.kat.files,all.kat.files0)


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

res_df_sq <- res_df %>% dplyr::select(-adj.P.Val) %>% pivot_wider(names_from = file, values_from = t)
cor_sq <- cor(res_df_sq[,-1], method = "spearman", 
              use = "pairwise.complete.obs")
cor_sg <- psych::corr.test(res_df_sq[,-1], method = "spearman", 
              use = "pairwise.complete.obs", adjust = "fdr")

print(cor_sg, short= FALSE)
dim(cor_sg$r)

cor_sg$stars

confidence_mat <- cor_sg$r
confidence_mat[!upper.tri(confidence_mat)] <- NA
confidence_mat[upper.tri(confidence_mat)] <- (cor_sg$p.adj)

replace(confidence_mat,upper.tri(confidence_mat),1:21 )

upper.tri(res_df_sq[,-1])

temp <- print(cor_sg, short= FALSE, digits = 20)

p.adjust(temp$raw.p, method = "fdr")


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
    
    
    colnames(res_df_sq) <- gsub("vsG2B2"," vs Ctrl",colnames(res_df_sq))
    colnames(res_df_sq) <- gsub("I45FvsI47F","A\u03B242-H vs A\u03B240-H",colnames(res_df_sq))
    colnames(res_df_sq) <- gsub("I45F","A\u03B242-H",colnames(res_df_sq))
    colnames(res_df_sq) <- gsub("I45","A\u03B242-H",colnames(res_df_sq))
     colnames(res_df_sq) <- gsub("I47F","A\u03B240-H",colnames(res_df_sq))
    #colnames(res_df_sq) <- gsub("v.*","",colnames(res_df_sq))
    colnames(res_df_sq) <- gsub("_"," ",colnames(res_df_sq))
    colnames(res_df_sq) <- gsub("Phospho","P-",colnames(res_df_sq))
    colnames(res_df_sq) <- gsub("Control","Ctrl",colnames(res_df_sq))
    

    colnames(res_df_sq) <- gsub("A\u03B242-H vs A\u03B240-H vs Ctrl","A\u03B242-H vs A\u03B240-H",colnames(res_df_sq))
    
    rownames(res_df_sq) <- colnames(res_df_sq)
    

    
    
    col_fun_2 <- circlize::colorRamp2(c(0,7),c("#f7f7f7", "#67001f") )
    col_fun2  <- c("coral","tomato")
    
    
    library(ComplexHeatmap)
    ha = HeatmapAnnotation(bar = sample(letters[1:3], 10, replace = TRUE),
                           col = list(bar = c("a" = "red", "b" = "green", "c" = "blue")))
    
    DE0 <- rownames(res_df_sq)
    
    DE0 <- data.frame(assay = DE0, source = "Transcriptome")
    DE0[grep("P-",DE0$assay),]$source <- "Phosphoproteome"
    
    
    
    row_ha = rowAnnotation("Source" = DE0$source, 
                           col = list(Source = c("Transcriptome" = "darkgreen",
                                                 "Phosphoproteome" = "darkmagenta")))
    
    
    

    set.seed(10)
    ht5 <- ComplexHeatmap::Heatmap(res_df_sq, 
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
                                   name = "correlation",
                                   width = unit(1.4, "cm")*ncol(res_df_sq), 
                                   height = unit(1.4, "cm")*nrow(res_df_sq),
                                   row_split = 2,
                                   column_split = 2,
                                   row_title = c("",  ""),
                                   column_title = c("",  ""),
                                   show_column_dend = F, 
                                   show_row_dend = F,
                                   cell_fun = function(j, i, x, y, width, height, fill) {
                                       grid.text(sprintf("%.2f", res_df_sq[i, j]), x, y, gp = gpar(fontsize = 10))}
    )
    
    
    set.seed(10)
    
    
    res_df_sq2 <- cor_sg$stars
    rownames(res_df_sq2) <- colnames(res_df_sq2) <- rownames(res_df_sq)
    
    ht6 <- ComplexHeatmap::Heatmap(res_df_sq, 
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
                                   name = "correlation",
                                   clustering_method_columns = "complete",
                                   width = unit(1.4, "cm")*ncol(res_df_sq), 
                                   height = unit(1.4, "cm")*nrow(res_df_sq),
                                   row_split = 2,
                                   column_split = 2,
                                   row_title = c("", ""),
                                   column_title = c("", ""),
                                   show_column_dend = F, 
                                   show_row_dend = F,
                                   cell_fun = function(j, i, x, y, width, height, fill) {
                                       grid.text(sprintf(res_df_sq2[i, j]), x, y, gp = gpar(fontsize = 8))}
    )
    
    
    draw(ht6)
    
    # p and adjust p are the same and already applied. For raw p-values get raw.pvalue from summary
    res_df_sq2 <- cor_sg$p

    rownames(res_df_sq2) <- colnames(res_df_sq2) <- rownames(res_df_sq)
    res_df_sq2 <- format(res_df_sq2,digits = 3)
    ht7 <- ComplexHeatmap::Heatmap(res_df_sq, 
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
                                   name = "correlation",
                                   clustering_method_columns = "complete",
                                   width = unit(1.4, "cm")*ncol(res_df_sq), 
                                   height = unit(1.4, "cm")*nrow(res_df_sq),
                                   row_split = 2,
                                   column_split = 2,
                                   row_title = c("", ""),
                                   column_title = c("", ""),
                                   show_column_dend = F, 
                                   show_row_dend = F,
                                   cell_fun = function(j, i, x, y, width, height, fill) {
                                       grid.text(sprintf((res_df_sq2[i, j])), x, y, gp = gpar(fontsize = 8))}
    )
    
    
    draw(ht7)
    
    
    analysis_name <- "pathway_correlation_"
    cairo_pdf(file = paste0("figures/Phospho_", analysis_name, ".pdf"),
              width = 10, height = 10, onefile = T)
    draw(ht5)
    draw(ht6)
    draw(ht7)
    dev.off()
    
    





