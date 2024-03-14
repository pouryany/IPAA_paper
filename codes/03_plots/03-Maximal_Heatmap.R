rm(list = ls())


AD.covariates <- readRDS("preprocessed/MSBB/COVARIATES.RDS")
AD.covariates$Tissue.Diagnosis <- paste0(AD.covariates$Tissue,".",AD.covariates$Diagnosis)
AD.covariates$Tissue.Diagnosis <- factor(AD.covariates$Tissue.Diagnosis)

de.paths <-  readRDS("Output/MSBB/MSBB_PHG_ADvsControl _DEP_AD_vs_Control.RDS")


phg_meta <- AD.covariates[grep("PHG.AD|PHG.CONTROL",AD.covariates$Tissue.Diagnosis),]
phg_residuals <- de.paths$PathwayResiduals[,rownames(phg_meta)]

library(gridExtra)
library(ggfortify)




sel.paths <- read.csv("tables/shared_pathways/45vs47_Mayo_MSBB_ROSMAP_2.csv")
sel.paths <- sel.paths[order(sel.paths$P.Value.x),]
sel.paths <- sel.paths[order(sel.paths$Direction.x,decreasing = T),]
var.paths <-  sel.paths$X.x




phg_meta <- AD.covariates[grep("PHG.AD|PHG.CONTROL",AD.covariates$Tissue.Diagnosis),]
phg_meta <- phg_meta[order(phg_meta$Diagnosis),]
phg_residuals <- de.paths$PathwayResiduals[,rownames(phg_meta)]

path.sums <- phg_residuals
rownames(path.sums) <- gsub("Pathway.","",rownames(path.sums))

paths.mat <- path.sums[var.paths,]


labs <- as.data.frame(phg_meta[,c("Tissue.Diagnosis")])
colnames(labs) <- "Diagnosis"
rownames(labs) <- rownames(phg_meta)

labs$Diagnosis <- as.character(labs$Diagnosis)
labs$Diagnosis <- factor(labs$Diagnosis)

# labs <- labs[order(labs$Diagnosis),]


col.pal <- RColorBrewer::brewer.pal(9, "Set3")

color.lab <- function(X){
    col.len <- length(X)
    col.pal <-  c("#F8766D" , "#00BFC4" )
    names(col.pal) <- X
    return(col.pal)
}


col.pal2 <- color.lab(unique(phg_meta$Tissue.Diagnosis))
col.pal2 <- 
annot_col <- list(
    Diagnosis = col.pal2)

# 
# save_pheatmap_pdf(res, "temp_report2/MSBB_heatmap.pdf",width = 3,height = 6)
# 
# save_pheatmap_pdf(res, "temp_report/Shared_all_MSBB_heatmap.pdf",width = 3,height = 6)
#save_pheatmap_pdf(res, "temp_report/MSBB_heatmap_all.pdf",width = 3,height = 6)


library(ComplexHeatmap)

scaled_mat = t(scale(t(paths.mat)))



ht2 <- Heatmap(scaled_mat, column_split = labs$Diagnosis,
               width =  unit(1.7,"inch"),
               height = unit(5,"inch"),
               row_title = "", column_title = "MSBB", 
               clustering_distance_columns =  "euclidean", clustering_method_columns = "ward.D2",
               col = circlize::colorRamp2(c(-2,0,2),c("#0571b0","#f7f7f7","#ca0020") ),
               top_annotation = HeatmapAnnotation(foo = labs$Diagnosis,
                                                  col = list(foo = annot_col$Diagnosis),
                                                  show_legend = F, show_annotation_name = F),
               show_row_names = F,
               show_column_names = F, cluster_rows = F,
               show_row_dend = F,
               show_column_dend = F,
               show_heatmap_legend = F,
               column_title_gp = gpar(fill = "white", col = "gray25", border = "white", cex = 1.5),
               column_title_side = "bottom")





ht_msbb <- ht2



### ROSMAP PLOTTING

ROS.covariates <- readRDS("preprocessed/ROSMAP/COVARIATES.RDS")
PathResiduals  <- readRDS("Output/ROSMAP/ROSMAP_DLPFC_ADvsControl _DEP_AD_vs_Control.RDS")

dlpfc_meta <- ROS.covariates[grep("AD|Control",ROS.covariates$Diagnosis),]
rownames(dlpfc_meta) 
dlpfc_residuals <- PathResiduals$PathwayResiduals[,rownames(dlpfc_meta)]





path.sums <- dlpfc_residuals
rownames(path.sums) <- gsub("Pathway.","",rownames(path.sums))

paths.mat <- path.sums[var.paths,]


labs <- as.data.frame(dlpfc_meta[,c("Diagnosis")])
colnames(labs) <- c("Diagnosis")
rownames(labs) <- rownames(dlpfc_meta)

labs$Diagnosis <- as.character(labs$Diagnosis)
labs$Diagnosis <- factor(labs$Diagnosis)

col.pal <- RColorBrewer::brewer.pal(12, "Set3")

# color.lab <- function(X){
#     col.len <- length(X)
#     col.pal <- col.pal[1:col.len]
#     names(col.pal) <- X
#     return(col.pal)
# }


color.lab <- function(X){
    col.len <- length(X)
    col.pal <-  c("#F8766D" , "#00BFC4" )
    names(col.pal) <- X
    return(col.pal)
}



col.pal2 <- color.lab(unique(dlpfc_meta$Diagnosis))
#col.pal1 <- color.lab(unique(dlpfc_meta$braaksc))

annot_col <- list(
    Diagnosis = col.pal2)


ord <- order(labs$Diagnosis)


scaled_mat = t(scale(t(paths.mat)))




ht2 <- Heatmap(scaled_mat, column_split = labs$Diagnosis,
               cluster_column_slices = F,
               width =  unit(1.7,"inch"),
               height = unit(5,"inch"),
               row_title = "", column_title = "ROSMAP", 
               clustering_distance_columns =  "euclidean", clustering_method_columns = "ward.D2",
               col = circlize::colorRamp2(c(-2,0,2),c("#0571b0","#f7f7f7","#ca0020") ),
               top_annotation = HeatmapAnnotation(foo = labs$Diagnosis,
                                                  col = list(foo = annot_col$Diagnosis),
                                                  show_legend = F, show_annotation_name = F),
               show_row_names = F,
               show_column_names = F, cluster_rows = F,
               show_row_dend = F,
               show_column_dend = F,
               show_heatmap_legend = F,
               column_title_gp = gpar(fill = "white", col = "gray25", border = "white", cex = 1.5),
               column_title_side = "bottom") 





ht_rosmap <- ht2




tcx_meta <- readRDS("preprocessed/MAYO/COVARIATES.RDS")
tcx_meta <- tcx_meta[grep("TCX.AD|TCX.CONTROL",tcx_meta$Tissue.Diagnosis),]

tcx_residuals <- readRDS("Output/MAYO/Mayo_ADvsControl _DEP_AD_vs_CONTROL.RDS")
tcx_residuals <- tcx_residuals$PathwayResiduals
tcx_residuals <- tcx_residuals[,rownames(tcx_meta)]


path.sums <- tcx_residuals
rownames(path.sums) <- gsub("Pathway.","",rownames(path.sums))

paths.mat <- path.sums[var.paths,]


labs <- as.data.frame(tcx_meta[,c("Tissue.Diagnosis")])
colnames(labs) <- "Diagnosis"
rownames(labs) <- rownames(tcx_meta)


labs$Diagnosis <- as.character(labs$Diagnosis)
labs$Diagnosis <- factor(labs$Diagnosis, levels = c( "TCX.CONTROL","TCX.AD"))




color.lab <- function(X){
    col.len <- length(X)
    col.pal <-  c("#F8766D" , "#00BFC4" )
    names(col.pal) <- X
    return(col.pal)
}


col.pal2 <- color.lab(unique(tcx_meta$Tissue.Diagnosis))
col.pal2 <- 
    annot_col <- list(
        Diagnosis = col.pal2)


scaled_mat = t(scale(t(paths.mat)))
#just to get the order correct
#scaled_mat <- scaled_mat[,rev(colnames(scaled_mat))]


labs$Diagnosis <- factor(labs$Diagnosis, levels = c("TCX.AD","TCX.CONTROL"))

ht2 <- Heatmap(scaled_mat, column_split = labs$Diagnosis,
               cluster_column_slices = F,
               width =  unit(1.7,"inch"),
               height = unit(5,"inch"),
               row_title = "", column_title = "Mayo", 
               clustering_distance_columns =  "euclidean", clustering_method_columns = "ward.D2",
               col = circlize::colorRamp2(c(-2,0,2),c("#0571b0","#f7f7f7","#ca0020") ),
               top_annotation = HeatmapAnnotation(foo = labs$Diagnosis,
                                                  col = list(foo = annot_col$Diagnosis),
                                                  show_legend = F, show_annotation_name = F),
               show_row_names = F,
               show_column_names = F, cluster_rows = F,
               show_row_dend = F,
               show_column_dend = F,
               show_heatmap_legend = F,
               column_title_gp = gpar(fill = "white", col = "gray25", border = "white", cex = 1.5),
               column_title_side = "bottom") 





scaled_mat2 <- scaled_mat



display_names    <- read.csv("preprocessed/MSigDB_display_names.csv")
rownames(display_names) <- display_names$STANDARD_NAME



rownames(scaled_mat2) <- display_names[rownames(scaled_mat2),]$display_name

# rownames(scaled_mat2) <- gsub("^([A-Z]*)_","\\1: ",rownames(scaled_mat2))
# rownames(scaled_mat2) <- gsub("_"," ",rownames(scaled_mat2))

library(stringr)
if(any(str_length(rownames(scaled_mat2)) > 70)){
    longInds <- str_length(rownames(scaled_mat2)) > 70
    rownames(scaled_mat2)[longInds] <- stringr::str_sub(rownames(scaled_mat2)[longInds],1,70)
    rownames(scaled_mat2)[longInds] <- paste0(rownames(scaled_mat2)[longInds],"*") 
}




ht4_mayo <- Heatmap(scaled_mat2, column_split = labs$Diagnosis,
               cluster_column_slices = F,
               width =  unit(1.7,"inch"),
               height = unit(5,"inch"),
               row_names_side = "left",
               row_title = "", column_title = "Mayo", 
               clustering_distance_columns =  "euclidean", clustering_method_columns = "ward.D2",
               col = circlize::colorRamp2(c(-2,0,2),c("#0571b0","#f7f7f7","#ca0020") ),
               top_annotation = HeatmapAnnotation(foo = labs$Diagnosis,
                                                  col = list(foo = annot_col$Diagnosis),
                                                  show_legend = T, show_annotation_name = F),
               show_row_names = T,
               show_column_names = F, cluster_rows = F,
               show_row_dend = F,
               show_column_dend = F,
               show_heatmap_legend = T,
               column_title_gp = gpar(fill = "white", col = "gray25", border = "white", cex = 1.5),
               column_title_side = "bottom",
               heatmap_legend_param = list(direction = "horizontal", title = "Pathway activity"))








ht_mayo <- ht2

# 45 47 stuff


library(tibble)
gene.address <- list.files("preprocessed/cell_models/",
                           pattern = "I45F_vs_I47F",
                           recursive = T,
                           full.names = T)

genes.counts2 <- read.delim(gene.address)
genes.counts2 <- genes.counts2 %>% dplyr::select(ID, contains("I45F"), contains("I47F")) %>% column_to_rownames("ID")


AD.covariates2 <- data.frame(Sample = colnames(genes.counts2),
                             Condition = sub("_[0-9].*", "",
                                             colnames(genes.counts2)))
rownames(AD.covariates2) <- AD.covariates2$Sample


cell_meta <- AD.covariates2

cell_residuals <- readRDS("Output/I45FvsI47F/I45FvsI47F _DEP_AD_vs_Control.RDS")
cell_residuals <- cell_residuals$pathwaySummaryStats



path.sums <- cell_residuals
rownames(path.sums) <- gsub("Pathway.","",rownames(path.sums))

paths.mat <- path.sums[var.paths,]


labs <- as.data.frame(cell_meta[,c("Condition")])
colnames(labs) <- "Diagnosis"
rownames(labs) <- rownames(cell_meta)


labs$Diagnosis <- as.character(labs$Diagnosis)
labs$Diagnosis <- factor(labs$Diagnosis)




color.lab <- function(X){
    col.len <- length(X)
    col.pal <-  c("#F8766D" , "#00BFC4" )
    names(col.pal) <- X
    return(col.pal)
}


col.pal2 <- color.lab(unique(cell_meta$Condition))
col.pal2 <- 
    annot_col <- list(
        Diagnosis = col.pal2)


scaled_mat = t(scale(t(paths.mat)))



display_names    <- read.csv("preprocessed/MSigDB_display_names.csv")
rownames(display_names) <- display_names$STANDARD_NAME


scaled_mat2 <- scaled_mat
rownames(scaled_mat2) <- display_names[rownames(scaled_mat2),]$display_name

# rownames(scaled_mat2) <- gsub("^([A-Z]*)_","\\1: ",rownames(scaled_mat2))
# rownames(scaled_mat2) <- gsub("_"," ",rownames(scaled_mat2))

library(stringr)
if(any(str_length(rownames(scaled_mat2)) > 70)){
    longInds <- str_length(rownames(scaled_mat2)) > 70
    rownames(scaled_mat2)[longInds] <- stringr::str_sub(rownames(scaled_mat2)[longInds],1,70)
    rownames(scaled_mat2)[longInds] <- paste0(rownames(scaled_mat2)[longInds],"*") 
}


ht2 <- Heatmap(scaled_mat2, column_split = labs$Diagnosis,
               width =  unit(1.7,"inch"),
               height = unit(5,"inch"),
               row_title = "", column_title = "A\u03B242-H vs A\u03B240-H", 
               clustering_distance_columns =  "euclidean", clustering_method_columns = "ward.D2",
               col = circlize::colorRamp2(c(-2,0,2),c("#0571b0","#f7f7f7","#ca0020") ),
               top_annotation = HeatmapAnnotation(foo = labs$Diagnosis,
                                                  col = list(foo = annot_col$Diagnosis),
                                                  show_legend = F, show_annotation_name = F),
               show_row_names = T,
               row_names_side = "left",
               show_column_names = F, cluster_rows = F,
               show_row_dend = F,
               show_column_dend = F,
               show_heatmap_legend = F,
               column_title_gp = gpar(fill = "white", col = "gray25", border = "white", cex = 1.5),
               column_title_side = "bottom") 




ht3 <- Heatmap(scaled_mat, column_split = labs$Diagnosis,
               width =  unit(1.7,"inch"),
               height = unit(5,"inch"),              
               row_names_side = "left",
               row_title = "", column_title = "A\u03B242-H vs A\u03B240-H", 
               clustering_distance_columns =  "euclidean", clustering_method_columns = "ward.D2",
               col = circlize::colorRamp2(c(-2,0,2),c("#0571b0","#f7f7f7","#ca0020") ),
               top_annotation = HeatmapAnnotation(foo = labs$Diagnosis,
                                                  col = list(foo = annot_col$Diagnosis),
                                                  show_legend = T, show_annotation_name = F),
               show_row_names = T,
               show_column_names = F, cluster_rows = F,
               show_row_dend = F,
               show_column_dend = F,
               show_heatmap_legend = T,
               column_title_gp = gpar(fill = "white", col = "gray25", border = "white", cex = 1.5),
               column_title_side = "bottom",
               heatmap_legend_param = list(direction = "horizontal")) 








ht_i45 <- ht2



ht_list <- ht_i45 + ht_mayo + ht_msbb + ht_rosmap


cairo_pdf(file = paste0("figures/revised_05_a_all.pdf"),
          width = 18, height = 10, onefile = T)
draw(ht_list, ht_gap =unit(0.5,"inch") )
ht3
ht4_mayo


dev.off()



