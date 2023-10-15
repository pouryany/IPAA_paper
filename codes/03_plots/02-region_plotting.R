rm(list = ls())


tcx_meta <- readRDS("preprocessed/MAYO/COVARIATES.RDS")
tcx_meta <- tcx_meta[grep("TCX.AD|TCX.CONTROL",tcx_meta$Tissue.Diagnosis),]
tcx_meta <- tcx_meta[order(tcx_meta$Tissue.Diagnosis),]

tcx_residuals <- readRDS(list.files("Output/MAYO/",pattern = "Residual",full.names = T))
tcx_residuals <- tcx_residuals[,rownames(tcx_meta)]





AD.covariates <- readRDS("preprocessed/MSBB/COVARIATES.RDS")
AD.covariates$Tissue.Diagnosis <- paste0(AD.covariates$Tissue,".",AD.covariates$Diagnosis)
AD.covariates$Tissue.Diagnosis <- factor(AD.covariates$Tissue.Diagnosis)

phg_meta <-  AD.covariates[grep("PHG.AD|PHG.CONTROL",AD.covariates$Tissue.Diagnosis),]
de.paths <-  readRDS("Output/MSBB/MSBB_PHG_ADvsControl _DEP_AD_vs_Control.RDS")
phg_residuals <- de.paths$PathwayResiduals[,rownames(phg_meta)]



library(gridExtra)
library(ggfortify)



sel.paths <- read.csv("tables/shared_pathways/01_shared_pathways_MSBB_Mayo.csv")

sel.paths <- sel.paths[order(sel.paths$P.Value.x),]
sel.paths <- sel.paths[order(sel.paths$Direction.x),]
var.paths <-  sel.paths$X


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
    col.pal <-  c("#F8766D","#00BFC4" )
    names(col.pal) <- X
    return(col.pal)
}


col.pal2 <- color.lab(unique(phg_meta$Tissue.Diagnosis))
col.pal2 <- 
    annot_col <- list(
        Diagnosis = col.pal2)


library(ComplexHeatmap)

scaled_mat = t(scale(t(paths.mat)))




scaled_mat_phg = t(scale(t(paths.mat)))
labs_phg       = labs
annot_col_phg  = annot_col



path.sums <- tcx_residuals
rownames(path.sums) <- gsub("Pathway.","",rownames(path.sums))

paths.mat <- path.sums[var.paths,]
scaled_mat = t(scale(t(paths.mat)))









labs <- as.data.frame(tcx_meta[,c("Tissue.Diagnosis")])
colnames(labs) <- "Diagnosis"
rownames(labs) <- rownames(tcx_meta)


labs$Diagnosis <- as.character(labs$Diagnosis)
labs$Diagnosis <- gsub("TCX.","",labs$Diagnosis)

labs$Diagnosis <- factor(labs$Diagnosis, levels = c("CONTROL","AD"))

labs2  <- labs$Diagnosis 
labs2 <- factor(labs2, levels = c("AD","CONTROL"))


color.lab <- function(X){
    col.len <- length(X)
    col.pal <-  c( "#F8766D", "#00BFC4")
    names(col.pal) <- X
    return(col.pal)
}


col.pal2 <- color.lab(unique(labs$Diagnosis))
col.pal2 <- 
    annot_col <- list(
        Diagnosis = rev(col.pal2))




scaled_mat = t(scale(t(paths.mat)))




#labs$Diagnosis <- factor(labs$Diagnosis, levels = c("TCX.AD", "TCX.CONTROL"))

#scaled_mat <- scaled_mat[,rev(colnames(scaled_mat))]


ht30 <- Heatmap(scaled_mat, column_split = labs,
               width = unit(45, "mm"),height = unit(200,"mm"),
               row_title = "", column_title = "Mayo", 
               clustering_distance_columns =  "manhattan",
               clustering_method_columns = "ward.D2",
               col = circlize::colorRamp2(c(-4,0,4),
                                          c("#0571b0","#f7f7f7","#ca0020") ),
               top_annotation = HeatmapAnnotation(foo = labs$Diagnosis,
                                                  col = list(foo = (annot_col$Diagnosis)),
                                                  show_legend = T, show_annotation_name = F),
               show_row_names = F,
               show_column_names = F,
               cluster_rows = F,
               show_row_dend = F,
               show_column_dend = F,
               show_heatmap_legend = F,
               column_title_gp = gpar(fill = "white", col = "gray35", border = "white", cex = 2),
               column_title_side = "bottom") 






sel.paths_up <- sel.paths[sel.paths$Direction.x == "UP",]
sel.paths_up <- sel.paths_up[1:10,]$X

sel.paths_dn <- sel.paths[sel.paths$Direction.x == "DOWN",]
sel.paths_dn <- sel.paths_dn[1:10,]$X

sel.paths_all <- c(sel.paths_dn,sel.paths_up)

r_labels <- rownames(scaled_mat)

r_labels[!(r_labels %in% sel.paths_all)] <- ""
r_labels <- gsub("_", " " ,r_labels)
r_labels <- gsub("^([A-Z]*)\\s", "\\1: " ,r_labels)

#r_labels[-grep(" ION|POTASSIUM|NEURO|P38|SYNAP|GABA|YAP1|ECADHERIN NA|IL3|ION REPAIR|NFKB P",r_labels)] <- ""




ht3 <- Heatmap(scaled_mat, column_split = labs$Diagnosis,
               width = unit(45, "mm"),height = unit(200,"mm"),
               row_title = "", column_title = "Mayo", 
               clustering_distance_columns =  "manhattan", clustering_method_columns = "ward.D2",
               col = circlize::colorRamp2(c(-4,0,4),c("#0571b0","#f7f7f7","#ca0020") ),
               top_annotation = HeatmapAnnotation(foo = labs$Diagnosis,
                                                  col = list(foo = annot_col$Diagnosis),
                                                  show_legend = F, show_annotation_name = F),
               show_row_names = F,
               show_column_names = F, cluster_rows = F,
               show_row_dend = F,
               show_column_dend = F,
               show_heatmap_legend = F,
               column_title_gp = gpar(fill = "white", col = "gray25", border = "white", cex = 2),
               column_title_side = "bottom") +
    rowAnnotation(
        link = anno_mark(at = which(r_labels != ""),labels  = r_labels[which(r_labels != "")],
                         labels_gp = gpar(fontsize = 14), padding = unit(2, "mm")))



# r_labels <- rownames(scaled_mat_phg)
# 
# r_labels <- gsub("_", " " ,r_labels)
# r_labels <- gsub("^([A-Z]*)\\s", "\\1: " ,r_labels)
# 
# r_labels[-grep(" ION|POTASSIUM|NEURO|P38|SYNAP|GABA|YAP1|ECADHERIN NA|IL3|ION REPAIR|NFKB P",r_labels)] <- ""
# 
# 

 

ht4 <- Heatmap(scaled_mat_phg, column_split = labs_phg$Diagnosis,
               width = unit(45, "mm"), height = unit(200,"mm"),
               row_title = "", column_title = "MSBB", 
               clustering_distance_columns =  "manhattan",
               clustering_method_columns = "ward.D2",
               col = circlize::colorRamp2(c(-4,0,4),c("#0571b0","#f7f7f7","#ca0020") ),
               top_annotation = HeatmapAnnotation(foo = labs_phg$Diagnosis,
                                                  col = list(foo = annot_col_phg$Diagnosis),
                                                  show_legend = F, show_annotation_name = F),
               show_row_names = F,
               show_column_names = F, cluster_rows = F,
               show_row_dend = F,
               show_column_dend = F,
               show_heatmap_legend = F,
               column_title_gp = gpar(fill = "white", col = "gray25", border = "white", cex = 2),
               column_title_side = "bottom") +
    rowAnnotation(
        link = anno_mark(at = which(r_labels != ""),labels  = r_labels[which(r_labels != "")],
                         labels_gp = gpar(fontsize = 14), padding = unit(2, "mm")))





cairo_pdf(file = paste0("figures/02_d_combined.pdf"),
          width = 20, height = 18, onefile = T)

draw(ht30+ ht4, gap = unit(8, "mm"))

dev.off()
