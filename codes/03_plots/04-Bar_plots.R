
rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ComplexHeatmap)

Core_paths0   <- read.csv("tables/shared_pathways/01_shared_pathways_I45FvsI47_MSBB_Mayo.csv", row.names = 1)
display_names <- read.csv("preprocessed/MSigDB_display_names.csv")
rownames(display_names) <- display_names$STANDARD_NAME

display_names <- display_names[Core_paths0$X,]

if(all(display_names$STANDARD_NAME == Core_paths0$X)){
  Core_paths0$X <- display_names$display_name
}





Core_paths0$I45FvsI47F <- sign(Core_paths0$logFC.x) * -log(Core_paths0$P.Value.x)
Core_paths0$Mayo <- sign(Core_paths0$logFC.y) * -log(Core_paths0$P.Value.y)
Core_paths0$MSBB <- sign(Core_paths0$logFC) * -log(Core_paths0$P.Value)

Core_paths0 <- Core_paths0[order(Core_paths0$Mayo),]
Core_paths0_up <- Core_paths0 %>% dplyr::filter(.,Direction == "UP") %>% dplyr::select(.,X, I45FvsI47F, Mayo, MSBB)

colnames(Core_paths0_up)[2] <- "A\u03B242-H vs A\u03B240-H"



rownames(Core_paths0_up) <- Core_paths0_up$X
Core_paths0_up$X <- NULL
Core_paths0_up <- Core_paths0_up[order(Core_paths0_up$Mayo,decreasing = T),c(2,3,1)]

set.seed(9)
ht5 <- ComplexHeatmap::Heatmap(Core_paths0_up, cluster_rows = F, cluster_columns = F,
                               col = circlize::colorRamp2(c(0,1,5,10, 30), 
                                                          colors = c("#fff5f0",
                                                                     "#fcbba1",
                                                                     "#fb6a4a",
                                                                     "#cb181d",
                                                                     "#67000d")),
                               # name = "-log10(pval)",
                               name = "Dysregulation\nscore",
                               width  = unit(0.6 * ncol(Core_paths0_up),"cm"),
                               height = unit(0.45 * nrow(Core_paths0_up),"cm"),
                               column_names_side = "bottom",
                               row_names_side = "left",
                               na_col = "grey95",
                               rect_gp = gpar(col = "white", lwd = 2),
                               border_gp = gpar(col = "black", lty = 1),
                               heatmap_legend_param = list(direction = "horizontal")
)





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
  

# 
# Core_paths0 <- read.csv("tables/shared_pathways/01_shared_pathways_I45FvsI47_MSBB_Mayo.csv", row.names = NULL)
# 
# Core_paths0$I45FvsI47F <- sign(Core_paths0$logFC.x) * -log(Core_paths0$P.Value.x)
# Core_paths0$Mayo <- sign(Core_paths0$logFC.y) * -log(Core_paths0$P.Value.y)
# Core_paths0$MSBB <- sign(Core_paths0$logFC) * -log(Core_paths0$P.Value)

Core_paths0 <- Core_paths0[order(Core_paths0$Mayo),]





Core_paths0_down <- Core_paths0 %>% dplyr::filter(.,Direction == "DOWN") %>%
  dplyr::select(.,X, I45FvsI47F, Mayo, MSBB)

colnames(Core_paths0_down)[2] <- "A\u03B242-H vs A\u03B240-H"




# Fix display names 




if(any(str_length(Core_paths0_down$X) > 70)){
  longInds <- str_length(Core_paths0_down$X) > 70
  Core_paths0_down[longInds,]$X <- stringr::str_sub(Core_paths0_down[longInds,]$X,1,70)
  Core_paths0_down[longInds,]$X <- paste0(Core_paths0_down[longInds,]$X,"*") 
}



rownames(Core_paths0_down) <- Core_paths0_down$X
Core_paths0_down$X <- NULL


Core_paths0_down <- Core_paths0_down[,c(2,3,1)]


set.seed(9)
ht6 <- ComplexHeatmap::Heatmap(Core_paths0_down, cluster_rows = F, cluster_columns = F,
                               col = circlize::colorRamp2(c(0,-1,-5,-10,-30), 
                                                          colors = c("#eff3ff",
                                                                     "#bdd7e7",
                                                                     "#6baed6",
                                                                     "#3182bd",
                                                                     "#08306b")),
                               # name = "-log10(pval)",
                               name = "Dysregulation\nscore",
                               width  = unit(0.6 * ncol(Core_paths0_down),"cm"),
                               height = unit(0.45 * nrow(Core_paths0_down),"cm"),
                               column_names_side = "bottom", row_names_side = "right",
                               na_col = "grey95",
                               rect_gp = gpar(col = "white", lwd = 2),
                               border_gp = gpar(col = "black", lty = 1),
                               heatmap_legend_param = list(direction = "horizontal")
)






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


cairo_pdf("figures/04_b_revised_shared_heatmap.pdf",
          width = 15,height = 15,onefile = T)
draw(ht5,heatmap_legend_side = "bottom")
draw(ht6,heatmap_legend_side = "bottom")
dev.off()





# PCxN Mapping.
library(PanomiR)
library(RColorBrewer)
library(igraph)
library(qvalue)



pcxn.dir   <- paste0("preprocessed/PCxN_MSigDB_withJaccard.RDS")
PCXN       <- readRDS(pcxn.dir)

pathways0  <- 
  read.csv("tables/shared_pathways/01_shared_pathways_I45FvsI47_MSBB_Mayo.csv",
           row.names = 1)

pathways    <- read.csv("Output/MAYO/Mayo_ADvsControl_diffPathways.csv",
                        row.names = 1)

pathways0$pathways   <- gsub("^.*: ","Pathway.",pathways0$path_dir)
pathways <- pathways[rownames(pathways) %in% pathways0$pathways,]


library(igraph)

igraph::cluster_edge_betweenness


correlationCutOff = 0.2 


temp <- PanomiR::mappingPathwaysClusters(PCXN, 
                             pathways, 
                             correlationCutOff = correlationCutOff,
                             clusteringFunction = "cluster_edge_betweenness",
                             pathwayFDR = 1,
                             outDir = "Output/PCxN/",
                             topClusters = 6)

temp$Clustering
PanomiR::clusterPlot
pathways <- pathways[temp$Clustering$Pathway,]

pathways$cluster <- temp$Clustering$cluster
pathways$pathway <- rownames(pathways)



pathways$pathway <- gsub("Pathway.","",pathways$pathway)


display_names    <- read.csv("preprocessed/MSigDB_display_names.csv")
rownames(display_names) <- display_names$STANDARD_NAME

display_names <- display_names[pathways$pathway,]

if(all(display_names$STANDARD_NAME == pathways$pathway)){
  pathways$display_names<- display_names$display_name
}

write.csv(pathways, "Output/PCxN/Figures/PCxN_table.csv")





AD_subnetwork <- PCXN[PCXN$Pathway.A %in% rownames(pathways) &
                      PCXN$Pathway.B %in% rownames(pathways) & 
                      PCXN$p.Adjust < 0.05 &
                      abs(PCXN$PathCor) > correlationCutOff, ]
  





write.csv(AD_subnetwork, "Output/PCxN/Figures/PCxN_network.csv")



display_names    <- read.csv("preprocessed/MSigDB_display_names.csv")
rownames(display_names) <- display_names$STANDARD_NAME

AD_subnetwork$Pathway.A <- gsub("Pathway.","",AD_subnetwork$Pathway.A)
AD_subnetwork$Pathway.B <- gsub("Pathway.","",AD_subnetwork$Pathway.B)

AD_subnetwork$Pathway_A_description <- display_names[AD_subnetwork$Pathway.A,]$display_name
AD_subnetwork$Pathway_B_description <- display_names[AD_subnetwork$Pathway.B,]$display_name




AD_subnetwork$PathCor  <- signif(AD_subnetwork$PathCor,3)
AD_subnetwork$p.Adjust  <- signif(AD_subnetwork$p.Adjust,3)
AD_subnetwork$Jaccard.Ind  <- signif(AD_subnetwork$Jaccard.Ind,3)
AD_subnetwork[AD_subnetwork$p.Adjust < 2.2e-16,]$p.Adjust  <- 2.2e-16
write.csv(AD_subnetwork, "Output/PCxN/Figures/cleaned_ PCxN_network.csv")





AD_PCxN <- temp$DE_PCXN
V(AD_PCxN)$shape




set.seed(4)

new_cluster_plot <- function (subNet, subplot = FALSE, topClusters = 2,
                              prefix = "", 
                              outDir = ".", plotSave = TRUE,asp = 0.5) 
{
  if (!dir.exists(outDir)) {
    stop("Output directory does not exist.")
  }
  figDir <- paste0(outDir, "/", prefix)
  legendCats <- data.frame(attr = c("Up-regulated", "Down-regulated"), 
                           shape = unique(igraph::V(subNet)$shape))
  cols <- RColorBrewer::brewer.pal(8, "Pastel2")
  clstMems <- igraph::V(subNet)$cluster
  nodeColors <- cols[clstMems]
  smallClust <- which(table(clstMems) <= table(clstMems)[5], 
                      useNames = TRUE)
  nodeColors[clstMems %in% smallClust] <- NA
  clusterPlotHelper(plotSave, figDir, subNet, nodeColors, 
                    legendCats,asp = asp)
  clusterSubPlotHelper(subplot, topClusters, clstMems, subNet, 
                       figDir, cols, legendCats,asp = asp)
  "%notin%" <- Negate("%in%")
  remove <- which(table(clstMems) < 5)
  remove <- which((clstMems %notin% remove))
  subNet2 <- igraph::induced_subgraph(subNet, remove)
  if (subplot == TRUE) {
    grDevices::pdf(paste0(figDir, "ConnectedPathways_PCxNCorGraph.pdf"), 
                   width = 18, height = 11)
    plot(subNet2, edge.width = 1.3, vertex.size = 4, vertex.label = NA, 
         vertex.color = cols[clstMems[remove]], legend = TRUE, 
         layout = igraph::layout_with_fr, asp = asp)
    graphics::legend(x = "bottomright", y = 200, legend = legendCats$attr, 
                     pch = c(0, 1), bty = "n", cex = 1.4)
    graphics::legend(x = "topleft", legend = c("Positive Cor", 
                                               "Negative Cor"),
                     col = c("#E41A1C", "#377EB8"), lty = 1, 
                     lwd = 2, cex = 1.4, bty = "n")
    grDevices::dev.off()
  }
}


clusterPlotHelper <- function(plotSave, figDir, subNet, nodeColors,
                              legendCats, asp = 0.5) {
  if (plotSave) {
    grDevices::pdf(paste0(figDir, "PCxNCorGraph.pdf"),
                   width = 18, height = 11)
  }
  plot(subNet, vertex.size = 4, vertex.label = NA,
       vertex.color = nodeColors, asp = asp)
  graphics::legend(x = "bottomleft", legend = legendCats$attr, pch = c(0, 1),
                   bty = "n", cex = 1.6)
  graphics::legend(x = "topleft", legend = c("Positive Cor", "Negative Cor"),
                   col = c("#E41A1C", "#377EB8"), lty = 1, lwd = 2,
                   cex = 1.6, bty = "n")
  if (plotSave)
    grDevices::dev.off()
}

clusterSubPlotHelper <- function(subplot, topClusters, clstMems, subNet,
                                 figDir, cols, legendCats, asp = 0.5) {
  if (subplot == TRUE) {
    for (k in seq_len(topClusters)) {
      keep <- which((clstMems) == k)
      subNet2 <- igraph::induced_subgraph(subNet, keep)
      if (length(igraph::V(subNet2)) < 2) next
      
      grDevices::pdf(paste0(figDir, "PCxNCorGraph_",
                            "Cluster_", k, ".pdf"))
      plot(subNet2, edge.width = 1.3, vertex.size = 4, vertex.label = NA,
           vertex.color = cols[clstMems[keep]], legend = TRUE,
           layout = igraph::layout_with_fr,
           asp = asp)
      graphics::legend(x = "bottomleft", legend = legendCats$attr,
                       pch = c(0, 1), bty = "n", cex = 1.4)
      graphics::legend(x = "topleft",
                       legend = c("Positive Cor", "Negative Cor"),
                       col = c("#E41A1C", "#377EB8"),
                       lty = 1, lwd = 2, cex = 1.4, bty = "n")
      grDevices::dev.off()
    }
  }
}





new_cluster_plot(AD_PCxN,subplot = T, plotSave = T, topClusters = 4,
                 outDir = "Output/PCxN", prefix = "customized")

# 
# 
# # Alternative stuff
# 
# V(AD_PCxN)$dir <- plyr::mapvalues(V(AD_PCxN)$shape,
#                                   c("circle","square"),
#                                   c("#08306b", "#67000d"))
# 
# 
# V(AD_PCxN)$cluster_shape <- plyr::mapvalues(V(AD_PCxN)$cluster,
#                                             c(1:4),
#                                             c("circle", "square","triangle","diamond"))
# 
# 
# 
# "%notin%" <- Negate("%in%")
# remove <- which(table(V(AD_PCxN)$cluster) < 5)
# remove <- which(V(AD_PCxN)$cluster %notin% remove)
# subNet2 <- igraph::induced_subgraph(AD_PCxN, remove)
# 
# 
# 
# legendCats <- data.frame(attr =  c(1:2), 
#                          shape = unique(igraph::V(subNet2)$cluster_shape))
# 
# redBu <- colorRampPalette(c("#08306b", "#67000d"))
# 
# colFun <- circlize::colorRamp2(c(-0.65,0,.5), 
#                      colors = c("#08306b","grey65", "#67000d"))
# 
# 
# 
# 
# 
# plot(subNet2, edge.width = 1.3, vertex.size = 5,
#      vertex.label = NA, vertex.shape =  V(subNet2)$cluster_shape,
#      vertex.color = V(subNet2)$dir, 
#      legend = TRUE,
#      #edge.color = colFun(as.numeric(E(subNet2)$PathCor)),
#      layout = igraph::layout_components)
# 
# graphics::legend(x = "bottomright", y = 200, legend = legendCats$attr, 
#                  pch = c(0, 1), bty = "n", cex = 1.4)
# graphics::legend(x = "topleft", legend = c("Positive Cor", 
#                                            "Negative Cor"),
#                  col = c("#E41A1C", "#377EB8"), lty = 1, 
#                  lwd = 2, cex = 1.4, bty = "n")
# 
# 
# 
# 
# new_cluster_plot2 <- function (subNet, subplot = FALSE, topClusters = 2,
#                               prefix = "", 
#                               outDir = ".", plotSave = TRUE) 
# {
#   if (!dir.exists(outDir)) {
#     stop("Output directory does not exist.")
#   }
#   figDir <- paste0(outDir, "/", prefix)
#   legendCats <- data.frame(attr =  c(1:4), 
#                            shape = unique(igraph::V(subNet)$cluster_shape))
#   cols <- RColorBrewer::brewer.pal(8, "Pastel2")
#   clstMems <- igraph::V(subNet)$cluster
#   nodeColors <- V(AD_PCxN)$dir
#   smallClust <- which(table(clstMems) <= table(clstMems)[5], 
#                       useNames = TRUE)
#   #nodeColors[clstMems %in% smallClust] <- NA
#   clusterPlotHelper(plotSave, figDir, subNet, nodeColors, 
#                     legendCats)
#   clusterSubPlotHelper(subplot, topClusters, clstMems, subNet, 
#                        figDir, cols, legendCats)
#   "%notin%" <- Negate("%in%")
#   remove <- which(table(clstMems) < 5)
#   remove <- which((clstMems %notin% remove))
#   subNet2 <- igraph::induced_subgraph(subNet, remove)
#   if (subplot == TRUE) {
#     grDevices::pdf(paste0(figDir, "ConnectedPathways_PCxNCorGraph.pdf"), 
#                    width = 18, height = 11)
#     plot(subNet2, edge.width = 1.3, vertex.size = 5, vertex.label = NA, 
#          vertex.color = cols[clstMems[remove]], legend = TRUE, 
#          layout = igraph::layout_components)
#     graphics::legend(x = "bottomright", y = 200, legend = legendCats$attr, 
#                      pch = c(0, 1), bty = "n", cex = 1.4)
#     graphics::legend(x = "topleft", legend = c("Positive Cor", 
#                                                "Negative Cor"),
#                      col = c("#E41A1C", "#377EB8"), lty = 1, 
#                      lwd = 2, cex = 1.4, bty = "n")
#     grDevices::dev.off()
#   }
# }
# 
# 
# clusterPlotHelper <- function(plotSave, figDir, subNet, nodeColors,
#                               legendCats) {
#   if (plotSave) {
#     grDevices::pdf(paste0(figDir, "PCxNCorGraph.pdf"),
#                    width = 18, height = 11)
#   }
#   plot(subNet, vertex.size = 5, vertex.label = NA, vertex.color = nodeColors)
#   graphics::legend(x = "bottomleft", legend = legendCats$attr, pch = c(0, 1),
#                    bty = "n", cex = 1.6)
#   graphics::legend(x = "topleft", legend = c("Positive Cor", "Negative Cor"),
#                    col = c("#E41A1C", "#377EB8"), lty = 1, lwd = 2,
#                    cex = 1.6, bty = "n")
#   if (plotSave)
#     grDevices::dev.off()
# }
# 
# clusterSubPlotHelper <- function(subplot, topClusters, clstMems, subNet,
#                                  figDir, cols, legendCats) {
#   if (subplot == TRUE) {
#     for (k in seq_len(topClusters)) {
#       keep <- which((clstMems) == k)
#       subNet2 <- igraph::induced_subgraph(subNet, keep)
#       if (length(igraph::V(subNet2)) < 2) next
#       
#       grDevices::pdf(paste0(figDir, "PCxNCorGraph_",
#                             "Cluster_", k, ".pdf"))
#       plot(subNet2, edge.width = 1.3, vertex.size = 5, vertex.label = NA,
#            vertex.color = cols[clstMems[keep]], legend = TRUE,
#            layout = igraph::layout.fruchterman.reingold)
#       graphics::legend(x = "bottomleft", legend = legendCats$attr,
#                        pch = c(0, 1), bty = "n", cex = 1.4)
#       graphics::legend(x = "topleft",
#                        legend = c("Positive Cor", "Negative Cor"),
#                        col = c("#E41A1C", "#377EB8"),
#                        lty = 1, lwd = 2, cex = 1.4, bty = "n")
#       grDevices::dev.off()
#     }
#   }
# }
# 
# 
# 
# new_cluster_plot2(AD_PCxN,subplot = T, plotSave = T, topClusters = 4,
#                  outDir = "Output/PCxN", prefix = "customized")
# 
# 
# 
# 
# 
