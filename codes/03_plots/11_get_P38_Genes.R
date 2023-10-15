# Are there differentially expressed genes in P38  activation
# Are they shared between assays?
# Can you make bar plots?
# Can you make heatmap?

# Let's collect all de genes 



rm(list = ls())


library(tidyverse)
all.kat.files <- list.files("tables/DE_genes_filtered_Mjd", full.names = T, recursive = T, pattern = "\\.csv|\\.txt")

all.kat.files <- setdiff(all.kat.files,grep("-",all.kat.files,value = T))
all.kat.files <- setdiff(all.kat.files,grep("G2B2",all.kat.files,value = T))

keep_genes <- readRDS("preprocessed/gene_meta_data.RDS")

head(keep_genes)




# The following codes are not provided

# Directory to the pathway gene membership
pathways2 <- readRDS('preprocessed/MSigDBPathGeneTab.RDS')


p38inds  <- grep("ACTIVATED_TAK1", pathways2$Pathway)
p38genes <- pathways2[p38inds,]
p38genes <- keep_genes[keep_genes$ID %in% p38genes$ENSEMBL,]

p38jubilee <- list()
for(de_file in all.kat.files){
    pour.paths <- read_csv(de_file)
    
    nameTag    <- tail(unlist(str_split(de_file, "/")),1)
    nameTag    <- gsub(".csv", "",nameTag)
    nameTag    <- gsub("_.*", "",nameTag)
    
    these_p38s <- pour.paths[pour.paths$ID %in% p38genes$ID,]
    these_p38s$tag <- nameTag
    p38jubilee[[nameTag]] <- these_p38s
    
}

p38jubilee <- Reduce(rbind, p38jubilee)

p38jubilee$score <- -log10(p38jubilee$FDR) * sign(p38jubilee$logFC)


p38jubilee <- dplyr::select(p38jubilee,Name,tag, score)


p38jubilee <- tidyr::pivot_wider(p38jubilee,names_from = tag, values_from = score)

p38jubilee <- as.data.frame(p38jubilee)
rownames(p38jubilee) <- p38jubilee$Name
p38jubilee$Name <- NULL
colnames(p38jubilee) <- c("A\u03B242-H vs A\u03B240-H",
                          "Mayo AD vs Ctrl",
                          "MSBB AD vs Ctrl",
                          "ROSMAP AD vs Ctrl")


set.seed(10)
cn <- colnames(p38jubilee)

library(ComplexHeatmap)
ht5 <- ComplexHeatmap::Heatmap(p38jubilee, 
                               col = circlize::colorRamp2(c(-20,-4,-2,0,2,4,20),
                                                          c("#053061",
                                                            "#0571b0",
                                                            "#4393c3",
                                                            "#f7f7f7",
                                                            "#d6604d",
                                                            "#ca0020",
                                                            "#67001f") ),
                               # col = circlize::colorRamp2(c(-10,-3,-0.7,0,0.7,3,10), 
                               #                            # c(-6,
                               #                            #                          -3,
                               #                            #                          -1,
                               #                            #                          0,
                               #                            #                          1,
                               #                            #                          3,
                               #                            #                          6),
                               #                            # c(-20, -5,0, 10, 20),
                               #                            colors = c("darkblue",
                               #                                       "cornflowerblue",
                               #                                       "grey95",
                               #                                       "grey95",
                               #                                       "grey95",
                               #                                       "coral",
                               #                                       "tomato2")),
                               # name = "-log10(pval)",
                               name = "Differential\nexpression",
                               width = unit(1, "cm")*ncol(p38jubilee), 
                               height = unit(0.8, "cm")*nrow(p38jubilee),
                               row_title = c(""),
                               column_title = c(""),
                               show_column_dend = F, 
                               show_row_dend = F,
                               show_column_names = F,
                               bottom_annotation = HeatmapAnnotation(
                                   text = anno_text(cn, rot = 45,
                                                    offset = unit(1, "npc"),
                                                    just = "right")),
                               
                               heatmap_legend_param = list(direction = "horizontal"),
                               cell_fun = function(j, i, x, y, w, h, fill) {
                                   if(abs(p38jubilee[i, j]) > 2) {
                                       grid.text("*", x, y, 
                                                 gp = gpar(fontsize = 24, fontface = "bold"))
                                   }}
                               )
draw(ht5)


cairo_pdf(file = paste0("figures/05_a_P38_genes.pdf"),
          width = 10, height = 10, onefile = T)

draw(ht5, heatmap_legend_side = "bottom")

dev.off()






