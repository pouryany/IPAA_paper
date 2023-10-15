rm(list = ls())
library(stringr)
library(dplyr)
library(ggpubr)






# Selecting all DEP from the I45F/I47F assay

temp        <- list.files("tables/DE_pathways2",recursive = T,
                          full.names = T)

temp        <- grep("ROS|Mayo|MSBB|I45FvsI",temp,value = T)


meta.labels <- list()
big.exp     <- list()

fdr.thresh   <- 0.01
logFC.thresh <- 0

plot.list    <- list()
big.tab.list <- list()
big.tab.list2 <- list()

for(i in temp) {
    base.ind <- which(temp == i)
    name.tag <- unlist(str_split(i,"/"))
    name.tag <- tail(name.tag,1)
    name.tag <- unlist(str_split(name.tag,pattern = "_"))
    name.tag <- name.tag[1]
    t.tab    <- read.csv(i)
    t.tab2   <- t.tab %>% 
                    dplyr::select(.,X,Direction,logFC,P.Value) %>%
                    mutate(.,assay = name.tag) %>%
                    filter(.,Direction != "NONE") %>%
                    mutate(.,path_dir = paste0(Direction,": ", X))
                    
    big.tab.list[[i]] <- t.tab2
    big.tab.list2[[i]] <- t.tab
}



Core_paths0 <- Reduce(function(x, y) merge(x, y, by = c("X"),all = T),
                     big.tab.list[2:3])


sign(Core_paths0$logFC.x) * -log(Core_paths0$P.Value.x)

Core_paths0$Mayo <- sign(Core_paths0$logFC.x) * -log(Core_paths0$P.Value.x)
Core_paths0$MSBB <- sign(Core_paths0$logFC.y) * -log(Core_paths0$P.Value.y)

Core_paths0  <- Core_paths0[,c("X","Mayo", "MSBB")]
colnames(Core_paths0) <- c("X","TCX", "PHG")
Core_paths0[is.na(Core_paths0)] <- 0
rownames(Core_paths0) <- Core_paths0$X


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
                               heatmap_legend_param = list(direction = "horizontal",title_gp = gpar(fontsize = 10, fontface = "bold"))
)


 cairo_pdf(file = "figures/02_c_heatmap_brains.pdf",
           width = 2, height = 10, onefile = T)
 ComplexHeatmap::draw(ht1, heatmap_legend_side = "bottom")
 dev.off()
 
 
 
 Core_paths0 <- Reduce(function(x, y) merge(x, y, by = c("X"),all = F),
                       big.tab.list[2:3])
 write.csv(Core_paths0,"tables/shared_pathways/01_shared_pathways_MSBB_Mayo.csv")
 
 Core_paths0$X<- gsub("_", " " ,Core_paths0$X)
 Core_paths0$X<- gsub("^([A-Z]*)\\s", "\\1: " ,Core_paths0$X)
 write.csv(Core_paths0,"tables/shared_pathways/01_shared_pathways_MSBB_Mayo_clean.csv")
 
 
 
 Core_paths0 <- Reduce(function(x, y) merge(x, y, by = c("X", "Direction"),all = F),
                       big.tab.list[1:2])
 write.csv(Core_paths0,"tables/shared_pathways/01_shared_pathways_I45vsI47_Mayo.csv")
 
 Core_paths0$X<- gsub("_", " " ,Core_paths0$X)
 Core_paths0$X<- gsub("^([A-Z]*)\\s", "\\1: " ,Core_paths0$X)
 write.csv(Core_paths0,"tables/shared_pathways/01_shared_pathways_I45vsI47_Mayo_clean.csv")
 
 
 
 
 
 Core_paths0 <- Reduce(function(x, y) merge(x, y, by = c("X","Direction"),all = F),
                       big.tab.list[1:3])
 write.csv(Core_paths0,"tables/shared_pathways/01_shared_pathways_I45FvsI47_MSBB_Mayo.csv")
 
 
 
 Core_paths0$X<- gsub("_", " " ,Core_paths0$X)
 Core_paths0$X<- gsub("^([A-Z]*)\\s", "\\1: " ,Core_paths0$X)
 
 
 write.csv(Core_paths0,"tables/shared_pathways/01_shared_pathways_I45FvsI47_MSBB_Mayo.csv")
 
 
 
 
 
Core_paths <- Reduce(function(x, y) merge(x, y, by = c("X","Direction"),all = F),
                big.tab.list[2:3])

Core_paths_all <- Reduce(function(x, y) merge(x, y, by = c("path_dir"),all = FALSE),
                     big.tab.list)


write.csv(Core_paths_all,file = "tables/shared_pathways/45vs47_Mayo_MSBB_ROSMAP_2.csv")

Core_paths_all$X.x <- gsub("_", " " ,Core_paths_all$X.x)
Core_paths_all$X.x <- gsub("^([A-Z]*)\\s", "\\1: " ,Core_paths_all$X.x)
write.csv(Core_paths_all,file = "tables/shared_pathways/45vs47_Mayo_MSBB_ROSMAP_2_clean.csv")


Core_paths_brain <- Reduce(function(x, y) merge(x, y, by = c("path_dir"),all = FALSE),
                         big.tab.list[2:4])



Core_paths_brain$X.x <- gsub("_", " " ,Core_paths_brain$X.x)
Core_paths_brain$X.x <- gsub("^([A-Z]*)\\s", "\\1: " ,Core_paths_brain$X.x)
write.csv(Core_paths_brain,file = "tables/shared_pathways/brain_Mayo_MSBB_ROSMAP_2_clean.csv")



temp <- merge(big.tab.list2[[2]],
              big.tab.list2[[3]],
              by = "X")


temp$dirx <- plyr::mapvalues(temp$Direction.x,c("UP","DOWN","NONE"), c(1,-1,0))
temp$diry <- plyr::mapvalues(temp$Direction.y,c("UP","DOWN","NONE"), c(1,-1,0))

temp$Direction.both <- as.factor(as.numeric(temp$dirx) * as.numeric(temp$diry))

temp$Direction.both <- plyr::mapvalues(temp$Direction.both, c(0,1), c("NONE","Significant in both\nMayo and MSBB"))

# temp    <- temp  %>% mutate(.,t.x = sign(logFC.x) * (-log(P.Value.x)))
# temp    <- temp  %>% mutate(.,t.y = sign(logFC.y) * (-log(P.Value.y)))

# 
# temp    <- temp  %>% mutate(.,t.x = ifelse(sign(logFC.x) > 0,  qnorm(1-P.Value.x/2,lower.tail = TRUE),
#                                            qnorm(1-P.Value.x/2,lower.tail = FALSE)))
# temp    <- temp  %>% mutate(.,t.y = ifelse(sign(logFC.y) > 0,  qnorm(1-P.Value.y/2,lower.tail = TRUE),
#                                            qnorm(1-P.Value.y/2,lower.tail = FALSE)))
# 


corrs <- signif(cor(temp$t.x,temp$t.y,method = "spearman"),2)
 


ggplot(data = temp, aes(t.x, t.y)) +
  geom_point(alpha = 0.5, size = 4) + 
  theme_classic() +
  labs (x = " Mayo AD vs control (t-score)",
        y = "MSBB AD vs control (t-score)",
        title = "Pathway dysregulation" ) +
  scale_color_manual(values = c("black", "#b2182b")) +
  theme(plot.title = element_text(size = 32 , hjust = 0.5, face = "plain"),
        axis.title = element_text(size = 24,
                                    angle = 0,
                                    vjust = 0.5, hjust = 0.5,
                                    face = "plain"),
        axis.text = element_text(size = 20), 
        legend.text = element_text(face = "plain"),
        legend.position = "none") +
  annotate("text",x=  -4, y = 5.5,
           label = paste0("Spearman's rho = ", corrs), size = 8) 

ggsave("figures/02_f__correlation.pdf",width = 8, height = 8)



temp2 <- table(temp$Direction.x , temp$Direction.y)




library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(4, "Pastel2")

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

a <- venn.diagram(
  x = list( big.tab.list[[2]]$path_dir,
            big.tab.list[[3]]$path_dir),
  category.names = c("Mayo" , "MSBB"),
  filename = NULL,
  output=F,
  # Numbers,
  col  = c("#FFB000", "#80B1D3"),#myCol[1:3],
  fill = c("#FFB000", "#80B1D3"),#myCol[1:3],
  cex = 6,
  fontface = "plain",
  fontfamily = "sans",
  cat.cex = 6,
  # Set names
  cat.fontface = "plain",
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  main.cex = 2.0,
  main.fontfamily = "sans",
  cat.dist = c(0.035, 0.035) ,
  cat.pos = c(-30, 30) # Modified

)


ggsave("figures/02_c_assays_venn.PDF",plot =a, width = 12, height = 12)


b <- venn.diagram(
  x = list( big.tab.list[[2]]$path_dir,
            big.tab.list[[3]]$path_dir,
            big.tab.list[[1]]$path_dir),
  category.names = c("Mayo" , "MSBB", "A\u03B242-H vs A\u03B240-H"),
  filename = NULL,
  output=F,
  # Numbers,
  col  = c("#FFB000", "#80B1D3", "#beaed4"),
  fill = c("#FFB000", "#80B1D3","#beaed4"),
  cex = 5,
  fontface = "plain",
  fontfamily = "sans",
  cat.cex = 5,
  # Set names
  cat.fontface = "plain",
  cat.default.pos = "outer",
  cat.fontfamily = "sans",
  main.cex = 2.0,
  main.fontfamily = "sans"#,
  # cat.dist = c(0.035, 0.035) ,
  # cat.pos = c(-30, 30) # Modified
  
)




cairo_pdf(file = "figures/04_a_assays__models_venn.PDF",
          width = 14,
          height = 14,
          onefile = T)
grid.draw(b)
dev.off()

