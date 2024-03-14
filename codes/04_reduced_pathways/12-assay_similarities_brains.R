rm(list = ls())
library(stringr)
library(dplyr)
library(ggpubr)






# Selecting all DEP from the I45F/I47F assay

temp        <- list.files("tables/DE_pathways_reduced/",recursive = T,
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
                    dplyr::select(.,X,Direction,logFC,P.Value, adj.P.Val) %>%
                    mutate(.,assay = name.tag) %>%
                    filter(.,Direction != "NONE") %>%
                    mutate(.,path_dir = paste0(Direction,": ", X))
                    
    big.tab.list[[i]] <- t.tab2
    big.tab.list2[[i]] <- t.tab
}


 
 display_names    <- read.csv("preprocessed/MSigDB_display_names.csv")
 rownames(display_names) <- display_names$STANDARD_NAME
 
 
 
 Core_paths0 <- Reduce(function(x, y) merge(x, y, by = c("X","Direction"),all = F),
                       big.tab.list[2:3])
 
 Core_paths0$display_name <- display_names[Core_paths0$X,]$display_name
 
 write.csv(Core_paths0,"tables/shared_pathways_reduced/01_shared_pathways_MSBB_Mayo.csv")
 
 Core_paths0$X<- gsub("_", " " ,Core_paths0$X)
 Core_paths0$X<- gsub("^([A-Z]*)\\s", "\\1: " ,Core_paths0$X)
 write.csv(Core_paths0,"tables/shared_pathways_reduced/01_shared_pathways_MSBB_Mayo_clean.csv")
 
 
 
 Core_paths0 <- Reduce(function(x, y) merge(x, y, by = c("X", "Direction"),all = F),
                       big.tab.list[1:2])
 Core_paths0$display_name <- display_names[Core_paths0$X,]$display_name
 
 write.csv(Core_paths0,"tables/shared_pathways_reduced/01_shared_pathways_I45vsI47_Mayo.csv")
 
 Core_paths0$X<- gsub("_", " " ,Core_paths0$X)
 Core_paths0$X<- gsub("^([A-Z]*)\\s", "\\1: " ,Core_paths0$X)
 write.csv(Core_paths0,"tables/shared_pathways_reduced/01_shared_pathways_I45vsI47_Mayo_clean.csv")
 
 
 
 
 
 Core_paths0 <- Reduce(function(x, y) merge(x, y, by = c("X","Direction"),all = F),
                       big.tab.list[1:3])
 Core_paths0$display_name <- display_names[Core_paths0$X,]$display_name
 
 write.csv(Core_paths0,"tables/shared_pathways_reduced/01_shared_pathways_I45FvsI47_MSBB_Mayo.csv")
 
 
 
 
 Core_paths0$X<- gsub("_", " " ,Core_paths0$X)
 Core_paths0$X<- gsub("^([A-Z]*)\\s", "\\1: " ,Core_paths0$X)
 
 write.csv(Core_paths0,"tables/shared_pathways_reduced/01_shared_pathways_I45FvsI47_MSBB_Mayo_clean.csv")
 
 
 
 


Core_paths_all <- Reduce(function(x, y) merge(x, y, by = c("path_dir"),all = FALSE),
                     big.tab.list)

Core_paths_all$display_name <- display_names[Core_paths_all$X.x,]$display_name

write.csv(Core_paths_all,file = "tables/shared_pathways_reduced/45vs47_Mayo_MSBB_ROSMAP_2.csv")

Core_paths_all$X.x <- gsub("_", " " ,Core_paths_all$X.x)
Core_paths_all$X.x <- gsub("^([A-Z]*)\\s", "\\1: " ,Core_paths_all$X.x)
write.csv(Core_paths_all,file = "tables/shared_pathways_reduced/45vs47_Mayo_MSBB_ROSMAP_2_clean.csv")


Core_paths_brain <- Reduce(function(x, y) merge(x, y, by = c("path_dir"),all = FALSE),
                         big.tab.list[2:4])



Core_paths_brain$X.x <- gsub("_", " " ,Core_paths_brain$X.x)
Core_paths_brain$X.x <- gsub("^([A-Z]*)\\s", "\\1: " ,Core_paths_brain$X.x)
write.csv(Core_paths_brain,file = "tables/shared_pathways_reduced/brain_Mayo_MSBB_ROSMAP_2_clean.csv")

