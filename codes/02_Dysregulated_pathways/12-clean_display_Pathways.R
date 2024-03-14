rm(list = ls())



#*** Cleaning pathway tables ***# #### 
library(stringr)
temp <- list.files("Output",recursive = T,full.names = T,pattern = "diff")
temp <- grep(".csv", temp,value = T)

meta.labels <- list()
big.exp     <- list()

fdr.thresh   <- 0.01
logFC.thresh <- 0

plot.list    <- list()
big.tab.list <- list()

for(i in temp) {
    base.ind <- which(temp == i)
    name.tag <- unlist(str_split(i,"/"))
    name.tag <- tail(name.tag,1)
    name.tag <- unlist(str_split(name.tag,pattern = "_"))
    name.tag <- name.tag[1]
    t.tab    <- read.csv(i)
    t.tab2   <- t.tab %>% dplyr::select(.,X,logFC,t,P.Value,adj.P.Val) %>%
        arrange(.,P.Value) %>%
        dplyr::mutate(Direction = logFC/abs(logFC),
                      Direction = factor(Direction, c(-1,1), c('-1' = 'DOWN', '1' = 'UP')),
                      Direction = as.character(Direction))
    
    
    if(name.tag %in% c("ROSMAP", "UCI")){
        t.tab2$Direction[t.tab2$adj.P.Val > 0.1 | abs(t.tab2$logFC) <
                             logFC.thresh] = 'NONE'
    }else{
        t.tab2$Direction[t.tab2$adj.P.Val > fdr.thresh | abs(t.tab2$logFC) <
                             logFC.thresh] = 'NONE'
    }
    
    t.tab2$Direction <- factor(t.tab2$Direction)
    t.readout        <- as.data.frame.matrix(t(table(t.tab2$Direction)))
    t.readout        <- data.frame(t.readout)
    t.readout$assay  <- name.tag
    tryCatch(expr = {big.tab.list  <- rbind(big.tab.list,t.readout)},
             error = function(e){print(i);})
    
    
    t.address   <- unlist(strsplit(i, "/"))
    t.address   <- tail(t.address,1) 
    t.address   <- paste0("tables/DE_pathways_clean/",t.address)
    
    
    
    t.tab2$X <- gsub("Pathway.","",t.tab2$X)
    
    display_names    <- read.csv("preprocessed/MSigDB_display_names.csv")
    rownames(display_names) <- display_names$STANDARD_NAME
    

    t.tab2$display_name <- display_names[t.tab2$X,]$display_name
    
    
    t.tab3 <- mutate_if(t.tab2,is.numeric,signif,digits = 3)
    colnames(t.tab3) <- c("MSigDB name",
                          "logFC",
                          "t","P.Value",
                          "adj.P.Val",
                          "Direction (AD vs Control)",
                          "Pathway Description")

    t.tab3 <- t.tab3[, c("Pathway Description",
                         "MSigDB name",
                         "Direction (AD vs Control)",
                         "logFC",
                         "t","P.Value",
                         "adj.P.Val")]
    write.csv(t.tab3,t.address,row.names = F)
    print(t.readout)
    
}

