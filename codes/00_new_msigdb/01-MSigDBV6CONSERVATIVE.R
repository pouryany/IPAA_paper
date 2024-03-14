rm(list = ls())
library(msigdb)
library(ExperimentHub)
library(GSEABase)

library(foreach)
library(doParallel)
detectCores()

registerDoParallel(14)

pathways2 <- readRDS('preprocessed/MSigDB.RDS')

jaccard.ind <- function(set1,set2){
    uns <- length(union(set1,set2))
    ovr <- length(intersect(set1,set2))
    jac <- ovr/uns
}


msigdb_ids <- names(pathways2)

jacc_mat <- matrix(0,length(msigdb_ids),length(msigdb_ids))
rownames(jacc_mat) <- colnames(jacc_mat) <- (msigdb_ids)

jacTable <- data.frame(Pathway.A = NA,Pathway.B = NA, Jaccard.Ind = NA )
for(i in seq_along(msigdb_ids)){
    print(paste0("scanning ", msigdb_ids[i]))
    
    temp <- foreach(j = seq_along(msigdb_ids),
                    .combine = "rbind") %dopar%
        {
            
            temp_ind      <-  jaccard.ind(pathways2[[i]],pathways2[[j]]) 
            #jacc_mat[i,j] <- temp_ind
            
            # jacTable <- rbind(jacTable, c(names(pathways2)[i],
            #                               names(pathways2)[j], temp_ind ))
            # 
            return(c(names(pathways2)[i],names(pathways2)[j], temp_ind ))
        }
    
   temp  <- as.data.frame(temp)
   names(temp) <- names(jacTable)
   rownames(temp) <- NULL
   jacTable <- rbind(jacTable, temp)
}



jacTable <- jacTable[!is.na(jacTable$Pathway.A),]



library(reshape2)
jacc_mat <-acast(jacTable, Pathway.A~Pathway.B, value.var="Jaccard.Ind")

jacTable <- jacTable[!(jacTable$Pathway.A == jacTable$Pathway.B),]


write.csv(jacTable,file = "preprocessed/MSigDB_jaccard.csv")
saveRDS(list(jacc_mat,jacTable),"preprocessed/MSigDB_jaccard.RDS")



diag(jacc_mat) <- 0

temp2 <- (jacc_mat)
temp2[lower.tri(temp2)] <- 0
temp2 <- as.data.frame(which(temp2 > .5 ,arr.ind = T))

temp2$row_path <- rownames(jacc_mat)[temp2$row]
temp2$col_path <- colnames(jacc_mat)[temp2$col]

temp2$Jaccard <- 0
for(i in 1:(nrow(temp2))){
    temp2[i,]$Jaccard <- jacc_mat[temp2$row[i],temp2$col[i]]
}


sel_paths <- setdiff((msigdb_ids),unique(temp2$col_path))



pathways2 <- readRDS('preprocessed/MSigDBPathGeneTab.RDS')


msigdb_ids <- pathways2[pathways2$Pathway %in% sel_paths,]



saveRDS(msigdb_ids,"preprocessed/MSigDBPathGeneTabLite.RDS")

