rm(list = ls())
library("sva")
library("limma")
library("edgeR")
library("cqn")
library("dplyr")
library(CovariateAnalysis)
library(ggfortify)




exp.data    <- readRDS(file = "preprocessed/collapsed_samples_A5.RDS")
meta.info   <- readRDS(file = "preprocessed/collapsed_meta_data_A5.RDS")

rownames(exp.data) <- as.character(exp.data$ID)

exp.data <- exp.data[,as.character(meta.info$ID)]
lapply(exp.data, class)


meta.info[,2] <- as.character(meta.info[,2])
meta.info[1:3,2] <- paste0(meta.info[1:3,2],"1")
meta.info$batch  <- "Batch1"

colnames(exp.data) <- meta.info$ID




exp.data2    <- readRDS(file = "preprocessed/collapsed_samples_I45F.RDS")
meta.info2   <- readRDS(file = "preprocessed/collapsed_meta_data_I45F.RDS")
rownames(exp.data2) <- as.character(exp.data2$ID)

exp.data2 <- exp.data2[,as.character(meta.info2$ID)]

meta.info2[,2] <- as.character(meta.info2[,2])
meta.info2$batch  <- "Batch2"

colnames(exp.data2) <- meta.info2$ID



shared_exp <- intersect(rownames(exp.data), rownames(exp.data2))

shared_exp <- cbind(exp.data[shared_exp,], exp.data2[shared_exp,])
shared_meta <- rbind(meta.info, meta.info2)
rownames(shared_meta) <- shared_meta$ID

colnames(shared_exp)

shared_exp[is.na(shared_exp)] <- 0

 mod0 = model.matrix(~1,data = shared_meta)
 mod  = model.matrix(~(Type),data = shared_meta)
 
 modcombat = model.matrix(~1, data=shared_meta)
 combat_edata = ComBat(dat=(shared_exp), batch=(shared_meta$batch),
                       mod=mod,
                       par.prior=TRUE, prior.plots=FALSE)


  
  totSignals   <- t(as.matrix(combat_edata))
  totType      <- (shared_meta$Type)
  
  
  Tissue.Type <- totType
  class(totSignals) <- "numeric"
  tot.pca <- prcomp(totSignals, center = T, scale. = T)
  

  
  
  p1 <- autoplot(tot.pca, data = shared_meta,x = 1, y = 2,colour = 'Type', size = 8,alpha = 0.8) +
    theme_bw(base_size = 30)+
    scale_color_manual(values=c("#7fc97f","#beaed4","#fdc086",
                                "#386cb0", "#f0027f", "#bf5b17"))+
    theme(axis.text=element_blank(),
          axis.ticks = element_blank())
  
  ggsave(filename = "figures/S02_All_I45F_pca.pdf", p1, width = 6, height = 4)
  

