rm(list = ls())
library(stringr)
library(dplyr)
library(ggpubr)

# The code to harmonize DEG table representation across all datasets.



clean.names <- c("Name","ID","logFC","logCPM","PValue","FDR","Direction")




mayo.de <- read.csv("preprocessed/DE_genes/DE_genes_TCX.csv")
mayo.de <- mayo.de %>% 
     dplyr::select(., hgnc_symbol,ensembl_gene_id, logFC,AveExpr,P.Value,adj.P.Val,Direction) 
names(mayo.de) <- clean.names
write.csv(mayo.de,"tables/DE_genes_filtered_Mjd/Mayo_AD_vs_Control.csv",
          row.names = F)

mayo.de <- read.csv("preprocessed/DE_genes/DE_genes_MSBB_AD_vs_Control.csv")

mayo.de <- mayo.de %>% 
    dplyr::select(., hgnc_symbol,ensembl_gene_id,
                  logFC,AveExpr,P.Value,adj.P.Val,Direction) 
names(mayo.de) <- clean.names

write.csv(mayo.de,"tables/DE_genes_filtered_Mjd/MSBB_AD_vs_Control.csv",
          row.names = F)


ros.de <- read.csv("preprocessed/DE_genes/DE_genes_ROSMAP.csv")
ros.de <- ros.de %>% 
    dplyr::select(., hgnc_symbol,ensembl_gene_id,
                  logFC,AveExpr,P.Value,adj.P.Val,Direction) 
names(ros.de) <- clean.names
write.csv(ros.de,"tables/DE_genes_filtered_Mjd/ROSMAP_AD_vs_Control.csv",row.names = F)



#### With t-scores

mayo.de <- read.csv("preprocessed/DE_genes/DE_genes_TCX.csv")
mayo.de <- mayo.de %>% 
    dplyr::select(., hgnc_symbol,ensembl_gene_id,
                  logFC,AveExpr,P.Value,adj.P.Val,Direction,t) 
clean.names <- c("Name","ID","logFC","logCPM",
                 "PValue","FDR","Direction","t")
names(mayo.de) <- clean.names

write.csv(mayo.de,"tables/DE_genes_clean_2/Mayo_AD_vs_Control.csv",
          row.names = F)




ros.de <- read.csv("preprocessed/DE_genes/DE_genes_ROSMAP.csv")
ros.de <- ros.de %>% 
    dplyr::select(., hgnc_symbol,ensembl_gene_id,
                  logFC,AveExpr,P.Value,adj.P.Val,Direction,t)
names(ros.de) <- clean.names
write.csv(ros.de,"tables/DE_genes_clean_2/ROSMAP_AD_vs_Control.csv",
          row.names = F)




msbb.de <- read.csv("preprocessed/DE_genes/DE_genes_MSBB_AD_vs_Control.csv")
msbb.de <- msbb.de %>% 
    dplyr::select(., hgnc_symbol,ensembl_gene_id,
                  logFC,AveExpr,P.Value,adj.P.Val,Direction,t) 

names(msbb.de) <- clean.names
write.csv(msbb.de,"tables/DE_genes_clean_2/MSBB_AD_vs_Control.csv",row.names = F)



