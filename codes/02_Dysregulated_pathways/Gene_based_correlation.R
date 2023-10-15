rm(list = ls())
library(tidyverse)



DEG_directory <- "tables/DE_genes_filtered_Mjd"
all.kat.files <- list.files(DEG_directory, full.names = T, recursive = T,
                            pattern = "\\.csv|\\.txt")

keep_genes <- readRDS("preprocessed/gene_meta_data.RDS")


pour.p.val <- 0.05
kat.p.val <- 0.05
n_paths <- 10000
base_file <- "ROSMAP_AD_vs_Control.csv"
# base_file <- "Mayo_ADvsControl_QC_diffPathways2020-11-20.csv"


# rank correlation sorted by t statistic 
res_df <- data.frame(ID = NA, logFC = NA, FDR = NA, file.de = NA)
res_df <- res_df[-1,]
for (j in all.kat.files){
    if (endsWith(j, "csv")){
        pour.paths <- read_csv(j)
        if(any(colnames(pour.paths) %in% c("X1")))
            pour.paths$X1 <- NULL
    } else if (endsWith(j, "txt")) {
        pour.paths <- read.delim(j)
    }
    
   
    pour.paths <- pour.paths %>% filter(ID %in% keep_genes$ID)

    pour.ups      <- pour.paths[,c("ID","logFC", "PValue")]
    pour.ups$file.de <- j
    head(pour.ups)
    res_df <- rbind(res_df, pour.ups)
}
res_df    <- res_df %>% distinct()
res_df2    <- res_df  %>% mutate(.,t_like = sign(logFC) * (-log10(PValue)))
res_df_sq <- res_df2 %>% dplyr::select(-PValue, -logFC) %>%
  pivot_wider(names_from = file.de, values_from = t_like)



cor_sq <- cor(res_df_sq[,-1], method = "spearman", use = "pairwise.complete.obs")

write.csv(cor_sq, file = paste0("Output/gene_similarity/",
                                "spearman_cor_AD_model_sq_gene_",
                                Sys.Date(), ".csv"))

