rm(list = ls())



library(tidyverse)
all.kat.files <- list.files("Output", full.names = T,recursive = T,
                            pattern = "diff")


base_file <- tail(all.kat.files,1)


# limit Rosmap to the common pathways
keep.paths <- read_csv(head(all.kat.files,1))
colnames(keep.paths)[1] <- "pathway"




# rank correlation sorted by t statistic 
res_df <- data.frame(pathway = NA, t = NA, adj.P.Val = NA, file = NA)
res_df <- res_df[-1,]
for (j in all.kat.files){
    pour.paths <- read_csv(j)
    colnames(pour.paths)[1] <- "pathway"
    pour.paths <- pour.paths %>% filter(pathway %in% keep.paths$pathway)
    pour.ups      <- pour.paths[,c("pathway","t","adj.P.Val")]
    pour.ups$file <- j
    head(pour.ups)
    res_df <- rbind(res_df, pour.ups)
}

res_df_sq <- res_df %>% dplyr::select(-adj.P.Val) %>% pivot_wider(names_from = file, values_from = t)
cor_sq <- cor(res_df_sq[,-1], method = "spearman", 
              use = "pairwise.complete.obs")

write.csv(cor_sq, file = paste0("Output/pathway_similarity/spearman_cor_AD_model_sq",
                                Sys.Date(), ".csv"))

