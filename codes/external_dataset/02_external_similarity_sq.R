rm(list = ls())



library(tidyverse)
all.kat.files <- list.files("codes/external_dataset", full.names = T, pattern = "diff")

all.kat.files2 <- list.files("Output", full.names = T,
                             pattern = "(ROSMAP|Mayo|MSBB)",recursive = T )
all.kat.files2 <- grep("diff",all.kat.files2,value = T)
all.kat.files <- union(all.kat.files, all.kat.files2)

pour.p.val <- 0.01
kat.p.val <- 0.01
n_paths <- 2000
base_file <- "ROSMAP_DLPFC_ADvsControl_diffPathways2023-02-14.csv"

res_df <- data.frame(path = rep(all.kat.files, each = length(all.kat.files)), file = NA, n_paths = NA, pval = NA, n_overlap = NA, up_up = NA, up_dn = NA, dn_dn = NA, dn_up = NA, n_up = NA, n_dn = NA, base_file = rep(all.kat.files, times = length(all.kat.files)), base_pval = pour.p.val, kat_pval = kat.p.val)

# limit Rosmap to the common pathways
keep.paths <- read_csv("Output/A5vsG2B2/A5vsG2B2_diffPathways2023-02-14.csv")
colnames(keep.paths)[1] <- "pathway"
static_paths  <- grep("Static_Module",keep.paths$pathway,value = T)
keep.paths     <- keep.paths[!(keep.paths$pathway %in% static_paths),]



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
cor_sq <- cor(res_df_sq[,-1], method = "spearman", use = "pairwise.complete.obs")

# cor(res_df_sq[,2], res_df_sq[,3], method = "spearman", use = "pairwise.complete.obs")
write.csv(cor_sq, file = paste0("Output/similarity_external/spearman_cor_AD_model_sq", Sys.Date(), ".csv"))
