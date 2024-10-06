# Clean up the excel sheet for the proteomic data
rm(list=ls())
library(dplyr)
# The choice of Excel file by peptide lengths will not affect the outcomes
proteome <- readxl::read_xlsx("data/TMT_MGHSample_phosphorylation_2024_7_11_90__MasterProteinOFF_6000Items_SortByCoverage.xlsx")

colnames(proteome) <- gsub("-","_",colnames(proteome))
colnames(proteome) <- gsub("#","No",colnames(proteome))

missing_vals <- apply(dplyr::select(proteome, contains("Abundances (Grouped):")),1, FUN= function(x)sum(is.na(x)))

proteome$missing <- missing_vals

proteome_filtered <- proteome |> 
    dplyr::filter(Master == "Master Protein",
                  `Protein FDR Confidence: Combined` == "High",
                  `Coverage [%]` >=2,
                  `Score Sequest HT: Sequest HT` >=1,
                  grepl("Phospho",`Modifications`),
                  missing  < 7) |>
    dplyr::arrange(desc(`Score Sequest HT: Sequest HT`)) 


proteome_filtered <- dplyr::filter(proteome_filtered,!duplicated(proteome_filtered$`Gene Symbol`))




library(clusterProfiler)
library(org.Hs.eg.db)
temp <- clusterProfiler::bitr(proteome_filtered$`Gene Symbol`,
                              fromType = "SYMBOL",toType = "ENSEMBL",OrgDb = org.Hs.eg.db)
temp <- temp[!(duplicated(temp$SYMBOL)),]
temp <- temp[!(duplicated(temp$ENSEMBL)),]

names(temp) <- c("Gene Symbol","ENSEMBL")
proteome_filtered <- dplyr::inner_join(proteome_filtered,temp)



trimmed_expression <- as.data.frame(dplyr::select(proteome_filtered, contains("Abundances (Grouped):")))
trimmed_expression <- as.data.frame(dplyr::select(trimmed_expression, -contains("2G11_Losma")))

row.names(trimmed_expression) <- proteome_filtered$ENSEMBL


saveRDS(trimmed_expression,"data/I45_trimmed_phospho_expression_2024.rds")

row.names(trimmed_expression) <- proteome_filtered$`Gene Symbol`
saveRDS(trimmed_expression,"data/I45_trimmed_phospho_expression_gene_names_2024.rds")
saveRDS(proteome_filtered,"data/phospho_proteome_filtered_2024.rds")

