rm(list = ls())
library(PanomiR)
#
# Directory to the pathway gene membership
pathways2 <- readRDS('preprocessed/MSigDBPathGeneTab.RDS')

# Use the code below if you want to run reduced overlap results instead
#pathways2 <- readRDS('preprocessed/MSigDBPathGeneTabLite.RDS')

# Directory to the normalized gene expression
# genes.counts2 <- readRDS('../../MSBB_data/preprocessed/NEWCOUNTS.RDS')


genes.counts2 <- readRDS("data/I45_trimmed_phospho_expression_2024.rds")


colnames(genes.counts2) <- gsub("Abundances \\(Grouped\\): ","",colnames(genes.counts2))
# genes.counts2 <- genes.counts2 %>%
#     dplyr::select(ID, contains("I45F"), contains("I47F")) %>%
#     column_to_rownames("ID")

# keep genes that are present in all data sets

# Directory to the normalized gene expression
# create AD.covariates from names
AD.covariates <- data.frame(Sample = colnames(genes.counts2), 
                            Condition = sub("_[0-9].*", "", colnames(genes.counts2)))
rownames(AD.covariates) <- AD.covariates$Sample

# Condition for finding DE pathways 
# should it be a column in AD.covariates
condition = 'Condition'

# Where to write the output data, DE pathways
out.dir = 'output_phospho/Phospho2024/'
if(!dir.exists(out.dir))
    dir.create(out.dir,recursive = T)

# These are the covariates that you may need to adjust for. 
# Can include both technical such as RIN or clinical such as Age at Death
# adj.covars <- readRDS("../../MSBB_data/preprocessed/POSTCOVAR.RDS")

# Not necessary if you are not using MSBB data remove this and all other instances
# could it be used as analysis name
tisss <- "PhosphoI45F_vs_Control"
#AD.covariates   <- AD.covariates[AD.covariates$Diagnosis %in% c("AD","CONTROL"),]
#AD.covariates   <- AD.covariates[AD.covariates$Tissue %in% tisss,]



# Making sure factor covariates are treated as such
# Any correction should be on either numeric or factor. No character type allowed
num.cols  <-  lapply(AD.covariates, is.numeric)
fact.cols <-  "Condition"
AD.covariates[,fact.cols] <-factor(as.character(AD.covariates[,fact.cols]))

levels(AD.covariates$Condition) <- c("Control","Abeta42","Losma")

# Making show gene expression and covariates data confirm
genes.counts2 <- genes.counts2[,rownames(AD.covariates)]

# In the following:
#   1- Gene counts should be corrected for length and GC content preferably
#   2- Modify outdir argument according to your own project
#   3- 
de.paths <- PanomiR::differentialPathwayAnalysis(geneCounts = genes.counts2,
                                                 pathways =  pathways2,
                                                 covariates = AD.covariates,
                                                 condition = condition,
                                                 outDir =  paste0(out.dir),
                                                 id = "ENSEMBL", trim = 0,
                                                 minPathSize = 5,
                                                 saveOutName = 
                                                     paste(tisss,"_DEP_AD_vs_Control.RDS"), 
                                                 contrastConds = "ConditionAbeta42-ConditionControl")

# contrast.conds = "DiagnosisAD-DiagnosisCONTROL")
temp <- de.paths$DEP
temp %>% 
    arrange(desc(logFC)) %>% 
    filter(adj.P.Val < 1) %>%
    write.csv(file = paste0(out.dir, tisss, "_diffPathways", ".csv"))


tisss <- "PhosphoLosmapimod_vs_I45F"

de.paths <- PanomiR::differentialPathwayAnalysis(geneCounts = genes.counts2,
                                                 pathways =  pathways2,
                                                 covariates = AD.covariates,
                                                 condition = condition,
                                                 outDir =  paste0(out.dir),
                                                 id = "ENSEMBL", trim = 0,
                                                 minPathSize = 5,
                                                 saveOutName = 
                                                     paste(tisss,"_DEP_AD_vs_Control.RDS"), 
                                                 contrastConds = "ConditionLosma-ConditionAbeta42")

# contrast.conds = "DiagnosisAD-DiagnosisCONTROL")
temp <- de.paths$DEP
temp %>% 
    arrange(desc(logFC)) %>% 
    filter(adj.P.Val < 1) %>%
    write.csv(file = paste0(out.dir, tisss, "_diffPathways", ".csv"))

