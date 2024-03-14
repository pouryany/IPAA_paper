rm(list = ls())

#
# Directory to the pathway gene membership
pathways2 <- readRDS('preprocessed/MSigDBPathGeneTab.RDS')

# Use the code below if you want to run reduced overlap results instead
#pathways2 <- readRDS('preprocessed/MSigDBPathGeneTabLite.RDS')

# Directory to the normalized gene expression
# genes.counts2 <- readRDS('../../MSBB_data/preprocessed/NEWCOUNTS.RDS')


gene.address <- list.files("preprocessed/cell_models/",
                           recursive = T,pattern = "I45F_vs_I47",full.names = T)
genes.counts2 <- read.delim(gene.address)

genes.counts2 <- genes.counts2 %>%
    dplyr::select(ID, contains("I45F"), contains("I47F")) %>%
    column_to_rownames("ID")

# keep genes that are present in all data sets
keep_genes <- readRDS("preprocessed/gene_meta_data.RDS")
genes.counts2 <- genes.counts2 %>% filter(rownames(genes.counts2) %in% keep_genes$ID)

# Directory to the normalized gene expression
# create AD.covariates from names
AD.covariates <- data.frame(Sample = colnames(genes.counts2), Condition = sub("_[0-9].*", "", colnames(genes.counts2)))
rownames(AD.covariates) <- AD.covariates$Sample

# Condition for finding DE pathways 
# should it be a column in AD.covariates
condition = 'Condition'

# Where to write the output data, DE pathways
out.dir = 'Output/I45FvsI47F/'
if(!dir.exists(out.dir))
    dir.create(out.dir)

# These are the covariates that you may need to adjust for. 
# Can include both technical such as RIN or clinical such as Age at Death
# adj.covars <- readRDS("../../MSBB_data/preprocessed/POSTCOVAR.RDS")

# Not necessary if you are not using MSBB data remove this and all other instances
# could it be used as analysis name
tisss <- "I45FvsI47F"
#AD.covariates   <- AD.covariates[AD.covariates$Diagnosis %in% c("AD","CONTROL"),]
#AD.covariates   <- AD.covariates[AD.covariates$Tissue %in% tisss,]



# Making sure factor covariates are treated as such
# Any correction should be on either numeric or factor. No character type allowed
num.cols  <-  lapply(AD.covariates, is.numeric)
fact.cols <-  names(which(num.cols == F))
AD.covariates[,fact.cols] <- lapply(AD.covariates[,fact.cols],
                                    function(X){factor(as.character(X))})


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
                                                 saveOutName = 
                                                     paste(tisss,"_DEP_AD_vs_Control.RDS"), 
                                                 contrastConds = "ConditionI45F-ConditionI47F")

                                        # contrast.conds = "DiagnosisAD-DiagnosisCONTROL")
temp <- de.paths$DEP
temp %>% 
    arrange(desc(logFC)) %>% 
    filter(adj.P.Val < 1) %>%
    write.csv(file = paste0(out.dir, tisss, "_diffPathways", ".csv"))

