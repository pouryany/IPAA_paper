rm(list = ls())

# Directory to the pathway gene membership
pathways2 <- readRDS('preprocessed/MSigDBPathGeneTabLite.RDS')

# Use the code below if you want to run reduced overlap results instead
#pathways2 <- readRDS('preprocessed/MSigDBPathGeneTabLite.RDS')


# Directory to the normalized gene expression
gene.address <- list.files("preprocessed/cell_models/",
                           recursive = T,pattern = "D4_vs_G2",full.names = T)
genes.counts2 <- read.delim(gene.address)

genes.counts2 <- genes.counts2 %>% 
    dplyr::select(ID, contains("D4"), contains("G2B2")) %>%
    column_to_rownames("ID")

# keep genes that are present in all data sets
keep_genes <- readRDS("preprocessed/gene_meta_data.RDS")
genes.counts2 <- genes.counts2 %>% filter(rownames(genes.counts2) %in% keep_genes$ID)

# create AD.covariates from names
AD.covariates <- data.frame(Sample = colnames(genes.counts2), Condition = sub("_[0-9]", "", colnames(genes.counts2)))
rownames(AD.covariates) <- AD.covariates$Sample

# Condition for finding DE pathways 
# should it be a column in AD.covariates
condition = 'Condition'

# Where to write the output data, DE pathways
out.dir = 'reduced_pathways_output/D4vsG2B2/'
if(!dir.exists(out.dir))
    dir.create(out.dir)

# These are the covariates that you may need to adjust for. 
# Can include both technical such as RIN or clinical such as Age at Death

# Not necessary if you are not using MSBB data remove this and all other instances
# could it be used as analysis name
tisss <- "D4vsG2B2"


# Making sure factor covariates are treated as such
# Any correction should be on either numeric or factor. No character type allowed
num.cols  <-  lapply(AD.covariates, is.numeric)
fact.cols <-  names(which(num.cols == F))
AD.covariates[,fact.cols] <- lapply(AD.covariates[,fact.cols],
                                    function(X){factor(as.character(X))})


# Making show gene expression and covariates data confirm
genes.counts2 <- genes.counts2[,rownames(AD.covariates)]
de.paths <- PanomiR::differentialPathwayAnalysis(geneCounts = genes.counts2,
                                                 pathways =  pathways2,
                                                 covariates = AD.covariates,
                                                 condition = condition,
                                                 outDir =  paste0(out.dir),
                                                 saveOutName = 
                                                     paste(tisss,"_DEP_AD_vs_Control.RDS"), 
                                                 contrastConds = "ConditionD4-ConditionG2B2")

temp <- de.paths$DEP
temp %>% 
    arrange(desc(logFC)) %>% 
    filter(adj.P.Val < 1) %>%
    write.csv(file = paste0(out.dir, tisss, "_diffPathways", ".csv"))

