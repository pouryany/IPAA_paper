rm(list = ls())

# Default directory assumed to be in place of the codes
library(CovariateAnalysis)
library(edgeR)
library(PanomiR)
library(dplyr)
library(ggplot2)



# Directory to the pathway gene membership

# Directory to the normalized gene expression
# genes.counts2 <- readRDS('../../MSBB_data/preprocessed/NEWCOUNTS.RDS')
genes.counts2 <- read.csv("codes/external_dataset/GSE128343_RNAseq_raw.txt",sep = "\t")
# genes.counts2 <- genes.counts2$

CQN_correction = T
if(CQN_correction){
    
    temp.exp <- genes.counts2
    rownames(temp.exp) <- temp.exp$GeneID
    
    # 
    # library(EDASeq)
    # temp.length <- getGeneLengthAndGCContent(id =rownames(temp.exp),  org=c( "hsa"))
    temp.exp <- temp.exp[,-c(1)]
    
    processed.counts = getGeneFilteredGeneExprMatrix(temp.exp,
                                                     MIN_GENE_CPM=1, 
                                                     MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM=0.2,
                                                     verbose = T)
    
    
    
    
    
    GENE.PARAM  <- read.csv("preprocessed/geneParameters.tsv",sep = "\t")
    GENE.LEN <-  dplyr::select(GENE.PARAM, ensembl_gene_id, gene.length) %>% 
        unique() 
    rownames(GENE.LEN) = GENE.LEN$ensembl_gene_id
    
    GENE.GC  <- dplyr::select(GENE.PARAM, ensembl_gene_id, percentage_gc_content) %>% 
        unique() 
    rownames(GENE.GC) = GENE.GC$ensembl_gene_id 
    
    
    genesToAnalyze = processed.counts$filteredExprMatrix$genes$genes
    
    genesToAnalyze = unlist(genesToAnalyze) %>% 
        unique() %>%
        intersect(GENE.GC$ensembl_gene_id[!is.na(GENE.GC$percentage_gc_content)]) %>%
        intersect(GENE.LEN$ensembl_gene_id[!is.na(GENE.LEN$gene.length)]) %>%
        setdiff(c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"))
    
    
    PROCESSED_COUNTS = getGeneFilteredGeneExprMatrix(temp.exp[genesToAnalyze, ], 
                                                     MIN_GENE_CPM=0, 
                                                     MIN_SAMPLE_PERCENT_WITH_MIN_GENE_CPM=0.0,
                                                     verbose = T)
    
    
    library(cqn)
    
    CQN.GENE_EXPRESSION = cqn(PROCESSED_COUNTS$filteredExprMatrix$counts, 
                              x = GENE.GC[PROCESSED_COUNTS$filteredExprMatrix$genes$genes, 'percentage_gc_content'],
                              lengths = GENE.LEN[PROCESSED_COUNTS$filteredExprMatrix$genes$genes, 'gene.length'],
                              lengthMethod = "smooth", 
                              verbose = FALSE)
    
    # CQN Normalization
    
    CQN.GENE_EXPRESSION$E = CQN.GENE_EXPRESSION$y + CQN.GENE_EXPRESSION$offset
    
    exp.data = (2^CQN.GENE_EXPRESSION$E) 
    
    
}


colnames(exp.data)
AD.covariates <- data.frame(ID = colnames(exp.data))
AD.covariates$Type <- gsub("RNAseq.","",AD.covariates$ID)
AD.covariates$Type <- gsub("[1-9]|_","",AD.covariates$Type)
AD.covariates$Type <- gsub("\\.","_",AD.covariates$Type)

# Directory to the normalized gene expression
# create AD.covariates from names

# Condition for finding DE pathways 
# should it be a column in AD.covariates
condition = 'Type'

# Where to write the output data, DE pathways
out.dir = 'codes/external_dataset/'

# These are the covariates that you may need to adjust for. 
# Can include both technical such as RIN or clinical such as Age at Death

# Not necessary if you are not using MSBB data remove this and all other instances
# could it be used as analysis name
tisss <- "GSE_dataset"
#AD.covariates   <- AD.covariates[AD.covariates$Diagnosis %in% c("AD","CONTROL"),]
#AD.covariates   <- AD.covariates[AD.covariates$Tissue %in% tisss,]



# Making sure factor covariates are treated as such
# Any correction should be on either numeric or factor. No character type allowed
num.cols  <-  lapply(AD.covariates, is.numeric)
fact.cols <-  names(which(num.cols == F))
AD.covariates[,fact.cols] <- lapply(AD.covariates[,fact.cols],
                                    function(X){factor(as.character(X))})


# Making show gene expression and covariates data confirm
genes.counts2 <- exp.data[,(AD.covariates$ID)]

# In the following:
#   1- Gene counts should be corrected for length and GC content preferably
#   2- Modify outdir argument according to your own project
#   3- 

pathways2 <- readRDS("preprocessed/MSigDBPathGeneTab.RDS")



de.paths <- PanomiR::differentialPathwayAnalysis(geneCounts = genes.counts2,
                                     pathways = pathways2,
                                     covariates = AD.covariates,
                                     condition = condition,
                                     outDir = out.dir,
                                     saveOutName = paste(tisss,"_DEP_APP_vs_wt.RDS"),
                                     contrastConds = "TypeAPP-Typewt")


temp <- de.paths$DEP
temp %>% 
    dplyr::arrange((P.Value)) %>% 
    write.csv(file = paste0(out.dir,
                            tisss, "APP_vs_wt_diffPathways", Sys.Date(), ".csv"),
              row.names = T)



de.paths <- PanomiR::differentialPathwayAnalysis(geneCounts = genes.counts2,
                                                 pathways = pathways2,
                                                 covariates = AD.covariates,
                                                 condition = condition,
                                                 outDir = out.dir,
                                                 saveOutName = paste(tisss,"_DEP_PSEN_vs_wt.RDS"),
                                                 contrastConds = "TypePSEN-Typewt")


temp <- de.paths$DEP
temp %>% 
    arrange((P.Value)) %>% 
    write.csv(file = paste0(out.dir,
                            tisss,
                            "PSEN_vs_wt_diffPathways",
                            Sys.Date(), ".csv"),
              row.names = T)


