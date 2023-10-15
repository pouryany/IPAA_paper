rm(list = ls())

library(PanomiR)
library(dplyr)


# Directory to the pathway gene membership
pathways2 <- readRDS('preprocessed/MSigDBPathGeneTab.RDS')

# Directory to the normalized gene expression
genes.counts2 <- readRDS("preprocessed/MAYO/NEWCOUNTS.RDS")
genes.counts2 <- as.data.frame(genes.counts2)


# Directory to the normalized gene expression
# TO DO: create AD.covariates from names
AD.covariates <- readRDS("preprocessed/MAYO/COVARIATES.RDS")
# Condition for finding DE pathways 
condition = 'Tissue.Diagnosis'

# Where to write the output data, DE pathways
out.dir0 = 'Output/MAYO/'

if(!dir.exists(out.dir0))
    dir.create(out.dir0)
# These are the covariates that you may need to adjust for. 
# Can include both technical such as RIN or clinical such as Age at Death
adj.covars <- readRDS("preprocessed/MAYO/POSTCOVAR.RDS")
tisss      <- "Mayo_ADvsControl"

# Making sure factor covariates are treated as such
# Any correction should be on either numeric or factor. No character type allowed

# Making show gene expression and covariates data confirm
genes.counts2 <- genes.counts2[,rownames(AD.covariates)]



de.paths <- PanomiR::differentialPathwayAnalysis(geneCounts = genes.counts2,
                                     pathways = pathways2,
                                     covariates = AD.covariates,
                                     condition = condition,
                                     outDir  = out.dir0,
                                     saveOutName =  paste(tisss,"_DEP_AD_vs_CONTROL.RDS"),
                                     adjustCovars =adj.covars,
                                     contrastConds = "Tissue.DiagnosisTCX.AD-Tissue.DiagnosisTCX.CONTROL")

temp <- de.paths$DEP
temp <- temp[order(temp$P.Value),]
temp %>% 
    write.csv(file = paste0(out.dir0, tisss, "_diffPathways", Sys.Date(), ".csv"),row.names = T)

temp <- tibble::rownames_to_column(temp)


de.paths$PathwayResiduals |> saveRDS(
                                     paste0(out.dir0,
                                            tisss, "_Residuals",
                                            Sys.Date(),
                                            ".RDS")
                                     )





