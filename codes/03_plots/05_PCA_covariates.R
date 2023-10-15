rm(list = ls())
library(CovariateAnalysis)


library(data.table)
library(plyr)
library(tidyverse)
library(psych)
library(limma)


msbb_dat      <- readRDS("Output/MSBB/MSBB_PHG_ADvsControl _DEP_AD_vs_Control.RDS")
AD.covariates <- readRDS("preprocessed/MSBB/COVARIATES.RDS")
adj.covars    <- readRDS("preprocessed/MSBB/POSTCOVAR.RDS")


path_express <- msbb_dat$pathwaySummaryStats


preAdjustedSigCovars = runPCAandPlotCorrelations(path_express, 
                                                 AD.covariates[,c("CDR", "Diagnosis",adj.covars)],
                                                 'uncorrected pathway activity in MSBB', 
                                                 isKeyPlot=TRUE, 
                                                 CORRELATION_TYPE = "spearman",
                                                 MIN_PVE_PCT_PC = 1)

path_express <- msbb_dat$PathwayResiduals


postAdjustedSigCovars = runPCAandPlotCorrelations(path_express, 
                                                 AD.covariates[,c("CDR","Diagnosis",adj.covars)],
                                                 'corrected pathway activity in MSBB', 
                                                 CORRELATION_TYPE = "spearman",
                                                 isKeyPlot=TRUE, 
                                                 MIN_PVE_PCT_PC = 1)

library(patchwork)


p1 <- preAdjustedSigCovars$PC_res[[2]]$plotData
p2 <- postAdjustedSigCovars$PC_res[[2]]$plotData


pdf("figures/S01_msbb_pathway_adjustment.pdf",
    width = 8,height = 12)
p1 / p2
dev.off()




MAYO_dat      <- readRDS("Output/MAYO/Mayo_ADvsControl _DEP_AD_vs_CONTROL.RDS")
AD.covariates <- readRDS("preprocessed/MAYO/COVARIATES.RDS")
adj.covars <- readRDS("preprocessed/MAYO/POSTCOVAR.RDS")


path_express <- MAYO_dat$pathwaySummaryStats


preAdjustedSigCovars = runPCAandPlotCorrelations(path_express, 
                                                 AD.covariates[,c("Tissue.Diagnosis",adj.covars)],
                                                 'uncorrected pathway activity in Mayo', 
                                                 isKeyPlot=TRUE, 
                                                 CORRELATION_TYPE = "spearman",
                                                 MIN_PVE_PCT_PC = 1)

path_express <- MAYO_dat$PathwayResiduals


postAdjustedSigCovars = runPCAandPlotCorrelations(path_express, 
                                                  AD.covariates[,c("Tissue.Diagnosis",adj.covars)],
                                                  'corrected pathway activity in Mayo', 
                                                  CORRELATION_TYPE = "spearman",
                                                  isKeyPlot=TRUE, 
                                                  MIN_PVE_PCT_PC = 1)

library(patchwork)


p1 <- preAdjustedSigCovars$PC_res[[2]]$plotData
p2 <- postAdjustedSigCovars$PC_res[[2]]$plotData


pdf("figures/S01_Mayo_pathway_adjustment.pdf",
    width = 8,height = 12)
p1 / p2
dev.off()
