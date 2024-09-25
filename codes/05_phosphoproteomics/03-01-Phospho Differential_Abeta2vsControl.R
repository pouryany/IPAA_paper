# perform differential expression of miRNAs
rm(list = ls())
library(limma)
library(CovariateAnalysis)
library(tidyr)


genes.counts <- readRDS("data/I45_trimmed_phospho_expression_gene_names_2024.rds")
colnames(genes.counts) <- gsub("Abundances \\(Grouped\\): ","",colnames(genes.counts))



library(preprocessCore)
genes.counts2 <- as.data.frame(normalize.quantiles(as.matrix(log2(genes.counts))))
rownames(genes.counts2) <- rownames(genes.counts)
colnames(genes.counts2) <- colnames(genes.counts)


# 
# prot.polish <- medpolish(as.matrix(log2(genes.counts), eps = 0.01, maxiter = 10,
#                          trace.iter = TRUE,na.rm = F))
# 
# genes.counts2 <- prot.polish$residuals
# rownames(genes.counts2) <- rownames(genes.counts)
# colnames(genes.counts2) <- colnames(genes.counts)

plotDensities(genes.counts2)


AD.covariates <- data.frame(Sample = colnames(genes.counts2), 
                            Condition = sub("_[0-9].*", "", colnames(genes.counts2)))
rownames(AD.covariates) <- AD.covariates$Sample


fact.cols <-  "Condition"
AD.covariates[,fact.cols] <-factor(as.character(AD.covariates[,fact.cols]))

levels(AD.covariates$Condition) <- c("Control","Abeta42","Losma")




DESIGN         <- getDesignMatrix(AD.covariates[,c("Condition"),
                                       drop = F], Intercept = F)



DESIGN$design  <- DESIGN$design[,linColumnFinder(DESIGN$design)$indepCols]
FIT            <- lmFit(genes.counts2, DESIGN$design)

contrast  <- makeContrasts(contrasts=c("ConditionAbeta42-(ConditionControl)"),
                           levels = colnames(FIT$coefficients))
FIT.CONTR <- contrasts.fit(FIT, contrasts=contrast)
FIT.CONTR <- eBayes(FIT.CONTR)
tT        <- topTable(FIT.CONTR, adjust="fdr", sort.by="p", number=Inf)


library(EnhancedVolcano)
library(tibble)
tT <- rownames_to_column(tT, var = "Gene Symbol")
# Here is the differential expression table


sel_labs <- (tT[1:50, ]$`Gene Symbol`)

sel_labs <- sel_labs[sel_labs!= ""]
keyvals <- ifelse(
        (tT$adj.P.Val > 0.01 | abs(tT$logFC) < 0.5), 'grey60',
    ifelse((tT$logFC) < 0, 'royalblue',
           'red3'))

names(keyvals)[keyvals == 'red3'] <- 'Up in A\u03B242-H'
names(keyvals)[keyvals == 'grey60'] <- ''
names(keyvals)[keyvals == 'royalblue'] <- 'Down in A\u03B242-H'


p1 <-   EnhancedVolcano(tT,
                        lab = tT$`Gene Symbol`, ylim = c(0,6),
                        selectLab = sel_labs,
                        x = "logFC",
                        y = "adj.P.Val",
                        caption = NULL,
                        xlab = bquote(~Log[2]~ "fold change"),
                        ylab = bquote(~-Log[10]~"Adj."~italic(p)~"value"),
                        title = NULL,
                        subtitle = "",
                        pCutoff = 0.01,
                        pointSize = 2.0,
                        labSize = 6,
                        labCol = 'black',
                        labFace = 'plain',
                        drawConnectors = T, 
                        arrowheads = F,
                        boxedLabels = FALSE,
                        widthConnectors = 1.2,
                        colConnectors = 'black',
                        legendPosition = "top",
                        colCustom = keyvals,
                        axisLabSize = 30,
                        legendLabSize = 20,
                        legendIconSize = 10,
                        FCcutoff = (0.5)
)

ggsave("figures/Volcano_phospho_proteome_Abeta42_vs_control.pdf",p1,
       width = 12, height = 9)




tT   <- tT %>% 
    arrange(.,P.Value) %>%
    mutate(Direction = logFC/abs(logFC),
           Direction = factor(Direction, c(-1,1),
                              c('-1' = 'DOWN', '1' = 'UP')),
           Direction = as.character(Direction))

tT$Direction[tT$adj.P.Val > 0.01 | abs(tT$logFC) < 0.5] = 'NONE'
table(tT$Direction)

full_data <- readRDS("data/phospho_proteome_filtered_2024.rds")

tT <- inner_join(tT,full_data)

write.csv(tT,"tables/DE_phospho/phospho_differential_Abeta42_control.csv")
