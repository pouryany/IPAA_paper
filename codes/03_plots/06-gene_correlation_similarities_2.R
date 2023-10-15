rm(list= ls())
library(stringr)
library(dplyr)
library(ggpubr)






# Selecting all DEP from the I45F/I47F assay

temp        <- list.files("tables/DE_genes_clean_2/",recursive = T,
                          full.names = T)
tab.rosmap  <- grep("Mayo_AD_|ROSMAP|MSBB", temp,value = T)
temp        <- c(tab.rosmap)
keep_genes <- readRDS("preprocessed/gene_meta_data.RDS")

meta.labels <- list()
big.exp     <- list()

fdr.thresh   <- 0.05
logFC.thresh <- 0

plot.list    <- list()
big.tab.list <- list()
big.tab.list2 <- list()

for(i in temp) {
    base.ind <- which(temp == i)
    name.tag <- unlist(str_split(i,"/"))
    name.tag <- tail(name.tag,1)
    name.tag <- unlist(str_split(name.tag,pattern = "_"))
    name.tag <- name.tag[1]
    t.tab    <- read.csv(i)
    t.tab    <- t.tab[t.tab$ID %in% as.character(keep_genes$ID),]
    t.tab2   <- t.tab %>% 
                    dplyr::select(.,Name,ID,Direction) %>%
                    mutate(.,assay = name.tag) %>%
                    filter(.,Direction != "NONE")
                    
    big.tab.list[[i]] <- t.tab2
    big.tab.list2[[i]] <- t.tab
    
}



  temp <- merge(big.tab.list2[[1]],
              big.tab.list2[[2]],
              by = "ID")




temp$dirx <- plyr::mapvalues(temp$Direction.x,c("UP","DOWN","NONE"), c(1,-1,0))
temp$diry <- plyr::mapvalues(temp$Direction.y,c("UP","DOWN","NONE"), c(1,-1,0))

temp$Direction.both <- as.factor(as.numeric(temp$dirx) * as.numeric(temp$diry))

temp$Direction.both <- plyr::mapvalues(temp$Direction.both, c(0,1), c("NONE","Significant in both\nMayo and MSBB"))
# 
# temp    <- temp  %>% mutate(.,t.x = sign(logFC.x) * (-log10(PValue.x)))
# temp    <- temp  %>% mutate(.,t.y = sign(logFC.y) * (-log10(PValue.y)))
# 
#  
# 
# temp    <- temp  %>% mutate(.,t.x = ifelse(sign(logFC.x) > 0,  qnorm(1-PValue.x/2,lower.tail = TRUE),
#                                            qnorm(1-PValue.x/2,lower.tail = FALSE)))
# temp    <- temp  %>% mutate(.,t.y = ifelse(sign(logFC.y) > 0,  qnorm(1-PValue.y/2,lower.tail = TRUE),
#                                            qnorm(1-PValue.y/2,lower.tail = FALSE)))
# 
# 
# 
# 
# 
# ifelse(sign(logFC.x) > 0,  qnorm(1-PValue.x,lower.tail = TRUE),
#        qnorm(1-PValue.x,lower.tail = FALSE))
# 
# 
# qnorm(0.99,lower.tail = F)
 corrs <- signif(cor(temp$t.x,temp$t.y,method = "spearman"),2)
# 
# ?pnorm

ggplot(data = temp, aes(t.x, t.y)) +
  geom_point(alpha = 0.5, size = 3) + 
  theme_classic() +
  labs (x = " Mayo AD vs control (t-score)",
        y = "MSBB AD vs control (t-score)",
        title = "Gene dysregulation" ) +
  scale_color_manual(values = c("black", "#b2182b", "grey")) +
  theme(plot.title = element_text(size = 32 , face = "bold", hjust = 0.5),
        axis.title = element_text(size = 20,
                                  angle = 0,
                                  vjust = 0.5, hjust = 0.5,
                                  face = "plain"),
        axis.text = element_text(size = 16),
        legend.text = element_text(face = "plain"),
        legend.position = "none") +
  annotate("text",x=  -4, y = 5.5,
           label = paste0("Spearman's rho = ", corrs), size = 8) 

ggsave("figures/02_f_correlation_genes.pdf",width = 8, height = 8)




