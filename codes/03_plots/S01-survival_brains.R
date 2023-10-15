rm(list = ls())
library(stringr)
library(dplyr)
library(ggpubr)






# Selecting all DEP from the I45F/I47F assay

temp        <- list.files("tables/DE_pathways2",recursive = T,
                          full.names = T)

temp        <- grep("ROS|Mayo|MSBB",temp,value = T)


meta.labels <- list()
big.exp     <- list()

fdr.thresh   <- 0.01
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
    t.tab$cohort <- name.tag
    t.tab2   <- t.tab %>% 
                    dplyr::select(.,X,Direction,logFC,P.Value) %>%
                    mutate(.,assay = name.tag) %>%
                    filter(.,Direction != "NONE") %>%
                    mutate(.,path_dir = paste0(Direction,": ", X))
                    
    big.tab.list[[i]] <- t.tab2
    big.tab.list2[[i]] <- t.tab
}


temp <- Reduce(rbind,big.tab.list2)
temp$score <-  -log(temp$P.Value)

ggplot(temp,aes(x = score, y = (1-..y..), color = cohort)) +
  stat_ecdf(geom = "line", size = 1.5,alpha=0.7) + xlim(0, 15) +
  theme_classic()+
  scale_y_continuous(labels = scales::percent) +
  labs(x = "-log(p-value)", y=  "Percentage of Pathways Dysregulated")+
  guides(colour = guide_legend(override.aes = list(size=15),show.legend = "line")) +
  theme(strip.text = element_text(face="bold", size=20),
        plot.title = element_text(size = 30),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title=element_blank(),
        axis.text.y=element_text(size = 15),
        axis.text.x=element_text(size = 15),
        axis.ticks.y=element_blank(),
        legend.position=c(0.75, 0.75), legend.box = "horizontal") + 
  scale_colour_manual(values = c("#FFB000", "#80B1D3", "#4d9221"))

ggsave("figures/S02_pvalue_distribution_pathways.pdf",width = 7, height = 7)


