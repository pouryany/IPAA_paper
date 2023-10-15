rm(list = ls())
library(stringr)
library(dplyr)
library(ggpubr)





fix_header_bigTab <- function(big.tab.list){
  
  big_names <- names(big.tab.list)
  
  big_names <- gsub(pattern = ".*/|_QC.*|_diffPathways.*|2021.*","", big_names)
  
  
  big_names<- gsub("vsG2B2","",big_names)
  big_names <- gsub("I45FvsI47F","A\u03B242-H / A\u03B240-H",big_names)
  big_names <- gsub("I45F","A\u03B242-H",big_names)
  big_names <- gsub("I47F","A\u03B240-H",big_names)
  big_names <- gsub("v.*","",big_names)
  big_names <- gsub("_"," ",big_names)
  big_names <- gsub(" PHG","",big_names)
  big_names <- gsub(" DLPFC","",big_names)
  
  big_names <- paste0(big_names, " vs Ctrl")
  big_ind   <- grep("A\u03B242-H / A\u03B240-H",big_names)
  big_names[big_ind] <- "A\u03B242-H vs A\u03B240-H"
  return(big_names)
}


plot_agreement <- function(temp, name_tag, legend_pos = "none"){
  
  nums <- table(temp$Direction.x , temp$Direction.y)
  chi_test <- chisq.test(nums)
  chi_test <- signif(chi_test$p.value,3)
  agreement    <- diag(nums)
  disagreement <- sum(nums) - sum(diag(nums))
  agreement    <- as.data.frame(t(agreement))
  agreement$disagreement <- -1 * disagreement
  agreement <- t(agreement)
  agreement <- tibble::rownames_to_column(as.data.frame(agreement),"Direction")
  agreement$Direction <- c("Down in AD (Agree)",
                           "Up in AD (Agree)",
                           "Disagree")
  
  p1 <- ggplot(agreement, aes(Direction,V1, fill = Direction)) +
    geom_bar(position= "stack", stat="identity",width = 0.6) +
    theme_minimal() + geom_hline(yintercept = 0, color = "grey30") +
    scale_colour_manual(values = c("#8dd3c7", "#80b1d3", "#fb8072")) +
    scale_fill_manual(values = c("#8073ac", "#4393c3", "#d6604d")) +
    labs(x = name_tag) +
    ylim(c(-90,170)) + 
    geom_text(aes(label = V1, y = V1 + (10 * sign(V1-0.1))), size= 4) + 
    theme(panel.grid.major.x =  element_blank(),
          panel.grid.minor.x =  element_blank(),
          panel.border = element_rect(colour = "grey80", fill=NA, size=2),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 20,
                                      angle = 0,
                                      vjust = 0.5, hjust = 1,
                                      face = "plain"),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          legend.text = element_text(face = "plain"),
          legend.position = legend_pos) +
    annotate("text",x=  1, y = 20, hjust = 0,
             label = paste0("p-value = ", chi_test), size = 4) + 
    coord_flip()
  
  return(p1)
}










# Selecting all DEP from the I45F/I47F assay

temp        <- list.files("tables/DE_pathways2",recursive = T,
                          full.names = T)

temp        <- grep("ROS|Mayo|MSBB|I45FvsI|G2B2",temp,value = T)


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
    t.tab2   <- t.tab %>% 
                    dplyr::select(.,X,Direction,logFC,P.Value) %>%
                    mutate(.,assay = name.tag) %>%
                    filter(.,Direction != "NONE") %>%
                    mutate(.,path_dir = paste0(Direction,": ", X))
                    
    big.tab.list[[i]] <- t.tab2
    big.tab.list2[[i]] <- t.tab
}

names(big.tab.list) <- fix_header_bigTab((big.tab.list))


plot_list <- list()


order <- c("MSBB AD vs Ctrl",
      "ROSMAP AD vs Ctrl",
      "Aβ42-H vs Aβ40-H",
      "Aβ42-H vs Ctrl",
      "A5 vs Ctrl",
      "D4 vs Ctrl",
      "H10 vs Ctrl")

big.tab.list3 <- big.tab.list[order]
for (i in seq_along(big.tab.list3)){
  
  temp <- merge(big.tab.list[["Mayo AD vs Ctrl"]],
                big.tab.list3[[i]],
                by = "X")
  
  plot_list[[i]]  <- plot_agreement(temp, name_tag = NULL)
  
}


p1 <-wrap_plots(c((plot_list)), ncol = 1)

p1 <- p1 + plot_annotation(
  title = 'Mayo AD vs Ctrl'
) & 
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "plain"))

plot_list <- list()
#plot_list[[1]] <-  plot_spacer() 


order <- c("MSBB AD vs Ctrl",
           "ROSMAP AD vs Ctrl",
           "Aβ42-H vs Aβ40-H",
           "Aβ42-H vs Ctrl",
           "A5 vs Ctrl",
           "D4 vs Ctrl",
           "H10 vs Ctrl")

# 
# order <- c(
#            "ROSMAP AD vs Ctrl",
#            "Aβ 42H / Aβ 40H ",
#            "Aβ 42H vs Ctrl",
#            "A5 vs Ctrl",
#            "D4 vs Ctrl",
#            "H10 vs Ctrl")

big.tab.list3 <- big.tab.list[order]
for (i in seq_along(big.tab.list3)){
  
  temp <- merge(big.tab.list[["MSBB AD vs Ctrl"]],
                big.tab.list3[[i]],
                by = "X")
  
  plot_list[[i]]  <- plot_agreement(temp, name_tag = NULL)
  
}



p2 <-wrap_plots(c((plot_list)), ncol = 1)


p2 <- p2 + plot_annotation(
  title = 'MSBB AD vs Ctrl'
) & 
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "plain"))

plot_list <- list()
order <- c("MSBB AD vs Ctrl",
           "ROSMAP AD vs Ctrl",
           "Aβ42-H vs Aβ40-H",
           "Aβ42-H vs Ctrl",
           "A5 vs Ctrl",
           "D4 vs Ctrl",
           "H10 vs Ctrl")

big.tab.list3 <- big.tab.list[order]
for (i in seq_along(big.tab.list3)){
  
  temp <- merge(big.tab.list[["ROSMAP AD vs Ctrl"]],
                big.tab.list3[[i]],
                by = "X")
  
  plot_list[[i]]  <- plot_agreement(temp, name_tag = NULL)
  
}


p3 <-wrap_plots(c((plot_list)), ncol = 1)
p3 <- p3 + plot_annotation(
  title = 'ROSMAP AD vs Ctrl'
) & 
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "plain"))


# p3 <- wrap_elements(grid::textGrob(' '))  / p3



(p1 | p2 |p3) 



plot_list <- list()
order <- c("MSBB AD vs Ctrl",
           "ROSMAP AD vs Ctrl",
           "Aβ42-H vs Aβ40-H",
           "Aβ42-H vs Ctrl",
           "A5 vs Ctrl",
           "D4 vs Ctrl",
           "H10 vs Ctrl")

big.tab.list3 <- big.tab.list[order]
for (i in seq_along(big.tab.list3)){
  
  temp <- merge(big.tab.list[["ROSMAP AD vs Ctrl"]],
                big.tab.list3[[i]],
                by = "X")
  
  plot_list[[i]]  <- plot_agreement(temp, name_tag =  names(big.tab.list3)[i])
  
}


p4 <-wrap_plots(c((plot_list)), ncol = 1)
p4 <- p4 + plot_annotation(
  title = 'ROSMAP AD vs Ctrl'
) & 
  theme(plot.title = element_text(size = 20, hjust = 0.5, face = "plain"))



 cairo_pdf("figures/03_c_Msbb_and_rest_Chi.pdf", width = 3.4, height = 7*1.5, onefile = T)
 p1
 p2
 p3
 p4
 dev.off()
# 