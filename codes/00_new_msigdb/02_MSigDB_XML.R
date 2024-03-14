# Make a table to create clean names for MSigDB names.
rm(list = ls())
library(msigdb)
library(ExperimentHub)
library(GSEABase)


library(xml2)
library(dplyr)

pg <-read_xml("~/Downloads/msigdb_v6.2_files_to_download_locally/msigdb_v6.2.xml")


node<-xml_find_all(pg, xpath = "//GENESET")

temp <- xml_attrs(node[[1]])

library(tidyverse)

records_data <- list()

# Loop through each record and extract attributes
for (i in seq_along(node)) {
    record <- node[i]
    
    # Extract attributes from the record
    attrs <- xml_attrs(record)
    
    # Add the extracted attributes as a new list element
    records_data[[i]] <- as.list(attrs)
}

df <- bind_rows(records_data)

temp <- df[grep("CP",df$SUB_CATEGORY_CODE),]
temp <- temp[,c("STANDARD_NAME","CONTRIBUTOR","DESCRIPTION_BRIEF")]

temp$CONTRIBUTOR <- gsub("Pathway Interaction Database","PID",temp$CONTRIBUTOR)
temp$CONTRIBUTOR <- gsub("Alexandra Naba","Naba",temp$CONTRIBUTOR)

temp$display_name <- gsub("^Genes involved in ","",temp$DESCRIPTION_BRIEF)
temp$display_name <- paste0("(",temp$CONTRIBUTOR,") ",temp$display_name)

write.csv(temp,"preprocessed/MSigDB_display_names.csv",row.names = F)
