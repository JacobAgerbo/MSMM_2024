library(tidyverse)
library(readxl)
library(ggpubr)

setwd("Desktop/MSMM_2024/")


FA <- read_excel("Fatty_Acids.xlsx")

#pdf("Fatty_Acids.pdf")
FA %>%
  mutate(Treatment = as.factor(Treatment)) %>%
  select(Treatment, contains("%")) %>%
  reshape2::melt() %>%
  ggboxplot(x = "Treatment",
            y = "value",
            fill = "Treatment") +
  facet_wrap(~variable, scales = "free_y") +
  stat_compare_means(label.x = 1.5, label = "p.signif") +
  scale_fill_manual(values = c("#eb4034", "#28478f")) +
  theme(legend.position = "none") +
  ylab("Sum-Normalised Abundace")
#dev.off()

md <- read_excel("metadata.xlsx")

MT <- read.csv("MetaTranscriptome.csv",
               row.names = 1)







# # Set the directory path where your data frames are located
# directory_path <- "KOs_across_MAGs/"
# 
# # Get a list of file names in the directory
# file_names <- list.files(path = directory_path, pattern = "\\.txt$", full.names = TRUE)
# 
# names <- sub(directory_path, "", basename(file_names))
# names <- sub("-gene_calls.txt", "", basename(names))
# # Read all data frames from the files into a list of data frames
# dfs_list <- lapply(file_names, function(file) {
#   read.csv(file, sep = "\t")
# })
# 
# names(dfs_list) <- names
# # Subtract the column from each dataframe in the list
# # Subtract the KOfam column from each data frame
# modified_dataframes <- lapply(dfs_list, function(df) {
#   df <- df$KOfam..ACCESSION.
#   return(df)
# })
# 
# KO_df <- plyr::ldply(modified_dataframes, data.frame) %>%
#   rename(MAG = ".id", koFAM = "X..i..") %>% 
#   group_by(MAG) %>%
#   distinct(koFAM) %>%
#   mutate(Presence = ifelse(koFAM == "", 0, 1)) %>%
#   pivot_wider(names_from = MAG, id_cols = koFAM, values_from = Presence) %>% 
#   filter(koFAM != "") %>%
#   column_to_rownames(var = "koFAM") %>%
#   mutate_all(~ if_else(is.na(.), 0, .))
# 
# 
# 
# KO_df <- KO_df %>%
#   rownames_to_column(var = "koFAM") %>%
#   filter(koFAM %in% rownames(MT))
# 
# 
# 
# MT <- MT %>%
#   rownames_to_column(var = "koFAM") %>%
#   filter(koFAM %in% KO_df$koFAM)
# 
# 
# FUNCTIONS <- read.csv("FUNCTIONS.txt", sep = "\t") %>%
#   filter(source == "KOfam")
# 
# 
# gene_cov <- read.csv("FUNCTIONS-GENE-COVERAGES.txt", sep = "\t") %>%
#   filter(key %in% FUNCTIONS$gene_callers_id) %>%
#   rename(gene_callers_id = "key") %>%
#   full_join(FUNCTIONS, by = "gene_callers_id") %>%
#   select(accession,contains("ERR"))
# 
# KO_cov <- gene_cov %>%
#   reshape2::melt() %>%
#   dplyr::group_by(accession, variable) %>%
#   summarise(value = sum(value)) %>%
#   pivot_wider(names_from = variable, id_cols = accession) %>%
#   filter(accession %in% KO_df$koFAM)
# 


MG <- read.csv("SUMMARY/bins_across_samples/mean_coverage_Q2Q3.txt", sep = "\t")
