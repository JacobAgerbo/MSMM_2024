library(tidyverse)
library(readxl)
library(ggpubr)
library(Hmisc)
setwd("Desktop/MSMM_2024/")


FA <- read_excel("Fatty_Acids.xlsx")
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

md <- read_excel("metadata.xlsx")

MT <- read.csv("MetaTranscriptome.csv",
               row.names = 1)







# Set the directory path where your data frames are located
directory_path <- "KOs_across_MAGs/"

# Get a list of file names in the directory
file_names <- list.files(path = directory_path, pattern = "\\.txt$", full.names = TRUE)

names <- sub(directory_path, "", basename(file_names))
names <- sub("-gene_calls.txt", "", basename(names))
# Read all data frames from the files into a list of data frames
dfs_list <- lapply(file_names, function(file) {
  read.csv(file, sep = "\t")
})

names(dfs_list) <- names
# Subtract the column from each dataframe in the list
# Subtract the KOfam column from each data frame
modified_dataframes <- lapply(dfs_list, function(df) {
  df <- df$KOfam..ACCESSION.
  return(df)
})

KO_df <- plyr::ldply(modified_dataframes, data.frame) %>%
  rename(MAG = ".id", koFAM = "X..i..") %>%
  group_by(MAG) %>%
  distinct(koFAM) %>%
  mutate(Presence = ifelse(koFAM == "", 0, 1)) %>%
  pivot_wider(names_from = MAG, id_cols = koFAM, values_from = Presence) %>%
  filter(koFAM != "") %>%
  column_to_rownames(var = "koFAM") %>%
  mutate_all(~ if_else(is.na(.), 0, .))



KO_df <- KO_df %>%
  rownames_to_column(var = "koFAM") %>%
  filter(koFAM %in% rownames(MT))

MT <- MT %>%
  rownames_to_column(var = "koFAM") %>%
  filter(koFAM %in% KO_df$koFAM)
# 
# 
FUNCTIONS <- read.csv("FUNCTIONS.txt", sep = "\t") %>%
  filter(source == "KOfam")

gene_cov <- read.csv("FUNCTIONS-GENE-COVERAGES.txt", sep = "\t") %>%
  filter(key %in% FUNCTIONS$gene_callers_id) %>%
  rename(gene_callers_id = "key") %>%
  full_join(FUNCTIONS, by = "gene_callers_id") %>%
  select(accession,contains("ERR"))

KO_cov <- gene_cov %>%
  reshape2::melt() %>%
  dplyr::group_by(accession, variable) %>%
  summarise(value = sum(value)) %>%
  pivot_wider(names_from = variable, id_cols = accession) %>%
  filter(accession %in% KO_df$koFAM)
# 

KO_FUNCTIONS <- read.csv("FUNCTIONS.txt", sep = "\t") %>%
  filter(source == "KOfam") %>%
  filter(accession %in% KO_cov$accession) %>%
  rename("ID" = gene_callers_id,
         "KOfam" = function.) %>%
  select(ID, accession, KOfam)

COG_CAT <- read.csv("FUNCTIONS.txt", sep = "\t") %>%
  filter(gene_callers_id %in% KO_FUNCTIONS$ID) %>%
  filter(source == "COG20_CATEGORY") %>%
  rename("ID" = gene_callers_id) %>%
  select(ID, function.) %>%
  rename("COG_CATEGORY" = function.) %>%
  mutate(COG_CATEGORY = str_remove(COG_CATEGORY, "!!!.*")) %>%
  distinct(ID, .keep_all = TRUE) 

FUNCTIONS <- KO_FUNCTIONS %>%
  full_join(COG_CAT, by = "ID") %>%
    distinct(accession, .keep_all = TRUE) %>%
  select(-ID)


MT_df <- MT %>%
  dplyr::mutate(koFAM = paste("MT_", koFAM, sep ="")) %>%
  column_to_rownames(var = "koFAM")
colnames(MT_df) <- md$Genome

KO_cov_df <- KO_cov %>%
  dplyr::mutate(accession = paste("MG_", accession, sep ="")) %>%
  rename("koFAM" = accession) %>%
  column_to_rownames(var = "koFAM")

MG_df <- read.csv("SUMMARY/bins_across_samples/mean_coverage_Q2Q3.txt", sep = "\t") %>%
  column_to_rownames(var = "bins")

FA_df <- FA %>%
  mutate(Sample = md$Genome) %>%
  select(Sample, contains("%")) %>%
  reshape2::melt() %>%
  pivot_wider(names_from = Sample, id_cols = variable) %>%
  column_to_rownames(var = "variable")

desired_order <- colnames(MG_df)
# Reorder the columns in the data frame to match the desired order
FA_df <- FA_df %>%
  select(all_of(desired_order)) %>%
  hilldiv::tss()
MT_df <- MT_df %>%
  select(all_of(desired_order)) %>%
  hilldiv::tss()
KO_cov_df <- KO_cov_df %>%
  select(all_of(desired_order)) %>%
  hilldiv::tss()

MG_df <- MG_df %>%
  hilldiv::tss()

cor_df <- rbind(FA_df,KO_cov_df,MT_df,MG_df)

library(Hmisc)
cor_FA_x_MG_x_MT <- rcorr(t(cor_df))


df <- cor_FA_x_MG_x_MT$r %>%
  reshape2::melt()
colnames(df) <- c("Source", "Target", "Pearson_Rho")

rcorr_padjust <- function(x, method = "BH") {
  stopifnot(class(x) == "rcorr")
  x$P[upper.tri(x$P)] <- p.adjust(x$P[upper.tri(x$P)], method = method)
  x$P[lower.tri(x$P)] <- p.adjust(x$P[lower.tri(x$P)], method = method)
  return(x)}

df_p <- rcorr_padjust(cor_FA_x_MG_x_MT, method = "bonferroni")
df_adj_p <- df_p$P %>%
  reshape2::melt()
colnames(df_adj_p) <- c("Source", "Target", "adj_p_value")

df <- df %>%
  mutate(adj_p_value = df_adj_p$adj_p_value) %>%
  filter(Source != Target,
         adj_p_value < 0.05)
  
MG_FUNCTIONs <- FUNCTIONS %>%
  dplyr::mutate(accession = paste("MG_",accession, sep = ""))
MT_FUNCTIONs <- FUNCTIONS %>%
  dplyr::mutate(accession = paste("MT_",accession, sep = ""))

FUNCTIONS <- rbind(MG_FUNCTIONs,MT_FUNCTIONs) %>%
  rename("Source" = accession)

nodes <- df %>%
  select(Source) %>%
  mutate(Type = case_when(
    str_detect(Source, "MGY") ~ "MAG",
    str_detect(Source, "MG_") ~ "KO_coverage",
    str_detect(Source, "MT_") ~ "KO_Expression",
    str_detect(Source, "%") ~ "Fatty Acid",
    TRUE ~ "Other")) %>%
  distinct(Source, .keep_all = TRUE) %>%
  full_join(FUNCTIONS, by = "Source") %>%
  filter(!is.na(Type)) %>%
  rename("Id" = Source)


# Save the dataframe to a CSV file for Gephi
write.csv(df, "correlation_network_for_gephi.csv", row.names = FALSE)
write.csv(nodes, "nodes_for_gephi.csv", row.names = FALSE)
