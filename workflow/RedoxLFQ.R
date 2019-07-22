# RedoxLFQ
# Reversibly oxidized Cys-level label-free quantification workflow


# ProcessLFQ ----


# Packages
library(Biostrings)
library(devtools)
library(tidyverse)


# Or install packages
#install.packages("BiocManager"); BiocManager::install("Biostrings")
#install.packages("devtools")
#install.packages("tidyverse")


# Set the working directory.
# The example below is based on downloading the package from GitHub.
setwd("~/ProgenesisLFQ/data") # ???


# Functions
url <- "https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/"
source_url(paste0(url, "R/ProcessLFQ.R"))


# Load Peptide Measurements spreadsheet exported from Progenesis.
pepm <- "20190123_EWM_AZD1_R_rank-lessthan11-include_uniprot_pepm.csv" %>% # ???
  read_csv(., skip = 2, col_types = cols())


# Load Peptide Measurements spreadsheet exported from Progenesis.
protm <- "20190123_EWM_AZD1_R_rank-lessthan11-include_uniprot_protm.csv" %>% # ???
  read_csv(., skip = 2, col_types = cols()) %>%
  select(1:4) %>%
  separate_rows(., Accession, sep = ";")


# Load protein sequence database.
database <- "Cr_uniprot_crap_20190130.fasta" %>% # ???
  readAAStringSet(.)

names(database) <- str_split(names(database), " ", simplify = TRUE)[, 1]


# Define normalized abundance columns in pepm
samples <- 19:34 # ???


# Process data.
data <- pepm %>%
  filter(Score > 13) %>%
  filter(Description != "cRAP") %>%
  left_join(., protm, by = "Accession") %>%
  reduce_features() %>%
  filter_redox(., reduced = "Nethylmaleimide") %>%
  get_identifier_redox(., database, reduced = "Nethylmaleimide") %>%
  reduce_identifiers(., samples) %>%
  select(Identifier, samples) %>%
  data.frame()


# Write parsed data to file
#write_csv(data, "###_processed.csv") # ???


# StatLFQ ----


# Packages
library(imp4p)
library(broom)


# Functions
url <- "https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/"
source_url(paste0(url, "R/StatLFQ.R"))


# Workflow
# Define the column indeces for replicates in each condition.
a <- 2:5; b <- 6:9; c <- 10:13; d <- 14:17

group <- list("0" = a, "15" = b, "30" = c, "60" = d) # ???

group.compare <- list("0-15" = list(a, b),
                      "0-30" = list(a, c),
                      "0-60" = list(a, d)) # ???


# Rename the abundance columns in a simplified "Condition-Replicate" format.
data <- data %>%
  rename_columns(., group)


# Clean, transform, and impute abundance columns
data2 <- data %>%
  clean_min(., group, nonzero = 3) %>% # ???
  transform_data(., group, method = "log2") %>%
  impute_imp4p(., group)


# Pairwise *t*-test
data3 <- data2 %>%
  calculate_ttest(., group.compare, fdr = TRUE)


# One-way ANOVA
data4 <- data2 %>%
  calculate_1anova()


# Fold change
data3 <- data3 %>%
  calculate_fc(., group.compare, difference = TRUE) %>%
  add_fc_max()

data4 <- data4 %>%
  calculate_fc(., group.compare, difference = TRUE) %>%
  add_fc_max()


# AnnotateLFQ ----


# Functions
url <- "https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/"
source_url(paste0(url, "R/AnnotateLFQ.R"))


# UniProt
data5 <- data4 %>%
  add_missingness(., data, group)


# PlotLFQ ----


# Functions
url <- "https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/"
source_url(paste0(url, "R/PlotLFQ.R"))


# PCA
# Principal component analysis (PCA)
data2 %>% plot_pca(., group) +
  theme_custom()


# Volcano Plot
data3 %>%
  plot_volcano(.,
               group,
               group.compare,
               fdr = TRUE,
               threshold = 2,
               xlimit = 8,
               ylimit = 5) +
  theme_custom()


# Trend Profiles
data4 %>%
  filter(FDR < 0.05) %>% # ???
  filter(abs(`0-60_FC`) >= 1) %>% # ???
  plot_hclust(., group, k = 2) +
  theme_custom()
