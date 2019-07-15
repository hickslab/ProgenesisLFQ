# PhosphoLFQ --------------------------------------------------------------
# Phosphosite-level label-free quantification workflow


# Load packages
library(Biostrings)
library(devtools)
library(tidyverse)


# Or install packages
#install.packages("devtools")
#install.packages("tidyverse")
#install.packages("BiocManager"); BiocManager::install("Biostrings")


# Set the working directory.
# The example below is based on downloading the package from GitHub.
setwd("~/ProgenesisLFQ/data") # ???


# Load functions
url <- "https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/"
source_url(paste0(url, "R/ProcessLFQ.R"))


# Load Peptide Measurements spreadsheet exported from Progenesis.
pepm <- "20180513_EWM_QP1_Phos_5600_pepm.csv" %>% #???
  read_csv(., skip = 2, col_types = cols())


# Load Peptide Measurements spreadsheet exported from Progenesis.
protm <- "20180513_EWM_QP1_Phos_5600_protm.csv" %>% #???
  read_csv(., skip = 2, col_types = cols()) %>%
  select(1:4) %>%
  separate_rows(., Accession, sep = ";")


# Load protein sequence database.
database <- "Cr_v5.5_mt_chl_crap_20160502.fasta" %>% # ???
  readAAStringSet(.)

names(database) <- str_split(names(database), " ", simplify = TRUE)[, 1]


# Define normalized abundance columns in pepm
samples <- 17:18 # ???


# Process data
data <- pepm %>%
  filter(Score > 13) %>%
  filter(Description != "cRAP") %>%
  left_join(., protm, by = "Accession") %>%
  reduce_features() %>%
  filter(str_detect(Modifications, "Phospho")) %>%
  get_identifier(., database, mod = "Phospho") %>%
  reduce_identifiers(., samples) %>%
  select(Identifier, samples) %>%
  data.frame()


# Write parsed data to file
#write_csv(data, "###_processed.csv")