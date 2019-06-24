# ProgenesisLFQ --------------------------------------------------------------

# Load packages
library(devtools)
library(tidyverse)
library(Biostrings)


# Or install packages
#install.packages("devtools")
#install.packages("tidyverse")
#install.packages("BiocManager"); BiocManager::install("Biostrings")


# Set working directory
setwd("./") # ???


# Load functions 
source_url("https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/R/ProgenesisLFQ_Peptide.R")


# Load Progenesis peptide measurements
pepm <- "####.csv" %>% # ???
  load_data()


# Load Progenesis protein measurements
protm <- "###.csv" %>% # ???
  load_data()


# Load FASTA protein database
database <- "https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/data/Cr_v5.5_mt_chl_crap_20160502.fasta" %>%
  load_database()


# Define normalized abundance columns in pepm
samples <- 19:34 # ???


# Process data
data <- pepm %>%
  filter_score() %>%
  filter_contaminants() %>%
  map_protm(protm) %>%
  reduce_features() %>%
  filter_phospho() %>%
  get_identifier_phospho(database = database) %>%
  reduce_identifiers(group = samples) %>%
  select(Identifier, samples) %>%
  data.frame()


# Write parsed data to file
#write_csv(data, "###_filtered.csv")