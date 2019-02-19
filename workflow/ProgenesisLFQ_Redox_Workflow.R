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
  filter_redox(reduced = "NEM") %>%
  get_identifier_redox(database = database, reduced = "NEM") %>%
  reduce_identifiers(group = samples) %>%
  select(Identifier, samples) %>%
  data.frame()


# Write processed data to output file
#write_csv(data, "###.csv") # ???