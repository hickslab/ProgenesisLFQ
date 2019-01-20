# Load packages
library(devtools)
library(tidyverse)
library(Biostrings)


# Or install packages
#install.packages("devtools")

#install.packages("tidyverse")

#install.packages("BiocManager")
#BiocManager::install("Biostrings")


# Set working directory where raw data is located
setwd("###")


# Load functions 
source_url("https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/R/ProgenesisLFQ_Phospho.R")


# Load Progenesis peptide measurements
pepm <- "###.csv" %>%
  load_data()


# Load Progenesis protein measurements
protm <- "###.csv" %>%
  load_data()


# Load Protein sequence database
database <- "###.fasta" %>%
  load_database()


# Define normalized abundance columns in pepm
samples <- #:#


# Process data
data <- pepm %>%
  filter_score() %>%
  filter_contaminants() %>%
  map_protm(protm) %>%
  reduce_features() %>%
  filter_phospho() %>%
  get_identifier_phospho(database = database) %>%
  reduce_identifiers(group = samples) %>%
  simplify_cols(variable = "Identifier", group = samples) %>%
  remove_PACid(., variable = "Identifier")


# Write parsed data to file
#write_csv(data, "###_filtered.csv")