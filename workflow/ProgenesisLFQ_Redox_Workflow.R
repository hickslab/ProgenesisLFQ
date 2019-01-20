# Load packages
library(devtools)
library(tidyverse)
library(Biostrings)


# Or install packages
#install.packages("devtools")

#install.packages("tidyverse")

#install.packages("BiocManager")
#BiocManager::install("Biostrings")


# Set working directory
setwd("###")


# Load functions 
source_url("https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/R/EWM_ProgenesisLFQ_Redox.R")


# Load Progenesis peptide measurements
pepm <- "###" %>%
  load_data()


# Load Progenesis protein measurements
protm <- "###" %>%
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
  filter_redox(reduced = "IAM") %>%
  get_identifier_redox(database = database, reduced = "IAM") %>%
  reduce_identifiers(group = samples) %>%
  simplify_cols(variable = "Identifier", group = samples) %>%
  data.frame()


# Write parsed data to file
#write_csv(data, "###_filtered.csv")

