# Load packages
library(devtools)
library(tidyverse)


# Or install packages
#install.packages("devtools")
#install.packages("tidyverse")


# Set working directory where raw data is located
setwd("###")


# Load functions 
source_url("https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/EWM_ProgenesisLFQ_Global.R")


# Load raw data
protm <- "###.csv" %>%
  load_data()


# Define columns with normalized abundance values
samples <- #:#
  
  
  # Workflow
  data <- protm %>%
  filter_contaminants() %>%
  filter_peptides(peptides = 2, unique = 1) %>%
  split_group() %>%
  simplify_cols(variable = "Accession", group = samples)


# Write processed data to output file
#write_csv(data, "###.csv")
