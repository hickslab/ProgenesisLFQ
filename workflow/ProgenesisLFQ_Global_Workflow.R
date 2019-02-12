# Load packages
library(devtools)
library(tidyverse)


# Or install packages
#install.packages("devtools")
#install.packages("tidyverse")


# Set working directory where raw data is located
setwd("./") # ???


# Load functions 
source_url("https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/R/ProgenesisLFQ_Global.R")


# Load raw data
protm <- read_csv("https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/data/20180502_WOS52_Cr_UPS_protm.csv", col_types = cols(), skip = 2)

#protm <- "###.csv" %>% # ???
#  load_data()


# Define columns with normalized abundance values
samples <- 10:21 # ???
  
  
# Workflow
data <- protm %>%
  filter_contaminants() %>%
  filter_peptides(peptides = 2, unique = 1) %>%
  split_group() %>%
  select(Accession, samples) %>%
  data.frame()


# Write processed data to output file
#write_csv(data, "###.csv") # ???
