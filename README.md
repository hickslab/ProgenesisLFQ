# ProgenesisLFQ

```{r}

# Load packages
library(devtools)
library(tidyverse)
library(Biostrings)
library(imp4p)

# Load functions 
source_url("https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/EWM_ProgenesisLFQ_Global.R")

# Load raw data
protm <- "https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/20180502_WOS52_Cr_UPS_protm.csv" %>%
  load_data()

# Define normalized abundance columns
samples <- 10:21

# Workflow
data <- protm %>%
  filter_contaminants() %>%
  filter_peptides(peptides = 2, unique = 1) %>%
  split_group() %>%
  simplify_cols(variable = "Accession", group = samples)
  
# Write parsed data to file
#write_csv(data, "20180821_EWM_AGG3-2_Global_filtered_threshold2.csv")

```
