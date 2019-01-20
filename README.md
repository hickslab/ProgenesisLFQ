# `ProgenesisLFQ`

Last updated: 01-18-2019 (EWM)

**This repository contains functions to process Progenesis QI for proteomics v2.0 data:**
* Protein measurements
  + [Global proteomics](https://github.com/hickslab/ProgenesisLFQ#global-proteomics)
  
* Peptide measurements:
  + [Phosphoproteomics](https://github.com/hickslab/ProgenesisLFQ#phosphoproteomics)
  + [Redox proteomics](https://github.com/hickslab/ProgenesisLFQ#redox-proteomics)

## Global proteomics
* Must edit in local session:
  + Set working directory
  + Raw data filename in working directory
  + Column indeces with normalized abundance values

**Copy script below and paste into your local RStudio session:**

```{r}
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

```

## Phosphoproteomics
* Must edit in local session:
  + Set working directory
  + Raw data filename in working directory
  + Column indeces with normalized abundance values

**Copy script below and paste into your local RStudio session:**

```{r}



```

## Redox proteomics
* Must edit in local session:
  + Set working directory
  + Raw data filename in working directory
  + Column indeces with normalized abundance values

**Copy script below and paste into your local RStudio session:**

```{r}
###
```

