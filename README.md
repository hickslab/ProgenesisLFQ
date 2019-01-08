# ProgenesisLFQ

### Copy and paste into local RStudio session
```{r}
# Load packages
library(devtools)
library(tidyverse)

# Load functions 
source_url("https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/EWM_ProgenesisLFQ_Global.R")

# Load raw data
protm <- "https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/20180502_WOS52_Cr_UPS_protm.csv" %>%
  load_data()

# Define columns with normalized abundance values
samples <- 10:21

# Workflow
data <- protm %>%
  filter_contaminants() %>%
  filter_peptides(peptides = 2, unique = 1) %>%
  split_group() %>%
  simplify_cols(variable = "Accession", group = samples)
  
# Write processed data to output file
#write_csv(data, "###.csv")

```
