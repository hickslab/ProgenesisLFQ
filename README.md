# ProgenesisLFQ

### Load packages
```{r}
library(devtools)
library(tidyverse)
library(Biostrings)
library(imp4p)
```

### Load functions 
```{r}
source_url("https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/EWM_ProgenesisLFQ_Global.R")
```

### Load raw data
```{r}

```


### Basic workflow
```{r}
data <- protm %>%
  filter_contaminants() %>%
  filter_peptides(threshold = 2) %>%
  split_group() %>%
  simplify_cols(variable = "Accession", group = samples) %>%
  data.frame()
```

