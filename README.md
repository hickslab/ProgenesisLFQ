
# *`ProgenesisLFQ`*

**This repository contains functions to process Progenesis QI for proteomics v2.0 data.**

Select and copy a workflow file below. In RStudio, find '???' in the workflow and update with your particular filename/path before running.**

To run, you need the following data:

* Protein measurements:
  + [Global proteomics](https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/workflow/ProgenesisLFQ_Global_Workflow.R)
  
* Peptide measurements, protein measurements, and FASTA database:
  + [Phosphoproteomics](https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/workflow/ProgenesisLFQ_Phospho_Workflow.R)
  + [Redox proteomics](https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/workflow/ProgenesisLFQ_Redox_Workflow.R)


# *`StatLFQ`*

**Statistical workflow for differential analysis.**

Takes data in simple format (variable column followed by group columns).

Append needed code to your workflow and run:

```{r}
# Load packages
library(imp4p)


# Load functions 
source_url("https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/R/StatLFQ.R")


# Condition indeces in data
a <- 2:5; b <- 6:9
group <- list("0" = a, "60" = b)


# Process data
data2 <- data %>%
  rename_columns(., group) %>%
  clean_min(., group, nonzero = 3) %>%
  transform_data(., group, method = "log2") %>%
  impute_imp4p(., group)
```

* Pairwise *t*-test

```{r}
group.compare <- list(list(a, b))


# Process data
data3 <- data2 %>%
  calculate_ttest(., group.compare, fdr = TRUE) %>%
  calculate_fc(., group.compare, difference = TRUE)
```

* One-way ANOVA

```{r}
# StatLFQ (ANOVA) -----------------------------------------------------------

data4 <- data2 %>%
  calculate_1anova() %>%
  calculate_fc(., group.compare, difference = TRUE) %>%
  add_max_fc()
```


# *`AnnotateLFQ`*

source_url(https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/R/AnnotateLFQ.R)

* [Clean Accession/Identifier]

* [GO](https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/R/AnnotateLFQ.R)

* [KEGG](https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/R/AnnotateLFQ.R)

# *`PlotLFQ`*

source_url(https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/R/PlotLFQ.R)