
# *`ProgenesisLFQ`*

**This repository contains functions to process Progenesis QI for proteomics v2.0 data.**

Select and copy a workflow file below.

In RStudio, find '???' in the workflow and update with your particular filename/path before running.

To run, you need the following data on your computer:

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
# StatLFQ (Imputation) ----------------------------------------------------


# Load packages
library(imp4p)


# Load functions 
source_url("https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/R/StatLFQ.R")


# Condition indeces in data
a <- 2:5; b <- 6:9; c <- 10:13
group <- list("25" = a, "50" = b, "100" = c)


# Process data
data2 <- data %>%
  rename_columns(., group) %>%
  clean_min(., group, nonzero = 3) %>%
  transform_data(., group, method = "log2") %>%
  impute_imp4p(., group)
```

* Pairwise *t*-test

```{r}
# StatLFQ (t-test) ---------------------------------------------------------


group.compare <- list(list(a, b),
                      list(a, c),
                      list(b, c))


# Process data
data3 <- data2 %>%
  calculate_ttest(., group.compare, fdr = TRUE) %>%
  calculate_fc(., group.compare, difference = TRUE)
```

* One-way ANOVA

```{r}
# StatLFQ (ANOVA) -----------------------------------------------------------


# Load packages
library(broom)


data4 <- data2 %>%
  calculate_1anova() %>%
  calculate_fc(., group.compare, difference = TRUE) %>%
  add_max_fc()
```


# *`AnnotateLFQ`*

**Annotation**


```{r}
# AnnotateLFQ -------------------------------------------------------------

source_url(https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/R/AnnotateLFQ.R)


phytozome <- source_url("https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/data/Creinhardtii_281_v5.5.annotation_info.txt")


data5 <- data3 %>%
  add_missingness(., data, group) %>%
  remove_PACid() %>%
  split_identifier() %>%
  add_phytozome(., phytozome)
```


# *`PlotLFQ`*

### **Plotting**

```{r}
# PlotLFQ ----------------------------------------------------------------


# Load functions 
source_url("https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/R/PlotLFQ.R")


# PCA
plot_pca(data2, group)


# Volcano plot
plot_volcano(data3,
             group,
             fdr = TRUE,
             threshold = 2,
             xlimit = 10,
             ylimit = 8)
```
