# Load packages
library(imp4p)


# Or install packages
#install.packages("imp4p")


# Load functions 
#source_url("https://raw.githubusercontent.com/hickslab/ProgenesisLFQ/master/R/StatLFQ.R")
source("Z:/Lab_Members/Evan/Scripts/R/LFQ/EWM_StatLFQ_v1.R")


# Simplified data after processing with ProgenesisLFQ workflow
protm <- "20180821_EWM_AGG3-2_Global_filtered_threshold2.csv" %>%
  load_simple()


# Condition indeces in data
a <- 2:5; b <- 6:9; c <- 10:13

group <- list(a, b, c)

group.compare <- list(list(a, b), # A vs B
                      list(a, c), # A vs C
                      list(b, c)) # B vs C


# Rename columns based on design
colnames(protm) <- get_design(group = group) %>%
  c("Accession", .)


# Statistical analysis
data <- protm %>%
  clean_min(., group = group, nonzero = 3) %>%
  transform_data(., group = group, method = "log2") %>%
  impute_imp4p(., group = group, method = "log2") %>%
  calculate_ttest(., group.compare = group.compare, fdr = TRUE) %>%
  calculate_fc(., group.compare = group.compare, method = "log2", difference = TRUE)


# Write output file
#write_csv(data, "###_filtered_clean3_log2_imp4p.csv")
