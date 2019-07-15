load_data <- function(name){
  data <- read.csv(name, skip = 2, stringsAsFactors = FALSE)
  
}


load_database <- function(name, split = " "){
  suppressPackageStartupMessages(library(Biostrings))
  
  # Load database
  data <- readAAStringSet(name)
  
  # Split accession by pattern and keep LHS
  names(data) <- sapply(strsplit(names(data), split), head, 1)
  
  return(data)
  
}


filter_score <- function(data, threshold = 13){
  data[data$Score == "---", "Score"] <- NA
  
  data %>%
    mutate(Score = as.numeric(as.character(Score))) %>%
    filter(Score > threshold)
  
}


filter_contaminants <- function(data){
  data %>%
    filter(Description != "cRAP")
  
}


map_protm <- function(data, protm){
  protm %>%
    select(Accession, Peptide.count, Unique.peptides, Confidence.score) %>%
    separate_rows(., Accession, sep = ";") %>%
    left_join(data, ., by = "Accession")
  
}


reduce_features <- function(data){
  data %>%
    # Summarize features with identical peptides - different accessions
    group_by(X., Sequence, Modifications, Score) %>%
    top_n(1, Unique.peptides) %>%
    top_n(1, Confidence.score) %>%
    slice_(1) %>%
    
    # Summarize features with multiple peptides
    group_by(X.) %>%
    top_n(1, Score) %>%
    slice_(1) %>%
    ungroup()
  
}


filter_redox <- function(pepm, reduced = "IAM", oxidized = "Empty"){
  if (reduced == "IAM"){
    mod <- "Carbamidomethyl"
    
  } else if (reduced == "NEM"){
    mod <- "Nethylmaleimide"
    
  }
  
  # Filter dataset for Cys-containing peptides only
  temp.data <- pepm %>%
    filter(., grepl("C", Sequence))
  
  # Filter dataset for oxidized Cys-containing peptides only
  if (oxidized == "Empty"){
    temp.data %>%
      filter(., str_count(Sequence, "C") > str_count(Modifications, mod))
    
  }
}


get_identifier_redox <- function(data, database, reduced = "IAM"){
  if (reduced == "IAM"){
    mod <- "Carbamidomethyl"
    
  } else if (reduced == "NEM"){
    mod <- "Nethylmaleimide"
    
  }
  
  # Map integer position of modification
  temp.data <- data %>%
    mutate(Modified = str_split(as.character(Modifications), "\\|") %>%
             lapply(., str_subset, mod) %>%
             lapply(., str_extract, "(?<=\\[)(.*)(?=\\])") %>%
             lapply(., as.numeric))
  
  # Map position of residue to peptide
  temp.data <- temp.data %>%
    rowwise() %>%
    mutate(Residue = str_locate_all(Sequence, "C")[[1]][, 1] %>% list())
  
  # Map unmodified position of residue on peptide
  temp.data <- temp.data %>%
    mutate(Residue = setdiff(Residue, Modified) %>% list())
  
  # Map position of peptide to protein
  temp.data <- temp.data %>%
    mutate(Start = str_locate(database[Accession], Sequence)[[1]])
  
  # Build identifier
  temp.data <- temp.data %>%
    mutate(Identifier = (Residue + Start - 1) %>%
             paste("C", ., sep = "") %>%
             paste(., collapse = "-") %>%
             paste(Accession, ., sep = "--"))
  
  # Return clean dataframe
  temp.data <- temp.data %>%
    select(-Modified, -Residue, -Start) %>%
    ungroup()
  
}


get_identifier <- function(data, database, modification = "Phospho"){
  # Map integer position of modification
  temp.data <- data %>%
    mutate(Identifier = str_split(as.character(Modifications), "\\|") %>%
             lapply(., str_subset, modification) %>%
             lapply(., str_extract, "(?<=\\[)(.*)(?=\\])") %>%
             lapply(., as.numeric))
  
  # Map modified residue by position
  temp.data <- temp.data %>%
    rowwise() %>%
    mutate(Residue = str_sub(Sequence, start = Identifier, end = Identifier) %>%
             list())
  
  # Map position of peptide to protein
  temp.data <- temp.data %>%
    mutate(Start = str_locate_all(database[Accession], as.character(Sequence))[[1]][[1]] %>%
             list())
  
  # Build identifier
  temp.data <- temp.data %>%
    mutate(Identifier = (Identifier + Start - 1) %>%
             paste(Residue, ., sep = "") %>%
             paste(., collapse = "-") %>%
             paste(Accession, ., sep = "--"))
  
  # Return clean dataframe
  temp.data <- temp.data %>%
    select(-Residue, -Start) %>%
    ungroup()
  
}


reduce_identifiers <- function(data, group){
  # Reduce to unique identifiers represented by highest score
  temp.data <- data %>%
    group_by(Identifier) %>%
    top_n(1, Score) %>%
    slice_(1) %>%
    arrange(Identifier)
  
  # Column-wise sum by group
  temp.issue <- data %>%
    group_by(Identifier) %>%
    select(Identifier, unlist(group)) %>%
    summarize_all(., funs(sum))
  
  # Replace summed abundance in ordered dataset
  temp.data[, unlist(group)] <- temp.issue[, -1]
  
  temp.data %>%
    ungroup()
  
}

