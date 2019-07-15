reduce_features <- function(df){
  df %>%
    # Summarize features with identical peptides - different accessions
    group_by(`#`, Sequence, Modifications, Score) %>%
    top_n(1, `Unique peptides`) %>%
    top_n(1, `Confidence score`) %>%
    dplyr::slice(1) %>%
    
    # Summarize features with multiple peptides
    group_by(`#`) %>%
    top_n(1, Score) %>%
    dplyr::slice(1) %>%
    ungroup()
  
}


get_identifier <- function(df, database, mod = "Phospho"){
  # Separate to one modification per row
  temp.df <- df %>%
    select(1, Accession, Sequence, Modifications) %>%
    separate_rows(Modifications, sep = "\\|") %>%
    filter(str_detect(Modifications, mod))
  
  # Get the modification position and residue from the peptide
  temp.df <- temp.df %>%
    mutate(position = str_extract(Modifications, "(?<=\\[)(.*)(?=\\])") %>% as.numeric()) %>%
    mutate(residue = str_sub(Sequence, start = position, end = position))
  
  # Get the modification site from the protein
  temp.df <- temp.df %>%
    mutate(locate = str_locate(database[Accession], Sequence)[, 1]) %>%
    mutate(site = locate + position - 1)
  
  # Group by feature and build identifier
  temp.df <- temp.df %>%
    group_by(`#`) %>%
    mutate(Identifier = paste(residue, site, sep = "") %>%
             paste(., collapse = "-") %>%
             paste(Accession, ., sep = "--")) %>%
    ungroup()
  
  # Join identifier column to input data
  temp.df %>%
    select(1, Identifier) %>%
    distinct() %>%
    inner_join(df, ., by = "#")
  
}


reduce_identifiers <- function(df, group){
  # Reduce to unique identifiers represented by highest score
  temp.df <- df %>%
    group_by(Identifier) %>%
    top_n(1, Score) %>%
    slice_(1) %>%
    arrange(Identifier)
  
  # Column-wise sum by group
  temp.issue <- df %>%
    group_by(Identifier) %>%
    select(Identifier, unlist(group)) %>%
    summarize_all(., funs(sum))
  
  # Replace summed abundance in ordered dataset
  temp.df[, unlist(group)] <- temp.issue[, -1]
  
  temp.df %>%
    ungroup()
  
}


filter_redox <- function(df, reduced = "Nethylmaleimide"){
  # Default Mascot variable modifications:
  # IAM = Carbamidomethyl
  # NEM = Nethylmaleimide
  
  temp.df <- df %>%
    mutate(Modifications = replace_na(Modifications, "")) %>%
    filter(., str_count(Sequence, "C") > str_count(Modifications, reduced))
  
}


get_identifier_redox <- function(df, database, reduced = "Nethylmaleimide"){
  # Default Mascot variable modifications:
  # IAM = Carbamidomethyl
  # NEM = Nethylmaleimide
  
  # Map integer position of modification
  temp.df <- df %>%
    mutate(Modified = str_split(as.character(Modifications), "\\|") %>%
             lapply(., str_subset, mod) %>%
             lapply(., str_extract, "(?<=\\[)(.*)(?=\\])") %>%
             lapply(., as.numeric))
  
  # Map position of residue to peptide
  temp.df <- temp.df %>%
    rowwise() %>%
    mutate(Residue = str_locate_all(Sequence, "C")[[1]][, 1] %>% list())
  
  # Map unmodified position of residue on peptide
  temp.df <- temp.df %>%
    mutate(Residue = setdiff(Residue, Modified) %>% list())
  
  # Map position of peptide to protein
  temp.df <- temp.df %>%
    mutate(Start = str_locate(database[Accession], Sequence)[[1]])
  
  # Build identifier
  temp.df <- temp.df %>%
    mutate(Identifier = (Residue + Start - 1) %>%
             paste("C", ., sep = "") %>%
             paste(., collapse = "-") %>%
             paste(Accession, ., sep = "--"))
  
  # Return clean dataframe
  temp.df <- temp.df %>%
    select(-Modified, -Residue, -Start) %>%
    ungroup()
  
}



get_identifier_redox_tidy <- function(df, database, reduced = "Nethylmaleimide"){
  # Default Mascot variable modifications:
  # IAM = Carbamidomethyl
  # NEM = Nethylmaleimide
  
  # Map integer position of modification
  temp.df <- df %>%
    mutate(Modified = str_split(as.character(Modifications), "\\|") %>%
             lapply(., str_subset, mod) %>%
             lapply(., str_extract, "(?<=\\[)(.*)(?=\\])") %>%
             lapply(., as.numeric))
  
  # Map position of residue to peptide
  temp.df <- temp.df %>%
    rowwise() %>%
    mutate(Residue = str_locate_all(Sequence, "C")[[1]][, 1] %>% list())
  
  # Map unmodified position of residue on peptide
  temp.df <- temp.df %>%
    mutate(Residue = setdiff(Residue, Modified) %>% list())
  
  # Map position of peptide to protein
  temp.df <- temp.df %>%
    mutate(Start = str_locate(database[Accession], Sequence)[[1]])
  
  # Build identifier
  temp.df <- temp.df %>%
    mutate(Identifier = (Residue + Start - 1) %>%
             paste("C", ., sep = "") %>%
             paste(., collapse = "-") %>%
             paste(Accession, ., sep = "--"))
  
  # Return clean dataframe
  temp.df <- temp.df %>%
    select(-Modified, -Residue, -Start) %>%
    ungroup()
  
}
