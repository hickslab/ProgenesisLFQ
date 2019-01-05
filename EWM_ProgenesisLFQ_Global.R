load_data <- function(name){
  data <- read.csv(name, skip = 2, stringsAsFactors = FALSE)
  
}


filter_contaminants <- function(data){
  data %>%
    filter(Description != "cRAP")
  
}


filter_peptides <- function(data, peptides = 2, unique = 1){
  data %>%
    filter(Peptide.count >= peptides, Unique.peptides >= unique)
  
}


split_group <- function(data){
  data %>%
    separate(Accession,
             sep = ";",
             into = c("Accession", "Group.members"),
             extra = "merge",
             fill = "right") %>%
    select(-Group.members, Group.members)
  
}


simplify_cols <- function(data, variable, group){
  temp.data <- data %>%
    select(variable)
  
  data[, unlist(group)] %>%
    bind_cols(temp.data, .)
  
}


remove_PACid <- function(data, variable){
  if (variable == "Accession"){
    # Conditional split and return of first element
    temp.data <- data %>%
      rowwise() %>%
      mutate(Accession = if_else(condition = str_detect(Accession, "Cre"),
                                 true = str_split(Accession, pattern = "\\|")[[1]][1],
                                 false = Accession))
    
  } else if (variable == "Identifier"){
    # Separate accession from modification sites
    temp.data <- data %>%
      separate(Identifier,
               into = c("Accession", "Sites"),
               sep = "--")
    
    # Conditional split and return of first element
    temp.data <- temp.data %>%
      rowwise() %>%
      mutate(Accession = if_else(condition = str_detect(Accession, "Cre"),
                                 true = str_split(Accession, pattern = "\\|")[[1]][1],
                                 false = Accession))
    
    # Rebuild identifier
    temp.data <- temp.data %>%
      unite(col = "Identifier",
            Accession,
            Sites,
            sep = "--")
    
  }
  # Exit
  return(temp.data)
  
}