
exchange_gene <- function(data, variable, organism){
  # Details: https://www.uniprot.org/help/api_queries
  # Columns: https://www.uniprot.org/help/uniprotkb_column_names
  
  # Arabidopsis thaliana: UP000006548
  # Chlamydomonas reinhardtii: UP000006906
  # Oryza sativa subsp. japonica: UP000059680
  # Homo sapiens: UP000005640

  
  # Load UniProt protein information table
  url <- paste("https://www.uniprot.org/uniprot/?query=proteome:",
               organism,
               "&columns=id,genes(OLN)&format=tab",
               sep = "")

  # Load UniProt protein information table
  annot <- read.table(url,
                      sep = "\t",
                      header = TRUE,
                      stringsAsFactors = FALSE)

  # Clean UniProt table
  annot <- annot %>%
    rename(Accession = Entry,
           Gene = Gene.names...ordered.locus..) %>%
    separate(Gene,
             into = c("Gene", "Extra"),
             sep = "[/; ]",
             fill = "right",
             extra = "drop") %>%
    select(-Extra) %>%
    mutate(Gene = str_to_upper(Gene))

  if (variable == "Accession"){
    # Clean data accession from UniProt
    temp.data <- data %>%
      mutate(Accession = str_sub(Accession, start = 4)) %>%
      separate(Accession,
               into = c("Accession", "Name"),
               sep = "\\|",
               extra = "merge") %>%
      select(-Name)
    
    # Join UniProt gene onto data
    temp.data <- temp.data %>%
      inner_join(., annot, by = "Accession")
    
    # Exit
    return(temp.data %>%
             select(Gene))
    
  } else if (variable == "Identifier"){
    # Clean data accession from UniProt
    temp.data <- data %>%
      select(Identifier) %>%
      separate(Identifier,
               into = c("Accession", "Sites"),
               sep = "--",
               extra = "drop",
               remove = FALSE) %>%
      mutate(Accession = str_sub(Accession, start = 4)) %>%
      separate(Accession,
               into = c("Accession", "Name"),
               sep = "\\|",
               extra = "merge") %>%
      select(-Name)
    
    # Join UniProt gene onto data
    temp.data <- temp.data %>%
      inner_join(., annot, by = "Accession")
    
    # Rebuild identifier and exit
    return(temp.data %>%
             unite(., col = "Identifier", Gene, Sites, sep = "--") %>%
             select(Identifier))
    
  }
}


remove_PACid <- function(data){
  variable <- data %>%
    select(1) %>%
    names()
  
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


split_identifier <- function(df){
  temp.data <- df %>%
    separate(Identifier,
             into = c("Accession", "Sites"),
             sep = "--",
             remove = FALSE) %>%
    select(-Accession, -Sites, everything())
  
}


add_missingness <- function(data, raw, variable){
  if (variable == "Accession"){
    # Filter raw data for StatLFQ passing observations
    temp.raw <- raw %>%
      semi_join(., data, by = "Accession")
    
    # Convert to long-format and summarize for complete missingness
    temp.raw <- temp.raw %>%
      gather(sample, value, -1) %>%
      separate(sample,
               into = c("condition", "replicate"),
               sep = "[0-9]+") %>%
      select(-replicate) %>%
      group_by(Accession, condition) %>%
      summarize(sum = sum(value) != 0) %>%
      filter(sum == 0)
    
    # Add column with groups completely missing
    temp.raw <- temp.raw %>%
      group_by(Accession) %>%
      mutate(Absent = paste(condition, collapse = "")) %>%
      select(-condition, -sum)
    
    # Join missingness column onto data
    data <- data %>%
      left_join(., temp.raw, by = "Accession")
    
  } else if (variable == "Identifier"){
      # Filter raw data for StatLFQ passing observations
      temp.raw <- raw %>%
        semi_join(., data, by = "Identifier")
      
      # Convert to long-format and summarize for complete missingness
      temp.raw <- temp.raw %>%
        gather(sample, value, -1) %>%
        separate(sample,
                 into = c("condition", "replicate"),
                 sep = "[0-9]+") %>%
        select(-replicate) %>%
        group_by(Identifier, condition) %>%
        summarize(sum = sum(value) != 0) %>%
        filter(sum == 0)
      
      # Add column with groups completely missing
      temp.raw <- temp.raw %>%
        group_by(Identifier) %>%
        mutate(Absent = paste(condition, collapse = "")) %>%
        select(-condition, -sum)
      
      # Join missingness column onto data
      data <- data %>%
        left_join(., temp.raw, by = "Identifier")
    
  }
  # Exit
  return(data)
  
}


add_phytozome <- function(df, path){
  # Load Phytozome annotations
  temp.data <- read_delim(path, delim = "\t", col_types = cols()) %>%
    mutate(Accession = peptideName)
  
  # Join onto input data
  temp.data <- temp.data %>%
    left_join(df, ., by = "Accession")
  
}


