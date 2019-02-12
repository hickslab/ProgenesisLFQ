load_data <- function(name){
  data <- read.csv(name, skip = 2, stringsAsFactors = FALSE)
  
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