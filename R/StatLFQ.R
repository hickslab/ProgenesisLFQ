load_simple <- function(name, type = ".csv"){
  if (type == ".csv"){
    data <- read.csv(name, stringsAsFactors = FALSE)
    
  } else if (type == ".txt"){
    data <- read.table(name, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
  }
}


get_design <- function(group){
  num.cond <- letters[1:length(group)] %>%
    toupper()
  
  num.sample = 1:(length(unlist(group)) / length(group))
  
  temp.design <- outer(num.cond, num.sample, paste, sep = "") %>%
    as.vector() %>%
    sort()
  
}


clean_min <- function(data, group = group, nonzero = 3){
  # Create vector to store index
  idx <- c()
  i <- 1
  
  # Iterate by row
  for (x in 1:nrow(data)){
    row <- data[x,]
    
    # Initiate condition
    keep <- FALSE
    
    # Iterate by replicates for each sample
    for (y in group){
      row2 <- row[, y]
      
      # Threshold how many columns can be not equal to 0
      if (sum(row2 != 0) >= nonzero){
        keep <- TRUE
        
      }
    }
    # Check on condition state
    if(keep == TRUE){
      idx[i] <- rownames(row)
      
    }
    i <- i + 1
    
  }
  return(data[which(rownames(data) %in% idx), ])
  
}


transform_data <- function(data, group = group, method = "log2"){
  if (method == "log2"){
    data[, unlist(group)] <- log2(data[, unlist(group)])
    data[, unlist(group)][data[, unlist(group)] == -Inf] <- 0
    
  } else if (method == "asinh"){
    data[, unlist(group)] <- asinh(data[, unlist(group)])
    
  }
  return(data)
  
}

impute_imp4p <- function(data, group = group, method = "log2"){
  # Set RNG state
  set.seed(123)
  
  # Define condition and sample number based on group
  num.cond <- length(group)
  num.sample <- length(unlist(group)) / length(group)
  
  # Create vector showing replicates for each sample
  membership <- gen.cond(num.cond, num.sample)
  
  # Store raw data to manipulate
  temp.data <- data.matrix(data[, unlist(group)])
  
  # All 0 to NA for 'imp4p' requirement
  temp.data[temp.data == 0] <- NA
  
  # Impute rows with with nonzeros in every condition
  temp.data <- impute.rand(temp.data, membership)
  
  # Impute rows with only zeros in at least one condition
  temp.data <- impute.pa(temp.data, membership, q.norm = 0)
  
  # Package 'imp4p' intended for log2-transformed data and may impute negative values
  if (method == "asinh"){
    temp.data$tab.imp <- abs(temp.data$tab.imp)
    
  }
  # Write onto original data
  data[, unlist(group)] <- data.frame(temp.data$tab.imp)
  
  return(data)
  
}


calculate_ttest <- function(data, group.compare = group.compare, fdr = TRUE){
  # @TODO
  
  # Dataframe with only variable column to build upon
  #temp.data <- data %>%
  #  select(1)
  
  #for (x in group.compare){
  #data %>%
  #  select(unlist(x)) %>%
  #  rownames_to_column() %>%
  #  gather(sample, abundance, -rowname) %>%
  #  mutate(cond = str_sub(sample, start = 1, end = 1)) %>%
  #  group_by(rowname) %>%
  #  mutate(p = t.test(abundance, alternative = "two.sided", var.equal = TRUE)$p.value)
  
  # @TODO: allow FDR toggle on/off
  
  
  for (x in group.compare){
    # Row-wise t-test
    p.ttest <- c()
    for (i in 1:nrow(data)){
      row <- data[i, ]
      
      ttest.out <- t.test(row[, x[[1]]],
                          row[, x[[2]]],
                          alternative = "two.sided",
                          var.equal = TRUE)
      
      p.ttest[i] <- ttest.out$p.value
    }
    
    if (fdr == TRUE){
      # Benjamini-Hochberg FDR correction
      p.ttest <- p.adjust(p.ttest, method = "BH", n = length(p.ttest))
      
    }
    # Add to dataframe
    data <- cbind(data, temp.name = p.ttest)
    
    
    # Column name defined by group.compare variable and 'get_design' nomenclature
    temp1 <- data %>%
      select(x[[1]][1]) %>%
      names() %>%
      str_sub(., start = 1, end = 1)
    
    temp2 <- data %>%
      select(x[[2]][1]) %>%
      names() %>%
      str_sub(., start = 1, end = 1)
    
    if (fdr == TRUE){
      temp.name <- paste(temp1, temp2, "_FDR", sep = "")
      
    } else {
      temp.name <- paste(temp1, temp2, "_P", sep = "")
      
    }
    
    # Rename appended column with dynamic variable
    names(data)[names(data) == "temp.name"] <- temp.name
    
  }
  return(data)
  
}


calculate_fc <- function(data, group.compare = group.compare, method = "log2", difference = TRUE){
  # @param: data - column 1 is a variable with abundance columns following, must be transformed
  
  for (x in group.compare){
    temp.data <- data
    
    
    # Specify data scale based on 'method' argument
    if (method == "log2"){
      # Do nothing
      
    } else if (method == "asinh"){
      temp.data[, unlist(x)] <- sinh(temp.data[, unlist(x)])
      
    }
    # Specify calculation based on 'difference' argument
    if (difference == TRUE){
      temp.fc <- rowMeans(temp.data[x[[2]]]) - rowMeans(temp.data[x[[1]]])
      
    } else if (difference == FALSE){
      temp.fc <- rowMeans(temp.data[x[[2]]]) / rowMeans(temp.data[x[[1]]])
      
    }
    # Add to dataframe
    data <- cbind(data, temp.name = temp.fc)
    
    
    # Column name defined by group.compare variable and 'get_design' nomenclature
    temp1 <- data %>%
      select(x[[1]][1]) %>%
      names() %>%
      str_sub(., start = 1, end = 1)
    
    temp2 <- data %>%
      select(x[[2]][1]) %>%
      names() %>%
      str_sub(., start = 1, end = 1)
    
    temp.name <- paste(temp1, temp2, "_FC", sep = "")
    
    
    # Rename appended column with dynamic variable
    names(data)[names(data) == "temp.name"] <- temp.name
    
  }
  return(data)
  
}



