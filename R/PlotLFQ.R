plot_corr <- function(data, group = group) {
  # Load packages
  library(broom)
  
  data2 %>%
    select(-1) %>%
    cor() %>%
    tidy() %>%
    gather(replicate, cor, -.rownames) %>%
    ggplot(., aes(x = .rownames, y = replicate, fill = cor)) +
    geom_tile()
  
}


plot_dendrogram <- function(df, group){
  # Load packages
  library(dendextend)
  
  temp.data <- df %>%
    select(group %>% flatten_int()) %>%
    scale() %>%
    t() %>%
    dist() %>%
    hclust() %>%
    #cutree(., 3) %>%
    as.dendrogram()
  
  temp.data %>%
    set("branches_k_color", k = 3) %>%
    plot(xlab = "Distance", ylab = "Replicate")
  
  temp.data %>%
    rect.dendrogram(k = 3)
  
}


plot_pca <- function(df, group){
  temp.data <- df %>%
    select(group %>% flatten_int())
  
  temp.data %>%
    prcomp() %>%
    .$rotation %>%
    data.frame() %>%
    rownames_to_column() %>%
    separate(rowname, into = c("condition", "replicate"), sep = "-") %>%
    ggplot(., aes(x = PC1, y = PC2, fill = condition, label = replicate)) +
    geom_point(shape = 21, size = 8) +
    scale_color_discrete(limits = names(group)) +
    #theme_bw(base_size = 20) +
    theme(text = element_text(size = 20)) +
    geom_text(color = "black")
  
}


plot_volcano <- function(data3, group, group.compare, fdr = TRUE, threshold = 2, xlimit = 10, ylimit = 8){
  # Data preparation
  temp.data <- 	data3 %>%
	  select(-unlist(group)) %>%
	  gather(compare, value, -1) %>%
	  separate(compare,
	           sep = "_",
	           into = c("compare", "variable"),
	           extra = "merge",
	           fill = "right") %>%
	  spread(variable, value)
  
  # Check if FDR-adjustment was applied
  if (fdr == TRUE){
    temp.data <- temp.data %>%
      mutate(significance = FDR)
    
  } else {
    temp.data <- temp.data %>%
      mutate(significance = P)
    
  }
  
  # Set significance types
  temp.data <- temp.data %>%
    mutate(down = if_else(FC <= -log2(threshold) & significance < 0.05, 1, 0),
           up = if_else(FC >= log2(threshold) & significance < 0.05, 1, 0),
           type = if_else(down == 1, "down", if_else(up == 1, "up", "same")))
  
  # Build facet titles
  temp.data <- temp.data %>%
    group_by(compare) %>%
    mutate(down = sum(down), up = sum(up)) %>%
    mutate(compare_count = paste(compare, "\nDown ", sep = "", down, " / Up ", up))
	
  # Set facet order
  temp.data <- temp.data %>%
    ungroup() %>%
    mutate(compare = factor(compare, levels = names(group.compare)),
           type = factor(type, levels = c("same", "down", "up")))
  
  # 
  temp.label <- temp.data %>%
    group_by(compare, compare_count) %>%
    count() %>%
    data.frame()
  
	# Plot
	temp.data %>%
	  ggplot(., aes(x = FC, y = -log10(significance), color = type)) +
	  geom_point(alpha = 1, size = 3) +
	  scale_color_manual(values = c("same" = "grey70", "down" = "blue", "up" = "red")) +
	  coord_cartesian(xlim = c(-xlimit, xlimit), ylim = c(0, ylimit)) +
	  xlab(expression("log"[2]*"(fold change)")) +
	  ylab(if_else(fdr == TRUE,
	               expression("-log"[10]*"(FDR-adjusted "*italic(p)*"-value)"),
	               expression( "-log"[10]*"("*italic(p)*"-value)"))) +
	  facet_wrap(~ compare) +
	  geom_text(data = temp.label, aes(x = 0, y = Inf, label = compare_count), inherit.aes = FALSE) +
	  #theme_bw(base_size = 16) +
	  theme(text = element_text(size = 20)) +
	  guides(color = FALSE)

}


plot_hclust <- function(df, group, k = 5){
  # Select abundance columns
  temp.data <- df %>%
    select(1, group %>% flatten_int())
  
  temp.data <- temp.data %>%
    gather(replicate, abundance, -1) %>%
    separate(replicate, into = c("condition", "replicate"), sep = "-") %>%
    group_by(.[[1]], condition) %>%
    summarize(mean = mean(abundance)) %>%
    ungroup() %>%
    spread(., condition, mean)
  
    
  # Z-score rowwise normalization
  temp.data[-1] <- temp.data %>%
    select(-1) %>%
    apply(., 1, scale) %>%
    t()
  
  temp.data <- temp.data %>%
    #select(unlist(group)) %>%
    select(-1) %>%
    
    mutate(clustered = dist(.) %>%
             hclust() %>%
             cutree(k)) %>%
    
    mutate(clustered = LETTERS[clustered]) %>%
    
    group_by(clustered) %>%
    mutate(clustered_count = paste(clustered, "\n", sep = "", n())) %>%
    ungroup() %>%
    
    rownames_to_column() %>%
    
    gather(condition, scaled, -rowname, -clustered, -clustered_count)
  
  # Set condition order
  temp.data <- temp.data %>%
    mutate(condition = factor(condition, level = names(group)))
  
  # Plot
  temp.data %>%
    #ggplot(., aes(x = condition, y = scaled, fill = factor(clustered))) + geom_boxplot(outlier.shape = NA) + guides(fill = FALSE) + scale_x_discrete(limits = names(group)) + 
    ggplot(., aes(x = condition, y = scaled, group = rowname, color = factor(clustered))) + geom_line() + guides(color = FALSE) + #scale_color_brewer(palette = "Pastel1") +
    #ggplot(., aes(x = condition, y = scaled)) + geom_line(group = 2) +# geom_smooth(method = "loess") +
    
    facet_wrap(~ clustered_count) +
    
    labs(x = "Condition", y = "Z-score") +
    
    theme_bw(base_size = 20)
  
}
 

plot_heatmap <- function(df){

  temp.data %>%
    ggplot(., aes(variable, reorder(Identifier, -value))) +
    #ggplot(., aes(variable, Identifier)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient(low = "white", high = "red") +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    xlab("")
  
}


plot_GO_uniprot <- function(df, top = 5){
  # Select GO columns
  df %>%
    #filter(fdr < 0.05) %>%
    #sample_n(., 1000, replace = FALSE) %>%
    select(contains("Gene ontology")) %>%
    
    # Melt GO
    gather(column, term) %>%
    mutate(column = str_split(column, "\\(", simplify = TRUE)[, 2] %>% str_sub(., end = -2)) %>%
    
    # Clean terms
    filter(!is.na(term)) %>%
    separate_rows(term, sep = "; ") %>%
    #separate(term, into = c("term", "GO"), sep = "\\[", extra = "merge") %>% View
    mutate(term = str_sub(term, end = -14)) %>%
    #mutate(GO = str_sub(GO, end = -2)) %>%
    
    # Count terms
    group_by(column, term) %>%
    tally() %>%
    
    # Filter for top terms
    group_by(column) %>%
    top_n(top, n) %>%
    #slice(1:n())
    
    #ggplot(., aes(x = count, fill = column)) + geom_histogram(binwidth = 1) + coord_cartesian(xlim = c(0, 100))
    
    #ggplot(., aes(x = reorder(term, count), y = count, color = column)) + geom_point() + coord_flip()
    
    ggplot(., aes(x = reorder(term, n), y = n, fill = column)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_grid(column ~ ., scales = "free") +
    guides(fill = FALSE) +
    labs(x = NULL, y = "Proteins") +
    theme(text = element_text(size = 16))
  
}


plot_box <- function(data, group = group){
  temp.data <- data %>%
    select(-Identifier) %>%
    #apply(., 1, scale) %>%
    #t() %>%
    #data.frame() %>%
    log(., 2) %>%
    #mutate_all(funs(replace(., is.na(.), 0))) %>%
    mutate_all(funs(replace(., . == -Inf, 0))) %>%
    bind_cols(data.frame(Identifier = data$Identifier), .) %>%
    melt() %>%
    #melt(na.rm = TRUE) %>%
    mutate(group = str_sub(variable, end = -3))
  
  
  # Plot
  temp.data %>%
    ggplot(., aes(variable, value, fill = group)) +
    geom_boxplot(outlier.alpha = 0.1) +
    facet_wrap(~group, scales = "free_x") +
    scale_fill_discrete()
  
}


plot_violin <- function(data, group = group){
  temp.data <- data %>%
    select(-Identifier) %>%
    #apply(., 1, scale) %>%
    #t() %>%
    #data.frame() %>%
    log(., 2) %>%
    #mutate_all(funs(replace(., is.na(.), 0))) %>%
    mutate_all(funs(replace(., . == -Inf, 0))) %>%
    bind_cols(data.frame(Identifier = data$Identifier), .) %>%
    melt() %>%
    #melt(na.rm = TRUE) %>%
    mutate(group = str_sub(variable, end = -3))
  
  
  # Plot
  temp.data %>%
    ggplot(., aes(variable, value, fill = group)) +
    geom_violin(trim = TRUE) +
    #geom_jitter(alpha = 0.1) +
    scale_fill_discrete()
  
}


plot_joy <- function(data, group = group){
  temp.data <- data %>%
    select(-Identifier) %>%
    #apply(., 1, scale) %>%
    #t() %>%
    #data.frame() %>%
    log(., 2) %>%
    #mutate_all(funs(replace(., is.na(.), 0))) %>%
    mutate_all(funs(replace(., . == -Inf, 0))) %>%
    bind_cols(data.frame(Identifier = data$Identifier), .) %>%
    melt() %>%
    #melt(na.rm = TRUE) %>%
    mutate(group = str_sub(variable, end = -3))
  
  
  # Plot
  temp.data %>%
    ggplot(., aes(value, variable, height = ..density.., fill = group)) +
    geom_joy(stat = "density", scale = 3) +
    #theme_joy(base_size = 14) +
    theme(axis.text = element_text(vjust = 0, color = "black")) +
    scale_x_continuous(expand = c(0.01, 0)) +
    scale_y_discrete(expand = c(0.01, 0)) +
    xlab(expression("log"[2]*"(Abundance)")) +
    #facet_wrap(~group, scales = "free_y") +
    ylab("")
  
}
