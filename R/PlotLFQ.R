plot_corr <- function(data, group = group) {
  # load required packages
  require(GGally, quietly = TRUE)
  
  temp.data <- data %>%
    select(unlist(group))
  
  ggcorr(temp.data,
         label = TRUE)
  
}


plot_volcano <- function(data3, group = group, fdr = TRUE, threshold = 2, xlimit = 10, ylimit = 6){
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
      rename(significance = FDR)
    
  } else {
    temp.data <- temp.data %>%
      rename(significance = P)
    
  }
  
  # Set significance types
  temp.data <- temp.data %>%
    mutate(type = if_else((FC >= log2(threshold) & significance < 0.05),
                                  "UP",
                                  if_else((FC <= -log2(threshold) & significance < 0.05),
                                          "DOWN",
                                          "SAME")))
  
  # Build descriptive plot title
  temp.data <- temp.data %>%
    filter(type != "SAME") %>%
    group_by(compare, type) %>%
    summarize(n = n()) %>%
    spread(type, n) %>%
    mutate(compare_count = paste(compare, "\nDown ", sep = "", DOWN, " / Up ", UP)) %>%
    select(compare, compare_count) %>%
    inner_join(temp.data, ., by = "compare")
	
	# Plot
	temp.data %>%
	  ggplot(., aes(x = FC, y = -log10(significance), color = type)) +
	  geom_point() +
	  scale_color_manual(values = c("SAME" = "grey70", "DOWN" = "blue", "UP" = "red")) +
	  xlim(-xlimit, xlimit) +
	  ylim(0, ylimit) +
	  xlab(expression("log"[2]*"(fold change)")) +
	  ylab(if_else(fdr == TRUE,
	               expression("-log"[10]*"(FDR-adjusted "*italic(p)*"-value)"),
	               expression( "-log"[10]*"("*italic(p)*"-value)"))) +
	  facet_wrap(~ compare_count) +
	  theme_bw(base_size = 16) +
	  guides(color = FALSE)

}


plot_heatmap <- function(data, group = group){
  temp.data <- data
  
  temp.data[unlist(group)] <- temp.data[unlist(group)] %>%
    apply(., 1, scale) %>%
    t()
  
  temp.data <- data %>%
    select(unlist(group)) %>%
    mutate(clustered = dist(.) %>%
             hclust() %>%
             cutree(4)) %>%
    
    mutate(clustered = LETTERS[clustered]) %>%
    
    group_by(clustered) %>%
    mutate(clustered_count = paste(clustered, "\n", sep = "", n())) %>%
    ungroup() %>%
    
    rownames_to_column() %>%
    
    gather(cond, scaled, unlist(group))
  
  
  # Plot
  temp.data %>%
    #ggplot(., aes(x = cond, y = scaled)) + geom_boxplot(outlier.shape = NA) +
    ggplot(., aes(x = cond, y = scaled, group = rowname, color = factor(clustered))) + geom_line() +
    
    facet_wrap(~ clustered_count) +
    
    labs(x = "Condition", y = "Z-score") +
    
    #scale_x_discrete(labels = c("TAP", "-N", "-C", "-P")) +
    
    guides(color = FALSE) +
    
    theme_bw(base_size = 16)
  
  
  
  temp.data %>%
    ggplot(., aes(variable, reorder(Identifier, -value))) +
    #ggplot(., aes(variable, Identifier)) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient(low = "white", high = "red") +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    xlab("")
  
}


plot_pca <- function(data2){
  data2 %>%
    select(-1) %>%
    prcomp() %>%
    .$rotation %>%
    data.frame() %>%
    rownames_to_column() %>%
    mutate(condition = str_sub(rowname, end = -2)) %>%
    ggplot(., aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 5)
  
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


