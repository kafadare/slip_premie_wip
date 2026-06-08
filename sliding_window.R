sliding_window_eqdist <- function(df, 
                           phenotype, 
                           pred_var, 
                           age_var,
                           window_size, 
                           overlap,
                           covars = NULL,
                           group_var = NA, 
                           dist_plot = FALSE) {
  #parameters
  #window_size <- window_size # years
  #overlap <- overlap # years
  step <- window_size - overlap
  
  age_min <- floor(min(df[[age_var]], na.rm = TRUE))
  age_max <- ceiling(max(df[[age_var]], na.rm = TRUE))
  
  #sliding window
  results <- data.frame(start_age = numeric(),
                        end_age = numeric(),
                        beta = numeric(),
                        pval = numeric(), 
                        N = numeric(),
                        group = character())
  
  all_subsets <- list()  # store subsets for optional plotting
  
  if (!is.na(group_var)){
    for (group_level in unique(df[[group_var]])) {
      df_group <- df %>% filter(get(group_var) == group_level)
      
      for (start_age in seq(age_min, age_max - window_size, by = step)) {
        end_age <- start_age + window_size
        
        subset <- df_group %>%
          filter(get(age_var) >= start_age & get(age_var) < end_age) %>%
          filter(!is.na(get(pred_var)))
        
        if (nrow(subset) > 2) {
          if(is.null(covars)) {
            formula = reformulate(pred_var, response = phenotype)
          } else {
            full_pred <- c(pred_var, covars)
            formula = reformulate(full_pred, response = phenotype)
          }
          fit <- lm(formula, data = subset)
          coef_summary <- summary(fit)$coefficients
          beta <- coef_summary[pred_var, "Estimate"]
          pval <- coef_summary[pred_var, "Pr(>|t|)"]
          N <- nrow(subset)
          results <- rbind(results, data.frame(start_age, end_age, beta, pval, N, group = group_level))
          if (dist_plot) {
            subset$window <- paste0(start_age, "-", end_age)
            subset$group  <- group_level
            all_subsets[[length(all_subsets) + 1]] <- subset
          }
        }
      }
    }
  } else {
    for (start_age in seq(age_min, age_max - window_size, by = step)) {
      end_age <- start_age + window_size
      
      subset <- df %>%
        filter(get(age_var) >= start_age & get(age_var) < end_age) %>%
        filter(!is.na(get(pred_var)))
      
      if (nrow(subset) > 2) {
        if(is.null(covars)) {
          formula = reformulate(pred_var, response = phenotype)
        } else {
          full_pred <- c(pred_var, covars)
          formula = reformulate(full_pred, response = phenotype)
        }
        fit <- lm(formula, data = subset)
        coef_summary <- summary(fit)$coefficients
        beta <- coef_summary[pred_var, "Estimate"]
        pval <- coef_summary[pred_var, "Pr(>|t|)"]
        N <- nrow(subset)
        results <- rbind(results, data.frame(start_age, end_age, beta, pval, N))
        if (dist_plot) {
          subset$window <- paste0(start_age, "-", end_age)
          all_subsets[[length(all_subsets) + 1]] <- subset
        }
      }
    }
  }
  
  # add midpoint and significance flag
  results <- results %>%
    mutate(mid_age = (start_age + end_age) / 2,
           sig = ifelse(pval < 0.05, "Significant", "Not Significant"))
  
    if (dist_plot && length(all_subsets) > 0) {
      df_plot <- dplyr::bind_rows(all_subsets)
      
      # order facet levels by sliding window order of creation
      df_plot$window <- factor(df_plot$window, 
                               levels = unique(df_plot$window))
      
      p <- ggplot(df_plot, aes(x = .data[[pred_var]])) +
        geom_histogram(bins = 20, fill = "steelblue", color = "white") +
        facet_wrap(~window, scales = "free_y") +
        theme_minimal() +
        labs(title = paste("Distribution of", pred_var, "across age windows"),
             x = pred_var, y = "Count")
      
      if (!is.na(group_var)) {
        p <- p + facet_wrap(vars(window, group), scales = "free_y")
      }
      
      print(p)
    }
  
  return(results)
}

sliding_window_eqsize <- function(df, 
                                  phenotype, 
                                  pred_var, 
                                  age_var,
                                  window_size, 
                                  overlap,
                                  covars = NULL,
                                  group_var = NA, 
                                  dist_plot = FALSE) {
  
  step <- window_size - overlap
  
  results <- data.frame(start_age = numeric(),
                        end_age = numeric(),
                        beta = numeric(),
                        pval = numeric(), 
                        N = numeric(),
                        group = character())
  
  all_subsets <- list()
  
  if (!is.na(group_var)){
    for (group_level in unique(df[[group_var]])) {
      
      df_group <- df %>%
        filter(get(group_var) == group_level) %>%
          filter(!is.na(get(pred_var))) %>% 
        arrange(get(age_var))
      
      n_total <- nrow(df_group)
      
      for (start_idx in seq(1, n_total - window_size + 1, by = step)) {
        end_idx <- start_idx + window_size - 1
        
        subset <- df_group[start_idx:end_idx, ]
        
        if (nrow(subset) > 2) {
          
          start_age <- min(subset[[age_var]], na.rm = TRUE)
          end_age   <- max(subset[[age_var]], na.rm = TRUE)
          
          if(is.null(covars)) {
            formula = reformulate(pred_var, response = phenotype)
          } else {
            full_pred <- c(pred_var, covars)
            formula = reformulate(full_pred, response = phenotype)
          }
          
          fit <- lm(formula, data = subset)
          coef_summary <- summary(fit)$coefficients
          
          beta <- coef_summary[pred_var, "Estimate"]
          pval <- coef_summary[pred_var, "Pr(>|t|)"]
          N <- nrow(subset)
          
          results <- rbind(results, data.frame(start_age, end_age, beta, pval, N, group = group_level))
          
          if (dist_plot) {
            subset$window <- paste0(start_age, "-", end_age)
            subset$group  <- group_level
            all_subsets[[length(all_subsets) + 1]] <- subset
          }
        }
      }
    }
    
  } else {
    
    df <- df %>%
      filter(!is.na(get(pred_var))) %>% 
      arrange(get(age_var))
    
    n_total <- nrow(df)
    
    for (start_idx in seq(1, n_total - window_size + 1, by = step)) {
      end_idx <- start_idx + window_size - 1
      
      subset <- df[start_idx:end_idx, ]
      
      if (nrow(subset) > 2) {
        
        start_age <- min(subset[[age_var]], na.rm = TRUE)
        end_age   <- max(subset[[age_var]], na.rm = TRUE)
        
        if(is.null(covars)) {
          formula = reformulate(pred_var, response = phenotype)
        } else {
          full_pred <- c(pred_var, covars)
          formula = reformulate(full_pred, response = phenotype)
        }
        
        fit <- lm(formula, data = subset)
        coef_summary <- summary(fit)$coefficients
        
        beta <- coef_summary[pred_var, "Estimate"]
        pval <- coef_summary[pred_var, "Pr(>|t|)"]
        N <- nrow(subset)
        
        results <- rbind(results, data.frame(start_age, end_age, beta, pval, N))
        
        if (dist_plot) {
          subset$window <- paste0(start_age, "-", end_age)
          all_subsets[[length(all_subsets) + 1]] <- subset
        }
      }
    }
  }
  
  results <- results %>%
    mutate(mid_age = (start_age + end_age) / 2,
           sig = ifelse(pval < 0.05, "Significant", "Not Significant"))
  
  if (dist_plot && length(all_subsets) > 0) {
    df_plot <- dplyr::bind_rows(all_subsets)
    
    df_plot$window <- factor(df_plot$window, 
                             levels = unique(df_plot$window))
    
    p <- ggplot(df_plot, aes(x = .data[[pred_var]])) +
      geom_histogram(bins = 20, fill = "steelblue", color = "white") +
      facet_wrap(~window, scales = "free_y") +
      theme_minimal() +
      labs(title = paste("Distribution of", pred_var, "across age windows"),
           x = pred_var, y = "Count")
    
    if (!is.na(group_var)) {
      p <- p + facet_wrap(vars(window, group), scales = "free_y")
    }
    
    print(p)
  }
  
  return(results)
}