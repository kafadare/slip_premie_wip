# Setup ----

library(tidyverse)
library(glue)
library(ggseg)
library(gridExtra)
library(ggsegTissue)
library(kableExtra)
library(table1)
library(flextable)
library(broom)
library(broom.mixed)
library(effectsize)
library(ggpubr)
library(cowplot)
library(reticulate)

rm(list=ls())
set.seed(2514)

# Define default color scheme
color_scheme <- c("royalblue","white","firebrick")

# Define global measures
global_measures <- c("eTIV", "GMV", "sGMV", "WMV", "Ventricles", "Cerebellum.all", 
                     "totalSA", "meanCT")

global_measure_names <- c("ICV", "GMV", "sGMV", "WMV", "Ventricles", "Cerebellum", 
                          "SA", "CT")


# Setup Atlas
# Some labels in each atlas are NA, replace with the ROI
dk_plot <- as_tibble(ggseg::dk)
dk_plot$label[is.na(dk_plot$label)] <- dk_plot$roi[is.na(dk_plot$label)]

aseg_plot <- as_tibble(ggseg::aseg) %>%
  mutate(label = str_replace_all(label, "-","_"))
aseg_plot$roi <- sprintf("%04d", as.numeric(rownames(aseg_plot)))
aseg_plot$label[is.na(aseg_plot$label)] <- aseg_plot$roi[is.na(aseg_plot$label)]


# Table Functions ---------------------
pvalue <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    if (length(unique(g)) == 2) {
      # For numeric variables, perform a standard 2-sample t-test
      p <- t.test(y ~ g)$p.value
      
    } else {
      # For numeric variables, perform a standard ANOVA
      p <- summary(aov(y ~ g))[[1]]["g","Pr(>F)"]
    }
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  # c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
  c("", format.pval(p, digits=3, eps=0.001))
}

pvalue_1_3 <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  y <- y[g %in% c(1,3)]
  g <- g[g %in% c(1,3)]
  g <- factor(g)
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  # c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
  c("", format.pval(p, digits=3, eps=0.001))
}

pvalue_1_2 <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  y <- y[g %in% c(1,2)]
  g <- g[g %in% c(1,2)]
  g <- factor(g)
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  # c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
  c("", format.pval(p, digits=3, eps=0.001))
}
pvalue_3_4 <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  y <- y[g %in% c(3,4)]
  g <- g[g %in% c(3,4)]
  g <- factor(g)
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  # c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
  c("", format.pval(p, digits=3, eps=0.001))
}

pvalue_2_4 <- function(x, ...) {
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  y <- y[g %in% c(2,4)]
  g <- g[g %in% c(2,4)]
  g <- factor(g)
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  # c("", sub("<", "&lt;", format.pval(p, digits=3, eps=0.001)))
  c("", format.pval(p, digits=3, eps=0.001))
}

# Primary Case-Control Model ---------------------------------
run_model <- function(data, cov_str, features, predictor = "dx", group_control = "CN", group_case = "22q11DS",rand_str = "", model_type = "centile",
                      mpr_only_cols = NA, simple = F) {
  
  # Define independent predictor
  data$predictor <- factor(data[[predictor]], 
                           levels = c(group_control,group_case), 
                           labels = c(FALSE, TRUE))
  # Drop NA predictors
  data <- data[!is.na(data$predictor),]
  # prop_thresh_low <- 0.05
  # prop_thresh_high <- 1 - prop_thresh_low
  # if (model_type == "zscore") {
  #   prop_thresh_low <- qnorm(prop_thresh_low)
  #   prop_thresh_high <- qnorm(prop_thresh_high)
  # }
  
  for (feature in features) {
    # dplyr::select centile feature
    if (model_type == "centile") {
      measure <- glue("{feature}Transformed.q.wre")
    } else if (model_type == "zscore") {
      measure <- glue("{feature}Transformed.q.zscore")
    } else if (model_type == "volume") {
      measure <- glue("{feature}Transformed.normalised")
    }
    # print(measure)
    if (measure %in% colnames(data)) {
      if (feature %in% mpr_only_cols) {
        cov_str_tmp <- str_replace_all(cov_str,"euler","euler_mpr")
      } else {
        cov_str_tmp <- cov_str
      }
      # Define final model, removing ICV if it is the dependent var
      form <- as.formula(glue("{measure} ~ predictor{cov_str_tmp}{rand_str}"))
      # print(form)
      data_mod <- data[!is.na(data[[measure]]) & !is.infinite(data[[measure]]),]
      #Skip any model that doesn't have any data in one of the levels
      if (any(table(data_mod$predictor) == 0)) {
        next
      }
      
      # Run lm or lmer
      if (rand_str == "") {
        res <- lm(form, data = data_mod)
      } else {
        res <- lmer(form, data = data_mod)
      }
      
      # Collect results
      if (simple) {
        res_tidy <- tidy(res) %>%
          filter(term == "predictorTRUE") %>%
          mutate(label = feature)
      } else {
        res_tidy <- tidy(res) %>%
          mutate(label = feature,
                 Predictor = glue("{predictor}_{group_control}_vs_{group_case}"),
                 Case = sum(data_mod$predictor == TRUE),
                 Control = sum(data_mod$predictor == FALSE))
        
        if (rand_str != "") {
          res_tidy$isSingular <- isSingular(res)
        }
      }
      
      #Cohen's D
      d <- as_tibble(t_to_d(res_tidy$statistic,df.residual(res))) %>%
        dplyr::select(Cohens_d = d,
                      Cohens_d_ci_lo = CI_low,
                      Cohens_d_ci_hi = CI_high)
      
      res_tidy <- bind_cols(res_tidy, d)
      
      # res_tidy$BIC <- BIC(res)
      
      if (feature == features[1]) {
        res_all <- res_tidy
      } else {
        res_all <- rbind(res_all,res_tidy)
      }
    }
  }
  
  
  # Final ordering
  res_all <- res_all %>% 
    relocate(p.value, estimate, std.error, statistic, 
             .after = label) %>%
    relocate(Cohens_d,Cohens_d_ci_lo, Cohens_d_ci_hi, .before = p.value)
  
  return(res_all)
}


# Create wrapper function to run linear regression with interactive term for a continous predictor
run_model_interactive_cont <- function(data, cov_str, features, rand_str = "",
                                       predictor = "dx", group_control = "CN", 
                                       group_case = "22q11DS", interactor = "age_days",  
                                       model_type = "centile", gam = F,
                                       mpr_only_cols = NA) {
  
  # Define independent predictor
  data$predictor <- factor(data[[predictor]], 
                           levels = c(group_control,group_case), 
                           labels = c(group_control, group_case))
  
  data$interactor <- data[[interactor]]
  
  # Drop NA predictors
  data <- data[!is.na(data$predictor) & !is.na(data$interactor),]
  
  prop_thresh_low <- 0.05
  prop_thresh_high <- 1 - prop_thresh_low
  
  if (model_type == "zscore") {
    prop_thresh_low <- qnorm(prop_thresh_low)
    prop_thresh_high <- qnorm(prop_thresh_high)
  }
  
  for (feature in features) {
    # dplyr::select centile feature
    if (model_type == "centile") {
      measure <- glue("{feature}Transformed.q.wre")
    } else if (model_type == "zscore") {
      measure <- glue("{feature}Transformed.q.zscore")
    } else if (model_type == "volume") {
      measure <- glue("{feature}Transformed.normalised")
    }
    
    if (feature %in% mpr_only_cols) {
      cov_str_tmp <- str_replace_all(cov_str,"euler","euler_mpr")
    } else {
      cov_str_tmp <- cov_str
    }
    
    # Define model predictors + covs + rand effects
    model_str <- glue("~ predictor*interactor{cov_str_tmp}{rand_str}")
    
    # print(measure)
    if (measure %in% colnames(data)) {
      form <- as.formula(glue("{measure}{model_str}"))
      print(form)
      data_mod <- data[!is.na(data[[measure]]),]
      #Skip any model that doesn't have any data in one of the levels
      if (any(table(data_mod$predictor) == 0)) {
        next
      }
      # Run lm or lmer
      if (rand_str == "") {
        res <- lm(form, data = data_mod)
      } else if (gam) {
        res <- gam(form, data = data_mod)
      }else {
        res <- lmer(form, data = data_mod)
      }
      
      # Collect results
      res_tidy <- tidy(res)
      
      #Cohen's D
      d <- as_tibble(t_to_d(res_tidy$statistic,res$df.residual)) %>%
        dplyr::select(Cohens_d = d,
                      Cohens_d_ci_lo = CI_low,
                      Cohens_d_ci_hi = CI_high)
      
      res_tidy <- bind_cols(res_tidy, d)
      res_tidy <- res_tidy %>% relocate(estimate, std.error, statistic, p.value, .after = Cohens_d_ci_hi)
      
      res_tidy$BIC <- BIC(res)
      res_tidy$label <- feature
      res_tidy <- res_tidy %>% relocate(label, .before = term)
      
      if (feature == features[1]) {
        res_all <- res_tidy
      } else {
        res_all <- rbind(res_all,res_tidy)
      }
    }
  }
  return(res_all)
}


summarize_results <- function(res, term) {
  res_final <- res[res$term == term,] %>%
    mutate(p.fdr = p.adjust(p.value, method = "fdr"),
           sig = factor(p.fdr <= 0.05, 
                        levels = c(FALSE, TRUE), 
                        labels = c("n.s.","p.fdr < 0.05")),
           sig_label = case_when(
             p.fdr <= 0.001 ~ "***",
             p.fdr <= 0.01 ~ "**",
             p.fdr <= 0.05 ~ "*",
             .default = NA)) %>%
    relocate(p.fdr, sig, .after = p.value)
  return(res_final)
}

fmt_table <- function(res) {
  res_fmt <- res %>% 
    select(label, contains("Cohens_d"),statistic, std.error,
           starts_with("p."), Case, Control) %>%
    mutate(label = case_match(label,
                              "Cerebellum.all" ~ "Cerebellum",
                              "SUBC.Left_Accumbens.area" ~ "SUBC.Left_Accumbens",
                              "SUBC.Right_Accumbens.area" ~ "SUBC.Right_Accumbens",
                              .default = label)) %>%
    separate_wider_delim(label, delim = ".", names = c("IDP Group", "IDP"),
                         too_few = "align_end") %>%
    mutate(`IDP Group` = ifelse(is.na(`IDP Group`), "Global", `IDP Group`),
           `IDP Group` = factor(`IDP Group`, levels = c("Global","SUBC","GM","SA","CT")),
           across(contains("Cohens"), ~ round(.x, 2)),
           across(starts_with("p."), ~ signif(.x, digits = 3)),
           statistic = round(statistic, 2),
           std.error = round(std.error, 2),
           `Cohen's d (95% CI)` = glue("[{Cohens_d_ci_lo} to {Cohens_d_ci_hi}]")) %>%
    select(-c(Cohens_d_ci_lo,Cohens_d_ci_hi)) %>%
    relocate(`Cohen's d (95% CI)`, .after = Cohens_d) %>%
    rename(`Cohen's d` = Cohens_d,
           Statistic = statistic,
           `Std Error` = std.error,
           p = p.value,
           `p (FDR)` = p.fdr,
           `N case` = Case,
           `N control` = Control) %>%
    arrange(`IDP Group`, IDP)
  return(res_fmt)
}

# Correlation Functions ---------
freedman_lane_perm <- function(df_perm, perm_features, n_perm = 1000, 
                               n_core = 32, suf = "Transformed.q.zscore",
                               covs = "sex + age_days + euler", mixed = F) {
  
  df_perm <- df_perm %>% select(participant, sex, dx, age_days, age_days_sq, euler, 
                                any_of(glue("{perm_features}{suf}")),
                                contains("eTIV"), site) %>%
    na.omit()
  
  df_perm[sapply(df_perm, is.infinite)] <- NA
  perms <- replicate(n_perm, sample(1:nrow(df_perm)))
  out_perms <- as_tibble(matrix(NA, nrow = length(perm_features),
                                ncol = n_perm + 1), .name_repair = "minimal")
  colnames(out_perms) <- c("label",glue("perm_{1:n_perm}"))
  out_perms$label <- perm_features
  cluster <- makeCluster(n_cores)
  registerDoParallel(cluster)
  
  for (n in 1:nrow(out_perms)) {
    results <- list()
    feature <- perm_features[n]
    # print(feature)
    # Fit the reduced model
    if (mixed) {
      form <- as.formula(glue("{feature}{suf} ~ {covs} + (1 | site)"))
      reduced_model <- lmer(form, df_perm)
    } else {
      form <- as.formula(glue("{feature}{suf} ~ {covs}"))
      reduced_model <- lm(form, df_perm)
    }
    
    
    # Permute the residuals and combine with fitted effects
    results <- foreach(i = 1:ncol(perms), 
                       .packages = c("tidyverse","glue","broom","effectsize"),
                       .export = c("feature", "suf","covs")) %dorng% {
                         # Calculate permuted outcome
                         permuted_resid <- residuals(reduced_model)[perms[,i]] + fitted.values(reduced_model)
                         df_perm$perm_outcome <- permuted_resid
                         
                         # Refit model to permuted outcome
                         if (mixed) {
                           form <- as.formula(glue("perm_outcome ~ dx + {covs} + (1 | site)"))
                           permuted_model <- lmer(form, df_perm)
                         } else {
                           form <- as.formula(glue("perm_outcome ~ dx + {covs}"))
                           permuted_model <- lm(form, df_perm)
                         }
                         
                         # Calculate permuted Cohen's d
                         results[i] <- t_to_d(coef(summary(permuted_model))["dx22q11DS","Pr(>|t|)"], 
                                              df.residual(permuted_model))$d
                         
                         # results[i] <- t_to_d(summary(permuted_model)$coefficients["dx22q11DS","t value"], 
                         #                      df.residual(permuted_model))$d
                       }
    out_perms[n,glue("perm_{1:n_perm}")] <- as.list(unlist(results))
  }
  
  stopCluster(cluster)
  return(out_perms)
}

fl_corr_test <- function(true_brain_map, perms, comp_brain_map, 
                         term_1 = "Cohens_d", term_2 = "Cohens_d",
                         cats = c("GM","SA","CT"), method = "pearson",
                         direction = "either") {
  out_res <- tibble(measure = cats, corr = NA, pval = NA)
  for (cat in cats) {
    # Extract sub-measures for the true map
    true_map_cat <- true_brain_map %>% filter(grepl(glue("{cat}\\."),label))
    
    # Extract sub-measures for the comparison map
    comp_map_cat <- comp_brain_map %>% 
      filter(grepl(glue("{cat}\\."),label)) %>%
      arrange(match(label, true_map_cat$label))
    
    # Extract sub-measures for the null maps
    null_map_cat <- perms %>% 
      filter(grepl(glue("{cat}\\."),label)) %>%
      arrange(match(label, true_map_cat$label)) %>%
      select(-label)
    
    # Calculate the true correlation
    true_corr <- cor(true_map_cat %>% pull(term_1),
                     comp_map_cat %>% pull(term_2),
                     method = method)
    
    # Calculate the null correlations
    null_corrs <- cor(null_map_cat,
                      comp_map_cat %>% pull(term_2),
                      method = method)
    
    # if (sign(true_corr) < 0) {
    #   pval <- sum(null_corrs <= true_corr)/nrow(null_corrs)
    # } else if (sign(true_corr) > 0) {
    #   pval <- sum(null_corrs >= true_corr)/nrow(null_corrs)
    # }
    if (direction == "either") {
      pval <- sum(abs(null_corrs) >= abs(true_corr))/nrow(null_corrs)
    } else if (direction == "greater") {
      pval <- sum(null_corrs >= true_corr)/nrow(null_corrs)
    }  else if (direction == "lesser") {
      pval <- sum(null_corrs <= true_corr)/nrow(null_corrs)
    } 
    
    # Save results
    out_res[out_res$measure == cat,"corr"] <- true_corr 
    out_res[out_res$measure == cat,"pval"] <- pval 
  }
  return(out_res)
}

# Plotting functions -------------------

plot_atlas <- function(data, atlas, fill_col = "Cohens_d",  color_scheme = c("royalblue","white","firebrick"), 
                       plot_lims = c(-10,1), colorbar_title = "",sig_col = NA, p_col = "p.value", 
                       position = "identity", lh_only = F, skip_sagittal = T, percent = F) {
  
  #If seeking to plot results by significance, columns of atlas need to be ordered 
  # such that the significant borders take precedence over the n.s. border
  atlas <- merge(atlas, data, by = "label", all.x = T)
  if (!is.na(sig_col)) {
    atlas <- atlas[order(!is.na(atlas[[p_col]]),atlas[[sig_col]], -atlas[[p_col]]),]
    # print(atlas$p.value)
    if (sum(atlas[[sig_col]] == "n.s.", na.rm = T) == 0) {
      fill_colors <- c("black")
    } else {
      fill_colors <- c("grey","black")
    }
  } else {
    atlas$sig_tmp <- FALSE
    sig_col <- "sig_tmp"
    fill_colors <- c("black","grey")
  }
  if (skip_sagittal) {
    atlas <- atlas %>% filter(side != "sagittal")
  }
  
  #Generate ggplot
  if (lh_only) {
    plot_out <- ggplot() + 
      geom_brain(atlas=atlas, size = 0.5, hemi = "left",
                 position = position_brain(hemi + side ~ .),
                 mapping=aes(fill  = get(fill_col),
                             color = get(sig_col)))
  } else if (position == "identity") {
    plot_out <- ggplot() + 
      geom_brain(atlas=atlas, size = 0.5, 
                 mapping=aes(fill  = get(fill_col),
                             color = get(sig_col)))
  } else if (position == "stacked") {
    plot_out <- ggplot() + 
      geom_brain(atlas=atlas, size = 0.5, 
                 position = position_brain(hemi ~ side),
                 mapping=aes(fill  = get(fill_col),
                             color = get(sig_col)))
  } else if (position == "stacked2") {
    plot_out <- ggplot() + 
      geom_brain(atlas=atlas, size = 0.5, 
                 position = position_brain(side ~ hemi),
                 mapping=aes(fill  = get(fill_col),
                             color = get(sig_col)))
  }
  plot_out <- plot_out +
    scale_color_manual(, values = fill_colors, na.value="grey", 
                       guide = "none",
                       na.translate=FALSE) +
    theme_brain2(text.colour = "black", text.family = "")  + 
    guides(fill=guide_colourbar(title=colorbar_title,
                                order = 1),
           color = guide_legend(title = "",
                                order = 2,
                                override.aes = list(fill = NA)))
  if (percent) {
    plot_out <- plot_out +
      scale_fill_gradientn(colours = color_scheme, na.value="grey",
                           limits = plot_lims, breaks = c(plot_lims[1],0,plot_lims[2]),
                           labels = scales::label_percent())
  } else {
    plot_out <- plot_out +
      scale_fill_gradientn(colours = color_scheme, na.value="grey",
                           limits = plot_lims, breaks = c(plot_lims[1],0,plot_lims[2]))
    
  }
  
  
  # Save legend
  plot_legend <- get_legend(plot_out)
  plot_out <- plot_out  + theme(legend.position="none")
  
  return(list(plot = plot_out, legend = plot_legend))
}


plot_global_measures <- function(data, raw_data, plot_type, comp_col, comp_col_1, comp_col_2) {
  
  # Plot Global
  if (plot_type == "zscore") {
    raw_data_suf <- "Transformed.q.zscore"
    sig_pos <- 4.7
    plot_ylab <- "Deviation Score"
  } else if (plot_type == "centile") {
    raw_data_suf <- "Transformed.q.wre"
    sig_pos <- 1.0
    plot_ylab <- "Centile"
  } else if (plot_type == "volume") {
    raw_data_suf <- "Transformed.normalised"
    sig_pos <- 75
    plot_ylab <- "Normalized Value"
  }
  
  tmp <- raw_data %>% 
    select(participant, any_of(comp_col), sex, age_days, (ends_with(raw_data_suf) &
                                                            starts_with(global_measures))) %>%
    pivot_longer(cols = contains(raw_data_suf), names_to = "Region", values_to = "Value") %>%
    mutate(Region = str_remove_all(Region,raw_data_suf),
           Region = factor(Region, 
                           levels = global_measures,
                           labels = global_measure_names)) %>%
    na.omit()
  
  sig_tab <- data %>% 
    select(label, Cohens_d, p.value, p.fdr) %>% 
    filter(label %in% global_measures) %>% 
    mutate(group1 = comp_col_1, group2 = comp_col_2,
           .y. = "Value") %>% 
    select(Region = label,.y., group1, group2, p.value, p.fdr) %>%
    mutate(p.signif = case_when(p.fdr <= 0.0001 ~ "****",
                                p.fdr <= 0.001 ~ "***",
                                p.fdr <= 0.01 ~ "**",
                                p.fdr <= 0.05 ~ "*"),
           Region = reduce2(global_measures, global_measure_names,  .init = Region,
                            ~ str_replace(..1, ..2, ..3)),
           # Region = str_replace_all(Region,"eTIV","ICV"),
           # Region = str_replace_all(Region,"meanCT","Mean CT"),
           # Region = str_replace_all(Region,"totalSA","Total SA"),
           y.position = sig_pos)
  
  # for(i in seq_along(words_to_replace)){
  #   strings <- str_replace_all(strings, words_to_replace[i], replacements[i])
  # }
  
  plot_global <- ggplot(tmp, aes(x = Region, color = .data[[comp_col]], y = Value)) +
    geom_point(position = position_jitterdodge(), size = 0.5) +
    geom_boxplot(outliers = F, alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    ylab(plot_ylab) +
    theme_pubclean() +
    scale_color_manual(values = c("darkgrey","firebrick"), drop=FALSE) +
    xlab("") + 
    stat_pvalue_manual(sig_tab, x = "Region", label = "p.signif") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.position="bottom")  + 
    guides(fill=guide_legend(title=""))
  
  return(plot_global)
}


plot_effect_map <- function(data, raw_data, fill_col, color_scheme, plot_lims,
                            plot_type = "zscore", comp_col = "dx", comp_col_1 = "CN", comp_col_2 = "22q11DS",
                            colorbar_title = "",sig_col = NA, cats = c("GM","SA","CT"),
                            cort_only = FALSE, formatted = TRUE, bar_global = F) {
  
  # Make plots for GM, CT, and SA
  plots_cortical <- list("Blank" = NULL)
  
  for (cat in cats) {
    cat_name <- case_match(cat,
                           "GM" ~ "Gray Matter Volume",
                           "SA" ~ "Surface Area",
                           "CT" ~ "Cortical Thickness")
    
    data_plot <- data %>% 
      filter(grepl(glue("{cat}\\."),label)) %>%
      mutate(label = str_remove_all(label, glue("{cat}\\.")))
    
    # Cortical Plot
    plot_cortical   <- plot_atlas(data_plot, dk_plot, fill_col, color_scheme, plot_lims,
                                  colorbar_title = colorbar_title,sig_col = sig_col,
                                  position = "stacked")
    
    plots_cortical[[cat]] <- plot_cortical$plot + 
      ggtitle(cat_name) + 
      theme(legend.position = "bottom",
            plot.title.position = "plot",
            legend.text=element_text(size=7))
    
    leg <- as_ggplot(get_legend(plots_cortical[[cat]]))
    
  }
  
  p_cort <- ggarrange(plotlist = plots_cortical, ncol = 4,
                      common.legend = T, widths = c(0.15,1,1,1),
                      legend = "none", labels = c("C)","","",""))
  
  # Add subcortical plot
  data_plot <- data %>% 
    filter(grepl(glue("SUBC\\."),label)) %>%
    mutate(label = str_remove_all(label, glue("SUBC\\.")))
  plot_subcortical   <- plot_atlas(data_plot, aseg_plot, fill_col, color_scheme, plot_lims,
                                   colorbar_title = colorbar_title,sig_col = sig_col)
  
  plot_subcortical <- plot_subcortical$plot + 
    ggtitle("Subcortical GMV") + 
    theme(legend.position = "none",
          plot.title.position = "plot",
          legend.text=element_text(size=7),
          plot.margin = margin(l = 25, r = 5.5, t = 15, b = 5.5, unit = "pt"))
  
  # Plot Global
  if (plot_type == "zscore") {
    raw_data_suf <- "Transformed.q.zscore"
    sig_pos <- 4.7
    plot_ylab <- "Deviation Score"
  } else if (plot_type == "centile") {
    raw_data_suf <- "Transformed.q.wre"
    sig_pos <- 1.0
    plot_ylab <- "Centile"
  } else if (plot_type == "volume") {
    raw_data_suf <- "Transformed.normalised"
    sig_pos <- 75
    plot_ylab <- "Normalized Value"
  }
  # print("GLOB")
  if (bar_global == F) {
    plot_global <- plot_global_measures(data, raw_data, plot_type, comp_col, comp_col_1, comp_col_2) +
      theme(plot.margin = margin(l = 30, r = 5.5, t = 10, b = 5.5, unit = "pt"),
            legend.margin = margin(0,0,0,0, unit = "pt"),
            legend.title=element_blank())
  } else {    
    tmp <- data %>% 
      mutate(label = str_replace_all(label,"\\.all",""),
             label = str_replace_all(label,"_Cortex",""),
             label = str_replace_all(label,"_"," "),
             label = str_replace_all(label,"eTIV","ICV"),
             label = str_replace_all(label,"totalSA","SA"),
             label = str_replace_all(label,"meanCT","CT"),
             label = factor(label, levels = global_measure_names)) %>%
      filter(!is.na(label))
    
    plot_global <- ggplot(tmp, aes(x = label, y = .data[[fill_col]], 
                                   ymin = .data[[ci_lo_col]], 
                                   ymax = .data[[ci_hi_col]])) +
      geom_col(aes(fill = .data[[sig_col]]), color = "black") +
      geom_errorbar(width = 0.5) + 
      geom_hline(yintercept = 0, linetype = "dashed") +
      ylab(colorbar_title) +
      theme_pubclean() +
      scale_fill_manual(values = c("white","grey"), drop=FALSE) +
      xlab("") + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.position="bottom")  + 
      guides(fill=guide_legend(title=""))
  }
  
  # print("SUBC")
  
  if (cort_only) {
    return(p_cort)
  } else if (!formatted) {
    return(list("Global" = plot_global,
                "Subcortical" = p_cort, 
                "Cortical" = plots_cortical,
                "Legend" = leg))
  } else {
    p_subcort <- ggarrange(plot_global,plot_subcortical,  ncol = 2, labels = c("A)","B)"),
                           widths = c(0.65, 0.35))
    
    # p <- ggarrange(p_subcort, p_cort, nrow = 2)
    p <- ggarrange(p_subcort, p_cort, leg, nrow = 3, heights = c(1.2, 1, 0.2))
    return(p)
  }
}

plot_outliers_single <- function(data, fill_col, color_scheme, plot_lims, fill_name,
                                 colorbar_title = "",sig_col = NA, cats = c("GM","SUBC","SA","CT"),
                                 cort_only = F, add_colorbar = T, add_fillbar = T, exclude_y_label = F,
                                 formatted = T, lh_only = F, p_col = "p.value",
                                 percent = T) {
  # Make plots for GM, CT, and SA
  if (add_colorbar) {
    leg <- "bottom"
  } else {
    leg <- "none"
  }
  plots <- list()
  for (m_cat in cats) {
    cat_name <- case_match(m_cat,
                           "GM" ~ "Gray Matter Volume",
                           "SUBC" ~ "Subcortical Volume",
                           "SA" ~ "Surface Area",
                           "CT" ~ "Cortical Thickness")
    
    data_plot <- data %>% 
      filter(grepl(glue("{m_cat}\\."),label)) %>%
      mutate(label = str_remove_all(label, glue("{m_cat}\\."))) 
    
    if (!is.na(sig_col)) {
      data_plot <- data_plot %>%
        arrange(get(sig_col))
    }
    if (m_cat == "SUBC") {
      plot   <- plot_atlas(data_plot, aseg_plot, fill_col, color_scheme, plot_lims,
                           colorbar_title = colorbar_title,sig_col = sig_col,
                           p_col = p_col, percent = percent)
    } else {
      plot   <- plot_atlas(data_plot, dk_plot, fill_col, color_scheme, plot_lims,
                           colorbar_title = colorbar_title,sig_col = sig_col,
                           position = "stacked", lh_only = lh_only,
                           p_col = p_col, percent = percent)
    }
    
    # Cortical Plot
    
    
    plots[[glue("{m_cat} {fill_col}")]] <- plot$plot + 
      ggtitle(" ") + 
      theme(legend.position = "none",
            plot.title.position = "plot",
            legend.text=element_text(size=7),
            legend.title=element_text(size=10),
            axis.title.y = element_text(size=10),
            plot.title = element_text(size=10, hjust = 0.15)) + ylab(" ") 
    if (!add_colorbar) {
      plots[[glue("{m_cat} {fill_col}")]] <- plots[[glue("{m_cat} {fill_col}")]] + guides(color = "none")
    }
    if (!add_fillbar) {
      plots[[glue("{m_cat} {fill_col}")]] <- plots[[glue("{m_cat} {fill_col}")]] + guides(fill = "none")
    }
    # scale_fill_stepsn(colors = c("white","firebrick"), 
    #                   breaks = seq(plot_lims[1],plot_lims[2],
    #                                (plot_lims[2]-plot_lims[1])/5),
    #                   limits = plot_lims)
    
    if (m_cat == cats[1] & exclude_y_label == FALSE) {
      plots[[glue("{m_cat} {fill_col}")]] <- plots[[glue("{m_cat} {fill_col}")]] +
        ylab(fill_name)
    }
    plots[[glue("{m_cat} {fill_col}")]] <- plots[[glue("{m_cat} {fill_col}")]] +
      ggtitle(cat_name)
  }
  
  if (formatted) {
    if ("SUBC" %in% cats) {
      if (length(cats) == 4) {
        plot_n_row <- 1
        plot_n_col <- 4
        w = c(1,1,1,0.9)
      } else if (length(cats) == 2) {
        plot_n_row <- 1
        plot_n_col <- 2
        w = c(1,0.75)
      } else if (length(cats) == 3) {
        plot_n_row <- 1
        plot_n_col <- 3
        w = c(1,1,0.75)
      }
    } else {
      plot_n_row <- 1
      plot_n_col <- 3
      w = 1
    }
    p <- ggarrange(plotlist =  plots, nrow = plot_n_row, ncol = plot_n_col, 
                   common.legend = TRUE, legend = "bottom", widths = w,
                   align = "v")
  } else {
    return(plots)
  }
}
# LMSz Functions ---------
backtransform_centiles <- function(NEWData, FITParam,
                                   Reference.Holder=NULL,
                                   MissingToZero=TRUE, NAToZero=TRUE, Prefix="",
                                   Pred.Name="PRED", Only.New=FALSE,
                                   verbose=FALSE ) {
  
  ##
  ## check arguments
  if( is.null(NEWData) || is.null(FITParam) ) { stop("Must supply valid NEWData and FITParam objects") }
  
  if( !is.null(Reference.Holder) ) {
    NEWData <- ValidateCleanInput(IN=NEWData,
                                  Reference.Subset=Reference.Holder$SUBSET,
                                  Reference.Model=Reference.Holder$MODEL,
                                  Reference.Param=Reference.Holder$FIT.EXTRACT$param )
    ## above line checks that input conforms to reference
    ##
  }
  
  ##
  ## vvv not sure this is necessary, but a precaution against mixing with model.frame(), model.matrix(), etc.
  Saved.Attributes <- attributes(NEWData)
  attributes(NEWData)[ !names(attributes(NEWData)) %in% c("names","row.names","class") ]  <- NULL
  
  if( Only.New ) {
    stop("This isn't quite right at the moment, need to fix this to perhaps Only.Modified?")
    IN.NAMES <- names(NEWData)
  }
  
  In.Model <- attr(FITParam,"model")
  
  Model.Columns <- unlist(In.Model$covariates[c("X","BY","OTHER","RANEF")])
  
  FAMILY <- get(FITParam$family)()
  
  Model.Parameters <- names(FAMILY$parameters)
  
  Any.Missing.Columns <- Model.Columns[!c(Model.Columns %in% names(NEWData))]
  if( length(Any.Missing.Columns) > 0 ) {
    NEWData[ , Any.Missing.Columns ] <- NA
    if( verbose ) { cat("Missing columns in NEWData from within Model. Set to NA. [", Any.Missing.Columns,"]\n") }
  }
  
  for( wPARAM in Model.Parameters ) {
    Model.Formula <- as.formula(sprintf(" ~ %s",FITParam[[ wPARAM ]]$equ$fixef))
    Model.Frame <- model.frame(formula=Model.Formula, data=NEWData, na.action=na.pass )
    Model.Matrix <- model.matrix( Model.Formula, Model.Frame, contrasts.arg=FITParam$contrasts, na.action=na.pass )
    
    if( MissingToZero==TRUE ) {
      lASSIGN <- attr(Model.Matrix,"assign")
      lAllNA <- apply(Model.Matrix,2,function(X){all(is.na(X))})
      if( length(lASSIGN)!=length(lAllNA) ) {stop("FATAL ERROR: lASSIGN different length from lAllNA")}
      if( sum(lAllNA) > 0 ) {
        lLABELS <- paste( attr(terms(Model.Formula),"term.labels")[unique(lASSIGN[lAllNA])], collapse=", " )
        if( verbose ) { 
          cat(wPARAM,"Some columns of the MODEL MATRIX are all NAs (",lLABELS,"), setting to zero in MODEL MATRIX.\n")
          cat(wPARAM,"NOTE: This is mainly acceptable for factors using the contr.sum() encoding, otherwise it may give unexpected results.\n")
        }
        Model.Matrix[ , which(lAllNA) ]  <- 0
      }
    }
    if( NAToZero==TRUE ) {
      local.sum <- sum(apply(Model.Matrix,1,function(X){any(is.na(X))}))
      if( local.sum > 0 ) {
        ARR.IND <- which( is.na(Model.Matrix), arr.ind=TRUE )
        
        lCOLUMNS <- paste( attr(terms(Model.Formula),"term.labels")[unique(attr(Model.Matrix,"assign")[unique(ARR.IND[,"col"])])], collapse=", " )
        STR <- sprintf("MODEL MATRIX has %i rows (across columns: %s ) with NAs, all set to zero",local.sum,lCOLUMNS)
        if( verbose ) { 
          cat(wPARAM,STR,"\n")
          cat(wPARAM,"NOTE: This is mainly acceptable for factors using the contr.sum() encoding, otherwise it may give unexpected results.\n")
        }
        Model.Matrix[ ARR.IND ]  <- 0
      }
    }
    if( !all( colnames(Model.Matrix) %in% names(FITParam[[ wPARAM ]]$fixef) ) ) {stop("FITParam inconsistent with INPUT")}
    Fit.Fixef <- matrix( FITParam[[ wPARAM ]]$fixef[colnames(Model.Matrix)], ncol=1, dimnames=list(colnames(Model.Matrix),"Beta") )
    
    NEWData[,sprintf("%s%s.pop",Prefix,wPARAM)] <- as.vector(Model.Matrix %*% Fit.Fixef)
    if( !is.na( FITParam[[ wPARAM ]]$equ$ranef ) ) {
      ##
      ## In this block we find the matching random-effect (from FITParam) and insert it into NEWData
      ## NOTE: This code block assumes only one random-effect (per gamlss-component) [ENHANCEMENT PLAN: expand to multiple random-effects]
      ## 
      POSITION.IN.MAP <- as.numeric(NEWData[, FITParam[[ wPARAM ]]$equ$ranef ])
      MAPPING <- match( levels(NEWData[, FITParam[[ wPARAM ]]$equ$ranef ]), names(FITParam[[ wPARAM ]]$ranef) )
      NEWData[,sprintf("%s%s.ranef",Prefix,wPARAM)] <- FITParam[[ wPARAM ]]$ranef[ MAPPING[POSITION.IN.MAP] ]
      
      NEWData[,sprintf("%s%s.wre",Prefix,wPARAM)] <- NEWData[,sprintf("%s%s.pop",Prefix,wPARAM)] + NEWData[,sprintf("%s%s.ranef",Prefix,wPARAM)]
    } else {
      NEWData[,sprintf("%s%s.ranef",Prefix,wPARAM)] <- Inf ## to differentiate from missing random-effect levels (NAs), use Inf if no ranef at all
      NEWData[,sprintf("%s%s.wre",Prefix,wPARAM)] <- NEWData[,sprintf("%s%s.pop",Prefix,wPARAM)]
    }
  }
  
  ARGUMENTS <- list()
  for( lTYPE in c("pop","wre") ) {
    ARGUMENTS[[ lTYPE ]] <- list()
    for( wPARAM in Model.Parameters ) {
      LinkInvFun <- FAMILY[[sprintf(sprintf("%s.linkinv",wPARAM))]]
      ComponentValues <- NEWData[,sprintf("%s%s.%s",Prefix,wPARAM,lTYPE)]
      
      ARGUMENTS[[ lTYPE ]][[ wPARAM ]] <- LinkInvFun( ComponentValues )
    }
  }
  
  ##
  ## Select rows (within arguments) that are valid, ie. not NAs
  KEEP <- lapply( ARGUMENTS, function(Y){ as.vector( Reduce(f=`&`,x=lapply(Y,function(X){!is.na(X)})), mode="logical") } )
  SHORT <- setNames(lapply(1:length(ARGUMENTS),function(IDX){ lapply(ARGUMENTS[[IDX]],function(X){ as.vector(X[ KEEP[[IDX]] ]) } ) } ), names(KEEP) )
  WHICH <- lapply( KEEP, which )
  
  ##
  ## in below we make versions that account for NAs in outcome (OUT.NAME) column
  OUT.NAME <- In.Model$covariates$Y
  if( (OUT.NAME %in% names(NEWData)) && any(!is.na(NEWData[,OUT.NAME])) ) {
    out.KEEP.ARGS <- lapply( ARGUMENTS, function(Y){ as.vector( Reduce(f=`&`,x=lapply(Y,function(X){!is.na(X)})), mode="logical") } )
    out.KEEP <- lapply( out.KEEP.ARGS, function(X) { (X) & (!is.na(NEWData[,OUT.NAME])) } )
    out.SHORT <- setNames(lapply(1:length(ARGUMENTS),function(IDX){ lapply(ARGUMENTS[[IDX]],function(X){ as.vector(X[ out.KEEP[[IDX]] ]) } ) } ), names(out.KEEP) )
    out.WHICH <- lapply( out.KEEP, which )
    
    KEEP.C <- Reduce(`&`,out.KEEP)
    SHORT.C <- lapply( ARGUMENTS, function(Y){ lapply(Y,function(X){ X[KEEP.C] } ) } )
    WHICH.C <- which( KEEP.C )
  }
  
  for( lTYPE in c("pop","wre") ) {
    if( length(WHICH[[lTYPE]])>0 ) {
      
      NEWData[WHICH[[lTYPE]],sprintf("%s%s.%s",Prefix,Pred.Name,lTYPE)] <- do.call(what=get(paste0("q",FAMILY$family[1])),
                                                                                   args=c(SHORT[[lTYPE]],list(p=NEWData$value)))
      if( (OUT.NAME %in% names(NEWData)) && any(!is.na(NEWData[,OUT.NAME])) ) {
        
        NEWData[out.WHICH[[lTYPE]],sprintf("%s%s.q.%s",Prefix,OUT.NAME,lTYPE)] <- do.call(what=get(paste0("p",FAMILY$family[1])),
                                                                                          args=c(out.SHORT[[lTYPE]],
                                                                                                 list(q=NEWData[out.WHICH[[lTYPE]],OUT.NAME,drop=TRUE])))
      }
    }
  }
  
  if( Only.New ) {
    NEWData <- NEWData[,names(NEWData)[!names(NEWData)%in%IN.NAMES]]
  }
  attributes(NEWData) <- c( attributes(NEWData), Saved.Attributes[!names(Saved.Attributes) %in% c("names","row.names","class")] )
  
  return( NEWData )
}

sim_list_slip <- list()
sim_list_slip$Male <- expand.grid(AgeTransformed = seq(log(280),log(365.25*21+280),0.01),
                                  sex = c("Male"))
sim_list_slip$Female <- expand.grid(AgeTransformed = seq(log(280),log(365.25*21+280),0.01),
                                    sex = c("Female"))


backtransform_centile_fans <- function(pheno) {
  # Load data
  pheno_name <- glue("{pheno}Transformed.q.wre")
  norm_name <- glue("{pheno}Transformed.normalised")
  zscore_name <- glue("{pheno}Transformed.q.zscore")
  
  old_pheno <- slip_features$old_pheno[slip_features$new_pheno == pheno]
  if (length(old_pheno) == 0) {
    return()
  }
  growthChartModel <- readRDS(glue("lmsz/lmsz_models_slip/SLIP_lmsz_{old_pheno}Transformed.rds"))
  sim_list <- sim_list_slip
  # Load SLIP fit model
  if (grepl("CT",pheno)) {
    fit_file <- glue("/mnt/isilon/bgdlab_processing/braincharts/SLIP/recon-all-clinical/gamlss/RDS/mpr-{old_pheno}Transformed/FIT.EXTRACT.rds")
  } else {
    fit_file <- glue("/mnt/isilon/bgdlab_processing/braincharts/SLIP/recon-all-clinical/gamlss/RDS/rac-{old_pheno}Transformed/FIT.EXTRACT.rds")
  }
  
  orig_fit <- readRDS(fit_file)
  
  df_tmp <- df %>%
    select(participant, AgeTransformed, sex, !!pheno_name, !!norm_name) %>%
    rename(pheno = !!pheno_name,
           normalised_value = !!norm_name) %>%
    mutate(pheno = !!pheno_name,
           zscore = !!zscore_name) %>%
    na.omit()
  
  # Simulate centiles
  lmsz_sim <- centile_predict(gamlssModel = growthChartModel, 
                              sim_df_list = sim_list, 
                              x_var = "AgeTransformed", 
                              desiredCentiles = c(0.025,0.250,0.50,0.750,0.975),
                              df = df_tmp,
                              average_over = FALSE,
                              resid_terms = FALSE)
  
  # Reformat simulated values
  lmsz_sim$fanCentiles_Female <- lmsz_sim$fanCentiles_Female %>%
    pivot_longer(-c(AgeTransformed),
                 names_to = "centile_cat",
                 values_to = "value") %>%
    mutate(sex = "Female")
  
  lmsz_sim$fanCentiles_Male <- lmsz_sim$fanCentiles_Male %>%
    pivot_longer(-c(AgeTransformed),
                 names_to = "centile_cat",
                 values_to = "value") %>%
    mutate(sex = "Male")
  
  lmsz_df <- bind_rows(lmsz_sim$fanCentiles_Female,
                       lmsz_sim$fanCentiles_Male) %>%
    mutate(site = NA,
           sex = factor(sex, levels = c("Female","Male")),
           LMSz_input_fans = value,
           LMSz_input_centile = pnorm(value)) %>%
    mutate(line_size = case_match(centile_cat,
                                  c("cent_0.025","cent_0.975") ~ 0.3,
                                  c("cent_0.25","cent_0.75") ~ 0.3,
                                  c("cent_0.5") ~ 0.4),
           line_type = case_match(centile_cat,
                                  c("cent_0.025","cent_0.975") ~ "dotted",
                                  c("cent_0.25","cent_0.75") ~ "longdash",
                                  c("cent_0.5") ~ "solid"),
           .after = sex) %>%
    select(-c(value))
  
  # Plot original SLIP Model
  NEWData <- lmsz_df %>%
    select(AgeTransformed, sex) %>% 
    distinct() %>%
    as.data.frame()
  
  original_fans <- Apply.Param(NEWData, orig_fit$param,
                               Pred.Set=c("l025"=0.025,"l250"=0.250,
                                          "m500"=0.5,"u750"=0.750,
                                          "u975"=0.975)) %>%
    pivot_longer(c(PRED.l025.pop,PRED.l250.pop,
                   PRED.m500.pop,PRED.u750.pop,
                   PRED.u975.pop),
                 names_to = "centile_cat", values_to = "SLIP_input_fans") %>%
    mutate(centile_cat = case_match(centile_cat,
                                    "PRED.l025.pop" ~ "cent_0.025",
                                    "PRED.l250.pop" ~ "cent_0.25",
                                    "PRED.m500.pop" ~ "cent_0.5",
                                    "PRED.u750.pop" ~ "cent_0.75",
                                    "PRED.u975.pop" ~ "cent_0.975")) %>%
    select(AgeTransformed, sex, centile_cat, SLIP_input_fans)
  
  # Plot backtransformed centiles
  NEWData <- lmsz_df %>%
    select(AgeTransformed, sex,
           centile_cat,
           value = LMSz_input_centile) %>% 
    distinct() %>%
    as.data.frame()
  
  backtransformed_fans <- backtransform_centiles(NEWData, orig_fit$param)
  
  backtransformed_fans <- backtransformed_fans %>%
    select(AgeTransformed, sex, centile_cat, 
           LMSz_fans = PRED.pop)
  
  # Join all results
  out_df <- lmsz_df %>%
    inner_join(original_fans, by= join_by(AgeTransformed, centile_cat, sex)) %>%
    inner_join(backtransformed_fans, by= join_by(AgeTransformed, centile_cat, sex)) %>%
    mutate(pheno = pheno, model = "SLIP", .before = AgeTransformed) %>%
    mutate(fit_file = fit_file)
  
  return(out_df)
}


pred_og_centiles <- function(gamlssModel, og.data, get.std.scores = FALSE, new.data=NULL){
  pheno <- gamlssModel$mu.terms[[2]]
  
  #subset df cols just to predictors from model
  predictor_list <- list_predictors(gamlssModel)
  stopifnot("Dataframe columns and model covariates don't match" = 
              predictor_list %in% names(og.data))
  if (is.null(new.data)) {
    newData <- og.data
    predict_me <- og.data
  } else {
    stopifnot("Dataframe columns and model covariates don't match" = 
                predictor_list %in% names(new.data))
    newData <- new.data
    predict_me <- new.data
    #make sure all vals are within range of those originally modeled
    check_range(subset(og.data, select = predictor_list), newData)
  }
  
  #predict
  predModel <- predictAll(gamlssModel, newdata=newData, data=og.data, type= "response")
  
  #get dist type (e.g. GG, BCCG) and write out function
  fname <- gamlssModel$family[1]
  pfun <- paste0("p", fname)
  
  #look for moments
  has_sigma <- "sigma" %in% gamlssModel[[2]]
  has_nu <- "nu" %in% gamlssModel[[2]]
  has_tau <- "tau" %in% gamlssModel[[2]]
  
  centiles <- c()
  #iterate through participants
  for (i in 1:nrow(predict_me)){
    cent_args <- list(predict_me[[pheno]][[i]], predModel$mu[[i]])
    
    if (has_sigma){
      cent_args$sigma <- predModel$sigma[[i]]
    }
    if (has_nu){
      cent_args$nu <- predModel$nu[[i]]
    } 
    if (has_tau){
      cent_args$tau <- predModel$tau[[i]]
    } 
    
    centiles[i] <- do.call(pfun, cent_args)
    
    #don't let centile = 1 (for z-scores)!
    if (centiles[i] == 1) {
      centiles[i] <- 0.9999999999999999
    }
    #don't let centile = 0 (for z-scores)!
    if (centiles[i] == 0) {
      centiles[i] <- 0.0000000000000001 #25 dec places
    }
    
  }
  if (get.std.scores == FALSE){
    return(centiles)
  } else {
    #get 'z scores' from normed centiles - how z.score() does it
    rqres <- qnorm(centiles)
    
    #return dataframe
    df <- data.frame("centile" = centiles,
                     "std_score" = rqres)
    return(df)
  } 
  
}
