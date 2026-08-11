library(gratia)

# Functions to extract info from models 
extract_bam_table <- function(models, predictor){
  tab <- as.data.frame(t(sapply(models, function(m){
    s_full <- summary(m$full)
    s_lin  <- summary(m$nolin)
    
    full_s <- s_full$s.table[1,]
    lin_s  <- s_lin$s.table[1,]
    
    # linear effect beta from parametric table
    lin_p <- s_lin$p.table[grep(predictor, rownames(s_lin$p.table)), ]
    
    c(
      full_edf = unname(full_s[1]),
      full_p = unname(full_s[4]),
      
      nolin_edf = unname(lin_s[1]),
      nolin_p = unname(lin_s[4]),
      
      lin_beta = unname(lin_p[1]),
      lin_beta_p = unname(lin_p[4])
    )
  })))
  
  tab$full_p_adj  <- p.adjust(tab$full_p, "BH")
  tab$nolin_p_adj <- p.adjust(tab$nolin_p, "BH")
  tab$lin_beta_p_adj <- p.adjust(tab$lin_beta_p, "BH")
  
  invisible(tab)
}

extract_bam_features_main <- function(models, data, predictor, covariates = NULL, random_effect = NULL){
  
  tab <- extract_bam_table(models, predictor)
  
  feats <- lapply(names(models), function(ph){
    
    m_full  <- models[[ph]]$full
    m_nolin <- models[[ph]]$nolin
    
    # rebuild reduced model WITHOUT predictor entirely
    covs <- if(!is.null(covariates)) paste(covariates, collapse = " + ") else ""
    re <- if(!is.null(random_effect)) paste0("s(", random_effect, ", bs='re')") else ""
    
    rhs <- paste(c(covs, re)[c(covs, re) != ""], collapse = " + ")

    f_red <- as.formula(paste(ph, "~ 1", 
                              if(rhs != "") paste("+", rhs)))
    
    f_onlylin <- as.formula(paste(ph, "~ 1 + ", predictor,
                                  if(rhs != "") paste("+", rhs)))
    
    # subset data to only include predictor is not NA to fit reduced model
    data_fit <- data %>% filter(!is.na(get(predictor)) & !is.na(get(ph)))
    m_red <- mgcv::bam(f_red, data = data_fit, method = "fREML", discrete = TRUE)
    m_onlylin <- mgcv::bam(f_onlylin, data = data_fit, method = "fREML", discrete = TRUE)
    
    # partial R2: full predictor contribution vs no predictor model
    r2_full_vs_none <- 1 - deviance(m_full)/deviance(m_red)
    #r2_full_vs_none <- (summary(m_full)$r.sq - summary(m_red)$r.sq) /  (1 - summary(m_red)$r.sq)
    
    # partial R2: spline contribution over linear model
    r2_spline_vs_linear <- 1 - deviance(m_full)/deviance(m_onlylin)
    #r2_spline_vs_linear <- (summary(m_full)$r.sq - summary(m_onlylin)$r.sq) /  (1 - summary(m_onlylin)$r.sq)
    
    # average first derivative
    d <- derivatives(m_full)
    avg_deriv <- mean(d$.derivative)
    
    # region of effect (CI crossing rule)
    d <- d %>%
      mutate (sig = !(0 >.lower_ci & 0 < .upper_ci))
    region <- if(sum(d$sig) > 0) c(min(d[d$sig == T, predictor]), max(d[d$sig == T, predictor])) else c(NA, NA)
    
    c(
      r2_full_vs_none = r2_full_vs_none,
      r2_spline_vs_linear = r2_spline_vs_linear,
      avg_first_deriv = avg_deriv,
      region_start = region[1],
      region_end = region[2]
    )
  })
  
  cbind(tab, do.call(rbind, feats))
}

extract_bam_interaction_table <- function(models,
                                          interaction_type = c("continuous", "ordered")){
  
  interaction_type <- match.arg(interaction_type)
  
  if(interaction_type == "continuous"){
    
    tab <- as.data.frame(t(sapply(models, function(m){
      
      s <- summary(m)$s.table
      
      c(
        pred_edf = s[1, 1],
        pred_p   = s[1, 4],
        int_edf  = s[3, 1],
        int_p    = s[3, 4]
      )
    })))
    
  } else {
    
    # determine maximum number of interaction smooths
    max_int <- max(sapply(models, function(m) nrow(summary(m)$s.table) - 1))
    
    tab <- t(sapply(models, function(m){
      
      s <- summary(m)$s.table
      
      out <- c(
        pred_edf = s[1, 1],
        pred_p   = s[1, 4]
      )
      
      for(i in seq_len(max_int)){
        
        if(i + 1 <= nrow(s)){
          out[paste0("int", i, "_edf")] <- s[i + 1, 1]
          out[paste0("int", i, "_p")]   <- s[i + 1, 4]
        } else {
          out[paste0("int", i, "_edf")] <- NA
          out[paste0("int", i, "_p")]   <- NA
        }
        
      }
      
      out
      
    }))
    
    tab <- as.data.frame(tab, check.names = FALSE)
    
  }
  
  p_cols <- grep("_p$", names(tab), value = TRUE)
  
  for(col in p_cols){
    tab[[paste0(col, "_adj")]] <- p.adjust(tab[[col]], method = "BH")
  }
  
  invisible(tab)
}

extract_bam_features_interaction <- function(models,
                                             data,
                                             predictor,
                                             moderator,
                                             covariates = NULL,
                                             random_effect = NULL,
                                             interaction_type = c("continuous",
                                                                  "ordered")){
  
  interaction_type <- match.arg(interaction_type)
  
  tab <- extract_bam_interaction_table(
    models,
    interaction_type = interaction_type
  )
  
  feats <- lapply(names(models), function(ph){
    
    m_full <- models[[ph]]
    
    covs <- if(!is.null(covariates))
      paste(covariates, collapse = " + ")
    else
      ""
    
    re <- if(!is.null(random_effect))
      paste0("s(", random_effect, ", bs='re')")
    else
      ""
    
    rhs_base <- paste(c(covs, re)[c(covs, re) != ""], collapse = " + ")
    
    if(interaction_type == "continuous"){
      
      f_red <- as.formula(
        paste(
          ph,
          "~ s(", predictor, ") + s(", moderator, ")",
          if(rhs_base != "") paste("+", rhs_base)
        )
      )
      
      data_fit <- data |>
        dplyr::filter(!is.na(.data[[predictor]]) &
                        !is.na(.data[[moderator]]))
      
    } else {
      
      f_red <- as.formula(
        paste(
          ph,
          "~ s(", predictor, ") +", moderator,
          if(rhs_base != "") paste("+", rhs_base)
        )
      )
      
      data_fit <- data |>
        dplyr::filter(!is.na(.data[[predictor]]) &
                        !is.na(.data[[moderator]]))
      
    }
    
    m_red <- mgcv::bam(
      f_red,
      data = data_fit,
      method = "fREML",
      discrete = TRUE
    )
    
    c(
      partial_r2 = 1 - deviance(m_full) / deviance(m_red)
    )
    
  })
  
  cbind(tab, do.call(rbind, feats))
}

extract_bam_resi <- function(models, term) {
  
  tab <- as.data.frame(
    t(
      sapply(models, function(m) {
        
        s <- summary(m$full)
        
        # Require exact match to the requested smooth term
        if (!term %in% rownames(s$s.table)) {
          stop(
            "Could not find exact smooth term: ",
            term
          )
        }
        
        smooth <- s$s.table[term, ]
        
        # Extract values needed for RESI
        F_stat <- unname(smooth["F"])
        edf    <- unname(smooth["edf"])
        ref_df <- unname(smooth["Ref.df"])
        p      <- unname(smooth["p-value"])
        
        rdf <- m$full$df.residual
        n <- nobs(m$full)
        
        # Calculate RESI
        resi <- RESI::f2S(
          f = F_stat,
          df = ref_df,
          rdf = rdf,
          n = n
        )
        
        c(
          term = term,
          F = F_stat,
          edf = edf,
          Ref.df = ref_df,
          rdf = rdf,
          n = n,
          RESI = resi,
          p = p
        )
      })
    ),
    stringsAsFactors = FALSE
  )
  
  # Convert numeric columns back to numeric
  numeric_cols <- c(
    "F",
    "edf",
    "Ref.df",
    "rdf",
    "n",
    "RESI",
    "p"
  )
  
  tab[numeric_cols] <- lapply(
    tab[numeric_cols],
    as.numeric
  )
  
  # BH adjustment across models
  tab$p_adj <- p.adjust(tab$p, method = "BH")
  
  invisible(tab)
}

test_interaction_lrt <- function(data, 
                                 phenotypes,
                                 predictor,
                                 moderator,
                                 interaction = "continuous",
                                 covariates = NULL,
                                 save = FALSE,
                                 method = "ML") {
  
  m1 <- fit_bam_interactions(
    data = data,
    phenotypes = phenotypes,
    predictor = predictor,
    moderator = moderator,
    interaction = interaction,
    covariates = covariates,
    save = save,
    method = method
  )
  
  m0 <- fit_bam_models(
    data = data,
    phenotypes = phenotypes,
    predictor = predictor,
    covariates = covariates,
    save = save,
    method = method
  )
  
  # Compare corresponding models
  results <- Map(
    function(mod1, mod0) {
      
      aov <- anova(
        mod0$full,
        mod1,
        test = "Chisq"
      )
      
      data.frame(
        p = aov$`Pr(>Chi)`[2]
      )
    },
    m1,
    m0
  )
  
  results <- dplyr::bind_rows(results, .id = "phenotype")
  
  results$p_adj <- p.adjust(
    results$p,
    method = "BH"
  )
  
  return(results)
}

# extract_bam_interaction_table <- function(models){
#   tab <- as.data.frame(t(sapply(models, function(m){
#     s <- summary(m)$s.table
#       c(pred_edf = s[1,1], pred_p = s[1,4],
#         mod_edf = s[2,1], mod_p = s[2,4],
#         int_edf = s[3,1], int_p = s[3,4])
#   })))
#   
#   num_cols <- names(tab)[grepl("_p$", names(tab))]
#   for(nm in num_cols) tab[[paste0(nm, "_adj")]] <- p.adjust(tab[[nm]], "BH")
#   invisible(tab)
# }
# 
# extract_bam_features_interaction <- function(models, data, predictor, moderator, 
#                                              covariates = NULL, random_effect = NULL){
#   
#   tab <- extract_bam_interaction_table(models)
#   
#   feats <- lapply(names(models), function(ph){
#       m <- models[[ph]]
#     
#       m_full <- m
#       
#       # reduced = drop interaction term
#       covs <- if(!is.null(covariates)) paste(covariates, collapse = " + ") else ""
#       re <- if(!is.null(random_effect)) paste0("s(", random_effect, ", bs='re')") else ""
#       rhs_base <- paste(c(covs, re)[c(covs, re) != ""], collapse = " + ")
#       
#       f_red <- as.formula(paste(ph, "~ s(", predictor, ") + s(", moderator, ")",
#                             if(rhs_base != "") paste("+", rhs_base)))
#       
#       data_fit <- data %>% filter(!is.na(get(predictor)) & !is.na(get(moderator)))
#       m_red <- mgcv::bam(f_red, data = data_fit, method = "fREML", discrete = TRUE)
#       
#       r2 <- 1 - deviance(m_full)/deviance(m_red)
#     
#       c(partial_r2 = r2)
#   })
#   
#   cbind(tab, do.call(rbind, feats))
# }

append_bam_summary <- function(table, model_name, predictor, covariates, test_term, output_file){
  n_sig <- sum(table[[test_term]] < 0.05, na.rm = TRUE)
  sig <- paste(rownames(table)[table[[test_term]] < 0.05], collapse = "; ")
  
  row <- data.frame(
    date = Sys.Date(),
    model_name = model_name,
    predictor = predictor,
    covariates = if(is.null(covariates)) "" else paste(covariates, collapse = "; "),
    test_term = test_term,
    n_sig = n_sig,
    sig_phenotypes = sig
  )
  
  out <- if(file.exists(output_file)) rbind(read.csv(output_file), row) else row
  write.csv(out, output_file, row.names = FALSE)
  
  print(paste("Number significant:", n_sig))
  print(sig)
  
  invisible(out)
}
