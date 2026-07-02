# Functions to fit models

fit_bam_models <- function(data, phenotypes, predictor, covariates = NULL, 
                           random_effect = NULL, save = FALSE, save_dir = NULL, 
                           save_name = NULL){
  
  if(save && (is.null(save_dir) || is.null(save_name))) stop("save_dir and save_name required when save = TRUE")
  
  models <- lapply(phenotypes, function(ph){
    covs <- if(!is.null(covariates)) paste(covariates, collapse = " + ") else ""
    re <- if(!is.null(random_effect)) paste0("s(", random_effect, ", bs='re')") else ""
    
    rhs <- paste(c(covs, re)[c(covs, re) != ""], collapse = " + ")
    
    f_full <- as.formula(paste(ph, "~ s(", predictor, ")",
                               if(rhs != "") paste("+", rhs)))
    f_lin  <- as.formula(paste(ph, "~ s(", predictor, ", m=c(2,0)) +", predictor,
                               if(rhs != "") paste("+", rhs)))
    
    list(
      full = mgcv::bam(f_full, data = data, method = "fREML", discrete = TRUE),
      nolin = mgcv::bam(f_lin, data = data, method = "fREML", discrete = TRUE)
    )
  })
  
  names(models) <- phenotypes
  attr(models, "predictor") <- predictor
  attr(models, "covariates") <- covariates
  attr(models, "random_effect") <- random_effect
  attr(models, "type") <- "main_effect"
  
  if(save) saveRDS(models, file.path(save_dir, save_name))
  invisible(models)
}

fit_bam_interactions <- function(data, phenotypes, predictor, moderator, 
                                 covariates = NULL, random_effect = NULL, save = FALSE, save_dir = NULL, 
                                 save_name = NULL){
  if(save && (is.null(save_dir) || is.null(save_name))) stop("save_dir and save_name required when save = TRUE")
  
  models <- lapply(phenotypes, function(ph){
    covs <- if(!is.null(covariates)) paste(covariates, collapse = " + ") else ""
    re <- if(!is.null(random_effect)) paste0("s(", random_effect, ", bs='re')") else ""
    rhs_base <- paste(c(covs, re)[c(covs, re) != ""], collapse = " + ")
    
    f <- as.formula(paste(ph, "~ s(", predictor, ") + s(", moderator, ") + ti(", predictor, ",", moderator, ")",
                            if(rhs_base != "") paste("+", rhs_base)))
      mod <- mgcv::bam(f, data = data, method = "fREML", discrete = TRUE)
      mod
  })
  
  names(models) <- phenotypes
  attr(models, "predictor") <- predictor
  attr(models, "moderator") <- moderator
  attr(models, "type") <- "interaction"
  
  if(save) saveRDS(models, file.path(save_dir, save_name))
  invisible(models)
}


bootstrap_partial_r2 <- function(models, data, predictor,
                                 covariates = NULL,
                                 random_effect = NULL,
                                 N = 100,
                                 seed = 1,
                                 cores = 1){
  
  set.seed(seed)
  
  get_r2 <- function(models_b, data_b){
    feats <- lapply(names(models_b), function(ph){
      
      m_full <- models_b[[ph]]$full
      m_nolin <- models_b[[ph]]$nolin
      
      covs <- if(!is.null(covariates)) paste(covariates, collapse = " + ") else ""
      re <- if(!is.null(random_effect)) paste0("s(", random_effect, ", bs='re')") else ""
      rhs <- paste(c(covs, re)[c(covs, re) != ""], collapse = " + ")
      
      f_red <- as.formula(paste(ph, "~", rhs))
      f_onlylin <- as.formula(paste(ph, "~", predictor,
                                    if(rhs != "") paste("+", rhs)))
      
      data_fit <- data_b %>% dplyr::filter(!is.na(.data[[predictor]]))
      
      m_red <- mgcv::bam(f_red, data = data_fit, method = "fREML", discrete = TRUE)
      m_onlylin <- mgcv::bam(f_onlylin, data = data_fit, method = "fREML", discrete = TRUE)
      
      r2_full_vs_none <- 1 - deviance(m_full) / deviance(m_red)
      r2_spline_vs_linear <- 1 - deviance(m_full) / deviance(m_onlylin)
      
      c(r2_full_vs_none = r2_full_vs_none,
        r2_spline_vs_linear = r2_spline_vs_linear)
    })
    
    do.call(rbind, feats)
  }
  
  boot_mat <- lapply(seq_len(N), function(i){
    
    idx <- sample(seq_len(nrow(data)), replace = TRUE)
    data_b <- data[idx, ]
    
    models_b <- lapply(models, function(m){
      list(full = m$full, nolin = m$nolin)
    })
    
    get_r2(models_b, data_b)
  })
  
  boot_arr <- simplify2array(boot_mat)
  
  summarize <- function(x){
    c(
      mean = mean(x, na.rm = TRUE),
      sd = sd(x, na.rm = TRUE),
      p2.5 = quantile(x, 0.025, na.rm = TRUE),
      p97.5 = quantile(x, 0.975, na.rm = TRUE)
    )
  }
  
  out <- apply(boot_arr, c(1,2), summarize)
  
  # reshape to table
  res <- lapply(dimnames(out)[[1]], function(ph){
    data.frame(
      phenotype = ph,
      metric = dimnames(out)[[2]],
      t(out[ph,,])
    )
  })
  
  do.call(rbind, res)
}
