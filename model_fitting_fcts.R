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

fit_bam_interactions <- function(data, main_models, phenotypes, predictor, moderator, 
                                 covariates = NULL, random_effect = NULL, save = FALSE, save_dir = NULL, 
                                 save_name = NULL){
  if(save && (is.null(save_dir) || is.null(save_name))) stop("save_dir and save_name required when save = TRUE")
  
  tab <- extract_bam_table(main_models)
  sig <- rownames(tab)[tab$full_p_adj < 0.05 | tab$nolin_p_adj < 0.05]
  
  models <- lapply(sig, function(ph){
    covs <- if(!is.null(covariates)) paste(covariates, collapse = " + ") else ""
    re <- if(!is.null(random_effect)) paste0("s(", random_effect, ", bs='re')") else ""
    rhs_base <- paste(c(covs, re)[c(covs, re) != ""], collapse = " + ")
    
    nonlinear <- tab[ph, "full_p_adj"] < 0.05 && tab[ph, "nolin_p_adj"] < 0.05
    
    if(nonlinear){
      f <- as.formula(paste(ph, "~ s(", predictor, ") + s(", moderator, ") + ti(", predictor, ",", moderator, ")",
                            if(rhs_base != "") paste("+", rhs_base)))
      mod <- mgcv::bam(f, data = data, method = "fREML", discrete = TRUE)
      attr(mod, "interaction") <- "nonlinear"
    } else {
      f <- as.formula(paste(ph, "~", predictor, "*", moderator,
                            if(rhs_base != "") paste("+", rhs_base)))
      mod <- mgcv::bam(f, data = data, method = "fREML", discrete = TRUE)
      attr(mod, "interaction") <- "linear"
    }
    
    mod
  })
  
  names(models) <- sig
  attr(models, "predictor") <- predictor
  attr(models, "moderator") <- moderator
  attr(models, "type") <- "interaction"
  
  if(save) saveRDS(models, file.path(save_dir, save_name))
  invisible(models)
}
