# Functions to fit models

# Remove exact duplicate top-level terms from a model formula built as a
# pasted string (e.g. when a moderator is also passed in `covariates`, so it
# ends up both as its own term and inside the covariate string). Splits the
# RHS on '+' at paren-depth 0, drops repeats (keeping the first occurrence),
# and prints the original and de-duplicated formulas whenever it changes
# anything, so silent formula edits are never silent.
dedupe_formula_terms <- function(formula_str) {
  parts <- strsplit(formula_str, "~", fixed = TRUE)[[1]]
  lhs <- trimws(parts[1])
  rhs <- trimws(paste(parts[-1], collapse = "~"))

  chars <- strsplit(rhs, "")[[1]]
  depth <- 0
  start <- 1
  terms_raw <- character(0)
  for (i in seq_along(chars)) {
    ch <- chars[i]
    if (ch == "(") depth <- depth + 1
    if (ch == ")") depth <- depth - 1
    if (ch == "+" && depth == 0) {
      terms_raw <- c(terms_raw, substr(rhs, start, i - 1))
      start <- i + 1
    }
  }
  terms_raw <- c(terms_raw, substr(rhs, start, nchar(rhs)))
  terms_trimmed <- trimws(terms_raw)

  is_dup <- duplicated(terms_trimmed)
  final_str <- paste0(lhs, " ~ ", paste(terms_trimmed[!is_dup], collapse = " + "))

  if (any(is_dup)) {
    message("Repeated term(s) removed from formula: ", paste(terms_trimmed[is_dup], collapse = ", "))
    message("  built formula:  ", formula_str)
    message("  formula used:   ", final_str)
  }

  final_str
}

fit_bam_models <- function(data, phenotypes, predictor, covariates = NULL,
                           random_effect = NULL, save = FALSE, save_dir = NULL, 
                           save_name = NULL,
                           method = "fREML",
                           fx = FALSE,
                           k = 10) {
  
  if(save && (is.null(save_dir) || is.null(save_name))) {
    stop("save_dir and save_name required when save = TRUE")
  }
  
  models <- lapply(phenotypes, function(ph) {
    
    covs <- if(!is.null(covariates)) {
      paste(covariates, collapse = " + ")
    } else ""
    
    re <- if(!is.null(random_effect)) {
      paste0("s(", random_effect, ", bs='re')")
    } else ""
    
    rhs <- paste(c(covs, re)[c(covs, re) != ""], collapse = " + ")
    
    fx_arg <- if(fx) ", fx = TRUE" else ""
    
    f_full <- as.formula(
      paste0(
        ph, " ~ s(", predictor, ", k = ", k, fx_arg, ")",
        if(rhs != "") paste0(" + ", rhs)
      )
    )
    
    f_lin <- as.formula(
      paste0(
        ph, " ~ s(", predictor, ", m = c(2,0)) + ", predictor,
        if(rhs != "") paste0(" + ", rhs)
      )
    )
    
    list(
      full = mgcv::bam(
        f_full, data = data, method = method, discrete = TRUE
      ),
      nolin = mgcv::bam(
        f_lin, data = data, method = method, discrete = TRUE
      )
    )
  })
  
  names(models) <- phenotypes
  attr(models, "predictor") <- predictor
  attr(models, "covariates") <- covariates
  attr(models, "random_effect") <- random_effect
  attr(models, "fx") <- fx
  attr(models, "k") <- k
  attr(models, "type") <- "main_effect"
  
  if(save) saveRDS(models, file.path(save_dir, save_name))
  invisible(models)
}


fit_bam_interactions <- function(data, 
                                 phenotypes, 
                                 predictor, 
                                 moderator, 
                                 interaction = "continuous",
                                 covariates = NULL, 
                                 random_effect = NULL, 
                                 save = FALSE, 
                                 save_dir = NULL, 
                                 save_name = NULL,
                                 method = "fREML",
                                 fx = FALSE,
                                 k = 10) {
  
  if(save && (is.null(save_dir) || is.null(save_name))) {
    stop("save_dir and save_name required when save = TRUE")
  }
  
  models <- lapply(phenotypes, function(ph) {
    
    covs <- if(!is.null(covariates)) {
      paste(covariates, collapse = " + ")
    } else ""
    
    re <- if(!is.null(random_effect)) {
      paste0("s(", random_effect, ", bs='re')")
    } else ""
    
    rhs_base <- paste(
      c(covs, re)[c(covs, re) != ""],
      collapse = " + "
    )
    
    if (interaction == "continuous") {

      f_str <- paste0(
        ph,
        " ~ s(", predictor, ", k = ", k, ", fx = ", fx, ")",
        " + s(", moderator, ", k = ", k, ", fx = ", fx, ")",
        " + ti(", predictor, ", ", moderator,
        ", k = c(", k, ", ", k, "), fx = ", fx, ")",
        if(rhs_base != "") paste0(" + ", rhs_base)
      )

    } else if(interaction == "categorical") {

      f_str <- paste0(
        ph,
        " ~ ", moderator,
        " + s(", predictor, ", k = ", k,
        ", by = ", moderator, ", fx = ", fx, ")",
        if(rhs_base != "") paste0(" + ", rhs_base)
      )

    } else if (interaction == "ordered") {

      f_str <- paste0(
        ph,
        " ~ s(", predictor, ", k = ", k, ", fx = ", fx, ")",
        " + ", moderator,
        " + s(", predictor, ", k = ", k,
        ", by = ", moderator, ", fx = ", fx, ")",
        if(rhs_base != "") paste0(" + ", rhs_base)
      )

    } else {
      stop("interaction must be 'continuous', 'categorical', or 'ordered'")
    }

    f <- as.formula(dedupe_formula_terms(f_str))

    mgcv::bam(
      f,
      data = data,
      method = method,
      discrete = TRUE
    )
  })
  
  names(models) <- phenotypes
  attr(models, "predictor") <- predictor
  attr(models, "moderator") <- moderator
  attr(models, "fx") <- fx
  attr(models, "k") <- k
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
