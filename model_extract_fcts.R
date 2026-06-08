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

# extract_bam_features_main <- function(models, data, predictor, outcome_col_map = NULL){
#   tab <- extract_bam_table(models)
#   
#   feats <- lapply(names(models), function(ph){
#     m_full <- models[[ph]]$full
#     m_nolin <- models[[ph]]$nolin
#     
#     # partial R2 (full vs reduced)
#     r2 <- 1 - deviance(m_full)/deviance(m_nolin)
#     
#     # average first derivative (finite difference on smooth)
#     nd <- 50
#     rng <- range(data[[predictor]], na.rm = TRUE)
#     grid <- seq(rng[1], rng[2], length.out = nd)
#     
#     pd <- predict(m_full, newdata = setNames(data.frame(grid), predictor), type = "lpmatrix")
#     co <- coef(m_full)
#     
#     fit <- as.numeric(pd %*% co)
#     d1 <- diff(fit)/diff(grid)
#     avg_deriv <- mean(d1, na.rm = TRUE)
#     
#     # CI-based region of effect (approx via SE of smooth)
#     pr <- predict(m_full, newdata = setNames(data.frame(grid), predictor), se.fit = TRUE)
#     lower <- pr$fit - 1.96 * pr$se.fit
#     upper <- pr$fit + 1.96 * pr$se.fit
#     
#     sig <- !(lower < 0 & upper > 0)
#     
#     if(any(sig)){
#       region <- c(min(grid[sig]), max(grid[sig]))
#     } else region <- c(NA, NA)
#     
#     c(partial_r2 = r2,
#       avg_first_deriv = avg_deriv,
#       region_start = region[1],
#       region_end = region[2])
#   })
#   
#   invisible(cbind(tab, do.call(rbind, feats)))
# }

extract_bam_interaction_table <- function(models){
  tab <- as.data.frame(t(sapply(models, function(m){
    s <- summary(m)$s.table
    if(attr(m, "interaction") == "nonlinear"){
      c(pred_edf = s[1,1], pred_p = s[1,4],
        mod_edf = s[2,1], mod_p = s[2,4],
        int_edf = s[3,1], int_p = s[3,4])
    } else {
      co <- summary(m)$p.table
      c(pred_beta = co[2,1], pred_p = co[2,4],
        mod_beta = co[3,1], mod_p = co[3,4],
        int_beta = co[4,1], int_p = co[4,4])
    }
  })))
  
  num_cols <- names(tab)[grepl("_p$", names(tab))]
  for(nm in num_cols) tab[[paste0(nm, "_adj")]] <- p.adjust(tab[[nm]], "BH")
  invisible(tab)
}

extract_bam_features_main <- function(models, data, predictor, covariates = NULL, random_effect = NULL){
  
  tab <- extract_bam_table(models, predictor)
  
  feats <- lapply(names(models), function(ph){
    
    m_full  <- models[[ph]]$full
    m_nolin <- models[[ph]]$nolin
    
    # rebuild reduced model WITHOUT predictor entirely
    rhs_cov <- c(covariates,
                 if(!is.null(random_effect)) paste0("s(", random_effect, ", bs='re')") else NULL)
    rhs_cov <- paste(rhs_cov[!is.na(rhs_cov)], collapse = " + ")
    
    f_red <- as.formula(paste(ph, "~", rhs_cov))
    
    # subset data to only include predictor is not NA to fit reduced model
    data_fit <- data %>% filter(!is.na(get(predictor)))
    m_red <- mgcv::bam(f_red, data = data_fit, method = "fREML", discrete = TRUE)
    
    # partial R2: full predictor contribution vs no predictor model
    r2_full_vs_none <- 1 - deviance(m_full)/deviance(m_red)

    # partial R2: spline contribution over linear model
    r2_spline_vs_linear <- 1 - deviance(m_full)/deviance(m_nolin)
    
    # average first derivative
    rng <- range(data[[predictor]], na.rm = TRUE)
    grid <- seq(rng[1], rng[2], length.out = 50)
    newd <- data.frame(tmp = grid)
    names(newd)[1] <- predictor
    
    # covariates fixed at typical values
    if(!is.null(covariates)){
      for(v in covariates){
        if(is.numeric(data[[v]])){
          newd[[v]] <- mean(data[[v]], na.rm = TRUE)
        } else {
          newd[[v]] <- as.character(stats::na.omit(data[[v]])[1])
        }
      }
    }
    
    # random effect fixed at a real observed level
    if(!is.null(random_effect)){
      newd[[random_effect]] <- sample(unique(data[[random_effect]]), 1)
    }
    
    
    fit <- predict(m_full, newdata = newd, type = "response")
    deriv <- diff(fit)/diff(grid)
    avg_deriv <- mean(deriv, na.rm = TRUE)
    
    # region of effect (CI crossing rule)
    pr <- predict(m_full, newdata = newd, se.fit = TRUE)
    lower <- pr$fit - 1.96 * pr$se.fit
    upper <- pr$fit + 1.96 * pr$se.fit
    
    sig <- !(lower < 0 & upper > 0)
    region <- if(any(sig)) c(min(grid[sig]), max(grid[sig])) else c(NA, NA)
    
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


extract_bam_features_interaction <- function(models, data, predictor, moderator){
  tab <- extract_bam_interaction_table(models)
  
  feats <- lapply(names(models), function(ph){
    m <- models[[ph]]
    is_nl <- attr(m, "interaction") == "nonlinear"
    
    if(is_nl){
      m_full <- m
      
      # reduced = drop interaction term (re-fit via update)
      m_red <- update(m_full, . ~ . - ti(get(predictor), get(moderator)))
      
      r2 <- 1 - deviance(m_full)/deviance(m_red)
      
      # interaction surface evaluated on grid
      gx <- seq(-2, 2, length.out = 20)
      gy <- seq(-2, 2, length.out = 20)
      
      grid <- expand.grid(x = gx, y = gy)
      names(grid) <- c(predictor, moderator)
      
      pr <- predict(m_full, newdata = grid, se.fit = TRUE)
      
      lower <- pr$fit - 1.96 * pr$se.fit
      upper <- pr$fit + 1.96 * pr$se.fit
      
      sig <- !(lower < 0 & upper > 0)
      
      region <- if(any(sig)){
        c(min(grid[[predictor]][sig]), max(grid[[predictor]][sig]),
          min(grid[[moderator]][sig]), max(grid[[moderator]][sig]))
      } else rep(NA, 4)
      
      c(partial_r2 = r2,
        region_x_start = region[1],
        region_x_end = region[2],
        region_y_start = region[3],
        region_y_end = region[4])
      
    } else {
      
      m_full <- m
      m_red <- update(m_full, . ~ . - get(predictor):get(moderator))
      
      r2 <- 1 - deviance(m_full)/deviance(m_red)
      
      # slopes at -1SD / mean / +1SD of moderator
      mod_vals <- scale(model.frame(m_full)[[moderator]], scale = TRUE)
      ref <- c(-1, 0, 1)
      
      slopes <- sapply(ref, function(z){
        coef(m_full)[predictor] + coef(m_full)[paste0(predictor, ":", moderator)] * z
      })
      
      c(partial_r2 = r2,
        slope_m1sd = slopes[1],
        slope_mean = slopes[2],
        slope_p1sd = slopes[3])
    }
  })
  
  cbind(tab, do.call(rbind, feats))
}

extract_bam_features_interaction_ifsig <- function(models, table, alpha = 0.05, predictor, moderator){
  
  feats <- lapply(names(models), function(ph){
    
    # interaction p (robust lookup)
    row <- table[ph, , drop = FALSE]
    p_int <- row[[grep("_p_adj$", names(row))]]
    if(length(p_int) == 0) p_int <- NA
    
    # gate: only proceed if significant
    if(is.na(p_int) || p_int >= alpha){
      return(c(partial_r2 = NA,
               slope_m1sd = NA,
               slope_mean = NA,
               slope_p1sd = NA,
               region_x_start = NA,
               region_x_end = NA,
               region_y_start = NA,
               region_y_end = NA))
    }
    
    m <- models[[ph]]
    
    is_nl <- attr(m, "interaction") == "nonlinear"
    
    if(is_nl){
      
      m_full <- m
      m_red <- update(m_full, . ~ . - ti(get(predictor), get(moderator)))
      
      r2 <- 1 - deviance(m_full)/deviance(m_red)
      
      gx <- seq(-2, 2, length.out = 20)
      gy <- seq(-2, 2, length.out = 20)
      
      grid <- expand.grid(x = gx, y = gy)
      names(grid) <- c(predictor, moderator)
      
      pr <- predict(m_full, newdata = grid, se.fit = TRUE)
      
      lower <- pr$fit - 1.96 * pr$se.fit
      upper <- pr$fit + 1.96 * pr$se.fit
      
      sig <- !(lower < 0 & upper > 0)
      
      region <- if(any(sig)){
        c(min(grid[[predictor]][sig]), max(grid[[predictor]][sig]),
          min(grid[[moderator]][sig]), max(grid[[moderator]][sig]))
      } else rep(NA, 4)
      
      return(c(partial_r2 = r2,
               region_x_start = region[1],
               region_x_end = region[2],
               region_y_start = region[3],
               region_y_end = region[4]))
    }
    
    # linear interaction case
    m_full <- m
    m_red <- update(m_full, . ~ . - get(predictor):get(moderator))
    
    r2 <- 1 - deviance(m_full)/deviance(m_red)
    
    mod_frame <- model.frame(m_full)
    z <- scale(mod_frame[[moderator]], scale = TRUE)
    
    slopes <- sapply(c(-1,0,1), function(v){
      coef(m_full)[predictor] + coef(m_full)[paste0(predictor, ":", moderator)] * v
    })
    
    c(partial_r2 = r2,
      slope_m1sd = slopes[1],
      slope_mean = slopes[2],
      slope_p1sd = slopes[3])
  })
  
  do.call(rbind, feats) |> as.data.frame()
}

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
