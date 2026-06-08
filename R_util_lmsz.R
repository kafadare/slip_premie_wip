## LMSz Utility functions written by Ben Jung, PhD
library(dplyr)
library(magrittr)
library(readr)
library(stringr)
library(tidyr)
library(tibble)
library(glue)
library(ggseg)
library(gridExtra)
#library(ggsegTissue)
#library(kableExtra)
library(table1)
#library(flextable)
#library(broom)
#library(broom.mixed)
#library(effectsize)
#library(ggpubr)
library(cowplot)
#library(reticulate)

set.seed(42)

source("/mnt/isilon/bgdlab_processing/braincharts/SLIP/2025_03/code/gamlss/102.gamlss-recode.r")
source("/mnt/isilon/bgdlab_processing/braincharts/SLIP/2025_03/code/gamlss/301.functions.r")

# Define default color scheme
color_scheme <- c("royalblue","white","firebrick")

# Define global measures
global_measures <- c("eTIV", "GMV", "sGMV", "WMV", "Ventricles", "Cerebellum.all", 
                     "totalSA", "meanCT")

global_measure_names <- c("ICV", "GMV", "sGMV", "WMV", "Ventricles", "Cerebellum", 
                          "SA", "CT")


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


#function adapted for gestational age inclusion in LMSz EK March 2026
backtransform_centile_fans <- function(pheno, lmsz_models_path, sim_list, df) {
  # Load data
  pheno_name <- pheno
  raw_name <- glue("{pheno}.raw")
  zscore_name <- glue("{pheno}.zscore")
  
  #old_pheno <- slip_features$old_pheno[slip_features$new_pheno == pheno]
  #if (length(old_pheno) == 0) {
  #return()
  #}
  growthChartModel <- readRDS(glue("{lmsz_models_path}lmsz_{pheno}.rds"))
  # Load SLIP fit model
  if (grepl("Thick",pheno)) {
    fit_file <- glue("{models_path}mpr-{pheno}Transformed/FIT.EXTRACT.rds")
  } else {
    fit_file <- glue("{models_path}median-{pheno}Transformed/FIT.EXTRACT.rds")
  }
  
  orig_fit <- readRDS(fit_file)
  
  df_tmp <- df %>%
    select(participant_id, adjusted_age_days_log, sex, gestational_age, !!pheno, !!raw_name) %>%
    rename(pheno = !!pheno,
           normalised_value = !!raw_name) %>%
    na.omit()
  
  # Simulate centiles
  lmsz_sim <- centile_predict(gamlssModel = growthChartModel, 
                              sim_df_list = sim_list, 
                              x_var = "adjusted_age_days_log", 
                              desiredCentiles = c(0.025,0.250,0.50,0.750,0.975),
                              df = df_tmp,
                              average_over = FALSE)
  
  # Reformat simulated values
  
  lmsz_sim <- lapply(names(lmsz_sim), function(x){ 
    out <- lmsz_sim[[x]] %>% pivot_longer(-c(adjusted_age_days_log),
                                            names_to = "centile_cat",
                                            values_to = "value")
    out$sex <- ifelse(grepl("Female", x), "Female", "Male")
    out$gestational_age <- as.numeric(gsub("\\D+", "", x))
    return(out)
  })
  
  lmsz_df <- bind_rows(lmsz_sim) %>%
    mutate(site = NA,
           sex = factor(sex, levels = c("Female","Male")),
           gestational_age = as.numeric(gestational_age),
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
           .after = site) %>%
    select(-c(value))
  
  # Plot original SLIP Model
  NEWData <- lmsz_df %>%
    select(adjusted_age_days_log, sex) %>% 
    rename(AgeTransformed = adjusted_age_days_log) %>%
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
    select(AgeTransformed, sex, centile_cat, SLIP_input_fans) %>%
    rename(adjusted_age_days_log = AgeTransformed)
  
  # Plot backtransformed centiles
  NEWData <- lmsz_df %>%
    select(adjusted_age_days_log, sex, gestational_age,
           centile_cat,
           value = LMSz_input_centile) %>% 
    rename(AgeTransformed = adjusted_age_days_log) %>%
    distinct() %>%
    as.data.frame()
  
  backtransformed_fans <- backtransform_centiles(NEWData, orig_fit$param)
  
  backtransformed_fans <- backtransformed_fans %>%
    rename(adjusted_age_days_log = AgeTransformed) %>%
    select(adjusted_age_days_log, sex, centile_cat, gestational_age,
           LMSz_fans = PRED.pop)
  
  # Join all results
  out_df <- lmsz_df %>%
    inner_join(original_fans, by= join_by(adjusted_age_days_log, centile_cat, sex)) %>%
    inner_join(backtransformed_fans, by= join_by(adjusted_age_days_log, centile_cat, sex, gestational_age)) %>%
    mutate(pheno = pheno, model = "SLIP", .before = adjusted_age_days_log) %>%
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


lmsz_modeling <- function(df_subset, m_models, s_models, vars, models_path, save_path, save_model = TRUE, reweight = FALSE) {
  if (!dir.exists(save_path)) dir.create(save_path)
  source("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/wp_taki_EK.R")
  
  #create empty list for wormplots
  worm_plots <- list()
  # Create tibble to hold best model
  models_df <- tibble(feature = full_vars,
                      m_model = "",
                      s_model = "",
                      AIC = NA,
                      BIC = 999999999999,
                      age_coef = NA,
                      sex_coef = NA,
                      intercept_coef = NA,
                      interaction_coef = NA,
                      age_sigma_coef = NA,
                      sex_sigma_coef = NA,
                      intercept_sigma_coef = NA,
                      interaction_sigma_coef = NA)
  
  # Build the growth chart model
  p <- "zscore" # specify the phenotype
  for (pheno in full_vars) {
    # Define phenotype
    print(pheno)
    #pheno_name <- glue("{pheno}Transformed")
    zscore_name <- glue("{pheno}.zscore")
    # Load reference SLIP fit model
    if (grepl("Thick",pheno)) {
      orig_fit <- readRDS(glue("{models_path}mpr-{pheno}Transformed/FIT.EXTRACT.rds"))
    } else {
      orig_fit <- readRDS(glue("{models_path}median-{pheno}Transformed/FIT.EXTRACT.rds"))
    }
    
    # Make subset
    df_tmp <- df_subset %>%
      select(participant_id, adjusted_age_days_log, sex, 
             adjusted_age_years,
             !!pheno, !!zscore_name) %>%
      rename(zscore = !!zscore_name) %>%
      na.omit()
    
    # Loop through each mu and sigma model
    mod <- models_df$feature == pheno
    for (m in m_models) {
      for (s in s_models) {
        print(m)
        print(s)
        
        # Run preliminary model - all subjects
        growthChartModel_all <- gamlss(formula = as.formula(glue("zscore {m}")),
                                       sigma.formula = as.formula(glue("zscore {s}")),
                                       family = NO,
                                       data = df_tmp,
                                       control = gamlss.control(n.cyc = 200),  # See (2)
                                       trace = F)
        if(reweight == TRUE){
          # Re-weight observations to exclude individuals with large residuals (> 3.5)
          df_tmp$weight_tmp <- rep(1, dim(df_tmp)[1])
          df_tmp$weight_tmp[which(abs(resid(growthChartModel_all)) > 3.5)] <- 0
          
          print(glue("The Number of Points Weighted 0 due to large residuals > 3.5 for {pheno} is {sum(df_tmp$weight_tmp == 0)}"))
          
          # Run refined model
          growthChartModel <- gamlss(formula = as.formula(glue("zscore {m}")),
                                     sigma.formula = as.formula(glue("zscore {s}")),
                                     family = NO,
                                     data = df_tmp,
                                     weights = df_tmp$weight_tmp,
                                     control = gamlss.control(n.cyc = 500),
                                     trace = F)
        } else{
          growthChartModel = growthChartModel_all
        }
        # Check if the model has the lowest BIC and record it
        if (growthChartModel$sbc <= models_df$BIC[mod]) {
          fin_model <- growthChartModel
          models_df$m_model[mod] <- m
          models_df$s_model[mod] <- s
          models_df$BIC[mod] <- growthChartModel$sbc
        }
        print("--------------")
      }
    }
    
    # Refine formatting of final models_df
    models_df$m_model[mod] <- format(fin_model$mu.formula)
    models_df$s_model[mod] <- format(fin_model$sigma.formula)
    
    # Extract model coefficients - mu
    models_df$intercept_coef[mod]   <- fin_model$mu.coefficients["(Intercept)"]
    models_df$age_coef[mod]         <- fin_model$mu.coefficients["adjusted_age_days_log"]
    models_df$sex_coef[mod]         <- fin_model$mu.coefficients["sexMale"]
    models_df$interaction_coef[mod] <- fin_model$mu.coefficients["adjusted_age_days_log:sexMale"]
    # Extract model coefficients - sigma
    models_df$intercept_sigma_coef[mod]   <- fin_model$sigma.coefficients["(Intercept)"]
    models_df$age_sigma_coef[mod]         <- fin_model$sigma.coefficients["adjusted_age_days_log"]
    models_df$sex_sigma_coef[mod]         <- fin_model$sigma.coefficients["sexMale"]
    models_df$interaction_sigma_coef[mod] <- fin_model$sigma.coefficients["adjusted_age_days_log:sexMale"]
    
    # Extract final AIC/BIC
    models_df$AIC[mod] <- fin_model$aic
    models_df$BIC[mod] <- fin_model$sbc
    
    # Save final model
    if(save_model == TRUE){
      saveRDS(fin_model,glue("{save_path}lmsz_{pheno}.rds"))
    }
    # Assign worm plot
    worm_plots[[pheno]] <- wp.taki_EK(object = fin_model, df_tmp = df_tmp) +
      ggtitle(paste0(pheno," LMSz")) +
      theme_classic(base_size = 12) +
      theme(plot.title = element_text(size = 12),
            axis.title = element_text(size = 12))
    
    # Update LMSz quantiles
    df_tmp <- df_tmp %>%
      add_column(!!glue("{pheno}.lmsz") := pred_og_centiles(fin_model,
                                                            df_tmp)) %>%
      select(-c(pheno, zscore, adjusted_age_years))
    
    # Add LMSz centile to final dataframe
    df_subset <- df_subset %>% dplyr::left_join(df_tmp, by = join_by(participant_id, sex, adjusted_age_days_log))
    
    
    print("-------------------------")
  }
  return_list <- list(models = models_df, df = df_subset, wp = worm_plots)
  return(return_list)
}
