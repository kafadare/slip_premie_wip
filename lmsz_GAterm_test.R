# LMSz for full sample test GA term
# Load libraries
require(dplyr)
require(magrittr)
require(readr)
require(stringr)
require(tidyr)
require(gamlssTools) # Margaret's package
require(gamlss.add)
require(tibble)
require(glue)
# For figures/plotting
require(ggplot2)
require(patchwork)
# For figures/plotting
require(gridExtra)
require(ggsegDKT)
require(ggseg)
# require(ggpubr)
require(ggcorrplot)
require(cowplot)
require(egg)
require(interactions)
require(growthstandards)
require(gtsummary)
require(svglite)

set.seed(42)
params <- list(
  output_folder = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/lmsz_median_Th-MPR_25.11_GAtest_splines/",
  centiles_data_path = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/median-Th-MPR-allgrades-25.11/data_clean/centiles_distinct_2026-03-18.csv", 
  raw_data_path = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11/raw_csv/slip_median_th_mpr_2026-03-17.csv",
  models_path = "/mnt/isilon/bgdlab_processing/braincharts/SLIP/2025_03/code/gamlss/RDS/"
)
save_model = TRUE
reweight = FALSE

# set folder paths
code_path <- "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/"
vars_path <-  "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/vars/"
output_folder <- params$output_folder
centiles_data_path <- params$centiles_data_path
raw_data_path <- params$raw_data_path
models_path <- params$models_path
# output folder
save_path <- params$output_folder
if (!dir.exists(save_path)) dir.create(save_path)

##  output subfolder paths
csv_save_path <-  paste0(save_path, "results_csv/")
if (!dir.exists(csv_save_path)) dir.create(csv_save_path)
plot_save_path <- paste0(save_path, "plots/")
if (!dir.exists(plot_save_path)) dir.create(plot_save_path)

panel_save <- function(plot, filename, extension, folder, plot_width = 8, plot_height = 10){
  save_name <- paste0(folder, filename, ".", extension)
  ggsave(save_name, plot = plot, dpi = 300, width = plot_width, height = plot_height)
}

source(glue("{code_path}R_util_lmsz.R"))
source("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/wp_taki_EK.R")


full_vars <- read.csv(paste0(vars_path,"all_vars.csv")) %>% pull(x)
global_vars <- read.csv(paste0(vars_path,"global_vars.csv")) %>% pull(x)
vol_vars <- read.csv(paste0(vars_path,"vol_vars.csv")) %>% pull(x)
sa_vars <- read.csv(paste0(vars_path,"sa_vars.csv")) %>% pull(x)
th_vars <- read.csv(paste0(vars_path,"th_vars.csv")) %>% pull(x)

models_save_path <- paste0(save_path,"lmsz_models/")
if (!dir.exists(models_save_path)) dir.create(models_save_path)

# *******

# Read in centiles and raw data
df <- read.csv(centiles_data_path) %>% select(-"X")
df_raw <- read.csv(raw_data_path)%>% select(any_of(c(full_vars, "participant_id", "session_id", "sex", "age_days",
                                                      "site", "avg_grade", "age_at_scan"))) %>%
rename_with(~paste0(.x, ".raw"), all_of(full_vars))

shared_names <- intersect(names(df_raw), names(df))

df <- merge(df, df_raw, by = shared_names, all.x = FALSE, all.y = FALSE)

df <- df %>% filter(adjusted_age_in_years <= 21)

## Set 0/1 centiles to min/max value to avoid Inf/-Inf
median_th_mpr_flag <- any(grepl(".median", names(df)))

if(median_th_mpr_flag){
th_vars_median <- paste0(th_vars, ".median")

complete_centile_vars <- c(full_vars, th_vars_median)
} else{
  complete_centile_vars <- c(full_vars)
}

for (x in complete_centile_vars){
  df[df[,x] == 0 & !is.na(df[,x]), x] <- 0.0000000000000000000000001
  df[df[,x] == 1 & !is.na(df[,x]), x] <- 0.99999999999999994
}

#zscore centiles and other variables of interest
df <- df %>%
  mutate(across(all_of(full_vars), qnorm,
                .names = "{.col}.zscore"))
#set sex and preterm as factor
df$sex <- as.factor(df$sex)
df$preterm <- as.factor(df$preterm)
zscore_vars_names <- paste0(full_vars, ".zscore")

## find Inf cells, make a table of N per phenotype, should be 0
inf_nan_table <- as.data.frame(t(sapply(df[zscore_vars_names], function(x){
  c(inf = sum(is.infinite(x)), nan = sum(is.nan(x)))
})))
max(inf_nan_table$inf)
max(inf_nan_table$nan)

inf_nan_table %>% filter(inf > 0)

# Specify Terms and create all combinations for models to test
lmsz_vars <- list(
  ga = c(NA, "gestational_age", "ns(gestational_age, 2)", "ns(gestational_age, 3)"),
  age = c(NA, "adjusted_age_in_days", "ns(adjusted_age_in_days, 2)", "ns(adjusted_age_in_days, 3)"),
  sex = c(NA, "sex")
)
  
  # all combinations of term versions
  main_grid <- expand.grid(lmsz_vars, stringsAsFactors = FALSE)
  
  formulas <- list()
  
  for (i in seq_len(nrow(main_grid))) {
    
    row <- as.list(main_grid[i, ])
    terms <- row[!is.na(row)]
    main_terms <- unname(terms)
    var_names  <- names(terms)
    
    base <- paste(main_terms, collapse = " + ")
    formulas[[length(formulas) + 1]] <- paste0("~ ", base)
    
    # two-way interactions
    if (length(terms) >= 2) {
      pair_ids <- combn(seq_along(terms), 2, simplify = FALSE)
      
      # get each valid combination of two vars in the model
      for (k in seq_along(pair_ids)) {
        pair_subsets <- combn(pair_ids, k, simplify = FALSE)
        
        # choose 1..n (as many as possible) interaction terms
        for (subset_pairs in pair_subsets) {
          int_terms <- vapply(subset_pairs,
            function(x) paste0(main_terms[x[1]], ":", main_terms[x[2]]),
            character(1) )
          mod <- paste(c(main_terms, int_terms), collapse = " + ")
          formulas[[length(formulas) + 1]] <- paste("~ ", mod)
        }
      }
    }
  }
  names(formulas) <- paste0("mod", seq_along(formulas))
  
  formulas_df <- as.data.frame(do.call(rbind,formulas)) %>%
    rename(formula = V1)
  
  formulas_df[1,1] <- '~ 0'
  formulas_df[nrow(formulas_df)+1, 1] <- '~ 1'
  

#' # Specify LMSz model search space
#' m_models = 
#' 
#' # Variance (sigma) models
#' s_models = 

  
df_subset <- df %>% filter(!is.na(gestational_age))
 
 #create empty list for wormplots
 worm_plots_ga <- list()
 # Create tibble to hold best model
 models_ga <- tibble(feature = full_vars,
                     m_model = "",
                     s_model = "",
                     AIC = NA,
                     BIC = 999999999999,
                     # mu coefs
                     intercept_mu_coef = NA,
                     ga_mu_linear_coef = NA,
                     ga_mu_ns2_edf = NA,
                     ga_mu_ns3_edf = NA,
                     age_mu_ns2_edf = NA,
                     age_mu_ns3_edf = NA,
                     ga.sex_mu_interaction = NA,
                     ga.age_mu_interaction = NA,
                     # sigma coefs
                     intercept_sigma_coef = NA,
                     ga_sigma_linear_coef = NA,
                     ga_sigma_ns2_edf = NA,
                     ga_sigma_ns3_edf = NA,
                     age_sigma_ns2_edf = NA,
                     age_sigma_ns3_edf = NA,
                     ga.sex_sigma_interaction = NA,
                     ga.age_sigma_interaction = NA
                     )
 
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
            adjusted_age_in_years, gestational_age,
            !!pheno, !!zscore_name) %>%
     rename(zscore = !!zscore_name) %>%
     na.omit() %>% 
     mutate(gestational_age = scale(gestational_age))
  
   # Loop through each mu and sigma model
   mod <- models_ga$feature == pheno
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
       if (growthChartModel$sbc <= models_ga$BIC[mod]) {
         fin_model <- growthChartModel
         models_ga$m_model[mod] <- m
         models_ga$s_model[mod] <- s
         models_ga$BIC[mod] <- growthChartModel$sbc
       }
      print("--------------")
     }
   }
   
   
   # Refine formatting of final models_ga
   models_ga$m_model[mod] <- format(fin_model$mu.formula)
   models_ga$s_model[mod] <- format(fin_model$sigma.formula)
   
   # Extract model coefficients - mu
   models_ga$intercept_coef[mod]   <- fin_model$mu.coefficients["(Intercept)"]
   models_ga$ga_coef[mod]         <- fin_model$mu.coefficients["gestational_age"]
   models_ga$ga_ns2_coef1[mod] <- fin_model$mu.coefficients["ns(gestational_age, df=2)"][1]
   models_ga$ga_ns2_coef2[mod] <- fin_model$mu.coefficients["ns(gestational_age, df=2)"][2]
   models_ga$ga_ns3_coef1[mod] <- fin_model$mu.coefficients["ns(gestational_age, df=3)"][1]
   models_ga$ga_ns3_coef2[mod] <- fin_model$mu.coefficients["ns(gestational_age, df=3)"][2]
   models_ga$ga_ns3_coef3[mod] <- fin_model$mu.coefficients["ns(gestational_age, df=3)"][3]
   models_ga$ga.sex_interaction_coef[mod] <- fin_model$mu.coefficients["gestational_age:sexMale"]
   models_ga$ga.age_interaction_coef[mod] <- fin_model$mu.coefficients["gestational_age:adjusted_age_days_log"]
   # Extract model coefficients - sigma
   models_ga$intercept_sigma_coef[mod]   <- fin_model$sigma.coefficients["(Intercept)"]
   models_ga$ga.sex_interaction_sigma_coef[mod] <- fin_model$sigma.coefficients["gestational_age:sexMale"]
   models_ga$ga.age_interaction_sigma_coef[mod] <- fin_model$sigma.coefficients["gestational_age:adjusted_age_days_log"]
   
   # Extract final AIC/BIC
   models_ga$AIC[mod] <- fin_model$aic
   models_ga$BIC[mod] <- fin_model$sbc
 
   
   # Save final model
   if(save_model == TRUE){
     saveRDS(fin_model,glue("{models_save_path}lmsz_{pheno}.rds"))
   }
   # Assign worm plot
   worm_plots_ga[[pheno]] <- wp.taki_EK(object = fin_model, df_tmp = df_tmp, xlim.worm = 4) +
     ggtitle(paste0(pheno," LMSz")) +
     theme_classic(base_size = 12) +
     theme(plot.title = element_text(size = 12),
           axis.title = element_text(size = 12))
  
   # Update LMSz quantiles
   df_tmp <- df_tmp %>%
    add_column(!!glue("{pheno}.lmsz") := pred_og_centiles(fin_model,
                                                           df_tmp)) %>%
     select(-c(pheno, zscore, adjusted_age_in_years, gestational_age))
   
   # Add LMSz centile to final dataframe
   df_subset <- df_subset %>% dplyr::left_join(df_tmp, 
                                               by = join_by(participant_id, sex, 
                                                            adjusted_age_days_log))
  
   
   print("-------------------------")
 }

 df_ga <- df_subset
 
write_csv(models_ga, glue("{csv_save_path}lmsz_models_ga.csv"))
write_csv(df_ga, glue("{csv_save_path}lmsz_adjusted_centiles_ga.csv"))
saveRDS(worm_plots_ga,glue("{csv_save_path}worm_plots_ga.rds"))


# Calculate Centile Fans for both SLIP (original) and LMSz-GA

#out_fans_slip <- read.csv(glue("{csv_save_path}lmsz_centfans_ga_uncorrected.csv"))
#df_ga <-  read.csv(glue("{csv_save_path}lmsz_adjusted_centiles_ga.csv"))
#models_ga <- read.csv(glue("{csv_save_path}lmsz_models_ga.csv"))
#readRDS(glue("{csv_save_path}worm_plots_ga.rds"))

##simulation frames, log days age, sex, and ga
#between 280 days (40 weeks gestation) and 21 years
sim_list_sexga <- list()

sim_list_sexga$Male22week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                         sex = c("Male"), gestational_age = 22)
sim_list_sexga$Female22week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                           sex = c("Female"), gestational_age = 22)

sim_list_sexga$Male23week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                         sex = c("Male"), gestational_age = 23)
sim_list_sexga$Female23week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                           sex = c("Female"), gestational_age = 23)

sim_list_sexga$Male24week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                  sex = c("Male"), gestational_age = 24)
sim_list_sexga$Female24week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                    sex = c("Female"), gestational_age = 24)

sim_list_sexga$Male25week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                         sex = c("Male"), gestational_age = 25)
sim_list_sexga$Female25week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                           sex = c("Female"), gestational_age = 25)

sim_list_sexga$Male26week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                         sex = c("Male"), gestational_age = 26)
sim_list_sexga$Female26week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                           sex = c("Female"), gestational_age = 26)

sim_list_sexga$Male27week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                         sex = c("Male"), gestational_age = 27)
sim_list_sexga$Female27week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                           sex = c("Female"), gestational_age = 27)

sim_list_sexga$Male28week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                        sex = c("Male"), gestational_age = 28)
sim_list_sexga$Female28week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                          sex = c("Female"), gestational_age = 28)

sim_list_sexga$Male29week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                         sex = c("Male"), gestational_age = 29)
sim_list_sexga$Female29week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                           sex = c("Female"), gestational_age = 29)

sim_list_sexga$Male30week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                         sex = c("Male"), gestational_age = 30)
sim_list_sexga$Female30week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                           sex = c("Female"), gestational_age = 30)

sim_list_sexga$Male31week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                         sex = c("Male"), gestational_age = 31)
sim_list_sexga$Female31week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                           sex = c("Female"), gestational_age = 31)

sim_list_sexga$Male32week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                        sex = c("Male"), gestational_age = 32)
sim_list_sexga$Female32week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                          sex = c("Female"), gestational_age = 32)

sim_list_sexga$Male33week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                         sex = c("Male"), gestational_age = 33)
sim_list_sexga$Female33week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                           sex = c("Female"), gestational_age = 33)

sim_list_sexga$Male34week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                         sex = c("Male"), gestational_age = 34)
sim_list_sexga$Female34week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                           sex = c("Female"), gestational_age = 34)

sim_list_sexga$Male35week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                         sex = c("Male"), gestational_age = 35)
sim_list_sexga$Female35week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                           sex = c("Female"), gestational_age = 35)

sim_list_sexga$Male36week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                        sex = c("Male"), gestational_age = 36)
sim_list_sexga$Female36week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                          sex = c("Female"), gestational_age = 36)

sim_list_sexga$Male37week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                         sex = c("Male"), gestational_age = 37)
sim_list_sexga$Female37week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                           sex = c("Female"), gestational_age = 37)

sim_list_sexga$Male38week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                         sex = c("Male"), gestational_age = 38)
sim_list_sexga$Female38week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                           sex = c("Female"), gestational_age = 38)

sim_list_sexga$Male39week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                         sex = c("Male"), gestational_age = 39)
sim_list_sexga$Female39week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                           sex = c("Female"), gestational_age = 39)

sim_list_sexga$Male40week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                        sex = c("Male"), gestational_age = 40)
sim_list_sexga$Female40week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                          sex = c("Female"), gestational_age = 40)

sim_list_sexga$Male41week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                         sex = c("Male"), gestational_age = 41)
sim_list_sexga$Female41week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                           sex = c("Female"), gestational_age = 41)

sim_list_sexga$Male42week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                         sex = c("Male"), gestational_age = 42)
sim_list_sexga$Female42week <- expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                                           sex = c("Female"), gestational_age = 42)

for (pheno in full_vars) {
  print(pheno)
  tmp <-  backtransform_centile_fans(pheno, 
                                     models_save_path, 
                                     sim_list_sexga, 
                                     df_ga)
  if (pheno == full_vars[1]) {
    out_fans_slip <- tmp
  } else {
    out_fans_slip <- bind_rows(out_fans_slip, tmp)
  }
}

write_csv(out_fans_slip, glue("{csv_save_path}lmsz_centfans_ga_uncorrected.csv"))
  
# **************

# Visualization 

out_fans_slip <- read.csv(glue("{csv_save_path}lmsz_centfans_ga_uncorrected.csv"))
df_ga <-  read.csv(glue("{csv_save_path}lmsz_adjusted_centiles_ga.csv"))
models_ga <- read.csv(glue("{csv_save_path}lmsz_models_ga.csv"))
#readRDS(glue("{csv_save_path}worm_plots_ga.rds"))

low_color = "#4878CF"
mid_color = "white"
high_color = "#D65F5F"

#add label column to be able to plot with ggseg
models_ga$domain <- ifelse(grepl("global", models_ga$feature), "global",
                            ifelse(grepl("GrayVol", models_ga$feature), "vol",
                                   ifelse(grepl("SurfArea", models_ga$feature), "sa",
                                          ifelse(grepl("ThickAvg", models_ga$feature), "th",
                                                 ifelse(grepl("aseg", models_ga$feature), "aseg", NA)))))
cols_edit <- grepl("aseg(?!_CC)", models_ga$feature, perl = TRUE)

models_ga$label <- models_ga$feature %>% 
  sub("aparc_GrayVol_", "", .) %>%
  sub("aparc_SurfArea_", "", .) %>%
  sub("aparc_ThickAvg_", "", .) %>%
  sub("aseg_", "", .) %>%
  sub("X", "x", .) %>%
  sub("rd_Ventricle", "rd_ventricle", .) %>%
  sub("th_Ventricle", "th_ventricle", .) %>%
  sub("Thalamus", "Thalamus-Proper", .)

models_ga$label[cols_edit] <- gsub("_", "-", models_ga$label[cols_edit])

# Look at which features have gestational age ns() in them


# Plot Coefficient of Gestational Age in Mu
# merge with atlas data
plot_data_dkt <- merge(models_ga, ggseg::dk, by = "label")
plot_data_aseg <- merge(models_ga, ggseg::aseg, by = "label")

# Regions not matched in atlases
#print(models_ga_noplot <- models_ga[!(models_ga$label %in% plot_data_dkt$label) & !(models_ga$label %in% plot_data_aseg$label),"label"], n = 100)

# plots
(ga_vol_mu_age_dk_plot <- ggplot(plot_data_dkt %>% filter(domain == "vol")) +
    geom_brain(atlas = dk, 
               position = position_brain(hemi ~ side),
               aes(fill = age_coef)) +
    scale_fill_gradient2(low = low_color, mid = mid_color, high = high_color, limits = c(-0.35, 0.35)) +
    theme_void()+
    labs(title = "ga Volume Mu GA"))

(ga_vol_mu_age_aseg_plot <- ggplot(plot_data_aseg %>% filter(domain == "aseg")) +
    geom_brain(atlas = aseg,
               aes(fill = age_coef)) +
    scale_fill_gradient2(low = low_color, mid = mid_color, high = high_color, limits = c(-0.35, 0.35)) +
    theme_void() +
    labs(title = "ga Volume Mu GA"))

(ga_sa_mu_age_dk_plot <- ggplot(plot_data_dkt %>% filter(domain == "sa")) +
    geom_brain(atlas = dk, 
               position = position_brain(hemi ~ side),
               aes(fill = age_coef)) +
    scale_fill_gradient2(low = low_color, mid = mid_color, high = high_color, limits = c(-0.35, 0.35)) +
    theme_void()+
    labs(title = "ga SA Mu age"))

(ga_th_mu_age_dk_plot <- ggplot(plot_data_dkt %>% filter(domain == "th")) +
    geom_brain(atlas = dk, 
               position = position_brain(hemi ~ side),
               aes(fill = age_coef)) +
    scale_fill_gradient2(low = low_color, mid = mid_color, high = high_color, limits = c(-0.35, 0.35)) +
    theme_void()+
    labs(title = "ga Thickness Mu age"))

if(save_figs == TRUE){
  panel_save(ga_vol_mu_age_dk_plot, "ga-vol-mu-age-dk", "pdf", plot_save_path)
  panel_save(ga_vol_mu_age_aseg_plot, "ga-vol-mu-age-aseg", "pdf", plot_save_path)
  panel_save(ga_sa_mu_age_dk_plot, "ga-sa-mu-age-dk", "pdf", plot_save_path)
  panel_save(ga_th_mu_age_dk_plot, "ga-th-mu-age-dk", "pdf", plot_save_path)
}

# Trajectory Plots Prep

## Fix scale of out_fans, scale everything back by 1000
#Scale all phenotypes by 1000
out_fans_slip <- out_fans_slip %>% 
  mutate(SLIP_input_fans = SLIP_input_fans*1000,
         LMSz_fans = LMSz_fans*1000) %>%
  mutate(sex_ga = paste0(sex, gestational_age))

## Define x axis marks
x_breaks <- log(c(280, 365.25*1+280, 365.25*2+280, 
                  365.25*6+280, 365.25*21+280))
x_labs <- c("0y","1y","2y","6y","21y")

## plotting loop function
separate_plot_loop <- function(pheno_list, df_plot, out_fans_df, x_breaks, x_labs){
  slip_plots <- list()
  lmsz_plots <- list()
  for (pheno_of_int in pheno_list) {
    df_plot$pheno <- df_plot[[glue("{pheno_of_int}.raw")]]
    df_plot$pheno[is.na(df_plot[[glue("{pheno_of_int}.lmsz")]])] <- NA
    if(grepl("global", pheno_of_int)){
      ylab_text = ""
      title_text = sub("global_", "", pheno_of_int)
    } else if(grepl("GrayVol", pheno_of_int))  {
      ylab_text = "Volume"
      title_text = sub("aparc_GrayVol_", "", pheno_of_int)
    } else if(grepl("SurfArea", pheno_of_int))  {
      ylab_text = "Surface Area"
      title_text = sub("aparc_SurfArea_", "", pheno_of_int)
    } else if(grepl("ThickAvg", pheno_of_int))  {
      ylab_text = "Thickness"
      title_text = sub("aparc_ThickAvg_", "", pheno_of_int)
    } else if(grepl("aseg", pheno_of_int))  {
      ylab_text = "Aseg/Subcortical"
      title_text = sub("aseg_", "", pheno_of_int)
    } else{
      ylab_text = pheno_of_int
      title_text = pheno_of_int
    }
    out_fans_plot <- out_fans_df %>% 
      filter(pheno == pheno_of_int,
             !(centile_cat %in% c("cent_0.25", "cent_0.75")))
    
    # CT data is questionable at 0y. Limit to 2 years
    if (grepl("Thick",pheno_of_int)) {
      df_plot <- df_plot %>% filter(adjusted_age_in_days >= 2*365.25+280)
        
      out_fans_plot <- out_fans_plot %>% filter(adjusted_age_days_log >= log(2*365.25+280))
    }
    
    df_plot$sex_ga <- paste0(df_plot$sex, df_plot$gestational_age)
    
    slip_plots[[pheno_of_int]] <- ggplot(out_fans_plot) +
      geom_point(data = df_plot, aes(y = pheno,
                                     x = adjusted_age_days_log, fill=sex), 
                 alpha=0.5, shape=21, size = 0.5) +
      geom_line(aes(x = adjusted_age_days_log, y = SLIP_input_fans,
                    group = interaction(centile_cat, sex),
                    color = sex,
                    linewidth = line_size,
                    linetype = line_type)) +
      scale_linewidth_identity() +
      scale_linetype_identity() +
      scale_colour_manual(values = c("firebrick","royalblue")) +
      scale_x_continuous(breaks = x_breaks, labels = x_labs) +
      ggtitle(title_text) + theme_classic() +
      xlab("") + ylab(ylab_text) +
      theme(plot.title = element_text(size = 12),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 8))
    
    lmsz_plots[[pheno_of_int]] <- ggplot(out_fans_plot %>%
                                           filter(gestational_age %in% c(24, 32, 36, 40, 42))) +
      geom_point(data = df_plot %>% 
                   filter(gestational_age %in% c(24, 32, 36, 40, 42)), aes(y = pheno,
                                     x = adjusted_age_days_log, color = sex_ga), 
                 alpha=0.5, shape=21, size = 0.5) +
      geom_line(aes(x = adjusted_age_days_log, y = LMSz_fans,
                    group = interaction(centile_cat, sex_ga),
                    color = sex_ga,
                    linewidth = line_size,
                    linetype = line_type)) +
      scale_linewidth_identity() +
      scale_linetype_identity() +
      scale_colour_manual(values = c(
        "#E3AEB7", # lightest pink
        "#D07C8F", # lighter pink
        "#B85C74", # medium light rose
        "#A64A63", # medium dark rose
        "#8A2F4D", # deeper red
        "#7FA6C9", # lightest blue
        "#4F85B5", # lighter blue
        "#2F5F8F", # medium light blue
        "#4C78A8",  # medium dark blue
        "#3F51B5" )) +  # indigo blue
      scale_x_continuous(breaks = x_breaks, labels = x_labs) +
      ggtitle(title_text) + theme_classic() +
      xlab("") + ylab(ylab_text) +
      theme(plot.title = element_text(size = 12),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 8))
    
    # Equalize Axes
    p_y_min <- min(layer_scales(slip_plots[[pheno_of_int]])$y$range$range[1],
                   layer_scales(lmsz_plots[[pheno_of_int]])$y$range$range[1])
    p_y_max <- max(layer_scales(slip_plots[[pheno_of_int]])$y$range$range[2],
                   layer_scales(lmsz_plots[[pheno_of_int]])$y$range$range[2])
    
    slip_plots[[pheno_of_int]] <- slip_plots[[pheno_of_int]] + ylim(p_y_min,p_y_max)
    lmsz_plots[[pheno_of_int]] <- lmsz_plots[[pheno_of_int]] + ylim(p_y_min,p_y_max)
  }
  all_plots_list <- list()
  all_plots_list$slip_plots <- slip_plots
  all_plots_list$lmsz_plots <- lmsz_plots
  return(all_plots_list)
}

overlaid_plot_loop <- function(pheno_list, df_plot, out_fans_df, x_breaks, x_labs){
  overlap_plots_male <- list()
  overlap_plots_female <- list()
  overlap_plots_avg <- list()
  
  pheno_of_int <- "global_SurfaceArea"
  
  for (pheno_of_int in pheno_list) {
    df_plot$pheno <- df_plot[[glue("{pheno_of_int}.raw")]]
    df_plot$pheno[is.na(df_plot[[glue("{pheno_of_int}.lmsz")]])] <- NA
    if(grepl("global", pheno_of_int)){
      ylab_text = ""
      title_text = sub("global_", "", pheno_of_int)
    } else if(grepl("GrayVol", pheno_of_int))  {
      ylab_text = "Volume"
      title_text = sub("aparc_GrayVol_", "", pheno_of_int)
    } else if(grepl("SurfArea", pheno_of_int))  {
      ylab_text = "Surface Area"
      title_text = sub("aparc_SurfArea_", "", pheno_of_int)
    } else if(grepl("ThickAvg", pheno_of_int))  {
      ylab_text = "Thickness"
      title_text = sub("aparc_ThickAvg_", "", pheno_of_int)
    } else if(grepl("aseg", pheno_of_int))  {
      ylab_text = "Aseg/Subcortical"
      title_text = sub("aseg_", "", pheno_of_int)
    } else{
      ylab_text = pheno_of_int
      title_text = pheno_of_int
    }
    
    out_fans_plot <- out_fans_df %>% 
      filter(pheno == pheno_of_int,
             !(centile_cat %in% c("cent_0.25", "cent_0.75")))
    
    if (grepl("Thick",pheno_of_int)) {
      df_plot <- df_plot %>% filter(adjusted_age_in_days >= 2*365.25+280)
      out_fans_plot <- out_fans_plot %>% filter(adjusted_age_days_log >= log(2*365.25+280))
    }
    
    
    overlap_plots_male[[pheno_of_int]] <- ggplot(out_fans_plot %>% filter(sex == "Male") %>% 
                                                   filter(gestational_age %in% c(24, 32, 36, 40, 42))) +
      geom_point(data = df_plot %>% filter(sex == "Male"), aes(y = pheno,
                                                            x = adjusted_age_days_log), 
                 alpha=1, size = 0.5, color = "gray") +
      geom_line(aes(x = adjusted_age_days_log, y = SLIP_input_fans,
                    group = centile_cat,
                    color = "SLIP",
                    linewidth = line_size,
                    linetype = line_type)) +
      geom_line(aes(x = adjusted_age_days_log, y = LMSz_fans,
                    group = interaction(centile_cat, gestational_age),
                    color = as.factor(gestational_age),
                    linewidth = line_size,
                    linetype = line_type)) +
      scale_linewidth_identity() +
      scale_linetype_identity() +
      scale_x_continuous(breaks = x_breaks, labels = x_labs) +
      ggtitle(paste0(title_text)) + theme_classic() +
      xlab("") + ylab(ylab_text) +
      scale_colour_manual(values = c(
        "#7FA6C9", # lightest blue
        "#4F85B5", # lighter blue
        "#2F5F8F", # medium light blue
        "#4C78A8",  # medium dark blue
        "#3F51B5",  # indigo blue
        "#2F2F2F")) + # black-ish gray
      theme(plot.title = element_text(size = 12),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 8))
    
    overlap_plots_female[[pheno_of_int]] <- ggplot(out_fans_plot %>% filter(sex == "Female") %>% 
                                                     filter(gestational_age %in% c(24, 32, 36, 40, 42))) +
      geom_point(data = df_plot %>% filter(sex == "Female"), aes(y = pheno,
                                                               x = adjusted_age_days_log), 
                 alpha=1, size = 0.5, color = "gray") +
      geom_line(aes(x = adjusted_age_days_log, y = SLIP_input_fans,
                    group = centile_cat,
                    color = "SLIP",
                    linewidth = line_size,
                    linetype = line_type)) +
      geom_line(aes(x = adjusted_age_days_log, y = LMSz_fans,
                    group = interaction(centile_cat, gestational_age),
                    color = as.factor(gestational_age),
                    linewidth = line_size,
                    linetype = line_type)) +
      scale_linewidth_identity() +
      scale_linetype_identity() +
      scale_x_continuous(breaks = x_breaks, labels = x_labs) +
      ggtitle(paste0(title_text)) + theme_classic() +
      xlab("") + ylab(ylab_text) +
      scale_colour_manual(values = c(
        "#E3AEB7", # lightest pink
        "#D07C8F", # lighter pink
        "#B85C74", # medium light rose
        "#A64A63", # medium dark rose
        "#8A2F4D",  # deeper red
        "#2F2F2F")) + # black-ish gray
      theme(plot.title = element_text(size = 12),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 8))
    
    
    # average over 2 sexes
    out_fans_plot_avg <-   out_fans_plot %>% 
      filter(grepl(str_replace_all(pheno_of_int,"SA.","SA.*"),pheno),
             !(centile_cat %in% c("cent_0.25", "cent_0.75"))) %>%
      mutate(SLIP_input_fans = SLIP_input_fans,
             LMSz_fans = LMSz_fans) %>%
      group_by(model, adjusted_age_days_log, gestational_age,
               centile_cat, line_size, line_type) %>%
      summarize(across(contains("fans"), ~mean(.x, na.rm = T)))
    
    overlap_plots_avg[[pheno_of_int]] <- ggplot(out_fans_plot_avg %>% 
                                                  filter(gestational_age %in% c(24, 32, 36, 40, 42))) +
      geom_point(data = df_plot, aes(y = pheno,
                                                                 x = adjusted_age_days_log), 
                 alpha=1, size = 0.5, color = "gray") +
      geom_line(aes(x = adjusted_age_days_log, y = SLIP_input_fans,
                    group = centile_cat,
                    color = "SLIP",
                    linewidth = line_size,
                    linetype = line_type)) +
      geom_line(aes(x = adjusted_age_days_log, y = LMSz_fans,
                    group = interaction(centile_cat, gestational_age),
                    color = as.factor(gestational_age),
                    linewidth = line_size,
                    linetype = line_type)) +
      scale_linewidth_identity() +
      scale_linetype_identity() +
      scale_x_continuous(breaks = x_breaks, labels = x_labs) +
      ggtitle(paste0(title_text)) + theme_classic() +
      xlab("") + ylab(ylab_text) +
      scale_colour_manual(values = c(
        "#AFC7B5", # lightest green
        "#7FA68F", # lighter green
        "#5F8F73", # medium light green
        "#4C7A5F", # medium dark green
        "#2F5A45",  # deeper green
        "#2F2F2F")) + # black-ish gray
      theme(plot.title = element_text(size = 12),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 8))
    
  }
  all_plots_list <- list()
  all_plots_list$overlap_male <- overlap_plots_male
  all_plots_list$overlap_female<- overlap_plots_female
  all_plots_list$overlap_avg <- overlap_plots_avg
  return(all_plots_list)
}

## Plot global trajectories
global_pheno_to_plot <- c("global_eTIV","global_VentricleChoroidVol","global_CortexVol",
                          "global_SubCortGrayVol","global_SurfaceArea","global_CerebralWhiteMatterVol",
                          "global_MeanThickness")

#Overlaid
global_plots_overlaid <- list()
global_plots_overlaid <- overlaid_plot_loop(global_pheno_to_plot, 
                                            df_ga, 
                                            out_fans_slip, 
                                            x_breaks, 
                                            x_labs)

#overlaid plots by sex
p_male_title <- ggdraw() + draw_label("Males", angle = 90,
                                      x = 0.5, hjust = 0.5, fontface = "bold")
p_male <- ggarrange(plots =  global_plots_overlaid$overlap_male, 
                    nrow = 2, legend = "none", common.legend = T)
p_female_title <- ggdraw() + draw_label("Females", angle = 90,
                                        x = 0.5, hjust = 0.5, fontface = "bold")
p_female <- ggarrange(plots =  global_plots_overlaid$overlap_female, 
                      nrow = 2, legend = "none", common.legend = T)

p_fits_global_overlaid_sex <- (p_male_title + p_male +
                                     p_female_title + p_female) + 
  plot_layout(ncol = 2, 
              nrow = 2, 
              widths = c(0.07, 1)) +
  plot_annotation(title = "Global, by GA")

p_fits_global_overlaid_sex

#overlaid plots averaged over sex
p_traj_sexavg <- ggarrange(plots =  global_plots_overlaid$overlap_avg, 
                           nrow = 2, legend = "bottom", common.legend = T)

p_traj_collapsed_title <- ggdraw() + draw_label("By GA, Sex Avg", angle = 90, x = 0.5, hjust = 0.5, fontface = "bold")


p_fits_global_overlaid_avg <- (p_traj_collapsed_title  + p_traj_sexavg) + 
  plot_layout(ncol = 2, 
              nrow = 1, 
              widths = c(0.07, 1)) +
  plot_annotation(title = "Global")


#Separate
global_plots <- list()
global_plots <- separate_plot_loop(global_pheno_to_plot, 
                                   df_ga, 
                                   out_fans_slip, 
                                   x_breaks, 
                                   x_labs)


p_slip_title <- ggdraw() + draw_label("Original (SLIP)", angle = 90,
                                      x = 0.5, hjust = 0.5, fontface = "bold")

p_slip <- ggarrange(plots =  global_plots$slip_plots, nrow = 2, legend = "none", common.legend = T)

p_lmsz_title <- ggdraw() + draw_label("GA LMSz", angle = 90,
                                      x = 0.5, hjust = 0.5, fontface = "bold")
p_lmsz <- ggarrange(plots =  global_plots$lmsz_plots, nrow = 2, legend = "bottom", common.legend = T)

p_fits_global_separate <- (p_slip_title + p_slip +
                                 p_lmsz_title + p_lmsz) + 
  plot_layout(ncol = 2, 
              nrow = 2, 
              widths = c(0.07, 1)) +
  plot_annotation(title = "Global")

panel_save(p_fits_global_separate, "global_traj_separate", "pdf", plot_save_path, plot_width = 24, plot_height = 16)
panel_save(p_fits_global_overlaid_sex, "global_traj_overlaid_sex", "pdf", plot_save_path, plot_width = 24, plot_height = 16)
panel_save(p_fits_global_overlaid_avg, "global_traj_overlaid_avg", "pdf", plot_save_path, plot_width = 24, plot_height = 16)
  