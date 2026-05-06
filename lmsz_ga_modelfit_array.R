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

# Command Line Args in Order
# --array_id=1 
# --output_folder=... 
# --centiles_data_path=... 
# --raw_data_path=... 
# --models_path=... 
# --save_model=logical 
# --reweight=logical
# paths should be FULL paths to avoid confusion

args <- commandArgs(trailingOnly = TRUE)

array_id <- as.integer(args[1])

print(args)

params <- list(
  output_folder = args[2],
  centiles_data_path = args[3],
  raw_data_path = args[4],
  models_path = args[5],
  save_model <- as.logical(as.numeric(args[6])),
  reweight   <- as.logical(as.numeric(args[7]))
)

print(params)

# set folder paths
code_path <- "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/"
vars_path <-  "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/vars/"
output_folder <- params$output_folder
centiles_data_path <- params$centiles_data_path
raw_data_path <- params$raw_data_path
models_path <- params$models_path

# output folder
save_path <- params$output_folder
if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE)

## output subfolder paths
csv_save_path <-  paste0(save_path, "results_csv/")
if (!dir.exists(csv_save_path)) dir.create(csv_save_path)

plot_save_path <- paste0(save_path, "plots/")
if (!dir.exists(plot_save_path)) dir.create(plot_save_path)

wormplot_save_path <- paste0(plot_save_path, "worm_plots/")
if (!dir.exists(wormplot_save_path)) dir.create(wormplot_save_path)

models_save_path <- paste0(save_path,"lmsz_models/")
if (!dir.exists(models_save_path)) dir.create(models_save_path)

# new output folders for array pieces
tmp_save_path <- paste0(csv_save_path, "tmp/")
if (!dir.exists(tmp_save_path)) dir.create(tmp_save_path)

tmp_models_path <- paste0(tmp_save_path, "models_ga/")
if (!dir.exists(tmp_models_path)) dir.create(tmp_models_path)

tmp_dfga_path <- paste0(tmp_save_path, "df_ga/")
if (!dir.exists(tmp_dfga_path)) dir.create(tmp_dfga_path)

fans_save_path <- paste0(csv_save_path, "out_fans_slip/")
if (!dir.exists(fans_save_path)) dir.create(fans_save_path)

source(glue("{code_path}R_util_lmsz.R"))
source("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/wp_taki_EK.R")

full_vars <- read.csv(paste0(vars_path,"all_vars.csv")) %>% pull(x)
global_vars <- read.csv(paste0(vars_path,"global_vars.csv")) %>% pull(x)
vol_vars <- read.csv(paste0(vars_path,"vol_vars.csv")) %>% pull(x)
sa_vars <- read.csv(paste0(vars_path,"sa_vars.csv")) %>% pull(x)
th_vars <- read.csv(paste0(vars_path,"th_vars.csv")) %>% pull(x)

if (is.na(array_id) || array_id < 1 || array_id > length(full_vars)) {
  stop("array_id / SLURM_ARRAY_TASK_ID is missing or out of bounds.")
}

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
  age = c(NA, "adjusted_age_days_log", "ns(adjusted_age_days_log, 2)", "ns(adjusted_age_days_log, 3)"),
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
                    intercept_mu = NA,
                    ga_term_mu = FALSE,
                    ga_linear_mu = FALSE,
                    ga_spline_2_mu = FALSE,
                    ga_spline_3_mu = FALSE,
                    sex_term_mu = FALSE,
                    ga.sex_mu_interaction = FALSE,
                    age_term_mu = FALSE,
                    age_linear_mu = FALSE,
                    age_spline_2_mu = FALSE,
                    age_spline_3_mu = FALSE,
                    ga.age_mu_interaction = FALSE,
                    age.sex_mu_interaction = FALSE,
                    # sigma coefs
                    intercept_sigma = NA,
                    ga_term_sigma = FALSE,
                    ga_linear_sigma = FALSE,
                    ga_spline_2_sigma = FALSE,
                    ga_spline_3_sigma = FALSE,
                    sex_term_sigma = FALSE,
                    ga.sex_sigma_interaction = FALSE,
                    age_term_sigma = FALSE,
                    age_linear_sigma = FALSE,
                    age_spline_2_sigma = FALSE,
                    age_spline_3_sigma = FALSE,
                    ga.age_sigma_interaction = FALSE,
                    age.sex_sigma_interaction = FALSE)

pheno <- full_vars[array_id]

# Build the growth chart model
p <- "zscore" # specify the phenotype

# Define phenotype
print(pheno)
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
for (m in formulas_df$formula) {
  for (s in formulas_df$formula) {
    print(m)
    print(s)
    m_terms <- unlist(strsplit(sub("~", "", m), "+", fixed = TRUE)) %>% gsub(" ", "", .)
    s_terms <- unlist(strsplit(sub("~", "", s), "+", fixed = TRUE)) %>% gsub(" ", "", .)
    
    if(sum(s_terms %in% m_terms) < length(s_terms)) next
    
    # Run preliminary model - all subjects
    growthChartModel_all <- gamlss(formula = as.formula(glue("zscore {m}")),
                                   sigma.formula = as.formula(glue("zscore {s}")),
                                   family = NO,
                                   data = df_tmp,
                                   control = gamlss.control(n.cyc = 200),
                                   trace = F)
    if(reweight == TRUE){
      # Re-weight observations to exclude individuals with large residuals (> 3.5)
      df_tmp$weight_tmp <- rep(1, dim(df_tmp)[1])
      df_tmp$weight_tmp[which(abs(resid(growthChartModel_all)) > 3.5)] <- 0
      
      print(glue("The Number of Points Weighted 0 due to large residuals > 3.5 for {pheno} is {sum(df_tmp$weight_tmp == 0)}"))
      
      N_reweight <- sum(df_tmp$weight_tmp == 0)
      
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

# Assign TRUE/FALSE based on what terms are in the final model + intercept coef if applicable

## mu
models_ga$intercept_mu[mod] <- fin_model$mu.coefficients["(Intercept)"]
models_ga$ga_term_mu[mod]   <- any(grepl("gestational_age", models_ga$m_model[mod]))
models_ga$ga_linear_mu[mod] <- any(grepl("gestational_age", models_ga$m_model[mod]))  & !any(grepl("ns\\(ge", models_ga$m_model[mod]))
models_ga$ga_spline_2_mu[mod] <- any(grepl("gestational_age", models_ga$m_model[mod]))  & any(grepl("ns\\(gestational_age,2)", models_ga$m_model[mod]))
models_ga$ga_spline_3_mu[mod] <- any(grepl("gestational_age", models_ga$m_model[mod]))  & any(grepl("ns\\(gestational_age,3)", models_ga$m_model[mod]))
models_ga$sex_term_mu[mod]  <- any(grepl("sex", models_ga$m_model[mod]))
ga.sex_pairs <- as.vector(outer(lmsz_vars$ga[!is.na(lmsz_vars$ga)], lmsz_vars$sex[!is.na(lmsz_vars$sex)], paste, sep = ":")) %>%
  gsub(" ", "", .) %>% gsub("\\(", "\\\\(", .)
models_ga$ga.sex_mu_interaction[mod] <- any(grepl(paste(ga.sex_pairs, collapse = "|"),
                                                  models_ga$m_model[mod]))
models_ga$age_term_mu[mod]   <- any(grepl("adjusted_age_days_log", models_ga$m_model[mod]))
models_ga$age_linear_mu[mod] <- any(grepl("adjusted_age_days_log", models_ga$m_model[mod]))  & !any(grepl("ns\\(adj", models_ga$m_model[mod]))
models_ga$age_spline_2_mu[mod] <- any(grepl("adjusted_age_days_log", models_ga$m_model[mod]))  & any(grepl("ns\\(adjusted_age_days_log,2)", models_ga$m_model[mod]))
models_ga$age_spline_3_mu[mod] <- any(grepl("adjusted_age_days_log", models_ga$m_model[mod]))  & any(grepl("ns\\(adjusted_age_days_log,3)", models_ga$m_model[mod]))
ga.age_pairs <- as.vector(outer(lmsz_vars$ga[!is.na(lmsz_vars$ga)], lmsz_vars$age[!is.na(lmsz_vars$age)], paste, sep = ":")) %>%
  gsub(" ", "", .) %>% gsub("\\(", "\\\\(", .)
models_ga$ga.age_mu_interaction[mod] <- any(grepl(paste(ga.age_pairs, collapse = "|"),
                                                  models_ga$m_model[mod]))
age.sex_pairs <- as.vector(outer(lmsz_vars$age[!is.na(lmsz_vars$age)], lmsz_vars$sex[!is.na(lmsz_vars$sex)], paste, sep = ":")) %>%
  gsub(" ", "", .) %>% gsub("\\(", "\\\\(", .)
models_ga$age.sex_mu_interaction[mod] <- any(grepl(paste(age.sex_pairs, collapse = "|"),
                                                   models_ga$m_model[mod]))

## sigma
models_ga$intercept_sigma[mod] <- fin_model$sigma.coefficients["(Intercept)"]
models_ga$ga_term_sigma[mod]   <- any(grepl("gestational_age", models_ga$s_model[mod]))
models_ga$ga_linear_sigma[mod] <- any(grepl("gestational_age", models_ga$s_model[mod]))  & !any(grepl("ns\\(ge", models_ga$s_model[mod])) < 0
models_ga$ga_spline_2_sigma[mod] <- any(grepl("gestational_age", models_ga$s_model[mod]))  & any(grepl("ns\\(gestational_age,2)", models_ga$s_model[mod]))
models_ga$ga_spline_3_sigma[mod] <- any(grepl("gestational_age", models_ga$s_model[mod]))  & any(grepl("ns\\(gestational_age,3)", models_ga$s_model[mod]))
models_ga$sex_term_sigma[mod]  <- any(grepl("sex", models_ga$s_model[mod]))
ga.sex_pairs <- as.vector(outer(lmsz_vars$ga[!is.na(lmsz_vars$ga)], lmsz_vars$sex[!is.na(lmsz_vars$sex)], paste, sep = ":")) %>%
  gsub(" ", "", .) %>% gsub("\\(", "\\\\(", .)
models_ga$ga.sex_sigma_interaction[mod] <- any(grepl(paste(ga.sex_pairs, collapse = "|"),
                                                     models_ga$s_model[mod]))
models_ga$age_term_sigma[mod]   <- any(grepl("adjusted_age_days_log", models_ga$s_model[mod]))
models_ga$age_linear_sigma[mod] <- any(grepl("adjusted_age_days_log", models_ga$s_model[mod]))  & any(grepl("ns\\(adj", models_ga$s_model[mod])) < 0
models_ga$age_spline_2_sigma[mod] <- any(grepl("adjusted_age_days_log", models_ga$s_model[mod]))  & any(grepl("ns\\(adjusted_age_days_log,2)", models_ga$s_model[mod]))
models_ga$age_spline_3_sigma[mod] <- any(grepl("adjusted_age_days_log", models_ga$s_model[mod]))  & any(grepl("ns\\(adjusted_age_days_log,3)", models_ga$s_model[mod]))
ga.age_pairs <- as.vector(outer(lmsz_vars$ga[!is.na(lmsz_vars$ga)], lmsz_vars$age[!is.na(lmsz_vars$age)], paste, sep = ":")) %>%
  gsub(" ", "", .) %>% gsub("\\(", "\\\\(", .)
models_ga$ga.age_sigma_interaction[mod] <- any(grepl(paste(ga.age_pairs, collapse = "|"),
                                                     models_ga$s_model[mod]))
age.sex_pairs <- as.vector(outer(lmsz_vars$age[!is.na(lmsz_vars$age)], lmsz_vars$sex[!is.na(lmsz_vars$sex)], paste, sep = ":")) %>%
  gsub(" ", "", .) %>% gsub("\\(", "\\\\(", .)
models_ga$age.sex_sigma_interaction[mod] <- any(grepl(paste(age.sex_pairs, collapse = "|"),
                                                      models_ga$s_model[mod]))

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

# save worm plot separately
saveRDS(worm_plots_ga[[pheno]], glue("{wormplot_save_path}worm_plot_{pheno}.rds"))

# Update LMSz quantiles
df_tmp <- df_tmp %>%
  add_column(!!glue("{pheno}.lmsz") := pred_og_centiles(fin_model,
                                                        df_tmp)) %>%
  select(-c(pheno, zscore, adjusted_age_in_years, gestational_age))

# Add LMSz centile to final dataframe
df_subset <- df_subset %>% dplyr::left_join(df_tmp,
                                            by = join_by(participant_id, sex,
                                                         adjusted_age_days_log))

df_ga <- df_subset

# save only this phenotype row for models_ga
write_csv(models_ga[mod, ], glue("{tmp_models_path}lmsz_models_ga_{pheno}.csv"))


write_csv(df_ga, glue("{tmp_dfga_path}lmsz_adjusted_centiles_ga_{pheno}.csv"))

# Calculate Centile Fans for this phenotype only

##simulation frames, log days age, sex, and ga
#between 280 days (40 weeks gestation) and 21 years

sim_df_sexga <- expand.grid(
  adjusted_age_days_log = seq(log(280), log(365.25 * 21 + 280), 0.01),
  sex = c("Male", "Female"),
  gestational_age = 22:42
)

sim_list_sexga <- split(sim_df_sexga, paste0(sim_df_sexga$sex, sim_df_sexga$gestational_age, "week"))

tmp <-  backtransform_centile_fans(pheno,
                                   models_save_path,
                                   sim_list_sexga,
                                   df_ga)

write_csv(tmp, glue("{fans_save_path}lmsz_centfans_ga_uncorrected_{pheno}.csv"))