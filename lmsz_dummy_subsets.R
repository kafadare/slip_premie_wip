# Looping R script to fit random subsets of same N as VPM to LMSz models, and save the lmsz-adjusted centiles output.
# Will use this to compare ICC values between these "fake" subsets and the TRUE VPM subset.
#usage(one argument per line)
# lmsz_dummy_subsets.R 
# iter_no : array index number from array job
# centiles_data_path: "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/median-Th-MPR-distinct-allgrades-25.11/data_clean/median_Th-MPR_2026-02-04.csv"
# raw data path: "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11/raw_csv/median_Th-MPR_2026-02-04.csv" 
# path for original gamlss models: "/mnt/isilon/bgdlab_processing/braincharts/SLIP/2025_03/code/gamlss/RDS/"
# main output folder: "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/lmsz_median_Th-MPR_25.11/"
# n_subset : 181 (number of participants in each dummy subset)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
iter_no <- as.numeric(args[1])
centiles_data_path <- args[2]
raw_data_path <- args[3]
models_path <- args[4]
output_folder <- args[5]
n_subset <- as.numeric(args[6])

#iter_no = 1
#centiles_data_path = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/median-Th-MPR-distinct-allgrades-25.11/data_clean/median_Th-MPR_2026-02-04.csv"
#raw_data_path = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11/raw_csv/median_Th-MPR_2026-02-04.csv" 
#models_path = "/mnt/isilon/bgdlab_processing/braincharts/SLIP/2025_03/code/gamlss/RDS/"
#output_folder = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/lmsz_median_Th-MPR_25.11/"
#n_subset = 181

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

seed_start = 42

# set folder paths
code_path <- "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/"
vars_path <-  "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/vars/"

# output folder
save_path <- output_folder
if (!dir.exists(save_path)) dir.create(save_path)

##  output subfolder paths
csv_save_path <-  paste0(save_path, "results_csv/")
if (!dir.exists(csv_save_path)) dir.create(csv_save_path)

source(glue("{code_path}R_util_lmsz.R"))

# Load Phenotype Variable Names
full_vars <- read.csv(paste0(vars_path,"all_vars.csv")) %>% pull(x)
global_vars <- read.csv(paste0(vars_path,"global_vars.csv")) %>% pull(x)
vol_vars <- read.csv(paste0(vars_path,"vol_vars.csv")) %>% pull(x)
sa_vars <- read.csv(paste0(vars_path,"sa_vars.csv")) %>% pull(x)
th_vars <- read.csv(paste0(vars_path,"th_vars.csv")) %>% pull(x)

## remove aseg_CC_Central and "aseg_CC_Mid_Posterior" from the features because there is one NaN value in the centiles causing issues. In the future will edit the functions to deal with this and/or figure out why there is this NaN value.
full_vars <- full_vars[!(full_vars %in% c("aseg_CC_Central", "aseg_CC_Mid_Posterior"))]

# Load Clean Data
df <- read.csv(centiles_data_path) %>% select(-"X") #%>%
#rename_with(~paste0(.x, ".cent"), all_of(full_vars))

# Load Raw Data and merge
df_raw <- read.csv(raw_data_path) %>% select(all_of(c(full_vars, "participant_id", "session_id", "sex", "age_days",
                                                      "site", "avg_grade", "age_at_scan"))) %>%
  #mutate(across(all_of(full_vars), ~c(scale(.)), .names = "{.col}.raw_normalised")) %>%
  rename_with(~paste0(.x, ".raw"), all_of(full_vars))

shared_names <- intersect(names(df_raw), names(df))

df <- merge(df, df_raw, by = shared_names, all.x = FALSE, all.y = FALSE)

# Remove Age > 21
df <- df %>% filter(adjusted_age_years <= 21)

# Z-score Centiles, remove Inf Cells

## zscore centiles (new df)
#zscore centiles and other variables of interest
df <- df %>%
  mutate(across(all_of(full_vars), qnorm,
                .names = "{.col}.zscore"))
#set sex and preterm as factor
df$sex <- as.factor(df$sex)
df$preterm <- as.factor(df$preterm)
zscore_vars_names <- paste0(full_vars, ".zscore")

## find Inf cells, make a table of N per phenotype, then convert to NaN
inf_nan_table <- as.data.frame(t(sapply(df[zscore_vars_names], function(x){
  c(inf = sum(is.infinite(x)), nan = sum(is.nan(x)))
})))
max(inf_nan_table$inf)
max(inf_nan_table$nan)

inf_nan_table %>% filter(inf > 0)

## set Inf cells to NaN
df <- df %>% 
  mutate(across(all_of(zscore_vars_names), ~ replace(., is.infinite(.), NaN)))

# Define LMSz models to process
# Mean (mu) models
m_models = c('~ adjusted_age_days_log + sex + adjusted_age_days_log*sex',
             '~ adjusted_age_days_log + sex',
             '~ adjusted_age_days_log',
             '~ ns(adjusted_age_days_log, df=2) + sex + ns(adjusted_age_days_log, df=2)*sex',
             '~ ns(adjusted_age_days_log, df=2) + sex',
             '~ ns(adjusted_age_days_log, df=2)',
             '~ ns(adjusted_age_days_log, df=3) + sex + ns(adjusted_age_days_log, df=3)*sex',
             '~ ns(adjusted_age_days_log, df=3) + sex',
             '~ ns(adjusted_age_days_log, df=3)',
             '~ sex',
             '~ 1',
             '~ 0')

# Variance (sigma) models
s_models = c('~ adjusted_age_days_log + sex',
             '~ adjusted_age_days_log',
             '~ sex',
             '~ 0')

# Choose a random subset of a specified number of participants
dummy_save_path <- paste0(save_path,"lmsz_dummy_models/")
seed = seed_start + iter_no
  set.seed(seed)
  df_subset <- df %>% slice_sample(n = n_subset)
  table(df_subset$preterm)
  lmsz <- lmsz_modeling(df_subset, 
                        m_models, 
                        s_models, 
                        full_vars, 
                        models_path, 
                        dummy_save_path, 
                        save_model = FALSE,
                        reweight = FALSE)
  
  models_df_subset <- lmsz$models
  df_out <- lmsz$df
  worm_plots<- lmsz$wp
  rm(lmsz)
  seed = seed + 1
  models_filename = glue("{dummy_save_path}models_set{iter_no}.csv")
  cent_filename = glue("{dummy_save_path}lmsz_adj_centiles_set{iter_no}.csv")
  write_csv(models_df_subset, models_filename)
  write_csv(df_out, cent_filename)

