#Script to calculate centile scores for the new SLIP releases using specified models
#usage(one argument per line)
# get_centiles.R 
# N which pheno (array number): pheno_i
# folder for gamlss fits: "/mnt/isilon/bgdlab_processing/braincharts/SLIP/2025_03/code/gamlss/RDS/" 
# file with raw values: "mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11/CSV/slip_median.csv"
# where to save output centiles: "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11/CSV/centile_csv/slip-median-all/"
# which type of model (should ideally match the file above) "median"
# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
ph_no <- as.numeric(args[1])
fits_folder <- args[2]
raw_data_path <- args[3]
out_path <- args [4]
model_case <- args[5]

#Load packages
require(dplyr)
require(glue)
require(magrittr)

##Testing with default arguments
#ph_no <- 2
#fits_folder <- "/mnt/isilon/bgdlab_processing/braincharts/SLIP/2025_03/code/gamlss/RDS/"
raw_data_path <- "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11/CSV/slip_median.csv"
#out_path <- "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11/CSV/centile_csv/slip-median-all/"
#model_case <- "median"
source("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/940.calc-novel-wo-subset-function.r")

if(!dir.exists(out_path)) dir.create(out_path)

raw_data <- read.csv(raw_data_path)

regions <- names(raw_data)[grep("global|aparc|aseg", names(raw_data))]

#define correct FIT.EXTRACT object
pheno <- regions[ph_no]
fit_name <- glue("{fits_folder}{model_case}-{pheno}Transformed/FIT.EXTRACT.rds")
message("Using fit: ", fit_name)
fit <- readRDS(fit_name)
#select df with covars and appropriate phenotype
df <- raw_data %>% select(all_of(c("participant_id", "session_id", "site", "sex", "age_days", "avg_grade", "dx", pheno)))
#calculate centile scores
out <- Calc.Novel(df, fit)
centile_var_name <- glue("{pheno}Transformed.q.wre")
out_df <- out$data %>% 
  select(all_of(c("participant_id", "session_id", "site", "sex", "age_days", "avg_grade", "dx", centile_var_name))) %>%
  rename(!!pheno := !!centile_var_name)
csv_name <- glue("{out_path}{pheno}_centiles.csv")
write.csv(out_df,csv_name)


