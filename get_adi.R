# Script to get median ADI before scan for 25.11 MRI cohort, should be matched to scan date, not participant


#Load packages. Install these packages with install.packages("pkg_name") if they are not installed.
library(bigrquery) #package to import SQL tables into R
library(tidyverse) #package suite for data manipulation, reading data into R, and plotting (ggplot2)
library(glue) #package to use expressions with strings
library(ggseg)
require(gridExtra)
library(interactions)
library(boot)
library(parameters)

panel_save <- function(plot, filename, extension, folder){
  save_name <- paste0(folder, filename, ".", extension)
  ggsave(save_name, plot = plot, dpi = 300)
}

bq_auth()
# proj_id <- if (exists('proj_id')) proj_id else bq_projects()[1]
proj_id = "scit605-healthymri-2d73d69f"

# read in centiles to subset release
cent_path <- "/mnt/arcus/lab/users/kafadare/slip_centiles_processed_tables/centiles_all_median_th_mpr_2026-04-24.csv"
cent <- read.csv(cent_path) %>%
  mutate(proc_ord_id_wrong = proc_ord_id, 
         proc_ord_id = sub(".*ses-([0-9]+)procId.*", "\\1", session_id))

pat_id_list <- paste0("'", sapply(strsplit(cent$pat_id, ","),trimws), "'", collapse = ", ")

env <- readRDS("/mnt/arcus/lab/shared/env_data.rds") %>% filter(pat_id %in% cent$pat_id)

participant_query <- glue('SELECT 
              pat_id, 
              dob, 
              dob_year 
              FROM arcus.patient
              WHERE pat_id IN ({pat_id_list});')
participants <- bq_project_query(proj_id, participant_query) %>% bq_table_download(., page_size=3500)

# Add DOB to env table to calculate age in days for start_date and end_data of env variables
env <- merge(env, participants %>% select(pat_id, dob, dob_year), by = "pat_id", all.x = TRUE)

sum(is.na(env$dob)) #14 entries have no DOB, all same patient id. does not appear in arcus.patient table

env <- env %>% filter(!is.na(dob))

#Calculate unadjusted age in days using dob and end date for a particular entry
env$unadjusted_age_in_days_at_end_date <- as.numeric(difftime(env$end_date, as.Date(env$dob), units = "days"))
env$unadjusted_age_in_days_at_start_date <- as.numeric(difftime(env$start_date, as.Date(env$dob), units = "days"))

env_vars<- c("ADI3_EDU", "ADI3_ECON", "ADI3_FINS", "ADI")

cent_env <- merge(cent %>% mutate(pat_id = sub("sub-", "", participant_id)), env, by = "pat_id", all.x = TRUE)

cent_env$scan_adi_start_age_diff <- cent_env$age_at_scan - cent_env$unadjusted_age_in_days_at_start_date
cent_env$scan_adi_end_age_diff <- cent_env$age_at_scan - cent_env$unadjusted_age_in_days_at_end_date

hist(cent_env$scan_adi_end_age_diff)
hist(cent_env$scan_adi_start_age_diff)

cent_env_median <- cent_env %>%
  group_by(proc_ord_id) %>%
  mutate(env_limit_row = min(scan_adi_start_age_diff[scan_adi_start_age_diff >= 0])) %>%
  filter(scan_adi_start_age_diff >= env_limit_row) %>%
  summarise(
    n_rows = n(),
    across(any_of(names(cent)),
           ~ first(.x),
           .names = "{.col}"),
    across(any_of(env_vars),
           list(
             median = ~ median(.x, na.rm = TRUE),
             n_val = ~ sum(!is.na(.x))
           ),
           .names = "{col}.{fn}")
  )

# Median ADI
hist(cent_env_median$ADI.median)
hist(cent_env_median$ADI.n_val, breaks = 20)

# Get z-scores for ADI measures
cent_env_median.zscore <- cent_env_median %>%
  mutate(across(all_of(c("ADI.median", "ADI3_ECON.median", 
                         "ADI3_FINS.median", "ADI3_EDU.median")), ~ as.numeric(scale(.))))

hist(cent_env_median.zscore$ADI.median)

# Add ADI median z-scores to centile zscore df and re-save
df.zscore_distinct_path <- "/mnt/arcus/lab/users/kafadare/slip_centiles_processed_tables/df_zscore_distinct.csv"
df.zscore.distinct <- read.csv(df.zscore_distinct_path) %>%
  mutate(proc_ord_id_wrong = proc_ord_id, 
         proc_ord_id = sub(".*ses-([0-9]+)procId.*", "\\1", session_id))

df.zscore_all_path <- "/mnt/arcus/lab/users/kafadare/slip_centiles_processed_tables/df_zscore_all.csv"
df.zscore.all <- read.csv(df.zscore_all_path) %>%
  mutate(proc_ord_id_wrong = proc_ord_id, 
         proc_ord_id = sub(".*ses-([0-9]+)procId.*", "\\1", session_id))

df.zscore.distinct <- merge(df.zscore.distinct, cent_env_median.zscore %>% 
                              select(proc_ord_id, ADI.median, ADI3_ECON.median, 
                                     ADI3_FINS.median, ADI3_EDU.median,
                                     ADI.n_val, ADI3_ECON.n_val, 
                                     ADI3_FINS.n_val, ADI3_EDU.n_val),
                            by = "proc_ord_id", all.x = TRUE)
sum(!is.na(df.zscore.distinct$ADI.median))

df.zscore.all <- merge(df.zscore.all, cent_env_median.zscore %>% 
                         select(proc_ord_id, ADI.median, ADI3_ECON.median, 
                                ADI3_FINS.median, ADI3_EDU.median,
                                ADI.n_val, ADI3_ECON.n_val, 
                                ADI3_FINS.n_val, ADI3_EDU.n_val),
                       by = "proc_ord_id", all.x = TRUE)
sum(!is.na(df.zscore.all$ADI.median))


write.csv(df.zscore.distinct, df.zscore_distinct_path)
write.csv(df.zscore.all, df.zscore_all_path)