# Script to get median weight/height within 90 days and 180 days of scan for 25.11 MRI cohort, should be matched to scan date, not participant


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

# Download 2025 11 release table
query <- glue("SELECT * FROM `lab.release_2025_11`") #the string here is the SQL query
release <- bq_project_query(proj_id, query) %>% bq_table_download(., page_size=3500)
print(paste0("Arcus lab.release_2511 table, dimensions: ", dim(release)))

# read in centiles to subset release
cent_path <- "/mnt/arcus/lab/users/kafadare/slip_centiles_processed_tables/centiles_all_median_th_mpr_2026-04-24.csv"
cent <- read.csv(cent_path) %>%
  mutate(proc_ord_id_wrong = proc_ord_id, 
         proc_ord_id = sub(".*ses-([0-9]+)procId.*", "\\1", session_id))

release <- release %>% filter(proc_ord_id %in% cent$proc_ord_id)

# Download proc-ord table
query <- glue("SELECT proc_ord_id, proc_ord_datetime, encounter_id FROM arcus.procedure_order AS proc 
              WHERE proc.proc_ord_id IN (SELECT proc_ord_id FROM lab.release_2025_11) ") #the string here is the SQL query
procedure <- bq_project_query(proj_id, query) %>% bq_table_download(., page_size=3500) %>%
  mutate(proc_ord_id = as.numeric(proc_ord_id))
print(paste0("Arcus arcus.procedure_order table with 25.11 proc_ord_ids, dimensions: ", dim(procedure)))

release[!(release$proc_ord_id %in% procedure$proc_ord_id), "study_id"]

# Download proc-ord table from 2023-02
query <- glue("SELECT proc_ord_id, proc_ord_datetime, encounter_id FROM arcus_2023_02_01.procedure_order AS proc 
              WHERE proc.proc_ord_id IN (SELECT proc_ord_id FROM lab.release_2025_11)") #the string here is the SQL query
procedure_2302 <- bq_project_query(proj_id, query) %>% bq_table_download(., page_size=3500)
print(paste0("Arcus arcus.procedure_ordern from 2023-02-01 table with 25.11 proc_ord_ids, dimensions: ", dim(procedure_2302)))

sum(!(release$proc_ord_id %in% procedure$proc_ord_id))/nrow(release)
release[!(release$proc_ord_id %in% procedure$proc_ord_id), "study_id"]

test <- merge(cent, procedure, by = "proc_ord_id", all.x = TRUE, all.y = FALSE)
sum(is.na(release$proc_ord_datetime))




cent <- merge(cent, release %>% select(c("proc_ord_id", "proc_ord_datetime", "encounter_id", "pat_id_age")), 
              by = "pat_id_age", all.x = TRUE, all.y = FALSE)

# Download Encounter table
query <- glue("SELECT pat_id, encounter_id, effective_datetime, effective_age,
weight_kg, height_cm,
coverage_id, 
hosp_admsn_age, hosp_admsn_datetime, 
hosp_disch_age, hosp_disch_datetime,  
inp_adm_age, inp_adm_datetime, length_of_stay FROM arcus.encounter AS encounter WHERE encounter.pat_id IN (SELECT pat_id FROM lab.release_2025_11) ") #the string here is the SQL query
encounter <- bq_project_query(proj_id, query) %>% bq_table_download(., page_size=3500) #assign the name of the dataframe in R here
print(paste0("Arcus arcus.encounter table with 25.11 participants, dimensions: ", dim(encounter))) #pastes the dimensions of the table downloaded

# LATER: get encounter dx through the lab.phecodes_20260228 table


# Merge encounter table with release table, keeping all encounters and selecting weight/height as well as admission data
encounter <- merge(encounter, release %>% select(-c("encounter_id")), by = c("pat_id"), all = TRUE) %>%
  mutate(mri_enc_timediff_days = difftime(effective_datetime, proc_ord_datetime, units = "days"))

range(encounter$mri_enc_timediff_days)

# Make a table of weight/height measures per participant for 90 and 180 days
wh_median_table <- encounter %>%
  group_by(proc_ord_id) %>%
 # mutate(env_limit_row = min(swyc_start_date_age_diff[swyc_start_date_age_diff >= 0])) %>%
  filter(mri_enc_timediff_days <= 180) %>%
  summarise(
    n_rows_wt = sum(!is.na(weight_kg)),
    n_rows_ht = sum(!is.na(height_cm)),
    pat_id = first(pat_id)
    proc_ord_datetime = first(proc_ord_datetime),
    age_at_scan = first(age_at_scan),
    coverage_id = first(coverage_id),
    swyc_start_date_age_diff   = first(swyc_start_date_age_diff),
    swyc_end_date_age_diff   = first(swyc_end_date_age_diff),
    dob = first(dob),
    dob_year = first(dob_year),
    pat_form = first(pat_form),
    across(any_of(env_vars_of_int),
           list(
             median = ~ median(.x, na.rm = TRUE),
             n_val = ~ sum(!is.na(.x))
           ),
           .names = "{col}.{fn}"),
    .groups = "drop"
  )


# Take a median of 90 days +/- scan and 180 days +/- scan




# For hospital/inpatient stuff, look at individuals with inp_admission_datetime in first months of life (check if there is a hosp_discharge datetime too)


# Also need to get diagnosis from procedure


# Coverage ID might be for the medicaid stuff (might need to use based on the scan)






