# Script to get median weight/height within 180 days of scan for 25.11 MRI cohort, should be matched to scan date, not participant


#Load packages. Install these packages with install.packages("pkg_name") if they are not installed.
library(bigrquery) #package to import SQL tables into R
library(tidyverse) #package suite for data manipulation, reading data into R, and plotting (ggplot2)
library(glue) #package to use expressions with strings
library(ggseg)
require(gridExtra)
library(interactions)
library(boot)
library(parameters)
library(childsds)

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

# Download flowsheets table for weight/height for patients in the centiles
flowsheet_query <- glue('SELECT
  pat_id,
  recorded_age,
  recorded_datetime,
  recorded_year,
  flowsheet_name,
  flowsheet_disp_name,
  flowsheet_value,
  flowsheet_unit,
  inpatient_data_id
  FROM arcus.flowsheet_measure
  WHERE flowsheet_name IN ("HEIGHT", "WEIGHT/SCALE")
                AND pat_id IN ({pat_id_list});')
flowsheet <- bq_project_query(proj_id, flowsheet_query) %>% bq_table_download(., page_size=3500)
print(paste0("arcus.flowsheet_measure table with 25.11 pat_ids, dimensions: ", dim(flowsheet)))

encounter_query <- glue('SELECT
  pat_id,
  encounter_id,
  primary_care_ind,
  specialty_care_ind,
  inpatient_ind,
  ed_ind,
  urgent_care_ind,
  enc_type_name,
  effective_datetime,
  effective_age,
  effective_year,
  length_of_stay,
  inpatient_data_id,
  height_cm,
  weight_kg,
  height_raw,
  weight_oz,
  hosp_admsn_age, 
  hosp_admsn_datetime, 
  hosp_disch_age, 
  hosp_disch_datetime,  
  inp_adm_age, 
  inp_adm_datetime,
  coverage_id
  FROM arcus.encounter
  WHERE pat_id IN ({pat_id_list});')
encounter <- bq_project_query(proj_id, encounter_query) %>% bq_table_download(., page_size=3500)
print(paste0("arcus.encounter table with 25.11 pat_ids, dimensions: ", dim(encounter)))

# Check the nature of weight/height in the flowsheet table
head(flowsheet[flowsheet$flowsheet_name == "HEIGHT",]$flowsheet_value) #appears to be in inches
head(flowsheet[flowsheet$flowsheet_name == "WEIGHT/SCALE",]$flowsheet_value) #appears to be in lbs

index_vars <- c("inpatient_ind", "primary_care_ind", "ed_ind", "specialty_care_ind", "urgent_care_ind")

# WEIGHT
## Get median weights for each subject, compare using both, only encounter table, only flow table weights
weight_180 <-  merge(encounter %>% 
                       select(c("pat_id", "weight_oz", "inpatient_data_id", "effective_age",
                                "inpatient_ind", "primary_care_ind", "ed_ind", "specialty_care_ind", "urgent_care_ind", "length_of_stay")), 
                     flowsheet %>% filter(flowsheet_name == "WEIGHT/SCALE") %>%
                       select(c("pat_id", "flowsheet_value", "inpatient_data_id")) %>%
                       mutate(flowsheet_value = as.numeric(flowsheet_value)),
                     by = c("pat_id", "inpatient_data_id"), all = TRUE) %>%
  pivot_longer(cols = c("weight_oz", "flowsheet_value"), 
               names_to = "weight_from", 
               values_to = "weight_oz",
               values_drop_na = TRUE) %>%
  mutate(weight_from = sub("weight_oz", "encounter_value", weight_from)) %>%
  filter(effective_age >= 0) %>%
  inner_join(cent %>% select(pat_id, age_at_scan),
             by = "pat_id") %>% 
  filter(abs(effective_age - age_at_scan) <= 180) %>%
  group_by(pat_id, age_at_scan) %>%
  summarise(across(all_of(index_vars),
                   list(n_avail = ~sum(.x == "1")),
                        .names = "{.col}_n"),
                   median_weight_oz = median(weight_oz),
                   .groups = "drop") %>%
  mutate(median_weight_kg = median_weight_oz*0.0283495, 
         age_at_scan_years = age_at_scan/365.25)

weight_180_onlyflow <-  merge(encounter %>% 
                                select(c("pat_id",  "inpatient_data_id", "effective_age")), 
                              flowsheet %>% filter(flowsheet_name == "WEIGHT/SCALE") %>%
                                select(c("pat_id", "flowsheet_value", "inpatient_data_id")) %>%
                                mutate(flowsheet_value = as.numeric(flowsheet_value)),
                              by = c("pat_id", "inpatient_data_id"), all = TRUE) %>%
  filter(effective_age >= 0) %>%
  inner_join(cent %>% select(pat_id, age_at_scan),
             by = "pat_id") %>% 
  filter(abs(effective_age - age_at_scan) <= 180) %>%
  group_by(pat_id, age_at_scan) %>%
  summarise(median_weight_oz = median(flowsheet_value, na.rm = TRUE),
            .groups = "drop") %>%
  filter(!is.na(median_weight_oz)) %>%
  mutate(median_weight_kg = median_weight_oz*0.0283495, 
         age_at_scan_years = age_at_scan/365.25)

weight_180_onlyenc <-  encounter %>%
  filter(effective_age >= 0) %>%
  inner_join(cent %>% select(pat_id, age_at_scan),
             by = "pat_id") %>% 
  filter(abs(effective_age - age_at_scan) <= 180) %>%
  group_by(pat_id, age_at_scan) %>%
  summarise(across(all_of(index_vars),
                   list(n_avail = ~sum(.x == "1")),
                   .names = "{.col}_n"),
            median_weight_oz = median(weight_oz, na.rm = TRUE),
            .groups = "drop") %>%
  filter(!is.na(median_weight_oz)) %>%
  mutate(median_weight_kg = median_weight_oz*0.0283495, 
         age_at_scan_years = age_at_scan/365.25)

weight_diff <- merge(weight_180_onlyenc %>% select(pat_id, age_at_scan, median_weight_oz) %>% rename(median_weight_oz_enc = median_weight_oz),
                     weight_180_onlyflow %>% select(pat_id, age_at_scan, median_weight_oz) %>% rename(median_weight_oz_flow = median_weight_oz),
                     by = c("pat_id", "age_at_scan"), all.x = FALSE, all.y = FALSE) %>%
  mutate(diff_enc_flow = median_weight_oz_enc - median_weight_oz_flow) %>%
  filter(diff_enc_flow != 0) %>%
  pivot_longer(cols = c("median_weight_oz_enc", "median_weight_oz_flow"),
               names_to = "median_from",
               values_to = "median_weight") %>%
  mutate(median_weight_kg = median_weight*0.0283495,
         diff_enc_flow_kg = diff_enc_flow*0.0283495,
         age_at_scan_years = age_at_scan/365.25)

hist(weight_diff$diff_enc_flow_kg, breaks = 50)

ggplot(weight_diff %>% filter(diff_enc_flow_kg > 0.5), 
       aes(x = age_at_scan_years, y = median_weight_kg, color = diff_enc_flow_kg)) +
  geom_point() +
  labs(x = "Age At Scan (years)", 
       y = "Median Weight (kg)",
       title = "Enc and Flow Median Diff > 0.5kg")


ggplot(weight_diff, aes(age_at_scan_years, y = diff_enc_flow_kg)) +
  geom_point() +
  labs(x = "Age At Scan (years)", 
       y = "Enc Median - Flow Median (kg)")

ggplot(weight_180, 
       aes(x = age_at_scan_years, y = median_weight_kg)) +
  geom_point() +
  labs(x = "Age At Scan (years)", 
       y = "Median Weight (kg)",
       title = glue("Using both Enc and Flow, N = {nrow(weight_180)}"))

ggplot(weight_180_onlyflow, 
       aes(x = age_at_scan_years, y = median_weight_kg)) +
  geom_point() +
  labs(x = "Age At Scan (years)", 
       y = "Median Weight (kg)",
       title = glue("Using only Flow, N = {nrow(weight_180_onlyflow)}"))

ggplot(weight_180_onlyenc, 
       aes(x = age_at_scan_years, y = median_weight_kg)) +
  geom_point() +
  labs(x = "Age At Scan (years)", 
       y = "Median Weight (kg)",
       title = glue("Using only Enc, N = {nrow(weight_180_onlyenc)}"))



# HEIGHT
## Get median heights for each subject, compare using both, only encounter table, only flow table heights
height_180 <-  merge(encounter %>% 
                       select(c("pat_id", "height_cm", "inpatient_data_id", "effective_age",
                                "inpatient_ind", "primary_care_ind", "ed_ind", "specialty_care_ind", "urgent_care_ind", "length_of_stay")), 
                     flowsheet %>% filter(flowsheet_name == "HEIGHT") %>%
                       select(c("pat_id", "flowsheet_value", "inpatient_data_id")) %>%
                       mutate(flowsheet_value = as.numeric(flowsheet_value)*2.54),
                     by = c("pat_id", "inpatient_data_id"), all = TRUE) %>%
  pivot_longer(cols = c("height_cm", "flowsheet_value"), 
               names_to = "height_from", 
               values_to = "height_cm",
               values_drop_na = TRUE) %>%
  mutate(height_from = sub("height_cm", "encounter_value", height_from)) %>%
  filter(effective_age >= 0) %>%
  inner_join(cent %>% select(pat_id, age_at_scan),
             by = "pat_id") %>% 
  filter(abs(effective_age - age_at_scan) <= 180) %>%
  group_by(pat_id, age_at_scan) %>%
  summarise(across(all_of(index_vars),
                   list(n_avail = ~sum(.x == "1")),
                        .names = "{.col}_n"),
                   median_height_cm = median(height_cm),
                   .groups = "drop") %>%
  mutate(age_at_scan_years = age_at_scan/365.25)

height_180_onlyflow <-  merge(encounter %>% 
                                select(c("pat_id",  "inpatient_data_id", "effective_age")), 
                              flowsheet %>% filter(flowsheet_name == "HEIGHT") %>%
                                select(c("pat_id", "flowsheet_value", "inpatient_data_id")) %>%
                                mutate(flowsheet_value = as.numeric(flowsheet_value)*2.54),
                              by = c("pat_id", "inpatient_data_id"), all = TRUE) %>%
  filter(effective_age >= 0) %>%
  inner_join(cent %>% select(pat_id, age_at_scan),
             by = "pat_id") %>% 
  filter(abs(effective_age - age_at_scan) <= 180) %>%
  group_by(pat_id, age_at_scan) %>%
  summarise(median_height_cm = median(flowsheet_value, na.rm = TRUE),
            .groups = "drop") %>%
  filter(!is.na(median_height_cm)) %>%
  mutate(age_at_scan_years = age_at_scan/365.25)

height_180_onlyenc <-  encounter %>%
  filter(effective_age >= 0) %>%
  inner_join(cent %>% select(pat_id, age_at_scan),
             by = "pat_id") %>% 
  filter(abs(effective_age - age_at_scan) <= 180) %>%
  group_by(pat_id, age_at_scan) %>%
  summarise(across(all_of(index_vars),
                   list(n_avail = ~sum(.x == "1")),
                   .names = "{.col}_n"),
            median_height_cm = median(height_cm, na.rm = TRUE),
            .groups = "drop") %>%
  filter(!is.na(median_height_cm)) %>%
  mutate(age_at_scan_years = age_at_scan/365.25)

height_diff <- merge(height_180_onlyenc %>% select(pat_id, age_at_scan, median_height_cm) %>% rename(median_height_cm_enc = median_height_cm),
                     height_180_onlyflow %>% select(pat_id, age_at_scan, median_height_cm) %>% rename(median_height_cm_flow = median_height_cm),
                     by = c("pat_id", "age_at_scan"), all.x = FALSE, all.y = FALSE) %>%
  mutate(diff_enc_flow_cm = median_height_cm_enc - median_height_cm_flow) %>%
  filter(diff_enc_flow_cm != 0) %>%
  pivot_longer(cols = c("median_height_cm_enc", "median_height_cm_flow"),
               names_to = "median_from",
               values_to = "median_height_cm") %>%
  mutate(age_at_scan_years = age_at_scan/365.25)

hist(height_diff$diff_enc_flow_cm, breaks = 50)

ggplot(height_diff %>% filter(diff_enc_flow_cm > 1), 
       aes(x = age_at_scan_years, y = median_height_cm, color = diff_enc_flow_cm)) +
  geom_point() +
  labs(x = "Age At Scan (years)", 
       y = "Median Height (cm)",
       title = "Enc and Flow Median Diff > 1cm")


ggplot(height_diff, aes(age_at_scan_years, y = diff_enc_flow_cm)) +
  geom_point() +
  labs(x = "Age At Scan (years)", 
       y = "Enc Median - Flow Median (cm)")

ggplot(height_180, 
       aes(x = age_at_scan_years, y = median_height_cm)) +
  geom_point() +
  labs(x = "Age At Scan (years)", 
       y = "Median Height (cm)",
       title = glue("Using both Enc and Flow, N = {nrow(height_180)}"))

ggplot(height_180_onlyflow, 
       aes(x = age_at_scan_years, y = median_height_cm)) +
  geom_point() +
  labs(x = "Age At Scan (years)", 
       y = "Median Height (cm)",
       title = glue("Using only Flow, N = {nrow(height_180_onlyflow)}"))

ggplot(height_180_onlyenc, 
       aes(x = age_at_scan_years, y = median_height_cm)) +
  geom_point() +
  labs(x = "Age At Scan (years)", 
       y = "Median Height (cm)",
       title = glue("Using only Enc, N = {nrow(height_180_onlyenc)}"))


# Calculate percentile for height/weight using packages and save tables

## weight

weight_180 <- merge(weight_180, cent %>% 
                      distinct(pat_id, .keep_all = TRUE) %>% 
                      select(c("pat_id", "sex")), by = "pat_id", all.x = TRUE, all.y = FALSE) %>%
  mutate(weight_zscore = case_when(
    age_at_scan_years <= 2 ~ sds(value = median_weight_kg,
                                 age = age_at_scan_years,
                                 sex = sex,
                                 item = "weight", 
                                 ref = who.ref,
                                 type = "SDS",
                                 male = "Male", 
                                 female = "Female"),
   TRUE ~ sds(value = median_weight_kg,
             age = age_at_scan_years,
             sex = sex,
             item = "weight2_20", 
             ref = cdc.ref,
             type = "SDS",
             male = "Male", 
             female = "Female")
   ) ) %>% 
    mutate(weight_cent = pnorm(weight_zscore))

sum(is.na(weight_180$weight_zscore))
weight_180[is.na(weight_180$weight_zscore), "age_at_scan_years"]


ggplot(weight_180, 
       aes(x = age_at_scan_years, y = weight_zscore)) +
  geom_point() +
  labs(x = "Age At Scan (years)", 
       y = "Weight Zscore (cm)",
       title = glue("Avail Up to 20 yrs, N = {sum(!is.na(weight_180$weight_zscore))}")) +
  theme_classic()


height_180 <- merge(height_180, cent %>% 
                      distinct(pat_id, .keep_all = TRUE) %>% 
                      select(c("pat_id", "sex")), by = "pat_id", all.x = TRUE, all.y = FALSE) %>%
  mutate(height_zscore = case_when(
    age_at_scan_years <= 2 ~ sds(value = median_height_cm,
                                 age = age_at_scan_years,
                                 sex = sex,
                                 item = "height", 
                                 ref = who.ref,
                                 type = "SDS",
                                 male = "Male", 
                                 female = "Female"),
    TRUE ~ sds(value = median_height_cm,
               age = age_at_scan_years,
               sex = sex,
               item = "height2_20", 
               ref = cdc.ref,
               type = "SDS",
               male = "Male", 
               female = "Female")
  ) ) %>% 
  mutate(height_cent = pnorm(height_zscore))

sum(is.na(height_180$height_zscore))
height_180[is.na(height_180$height_zscore), "age_at_scan_years"]


ggplot(height_180, 
       aes(x = age_at_scan_years, y = height_zscore)) +
  geom_point() +
  labs(x = "Age At Scan (years)", 
       y = "height Zscore (cm)",
       title = glue("Avail Up to 20 yrs, N = {sum(!is.na(height_180$height_zscore))}")) +
  theme_classic()

wh <- merge(weight_180 %>% 
              select(pat_id, age_at_scan, age_at_scan_years, sex, median_weight_kg, weight_zscore, weight_cent),
            height_180 %>% 
              select(pat_id, age_at_scan, age_at_scan_years, sex, median_height_cm, height_zscore, height_cent), 
            by = c("pat_id", "age_at_scan", "age_at_scan_years", "sex"), all = TRUE)

write.csv(weight_180, "/mnt/arcus/lab/users/kafadare/slip_download_tables/weight_zscores.csv")
write.csv(height_180, "/mnt/arcus/lab/users/kafadare/slip_download_tables/height_zscores.csv")
write.csv(wh, "/mnt/arcus/lab/users/kafadare/slip_download_tables/wh_zscores.csv")
