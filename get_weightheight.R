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

# Check the nature of weight/height in the flowsheet table
head(flowsheet[flowsheet$flowsheet_name == "HEIGHT",]$flowsheet_value) #appears to be in inches
head(flowsheet[flowsheet$flowsheet_name == "WEIGHT/SCALE",]$flowsheet_value) #appears to be in lbs

index_vars <- c("inpatient_ind", "primary_care_ind", "ed_ind", "specialty_care_ind", "urgent_care_ind")

# WEIGHT
## Get median weights for each subject using flowsheet
weight_180_onlyflow <-  flowsheet %>% 
  filter(flowsheet_name == "WEIGHT/SCALE") %>%
  select(c("pat_id", "flowsheet_value", "inpatient_data_id", "recorded_age")) %>%
  mutate(flowsheet_value = as.numeric(flowsheet_value)) %>%
  rename(age_at_anthro = recorded_age) %>%
  filter(age_at_anthro >= 0) %>%
  inner_join(cent %>% select(pat_id, age_at_scan),
             by = "pat_id") %>% 
  filter(abs(age_at_anthro - age_at_scan) <= 180) %>%
  group_by(pat_id, age_at_scan) %>%
  summarise(median_weight_oz = median(flowsheet_value, na.rm = TRUE),
            .groups = "drop") %>%
  filter(!is.na(median_weight_oz)) %>%
  mutate(median_weight_kg = median_weight_oz*0.0283495, 
         age_at_scan_years = age_at_scan/365.25)

ggplot(weight_180, 
       aes(x = age_at_scan_years, y = median_weight_kg)) +
  geom_point() +
  labs(x = "Age At Scan (years)", 
       y = "Median Weight (kg)",
       title = glue("Using both Enc and Flow, N = {nrow(weight_180)}"))


# HEIGHT
## Get median heights for each subject using flowsheet
height_180 <- flowsheet %>% 
  filter(flowsheet_name == "HEIGHT") %>%
  select(c("pat_id", "flowsheet_value", "inpatient_data_id", "recorded_age")) %>%
  mutate(flowsheet_value = as.numeric(flowsheet_value)) %>%
  rename(age_at_anthro = recorded_age) %>%
  filter(age_at_anthro >= 0) %>%
  inner_join(cent %>% select(pat_id, age_at_scan),
             by = "pat_id") %>% 
  filter(abs(age_at_anthro - age_at_scan) <= 180) %>%
  group_by(pat_id, age_at_scan) %>%
  summarise(median_height_in = median(flowsheet_value, na.rm = TRUE),
            .groups = "drop") %>%
  filter(!is.na(median_height_in)) %>%
  mutate(median_height_cm = median_height_in*2.54,
    age_at_scan_years = age_at_scan/365.25)

ggplot(height_180, 
       aes(x = age_at_scan_years, y = median_height_cm)) +
  geom_point() +
  labs(x = "Age At Scan (years)", 
       y = "Median Height (cm)",
       title = glue("Using both Enc and Flow, N = {nrow(height_180)}"))


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

# Save before removing for QC
write.csv(weight_180, "/mnt/arcus/lab/users/kafadare/slip_download_tables/weight_zscores.csv")
write.csv(height_180, "/mnt/arcus/lab/users/kafadare/slip_download_tables/height_zscores.csv")

# Remove weight and height zscores > 5SDs, and infinites
## Document how many removed
weight_180_qc <- weight_180 %>%
  filter(abs(weight_zscore) < 5)
height_180_qc <- height_180 %>%
  filter(abs(height_zscore) < 5)

wh <- merge(weight_180_qc %>% 
              select(pat_id, age_at_scan, age_at_scan_years, sex, median_weight_kg, weight_zscore, weight_cent),
            height_180_qc %>% 
              select(pat_id, age_at_scan, age_at_scan_years, sex, median_height_cm, height_zscore, height_cent), 
            by = c("pat_id", "age_at_scan", "age_at_scan_years", "sex"), all = TRUE)

wh <- merge(wh, cent %>% select(proc_ord_id, pat_id, age_at_scan),
            by = c("pat_id", "age_at_scan"), all.x = TRUE, all.y = FALSE)

write.csv(wh, "/mnt/arcus/lab/users/kafadare/slip_download_tables/wh_zscores_QC.csv")


# Add weight/height zscores to centile zscore df and save full df
df.zscore_distinct_path <- "/mnt/arcus/lab/users/kafadare/slip_centiles_processed_tables/df_zscore_distinct.csv"
df.zscore.distinct <- read.csv(df.zscore_distinct_path) %>%
  mutate(proc_ord_id_wrong = proc_ord_id, 
         proc_ord_id = sub(".*ses-([0-9]+)procId.*", "\\1", session_id))

df.zscore_all_path <- "/mnt/arcus/lab/users/kafadare/slip_centiles_processed_tables/df_zscore_all.csv"
df.zscore.all <- read.csv(df.zscore_all_path) %>%
  mutate(proc_ord_id_wrong = proc_ord_id, 
         proc_ord_id = sub(".*ses-([0-9]+)procId.*", "\\1", session_id))

df.zscore.distinct <- merge(df.zscore.distinct, wh %>% 
                                   select(proc_ord_id, weight_zscore, height_zscore),
                                 by = "proc_ord_id", all.x = TRUE)
sum(!is.na(df.zscore.distinct$weight_zscore))
sum(!is.na(df.zscore.distinct$height_zscore))

df.zscore.all <- merge(df.zscore.all, wh %>% 
                                   select(proc_ord_id, weight_zscore, height_zscore),
                                 by = "proc_ord_id", all.x = TRUE)
sum(!is.na(df.zscore.all$weight_zscore))
sum(!is.na(df.zscore.all$height_zscore))

write.csv(df.zscore.distinct, df.zscore_distinct_path)
write.csv(df.zscore.all, df.zscore_all_path)
