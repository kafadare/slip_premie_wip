# make median except for CT (mpr): 
# provide paths for both centiles and raw segmentation outputs

#Load packages
require(dplyr)
require(glue)
require(magrittr)

`%nin%` = Negate(`%in%`)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
  median_cent_distinct = args[1]
  mpr_cent_distinct = args[2]
  median_cent_all = args[3]
  mpr_cent_all = args[4]
  median_raw_file = args[5]
  mpr_raw_file = args[6]
  cent_save_path = args[7]
  raw_save_path = args[8]
} else {
median_cent_distinct = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11_new_models/data_clean/csv/centiles_distinct_median_2026-04-24.csv"
mpr_cent_distinct = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11_new_models/data_clean/csv/centiles_distinct_mpr_2026-04-24.csv"
median_cent_all = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11_new_models/data_clean/csv/centiles_all_median_2026-04-24.csv"
mpr_cent_all = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11_new_models/data_clean/csv/centiles_all_mpr_2026-04-24.csv"
median_raw_file = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11_new_models/raw_csv/slip_median_041526.csv"
mpr_raw_file = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11_new_models/raw_csv/slip_mpr_041526.csv"
cent_save_path = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11_new_models/data_clean/csv/"
raw_save_path = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11_new_models/raw_csv/"
}

if(!dir.exists(cent_save_path)) dir.create(cent_save_path)
if(!dir.exists(raw_save_path)) dir.create(raw_save_path)

median_distinct <- read.csv(median_cent_distinct) %>% select(-("X"))
mpr_distinct <- read.csv(mpr_cent_distinct) %>% select(-("X"))
median_all <- read.csv(median_cent_all) %>% select(-("X"))
mpr_all <- read.csv(mpr_cent_all) %>% select(-("X"))
median_raw <- read.csv(median_raw_file)
mpr_raw <- read.csv(mpr_raw_file)

cth_names <- names(median_distinct)[grepl("Thick", names(median_distinct))]

merged_distinct <- merge(median_distinct, 
                         mpr_distinct %>% select(all_of(c("participant_id", "session_id", "scan_id", "euler_mean", cth_names))) %>% 
                           rename(euler_mean_mpr = euler_mean, scan_id_mpr = scan_id),
                         by = c("participant_id", "session_id"), suffixes = c(".median", ""), all.x = TRUE)


filename <- paste0(cent_save_path, "centiles_distinct_median_th_mpr_", Sys.Date(), ".csv")
write.csv(merged_distinct, filename)


cth_names <- names(median_all)[grepl("Thick", names(median_all))]

merged_all <- merge(median_all, 
                         mpr_all %>% select(all_of(c("participant_id", "session_id", "scan_id", "euler_mean", cth_names))) %>% 
                           rename(euler_mean_mpr = euler_mean, scan_id_mpr = scan_id),
                         by = c("participant_id", "session_id"), suffixes = c(".median", ""), all.x = TRUE)


filename <- paste0(cent_save_path, "centiles_all_median_th_mpr_", Sys.Date(), ".csv")
write.csv(merged_all, filename)


cth_names <- names(median_raw)[grepl("Thick", names(median_raw))]

merged_raw <- merge(median_raw, 
                         mpr_raw %>% select(all_of(c("participant_id", "session_id", "scan_id", cth_names))) %>% 
                           rename(scan_id_mpr = scan_id),
                         by = c("participant_id", "session_id"), suffixes = c(".median", ""), all.x = TRUE)


filename <- paste0(raw_save_path, "slip_median_th_mpr_", Sys.Date(), ".csv")
write.csv(merged_raw, filename)
