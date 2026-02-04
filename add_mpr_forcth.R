#make median except for CT (mpr) centiles

#Load packages
require(dplyr)
require(glue)
require(magrittr)

args <- commandArgs(trailingOnly = TRUE)
median_file <- as.numeric(args[1])
mpr_file <- args[2]
save_path <- args [3]

median_file = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/median-distinct-allgrades-25.11/data_clean/centiles_distinct_2026-01-29.csv"
mpr_file = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/mpr-distinct-allgrades-25.11/data_clean/centiles_distinct_2026-02-04.csv"
save_path = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/median-ThMPR-distinct-allgrades-25.11/data_clean/"

median <- read.csv(median_file)
mpr <- read.csv(mpr_file)

cth_names <- names(median)[grepl("Thick", names(median))]

merged <- merge(median, 
                mpr %>% select(all_of(c("participant_id", "session_id", cth_names))),
                by = c("participant_id", "session_id"), suffixes = c(".median", ""))


if(!dir.exists(save_path)) dir.create(save_path)
filename <- paste0(save_path, "median_Th-MPR_", Sys.Date(), ".csv")
write.csv(merged, filename)



