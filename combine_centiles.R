##Script to put together all centile csvs of different phenotypes from one run of get_centiles.R
# N which pheno (array number): pheno_i
# folder for gamlss fits: "/mnt/isilon/bgdlab_processing/braincharts/SLIP/2025_03/code/gamlss/RDS/" 
# file with raw values: "mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11/CSV/slip_median.csv"
# where to save output centiles: "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11/CSV/centile_csv/slip-median-all/"
# which type of model (should ideally match the file above) "median"
# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
out_path <- args[1]
model_type <- args[2]

#Load packages
require(dplyr)
require(glue)
require(magrittr)

#List files
files_list <- list.files(out_path, full.names = TRUE)
all_csv <- lapply(files_list, read.csv)
n_rows <- sapply(all_csv, nrow) %>% unique()
if(length(n_rows) > 1){
  print(glue("Unique number of rows across csvs > 1 (to catch if any CSVs are missing rows: {n_rows}"))
  stop()
}

#check if all participant and session ids are identical across CSVs -- commented out when running to merge for efficiency,
#but might be useful in the future for identifying potential issues
# id_check <- Reduce(
#   function(x, y) {
#     identical(
#       x[, c("participant_id", "session_id")],
#       y[, c("participant_id", "session_id")]
#     )
#   },
#   all_csv
# )

#merge all centile csvs
merged_df <- Reduce(
  function(x, y) merge(x, y, by = intersect(names(x), names(y))),
  all_csv
) %>% select(-c("X"))

#save one merged csv
csv_name <- glue("{out_path}/../slip_{model_type}_centiles.csv")
write.csv(merged_df,csv_name)

#delete all individual csvs and folder
file.remove(files_list)
unlink(out_path, recursive = TRUE)
