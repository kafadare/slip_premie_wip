library(readr)
library(dplyr)

#args <- commandArgs(trailingOnly = TRUE)

#output_folder = args[1]
output_folder = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/lmsz_median_Th-MPR_25.11_GA_splines_v2/"

tmp_models_path <- paste0(output_folder, "results_csv/tmp/models_ga/")
tmp_dfga_path   <- paste0(output_folder, "results_csv/tmp/df_ga/")
csv_save_path   <- paste0(output_folder, "results_csv/")

# # combine all models_ga tmp files, by row (each row is a pheno)
# model_files <- list.files(
#   tmp_models_path,
#   full.names = TRUE,
#   pattern = "\\.csv$"
# )
# 
# models_ga <- bind_rows(lapply(model_files, read_csv, show_col_types = FALSE))
# 
# filename_models <- paste0("lmsz_models_ga_", Sys.Date(),".csv")
# 
# write_csv(models_ga, paste0(csv_save_path, filename_models))

# combine all df_ga tmp files, by column
df_files <- list.files(
  tmp_dfga_path,
  full.names = TRUE,
  pattern = "\\.csv$"
)

df_list <- lapply(df_files, read_csv, show_col_types = FALSE)


df_ga <- df_list[[1]]

for (i in 2:length(df_list)) {
  join_cols <- intersect(names(df_ga), names(df_list[[i]]))
  df_ga <- full_join(df_ga, df_list[[i]], by = join_cols)
}

filename_df <- paste0("lmsz_adjusted_centiles_ga_", Sys.Date(),".csv")
write_csv(df_ga, paste0(csv_save_path, filename_df))