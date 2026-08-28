## Merges the per-variant outputs written by run_spline_analysis_variant.R
## (one call per covariate-set variant, run in parallel) into the combined
## long-format tables used for visualization/comparison across variants.
## No p-value re-adjustment happens here -- each per-variant CSV's *_p_adj
## column was already computed within its own (predictor x phenotype_group[ x
## moderator]) batch when it was written, which is the correction scope that
## was requested; merging here is a plain row-bind.
##
## Command line args, in order:
##   1. output_folder -- same folder passed to run_spline_analysis_variant.R
##                        (reads from <output_folder>/models/ and /results_csv/,
##                        writes the merged files back into the same folders)

require(dplyr)

#args <- commandArgs(trailingOnly = TRUE)
#output_folder <- args[1]
output_folder <- "/mnt/arcus/lab/users/kafadare/slip_premie_results_fxfalse/"

# models_folder <- paste0(output_folder, "models/")
stats_folder <- paste0(output_folder, "results_csv/")

read_variant_csvs <- function(suffix) {
  files <- list.files(stats_folder, pattern = paste0("_", suffix, "\\.csv$"), full.names = TRUE)
  if (length(files) == 0) stop("No files found matching *_", suffix, ".csv in ", stats_folder)
  bind_rows(lapply(files, read.csv))
}

main_stats_long <- read_variant_csvs("main_effects_stats")
resi_stats_long <- read_variant_csvs("resi_stats")
int_stats_long  <- read_variant_csvs("interaction_stats")

write.csv(main_stats_long, file.path(stats_folder, "main_effects_stats_all_variants.csv"), row.names = FALSE)
write.csv(resi_stats_long, file.path(stats_folder, "resi_stats_all_variants.csv"), row.names = FALSE)
write.csv(int_stats_long, file.path(stats_folder, "interaction_stats_all_variants.csv"), row.names = FALSE)

# comment out merging all models because unnecessary and object becomes too big.

# model_files <- list.files(models_folder, pattern = "_models\\.rds$", full.names = TRUE)
# if (length(model_files) == 0) stop("No files found matching *_models.rds in ", models_folder)
# 
# variant_names <- sub("_models\\.rds$", "", basename(model_files))
# all_models <- setNames(lapply(model_files, readRDS), variant_names)
# 
# saveRDS(all_models, file.path(models_folder, "all_variants_models.rds"))
# 
# message("Merged ", length(variant_names), " variant(s): ", paste(variant_names, collapse = ", "))
