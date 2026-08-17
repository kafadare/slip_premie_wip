## Runs the GA/BWP spline analyses (main effects + interactions) for one or
## more covariate-set variants, and saves long-format stat tables + models
## per variant for downstream visualization. Replaces the four separate
## slipga_spline_analysis_*.Rmd files for the model-fitting/stat-extraction
## step; plotting stays in a separate script/Rmd that reads the CSVs saved
## here.
##
## Each variant writes its own files (<variant>_main_effects_stats.csv,
## <variant>_resi_stats.csv, <variant>_interaction_stats.csv,
## <variant>_models.rds), so multiple variants can be run at the same time
## as separate parallel processes (e.g. one Rscript call per variant, or one
## SLURM job per variant via rscript-submit.sh) without collisions, and
## merge_spline_analysis_variants.R combines them afterward.
##
## Command line args, in order:
##   1. code_path             -- folder containing model_fitting_fcts.R etc.
##   2. vars_path              -- folder containing vol_vars.csv/sa_vars.csv/th_vars.csv
##   3. qnorm_centiles_path    -- full path to df_zscore_distinct.csv
##   4. output_folder          -- results are written under <output_folder>/models/ and /results_csv/
##   5. predictors              -- comma-separated predictor names to run
##                                (default "gestational_age,bwp_fen" if omitted)
##   6. covariate_group        -- comma-separated variant names to run, or "all"
##                                (choices: nocovar, eTIV, ADI, eTIV_ADI, wh, eTIV_wh, neonatal_dx, eTIV_neonatal_dx)
##   7. fx                     -- 0/1, passed to fit_bam_models()/fit_bam_interactions() as fx = TRUE/FALSE
##   8. age_filter_thickness   -- 0/1, whether to restrict global_meanthick/th phenotype
##                                groups to age_over1_unadj (default 1/TRUE if omitted)
##   9. rm_sd                   -- outlier cutoff (in SD) passed to remove_outliers()
##                                (default 5 if omitted)
##
## Example (one variant per SLURM job, all in parallel):
##   Rscript run_spline_analysis_variant.R $CODE_PATH $VARS_PATH $CENTILES_PATH $OUT_FOLDER gestational_age,bwp_fen eTIV 1
##   Rscript run_spline_analysis_variant.R $CODE_PATH $VARS_PATH $CENTILES_PATH $OUT_FOLDER gestational_age,bwp_fen ADI 1
##   ... (one call per variant, or "all" in a single call to run them serially)

# ---- Packages ----
require(dplyr)
require(glue)
require(purrr)
require(tibble)
require(mgcv)
require(RESI)

# ---- Command line args ----
args <- commandArgs(trailingOnly = TRUE)

code_path <- args[1]
vars_path <- args[2]
qnorm_centiles_path <- args[3]
output_folder <- args[4]
predictors_arg <- if (length(args) >= 5) args[5] else "gestational_age,bwp_fen"
covariate_group_arg <- if (length(args) >= 6) args[6] else "all"
fx <- if (length(args) >= 7) as.logical(as.numeric(args[7])) else FALSE
age_filter_thickness <- if (length(args) >= 8) as.logical(as.numeric(args[8])) else TRUE
rm_sd <- if (length(args) >= 9) as.numeric(args[9]) else 5

predictors <- strsplit(predictors_arg, ",")[[1]]

print(list(code_path = code_path, vars_path = vars_path,
           qnorm_centiles_path = qnorm_centiles_path,
           output_folder = output_folder, predictors = predictors,
           covariate_group_arg = covariate_group_arg, fx = fx,
           age_filter_thickness = age_filter_thickness, rm_sd = rm_sd))

# ---- Paths ----
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

models_folder <- paste0(output_folder, "models/")
if (!dir.exists(models_folder)) dir.create(models_folder)

stats_folder <- paste0(output_folder, "results_csv/")
if (!dir.exists(stats_folder)) dir.create(stats_folder)

source(glue("{code_path}model_fitting_fcts.R"))
source(glue("{code_path}model_extract_fcts.R"))
source(glue("{code_path}slip_data_utils.R"))

set.seed(42)

# ---- Load data ----
vol_vars <- read.csv(paste0(vars_path, "vol_vars.csv")) %>% pull(x)
sa_vars  <- read.csv(paste0(vars_path, "sa_vars.csv")) %>% pull(x)
th_vars  <- read.csv(paste0(vars_path, "th_vars.csv")) %>% pull(x)

global_self_var  <- "global_eTIV"
global_other_vars <- c("global_SubCortGrayVol", "global_CortexVol",
                        "global_CerebralWhiteMatterVol", "global_VentricleChoroidVol",
                        "global_SurfaceArea")
global_meanthick_var <- "global_MeanThickness"

df.zscore_full <- read.csv(qnorm_centiles_path)

median_th_mpr_flag <- any(grepl(".median", names(df.zscore_full)))
th_qc_var <- if (median_th_mpr_flag) "euler_mean_mpr" else "qc_var"

df.zscore <- remove_outliers(
  df.zscore_full,
  cols = c(global_self_var, global_other_vars, global_meanthick_var,
           vol_vars, sa_vars, th_vars, "bwp_fen"),
  limit = rm_sd
)
rm(df.zscore_full)
gc()

# ---- Variant definitions: each variant fully and explicitly states its own
# covariates per phenotype group (no shared "base" that variants add to or
# subtract from), so combinations like ADI-with-eTIV vs ADI-without-eTIV are
# each a plain, self-contained entry. ----
variants <- list(
  nocovar = list(
    group_covariates = list(
      global_self = c("qc_var"), global_other = c("qc_var"), global_meanthick = c(th_qc_var),
      vol = c("qc_var"), sa = c("qc_var"), th = c(th_qc_var)
    ),
    moderators = c(age = "age_days_adj")
  ),
  eTIV = list(
    group_covariates = list(
      global_self = c("qc_var"),
      global_other = c("qc_var", "global_eTIV"),
      global_meanthick = c(th_qc_var, "global_eTIV"),
      vol = c("qc_var", "global_eTIV"),
      sa = c("qc_var", "global_SurfaceArea"),
      th = c(th_qc_var, "global_MeanThickness")
    ),
    moderators = c(age = "age_days_adj")
  ),
  ADI = list(
    group_covariates = list(
      global_self = c("qc_var", "ADI.median"),
      global_other = c("qc_var", "ADI.median"),
      global_meanthick = c(th_qc_var, "ADI.median"),
      vol = c("qc_var", "ADI.median"),
      sa = c("qc_var", "ADI.median"),
      th = c(th_qc_var, "ADI.median")
    ),
    moderators = c(age = "age_days_adj", adi = "ADI.median")
  ),
  eTIV_ADI = list(
    group_covariates = list(
      global_self = c("qc_var", "ADI.median"),
      global_other = c("qc_var", "global_eTIV", "ADI.median"),
      global_meanthick = c(th_qc_var, "global_eTIV", "ADI.median"),
      vol = c("qc_var", "global_eTIV", "ADI.median"),
      sa = c("qc_var", "global_SurfaceArea", "ADI.median"),
      th = c(th_qc_var, "global_MeanThickness", "ADI.median")
    ),
    moderators = c(age = "age_days_adj", adi = "ADI.median")
  ),
  wh = list(
    group_covariates = list(
      global_self = c("qc_var", "weight_zscore", "height_zscore"),
      global_other = c("qc_var", "weight_zscore", "height_zscore"),
      global_meanthick = c(th_qc_var, "weight_zscore", "height_zscore"),
      vol = c("qc_var", "weight_zscore", "height_zscore"),
      sa = c("qc_var", "weight_zscore", "height_zscore"),
      th = c(th_qc_var, "weight_zscore", "height_zscore")
    ),
    moderators = c(age = "age_days_adj")
  ),
  eTIV_wh = list(
    group_covariates = list(
      global_self = c("qc_var", "weight_zscore", "height_zscore"),
      global_other = c("qc_var", "global_eTIV", "weight_zscore", "height_zscore"),
      global_meanthick = c(th_qc_var, "global_eTIV", "weight_zscore", "height_zscore"),
      vol = c("qc_var", "global_eTIV", "weight_zscore", "height_zscore"),
      sa = c("qc_var", "global_SurfaceArea", "weight_zscore", "height_zscore"),
      th = c(th_qc_var, "global_MeanThickness", "weight_zscore", "height_zscore")
    ),
    moderators = c(age = "age_days_adj")
  ),
  eTIV_neonatal_dx = list(
    group_covariates = list(
      global_self = c("qc_var"),
      global_other = c("qc_var", "global_eTIV", "neonatal_dx"),
      global_meanthick = c(th_qc_var, "global_eTIV", "neonatal_dx"),
      vol = c("qc_var", "global_eTIV", "neonatal_dx"),
      sa = c("qc_var", "global_SurfaceArea", "neonatal_dx"),
      th = c(th_qc_var, "global_MeanThickness", "neonatal_dx")
    ),
    moderators = c(age = "age_days_adj", neo_dx = "neonatal_dx")
  ),
  neonatal_dx = list(
    group_covariates = list(
      global_self = c("qc_var"),
      global_other = c("qc_var", "neonatal_dx"),
      global_meanthick = c(th_qc_var, "neonatal_dx"),
      vol = c("qc_var", "neonatal_dx"),
      sa = c("qc_var", "neonatal_dx"),
      th = c(th_qc_var, "neonatal_dx")
    ),
    moderators = c(age = "age_days_adj", neo_dx = "neonatal_dx")
  )
)

covariate_groups_to_run <- if (identical(covariate_group_arg, "all")) {
  names(variants)
} else {
  strsplit(covariate_group_arg, ",")[[1]]
}
unknown <- setdiff(covariate_groups_to_run, names(variants))
if (length(unknown) > 0) stop("Unknown covariate_group(s): ", paste(unknown, collapse = ", "))

build_phenotype_groups <- function(variant_name) {
  cov <- variants[[variant_name]]$group_covariates
  list(
    global_self = list(phenotypes = global_self_var, covariates = cov$global_self, age_filter = FALSE),
    global_other = list(phenotypes = global_other_vars, covariates = cov$global_other, age_filter = FALSE),
    global_meanthick = list(phenotypes = global_meanthick_var, covariates = cov$global_meanthick, age_filter = age_filter_thickness),
    vol = list(phenotypes = vol_vars, covariates = cov$vol, age_filter = FALSE),
    sa = list(phenotypes = sa_vars, covariates = cov$sa, age_filter = FALSE),
    th = list(phenotypes = th_vars, covariates = cov$th, age_filter = age_filter_thickness)
  )
}

tag_run <- function(df, variant_name, predictor, group_id) {
  df$variant <- variant_name
  df$predictor <- predictor
  df$group <- group_id
  df
}

# ---- Narrow df.zscore down to only the columns a given (predictor, group)
# combination's models actually need before any fitting/filtering happens.
# fit_bam_models()/extract_bam_features_*() filter and copy whatever data
# frame they're given once per phenotype internally, so handing them the full
# wide table (every global/vol/sa/th column) for a 100+-phenotype group means
# 100+ full-width copies for models that each use 4-6 columns. Slicing once
# here keeps every downstream copy small instead. ----
slice_group_data <- function(g, predictor, moderators) {
  needed <- unique(c(g$phenotypes, predictor, g$covariates, moderators))
  if (g$age_filter) needed <- c(needed, "age_over1_unadj")
  d <- df.zscore[, intersect(names(df.zscore), needed), drop = FALSE]
  if (g$age_filter) d <- filter(d, age_over1_unadj == TRUE)
  d
}

# ---- Run one variant: fits every predictor x phenotype_group[ x moderator]
# combination for that variant, and returns the three long-format tables plus
# the nested models list, all scoped to this variant only. ----
run_variant <- function(variant_name) {

  groups <- build_phenotype_groups(variant_name)

  main_stats <- list()
  resi_stats <- list()
  int_stats  <- list()
  models     <- list()

  moderators <- unname(variants[[variant_name]]$moderators)

  for (predictor in predictors) {
    for (group_id in names(groups)) {

      g <- groups[[group_id]]
      data_fit <- slice_group_data(g, predictor, moderators)
      run_key <- paste(predictor, group_id, sep = "_")

      # -- main effects: BH correction (inside extract_bam_features_main/
      # extract_bam_table) is scoped exactly to this group's phenotypes --
      m_main <- fit_bam_models(data_fit, g$phenotypes, predictor, covariates = g$covariates, fx = fx)
      stats_main <- extract_bam_features_main(m_main, data_fit, predictor, covariates = g$covariates) %>%
        rownames_to_column("phenotype") %>% tag_run(variant_name, predictor, group_id)
      resi_main <- extract_bam_resi(m_main, term = paste0("s(", predictor, ")")) %>%
        rownames_to_column("phenotype") %>% tag_run(variant_name, predictor, group_id)

      main_stats[[run_key]] <- stats_main
      resi_stats[[run_key]] <- resi_main

      run_models <- list(main = m_main, interactions = list())

      # -- interactions: one BH-scoped batch per moderator. When the moderator
      # is itself one of this group's covariates (e.g. ADI.median as both a
      # sensitivity covariate and the interaction moderator), drop it from the
      # covariate list so it isn't entered twice (once as s(moderator), once
      # as a bare linear term). --
      for (mod_label in names(variants[[variant_name]]$moderators)) {
        moderator <- variants[[variant_name]]$moderators[[mod_label]]
        covs_int <- setdiff(g$covariates, moderator)

        m_int <- fit_bam_interactions(data_fit, g$phenotypes, predictor, moderator, covariates = covs_int, fx = fx)
        stats_int <- extract_bam_features_interaction(m_int, data_fit, predictor, moderator, covariates = covs_int) %>%
          rownames_to_column("phenotype") %>% tag_run(variant_name, predictor, group_id)
        stats_int$moderator_label <- mod_label
        stats_int$moderator <- moderator

        int_key <- paste(run_key, mod_label, sep = "_")
        int_stats[[int_key]] <- stats_int
        run_models$interactions[[mod_label]] <- m_int
      }

      models[[run_key]] <- run_models

      rm(data_fit, m_main)
      gc()
    }
  }

  list(
    main_stats = bind_rows(main_stats),
    resi_stats = bind_rows(resi_stats),
    int_stats  = bind_rows(int_stats),
    models     = models
  )
}

# ---- Run requested variant(s), writing each variant's output as soon as it's
# done (rather than accumulating in memory) so parallel invocations never
# write to the same file, and a long run doesn't lose earlier variants'
# results if a later one fails. ----
for (variant_name in covariate_groups_to_run) {

  message("Running variant: ", variant_name)
  result <- run_variant(variant_name)

  saveRDS(result$models, file.path(models_folder, paste0(variant_name, "_models.rds")))
  write.csv(result$main_stats, file.path(stats_folder, paste0(variant_name, "_main_effects_stats.csv")), row.names = FALSE)
  write.csv(result$resi_stats, file.path(stats_folder, paste0(variant_name, "_resi_stats.csv")), row.names = FALSE)
  write.csv(result$int_stats, file.path(stats_folder, paste0(variant_name, "_interaction_stats.csv")), row.names = FALSE)
}
