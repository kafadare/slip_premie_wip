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
## Writing is merge-safe, not overwrite-safe: each run's output for a given
## (predictor, phenotype_group[, moderator]) replaces just that slice of an
## existing file (if any) and leaves the rest untouched. That means you can
## e.g. run every covariate_group with phenotype_groups=cmd only, after
## already having run the other phenotype groups for those same variants,
## without erasing the earlier results.
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
##  10. phenotype_groups        -- comma-separated phenotype_group selectors to run, or "all"
##                                (choices: global, vol, sa, th, cmd; default "all" if omitted)
##                                "global" expands to every non-CMD global group
##                                (global_self, global_other, global_meanthick); "cmd" expands to
##                                every CMD group (global_cmd, vol_cmd, sa_cmd, th_cmd, subcort_cmd);
##                                "vol"/"sa"/"th" each mean just that one regional group.
##
## Example (one variant per SLURM job, all in parallel):
##   Rscript run_spline_analysis_variant.R $CODE_PATH $VARS_PATH $CENTILES_PATH $OUT_FOLDER gestational_age,bwp_fen eTIV 1
##   Rscript run_spline_analysis_variant.R $CODE_PATH $VARS_PATH $CENTILES_PATH $OUT_FOLDER gestational_age,bwp_fen ADI 1
##   ... (one call per variant, or "all" in a single call to run them serially)
##
## Example (backfilling only the new cmd phenotype groups for variants already run):
##   Rscript run_spline_analysis_variant.R $CODE_PATH $VARS_PATH $CENTILES_PATH $OUT_FOLDER gestational_age,bwp_fen all 0 1 5 cmd

# ---- Packages ----
require(dplyr)
require(glue)
require(purrr)
require(tibble)
require(mgcv)
require(RESI)

# ---- Command line args ----
args <- commandArgs(trailingOnly = TRUE)

# args <- c("/mnt/arcus/lab/users/kafadare/slip_premie_wip/",
#           "/mnt/arcus/lab/users/kafadare/slip_brain_vars/",
#           "/mnt/arcus/lab/users/kafadare/slip_centiles_processed_tables/df_zscore_distinct.csv",
#           "/mnt/arcus/lab/users/kafadare/slip_premie_results/",
#           "gestational_age,bwp_fen,birth_weight_kg",
#           "nocovar",
#           1, 1, 5, "cmd")

code_path <- args[1]
vars_path <- args[2]
qnorm_centiles_path <- args[3]
output_folder <- args[4]
predictors_arg <- if (length(args) >= 5) args[5] else "gestational_age,bwp_fen"
covariate_group_arg <- if (length(args) >= 6) args[6] else "all"
fx <- if (length(args) >= 7) as.logical(as.numeric(args[7])) else FALSE
age_filter_thickness <- if (length(args) >= 8) as.logical(as.numeric(args[8])) else TRUE
rm_sd <- if (length(args) >= 9) as.numeric(args[9]) else 5
phenotype_groups_arg <- if (length(args) >= 10) args[10] else "all"

predictors <- strsplit(predictors_arg, ",")[[1]]

print(list(code_path = code_path, vars_path = vars_path,
           qnorm_centiles_path = qnorm_centiles_path,
           output_folder = output_folder, predictors = predictors,
           covariate_group_arg = covariate_group_arg, fx = fx,
           age_filter_thickness = age_filter_thickness, rm_sd = rm_sd,
           phenotype_groups_arg = phenotype_groups_arg))

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
cmd_vars <- c("global.cmd", "vol.cmd", "sa.cmd", "th.cmd", "subcort.cmd")

df.zscore_full <- read.csv(qnorm_centiles_path)

median_th_mpr_flag <- any(grepl(".median", names(df.zscore_full)))
th_qc_var <- if (median_th_mpr_flag) "euler_mean_mpr" else "qc_var"

# This function removes outliers assuming input is already z-scored. CMD outliers have already been removed, so leave those out of this function.
df.zscore <- remove_outliers(
  df.zscore_full,
  cols = c(global_self_var, global_other_vars, global_meanthick_var,
           vol_vars, sa_vars, th_vars, "bwp_fen"),
  limit = rm_sd
)
# Save information on outliers removed
outlier_summary_path <- file.path(stats_folder, "outlier_removal_summary.csv")

if (!file.exists(outlier_summary_path)) {
  outlier_cols <- intersect(
    c(global_self_var, global_other_vars, global_meanthick_var,
      vol_vars, sa_vars, th_vars, "bwp_fen"),
    names(df.zscore_full)
  )
  
  outlier_summary <- data.frame(
    phenotype = outlier_cols,
    n_total   = sapply(outlier_cols, function(col) sum(!is.na(df.zscore_full[[col]]))),
    n_removed = sapply(outlier_cols, function(col) {
      sum(!is.na(df.zscore_full[[col]]) & is.na(df.zscore[[col]]))
    })
  )
  outlier_summary$pct_removed <- 100 * outlier_summary$n_removed / outlier_summary$n_total
  
  write.csv(outlier_summary, outlier_summary_path, row.names = FALSE)
}
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
    moderators = list(age = list(var = "age_days_adj", type = "continuous"))
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
    moderators = list(age = list(var = "age_days_adj", type = "continuous"))
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
    moderators = list(
      age = list(var = "age_days_adj", type = "continuous"),
      adi = list(var = "ADI.median", type = "continuous")
    )
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
    moderators = list(
      age = list(var = "age_days_adj", type = "continuous"),
      adi = list(var = "ADI.median", type = "continuous")
    )
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
    moderators = list(age = list(var = "age_days_adj", type = "continuous"))
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
    moderators = list(age = list(var = "age_days_adj", type = "continuous"))
  ),
  eTIV_neonatal_dx = list(
    group_covariates = list(
      global_self = c("qc_var", "neonatal_dx"),
      global_other = c("qc_var", "global_eTIV", "neonatal_dx"),
      global_meanthick = c(th_qc_var, "global_eTIV", "neonatal_dx"),
      vol = c("qc_var", "global_eTIV", "neonatal_dx"),
      sa = c("qc_var", "global_SurfaceArea", "neonatal_dx"),
      th = c(th_qc_var, "global_MeanThickness", "neonatal_dx")
    ),
    moderators = list(
      age = list(var = "age_days_adj", type = "continuous")
      #neo_dx = list(var = "neonatal_dx", type = "categorical")
    )
  ),
  neonatal_dx = list(
    group_covariates = list(
      global_self = c("qc_var" , "neonatal_dx"),
      global_other = c("qc_var", "neonatal_dx"),
      global_meanthick = c(th_qc_var, "neonatal_dx"),
      vol = c("qc_var", "neonatal_dx"),
      sa = c("qc_var", "neonatal_dx"),
      th = c(th_qc_var, "neonatal_dx")
    ),
    moderators = list(
      age = list(var = "age_days_adj", type = "continuous")
      #neo_dx = list(var = "neonatal_dx", type = "categorical")
    )
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
    th = list(phenotypes = th_vars, covariates = cov$th, age_filter = age_filter_thickness),
    # CMD = composite/aggregate deviation score per domain (see slipga_qnorm_df_cmd.Rmd);
    # each reuses that domain's own covariate treatment rather than a new one.
    global_cmd = list(phenotypes = "global.cmd", covariates = cov$global_other, age_filter = FALSE),
    vol_cmd = list(phenotypes = "vol.cmd", covariates = cov$vol, age_filter = FALSE),
    sa_cmd = list(phenotypes = "sa.cmd", covariates = cov$sa, age_filter = FALSE),
    th_cmd = list(phenotypes = "th.cmd", covariates = cov$th, age_filter = age_filter_thickness),
    subcort_cmd = list(phenotypes = "subcort.cmd", covariates = cov$vol, age_filter = FALSE)
  )
}

phenotype_group_aliases <- list(
  global = c("global_self", "global_other", "global_meanthick"),
  vol = "vol",
  sa = "sa",
  th = "th",
  cmd = c("global_cmd", "vol_cmd", "sa_cmd", "th_cmd", "subcort_cmd")
)

phenotype_selectors_to_run <- if (identical(phenotype_groups_arg, "all")) {
  names(phenotype_group_aliases)
} else {
  strsplit(phenotype_groups_arg, ",")[[1]]
}
unknown_groups <- setdiff(phenotype_selectors_to_run, names(phenotype_group_aliases))
if (length(unknown_groups) > 0) stop("Unknown phenotype_group(s): ", paste(unknown_groups, collapse = ", "))

phenotype_groups_to_run <- unique(unlist(phenotype_group_aliases[phenotype_selectors_to_run], use.names = FALSE))

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
  groups <- groups[intersect(names(groups), phenotype_groups_to_run)]

  main_stats <- list()
  resi_stats <- list()
  int_stats  <- list()
  models     <- list()

  moderators <- vapply(variants[[variant_name]]$moderators, function(m) m$var, character(1))

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
        mod_spec <- variants[[variant_name]]$moderators[[mod_label]]
        moderator <- mod_spec$var
        interaction_type <- mod_spec$type
        covs_int <- setdiff(g$covariates, moderator)

        m_int <- fit_bam_interactions(data_fit, g$phenotypes, predictor, moderator,
                                       interaction = interaction_type, covariates = covs_int, fx = fx)
        stats_int <- extract_bam_features_interaction(m_int, data_fit, predictor, moderator, covariates = covs_int,
                                                        interaction_type = interaction_type) %>%
          rownames_to_column("phenotype") %>% tag_run(variant_name, predictor, group_id)
        stats_int$moderator_label <- mod_label
        stats_int$moderator <- moderator
        stats_int$interaction_type <- interaction_type

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

# ---- Merge-safe writers: replace only the rows/entries for the
# (predictor, group[, moderator_label]) combinations just computed, leaving
# any other phenotype_groups already present in the file untouched. This is
# what lets you e.g. backfill phenotype_groups=cmd for variants whose other
# phenotype_groups were already run, without erasing those earlier results. --
merge_write_csv <- function(new_df, path, key_cols) {
  if (nrow(new_df) == 0) return(invisible(NULL))

  if (file.exists(path)) {
    existing <- read.csv(path)
    existing_key <- do.call(paste, c(existing[key_cols], sep = "\r"))
    new_key <- do.call(paste, c(new_df[key_cols], sep = "\r"))
    existing <- existing[!(existing_key %in% unique(new_key)), ]
    out <- bind_rows(existing, new_df)
  } else {
    out <- new_df
  }
  write.csv(out, path, row.names = FALSE)
}

merge_write_models <- function(new_models, path) {
  if (length(new_models) == 0) return(invisible(NULL))

  out <- if (file.exists(path)) readRDS(path) else list()
  out[names(new_models)] <- new_models
  saveRDS(out, path)
}

# ---- Run requested variant(s), writing each variant's output as soon as it's
# done (rather than accumulating in memory) so parallel invocations never
# write to the same file, and a long run doesn't lose earlier variants'
# results if a later one fails. ----
for (variant_name in covariate_groups_to_run) {

  message("Running variant: ", variant_name)
  result <- run_variant(variant_name)

  merge_write_models(result$models, file.path(models_folder, paste0(variant_name, "_models.rds")))
  merge_write_csv(result$main_stats, file.path(stats_folder, paste0(variant_name, "_main_effects_stats.csv")),
                  key_cols = c("predictor", "group"))
  merge_write_csv(result$resi_stats, file.path(stats_folder, paste0(variant_name, "_resi_stats.csv")),
                  key_cols = c("predictor", "group"))
  merge_write_csv(result$int_stats, file.path(stats_folder, paste0(variant_name, "_interaction_stats.csv")),
                  key_cols = c("predictor", "group", "moderator_label"))
}
