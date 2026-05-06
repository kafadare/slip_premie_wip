
output_path <- "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/lmsz_reweight_median_Th-MPR_25.11/"

models_vpm_list <- list.files(glue("{output_path}lmsz_models_vpm/"))

out_df_vpm <- data.frame(
  pheno = character(),
  reweight_N = numeric(),
  m_model = character(),
  s_model = character()
)

for (i in 1:length(models_vpm_list)){
  mod_name = models_vpm_list[i]
  mod <- readRDS(glue("{output_path}lmsz_models_vpm/{mod_name}"))
  print(mod_name)
  out_df_vpm[i, "pheno"] <- gsub("lmsz_|.rds", "", mod_name)
  out_df_vpm[i, "reweight_N"] <- sum(mod$weights == 0)
  out_df_vpm[i, "m_model"] <- paste(deparse(mod$mu.formula), collapse = "")
  out_df_vpm[i, "s_model"] <- paste(deparse(mod$sigma.formula), collapse = "")
}

models_lpm_list <- list.files(glue("{output_path}lmsz_models_lpm/"))

out_df_lpm <- data.frame(
  pheno = character(),
  reweight_N = numeric(),
  m_model = character(),
  s_model = character()
)

for (i in 1:length(models_lpm_list)){
  mod_name = models_lpm_list[i]
  mod <- readRDS(glue("{output_path}lmsz_models_lpm/{mod_name}"))
  print(mod_name)
  out_df_lpm[i, "pheno"] <- gsub("lmsz_|.rds", "", mod_name)
  out_df_lpm[i, "reweight_N"] <- sum(mod$weights == 0)
  out_df_lpm[i, "m_model"] <- paste(deparse(mod$mu.formula), collapse = "")
  out_df_lpm[i, "s_model"] <- paste(deparse(mod$sigma.formula), collapse = "")
}

hist(out_df_vpm$reweight_N)

hist(out_df_lpm$reweight_N)
sum(out_df_lpm$reweight_N > 0)
