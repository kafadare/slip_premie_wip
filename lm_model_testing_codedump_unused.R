# more efficient and clean way to test a bunch of lm models
# probably can be applied to lme (??)
# copied over from arcus
# can make this into a function
# eventually put all functions together ??

vars_totest <- c("global_eTIV", "global_SubCortGrayVol","global_CortexVol", "global_CerebralWhiteMatterVol", "global_VentricleChoroidVol","global_SurfaceArea", "global_MeanThickness") 
swyc_forms_totest <- c("Month-6", "Month-9", "Month-15", "Month-18", "Month-24", "Month-30")

grid <- expand.grid(
  pheno = vars_totest,
  swyc_form = swyc_forms_totest,
  ga_adj = c(TRUE, FALSE),
  stringsAsFactors = FALSE
)

df <- df.zscore.swyc.all.cent %>% 
  select(pat_id, ga_days_total_qc, birth_weight_kg_qc, vars_totest, 
         form_friendly_name, norm_devo_quotient, ga_adj_norm_devo_quotient_round, 
         median_euler_mean, euler_mean_mpr,complete_centile_vars) %>%
  rename(swyc_form = form_friendly_name) %>%
  mutate(swyc_form = recode(swyc_form,
                            "SWYC : 04 Month" = "Month-4",
                            "SWYC : 06 Month" = "Month-6",
                            "SWYC : 09 Month" = "Month-9",
                            "SWYC : 12 Month" = "Month-12",
                            "SWYC : 15 Month" = "Month-15",
                            "SWYC : 18 Month" = "Month-18",
                            "SWYC : 24 Month" = "Month-24",
                            "SWYC : 30 Month" = "Month-30",
                            "SWYC : 36 Month" = "Month-36",
                            "SWYC : 48 Month" = "Month-48",
                            "SWYC : 60 Month" = "Month-60")) %>%
  pivot_longer(cols = c(ga_adj_norm_devo_quotient_round, norm_devo_quotient),
               names_to = "ga_adj",
               values_to = "dq") %>%
  mutate(ga_adj = grepl("ga_adj", ga_adj))

results <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  
  pheno_i <- grid$pheno[i]
  form_i  <- grid$swyc_form[i]
  ga_i    <- grid$ga_adj[i]
  
  df_sub <- df %>%
    dplyr::filter(swyc_form == form_i, ga_adj == ga_i) %>%
    dplyr::select(all_of(c("dq", pheno_i, "median_euler_mean"))) %>%
    na.omit()
  
  if (nrow(df_sub) < 3) {
    return(data.frame(
      pheno = pheno_i,
      swyc_type = form_i,
      ga_adj = ga_i,
      N = nrow(df_sub),
      beta = NA,
      p = NA,
      se = NA
    ))
  }
  
  fit <- lm(reformulate(c(pheno_i, "median_euler_mean"), response = "dq"),
            data = df_sub)
  
  coef_tab <- summary(fit)$coefficients
  
  data.frame(
    pheno = pheno_i,
    swyc_form = form_i,
    ga_adj = ga_i,
    N = nrow(df_sub),
    beta = coef_tab[pheno_i, "Estimate"],
    p = coef_tab[pheno_i, "Pr(>|t|)"],
    se = coef_tab[pheno_i, "Std. Error"],
    p_adj_pheno = p.adjust(coef_tab[pheno_i, "Pr(>|t|)"], method = "bonferroni", n = length(vars_totest)),
    p_adj_phenoform = p.adjust(coef_tab[pheno_i, "Pr(>|t|)"], method = "bonferroni", n = length(vars_totest)*length(swyc_forms_totest))
  )
}))

# arrange by swyc form and add columns for plotting
results <- results %>%
  mutate(swyc_form = factor(swyc_form, levels = c("Month-6", "Month-9", "Month-15", "Month-18", "Month-24", "Month-30"))) %>% arrange(swyc_form) %>%
  mutate(sig = p_adj_phenoform < 0.05,
         color = case_when(
           sig & beta > 0 ~ "#c87a7a",   #  muted red
           sig & beta < 0 ~ "#7293b0",   #  muted blue
           TRUE ~ "gray70"),
         y_labels = paste(swyc_form, N, sep = ", N="))


(swyc_allforms_global_plot <- ggplot(results, aes(x = beta, y = swyc_form, 
                                                  color = ga_adj, alpha = p_adj_pheno <= 0.05)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbarh(aes(xmin = beta - se, xmax = beta + se),
                   height = 0.2,
                   position = position_dodge(width = 0.5)) +
    facet_wrap(~pheno, scales = "free_x") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#B0B0B0") +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.border = element_blank(),
      axis.text.y  = element_text(color = "black"),
      axis.line.x  = element_line(color = "black"),
      axis.line.y  = element_blank(),
      legend.justification = c("left", "top")
    ) +
    xlab("Std Beta DQ ~ Global Centiles, Across Forms > N = 40") +
    ylab(""))


# if(save_figs == TRUE){
# panel_save(swyc_bf_all_plot, "swyc-bf-all", "pdf", plot_output_folder)
# panel_save(swyc_bf_median_plot, "swyc-bf-median", "pdf", plot_output_folder)
# panel_save(swyc_bf_allforms_ga, "swyc-bf-allforms-ga", "pdf", plot_output_folder)
# }