# Regional phenotypes (non-linear models, with interaction term test)

## Volume

### GA

#### fit models
```{r, echo = FALSE}
gabyage_sp_models_vol <- lapply(vol_vars, function(x) {
  form <- as.formula(glue("{x} ~ s(gestational_age) + s(age_days_adj) + ti(gestational_age, age_days_adj) + qc_var"))
  gam_model <- gamm(as.formula(form),
              data = df.zscore)
  return(gam_model)
  })

gabyage_sp_models_vol <- setNames(gabyage_sp_models_vol, vol_vars)

gabyage_sp_gam_table_vol <- as.data.frame(t(sapply(gabyage_sp_models_vol,function(x) {
  summary <- summary(x$gam)
  c(ga_edf = summary$s.table[1,1], ga_p = summary$s.table[1,4],
    age_edf = summary$s.table[2,1], age_p = summary$s.table[2,4],
    gabyage_edf = summary$s.table[3,1], gabyage_p = summary$s.table[3,4])
  })))

rownames(gabyage_sp_gam_table_vol) <- vol_vars

gabyage_sp_gam_table_vol$gabyage_p_adj <- p.adjust(gabyage_sp_gam_table_vol$gabyage_p, method = "BH")

print(paste0("Number of significant IDPs for GA x Age- NO correction: ", dim(gabyage_sp_gam_table_vol %>% filter(gabyage_p < 0.05))[1]))
print(paste0("Number of significant IDPs for GA x Age- AFTER FDR correction: ", dim(gabyage_sp_gam_table_vol %>% filter(gabyage_p_adj < 0.05))[1]))
print(paste0("Significant IDPs for GA x Age- AFTER FDR correction: ", rownames(gabyage_sp_gam_table_vol %>% filter(gabyage_p_adj < 0.05)) ))

gabyage_sp_vol_sig_pheno <- rownames(gabyage_sp_gam_table_vol %>% filter(gabyage_p_adj < 0.05))
```

#### Visualize models
```{r}
gabyage_sp_vol_plot_pheno <- rownames(gabyage_sp_gam_table_vol %>% filter(ga_edf > 2))
font_size <- 16
theme_set(theme_classic(base_family = "sans",base_size = font_size))
line_size <- 1.5
point_size <- 2

mu_age <- mean(cent_all$age_days_adj, na.rm = TRUE)
sd_age <- sd(cent_all$age_days_adj, na.rm = TRUE)

mu_ga <- mean(cent_all$gestational_age, na.rm = TRUE)
sd_ga <- sd(cent_all$gestational_age, na.rm = TRUE)

gabyage_sp_plots_vol <- lapply(gabyage_sp_models_vol[gabyage_sp_vol_plot_pheno], function(x) {
  
  
  low_age_label <- (quantile(df.zscore[!is.na(df.zscore$gestational_age),
                                       "age_days_adj"],
                             probs = c(.10)) * sd_age + mu_age)/365.25
  
  mid_age_label <- (quantile(df.zscore[!is.na(df.zscore$gestational_age),
                                       "age_days_adj"],
                             probs = c(.50)) * sd_age + mu_age)/365.25

  high_age_label <- (quantile(df.zscore[!is.na(df.zscore$gestational_age),
                                      "age_days_adj"],
                            probs = c(.90)) * sd_age + mu_age)/365.25
  
  model_plot <- visualize_model(
  modobj = x,
  smooth_var = "gestational_age",
  int_var = "age_days_adj",
  plabels = "continuous interaction",
  derivative_plot = F)

  plot_colors <- ggplot_build(model_plot)$plot$scales$get_scales("fill")$palette(1)
  
  model_plot <- model_plot + 
    labs(x="Gestational Age (wks)",
         title = glue("{x$gam$terms[[2]]}"),
                  fill = "Age") +
    theme(axis.title.y = element_blank()) +
  scale_x_continuous(labels = function(x) round(x * sd_ga + mu_ga, 1)) +
  scale_fill_manual(labels = c(paste0("Low: ",
                                      round(low_age_label,1),"yrs"),
                               paste0("Mid: ",
                                      round(mid_age_label,1),"yrs"),
                               paste0("High: ",
                                      round(high_age_label,1),"yrs")), 
                    values = plot_colors) +
   scale_color_gradientn(colors = plot_colors,
    labels = function(x) round((x * sd_age + mu_age)/365.25, 1),
    name = "") +
  scale_y_continuous(labels = function(y) signif(pnorm(y)*100, 3))
  
  return(model_plot)
})

gabyage_sp_plots_vol <- setNames(gabyage_sp_plots_vol, gabyage_sp_vol_plot_pheno)

gabyage_sp_plots_vol_nopt <- lapply(gabyage_sp_plots_vol, function(x){
x$layers[1] <- NULL
x})

# wrap_plots(gabyage_sp_plots_vol_nopt, ncol = 2, guides = "collect") &
#   theme(legend.position = "right")
# 
# wrap_plots(gabyage_sp_plots_vol, ncol = 2, guides = "collect") &
#   theme(legend.position = "right")

wrapped_plot_n <- 4 

plot_batches <- split(gabyage_sp_plots_vol, ceiling(seq_along(gabyage_sp_plots_vol) / wrapped_plot_n))

 wrapped_list <- lapply(plot_batches, function(batch) {
   wrap_plots(batch, ncol = 2)
 })

wrapped_list

plot_batches_nopt <- split(gabyage_sp_plots_vol_nopt, ceiling(seq_along(gabyage_sp_plots_vol_nopt) / wrapped_plot_n))

 wrapped_list_nopt <- lapply(plot_batches_nopt, function(batch) {
   wrap_plots(batch, ncol = 2)
 })

wrapped_list_nopt
```

### BW

#### fit models
```{r, echo = FALSE}
bwbyage_sp_models_vol <- lapply(vol_vars, function(x) {
  form <- as.formula(glue("{x} ~ s(birth_weight_kg) + s(age_days_adj) + ti(birth_weight_kg, age_days_adj) + qc_var"))
  gam_model <- gamm(as.formula(form),
              data = df.zscore)
  return(gam_model)
  })

bwbyage_sp_models_vol <- setNames(bwbyage_sp_models_vol, vol_vars)

bwbyage_sp_gam_table_vol <- as.data.frame(t(sapply(bwbyage_sp_models_vol,function(x) {
  summary <- summary(x$gam)
  c(bw_edf = summary$s.table[1,1], bw_p = summary$s.table[1,4],
    age_edf = summary$s.table[2,1], age_p = summary$s.table[2,4],
    bwbyage_edf = summary$s.table[3,1], bwbyage_p = summary$s.table[3,4])
  })))

rownames(bwbyage_sp_gam_table_vol) <- vol_vars

bwbyage_sp_gam_table_vol$bwbyage_p_adj <- p.adjust(bwbyage_sp_gam_table_vol$bwbyage_p, method = "BH")

print(paste0("Number of significant IDPs for bw x Age- NO correction: ", dim(bwbyage_sp_gam_table_vol %>% filter(bwbyage_p < 0.05))[1]))
print(paste0("Number of significant IDPs for bw x Age- AFTER FDR correction: ", dim(bwbyage_sp_gam_table_vol %>% filter(bwbyage_p_adj < 0.05))[1]))
print(paste0("Significant IDPs for bw x Age- AFTER FDR correction: ", rownames(bwbyage_sp_gam_table_vol %>% filter(bwbyage_p_adj < 0.05)) ))

bwbyage_sp_vol_sig_pheno <- rownames(bwbyage_sp_gam_table_vol %>% filter(bwbyage_p_adj < 0.05))
```

#### Visualize models
```{r}
bwbyage_sp_vol_plot_pheno <- rownames(bwbyage_sp_gam_table_vol %>% filter(bw_edf > 2))
font_size <- 16
theme_set(theme_classic(base_family = "sans",base_size = font_size))
line_size <- 1.5
point_size <- 2
source("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/visualize_model_bartfct.R")

mu_age <- mean(cent_all$age_days_adj, na.rm = TRUE)
sd_age <- sd(cent_all$age_days_adj, na.rm = TRUE)

mu_bw <- mean(cent_all$birth_weight_kg, na.rm = TRUE)
sd_bw <- sd(cent_all$birth_weight_kg, na.rm = TRUE)

bwbyage_sp_plots_vol <- lapply(bwbyage_sp_models_vol[bwbyage_sp_vol_plot_pheno], function(x) {
  
  
  low_age_label <- (quantile(df.zscore[!is.na(df.zscore$birth_weight_kg),
                                       "age_days_adj"],
                             probs = c(.10)) * sd_age + mu_age)/365.25
  
  mid_age_label <- (quantile(df.zscore[!is.na(df.zscore$birth_weight_kg),
                                       "age_days_adj"],
                             probs = c(.50)) * sd_age + mu_age)/365.25

  high_age_label <- (quantile(df.zscore[!is.na(df.zscore$birth_weight_kg),
                                      "age_days_adj"],
                            probs = c(.90)) * sd_age + mu_age)/365.25
  
  model_plot <- visualize_model(
  modobj = x,
  smooth_var = "birth_weight_kg",
  int_var = "age_days_adj",
  plabels = "continuous interaction",
  derivative_plot = F)

  plot_colors <- ggplot_build(model_plot)$plot$scales$get_scales("fill")$palette(1)
  
  model_plot <- model_plot + 
    labs(x="Birth Weight (kg)",
         title = glue("{x$gam$terms[[2]]}"),
                  fill = "Age") +
    theme(axis.title.y = element_blank()) +
  scale_x_continuous(labels = function(x) round(x * sd_bw + mu_bw, 1)) +
  scale_fill_manual(labels = c(paste0("Low: ",
                                      round(low_age_label,1),"yrs"),
                               paste0("Mid: ",
                                      round(mid_age_label,1),"yrs"),
                               paste0("High: ",
                                      round(high_age_label,1),"yrs")), 
                    values = plot_colors) +
   scale_color_gradientn(colors = plot_colors,
    labels = function(x) round((x * sd_age + mu_age)/365.25, 1),
    name = "") +
  scale_y_continuous(labels = function(y) signif(pnorm(y)*100, 3))
  
  return(model_plot)
})

bwbyage_sp_plots_vol <- setNames(bwbyage_sp_plots_vol, bwbyage_sp_vol_plot_pheno)

bwbyage_sp_plots_vol_nopt <- lapply(bwbyage_sp_plots_vol, function(x){
x$layers[1] <- NULL
x})

wrapped_plot_n <- 4 

plot_batches <- split(bwbyage_sp_plots_vol, ceiling(seq_along(bwbyage_sp_plots_vol) / wrapped_plot_n))

 wrapped_list <- lapply(plot_batches, function(batch) {
   wrap_plots(batch, ncol = 2)
 })

wrapped_list

plot_batches_nopt <- split(bwbyage_sp_plots_vol_nopt, ceiling(seq_along(bwbyage_sp_plots_vol_nopt) / wrapped_plot_n))

 wrapped_list_nopt <- lapply(plot_batches_nopt, function(batch) {
   wrap_plots(batch, ncol = 2)
 })

wrapped_list_nopt
```


## Surface Area
### GA

#### fit models
```{r, echo = FALSE}
gabyage_sp_models_sa <- lapply(sa_vars, function(x) {
  form <- as.formula(glue("{x} ~ s(gestational_age) + s(age_days_adj) + ti(gestational_age, age_days_adj) + qc_var"))
  gam_model <- gamm(as.formula(form),
              data = df.zscore)
  return(gam_model)
  })

gabyage_sp_models_sa <- setNames(gabyage_sp_models_sa, sa_vars)

gabyage_sp_gam_table_sa <- as.data.frame(t(sapply(gabyage_sp_models_sa,function(x) {
  summary <- summary(x$gam)
  c(ga_edf = summary$s.table[1,1], ga_p = summary$s.table[1,4],
    age_edf = summary$s.table[2,1], age_p = summary$s.table[2,4],
    gabyage_edf = summary$s.table[3,1], gabyage_p = summary$s.table[3,4])
  })))

rownames(gabyage_sp_gam_table_sa) <- sa_vars

gabyage_sp_gam_table_sa$gabyage_p_adj <- p.adjust(gabyage_sp_gam_table_sa$gabyage_p, method = "BH")

print(paste0("Number of significant IDPs for GA x Age- NO correction: ", dim(gabyage_sp_gam_table_sa %>% filter(gabyage_p < 0.05))[1]))
print(paste0("Number of significant IDPs for GA x Age- AFTER FDR correction: ", dim(gabyage_sp_gam_table_sa %>% filter(gabyage_p_adj < 0.05))[1]))
print(paste0("Significant IDPs for GA x Age- AFTER FDR correction: ", rownames(gabyage_sp_gam_table_sa %>% filter(gabyage_p_adj < 0.05)) ))

gabyage_sp_sa_sig_pheno <- rownames(gabyage_sp_gam_table_sa %>% filter(gabyage_p_adj < 0.05))
```

#### Visualize models
```{r}
gabyage_sp_sa_plot_pheno <- rownames(gabyage_sp_gam_table_sa %>% filter(ga_edf > 2))
font_size <- 16
theme_set(theme_classic(base_family = "sans",base_size = font_size))
line_size <- 1.5
point_size <- 2
source("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/visualize_model_bartfct.R")

mu_age <- mean(cent_all$age_days_adj, na.rm = TRUE)
sd_age <- sd(cent_all$age_days_adj, na.rm = TRUE)

mu_ga <- mean(cent_all$gestational_age, na.rm = TRUE)
sd_ga <- sd(cent_all$gestational_age, na.rm = TRUE)

gabyage_sp_plots_sa <- lapply(gabyage_sp_models_sa[gabyage_sp_sa_plot_pheno], function(x) {
  
  
  low_age_label <- (quantile(df.zscore[!is.na(df.zscore$gestational_age),
                                       "age_days_adj"],
                             probs = c(.10)) * sd_age + mu_age)/365.25
  
  mid_age_label <- (quantile(df.zscore[!is.na(df.zscore$gestational_age),
                                       "age_days_adj"],
                             probs = c(.50)) * sd_age + mu_age)/365.25

  high_age_label <- (quantile(df.zscore[!is.na(df.zscore$gestational_age),
                                      "age_days_adj"],
                            probs = c(.90)) * sd_age + mu_age)/365.25
  
  model_plot <- visualize_model(
  modobj = x,
  smooth_var = "gestational_age",
  int_var = "age_days_adj",
  plabels = "continuous interaction",
  derivative_plot = F)

  plot_colors <- ggplot_build(model_plot)$plot$scales$get_scales("fill")$palette(1)
  
  model_plot <- model_plot + 
    labs(x="Gestational Age (wks)",
         title = glue("{x$gam$terms[[2]]}"),
                  fill = "Age") +
    theme(axis.title.y = element_blank()) +
  scale_x_continuous(labels = function(x) round(x * sd_ga + mu_ga, 1)) +
  scale_fill_manual(labels = c(paste0("Low: ",
                                      round(low_age_label,1),"yrs"),
                               paste0("Mid: ",
                                      round(mid_age_label,1),"yrs"),
                               paste0("High: ",
                                      round(high_age_label,1),"yrs")), 
                    values = plot_colors) +
   scale_color_gradientn(colors = plot_colors,
    labels = function(x) round((x * sd_age + mu_age)/365.25, 1),
    name = "") +
  scale_y_continuous(labels = function(y) signif(pnorm(y)*100, 3))
  
  return(model_plot)
})

gabyage_sp_plots_sa <- setNames(gabyage_sp_plots_sa, gabyage_sp_sa_plot_pheno)

gabyage_sp_plots_sa_nopt <- lapply(gabyage_sp_plots_sa, function(x){
x$layers[1] <- NULL
x})

wrapped_plot_n <- 4 

plot_batches <- split(gabyage_sp_plots_sa, ceiling(seq_along(gabyage_sp_plots_sa) / wrapped_plot_n))

 wrapped_list <- lapply(plot_batches, function(batch) {
   wrap_plots(batch, ncol = 2)
 })

wrapped_list

plot_batches_nopt <- split(gabyage_sp_plots_sa_nopt, ceiling(seq_along(gabyage_sp_plots_sa_nopt) / wrapped_plot_n))

 wrapped_list_nopt <- lapply(plot_batches_nopt, function(batch) {
   wrap_plots(batch, ncol = 2)
 })

wrapped_list_nopt

#wrapped_list
```

### BW

#### fit models
```{r, echo = FALSE}
bwbyage_sp_models_sa <- lapply(sa_vars, function(x) {
  form <- as.formula(glue("{x} ~ s(birth_weight_kg) + s(age_days_adj) + ti(birth_weight_kg, age_days_adj) + qc_var"))
  gam_model <- gamm(as.formula(form),
              data = df.zscore)
  return(gam_model)
  })

bwbyage_sp_models_sa <- setNames(bwbyage_sp_models_sa, sa_vars)

bwbyage_sp_gam_table_sa <- as.data.frame(t(sapply(bwbyage_sp_models_sa,function(x) {
  summary <- summary(x$gam)
  c(bw_edf = summary$s.table[1,1], bw_p = summary$s.table[1,4],
    age_edf = summary$s.table[2,1], age_p = summary$s.table[2,4],
    bwbyage_edf = summary$s.table[3,1], bwbyage_p = summary$s.table[3,4])
  })))

rownames(bwbyage_sp_gam_table_sa) <- sa_vars

bwbyage_sp_gam_table_sa$bwbyage_p_adj <- p.adjust(bwbyage_sp_gam_table_sa$bwbyage_p, method = "BH")

print(paste0("Number of significant IDPs for bw x Age- NO correction: ", dim(bwbyage_sp_gam_table_sa %>% filter(bwbyage_p < 0.05))[1]))
print(paste0("Number of significant IDPs for bw x Age- AFTER FDR correction: ", dim(bwbyage_sp_gam_table_sa %>% filter(bwbyage_p_adj < 0.05))[1]))
print(paste0("Significant IDPs for bw x Age- AFTER FDR correction: ", rownames(bwbyage_sp_gam_table_sa %>% filter(bwbyage_p_adj < 0.05)) ))

bwbyage_sp_sa_sig_pheno <- rownames(bwbyage_sp_gam_table_sa %>% filter(bwbyage_p_adj < 0.05))
```

#### Visualize models
```{r}
bwbyage_sp_sa_plot_pheno <- rownames(bwbyage_sp_gam_table_sa %>% filter(bw_edf > 2))
font_size <- 16
theme_set(theme_classic(base_family = "sans",base_size = font_size))
line_size <- 1.5
point_size <- 2
source("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/visualize_model_bartfct.R")

mu_age <- mean(cent_all$age_days_adj, na.rm = TRUE)
sd_age <- sd(cent_all$age_days_adj, na.rm = TRUE)

mu_bw <- mean(cent_all$birth_weight_kg, na.rm = TRUE)
sd_bw <- sd(cent_all$birth_weight_kg, na.rm = TRUE)

bwbyage_sp_plots_sa <- lapply(bwbyage_sp_models_sa[bwbyage_sp_sa_plot_pheno], function(x) {
  
  
  low_age_label <- (quantile(df.zscore[!is.na(df.zscore$birth_weight_kg),
                                       "age_days_adj"],
                             probs = c(.10)) * sd_age + mu_age)/365.25
  
  mid_age_label <- (quantile(df.zscore[!is.na(df.zscore$birth_weight_kg),
                                       "age_days_adj"],
                             probs = c(.50)) * sd_age + mu_age)/365.25

  high_age_label <- (quantile(df.zscore[!is.na(df.zscore$birth_weight_kg),
                                      "age_days_adj"],
                            probs = c(.90)) * sd_age + mu_age)/365.25
  
  model_plot <- visualize_model(
  modobj = x,
  smooth_var = "birth_weight_kg",
  int_var = "age_days_adj",
  plabels = "continuous interaction",
  derivative_plot = F)

  plot_colors <- ggplot_build(model_plot)$plot$scales$get_scales("fill")$palette(1)
  
  model_plot <- model_plot + 
    labs(x="Birth Weight (kg)",
         title = glue("{x$gam$terms[[2]]}"),
                  fill = "Age") +
    theme(axis.title.y = element_blank()) +
  scale_x_continuous(labels = function(x) round(x * sd_bw + mu_bw, 1)) +
  scale_fill_manual(labels = c(paste0("Low: ",
                                      round(low_age_label,1),"yrs"),
                               paste0("Mid: ",
                                      round(mid_age_label,1),"yrs"),
                               paste0("High: ",
                                      round(high_age_label,1),"yrs")), 
                    values = plot_colors) +
   scale_color_gradientn(colors = plot_colors,
    labels = function(x) round((x * sd_age + mu_age)/365.25, 1),
    name = "") +
  scale_y_continuous(labels = function(y) signif(pnorm(y)*100, 3))
  
  return(model_plot)
})

bwbyage_sp_plots_sa <- setNames(bwbyage_sp_plots_sa, bwbyage_sp_sa_plot_pheno)

bwbyage_sp_plots_sa_nopt <- lapply(bwbyage_sp_plots_sa, function(x){
x$layers[1] <- NULL
x})

wrapped_plot_n <- 4 

plot_batches <- split(bwbyage_sp_plots_sa, ceiling(seq_along(bwbyage_sp_plots_sa) / wrapped_plot_n))

 wrapped_list <- lapply(plot_batches, function(batch) {
   wrap_plots(batch, ncol = 2)
 })

wrapped_list

plot_batches_nopt <- split(bwbyage_sp_plots_sa_nopt, ceiling(seq_along(bwbyage_sp_plots_sa_nopt) / wrapped_plot_n))

 wrapped_list_nopt <- lapply(plot_batches_nopt, function(batch) {
   wrap_plots(batch, ncol = 2)
 })

wrapped_list_nopt
```


## Thickness median full sample
### GA

#### fit models
```{r, echo = FALSE}
if (median_th_mpr_flag == TRUE){
gabyage_sp_models_th_median <- lapply(th_vars_median, function(x) {
  form <- as.formula(glue("{x} ~ s(gestational_age) + s(age_days_adj) + ti(gestational_age, age_days_adj) + qc_var"))
  gam_model <- gamm(as.formula(form),
              data = df.zscore)
  return(gam_model)
  })

gabyage_sp_models_th_median <- setNames(gabyage_sp_models_th_median, th_vars)

gabyage_sp_gam_table_th_median <- as.data.frame(t(sapply(gabyage_sp_models_th_median,function(x) {
  summary <- summary(x$gam)
  c(ga_edf = summary$s.table[1,1], ga_p = summary$s.table[1,4],
    age_edf = summary$s.table[2,1], age_p = summary$s.table[2,4],
    gabyage_edf = summary$s.table[3,1], gabyage_p = summary$s.table[3,4])
  })))

rownames(gabyage_sp_gam_table_th_median) <- th_vars

gabyage_sp_gam_table_th_median$gabyage_p_adj <- p.adjust(gabyage_sp_gam_table_th_median$gabyage_p, method = "BH")

print(paste0("Number of significant IDPs for GA x Age- NO correction: ", dim(gabyage_sp_gam_table_th_median %>% filter(gabyage_p < 0.05))[1]))
print(paste0("Number of significant IDPs for GA x Age- AFTER FDR correction: ", dim(gabyage_sp_gam_table_th_median %>% filter(gabyage_p_adj < 0.05))[1]))
print(paste0("Significant IDPs for GA x Age- AFTER FDR correction: ", rownames(gabyage_sp_gam_table_th_median %>% filter(gabyage_p_adj < 0.05)) ))

gabyage_sp_th_median_sig_pheno <- rownames(gabyage_sp_gam_table_th_median %>% filter(gabyage_p_adj < 0.05))
}
```

#### Visualize models
```{r}
if (median_th_mpr_flag == TRUE){
gabyage_sp_th_median_plot_pheno <- rownames(gabyage_sp_gam_table_th_median %>% filter(ga_edf > 2))
font_size <- 16
theme_set(theme_classic(base_family = "sans",base_size = font_size))
line_size <- 1.5
point_size <- 2
source("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/visualize_model_bartfct.R")

mu_age <- mean(cent_all$age_days_adj, na.rm = TRUE)
sd_age <- sd(cent_all$age_days_adj, na.rm = TRUE)

mu_ga <- mean(cent_all$gestational_age, na.rm = TRUE)
sd_ga <- sd(cent_all$gestational_age, na.rm = TRUE)

gabyage_sp_plots_th_median <- lapply(gabyage_sp_models_th_median[gabyage_sp_th_median_sig_pheno], function(x) {
  
  
  low_age_label <- (quantile(df.zscore[!is.na(df.zscore$gestational_age),
                                       "age_days_adj"],
                             probs = c(.10)) * sd_age + mu_age)/365.25
  
  mid_age_label <- (quantile(df.zscore[!is.na(df.zscore$gestational_age),
                                       "age_days_adj"],
                             probs = c(.50)) * sd_age + mu_age)/365.25

  high_age_label <- (quantile(df.zscore[!is.na(df.zscore$gestational_age),
                                      "age_days_adj"],
                            probs = c(.90)) * sd_age + mu_age)/365.25
  
  model_plot <- visualize_model(
  modobj = x,
  smooth_var = "gestational_age",
  int_var = "age_days_adj",
  plabels = "continuous interaction",
  derivative_plot = F)

  plot_colors <- ggplot_build(model_plot)$plot$scales$get_scales("fill")$palette(1)
  
  model_plot <- model_plot + 
    labs(x="Gestational Age (wks)",
         title = glue("{x$gam$terms[[2]]}"),
                  fill = "Age") +
    theme(axis.title.y = element_blank()) +
  scale_x_continuous(labels = function(x) round(x * sd_ga + mu_ga, 1)) +
  scale_fill_manual(labels = c(paste0("Low: ",
                                      round(low_age_label,1),"yrs"),
                               paste0("Mid: ",
                                      round(mid_age_label,1),"yrs"),
                               paste0("High: ",
                                      round(high_age_label,1),"yrs")), 
                    values = plot_colors) +
   scale_color_gradientn(colors = plot_colors,
    labels = function(x) round((x * sd_age + mu_age)/365.25, 1),
    name = "") +
  scale_y_continuous(labels = function(y) signif(pnorm(y)*100, 3))
  
  return(model_plot)
})

gabyage_sp_plots_th_median <- setNames(gabyage_sp_plots_th_median, gabyage_sp_th_median_sig_pheno)

gabyage_sp_plots_th_median_nopt <- lapply(gabyage_sp_plots_th_median, function(x){
x$layers[1] <- NULL
x})

wrapped_plot_n <- 4 

plot_batches <- split(gabyage_sp_plots_th_median, ceiling(seq_along(gabyage_sp_plots_th_median) / wrapped_plot_n))

 wrapped_list <- lapply(plot_batches, function(batch) {
   wrap_plots(batch, ncol = 2)
 })

wrapped_list

plot_batches_nopt <- split(gabyage_sp_plots_th_median_nopt, ceiling(seq_along(gabyage_sp_plots_th_median_nopt) / wrapped_plot_n))

 wrapped_list_nopt <- lapply(plot_batches_nopt, function(batch) {
   wrap_plots(batch, ncol = 2)
 })

wrapped_list_nopt
}
```

### BW

#### fit models
```{r, echo = FALSE}
if (median_th_mpr_flag == TRUE){
bwbyage_sp_models_th_median <- lapply(th_vars_median, function(x) {
  form <- as.formula(glue("{x} ~ s(birth_weight_kg) + s(age_days_adj) + ti(birth_weight_kg, age_days_adj) + qc_var"))
  gam_model <- gamm(as.formula(form),
              data = df.zscore)
  return(gam_model)
  })

bwbyage_sp_models_th_median <- setNames(bwbyage_sp_models_th_median, th_vars)

bwbyage_sp_gam_table_th_median <- as.data.frame(t(sapply(bwbyage_sp_models_th_median,function(x) {
  summary <- summary(x$gam)
  c(bw_edf = summary$s.table[1,1], bw_p = summary$s.table[1,4],
    age_edf = summary$s.table[2,1], age_p = summary$s.table[2,4],
    bwbyage_edf = summary$s.table[3,1], bwbyage_p = summary$s.table[3,4])
  })))

rownames(bwbyage_sp_gam_table_th_median) <- th_vars

bwbyage_sp_gam_table_th_median$bwbyage_p_adj <- p.adjust(bwbyage_sp_gam_table_th_median$bwbyage_p, method = "BH")

print(paste0("Number of significant IDPs for bw x Age- NO correction: ", dim(bwbyage_sp_gam_table_th_median %>% filter(bwbyage_p < 0.05))[1]))
print(paste0("Number of significant IDPs for bw x Age- AFTER FDR correction: ", dim(bwbyage_sp_gam_table_th_median %>% filter(bwbyage_p_adj < 0.05))[1]))
print(paste0("Significant IDPs for bw x Age- AFTER FDR correction: ", rownames(bwbyage_sp_gam_table_th_median %>% filter(bwbyage_p_adj < 0.05)) ))

bwbyage_sp_th_median_sig_pheno <- rownames(bwbyage_sp_gam_table_th_median %>% filter(bwbyage_p_adj < 0.05))
}
```

#### Visualize models
```{r}
if (median_th_mpr_flag == TRUE){
bwbyage_sp_th_median_plot_pheno <- rownames(bwbyage_sp_gam_table_th_median %>% filter(bw_edf > 2))
font_size <- 16
theme_set(theme_classic(base_family = "sans",base_size = font_size))
line_size <- 1.5
point_size <- 2
source("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/visualize_model_bartfct.R")

mu_age <- mean(cent_all$age_days_adj, na.rm = TRUE)
sd_age <- sd(cent_all$age_days_adj, na.rm = TRUE)

mu_bw <- mean(cent_all$birth_weight_kg, na.rm = TRUE)
sd_bw <- sd(cent_all$birth_weight_kg, na.rm = TRUE)

bwbyage_sp_plots_th_median <- lapply(bwbyage_sp_models_th_median[bwbyage_sp_th_median_sig_pheno], function(x) {
  
  
  low_age_label <- (quantile(df.zscore[!is.na(df.zscore$birth_weight_kg),
                                       "age_days_adj"],
                             probs = c(.10)) * sd_age + mu_age)/365.25
  
  mid_age_label <- (quantile(df.zscore[!is.na(df.zscore$birth_weight_kg),
                                       "age_days_adj"],
                             probs = c(.50)) * sd_age + mu_age)/365.25

  high_age_label <- (quantile(df.zscore[!is.na(df.zscore$birth_weight_kg),
                                      "age_days_adj"],
                            probs = c(.90)) * sd_age + mu_age)/365.25
  
  model_plot <- visualize_model(
  modobj = x,
  smooth_var = "birth_weight_kg",
  int_var = "age_days_adj",
  plabels = "continuous interaction",
  derivative_plot = F)

  plot_colors <- ggplot_build(model_plot)$plot$scales$get_scales("fill")$palette(1)
  
  model_plot <- model_plot + 
    labs(x="Birth Weight (kg)",
         title = glue("{x$gam$terms[[2]]}"),
                  fill = "Age") +
    theme(axis.title.y = element_blank()) +
  scale_x_continuous(labels = function(x) round(x * sd_bw + mu_bw, 1)) +
  scale_fill_manual(labels = c(paste0("Low: ",
                                      round(low_age_label,1),"yrs"),
                               paste0("Mid: ",
                                      round(mid_age_label,1),"yrs"),
                               paste0("High: ",
                                      round(high_age_label,1),"yrs")), 
                    values = plot_colors) +
   scale_color_gradientn(colors = plot_colors,
    labels = function(x) round((x * sd_age + mu_age)/365.25, 1),
    name = "") +
  scale_y_continuous(labels = function(y) signif(pnorm(y)*100, 3))
  
  return(model_plot)
})

bwbyage_sp_plots_th_median <- setNames(bwbyage_sp_plots_th_median, bwbyage_sp_th_median_sig_pheno)

bwbyage_sp_plots_th_median_nopt <- lapply(bwbyage_sp_plots_th_median, function(x){
x$layers[1] <- NULL
x})

wrapped_plot_n <- 4 

plot_batches <- split(bwbyage_sp_plots_th_median, ceiling(seq_along(bwbyage_sp_plots_th_median) / wrapped_plot_n))

 wrapped_list <- lapply(plot_batches, function(batch) {
   wrap_plots(batch, ncol = 2)
 })

wrapped_list

plot_batches_nopt <- split(bwbyage_sp_plots_th_median_nopt, ceiling(seq_along(bwbyage_sp_plots_th_median_nopt) / wrapped_plot_n))

 wrapped_list_nopt <- lapply(plot_batches_nopt, function(batch) {
   wrap_plots(batch, ncol = 2)
 })

wrapped_list_nopt
}
```


## Thickness (mpr, df.zscore_overOne overOne sample)
### GA

#### fit models
```{r, echo = FALSE}
gabyage_sp_models_th <- lapply(th_vars, function(x) {
  form <- as.formula(glue("{x} ~ s(gestational_age) + s(age_days_adj) + ti(gestational_age, age_days_adj) + euler_mean_mpr"))
  gam_model <- gamm(as.formula(form),
              data = df.zscore_overOne)
  return(gam_model)
  })

gabyage_sp_models_th <- setNames(gabyage_sp_models_th, th_vars)

gabyage_sp_gam_table_th <- as.data.frame(t(sapply(gabyage_sp_models_th,function(x) {
  summary <- summary(x$gam)
  c(ga_edf = summary$s.table[1,1], ga_p = summary$s.table[1,4],
    age_edf = summary$s.table[2,1], age_p = summary$s.table[2,4],
    gabyage_edf = summary$s.table[3,1], gabyage_p = summary$s.table[3,4])
  })))

rownames(gabyage_sp_gam_table_th) <- th_vars

gabyage_sp_gam_table_th$gabyage_p_adj <- p.adjust(gabyage_sp_gam_table_th$gabyage_p, method = "BH")

print(paste0("Number of significant IDPs for GA x Age- NO correction: ", dim(gabyage_sp_gam_table_th %>% filter(gabyage_p < 0.05))[1]))
print(paste0("Number of significant IDPs for GA x Age- AFTER FDR correction: ", dim(gabyage_sp_gam_table_th %>% filter(gabyage_p_adj < 0.05))[1]))
print(paste0("Significant IDPs for GA x Age- AFTER FDR correction: ", rownames(gabyage_sp_gam_table_th %>% filter(gabyage_p_adj < 0.05)) ))

gabyage_sp_th_sig_pheno <- rownames(gabyage_sp_gam_table_th %>% filter(gabyage_p_adj < 0.05))
```

#### Visualize models
```{r}
gabyage_sp_th_plot_pheno <- rownames(gabyage_sp_gam_table_th %>% filter(ga_edf > 2))
font_size <- 16
theme_set(theme_classic(base_family = "sans",base_size = font_size))
line_size <- 1.5
point_size <- 2
source("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/visualize_model_bartfct.R")

mu_age <- mean(cent_all$age_days_adj, na.rm = TRUE)
sd_age <- sd(cent_all$age_days_adj, na.rm = TRUE)

mu_ga <- mean(cent_all$gestational_age, na.rm = TRUE)
sd_ga <- sd(cent_all$gestational_age, na.rm = TRUE)

gabyage_sp_plots_th <- lapply(gabyage_sp_models_th[gabyage_sp_th_sig_pheno], function(x) {
  
  
  low_age_label <- (quantile(df.zscore_overOne[!is.na(df.zscore_overOne$gestational_age),
                                       "age_days_adj"],
                             probs = c(.10)) * sd_age + mu_age)/365.25
  
  mid_age_label <- (quantile(df.zscore_overOne[!is.na(df.zscore_overOne$gestational_age),
                                       "age_days_adj"],
                             probs = c(.50)) * sd_age + mu_age)/365.25

  high_age_label <- (quantile(df.zscore_overOne[!is.na(df.zscore_overOne$gestational_age),
                                      "age_days_adj"],
                            probs = c(.90)) * sd_age + mu_age)/365.25
  
  model_plot <- visualize_model(
  modobj = x,
  smooth_var = "gestational_age",
  int_var = "age_days_adj",
  plabels = "continuous interaction",
  derivative_plot = F)

  plot_colors <- ggplot_build(model_plot)$plot$scales$get_scales("fill")$palette(1)
  
  model_plot <- model_plot + 
    labs(x="Gestational Age (wks)",
         title = glue("{x$gam$terms[[2]]}"),
                  fill = "Age") +
    theme(axis.title.y = element_blank()) +
  scale_x_continuous(labels = function(x) round(x * sd_ga + mu_ga, 1)) +
  scale_fill_manual(labels = c(paste0("Low: ",
                                      round(low_age_label,1),"yrs"),
                               paste0("Mid: ",
                                      round(mid_age_label,1),"yrs"),
                               paste0("High: ",
                                      round(high_age_label,1),"yrs")), 
                    values = plot_colors) +
   scale_color_gradientn(colors = plot_colors,
    labels = function(x) round((x * sd_age + mu_age)/365.25, 1),
    name = "") +
  scale_y_continuous(labels = function(y) signif(pnorm(y)*100, 3))
  
  return(model_plot)
})

gabyage_sp_plots_th <- setNames(gabyage_sp_plots_th, gabyage_sp_th_sig_pheno)

gabyage_sp_plots_th_nopt <- lapply(gabyage_sp_plots_th, function(x){
x$layers[1] <- NULL
x})

wrapped_plot_n <- 4 

plot_batches <- split(gabyage_sp_plots_th, ceiling(seq_along(gabyage_sp_plots_th) / wrapped_plot_n))

 wrapped_list <- lapply(plot_batches, function(batch) {
   wrap_plots(batch, ncol = 2)
 })

wrapped_list

plot_batches_nopt <- split(gabyage_sp_plots_th_nopt, ceiling(seq_along(gabyage_sp_plots_th_nopt) / wrapped_plot_n))

 wrapped_list_nopt <- lapply(plot_batches_nopt, function(batch) {
   wrap_plots(batch, ncol = 2)
 })

wrapped_list_nopt
```

### BW

#### fit models
```{r, echo = FALSE}
bwbyage_sp_models_th <- lapply(th_vars, function(x) {
  form <- as.formula(glue("{x} ~ s(birth_weight_kg) + s(age_days_adj) + ti(birth_weight_kg, age_days_adj) + euler_mean_mpr"))
  gam_model <- gamm(as.formula(form),
              data = df.zscore_overOne)
  return(gam_model)
  })

bwbyage_sp_models_th <- setNames(bwbyage_sp_models_th, th_vars)

bwbyage_sp_gam_table_th <- as.data.frame(t(sapply(bwbyage_sp_models_th,function(x) {
  summary <- summary(x$gam)
  c(bw_edf = summary$s.table[1,1], bw_p = summary$s.table[1,4],
    age_edf = summary$s.table[2,1], age_p = summary$s.table[2,4],
    bwbyage_edf = summary$s.table[3,1], bwbyage_p = summary$s.table[3,4])
  })))

rownames(bwbyage_sp_gam_table_th) <- th_vars

bwbyage_sp_gam_table_th$bwbyage_p_adj <- p.adjust(bwbyage_sp_gam_table_th$bwbyage_p, method = "BH")

print(paste0("Number of significant IDPs for bw x Age- NO correction: ", dim(bwbyage_sp_gam_table_th %>% filter(bwbyage_p < 0.05))[1]))
print(paste0("Number of significant IDPs for bw x Age- AFTER FDR correction: ", dim(bwbyage_sp_gam_table_th %>% filter(bwbyage_p_adj < 0.05))[1]))
print(paste0("Significant IDPs for bw x Age- AFTER FDR correction: ", rownames(bwbyage_sp_gam_table_th %>% filter(bwbyage_p_adj < 0.05)) ))

bwbyage_sp_th_sig_pheno <- rownames(bwbyage_sp_gam_table_th %>% filter(bwbyage_p_adj < 0.05))
```

#### Visualize models
```{r}
bwbyage_sp_th_plot_pheno <- rownames(bwbyage_sp_gam_table_th %>% filter(bw_edf > 2))
font_size <- 16
theme_set(theme_classic(base_family = "sans",base_size = font_size))
line_size <- 1.5
point_size <- 2
source("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/visualize_model_bartfct.R")

mu_age <- mean(cent_all$age_days_adj, na.rm = TRUE)
sd_age <- sd(cent_all$age_days_adj, na.rm = TRUE)

mu_bw <- mean(cent_all$birth_weight_kg, na.rm = TRUE)
sd_bw <- sd(cent_all$birth_weight_kg, na.rm = TRUE)

bwbyage_sp_plots_th <- lapply(bwbyage_sp_models_th[bwbyage_sp_th_sig_pheno], function(x) {
  
  
  low_age_label <- (quantile(df.zscore_overOne[!is.na(df.zscore_overOne$birth_weight_kg),
                                       "age_days_adj"],
                             probs = c(.10)) * sd_age + mu_age)/365.25
  
  mid_age_label <- (quantile(df.zscore_overOne[!is.na(df.zscore_overOne$birth_weight_kg),
                                       "age_days_adj"],
                             probs = c(.50)) * sd_age + mu_age)/365.25

  high_age_label <- (quantile(df.zscore_overOne[!is.na(df.zscore_overOne$birth_weight_kg),
                                      "age_days_adj"],
                            probs = c(.90)) * sd_age + mu_age)/365.25
  
  model_plot <- visualize_model(
  modobj = x,
  smooth_var = "birth_weight_kg",
  int_var = "age_days_adj",
  plabels = "continuous interaction",
  derivative_plot = F)

  plot_colors <- ggplot_build(model_plot)$plot$scales$get_scales("fill")$palette(1)
  
  model_plot <- model_plot + 
    labs(x="Birth Weight (kg)",
         title = glue("{x$gam$terms[[2]]}"),
                  fill = "Age") +
    theme(axis.title.y = element_blank()) +
  scale_x_continuous(labels = function(x) round(x * sd_bw + mu_bw, 1)) +
  scale_fill_manual(labels = c(paste0("Low: ",
                                      round(low_age_label,1),"yrs"),
                               paste0("Mid: ",
                                      round(mid_age_label,1),"yrs"),
                               paste0("High: ",
                                      round(high_age_label,1),"yrs")), 
                    values = plot_colors) +
   scale_color_gradientn(colors = plot_colors,
    labels = function(x) round((x * sd_age + mu_age)/365.25, 1),
    name = "") +
  scale_y_continuous(labels = function(y) signif(pnorm(y)*100, 3))
  
  return(model_plot)
})

bwbyage_sp_plots_th <- setNames(bwbyage_sp_plots_th, bwbyage_sp_th_sig_pheno)

bwbyage_sp_plots_th_nopt <- lapply(bwbyage_sp_plots_th, function(x){
x$layers[1] <- NULL
x})

wrapped_plot_n <- 4 

plot_batches <- split(bwbyage_sp_plots_th, ceiling(seq_along(bwbyage_sp_plots_th) / wrapped_plot_n))

 wrapped_list <- lapply(plot_batches, function(batch) {
   wrap_plots(batch, ncol = 2)
 })

wrapped_list

plot_batches_nopt <- split(bwbyage_sp_plots_th_nopt, ceiling(seq_along(bwbyage_sp_plots_th_nopt) / wrapped_plot_n))

 wrapped_list_nopt <- lapply(plot_batches_nopt, function(batch) {
   wrap_plots(batch, ncol = 2)
 })

wrapped_list_nopt
```

