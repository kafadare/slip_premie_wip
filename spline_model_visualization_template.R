gabyage_sp_plots_vol <- lapply(gabyage_sp_models_vol[gabyage_sp_vol_sig_pheno], function(x) {
  
  low_age_label <- (quantile(df.zscore[!is.na(df.zscore$gestational_age),
                                       "gestational_age"],
                             probs = c(.05)) * sd_age + mu_age)/365.25
  
  high_age_label <- (quantile(df.zscore[!is.na(df.zscore$gestational_age),
                                        "gestational_age"],
                              probs = c(.95)) * sd_age + mu_age)/365.25
  
  model_plot <- visualize_model(
    modobj = x,
    int_var = "gestational_age",
    smooth_var = "adjusted_age_in_days",
    plabels = "continuous interaction",
    derivative_plot = F)
  
  if(grepl("Thick",x$gam$pterms[[2]])) {
    low_age_label <- (quantile(df.zscore_overOne[!is.na(df.zscore_overOne$gestational_age),
                                                 "gestational_age"],
                               probs = c(.05)) * sd_ga + mu_ga)/365.25
    
    high_age_label <- (quantile(df.zscore_overOne[!is.na(df.zscore_overOne$gestational_age),
                                                  "gestational_age"],
                                probs = c(.95)) * sd_ga + mu_ga)/365.25
  }
  plot_colors <- ggplot_build(model_plot)$plot$scales$get_scales("fill")$palette(1)
  
  model_plot <- model_plot + 
    labs(x="Age",
         title = glue("{x$gam$terms[[2]]}"),
         fill = "Age") +
    theme(axis.title.y = element_blank()) +
    scale_x_continuous(labels = function(x) round(x * sd_age + mu_age, 1)) +
    scale_fill_manual(labels = c(paste0("Low: ",
                                        round(low_age_label,1),"yrs"),
                                 paste0("High: ",
                                        round(high_age_label,1),"yrs")), 
                      values = plot_colors) +
    scale_color_gradientn(colors = plot_colors,
                          labels = function(x) round((x * sd_ga + mu_ga), 1),
                          name = "") +
    scale_y_continuous(labels = function(y) signif(pnorm(y)*100, 3))
  
  return(model_plot)
})