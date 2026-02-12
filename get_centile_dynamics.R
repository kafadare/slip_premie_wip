require(dplyr)
require(glue)
require(magrittr)
require(gamlssTools)
# For figures/plotting
require(ggplot2)
require(ggseg)
require(egg)
require(patchwork)

# script to find the plateau point for eTIV (and whatever other phenotype of interest)
# set folder paths
code_path <- "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/"
source(glue("{code_path}R_util_lmsz.R"))

vars_path <-  "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/vars/"

centiles_data_path = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/median-Th-MPR-distinct-allgrades-25.11/data_clean/median_Th-MPR_2026-02-04.csv"
raw_data_path = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11/raw_csv/median_Th-MPR_2026-02-04.csv"

models_path = "/mnt/isilon/bgdlab_processing/braincharts/SLIP/2025_03/code/gamlss/RDS/"

csv_save_path = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/"

full_vars <- read.csv(paste0(vars_path,"all_vars.csv")) %>% pull(x)
global_vars <- read.csv(paste0(vars_path,"global_vars.csv")) %>% pull(x)
vol_vars <- read.csv(paste0(vars_path,"vol_vars.csv")) %>% pull(x)
sa_vars <- read.csv(paste0(vars_path,"sa_vars.csv")) %>% pull(x)
th_vars <- read.csv(paste0(vars_path,"th_vars.csv")) %>% pull(x)

df <- read.csv(centiles_data_path) %>% select(-"X")

df_raw <- read.csv(raw_data_path) %>% select(all_of(c(full_vars, "participant_id", "session_id", "sex", "age_days",
                                                      "site", "avg_grade", "age_at_scan"))) %>%
  rename_with(~paste0(.x, ".raw"), all_of(full_vars))

shared_names <- intersect(names(df_raw), names(df))

df <- merge(df, df_raw, by = shared_names, all.x = FALSE, all.y = FALSE)

sim_df <- rbind(expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                            sex = c("Male")),
                expand.grid(adjusted_age_days_log = seq(log(280),log(365.25*21+280),0.01),
                            sex = c("Female")))
NEWData <- sim_df %>%
  select(adjusted_age_days_log, sex) %>% 
  rename(AgeTransformed = adjusted_age_days_log) %>%
  mutate(sex = factor(sex, levels = c("Female","Male"))) %>%
  distinct() %>%
  as.data.frame()

out_df <- data.frame(
  pheno = character(),
  cent_curve = numeric(),
  cent_deriv_adj_age_days_log = numeric(),
  cent_deriv2_adj_age_days_log = numeric(),
  cent_deriv_adj_age_years = numeric(),
  cent_deriv2_adj_age_years = numeric(),
  adj_age_days_log_peak50 = numeric(),
  rank_peak50 = numeric(),
  rank_centderiv = numeric(),
  rank_centderiv2 = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:length(full_vars)) {
  pheno <- full_vars[i]
  print(pheno)
  out_df[i, "pheno"] <- pheno
  
  if (grepl("Thick",pheno)) {
    orig_fit <- readRDS(glue("{models_path}mpr-{pheno}Transformed/FIT.EXTRACT.rds"))
  } else {
    orig_fit <- readRDS(glue("{models_path}median-{pheno}Transformed/FIT.EXTRACT.rds"))
  }
  
  tmp <- Apply.Param(NEWData, orig_fit$param,
                               Pred.Set=c("l025"=0.025,"l250"=0.250,
                                          "m500"=0.5,"u750"=0.750,
                                          "u975"=0.975)) %>%
    pivot_longer(c(PRED.l025.pop,PRED.l250.pop,
                   PRED.m500.pop,PRED.u750.pop,
                   PRED.u975.pop),
                 names_to = "centile_cat", values_to = "SLIP_fans") %>%
    mutate(centile_cat = case_match(centile_cat,
                                    "PRED.l025.pop" ~ "cent_0.025",
                                    "PRED.l250.pop" ~ "cent_0.25",
                                    "PRED.m500.pop" ~ "cent_0.5",
                                    "PRED.u750.pop" ~ "cent_0.75",
                                    "PRED.u975.pop" ~ "cent_0.975")) %>%
    select(AgeTransformed, sex, centile_cat, SLIP_fans) %>%
    rename(adjusted_age_days_log = AgeTransformed)
  
  tmp_50 <- tmp %>% filter(centile_cat == "cent_0.5") %>%
    group_by(adjusted_age_days_log) %>%
    summarise(sex_average = mean(SLIP_fans)) %>%
    ungroup() %>%
    mutate(adjusted_age_years = exp(adjusted_age_days_log)/365.25)
  
  dx <- diff(tmp_50$adjusted_age_days_log)
  dy <- diff(tmp_50$sex_average)
  tmp_50$cent_deriv <- c(99999, (dy / dx))
  tmp_50$cent_deriv2 <- c(99999, 99999, diff(dy / dx))
  
  #plot(tmp_50$adjusted_age_years, tmp_50$sex_average)
  #plot(tmp_50[2:nrow(tmp_50),]$adjusted_age_years, tmp_50[2:nrow(tmp_50),]$cent_deriv)
  #plot(tmp_50[3:nrow(tmp_50),]$adjusted_age_years, tmp_50[3:nrow(tmp_50),]$cent_deriv2)
  
  out_df[i, "cent_deriv_adj_age_days_log"] <- tmp_50[which.min(tmp_50$cent_deriv), "adjusted_age_days_log"]
  out_df[i, "cent_deriv2_adj_age_days_log"] <- tmp_50[which.min(tmp_50$cent_deriv2), "adjusted_age_days_log"]
  out_df[i, "cent_deriv_adj_age_years"] <- tmp_50[which.min(tmp_50$cent_deriv), "adjusted_age_years"]
  out_df[i, "cent_deriv2_adj_age_years"] <- tmp_50[which.min(tmp_50$cent_deriv2), "adjusted_age_years"]
  
  age_at_peak <- tmp_50[which.max(tmp_50$sex_average), "adjusted_age_days_log"]
  out_df[i,"adj_age_days_log_peak50"] <- age_at_peak
}

out_df$cent_curve <- 0.5
out_df <- out_df %>%
  mutate(adj_age_years_peak50 = exp(adj_age_days_log_peak50)/365.25)

# add columns for spin testing and plotting
out_df$hemisphere <- ifelse(grepl("lh", out_df$pheno), "L", 
                            ifelse(grepl("rh", out_df$pheno), "R", NA))
out_df$domain <- ifelse(grepl("aparc_GrayVol", out_df$pheno), "vol", 
                        ifelse(grepl("aparc_SurfArea", out_df$pheno), "sa",
                               ifelse(grepl("aparc_ThickAvg", out_df$pheno), "th",
                                      ifelse(grepl("global", out_df$pheno), "global", 
                                             ifelse(grepl("aseg", out_df$pheno), "aseg",NA)))))

out_df$label <- sub("^(aparc_(GrayVol|SurfArea|ThickAvg)_(lh|rh)_)", "", out_df$pheno)


out_df %>% filter(domain == "vol") %>% frank(adj_age_days_log_peak50, ties.method = "dense")

out_df[which(out_df$domain == "vol"),]$rank_peak50 <- out_df %>% filter(domain == "vol") %>% frank(adj_age_days_log_peak50, ties.method = "dense")
out_df[which(out_df$domain == "sa"),]$rank_peak50 <- out_df %>% filter(domain == "sa") %>% frank(adj_age_days_log_peak50, ties.method = "dense")
out_df[which(out_df$domain == "th"),]$rank_peak50 <- out_df %>% filter(domain == "th") %>% frank(adj_age_days_log_peak50, ties.method = "dense")

out_df[which(out_df$domain == "vol"),]$rank_centderiv <- out_df %>% filter(domain == "vol") %>% frank(cent_deriv_adj_age_days_log, ties.method = "dense")
out_df[which(out_df$domain == "sa"),]$rank_centderiv <- out_df %>% filter(domain == "sa") %>% frank(cent_deriv_adj_age_days_log, ties.method = "dense")
out_df[which(out_df$domain == "th"),]$rank_centderiv <- out_df %>% filter(domain == "th") %>% frank(cent_deriv_adj_age_days_log, ties.method = "dense")

out_df[which(out_df$domain == "vol"),]$rank_centderiv2 <- out_df %>% filter(domain == "vol") %>% frank(cent_deriv2_adj_age_days_log, ties.method = "dense")
out_df[which(out_df$domain == "sa"),]$rank_centderiv2 <- out_df %>% filter(domain == "sa") %>% frank(cent_deriv2_adj_age_days_log, ties.method = "dense")
out_df[which(out_df$domain == "th"),]$rank_centderiv2 <- out_df %>% filter(domain == "th") %>% frank(cent_deriv2_adj_age_days_log, ties.method = "dense")


filename <- glue(("{csv_save_path}centile_peaks.csv"))
write.csv(out_df, filename)


#PLOTS
curves <- list()
brain_plot_df <- data.frame(
  pheno = full_vars,
  value = rep(NA, length(full_vars)),
  label = sub("(aparc_(GrayVol|SurfArea|ThickAvg)|aseg)_", "", full_vars),
  hemi = rep(NA, length(full_vars))
)

rows_edit <- grepl("aseg(?!_CC)", brain_plot_df$pheno, perl = TRUE)

brain_plot_df$label <- brain_plot_df$label %>% 
  sub("X", "x", .) %>%
  sub("rd_Ventricle", "rd_ventricle", .) %>%
  sub("th_Ventricle", "th_ventricle", .) %>% 
  sub("Thalamus", "Thalamus-Proper", .)

brain_plot_df$label[rows_edit] <- gsub("_", "-", brain_plot_df$label[rows_edit])

brain_plot_df$hemi <- ifelse(grepl("lh", brain_plot_df$pheno), "left",
                          ifelse(grepl("rh", brain_plot_df$pheno), "right", NA))

limits = c(0,1)
high_color = "#0F766E"

# plotting curves for inspection
for (i in 1:length(full_vars)) {
  pheno <- full_vars[i]
  print(pheno)
  pheno_title <- sub("(aparc|aseg|global)_", "", pheno)

  if (grepl("Thick",pheno)) {
    orig_fit <- readRDS(glue("{models_path}mpr-{pheno}Transformed/FIT.EXTRACT.rds"))
  } else {
    orig_fit <- readRDS(glue("{models_path}median-{pheno}Transformed/FIT.EXTRACT.rds"))
  }
  
  tmp <- Apply.Param(NEWData, orig_fit$param,
                     Pred.Set=c("l025"=0.025,"l250"=0.250,
                                "m500"=0.5,"u750"=0.750,
                                "u975"=0.975)) %>%
    pivot_longer(c(PRED.l025.pop,PRED.l250.pop,
                   PRED.m500.pop,PRED.u750.pop,
                   PRED.u975.pop),
                 names_to = "centile_cat", values_to = "SLIP_fans") %>%
    mutate(centile_cat = case_match(centile_cat,
                                    "PRED.l025.pop" ~ "cent_0.025",
                                    "PRED.l250.pop" ~ "cent_0.25",
                                    "PRED.m500.pop" ~ "cent_0.5",
                                    "PRED.u750.pop" ~ "cent_0.75",
                                    "PRED.u975.pop" ~ "cent_0.975")) %>%
    select(AgeTransformed, sex, centile_cat, SLIP_fans) %>%
    rename(adjusted_age_days_log = AgeTransformed)
  
  tmp_50 <- tmp %>% filter(centile_cat == "cent_0.5") %>%
    group_by(adjusted_age_days_log) %>%
    summarise(sex_average = mean(SLIP_fans)) %>%
    ungroup() %>%
    mutate(adjusted_age_years = exp(adjusted_age_days_log)/365.25)
  
  dx <- diff(tmp_50$adjusted_age_days_log)
  dy <- diff(tmp_50$sex_average)
  tmp_50$cent_deriv <- c(99999, (dy / dx))
  tmp_50$cent_deriv2 <- c(99999, 99999, diff(dy / dx))
  
  i_peak50   <- which.max(tmp_50$sex_average)
  i_d1    <- which.min(tmp_50$cent_deriv)
  i_d2    <- which.min(tmp_50$cent_deriv2)
  
  p1 <- ggplot(tmp_50, aes(x = adjusted_age_years, y = sex_average)) +
    geom_point(color = "#6B7C8E", alpha = 0.7, size = 1.5) +
      geom_point(data = tmp_50[i_peak50, ], aes(adjusted_age_years, sex_average, color = "50th centile peak"),
                 size = 2.5, shape = 8) +
      geom_point(data = tmp_50[i_d1, ], aes(adjusted_age_years, sex_average, color = "First Deriv Min"),
                  size = 2.5, shape = 8) +
      geom_point(data = tmp_50[i_d2, ], aes(adjusted_age_years, sex_average, color = "Second Deriv Min"),
                 size = 2.5, shape = 8) +
    scale_color_manual(
      values = c("50th centile peak" = "#1C1C1C",
                 "First Deriv Min"   = "#1F60B2",
                 "Second Deriv Min"  = "#1FAF9D")) +
      labs(y = "Sex Averaged 50th Centile") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
            axis.text.x  = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "none",
            axis.title.y = element_text(margin = margin(r = 13), hjust = 0.5))
  
  p2 <- ggplot(tmp_50[2:nrow(tmp_50),], aes(x = adjusted_age_years, y = cent_deriv)) +
    geom_point(color = "#4C72B0", alpha = 0.7, size = 1.5) +
      geom_point(data = tmp_50[i_peak50, ], aes(adjusted_age_years, cent_deriv, color = "50th centile peak"),
                 size = 2.5, shape = 8) +
      geom_point(data = tmp_50[i_d1, ], aes(adjusted_age_years, cent_deriv, color = "First Deriv Min"),
                 size = 2.5, shape = 8) +
      geom_point(data = tmp_50[i_d2, ], aes(adjusted_age_years, cent_deriv, color = "Second Deriv Min"),
                 size = 2.5, shape = 8) +
    scale_color_manual(
      values = c("50th centile peak" = "#1C1C1C",
                 "First Deriv Min"   = "#1F60B2",
                 "Second Deriv Min"  = "#1FAF9D")) +
    labs(y = "First Derivative") +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(margin = margin(r = 13), hjust = 0.5)) +
    theme(legend.position = "right")
  
  p3 <- ggplot(tmp_50[3:nrow(tmp_50),], aes(x = adjusted_age_years, y = cent_deriv2)) +
    geom_point(color = "#55A9A3", alpha = 0.7, size = 1.5) +
      geom_point(data = tmp_50[i_peak50, ], aes(adjusted_age_years, cent_deriv2, color = "50th centile peak"),
                 size = 2.5, shape = 8) +
      geom_point(data = tmp_50[i_d1, ], aes(adjusted_age_years, cent_deriv2, color = "First Deriv Min"),
                 size = 2.5, shape = 8) +
      geom_point(data = tmp_50[i_d2, ], aes(adjusted_age_years, cent_deriv2, color = "Second Deriv Min"),
                 size = 2.5, shape = 8) +
    scale_color_manual(
      values = c("50th centile peak" = "#1C1C1C",
                 "First Deriv Min"   = "#1F60B2",
                 "Second Deriv Min"  = "#1FAF9D")) +
    labs(x = "Adjusted Age (Years)",
      y = "Second Derivative") +
    theme_classic() +
    theme(axis.title.y = element_text(margin = margin(r = 10), hjust = 0.5),
          legend.position = "none")
  
  all_plots <- (p1 / p2 / p3) +
    plot_annotation(title = pheno_title) &
    theme(
      plot.title.position = "plot",
      plot.title = element_text(hjust = 0.5)
    )
  
  plot_df <- brain_plot_df %>% select(pheno, value, label)
  pheno_base <- sub("aparc_(GrayVol|SurfArea|ThickAvg)_", "", pheno)
  plot_df[grepl(pheno_base, plot_df$pheno),"value"] <- 1
  
  if(grepl("aparc", pheno)) {
    plot_data_dkt <- merge(plot_df, ggseg::dk, by = "label")
    hemi_plot <- plot_data_dkt[which.max(plot_data_dkt$value), "hemi"]
    side_plot <- plot_data_dkt[which.max(plot_data_dkt$value), "side"]
  
    brain_plot <- ggseg(plot_data_dkt, 
                         atlas = dk, 
                         colour = "black",
                         size = .2,
                   hemisphere = hemi_plot,
                   view = side_plot,
                   mapping = aes(fill = value)) +
        scale_fill_gradient2(high = high_color, limits = limits) +
        theme_void()+
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none")
    
    all_plots <- all_plots +  inset_element(brain_plot, left = 0.95, bottom = 0.2, right = 1.5, top = 1, align_to = "plot")
    
  
  } else if (grepl("aseg", pheno)){ 
    plot_data_aseg <- merge(plot_df, ggseg::aseg, by = "label")
    brain_plot <- ggplot(plot_data_aseg) +
       geom_brain(atlas = aseg,
                  aes(fill = value)) +
       scale_fill_gradient2(high = high_color, limits = limits) +
       theme_void()+
       theme(axis.title.x = element_blank(),
             axis.title.y = element_blank(),
             legend.position = "none")
    
    all_plots <- all_plots +  inset_element(brain_plot, left = 1.02, bottom = 0.2, right = 1.4, top = 1, align_to = "plot")
    
    # side_plot <- plot_data_aseg[which.max(plot_data_aseg$value), "side"]
    # 
    # (brain_plot <- ggseg(plot_data_aseg, 
    #                      atlas = aseg, 
    #                      colour = "black",
    #                      size = .2,
    #                      view = side_plot,
    #                      mapping = aes(fill = value)) +
    #     scale_fill_gradient2(high = high_color, limits = limits) +
    #     theme_void()+
    #     theme(axis.title.x = element_blank(),
    #           axis.title.y = element_blank(),
    #           legend.position = "none"))
  }
   
   curves[[pheno_title]] <- all_plots
   curves[[pheno_title]]
}

curves_filename <- glue(("{csv_save_path}centile_peaks_plots.rds"))
saveRDS(curves, curves_filename)
