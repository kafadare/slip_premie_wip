library(irr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr)
library(glue)
library(patchwork)

panel_save <- function(plot, filename, extension, folder, plot_width = 8, plot_height = 10){
  save_name <- paste0(folder, filename, ".", extension)
  ggsave(save_name, plot = plot, dpi = 300, width = plot_width, height = plot_height)
}

# Compare LMSz-centiles to original centiles

## VPM
lmsz_dummy_folder <- "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/lmsz_median_Th-MPR_25.11/lmsz_dummy_models/"
out_folder <- "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/lmsz_median_Th-MPR_25.11/"
plot_save_path <- paste0(out_folder, "plots/")
if (!dir.exists(plot_save_path)) dir.create(plot_save_path)

df_vpm_path <- glue("{out_folder}results_csv/lmsz_adjusted_centiles_vpm.csv")
df_lpm_path <- glue("{out_folder}results_csv/lmsz_adjusted_centiles_lpm.csv")                 
vars_path <-  "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/vars/"
full_vars <- read.csv(paste0(vars_path,"all_vars.csv")) %>% pull(x)

lmsz_vars <- full_vars[!(full_vars %in% c("aseg_CC_Central", "aseg_CC_Mid_Posterior"))]
lmsz_dummy_files <- list.files(lmsz_dummy_folder, pattern = "lmsz_adj_centiles", full.names = TRUE)
n = length(lmsz_dummy_files) + 2
lmsz_full_files <- (c(lmsz_dummy_files, df_vpm_path, df_lpm_path))

icc_consistency <- data.frame(matrix(nrow = length(lmsz_vars), ncol = n))
rownames(icc_consistency) <- lmsz_vars
colnames(icc_consistency) <- c(paste0("dummy_",1:100), "VPM", "LPM")

icc_agreement <- data.frame(matrix(nrow = length(lmsz_vars), ncol = n))
rownames(icc_agreement) <- lmsz_vars
colnames(icc_agreement) <- c(paste0("dummy_",1:100), "VPM", "LPM")

#could I do this part faster with apply() ???
for (i in 1:n) {
  df <- read_csv(lmsz_full_files[i])
  for (j in 1:length(lmsz_vars)) {
    pheno <- lmsz_vars[j]
    print(pheno)
    pheno_lmsz <- glue("{pheno}.lmsz")
    tmp <- df %>% select(participant_id, adjusted_age_in_days, sex, names(df)[grepl(pheno, names(df))])
    
    # tmp_long <- tmp %>% 
    #   pivot_longer(cols = all_of(c(pheno, glue("{pheno}.raw"), glue("{pheno}.zscore"), glue("{pheno}.lmsz"))),
    #                names_to = "key",
    #                values_to = "value") %>%
    #   mutate(column_key = case_when(
    #     !grepl("\\.", key) ~ "orig",
    #     grepl("\\.raw$", key) ~ "raw",
    #     grepl("\\.zscore$", key) ~ "orig_zscore",
    #     grepl("\\.lmsz$", key) ~ "lmsz"
    #   ))
    
    tmp <- tmp %>%
      mutate(across(c(pheno, pheno_lmsz),
                    qnorm,
                    .names = "{.col}_qnorm"))
    
    # ggplot(tmp, aes(x = !!sym(glue("{pheno}_qnorm")), y = !!sym(glue("{pheno_lmsz}_qnorm")))) +
    #   geom_point(alpha = 0.6, color = "darkgrey") + labs(x = "Original QNorm", y = "VPM LMSz QNorm", title = pheno) +
    #   theme_classic()
    
    tmp_icc <- tmp %>% select(all_of(c(pheno, pheno_lmsz)))
    icc_agree <- icc(tmp_icc, model="oneway", type="agreement", unit="single")$value
    icc_cons <- icc(tmp_icc, model="oneway", type="consistency", unit="single")$value
    icc_consistency[j, i] <- icc_agree
    icc_agreement[j, i] <- icc_cons
  }
}

write.csv(icc_agreement, glue("{out_folder}results_csv/icc_agreement_table.csv"))
write.csv(icc_consistency, glue("{out_folder}results_csv/icc_consistency_table.csv"))


#Get a summary table for plotting
icc_stats <- as.data.frame(t(rbind(icc_a_mean = colMeans(icc_agreement), 
                   icc_a_median = apply(icc_agreement, 2, median),
                   icc_a_sd = apply(icc_agreement, 2, sd),
                   icc_c_mean = colMeans(icc_consistency),
                   icc_c_median = apply(icc_consistency, 2, median),
                   icc_c_sd = apply(icc_consistency, 2, sd))))

#PLOTS

# colors
vpm_color <- "#E69F00"
lpm_color <- "#56B4E9"

# Plot across all phenotypes

## ICC Agreement plots
icc_a_mean_plot <- ggplot(icc_stats, aes(x = icc_a_mean)) + 
  geom_bar(color = "gray70", fill = "gray70") + 
  geom_point(data = icc_stats[rownames(icc_stats) == "VPM",], aes(x = icc_a_mean, y = 1.05), color = vpm_color,
             size = 2, shape = 8) +
  annotate("text", x = icc_stats[rownames(icc_stats) == "VPM", "icc_a_mean"], y = 1.1, color = vpm_color,
           size = 3, label = "VPM") +
  geom_point(data = icc_stats[rownames(icc_stats) == "LPM",], aes(x = icc_a_mean, y = 1.05), color = lpm_color,
             size = 2, shape = 8) +
  annotate("text", x = icc_stats[rownames(icc_stats) == "LPM", "icc_a_mean"], y = 1.1, color = lpm_color,
           size = 3, label = "LPM") +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  ggtitle("ICC-A Mean Across All Phenotypes") +
  theme(plot.title = element_text(hjust = 0.5))

icc_a_median_plot <- ggplot(icc_stats, aes(x = icc_a_median)) + 
  geom_bar(color = "gray70", fill = "gray70") + 
  geom_point(data = icc_stats[rownames(icc_stats) == "VPM",], aes(x = icc_a_median, y = 105), color = vpm_color,
             size = 2, shape = 8) +
  annotate("text", x = icc_stats[rownames(icc_stats) == "VPM", "icc_a_median"], y = 110, color = vpm_color,
           size = 3, label = "VPM") +
  geom_point(data = icc_stats[rownames(icc_stats) == "LPM",], aes(x = icc_a_median, y = 105), color = lpm_color,
             size = 2, shape = 8) +
  annotate("text", x = icc_stats[rownames(icc_stats) == "LPM", "icc_a_median"], y = 110, color = lpm_color,
           size = 3, label = "LPM") +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  ggtitle("ICC-A Median Across All Phenotypes") +
  theme(plot.title = element_text(hjust = 0.5))

icc_a_sd_plot <- ggplot(icc_stats, aes(x = icc_a_sd)) + 
  geom_bar(color = "gray70", fill = "gray70") + 
  geom_point(data = icc_stats[rownames(icc_stats) == "VPM",], aes(x = icc_a_sd, y = 1.05), color = vpm_color,
             size = 2, shape = 8) +
  annotate("text", x = icc_stats[rownames(icc_stats) == "VPM", "icc_a_sd"], y = 1.10, color = vpm_color,
           size = 3, label = "VPM") +
  geom_point(data = icc_stats[rownames(icc_stats) == "LPM",], aes(x = icc_a_sd, y = 1.05), color = lpm_color,
             size = 2, shape = 8) +
  annotate("text", x = icc_stats[rownames(icc_stats) == "LPM", "icc_a_sd"], y = 1.10, color = lpm_color,
           size = 3, label = "LPM") +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  ggtitle("ICC-A SD Across All Phenotypes") +
  theme(plot.title = element_text(hjust = 0.5))

## ICC Consistency Plots
  
icc_c_mean_plot <- ggplot(icc_stats, aes(x = icc_c_mean)) + 
  geom_bar(color = "gray70", fill = "gray70") + 
  geom_point(data = icc_stats[rownames(icc_stats) == "VPM",], aes(x = icc_c_mean, y = 1.05), color = vpm_color,
             size = 2, shape = 8) +
  annotate("text", x = icc_stats[rownames(icc_stats) == "VPM", "icc_c_mean"], y = 1.1, color = vpm_color,
           size = 3, label = "VPM") +
  geom_point(data = icc_stats[rownames(icc_stats) == "LPM",], aes(x = icc_c_mean, y = 1.05), color = lpm_color,
             size = 2, shape = 8) +
  annotate("text", x = icc_stats[rownames(icc_stats) == "LPM", "icc_c_mean"], y = 1.1, color = lpm_color,
           size = 3, label = "LPM") +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  ggtitle("ICC-C Mean Across All Phenotypes") +
  theme(plot.title = element_text(hjust = 0.5))

icc_c_median_plot <- ggplot(icc_stats, aes(x = icc_c_median)) + 
  geom_bar(color = "gray70", fill = "gray70") + 
  geom_point(data = icc_stats[rownames(icc_stats) == "VPM",], aes(x = icc_c_median, y = 105), color = vpm_color,
             size = 2, shape = 8) +
  annotate("text", x = icc_stats[rownames(icc_stats) == "VPM", "icc_c_median"], y = 110, color = vpm_color,
           size = 3, label = "VPM") +
  geom_point(data = icc_stats[rownames(icc_stats) == "LPM",], aes(x = icc_c_median, y = 105), color = lpm_color,
             size = 2, shape = 8) +
  annotate("text", x = icc_stats[rownames(icc_stats) == "LPM", "icc_c_median"], y = 110, color = lpm_color,
           size = 3, label = "LPM") +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  ggtitle("ICC-C Median Across All Phenotypes") +
  theme(plot.title = element_text(hjust = 0.5))

icc_c_sd_plot <- ggplot(icc_stats, aes(x = icc_c_sd)) + 
  geom_bar(color = "gray70", fill = "gray70") + 
  geom_point(data = icc_stats[rownames(icc_stats) == "VPM",], aes(x = icc_c_sd, y = 1.05), color = vpm_color,
             size = 2, shape = 8) +
  annotate("text", x = icc_stats[rownames(icc_stats) == "VPM", "icc_c_sd"], y = 1.10, color = vpm_color,
           size = 3, label = "VPM") +
  geom_point(data = icc_stats[rownames(icc_stats) == "LPM",], aes(x = icc_c_sd, y = 1.05), color = lpm_color,
             size = 2, shape = 8) +
  annotate("text", x = icc_stats[rownames(icc_stats) == "LPM", "icc_c_sd"], y = 1.10, color = lpm_color,
           size = 3, label = "LPM") +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  ggtitle("ICC-C SD Across All Phenotypes") +
  theme(plot.title = element_text(hjust = 0.5))

panel_save(icc_a_mean_plot, "icc-a-mean-plot", "pdf", plot_save_path)
panel_save(icc_a_median_plot, "icc-a-median-plot", "pdf", plot_save_path)
panel_save(icc_a_sd_plot, "icc-a-sd-plot", "pdf", plot_save_path)
panel_save(icc_c_mean_plot, "icc-c-mean-plot", "pdf", plot_save_path)
panel_save(icc_c_median_plot, "icc-c-median-plot", "pdf", plot_save_path)
panel_save(icc_c_sd_plot, "icc-c-sd-plot", "pdf", plot_save_path)
