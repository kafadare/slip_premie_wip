# get distinct DF

##Subset to one data point per participant, youngest scan available
# cent_df_distinct <- cent_df %>% 
#   group_by(participant_id) %>%
#   slice_min(age_days_adj, n = 1, with_ties = FALSE) %>% 
#   ungroup()

library(gamlssTools)
library(dplyr)
library(glue)

panel_save <- function(plot, filename, extension, folder, plot_width = 8, plot_height = 10){
  save_name <- paste0(folder, filename, ".", extension)
  ggsave(save_name, plot = plot, dpi = 300, width = plot_width, height = plot_height)
}

# Plots of models

model_folder <- "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11_new_models/gamlss_models_rds/"

ctx_model <- readRDS(glue("{model_folder}slip_median_global_CortexVol.rds"))

raw_data_folder <- "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11_new_models/raw_csv/"

slip_median <- read.csv(glue("{raw_data_folder}slip_median_041526.csv")) %>%
  mutate(site = as.factor(site))

(ctx_model_plot <- make_centile_fan(gamlssModel = ctx_model, 
                                    df = slip_median, 
                                    x_var = "log10_age_days", 
                                    desiredCentiles = c(0.05,  0.25, 0.5, 0.75,  0.95),
                                    color_var = "sex"
                                    ) +
    scale_x_continuous(limits = c(log10(0.94*365.25),log10(22.5*365.25)),
                       breaks = c(log10(1*365.25),log10(5*365.25),log10(10*365.25), log10(15*365.25),log10(20*365.25)),
                       labels = c(1,5,10,15,20)) +
    geom_point(alpha = 0.1, size = 5) +
    scale_color_manual(values = c("Male" = "#3F9EE8",
                                  "Female" = "#C486E3")) +
    theme_classic(base_size = 30) +
    labs(x = "Adjusted Age (Years)",
         y = "Cortical Gray Volume",
         title = "") +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        axis.line = element_line(linewidth = 0.4),
        legend.key.size = unit(2, "cm")) +
  guides(shape = guide_legend(override.aes = list(size = 5)),
         color = guide_legend(override.aes = list(size = 5))) )

                     
panel_save(ctx_model_plot, "plots/ctx-model-bysex", "svg", model_folder, plot_width = 14, plot_height = 12)  
