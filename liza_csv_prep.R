#comparing Liza's SLIP 25-11 sample and the OG 25-11 sample and modifying liza's output to fit with rest of this pipeline
require(tidyverse)
require(glue)
liza_folder <-  "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11_liza/"

liza_median <- read.csv(paste0(liza_folder,"slip_median_EulerQC_centiles_20260318.csv"))
liza_mpr <- read.csv(paste0(liza_folder,"slip_mpr_EulerQC_centiles_20260318.csv"))
liza_mpr_qc <- read.csv("/mnt/isilon/bgdlab_processing/Liza/SLIP_nnUNet_braincharts/2025_11/code/gamlss/CSV/slip_mpr_EulerQC.csv")
liza_median_qc <- read.csv("/mnt/isilon/bgdlab_processing/Liza/SLIP_nnUNet_braincharts/2025_11/code/gamlss/CSV/slip_median_EulerQC.csv")

slip_folder <-  "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11/"

slip_median <- read.csv(paste0(slip_folder,"centile_csv/slip_median_centiles.csv"))
slip_mpr <- read.csv(paste0(slip_folder,"centile_csv/slip_mpr_centiles.csv"))


#Compare distinct participant IDs
liza_extra_ids <- liza_median$participant_id[!(liza_median$participant_id %in% slip_median$participant_id)]
length(unique(liza_extra_ids))

reg_extra_ids <- slip_median$participant_id[!(slip_median$participant_id %in% liza_median$participant_id)]
length(unique(reg_extra_ids))

# Plot age ranges for extra IDs on both sides
ggplot(data = liza_median %>% filter(participant_id %in% liza_extra_ids) %>%
         distinct(participant_id, .keep_all = TRUE),
       aes(x = age_days/365.25)) +
      geom_histogram(binwidth = 0.1, fill = "grey70", color = "white") +
      labs(
        x = paste0("Age (years),  N = ", length(unique(liza_extra_ids))),
        y = NULL,
        title = "Extra IDs in Liza DF"
      ) +
      theme_minimal(base_size = 16) +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
      )

ggplot(data = slip_median %>% filter(participant_id %in% reg_extra_ids) %>%
         distinct(participant_id, .keep_all = TRUE),
       aes(x = adjusted_age_in_days/365.25)) +
  geom_histogram(binwidth = 0.2, fill = "grey70", color = "white") +
  labs(
    x = paste0("Adjusted Age (years),  N = ", length(unique(reg_extra_ids))),
    y = NULL,
    title = "Extra IDs in Reg DF"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

liza_median <- liza_median %>% rename(adjusted_age_in_days = age_days)
liza_mpr <- liza_mpr %>% rename(adjusted_age_in_days = age_days)

write.csv(liza_median, glue("{liza_folder}centile_csv/liza_median.csv"), row.names=F)
write.csv(liza_mpr, glue("{liza_folder}centile_csv/liza_mpr.csv"), row.names=F)
