#comparing Liza's SLIP sample and the 2025-03 SLIP sample participants
require(tidyverse)
vars_save_path <-  "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/data_processed_csv/"

liza_median <- read.csv(paste0(vars_save_path,"liza_median_all_clean_2026-01-08.csv")) %>% select(-("X"))
median <- read.csv(paste0(vars_save_path,"median_all_clean_2026-01-08.csv")) %>% select(-("X"))

liza_median_distinct <- read.csv(paste0(vars_save_path,"liza_median_distinct_clean_2026-01-08.csv")) %>% select(-("X"))
median_distinct <- read.csv(paste0(vars_save_path,"median_distinct_clean_2026-01-08.csv")) %>% select(-("X"))


#Compare distinct participant IDs
liza_extra_ids <- liza_median_distinct$participant_id[!(liza_median_distinct$participant_id %in% median_distinct$participant_id)]

reg_extra_ids <- median_distinct$participant_id[!(median_distinct$participant_id %in% liza_median_distinct$participant_id)]


ggplot(data = liza_median_distinct %>% filter(participant_id %in% liza_extra_ids),
       aes(x = adjusted_age_years)) +
      geom_histogram(binwidth = 0.1, fill = "grey70", color = "white") +
      labs(
        x = paste0("Adjusted Age (years),  N = ", length(liza_extra_ids)),
        y = NULL,
        title = "Extra IDs in Liza DF"
      ) +
      theme_minimal(base_size = 16) +
      theme(
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5)
      )

ggplot(data = median_distinct %>% filter(participant_id %in% reg_extra_ids),
       aes(x = adjusted_age_years)) +
  geom_histogram(binwidth = 0.1, fill = "grey70", color = "white") +
  labs(
    x = paste0("Adjusted Age (years),  N = ", length(reg_extra_ids)),
    y = NULL,
    title = "Extra IDs in Reg DF"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )