#Script for looking at data in the slip 25-11 release
require(tidyverse)
require(growthstandards)

code_path <- "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/"
source(paste0(code_path,"slip_data_utils.R"))

vars_save_path <-  "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/data_processed_csv/"
median_distinct <- read.csv(paste0(vars_save_path,"median_distinct_clean_2026-01-07.csv")) %>% select(-("X"))

main_folder <- "/mnt/isilon/bgdlab_processing/releases_clinical/2025_11_release/"
participants_path <- paste0(main_folder,"BIDS/participants.tsv")
participants <- read.table(participants_path,
                           sep = "\t",
                           header = TRUE,
                           stringsAsFactors = FALSE,
                           quote = "",
                           comment.char = "")
length(unique(participants$participant_id[!(participants$participant_id %in% median_distinct$participant_id)]))

#Data Cleanup Steps
##Convert birth variables to numerical where "n/a" is set to NA
participants <- participants %>%
  mutate(
    gestational_age = as.numeric(na_if(gestational_age, "n/a")),
    birth_weight_kg = as.numeric(na_if(birth_weight_kg, "n/a")),
    birth_length_cm = as.numeric(na_if(birth_length_cm, "n/a"))
    )
##See if there are discrepancies between session ids for the same participant for birth variables
birth_var_discrepancy_table <- participants %>%
  group_by(participant_id) %>%
  summarise(n_var_bl = n_distinct(birth_length_cm),
            n_var_bw = n_distinct(birth_weight_kg),
            n_var_ga = n_distinct(gestational_age))%>%
  filter(n_var_bl > 1 | n_var_bw > 1 | n_var_ga > 1)
participants %>% filter(participant_id %in% birth_var_discrepancy_table$participant_id) %>%
  select(c("participant_id", "session_id", "birth_length_cm", "birth_weight_kg", "gestational_age")) %>%
  arrange(participant_id)
##keep the available value, or mean of available values for all sessions for birth variables
participants <- participants %>%
  group_by(participant_id) %>%
  mutate(
    gestational_age = if (all(is.na(gestational_age))) NA_real_
    else mean(gestational_age, na.rm = TRUE),
    birth_weight_kg = if (all(is.na(birth_weight_kg))) NA_real_
    else mean(birth_weight_kg, na.rm = TRUE),
    birth_length_cm = if (all(is.na(birth_length_cm))) NA_real_
    else mean(birth_length_cm, na.rm = TRUE)
    ) %>%
  ungroup()

#Birthweight and GA QC
##Remove really high birthweight values to get a meaningful mean/sd value
bw_absurd <- sum(participants$birth_weight_kg  > 1000, na.rm = TRUE)
print(paste0("Number of data point set to NA for recorded BW > 1000kg: ", bw_absurd))
participants$birth_weight_kg[participants$birth_weight_kg>1000] <- NA
bw_mean <- mean(median_distinct$birth_weight_kg, na.rm = TRUE)
bw_sd <- sd(median_distinct$birth_weight_kg, na.rm = TRUE)
#apply QC function
participants <- slip_bwga_qc(cent_df = participants, bw_mean = bw_mean, bw_sd = bw_sd,  bw_sd_limit = 4)

#Read in QC File
qc_path <- paste0(main_folder, "derivatives/recon-all-clinical/fs7.4.1-clinical-reconall-qc.csv")
qc_metrics <- read.csv(qc_path)
qc_metrics_filtered <- qc_metrics %>% filter(euler_mean > -60 &
                                    qc_general_white_matter > 0.65 &
                                    qc_general_grey_matter > 0.65 &
                                    qc_general_csf > 0.65 &
                                    qc_cerebellum > 0.65 &
                                    qc_brainstem > 0.65 &
                                    qc_thalamus > 0.65 &
                                    qc_putamen_pallidum > 0.65 &
                                    qc_hippocampus_amygdala > 0.65)

participants_postqc <- participants %>% filter(participants$participant_id %in% qc_metrics_filtered$subject_id)
length(unique(participants_postqc$participant_id[!(participants_postqc$participant_id %in% median_distinct$participant_id)]))

#See distribution of variables in the new participants df (including the previously included IDs)
participants_postqc_distinct <- participants_postqc %>% distinct(participant_id, .keep_all = TRUE)
N_participants_full <- dim(participants_postqc_distinct)[1]

#Age
participants_postqc_distinct$age_years <- participants_postqc_distinct$age_at_scan/365
participants_postqc_distinct$adjusted_age_years <- participants_postqc_distinct$adjusted_age_in_days/365
##Cut < 21 years for visual comparison & continuity
participants_postqc_distinct_u21 <- participants_postqc_distinct %>% filter(adjusted_age_years < 21)
N_participants_u21 <- dim(participants_postqc_distinct_u21)[1]

##plot full age range
ggplot(participants_postqc_distinct, aes(x = adjusted_age_years)) +
  geom_histogram(binwidth = 0.1, fill = "grey70", color = "white") +
  labs(
    x = paste0("Adjusted Age (years),  N = ", N_participants_full),
    y = NULL,
    title = "25-11 all participants postQC"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )
##plot u21 age range
ggplot(participants_postqc_distinct_u21, aes(x = adjusted_age_years)) +
  geom_histogram(binwidth = 0.1, fill = "grey70", color = "white") +
  labs(
    x = paste0("Adjusted Age (years),  N = ", N_participants_u21),
    y = NULL,
    title = "25-11 all participants postQC"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

#Availability of Birth Variables Plots

##get numbers for birth variables availability
N_all <- sum(!is.na(participants_postqc_distinct_u21$adjusted_age_years))
N_ga <- sum(!is.na(participants_postqc_distinct_u21$gestational_age))
N_bw <- sum(!is.na(participants_postqc_distinct_u21$birth_weight_kg))
N_ga_only <- sum(!is.na(participants_postqc_distinct_u21$gestational_age) & is.na(participants_postqc_distinct_u21$birth_weight_kg))
N_bw_only <- sum(!is.na(participants_postqc_distinct_u21$birth_weight_kg) & is.na(participants_postqc_distinct_u21$gestational_age))
N_gabw <- sum(!is.na(participants_postqc_distinct_u21$gestational_age) & !is.na(participants_postqc_distinct_u21$birth_weight_kg))
N_none <- sum(is.na(participants_postqc_distinct_u21$gestational_age) & is.na(participants_postqc_distinct_u21$birth_weight_kg))

##create a data table with categories based on birth variables availability
df_age <- participants_postqc_distinct_u21 %>%
  mutate(
    ga_known  = !is.na(gestational_age),
    bw_known  = !is.na(birth_weight_kg),
    group = case_when(
      ga_known & bw_known ~ "GA & BW Known",
      ga_known & !bw_known ~ "Only GA Known",
      bw_known & !ga_known ~ "Only BW Known",
      TRUE ~ "None"
    )) %>%
  mutate(group = factor(group, levels = c(
    "Only GA Known",
    "Only BW Known",
    "GA & BW Known",
    "None"
  )))

##plot birth variable availability across age range
ggplot(df_age, aes(x = adjusted_age_years, fill = group)) +
    geom_histogram(binwidth = 0.1, color = "white") +
    scale_fill_manual(
      values = c(
        "Only GA Known" = "#9EC3D2",
        "Only BW Known" = "#D6B8A8",
        "GA & BW Known" = "#6A5A78",
        "None" = "#E0E0E0"
      ),
      labels = c(
        paste0("Only GA Known, N=", N_ga_only),
        paste0("Only BW Known, N=", N_bw_only),
        paste0("GA & BW Known, N=", N_gabw),
        paste0("Neither Known, N=", N_none)
      ),name = NULL) +
    labs(
      x = paste0("Adjusted Age (years),  N = ", N_all),
      y = NULL
    ) +
    theme_minimal(base_size = 16) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = element_text(size = 14),
      legend.key.width = unit(1.5, "cm"),
      legend.box.just = "center",
      panel.grid = element_blank()
    ) +
    guides(fill = guide_legend(nrow = 2))

#Birth Variable Distributions
##Assign preterm categories
participants_postqc_distinct_u21$preterm <- as.factor(cut(participants_postqc_distinct_u21$gestational_age, breaks = c(-Inf, 31.9, 36.9, Inf), labels = c("VPM", "LPM", "Term"), include.lowest = TRUE))
##Calculate birth weight percentile
participants_postqc_distinct_u21 <- participants_postqc_distinct_u21 %>%
  mutate(gestAge_days = gestational_age*7+3) %>%
  filter(sex != "Unknown")
participants_postqc_distinct_u21$bweight_percentile <- igb_wtkg2centile(participants_postqc_distinct_u21$gestAge_days, participants_postqc_distinct_u21$birth_weight_kg, sex = participants_postqc_distinct_u21$sex)/100

##Gestational Age Distribution
ga_ptb <- participants_postqc_distinct_u21 %>%
  filter(!is.na(gestational_age)) %>%
  mutate(group = ifelse(gestational_age < 37, "Preterm", "Term"))

pct_preterm <- round(mean(ga_ptb$group == "Preterm") * 100, 1)

(ga_preterm_dist_plot <- ggplot(ga_ptb, aes(x = gestational_age, fill = group)) +
    geom_histogram(binwidth = 0.5, color = "white") +
    scale_fill_manual(values = c("Preterm" = "#F4A475", "Term" = "#A7D9A0")) +
    labs(
      x = paste0("Gestational Age(weeks),  % Preterm(<37wks) = ", pct_preterm),
      y = NULL
    ) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    ))

participants_postqc_distinct_u21 <- participants_postqc_distinct_u21 %>%
  mutate(preterm = factor(preterm, levels = c("VPM", "LPM", "Term")))

(ga_preterm_age_plot <- ggplot(participants_postqc_distinct_u21 %>% filter(!is.na(preterm)), aes(x = adjusted_age_years, fill = preterm)) +
    geom_histogram(position = "fill", binwidth = 0.25, color = "white") +
    scale_fill_manual(
      values = c(
        "VPM"  = "#E8A0A8",
        "LPM"  = "#A8C4E8",
        "Term" = "#A7D9A0"
      )
    ) +
    labs(
      x = "Adjusted Age (years)",
      y = "Proportion",
      fill = "Preterm Status"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal"
    ))

##Birthweight Distribution
participants_postqc_distinct_u21$bw_category <- as.factor(cut(participants_postqc_distinct_u21$birth_weight_kg, breaks = c(-Inf, 1.5, 2.5, Inf), labels = c("VLBW", "LBW", "Normal BW"), include.lowest = TRUE))

bw_low <- participants_postqc_distinct_u21 %>%
  filter(!is.na(birth_weight_kg)) %>%
  mutate(group = ifelse(birth_weight_kg < 2.5, "Low BW", "Normal BW"))

pct_lowbw <- round(mean(bw_low$group == "Low BW") * 100, 1)

(bw_low_dist_plot <- ggplot(bw_low, aes(x = birth_weight_kg, fill = group)) +
    geom_histogram(
      binwidth = 0.1,
      boundary = 0,
      color = "white"
    ) +
    scale_fill_manual(values = c("Low BW" = "#C76E6A", "Normal BW" = "#7FAF8A")) +
    labs(
      x = paste0("Birthweight(kg),  % Low(<2.5kg) = ", pct_lowbw),
      y = NULL
    ) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    ))

(bw_low_age_plot <- ggplot(participants_postqc_distinct_u21 %>% filter(!is.na(bw_category)), aes(x = adjusted_age_years, fill = bw_category)) +
    geom_histogram(position = "fill", binwidth = 0.25, color = "white") +
    scale_fill_manual(
      values = c(
        "VLBW"  = "#B9836A",
        "LBW"  = "#6A8F9F",
        "Normal BW" = "#7FAF8A"
      )
    ) +
    labs(
      x = "Adjusted Age (years)",
      y = "Proportion",
      fill = "Birthweight"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal"
    ))

participants_postqc_distinct_u21$age_group <- cut(participants_postqc_distinct_u21$adjusted_age_years, breaks = 80, include.lowest = TRUE, labels = FALSE)

participants_postqc_distinct_u21 <- participants_postqc_distinct_u21 %>%
  group_by(age_group) %>%
  arrange(birth_weight_kg, .by_group = TRUE) %>%
  ungroup()

#for x axis labeling
age_labels <- c(0, 5, 10, 15, 20)
group_width <- (max(participants_postqc_distinct_u21$adjusted_age_years) - min(participants_postqc_distinct_u21$adjusted_age_years)) / 80
x_breaks <- floor((age_labels - min(participants_postqc_distinct_u21$adjusted_age_years)) / group_width) + 1
#fix position of 0
x_breaks[1] <- 1

(bw_cont_age_plot <- ggplot(participants_postqc_distinct_u21 %>% filter(!is.na(birth_weight_kg)), aes(x = age_group, y = 1, fill = birth_weight_kg)) + geom_bar(stat = "identity", width = 0.8) +
    scale_fill_gradient2(low = "#F2EAD3", mid = "#E6E2C8", high = "#7FAF8A", limits = range(participants_postqc_distinct_u21$birth_weight_kg, na.rm = TRUE)) +
    labs(
      x = "Adjusted Age (years)",
      y = NULL,
      fill = "Birthweight"
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = age_labels
    ) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(margin = margin(r = 10))
    ) +
    guides(fill = guide_colorbar(
      barwidth = 10,
      barheight = 1,
      title.position = "left",
      title.vjust = 1
    )))

##Birthweight Percentile Distribution
participants_postqc_distinct_u21 <- participants_postqc_distinct_u21 %>%
  group_by(age_group) %>%
  arrange(bweight_percentile, .by_group = TRUE) %>%
  ungroup()

#for x axis labeling
age_labels <- c(0, 5, 10, 15, 20)
group_width <- (max(participants_postqc_distinct_u21$adjusted_age_years) - min(participants_postqc_distinct_u21$adjusted_age_years)) / 80
x_breaks <- floor((age_labels - min(participants_postqc_distinct_u21$adjusted_age_years)) / group_width) + 1
#fix position of 0
x_breaks[1] <- 1

(bwp_cont_age_plot <- ggplot(participants_postqc_distinct_u21 %>% filter(!is.na(bweight_percentile)), aes(x = age_group, y = 1, fill = bweight_percentile)) + geom_bar(stat = "identity", width = 0.8) +
    scale_fill_gradient2(low = "#F2EAD3", mid = "#E6E2C8", high = "#7FAF8A", limits = c(0,1)) +
    labs(
      x = "Adjusted Age (years)",
      y = NULL,
      fill = "Birthweight Percentile"
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = age_labels
    ) +
    theme_minimal(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(margin = margin(r = 10))
    ) +
    guides(fill = guide_colorbar(
      barwidth = 10,
      barheight = 1,
      title.position = "left",
      title.vjust = 1
    )))

#Data paths for record keeping
#2. SynthSeg+ volumes: `derivatives/recon-all-clinical/fs7.4.1-clinical-reconall-ss.csv`
#4. glasser cortical metrics: `derivatives/recon-all-clinical/fs7.4.1-clinical-reconall-glasser.csv`
#5. global volumes: `derivatives/recon-all-clinical/fs7.4.1-clinical-reconall-global.csv`
#6. aparc cortical metrics: `derivatives/recon-all-clinical/fs7.4.1-clinical-reconall-aparc.csv`
#7. aseg subcortical volumes: `derivatives/recon-all-clinical/fs7.4.1-clinical-reconall-aseg.csv`

