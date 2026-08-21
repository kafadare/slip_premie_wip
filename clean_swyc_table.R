# Save median swyc table, and patient-form expanded tables.
# A lot of the code should be in slip_swyc.Rmd and slip_geocode.Rmd

# ??? do I re-do SWYC calculation to calculate the non-verified domain-specific SWYCs?

library(bigrquery) #package to import SQL tables into R
library(tidyverse) #package suite for data manipulation, reading data into R, and plotting (ggplot2)
library(glue) #package to use expressions with strings
library(ggseg)
require(gridExtra)
library(interactions)
library(boot)
library(parameters)

panel_save <- function(plot, filename, extension, folder, plot_width = 8, plot_height = 10){
  save_name <- paste0(folder, filename, ".", extension)
  ggsave(save_name, plot = plot, dpi = 300, width = plot_width, height = plot_height)
}

`%nin%` = Negate(`%in%`)

bq_auth()
proj_id = "scit605-healthymri-2d73d69f"

#Download SWYC table
query_swyc <- glue("SELECT * FROM `lab.swyc_normative_scores_2026-04-13`")
swyc <- bq_project_query(proj_id, query_swyc) %>% bq_table_download(., page_size=3500)
print(paste0("Arcus lab.SWYC table, dimensions: ", dim(swyc)))

##
##
##

## Add Date Range to SWYC ##
encounter_id_list <- paste0("'", sapply(strsplit(swyc$encounter_id, ","),trimws), "'", collapse = ", ")

#Download encounter table to get datetime
enc_query <- glue("SELECT encounter_id, effective_datetime, appt_time FROM 
              `arcus.encounter`
              WHERE encounter_id IN ({encounter_id_list});")
encounter <- bq_project_query(proj_id, enc_query) %>% bq_table_download(., page_size=3500)
print(paste0("Arcus encounter table, dimensions: ", dim(encounter)))

sum(swyc$encounter_id %nin% encounter$encounter_id) #0

swyc <- merge(swyc, encounter, by = "encounter_id", all.x = TRUE)
sum(is.na(swyc$effective_datetime)) # 0
sum(is.na(swyc$appt_time)) # 3

range(swyc$effective_datetime)

##
##
##

## Clean up duplicate patient-SWYC form pairs ##

swyc$pat_form <- paste0(swyc$pat_id, swyc$form_friendly_name)

dup_df <- swyc %>%
  group_by(pat_form) %>%
  filter(n() > 1)

dup_df_summary <- dup_df %>%
  group_by(pat_form) %>%
  summarise(
    across(c(n_answers, effective_age, age_mo, norm_devo_quotient),
           ~ sd(.x, na.rm = TRUE),
           .names = "sd_{.col}"),
    .groups = "drop"
  )

median(dup_df_summary$sd_n_answers)
median(dup_df_summary$sd_effective_age)
median(dup_df_summary$sd_age_mo)
median(dup_df_summary$sd_norm_devo_quotient)

max(dup_df_summary$sd_n_answers)
max(dup_df_summary$sd_effective_age)
max(dup_df_summary$sd_age_mo)
max(dup_df_summary$sd_norm_devo_quotient)

swyc_clean <- swyc %>%
  group_by(pat_form) %>%
  slice_max(order_by = n_answers, with_ties = TRUE) %>%
  slice_max(order_by = effective_age, with_ties = TRUE) %>%
  distinct(pat_form, .keep_all = TRUE)

# Table of SWYC Form Availability
tab_swyc_forms <- table(swyc_clean$form_friendly_name)
data.frame(form = names(tab_swyc_forms), count = as.vector(tab_swyc_forms))

##
##
##

## Get BW and GA from participants table for all of SWYC ##

pat_id_list <- paste0("'", sapply(strsplit(swyc_clean$pat_id, ","),trimws), "'", collapse = ", ")

#Download Patient Table with GA, BW, and DOB for subjects in the SWYC table
query_pat <- glue("SELECT 
              pat_id,
              sex,
              dob,
              birth_weight_kg,
              gestational_age
              FROM arcus.patient
              WHERE pat_id IN ({pat_id_list});")
pat <- bq_project_query(proj_id, query_pat) %>% bq_table_download(., page_size=3500) %>%
  rename(gestational_age_string = gestational_age)

# Convert different gestational_age_string formats and convert to both 1) number of days AND 2) number of weeks (rounded up/down) AND 3) number of COMPLETED weeks to the best of our knowledge
data <- pat$gestational_age_string

data <- sub("\\s", "", data)
complete_weeks <- as.numeric(sub("^(.{2}).*", "\\1", data))
indices_plus <- grep("\\+", data)
indices_fraction <- grep("/", data)
indices_decimal <- grep("\\.", data)
days <- rep(NA, length(data))
days[!is.na(data)] <- 3 #add three days by default if no other information is encoded for days
days[indices_fraction] <- as.numeric(sub("^.*?(\\d{1}/\\d{1})$", "\\1", data[indices_fraction]) %>% sub("^([0-9]{1}).*", "\\1", .))
days[indices_plus] <- 3 # add three days if gestational age is entered as [week number]+ (like 40+ of 38+)
days[indices_decimal] <- round(as.numeric(sub("^.*?(\\.\\d+)$", "\\1", data[indices_decimal]))*7)

days_total <- complete_weeks*7 + days
weeks_rounded <- round(days_total/7)

pat$ga_complete_weeks <- complete_weeks
pat$ga_days_total <- days_total
pat$ga_weeks_rounded <- weeks_rounded

# Deal with BW/GA QC
range(pat$birth_weight_kg, na.rm = T)
range(pat$ga_days_total/7, na.rm = T)
# GA: remove anyone encoded as > 44*7 ga days or < 22 weeks (there is a single 43.5 weeker, I don't really want to remove this entry)
# BW: remove anyone with birthweight recorded as > 10kg and < 0.3kg
pat$ga_days_total_qc <- ifelse(pat$ga_days_total < 22*7 | pat$ga_days_total >= 44*7, NA, pat$ga_days_total)
sum(!is.na(pat$ga_days_total) & is.na(pat$ga_days_total_qc)) # 0 removed
pat$birth_weight_kg_qc <- ifelse(pat$birth_weight_kg < 0.3 | pat$birth_weight_kg > 10, NA, pat$birth_weight_kg)
sum(!is.na(pat$birth_weight_kg) & is.na(pat$birth_weight_kg_qc)) # 5 removed


# Impute gestational age for entries with birth weight but no gestational age
bw_ga_table <- read.csv("/mnt/arcus/lab/users/kafadare/Fenton_2013_bw_LMS_upload.csv")

bwga_fct_female <- with(bw_ga_table %>% filter(Sex == "Female"), approxfun(x = M, y = as.numeric(wk)))
bwga_fct_male <- with(bw_ga_table %>% filter(Sex == "Male"), approxfun(x = M, y = as.numeric(wk)))

pat <- as.data.frame(pat)

# Predict GA from birthweight (post-QC)
pat[pat$sex == "Female", "ga_pred"] <- bwga_fct_female(pat[pat$sex == "Female","birth_weight_kg_qc"]*1000)
pat[pat$sex == "Male", "ga_pred"] <- bwga_fct_male(pat[pat$sex == "Male","birth_weight_kg_qc"]*1000)

# Examine distribution of "unmatched" birthweights. Appears mostly higher, but a couple on the lower end
hist(pat[is.na(pat$ga_pred) & !is.na(pat$birth_weight_kg_qc), "birth_weight_kg"])
range(pat[is.na(pat$ga_pred) & !is.na(pat$birth_weight_kg_qc), "birth_weight_kg"])

# set predicted ga for out of range birthweights to 41 for high and 22 for low out of range.
pat[is.na(pat$ga_pred) & !is.na(pat$birth_weight_kg_qc) & pat$birth_weight_kg_qc > 3.5, "ga_pred"] <- 41
pat[is.na(pat$ga_pred) & !is.na(pat$birth_weight_kg_qc) & pat$birth_weight_kg_qc < 0.55, "ga_pred"] <- 22
# This below should be 0 after setting floor and ceiling
sum(is.na(pat$ga_pred) & !is.na(pat$birth_weight_kg_qc)) # 0
pat$ga_pred_days_total <- pat$ga_pred*7

# Can get rid of the plot & cor.test below, it's for validation purposes
ggplot(pat %>% filter(!is.na(ga_pred_days_total) & !is.na(ga_days_total_qc)), aes(x = ga_pred_days_total, y = ga_days_total_qc)) +
  geom_point() +
  labs(x = "GA Pred Using Fenton2013 completed weeks w/ Interpolation (days)", y = "Ground Truth GA (days)",
       title = "For Whole Sample with Birthweight") +
  theme_classic()

cor.test(pat$ga_days_total_qc, pat$ga_pred_days_total)


# Create new column for imputed ga in days total, using known GA if available and predicted GA if not
# if predicted GA is not available, replace with median in the sample
pat$ga_imputed_days_total <- ifelse(!is.na(pat$ga_days_total_qc), pat$ga_days_total_qc,
                                   ifelse(!is.na(pat$ga_pred), pat$ga_pred,
                                          median(pat$ga_days_total_qc, na.rm = T)))

sum(is.na(pat$ga_imputed_days_total)) #0 -- should be 0

# Create new column for imputed ga in days total, using known GA if available and predicted GA if not
# if predicted GA is not available, keep NA
pat$ga_imputed_nomedian_days_total <- ifelse(!is.na(pat$ga_days_total_qc), pat$ga_days_total_qc,
                                            ifelse(!is.na(pat$ga_pred_days_total), pat$ga_pred_days_total,
                                                   NA))

sum(is.na(pat$ga_imputed_nomedian_days_total)) #7598
sum(is.na(pat$ga_days_total_qc)) #8392

##
##
##

## Calculate BWP ## (Also using Fenton Table Loaded Above for Imputation) ##

# Function to calculate z score from L, M, S values
zscore_fromLMS <- function(value, L, M, S){
  out <- ((((value/M)^L)-1)/(L*S))
  return(out)
}

# Interpolation for LMS functions
L_fct_female <- with(bw_ga_table %>% filter(Sex == "Female"), approxfun(x = as.numeric(wk), y = L))
L_fct_male <- with(bw_ga_table %>% filter(Sex == "Male"), approxfun(x = as.numeric(wk), y = L))
M_fct_female <- with(bw_ga_table %>% filter(Sex == "Female"), approxfun(x = as.numeric(wk), y = M))
M_fct_male <- with(bw_ga_table %>% filter(Sex == "Male"), approxfun(x = as.numeric(wk), y = M))
S_fct_female <- with(bw_ga_table %>% filter(Sex == "Female"), approxfun(x = as.numeric(wk), y = S))
S_fct_male <- with(bw_ga_table %>% filter(Sex == "Male"), approxfun(x = as.numeric(wk), y = S))

# Calculate birth weight percentile using LMS to z function and defining the correct LMS parameters from Fenton 2013 charts
pat <- pat %>%
  mutate(gestational_age_restrict = case_when(
    ga_days_total_qc/7 >= 42 ~ 42*7,
    TRUE ~ ga_days_total_qc) ) %>%
  mutate(bwz_fen = case_when(
    sex == "Female" ~ zscore_fromLMS(birth_weight_kg*1000, 
                                     L_fct_female(gestational_age_restrict/7),
                                     M_fct_female(gestational_age_restrict/7),
                                     S_fct_female(gestational_age_restrict/7)),
    sex == "Male" ~ zscore_fromLMS(birth_weight_kg*1000, 
                                   L_fct_male(gestational_age_restrict/7),
                                   M_fct_male(gestational_age_restrict/7),
                                   S_fct_male(gestational_age_restrict/7)),
    TRUE ~ NA)) %>%
  mutate(bwp_fen = pnorm(bwz_fen) ) %>%
  select(-c("gestational_age_restrict"))

sum(!is.na(pat$birth_weight_kg) & !is.na(pat$ga_days_total_qc))

sum(!is.na(pat$bwz_fen))
sum(!is.na(pat$bwp_fen))

#pat[!is.na(pat$birth_weight_kg) & !is.na(pat$ga_days_total_qc) & is.na(pat$bwp_fen), c("ga_days_total_qc", "birth_weight_kg")]

#hist(pat$bwp_fen)

# Save birth factors DF including imputed GA and pre-QC bw/ga values
write.csv(pat, "/mnt/arcus/lab/users/kafadare/swyc_premie_analysis_results/processed_tables/birth_factors.csv")

## Add Birth Factors (GA, BW, and BWP) to the SWYC Table

swyc_clean <- merge(swyc_clean, pat %>% 
                select(pat_id, dob, birth_weight_kg_qc, ga_days_total_qc, bwz_fen, bwp_fen,
                       ga_complete_weeks) %>%
                rename(birth_weight_kg = birth_weight_kg_qc,
                       gestational_age = ga_days_total_qc),
              by = "pat_id", all.x = TRUE)

##
##
##

## Get ADI values and add to SWYC. Median ADI up to each SWYC entry age/date

env <- readRDS("/mnt/arcus/lab/shared/env_data.rds") %>% filter(pat_id %in% swyc_clean$pat_id)

# Add DOB to env table to calculate age in days for start_date and end_data of env variables
env <- merge(env, pat %>% select(pat_id, dob), by = "pat_id", all.x = TRUE)

sum(is.na(env$dob)) #0

env <- env %>% filter(!is.na(dob))

#Calculate unadjusted age in days using dob and end date for a particular entry
env$unadjusted_age_in_days_at_end_date <- as.numeric(difftime(env$end_date, as.Date(env$dob), units = "days"))
env$unadjusted_age_in_days_at_start_date <- as.numeric(difftime(env$start_date, as.Date(env$dob), units = "days"))

env_vars<- c("ADI3_EDU", "ADI3_ECON", "ADI3_FINS", "ADI")

swyc_env <- merge(swyc_clean , env, by = "pat_id", all.x = TRUE)

swyc_env$swyc_adi_start_age_diff <- swyc_env$effective_age - swyc_env$unadjusted_age_in_days_at_start_date
swyc_env$swyc_adi_end_age_diff <- swyc_env$effective_age - swyc_env$unadjusted_age_in_days_at_end_date

hist(swyc_env$swyc_adi_end_age_diff)
hist(swyc_env$swyc_adi_start_age_diff)

swyc_env_median <- swyc_env %>%
  group_by(pat_form) %>%
  mutate(env_limit_row = min(swyc_adi_start_age_diff[swyc_adi_start_age_diff >= 0])) %>%
  filter(swyc_adi_start_age_diff >= env_limit_row) %>%
  summarise(
    n_rows = n(),
    across(any_of(names(swyc_clean)),
           ~ first(.x),
           .names = "{.col}"),
    across(any_of(env_vars),
           list(
             median = ~ median(.x, na.rm = TRUE),
             n_val = ~ sum(!is.na(.x))
           ),
           .names = "{col}.{fn}")
  )

swyc_clean <- merge(swyc_clean, swyc_env_median %>% 
                 select(pat_form, ADI.median, ADI3_ECON.median, 
                        ADI3_FINS.median, ADI3_EDU.median,
                        ADI.n_val, ADI3_ECON.n_val, 
                        ADI3_FINS.n_val, ADI3_EDU.n_val),
               by = "pat_form", all.x = TRUE)

sum(!is.na(swyc_clean$ADI.median))

##
##
##

## Calculate GA-Adjusted DQ for each swyc patient-form entry ##

median_ga <- median(swyc_clean %>% distinct(pat_id, .keep_all = TRUE) %>% select(gestational_age) %>% pull(.), 
                    na.rm = TRUE)
median_ga_term <- median(swyc_clean %>% distinct(pat_id, .keep_all = TRUE) %>% 
                           filter(gestational_age >= 7*37) %>% select(gestational_age) %>% pull(.), 
                         na.rm = TRUE)
swyc_clean$ga_diff_termmedian <- median_ga_term - swyc_clean$gestational_age

swyc_clean$ga_adj_effective_age <- swyc_clean$effective_age - swyc_clean$ga_diff_termmedian
swyc_clean$ga_adj_age_mo <- round(swyc_clean$ga_adj_effective_age/30.437)
#plot(swyc_clean$age_mo, swyc_clean$ga_adj_age_mo)
swyc_clean$ga_adj_norm_devo_quotient <- swyc_clean$norm_devo_age/swyc_clean$ga_adj_age_mo
swyc_clean$ga_adj_norm_devo_quotient_round <- round(swyc_clean$ga_adj_norm_devo_quotient, 5)
#plot(swyc_clean$norm_devo_quotient , swyc_clean$ga_adj_norm_devo_quotient)

ggplot(swyc_clean, aes(x = norm_devo_quotient, y = ga_adj_norm_devo_quotient_round, color = ga_complete_weeks)) +
  geom_point() +
  labs(x = "OG DQ",
       y = "GA-ADJ DQ",
       title = "All GA") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#B0B0B0")

cor.test(swyc_clean$ga_adj_norm_devo_quotient, swyc_clean$gestational_age)
cor.test(swyc_clean$norm_devo_quotient, swyc_clean$gestational_age)

sum((swyc_clean$gestational_age < median_ga_term) & (swyc_clean$ga_adj_norm_devo_quotient_round < swyc_clean$norm_devo_quotient) & !is.na(swyc_clean$gestational_age)) # should be 0

# swyc_clean[(swyc_clean$gestational_age < median_ga_term) & (swyc_clean$ga_adj_norm_devo_quotient_round < swyc_clean$norm_devo_quotient) & !is.na(swyc_clean$gestational_age),
#            c("gestational_age", "ga_complete_weeks", "ga_adj_norm_devo_quotient", "norm_devo_quotient", "ga_adj_age_mo", "age_mo", "norm_devo_age")]


write.csv(swyc_clean, "/mnt/arcus/lab/users/kafadare/swyc_premie_analysis_results/processed_tables/swyc_allforms.csv")

## 
##
##

## Get median SWYC DQ for each participant, keep ADI for the latest SWYC DQ ##

swyc_avg <- swyc_clean %>%
  group_by(pat_id) %>%
  summarise(n_forms = n(),
            dob = first(dob),
            sex = first(sex_abbr),
            birth_weight_kg = first(birth_weight_kg),
            bwp_fen = first(bwp_fen),
            bwz_fen = first(bwz_fen),
            ga_complete_weeks = first(ga_complete_weeks),
            gestational_age   = first(gestational_age),
            all_2s_sum = sum(all_2s == TRUE, na.rm = TRUE),
            all_0s_sum = sum(all_0s == TRUE, na.rm = TRUE),
            final_swyc_age = max(effective_age, na.rm = TRUE),
            final_swyc_age_ga_adj = max(ga_adj_effective_age),
            median_dq = median(norm_devo_quotient, na.rm = TRUE), 
            ga_adj_median_dq = median(ga_adj_norm_devo_quotient_round, na.rm = TRUE),
            final_swyc_date = effective_datetime[which.max(effective_age)],
            final_ADI = ADI.median[which.max(effective_age)],
            final_ADI3_FINS = ADI3_FINS.median[which.max(effective_age)],
            final_ADI3_ECON = ADI3_ECON.median[which.max(effective_age)],
            final_ADI3_EDU = ADI3_EDU.median[which.max(effective_age)],
            .groups = "drop")

hist(swyc_avg$median_dq)
hist(swyc_avg$ga_adj_median_dq)
hist(swyc_avg$final_swyc_age/30.437)
hist(swyc_avg$final_ADI)

write.csv(swyc_avg, "/mnt/arcus/lab/users/kafadare/swyc_premie_analysis_results/processed_tables/swyc_avg.csv")


##
##
##

## Add SWYC All Forms and Avg Separately to Centile DFs and Save ##

# Add swyc median values to df.zscore file and save
# Save several versions with different date-limits between SWYC and Centile Scores
### Find SWYCs after a scan and take the median among those (using centiles all)
# For all available forms, save under different name (because it will duplicate centile rows)
# Again, save several versions depending on time between scan and swyc

cent_all_path <- "/mnt/arcus/lab/users/kafadare/slip_centiles_processed_tables/centiles_all_median_th_mpr_2026-04-24.csv"
cent_all <- read.csv(cent_all_path)

cent_distinct_path <- "/mnt/arcus/lab/users/kafadare/slip_centiles_processed_tables/df_zscore_distinct.csv"
cent_distinct <- read.csv(cent_distinct_path)
load("/mnt/arcus/lab/users/kafadare/slip_centiles_processed_tables/lookup_stats_distinct.RData")
cent_distinct <- cent_distinct %>%
  mutate(age_at_scan_raw = map2_dbl(age_at_scan, "age_at_scan", ~ .x*sd_lookup_distinct[[.y]] + mean_lookup_distinct[[.y]]))

# cent_distinct <- cent_all %>% 
#   group_by(participant_id) %>%
#   slice_min(age_days_adj, n = 1, with_ties = FALSE) %>% 
#   ungroup()
# 
# # Add QC Variable Name
# if(!("median_euler_mean" %in% names(cent_distinct))){
#   cent_distinct$qc_var <- cent_distinct$euler_mean
# } else if("median_euler_mean" %in% names(cent_distinct)){
#   cent_distinct$qc_var <- cent_distinct$median_euler_mean
# }
# 
# # Remove BW Outliers
# mu_bw <- mean(cent_distinct$birth_weight_kg, na.rm = TRUE)
# sd_bw <- sd(cent_distinct$birth_weight_kg, na.rm = TRUE)
# 
# rm_dist = 5
# 
# nrow(cent_distinct %>% filter(
#   birth_weight_kg < mu_bw-sd_bw*rm_dist | birth_weight_kg > mu_bw+sd_bw*rm_dist) )
# 
# cent_distinct <- cent_distinct %>%
#   mutate(birth_weight_kg_inc_out = birth_weight_kg,
#          birth_weight_kg = case_when(
#            (birth_weight_kg < mu_bw-sd_bw*rm_dist | birth_weight_kg > mu_bw+sd_bw*rm_dist) ~ NA,
#            TRUE ~ birth_weight_kg),
#          bwp_fen_inc_out = bwp_fen,
#          bwp_fen = case_when(
#            is.na(birth_weight_kg) ~ NA,
#            TRUE ~ bwp_fen))
# 
# sum(!is.na(cent_distinct$birth_weight_kg))
# range(cent_distinct$birth_weight_kg, na.rm = TRUE)
# sum(!is.na(cent_distinct$birth_weight_kg_outrm))
# range(cent_distinct$birth_weight_kg_outrm, na.rm = TRUE)


# Calculate CMD


swyc_allforms_cent <- merge(cent_distinct, swyc_clean %>%
                                 select(pat_id, form_friendly_name, pat_form,
                                        norm_devo_quotient,
                                        ga_adj_norm_devo_quotient_round, effective_age,
                                        ga_adj_effective_age, sum_answers,
                                        all_0s, all_2s, effective_datetime) %>%
                              rename(ga_adj_norm_devo_quotient = ga_adj_norm_devo_quotient_round),
                               by = "pat_id", all.x = FALSE, all.y = FALSE)

# Plot Scan-SWYC timing for All Available SWYC Forms

(swyc_allforms_cent_timing_plot <- ggplot(swyc_allforms_cent,
                                              aes(y = reorder(pat_id, age_at_scan_raw))) +
    geom_point(aes(x = effective_age, color = "SWYC"), size = 0.5) +
    geom_point(aes(x = age_at_scan_raw, color = "Scan"), size = 0.5) +
    scale_x_continuous(labels = \(x) round(x/30.4)) +
    scale_color_manual( values = c("SWYC" = "blue",
                                   "Scan" = "darkgray")
    ) +
    labs(x = "Age (months)",
         y = "Patients",
         title = "All SWYC Forms for Patients with Scans (One Scan per Patient)") +
    theme_classic(base_family = "Helvetica", base_size = 14) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "bottom")
)

swyc_allforms_cent$scan_swyc_age_diff_days <- swyc_allforms_cent$age_at_scan_raw - swyc_allforms_cent$effective_age
ggplot(swyc_allforms_cent, aes(x = scan_swyc_age_diff_days)) + 
  geom_histogram(bins = 100) + 
  labs(x = "", y = "Count", 
       title = "Days Between Scan and SWYC Form") +
  theme_classic()
swyc_allforms_cent$scan_before_swyc <- swyc_allforms_cent$age_at_scan_raw < swyc_allforms_cent$effective_age
sum(swyc_allforms_cent$scan_before_swyc) # 872
length(unique(swyc_allforms_cent[swyc_allforms_cent$scan_before_swyc == T, "pat_id"])) #421

# if limited scan within two years of swyc
sum(swyc_allforms_cent$scan_swyc_age_diff_days <= (365.25*2)) # 1451 rows

length(unique(swyc_allforms_cent[swyc_allforms_cent$scan_swyc_age_diff_days <= (365.25*2), "pat_id"])) #591

swyc_avg_cent <- merge(cent_distinct, swyc_avg %>%
                         select(pat_id, all_0s_sum, all_2s_sum, 
                                median_dq, ga_adj_median_dq, 
                                final_swyc_age, final_swyc_age_ga_adj,
                                final_swyc_date),
                       by = "pat_id", all.x = FALSE, all.y = FALSE)

swyc_avg_cent$scan_swyc_age_diff_days <- swyc_avg_cent$age_at_scan_raw - swyc_avg_cent$final_swyc_age

swyc_avg_cent <- swyc_avg_cent %>%
  mutate(scan_swyc_twoyears = abs(scan_swyc_age_diff_days) <= (365.25*2))

sum(swyc_avg_cent$scan_swyc_twoyears) # 553

(swyc_avg_cent_overlap_plot <- ggplot(swyc_avg_cent, 
                                      aes(y = reorder(pat_id, age_at_scan_raw))) +
    geom_point(aes(x = final_swyc_age, color = "SWYC Final Age"), size = 0.5) +
    geom_point(aes(x = age_at_scan_raw, color = "Scan"), size = 0.5) +
    scale_x_continuous(labels = \(x) round(x/30.4)) +
    scale_color_manual( values = c("SWYC Final Age" = "blue",
                                   "Scan" = "darkgray")
    ) +
    labs(x = "Age (months)",
         y = "Patients",
         title = "Final SWYC Age With Scan (One Scan per Patient)") +
    theme_classic(base_family = "Helvetica", base_size = 14) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "bottom")
)

(swyc_avg_cent_overlap_twoyears_plot <- ggplot(swyc_avg_cent %>% 
                                        filter(scan_swyc_twoyears), 
                                      aes(y = reorder(pat_id, age_at_scan_raw))) +
    geom_point(aes(x = final_swyc_age, color = "SWYC Final Age"), size = 0.5) +
    geom_point(aes(x = age_at_scan_raw, color = "Scan"), size = 0.5) +
    scale_x_continuous(labels = \(x) round(x/30.4)) +
    scale_color_manual( values = c("SWYC Final Age" = "blue",
                                   "Scan" = "darkgray")
    ) +
    labs(x = "Age (months)",
         y = "Patients",
         title = "Final SWYC Age Within Two Years of Scan (One Scan Per Patient)") +
    theme_classic(base_family = "Helvetica", base_size = 14) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "bottom")
)

ggplot(swyc_avg_cent, aes(x = scan_swyc_age_diff_days)) + 
  geom_histogram(bins = 100) + 
  labs(x = "", y = "Count", 
       title = "Days Between Scan and Final Age Used for SWYC Median") +
  theme_classic()


write.csv(swyc_avg_cent, "/mnt/arcus/lab/users/kafadare/swyc_premie_analysis_results/processed_tables/swyc_avg_cent.csv")

## Take SWYC DQ Median Within SWYC forms that are <= 2 years from Scan
swyc_avg_forcent <- swyc_allforms_cent %>%
  filter(scan_swyc_age_diff_days <= (365.25*2)) %>%
  group_by(pat_id) %>%
  summarise(n_forms = n(),
            all_2s_sum = sum(all_2s == TRUE, na.rm = TRUE),
            all_0s_sum = sum(all_0s == TRUE, na.rm = TRUE),
            final_swyc_age = max(effective_age, na.rm = TRUE),
            final_swyc_age_ga_adj = max(ga_adj_effective_age),
            median_dq = median(norm_devo_quotient, na.rm = TRUE), 
            ga_adj_median_dq = median(ga_adj_norm_devo_quotient, na.rm = TRUE),
            final_swyc_date = effective_datetime[which.max(effective_age)],
            .groups = "drop")

write.csv(swyc_avg_forcent, "/mnt/arcus/lab/users/kafadare/swyc_premie_analysis_results/processed_tables/swyc_avg_twoyrswithinscan.csv")

swyc_avg_cent_twoyrs <- merge(cent_distinct, swyc_avg_forcent %>%
                         select(pat_id, all_0s_sum, all_2s_sum, 
                                median_dq, ga_adj_median_dq, 
                                final_swyc_age, final_swyc_age_ga_adj,
                                final_swyc_date),
                       by = "pat_id", all.x = FALSE, all.y = FALSE)

(swyc_avg_cent_overlap_plot <- ggplot(swyc_avg_cent_twoyrs, 
                                      aes(y = reorder(pat_id, age_at_scan_raw))) +
    geom_point(aes(x = final_swyc_age, color = "SWYC Final Age"), size = 0.5) +
    geom_point(aes(x = age_at_scan_raw, color = "Scan"), size = 0.5) +
    scale_x_continuous(labels = \(x) round(x/30.4)) +
    scale_color_manual( values = c("SWYC Final Age" = "blue",
                                   "Scan" = "darkgray")
    ) +
    labs(x = "Age (months)",
         y = "Patients",
         title = glue("Final SWYC Age With Scan, \n SWYC Avg Within Scans Two Yrs Apart \n (One Scan per Patient)")) +
    theme_classic(base_family = "Helvetica", base_size = 14) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "bottom")
)

# Save cent_distinct with SWYC Data

write.csv(swyc_allforms_cent,  "/mnt/arcus/lab/users/kafadare/swyc_premie_analysis_results/processed_tables/swyc_allforms_cent.csv")
write.csv(swyc_avg_cent_twoyrs,  "/mnt/arcus/lab/users/kafadare/swyc_premie_analysis_results/processed_tables/swyc_avgtwoyrs_cent.csv")

