# compare old bwp (using intergrowth) to new bwp (using Fenton L, M, S tables

library(tidyr)
library(dplyr)
library(growthstandards)
library(ggplot2)
library(ggcorrplot)

folder <-  "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11_new_models/"

fenton_table_path <- "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/data_files/Fenton_2013_bw_LMS_upload.csv"

df_vars_select <- c("participant_id", "sex", 
                              "ga_days_total_qc", "ga_imputed_nomedian_days_total", "ga_imputed_days_total",
                              "birth_weight_kg_qc", "birth_length_cm")

df <- read.csv(paste0(folder,"/participants_ga_imputed.csv")) %>%
  select(all_of(df_vars_select)) %>%
  rename(gestational_age = ga_days_total_qc,
         ga_bw_imputed = ga_imputed_nomedian_days_total,
         ga_full_imputed = ga_imputed_days_total,
         birth_weight_kg = birth_weight_kg_qc)

fenton_lms_table <- read.csv(fenton_table_path)
df <- df %>% filter(sex != "Unknown")

#Calculate birth weight percentile using growthstandards function
df <- df %>%
  mutate(sex_recode= recode(sex, "M" = "Male", "F" = "Female")) %>%
  mutate(bwp_ig = igb_wtkg2centile(gestational_age, birth_weight_kg, sex = sex_recode)/100,
         bwz_ig = igb_wtkg2zscore(gestational_age, birth_weight_kg, sex = sex_recode))

# Function to calculate z score from L, M, S values
zscore_fromLMS <- function(value, L, M, S){
  out <- ((((value/M)^L)-1)/(L*S))
  return(out)
}

# Interpolation for LMS functions
L_fct_female <- with(fenton_lms_table %>% filter(Sex == "Female"), approxfun(x = as.numeric(wk), y = L))
L_fct_male <- with(fenton_lms_table %>% filter(Sex == "Male"), approxfun(x = as.numeric(wk), y = L))
M_fct_female <- with(fenton_lms_table %>% filter(Sex == "Female"), approxfun(x = as.numeric(wk), y = M))
M_fct_male <- with(fenton_lms_table %>% filter(Sex == "Male"), approxfun(x = as.numeric(wk), y = M))
S_fct_female <- with(fenton_lms_table %>% filter(Sex == "Female"), approxfun(x = as.numeric(wk), y = S))
S_fct_male <- with(fenton_lms_table %>% filter(Sex == "Male"), approxfun(x = as.numeric(wk), y = S))

# Calculate birth weight percentile using LMS to z function and defining the correct LMS parameters from Fenton 2013 charts
df <- df %>%
  mutate(gestational_age_restrict = case_when(
    gestational_age/7 >= 42 ~ 42*7,
    TRUE ~ gestational_age) ) %>%
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

sum(!is.na(df$birth_weight_kg) & !is.na(df$gestational_age))

sum(!is.na(df$bwz_fen))
sum(!is.na(df$bwp_fen))

df[!is.na(df$birth_weight_kg) & !is.na(df$gestational_age) & is.na(df$bwp_fen), c("gestational_age", "birth_weight_kg", "bwp_ig")]

hist(df$bwp_fen)

cor.test(df$bwp_ig, df$bwp_fen)

ggplot(df, aes(x = bwp_ig, y = bwp_fen)) +
  geom_point(color = "gray", alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  theme_classic() +
  labs(
    x = "BW Percentile Intergrowth",
    y = "BW Percentile Fenton 2013, Linearly Interpolating L, M, S for a given GA"
  )

# Old method without interpolation

df <- df %>%
  mutate(gestational_age_weeks_rounded = round(gestational_age/7),
         gestational_age_weeks_floor = floor(gestational_age/7))

# df <- merge(df, fenton_lms_table %>%
#               rename(gestational_age_weeks_rounded = wk,
#                      sex = Sex), 
#               by = c("gestational_age_weeks_rounded", "sex"), all.x = TRUE, all.y = FALSE)

df <- merge(df, fenton_lms_table %>% 
              rename(gestational_age_weeks_floor = wk,
       sex = Sex), 
by = c("gestational_age_weeks_floor", "sex"), all.x = TRUE, all.y = FALSE)

sum(df$L == 0, na.rm = TRUE)
# nothing is 0 so we can just use the regular formula for all cases
df <- df %>%
  mutate(
    bwz_fen_old = zscore_fromLMS(birth_weight_kg*1000, L, M, S)
  ) %>%
  mutate(
    bwp_fen_old = pnorm(bwz_fen_old)
  )

sum(!is.na(df$bwz_fen_old))
sum(!is.na(df$bwp_fen_old))

hist(df$bwp_ig)
hist(df$bwp_fen_old)
hist(df$bwz_ig)
hist(df$bwz_fen_old)

cor.test(df$bwp_ig, df$bwp_fen_old)

ggplot(df, aes(x = bwp_fen, y = bwp_fen_old)) +
  geom_point(color = "gray", alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  theme_classic() +
  labs(
    x = "BW Percentile Fenton 2013, Linearly Interpolating L, M, S for a given GA",
    y = "BW Percentile Fenton 2013, Using Complete Weeks"
  )

ggplot(df, aes(x = bwp_ig, y = bwp_fen_old)) +
  geom_point(color = "gray", alpha = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  theme_classic() +
  labs(
    x = "BW Percentile Intergrowth",
    y = "BW Percentile Fenton 2013, Using Complete Weeks"
  )

# Correlation Plot with ggcorrplot
# correlations between CMDs across three centiles
cor_mat <- cor(df[,c("bwp_ig", "bwp_fen_old", "bwp_fen", "gestational_age")], use = "pairwise.complete.obs")

p_mat <- cor_pmat(df[,c("bwp_ig", "bwp_fen_old", "bwp_fen", "gestational_age")], use = "pairwise.complete.obs")

(bwp_ga_corrs_plot <- ggcorrplot(
  cor_mat,
  p.mat = p_mat,
  type = "lower",
  insig = "blank",
  lab = TRUE
))