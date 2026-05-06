#library(growthstandards)
library(dplyr)
library(magrittr)
library(ggplot2)


clean_centiles_path <- "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/median-Th-MPR-allgrades-25.11/data_clean/centiles_distinct_2026-03-18.csv"
df <- read.csv(clean_centiles_path) %>% select(-("X"))

bw_ga_table <- read.csv("/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/data_files/Fenton_2013_bw_LMS_upload.csv")

bwga_fct_female <- with(bw_ga_table %>% filter(Sex == "Female"), approxfun(x = M, y = as.numeric(wk)))
bwga_fct_male <- with(bw_ga_table %>% filter(Sex == "Male"), approxfun(x = M, y = as.numeric(wk)))


# Predict GA from birthweight

df[df$sex == "Female", "ga_pred"] <- bwga_fct_female(df[df$sex == "Female","birth_weight_kg"]*1000)
df[df$sex == "Male", "ga_pred"] <- bwga_fct_male(df[df$sex == "Male","birth_weight_kg"]*1000)
# set predicted ga for out of range birthweights to 41 for high and 22 for low out of range
df[is.na(df$ga_pred) & !is.na(df$birth_weight_kg) & df$birth_weight_kg > 3.6, "ga_pred"] <- 41
df[is.na(df$ga_pred) & !is.na(df$birth_weight_kg) & df$birth_weight_kg < 0.5, "ga_pred"] <- 22
# set any predicted ga > 41 to 41, do not want to inflate 42 weekers (very rare)
df[!is.na(df$ga_pred) & df$ga_pred > 41, "ga_pred"] <- 41
sum(is.na(df$ga_pred) & !is.na(df$birth_weight_kg)) # should be 0

ggplot(df %>% filter(!is.na(ga_pred) & !is.na(gestational_age)), aes(x = ga_pred, y = gestational_age)) +
  geom_point() +
  labs(x = "GA Pred Using Fenton2013 completed weeks w/ Interpolation (weeks)", y = "Ground Truth GA (weeks)",
       title = "For Whole Sample with Birthweight") +
  #scale_x_discrete(breaks = c(seq(22, 41, 1))) +
  #scale_y_discrete(breaks = c(seq(27, 42, 1))) +
  theme_classic()

cor.test(df$gestational_age, df$ga_pred)

df$ga_imputed <- ifelse(!is.na(df$gestational_age), df$gestational_age,
                        df$ga_pred)

sum(!is.na(df$ga_imputed))
sum(!is.na(df$gestational_age))