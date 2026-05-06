require(bigrquery)
require(tidyr)
require(magrittr)
require(glue)
require(growthstandards)
bq_auth()
proj_id <- "scit605-healthymri-2d73d69f"

#Download slip 2511 release table
query <- glue("SELECT * FROM lab.release_2025_11") #the string here is the SQL query
df <- bq_project_query(proj_id, query) %>% bq_table_download(., page_size=3500) #assign the name of the dataframe in R here
print(paste0("Arcus lab.release_2025_11 table, dimensions: ", dim(df))) #pastes the dimensions of the table downloaded

#Download full patients table where patient IDs are found in the slip release
query2 <- glue("SELECT * FROM arcus.patient AS p WHERE p.pat_id IN (SELECT DISTINCT pat_id FROM lab.release_2025_11)")
patients <- bq_project_query(proj_id, query2) %>% bq_table_download(., page_size=3500) #assign the name of the dataframe in R here
print(paste0("Arcus arcus.patients table intersection with lab.release_2025_11, dimensions: ", dim(patients))) #pastes the dimensions of the table downloaded


df <- merge(df %>% rename(gestational_age_old = gestational_age, birth_weight_kg_old = birth_weight_kg, birth_length_cm_old = birth_length_cm), 
            patients %>% select(c("pat_id", "gestational_age", "gestational_age_num", "birth_weight_kg", "birth_length_cm")) %>% 
              rename(gestational_age_string = gestational_age),
            by = "pat_id", all.x = TRUE, all.y = FALSE)


sum(is.na(df$gestational_age_old) & !is.na(df$gestational_age_num)) #117
sum(is.na(df$gestational_age_old) & !is.na(df$gestational_age_string)) #117
sum(is.na(df$birth_weight_kg_old) & !is.na(df$birth_weight_kg)) #104 
sum(is.na(df$birth_length_cm_old) & !is.na(df$birth_length_cm)) #12

# Check if there are mismatch entries for the same patient
## only birth_weight_kg_old and gestational_age_old has mismatches
df %>%
  group_by(pat_id) %>%
  summarise(n_var_bl = n_distinct(birth_length_cm_old)) %>%
  filter(n_var_bl > 1)
df %>%
  group_by(pat_id) %>%
  summarise(n_var_bl = n_distinct(birth_length_cm)) %>%
  filter(n_var_bl > 1)
df %>%
  group_by(pat_id) %>%
  summarise(n_var_bw = n_distinct(birth_weight_kg_old)) %>%
  filter(n_var_bw > 1)
df %>%
    group_by(pat_id) %>%
    summarise(n_var_bw = n_distinct(birth_weight_kg)) %>%
    filter(n_var_bw > 1)
df %>%
  group_by(pat_id) %>%
  summarise(n_var_ga = n_distinct(gestational_age_old)) %>%
  filter(n_var_ga > 1)
df %>%
  group_by(pat_id) %>%
  summarise(n_var_ga = n_distinct(gestational_age_num)) %>%
  filter(n_var_ga > 1)
df %>%
  group_by(pat_id) %>%
  summarise(n_var_ga = n_distinct(gestational_age_string)) %>%
  filter(n_var_ga > 1)

# Deal with missing cells -- in this case only the "old" (OG) columns have this issue so commenting it out and will only use the newly merged columns going ahead.
# df <- df %>%
#   group_by(pat_id) %>%
#   mutate(
#     gestational_age_old = if (all(is.na(gestational_age_old))) NA_real_
#     else mean(gestational_age_old, na.rm = TRUE),
#     birth_weight_kg_old = if (all(is.na(birth_weight_kg_old))) NA_real_
#     else mean(birth_weight_kg_old, na.rm = TRUE)
#   ) %>%
#   ungroup()

# How many new patients have GA/BW/BL info that didn't have it in the old column with this new merge
## This is mainly to keep track/see how many people would be affected, can be deleted from final code
length(unique(df[which(is.na(df$gestational_age_old) & !is.na(df$gestational_age_string)),"pat_id"])) #89
pat_id_newGA <- unique(df[which(is.na(df$gestational_age_old) & !is.na(df$gestational_age_string)),"pat_id"])
length(unique(df[which(is.na(df$birth_weight_kg_old) & !is.na(df$birth_weight_kg)),"pat_id"])) #63
pat_id_newBW <- unique(df[which(is.na(df$birth_weight_kg_old) & !is.na(df$birth_weight_kg)),"pat_id"])
length(unique(df[which(is.na(df$birth_length_cm_old) & !is.na(df$birth_length_cm)),"pat_id"])) #12
pat_id_newBL <- unique(df[which(is.na(df$birth_length_cm_old) & !is.na(df$birth_length_cm)),"pat_id"])

# no one that had info in an old column is missing info in new column
sum(!is.na(df$gestational_age_old) & is.na(df$gestational_age_string)) #0
sum(!is.na(df$birth_weight_kg_old) & is.na(df$birth_weight_kg)) #0
sum(!is.na(df$birth_length_cm_old) & is.na(df$birth_Length_cm)) #0

#df <- df %>% distinct(pat_id, .keep_all = TRUE)
# look at the distribution for the newly merged gestational age variables
#table(df$gestational_age_string)
#table(df$gestational_age_num)

# how many entries in this sample have fractional information (beyond whole numbers for weeks?)
# sum(grepl("\\+|/|\\.", df$gestational_age_string))/sum(!is.na(df$gestational_age_string))

# Convert different gestational_age_string formats and convert to both 1) number of days AND 2) number of weeks (rounded up/down) AND 3) number of COMPLETED weeks to the best of our knowledge
data <- df$gestational_age_string

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

df$ga_complete_weeks <- complete_weeks
df$ga_days_total <- days_total
df$ga_weeks_rounded <- weeks_rounded

# Checks
# 0.5 is rounded up in ga_weeks_rounded but its rounded to complete weeks in gestational_age_num (original field from arcus.patients)
# Seems that .7 and up are rounded up in gestational_age_num
df %>% filter(gestational_age_num != ga_weeks_rounded) %>% select(gestational_age_string, gestational_age_num, ga_weeks_rounded, ga_complete_weeks, ga_days_total)
df %>% filter(gestational_age_num != ga_complete_weeks) %>% select(gestational_age_string, gestational_age_num, ga_weeks_rounded, ga_complete_weeks, ga_days_total)

# Deal with BW/GA QC (can comment this out if you want to keep nonQC birthweight)
df %>% filter(gestational_age_num > 42) %>% select(gestational_age_old, gestational_age_num, gestational_age_string, ga_weeks_rounded)
# ^ This shows that gestational_age_string and gestational_age_num don't actually match up. 
# Will be using gestational_age_string and values derived from that field for all QC

range(df$birth_weight_kg, na.rm = T)
range(df$ga_days_total/7, na.rm = T)
# GA: remove anyone encoded as > 44*7 ga days or < 22 weeks (there is a single 43.5 weeker, I don't really want to remove this entry)
# BW: remove anyone with birthweight recorded as > 10kg and < 0.3kg
df$ga_days_total_qc <- ifelse(df$ga_days_total < 22*7 | df$ga_days_total >= 44*7, NA, df$ga_days_total)
sum(!is.na(df$ga_days_total) & is.na(df$ga_days_total_qc)) # 2 removed
df$birth_weight_kg_qc <- ifelse(df$birth_weight_kg < 0.3 | df$birth_weight_kg > 10, NA, df$birth_weight_kg)
sum(!is.na(df$birth_weight_kg) & is.na(df$birth_weight_kg_qc)) # 12 removed


# Impute gestational age for entries with birth weight but no gestational age
bw_ga_table <- read.csv("/mnt/arcus/lab/shared/Eren/Fenton_2013_bw_LMS_upload.csv")

bwga_fct_female <- with(bw_ga_table %>% filter(Sex == "Female"), approxfun(x = M, y = as.numeric(wk)))
bwga_fct_male <- with(bw_ga_table %>% filter(Sex == "Male"), approxfun(x = M, y = as.numeric(wk)))

# Predict GA from birthweight (post-QC)

df[df$sex == "Female", "ga_pred"] <- bwga_fct_female(df[df$sex == "Female","birth_weight_kg_qc"]*1000)
df[df$sex == "Male", "ga_pred"] <- bwga_fct_male(df[df$sex == "Male","birth_weight_kg_qc"]*1000)

# Examine distribution of "unmatched" birthweights. Appears mostly higher, but a couple on the lower end
hist(df[is.na(df$ga_pred) & !is.na(df$birth_weight_kg_qc), "birth_weight_kg"])
range(df[is.na(df$ga_pred) & !is.na(df$birth_weight_kg_qc), "birth_weight_kg"])

# set predicted ga for out of range birthweights to 41 for high and 22 for low out of range.
df[is.na(df$ga_pred) & !is.na(df$birth_weight_kg_qc) & df$birth_weight_kg_qc > 3.5, "ga_pred"] <- 41
df[is.na(df$ga_pred) & !is.na(df$birth_weight_kg_qc) & df$birth_weight_kg_qc < 0.55, "ga_pred"] <- 22
# This below should be 0 after setting floor and ceiling
sum(is.na(df$ga_pred) & !is.na(df$birth_weight_kg_qc)) # 0
df$ga_pred_days_total <- df$ga_pred*7

# Can get rid of the plot & cor.test below, it's for validation purposes
ggplot(df %>% filter(!is.na(ga_pred_days_total) & !is.na(ga_days_total_qc)), aes(x = ga_pred_days_total, y = ga_days_total_qc)) +
  geom_point() +
  labs(x = "GA Pred Using Fenton2013 completed weeks w/ Interpolation (days)", y = "Ground Truth GA (days)",
       title = "For Whole Sample with Birthweight") +
  theme_classic()

cor.test(df$ga_days_total_qc, df$ga_pred_days_total)



# Create new column for imputed ga in days total, using known GA if available and predicted GA if not
# if predicted GA is not available, replace with median in the sample
df$ga_imputed_days_total <- ifelse(!is.na(df$ga_days_total_qc), df$ga_days_total_qc,
                        ifelse(!is.na(df$ga_pred), df$ga_pred,
                               median(df$ga_days_total_qc, na.rm = T)))

sum(is.na(df$ga_imputed_days_total)) #0 -- should be 0

# Create new column for imputed ga in days total, using known GA if available and predicted GA if not
# if predicted GA is not available, keep NA
df$ga_imputed_nomedian_days_total <- ifelse(!is.na(df$ga_days_total_qc), df$ga_days_total_qc,
                                 ifelse(!is.na(df$ga_pred_days_total), df$ga_pred_days_total,
                                        NA))

sum(is.na(df$ga_imputed_nomedian_days_total)) #7598
sum(is.na(df$ga_days_total_qc)) #8392


### Note: you don't have to keep all columns I created in the final table, 
# especially the conversions between complete weeks and weeks rounded are easily achievable from the ga_days_total column
# there is very few QC-ed GA entries between ga_days_total and ga_days_total qc, I don't know if this is helpful to keep in the table for other people.
# ^ same for birth_weight_kg and birth_weight_kg_qc columns. I will need to report on these numbers though, but I can come back to this script.
# I have two imputed columns one that includes everyone with known ga + ga predicted from bw, 
# and one also includes using the median of the sample for the completely unknown ga. I think both could be useful.
# since the gestational_age_num column is not used to get the final columns, and the gestational_age_string column doesn't seem to be useful anymore,
# those could be removed (but if anyone wants to see the original column everything is derived from I would keep gestational_age_string)
# gestational_age_old and birth_weight_kg_old could be useful to keep for people whose analyses depend on these original
# columns, so they can replicate their results.