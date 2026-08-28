# Script to get phecodes and merge with centiles

# Script to get median weight/height within 180 days of scan for 25.11 MRI cohort, should be matched to scan date, not participant


#Load packages. Install these packages with install.packages("pkg_name") if they are not installed.
library(bigrquery) #package to import SQL tables into R
library(tidyverse) #package suite for data manipulation, reading data into R, and plotting (ggplot2)
library(glue) #package to use expressions with strings
library(ggseg)
require(gridExtra)
library(interactions)
library(boot)
library(parameters)
library(PheWAS)

panel_save <- function(plot, filename, extension, folder){
  save_name <- paste0(folder, filename, ".", extension)
  ggsave(save_name, plot = plot, dpi = 300)
}

bq_auth()
# proj_id <- if (exists('proj_id')) proj_id else bq_projects()[1]
proj_id = "scit605-healthymri-2d73d69f"

# read in centiles to subset release
cent_path <- "/mnt/arcus/lab/users/kafadare/slip_centiles_processed_tables/df_zscore_distinct.csv"
cent <- read.csv(cent_path) %>%
  mutate(proc_ord_id_wrong = proc_ord_id, 
         proc_ord_id = sub(".*ses-([0-9]+)procId.*", "\\1", session_id))

pat_id_list <- paste0("'", sapply(strsplit(cent$pat_id, ","),trimws), "'", collapse = ", ")


phecode_query <- glue('
  SELECT DISTINCT
    p.*,
  FROM lab.phecodes_20260228 AS p
  WHERE p.pat_id IN ({pat_id_list});
')

phecode <- bq_project_query(proj_id, phecode_query) %>%
  bq_table_download(., page_size = 3500)

print(paste0(
  "lab.phecodes_20260228 distinct rows, dimensions: ",
  dim(phecode)
))

phecode$pat_phecode <- paste(phecode$pat_id, phecode$phecode, sep = "_")
phecode <- addPhecodeInfo(phecode, groupcolors = TRUE, groupnums = TRUE)

phecode <- merge(phecode, cent %>%
                   select(pat_id, gestational_age, birth_weight_kg, bwp_fen,
                          global_eTIV, global_VentricleChoroidVol, global_SubCortGrayVol),
                 by = "pat_id", all.x = TRUE, all.y = FALSE)

sum(duplicated(phecode$pat_phecode))/nrow(phecode)

test <- phecode %>%
  group_by(pat_phecode) %>%
  slice_min(age)

sum(duplicated(test$pat_phecode))/nrow(test)

test_distinct <- phecode %>%
  distinct(pat_phecode, .keep_all = TRUE)

# Convert to wide df with 0/1 for each phecode-patient pair
phewas_phecode_df <- test_distinct %>%
  filter(!is.na(gestational_age)) %>%
  select(pat_id, phecode) %>%
  rename(id = pat_id) %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = phecode, 
    values_from = value, 
    values_fill = 0
  )

### Gestational Age

phewas_predictor_df <- test_distinct %>%
  distinct(pat_id, .keep_all = TRUE) %>%
  filter(!is.na(gestational_age)) %>%
  select(pat_id, gestational_age) %>%
  rename(id = pat_id)

results_ga <- phewas(
  phenotypes = phewas_phecode_df,
  genotypes = phewas_predictor_df,
  predictor = "gestational_age",
  cores = 4
)

phewasManhattan(
  results_ga %>%
    select(-OR) %>%
    rename(OR = beta), 
  title = "PheWAS: Gestational Age at Birth",
  OR.direction = TRUE,
  OR.size = TRUE 
) +
  scale_size_continuous(
    range = c(1, 3),
    limits = c(-0.03, 0.07),
    breaks = seq(-0.03, 0.07, length.out = 10)
  )


### BWP
phewas_predictor_df <- test_distinct %>%
  distinct(pat_id, .keep_all = TRUE) %>%
  filter(!is.na(bwp_fen)) %>%
  select(pat_id, bwp_fen) %>%
  rename(id = pat_id)

results_bwp <- phewas(
  phenotypes = phewas_phecode_df,
  genotypes = phewas_predictor_df,
  predictor = "bwp_fen",
  cores = 4
)

phewasManhattan(
  results_bwp %>%
    select(-OR) %>%
    rename(OR = beta), 
  title = "PheWAS: Birth Weight Percentile",
  OR.direction = TRUE,
  OR.size = TRUE 
)


### Brain Size (eTIV)
phewas_predictor_df <- test_distinct %>%
  distinct(pat_id, .keep_all = TRUE) %>%
  filter(!is.na(global_eTIV)) %>%
  select(pat_id, global_eTIV) %>%
  rename(id = pat_id)

results_etiv <- phewas(
  phenotypes = phewas_phecode_df,
  genotypes = phewas_predictor_df,
  predictor = "global_eTIV",
  cores = 4
)

phewasManhattan(
  results_etiv %>%
    select(-OR) %>%
    rename(OR = beta), 
  title = "PheWAS: eTIV centile",
  OR.direction = TRUE,
  OR.size = TRUE 
)


### CMDs
# volume
phewas_predictor_df <- cent %>%
  filter(!is.na(vol.cmd)) %>%
  select(pat_id, vol.cmd) %>%
  rename(id = pat_id)

results_volcmd <- phewas(
  phenotypes = phewas_phecode_df,
  genotypes = phewas_predictor_df,
  predictor = "vol.cmd",
  cores = 4
)

phewasManhattan(
  results_volcmd %>%
    select(-OR) %>%
    rename(OR = beta), 
  title = "PheWAS: Volume CMD",
  OR.direction = TRUE,
  OR.size = TRUE 
)

# sa
phewas_predictor_df <- cent %>%
  filter(!is.na(sa.cmd)) %>%
  select(pat_id, sa.cmd) %>%
  rename(id = pat_id)

results_sacmd <- phewas(
  phenotypes = phewas_phecode_df,
  genotypes = phewas_predictor_df,
  predictor = "sa.cmd",
  cores = 4
)

phewasManhattan(
  results_sacmd %>%
    select(-OR) %>%
    rename(OR = beta), 
  title = "PheWAS: SA CMD",
  OR.direction = TRUE,
  OR.size = TRUE 
)


# th
phewas_predictor_df <- cent %>%
  filter(!is.na(th.cmd)) %>%
  select(pat_id, th.cmd) %>%
  rename(id = pat_id)

results_thcmd <- phewas(
  phenotypes = phewas_phecode_df,
  genotypes = phewas_predictor_df,
  predictor = "th.cmd",
  cores = 4
)

phewasManhattan(
  results_thcmd %>%
    select(-OR) %>%
    rename(OR = beta), 
  title = "PheWAS: th CMD",
  OR.direction = TRUE,
  OR.size = TRUE 
)


# subcort
phewas_predictor_df <- cent %>%
  filter(!is.na(subcort.cmd)) %>%
  select(pat_id, subcort.cmd) %>%
  rename(id = pat_id)

results_subcortcmd <- phewas(
  phenotypes = phewas_phecode_df,
  genotypes = phewas_predictor_df,
  predictor = "subcort.cmd",
  cores = 4
)

phewasManhattan(
  results_subcortcmd %>%
    select(-OR) %>%
    rename(OR = beta), 
  title = "PheWAS: subcort CMD",
  OR.direction = TRUE,
  OR.size = TRUE 
)

phecode <- merge(phecode,PheWAS::pheinfo %>%
                           select(phecode, description, group, color) %>%
                           rename(disease_group = group, disease_color = color,
                                  disease_name = description) %>%
                           distinct(),
                         by = "phecode",
                         all.x = TRUE, all.y = FALSE
)


# Check neonatal diagnoses
phecode_neonate <- phecode %>% filter(age < 28)

table((phecode_neonate$disease_name))

table((phecode_neonate$disease_group))

length(unique(phecode_neonate$pat_id)) # 1695

pat_id_neonate_dx <- phecode_neonate %>% pull(pat_id)

pat_id_pregnancy_comp <- phecode_neonate %>% pull(pat_id)

cent_noneo_dx <- cent %>% filter(!(pat_id %in% pat_id_neonate_dx))

cor.test(cent_noneo_dx$global_eTIV, cent_noneo_dx$gestational_age)
cor.test(cent$global_eTIV, cent$gestational_age)

cor.test(cent_noneo_dx$global_eTIV, cent_noneo_dx$bwp_fen)
cor.test(cent$global_eTIV, cent$bwp_fen)


# Add true/false column of ANY neonatal diagnosis
cent$neonatal_dx <- ifelse(cent$pat_id %in% pat_id_neonate_dx, 1, 0)

# Save with neonatal diagnosis

write.csv(cent, cent_path)

# Make table of neuropsych diagnoses

phecode_psych <- phecode %>% filter(disease_group %in% c("mental disorders"))

table((phecode_psych$disease_name))

length(unique(phecode_psych$pat_id)) # 7895

sum(duplicated(phecode_psych$pat_phecode))/nrow(phecode_psych)

phecode_psych_minage <- phecode_psych %>%
  group_by(pat_phecode) %>%
  slice_min(age)

sum(duplicated(phecode_psych_minage$pat_phecode))/nrow(phecode_psych_minage)

phecode_psych_minage <- phecode_psych_minage %>%
  distinct(pat_phecode, .keep_all = TRUE)

cent_psych <- merge(cent,
                    phecode_psych_minage %>%
                      select(phecode, pat_id, encounter_id, age, disease_name,
                             disease_group, disease_color) %>%
                      rename(age_min_psych_dx = age), 
                    by = "pat_id", all.x = TRUE, all.y = FALSE)


write.csv(cent_psych, "/mnt/arcus/lab/users/kafadare/slip_centiles_processed_tables/df_zscore_distinct_psychdx.csv")
