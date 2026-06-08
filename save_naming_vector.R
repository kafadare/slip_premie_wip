# Script to define and save naming vector
require(dplyr)
require(glue)
require(magrittr)

`%nin%` = Negate(`%in%`)

vars_path <-  "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/vars_2511_new/"

#csv_save_path = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/"

full_vars <- read.csv(paste0(vars_path,"all_vars.csv")) %>% pull(x)
global_vars <- read.csv(paste0(vars_path,"global_vars.csv")) %>% pull(x)
vol_vars <- read.csv(paste0(vars_path,"vol_vars.csv")) %>% pull(x)
sa_vars <- read.csv(paste0(vars_path,"sa_vars.csv")) %>% pull(x)
th_vars <- read.csv(paste0(vars_path,"th_vars.csv")) %>% pull(x)

global_names <- global_vars %>%
  sub("global_", "", .) %>%
  sub("lh", "Left", .) %>%
  sub("rh", "Right", .) %>%
  sub("Vol", "Volume", .) %>%
  sub("SubCort", "Subcortical", .) %>%
  gsub("(?<!^)([A-Z])", " \\1", ., perl = TRUE) %>%
  sub("e T I V", "eTIV (Total Intracranial Volume)", .)

vol_names <- vol_vars %>%
  sub("aparc_GrayVol_", "", .) %>%
  sub("lh_", "L ", .) %>%
  sub("rh_", "R ", .) %>%
  sub("Right", "R", .) %>%
  sub("Left", "L", .) %>%
  sub("aseg_", "", .) %>%
  sub("Inf_Lat_Vent", "Inferior Lateral Ventricle", .) %>%
  gsub("_", " ", .)

sa_names <- sa_vars %>%
  sub("aparc_SurfArea_", "", .) %>%
  sub("lh_", "L ", .) %>%
  sub("rh_", "R ", .) %>%
  gsub("_", " ", .)

th_names <- th_vars %>%
  sub("aparc_ThickAvg_", "", .) %>%
  sub("lh_", "L ", .) %>%
  sub("rh_", "R ", .) %>%
  gsub("_", " ", .)


name_mapping <- data.frame(
  key = c(global_vars, vol_vars, sa_vars, th_vars),
  plot_names = c(global_names, vol_names, sa_names, th_names)
) 

filename <- glue("{vars_path}name_mapping.csv")

write.csv(name_mapping, filename)