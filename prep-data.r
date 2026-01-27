library(dplyr)
library(tidyr)
library(stringr)
setwd('/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/SLIP_25_11/')
main_data_folder <- "/mnt/isilon/bgdlab_processing/releases_clinical/2025_11_release/"

#Read in raw participants file
slip_participants <- read.csv(paste0(main_data_folder,"BIDS/participants.tsv"), sep='\t') 

#Check BW/GA missing values
slip_participants %>%
  group_by(participant_id) %>%
  summarise(n_var = n_distinct(birth_length_cm)) %>%
  filter(n_var > 1)
slip_participants %>%
  group_by(participant_id) %>%
  summarise(n_var = n_distinct(birth_weight_kg)) %>%
  filter(n_var > 1)
slip_participants %>%
  group_by(participant_id) %>%
  summarise(n_var = n_distinct(gestational_age)) %>%
  filter(n_var > 1)

#Replace missing bw/ga values for participants across sessions -- if more than one value recorded, then take the mean
slip_participants <- slip_participants %>%
  group_by(participant_id) %>%
  mutate(
    gestational_age = if (all(is.na(gestational_age))) NA_real_
    else mean(gestational_age, na.rm = TRUE)
  ) %>%
  ungroup()

slip_participants <- slip_participants %>%
  group_by(participant_id) %>%
  mutate(
    birth_weight_kg = if (all(is.na(birth_weight_kg))) NA_real_
    else mean(birth_weight_kg, na.rm = TRUE)
  ) %>%
  ungroup()

slip_participants <- slip_participants %>%
  group_by(participant_id) %>%
  mutate(
    birth_length_cm = if (all(is.na(birth_length_cm))) NA_real_
    else mean(birth_length_cm, na.rm = TRUE)
  ) %>%
  ungroup()

# mutate age column
slip_participants <- slip_participants %>%
    mutate(
        age_days = ifelse(is.na(gestational_age), age_at_scan + 280, age_at_scan + (gestational_age * 7))
    ) %>%
    filter(useable) %>%
    select(participant_id, subject_id, session_id, sex, age_days, device_serial_number, magnetic_field_strength, avg_grade)

#Read QC data
qc = read.csv(paste0(main_data_folder,"derivatives/recon-all-clinical/fs7.4.1-clinical-reconall-qc.csv")) %>%
  mutate(
    threshold = if_else(
      do.call(pmin, c(across(starts_with("qc")), list(na.rm = TRUE))) > 0.65 &
        euler_mean > -60,
      "pass",
      "fail"
    )) %>%
  filter(threshold=='pass') %>%
  select(subject_id, session_id, scan_id)

slip_qc_participants <- inner_join(slip_participants, qc)
slip_qc_scans <- unique(slip_qc_participants$scan_id)


synthseg = read.csv(paste0(main_data_folder,"derivatives/recon-all-clinical/fs7.4.1-clinical-reconall-ss.csv")) %>%
  rename_with(~ paste0("ss_", .x), .cols = -c(subject_id, session_id, scan_id)) %>%
  filter(scan_id %in% slip_qc_scans)
aseg = read.csv(paste0(main_data_folder,"derivatives/recon-all-clinical/fs7.4.1-clinical-reconall-aseg.csv")) %>%
  rename_with(~ paste0("aseg_", .x), .cols = -c(subject_id, session_id, scan_id)) %>%
  filter(scan_id %in% slip_qc_scans)
glasser = read.csv(paste0(main_data_folder,"derivatives/recon-all-clinical/fs7.4.1-clinical-reconall-glasser.csv")) %>%
  select(subject_id, session_id, scan_id, starts_with('GrayVol'), starts_with('SurfArea'),starts_with('ThickAvg')) %>%
  rename_with(~ paste0("glasser_", .x), .cols = -c(subject_id, session_id, scan_id)) %>%
  filter(scan_id %in% slip_qc_scans)
aparc = read.csv(paste0(main_data_folder,"derivatives/recon-all-clinical/fs7.4.1-clinical-reconall-aparc.csv")) %>%
  select(subject_id, session_id, scan_id, starts_with('GrayVol'), starts_with('SurfArea'),starts_with('ThickAvg'))%>%
  rename_with(~ paste0("aparc_", .x), .cols = -c(subject_id, session_id, scan_id)) %>%
  filter(scan_id %in% slip_qc_scans)
global = read.csv(paste0(main_data_folder,"derivatives/recon-all-clinical/fs7.4.1-clinical-reconall-global.csv")) %>%
  rename_with(~ paste0("global_", .x), .cols = -c(subject_id, session_id, scan_id)) %>%
  filter(scan_id %in% slip_qc_scans)

slip_scans <- slip_qc_participants %>%
  left_join(global) %>%
  left_join(synthseg) %>%
  left_join(glasser) %>%
  left_join(aparc) %>%
  left_join(aseg) %>%
  mutate(
    global_MeanThickness = rowMeans(select(., starts_with("aparc_ThickAvg")), na.rm = TRUE),
    global_SurfaceArea = rowSums(select(., starts_with("aparc_SurfArea")), na.rm = TRUE)
  )

regions <- colnames(slip_scans)[9:1455]
regions2 <- colnames(slip_scans)[1454:1455]

write.table(regions, 'data-prep/regions.txt', quote=F, row.names=F)
write.table(regions2, 'data-prep/regions2.txt', quote=F, row.names=F)


slip_median <- slip_scans %>%
  filter(sex != 'Unknown')%>%
  rename(site = device_serial_number) %>%
  group_by(participant_id, session_id, sex, age_days, site, avg_grade) %>%
  summarise(across(
    .cols = where(is.numeric), 
    .fns = ~median(.x, na.rm = TRUE)
  ), .groups = "drop")

slip_median <- slip_median %>%
  group_by(participant_id) %>%
  mutate(session= dense_rank(interaction(session_id, age_days)),
         dx='CN') %>%
  ungroup()


slip_mpr <- slip_scans %>%
  filter(str_detect(scan_id, "MPR")) %>%
  filter(sex != 'Unknown')%>%
  rename(site = device_serial_number) %>%
  group_by(participant_id, session_id) %>%
  slice(1) %>%
  ungroup() %>%
  group_by(participant_id) %>%
  mutate(session = dense_rank(interaction(session_id, age_days)),
         dx='CN')

slip_median2 <- slip_median %>%
  filter(avg_grade == 2)

slip_mpr2 <- slip_mpr %>%
  filter(avg_grade == 2)

write.csv(slip_median, 'CSV/slip_median.csv', row.names=F)
write.csv(slip_mpr, 'CSV/slip_mpr.csv', row.names=F)
write.csv(slip_median2, 'CSV/slip_median2.csv', row.names=F)
write.csv(slip_mpr2, 'CSV/slip_mpr2.csv', row.names=F)
