save_folder <- "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/vars/"
data_path <- "/mnt/isilon/bgdlab_processing/braincharts/SLIP/2025_03/code/gamlss/OUT/eren/"
centiles <- read.csv(paste0(data_path, "slip_median_centiles_20250912.csv"))

vol_vars <- c(names(centiles)[grepl("aparc_GrayVol", names(centiles))],
                names(centiles)[grepl("aseg", names(centiles))])
vol_vars <- vol_vars[!vol_vars %in% c("aseg_Left_Cerebral_Cortex", "aseg_Left_Cerebral_White_Matter",
                                        "aseg_Right_Cerebral_Cortex", "aseg_Right_Cerebral_White_Matter",
                                        "aseg_CSF" )]
sa_vars <- c(names(centiles)[grepl("aparc_SurfArea", names(centiles))])
th_vars <- c(names(centiles)[grepl("aparc_ThickAvg", names(centiles))])
global_vars <- c(names(centiles)[grepl("global", names(centiles))])
all_vars <- c(vol_vars, sa_vars, th_vars, global_vars)


if (!dir.exists(save_folder)) dir.create(save_folder)  
write.csv(all_vars, paste0(save_folder, "all_vars", ".csv"))
write.csv(global_vars, paste0(save_folder, "global_vars", ".csv"))
write.csv(vol_vars, paste0(save_folder, "vol_vars", ".csv"))
write.csv(sa_vars, paste0(save_folder, "sa_vars", ".csv"))
write.csv(th_vars, paste0(save_folder, "th_vars", ".csv"))
