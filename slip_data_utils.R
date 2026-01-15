slip_centdf_clean <- function (cent_df, participant_df, qc_df, type = c("median", "mpr")) {
  ##Merge centile tables with participant tables
  ##overall including repeat participant IDs
  cent_df <- merge(cent_df, participant_df %>% select(c("participant_id", 
                                                          "session_id", 
                                                          "age_at_scan", 
                                                          "adjusted_age_in_days",
                                                          "avg_grade",
                                                          "height_at_enc",
                                                          "weight_at_enc",
                                                          "year_of_scan")), by = c("participant_id", "session_id"), 
                     all.x = TRUE, all.y = FALSE)
  
  cent_df <- merge(cent_df, participant_df %>% select(c("participant_id",
                                                          "gestational_age",
                                                          "birth_weight_kg",
                                                          "birth_length_cm")), by = c("participant_id"), 
                     all.x = TRUE, all.y = FALSE)
  
  ##Merge centile tables with QC metrics
  ##group by session-id and take median of euler_mean
  ##merge centile tables with QC metrics if calculated from median of scans
  #MPRage centiles already have corresponding euler_mean for scan
  if (type == "median"){
     qc_df <- qc_df %>%
      group_by(session_id) %>%
      summarise(median_euler_mean = first(na.omit(median_euler_mean)),
                .groups = "drop")
    
    cent_df <- merge(cent_df, qc_df %>% select(c("session_id", 
                                                      "median_euler_mean")), by = c("session_id"), 
                     all.x = TRUE, all.y = FALSE)
  }
  
  ##New Variables
  
  ###Age Variables
  cent_df$age_years <- cent_df$age_days/365
  cent_df$adjusted_age_years <- cent_df$adjusted_age_in_days/365
  ###Log transform Age (ln)
  cent_df$age_days_log <- log(cent_df$age_days)
  cent_df$adjusted_age_days_log <- log(cent_df$adjusted_age_in_days)
  
  #Assign preterm categories
  cent_df$preterm <- as.factor(cut(cent_df$gestational_age, breaks = c(-Inf, 31.9, 36.9, Inf), labels = c("VPM", "LPM", "Term"), include.lowest = TRUE))
  #Calculate birth weight percentile
  cent_df <- cent_df %>%
    mutate(gestAge_days = gestational_age*7+3) %>%
    mutate(sex_recode= recode(sex, "M" = "Male", "F" = "Female"))
  cent_df$bweight_percentile <- igb_wtkg2centile(cent_df$gestAge_days, cent_df$birth_weight_kg, sex = cent_df$sex_recode)/100

  ##Subset to distinct participant IDS
  cent_df_distinct <- cent_df %>% distinct(participant_id, .keep_all = TRUE)
  
  return(list(distinct = cent_df_distinct, all = cent_df))
}

slip_bwga_qc <- function(cent_df, bw_mean, bw_sd, bw_sd_limit, ga_low_limit = 22, ga_high_limit = 42){
  bw_high_limit <- bw_mean + bw_sd_limit*bw_sd
  bw_low_limit <- bw_mean - bw_sd_limit*bw_sd
  print(paste0("BW HIGH LIMIT= ", bw_high_limit))
  print(paste0("BW LOW LIMIT= ", bw_low_limit))
  
  print(paste0("GA HIGH LIMIT= ", ga_high_limit))
  print(paste0("GA LOW LIMIT= ", ga_low_limit))
  
  #BW
  bw_high_removed <- sum(cent_df$birth_weight_kg > bw_high_limit, na.rm = TRUE)
  bw_low_removed <- sum(cent_df$birth_weight_kg < bw_low_limit, na.rm = TRUE)
  print(paste0("Number of participants set to NA because of BW ABOVE qc limit: ", bw_high_removed))
  print(paste0("Number of participants set to NA because of BW BELOW qc limit: ", bw_low_removed))
  cent_df$birth_weight_kg[cent_df$birth_weight_kg > bw_mean + bw_sd_limit*bw_sd] <- NA
  cent_df$birth_weight_kg[cent_df$birth_weight_kg < bw_mean - bw_sd_limit*bw_sd] <- NA
  
  #GA
  ga_high_removed <- sum(cent_df$gestational_age > ga_high_limit, na.rm = TRUE)
  ga_low_removed <- sum(cent_df$gestational_age < ga_low_limit, na.rm = TRUE)
  print(paste0("Number of participants set to NA because of GA ABOVE qc limit: ", ga_high_removed))
  print(paste0("Number of participants set to NA because of GA BELOW qc limit: ", ga_low_removed))
  cent_df$gestational_age[cent_df$gestational_age<ga_low_limit] <- NA
  cent_df$gestational_age[cent_df$gestational_age>ga_high_limit] <- NA
  
  ##set preterm category to NA if GA is now NA, and bweight percentile to NA if bw is now NA
  cent_df$preterm[is.na(cent_df$gestational_age)] <- NA
  cent_df$bweight_percentile[is.na(cent_df$gestational_age) | is.na(cent_df$birth_weight_kg)] <- NA
  
  return(cent_df)
}