slip_centdf_clean <- function (cent_df, participant_df, qc_df, 
                               type = c("median", "mpr"), add_qc = TRUE,
                               participants_vars_select,
                               cent_df_vars_select) {
  ##Merge centile tables with participant tables
  ##overall including repeat participant IDs
  cent_df <- merge(cent_df %>% select(any_of(cent_df_vars_select)), 
                   participant_df %>% select(all_of(participants_vars_select)), 
                   by = c("participant_id", "session_id", "sex"), 
                     all.x = TRUE, all.y = FALSE)
  
  ##Merge centile tables with QC metrics
  ##group by session_id and take median of euler_mean
  ##merge centile tables with QC metrics if calculated from median of scans
  #MPRage centiles already have corresponding euler_mean matched on scan_id
  if(add_qc == TRUE){
    if (type == "median"){
       qc_df <- qc_df %>%
        group_by(session_id) %>%
        summarise(median_euler_mean = first(na.omit(median_euler_mean)),
                  .groups = "drop")
      
      cent_df <- merge(cent_df, qc_df %>% select(c("session_id", 
                                                        "median_euler_mean")), 
                       by = c("session_id"), all.x = TRUE, all.y = FALSE)
      
    }  else if (type == "mpr" & !("euler_mean" %in% names(cent_df))) {
      cent_df <- merge(cent_df, qc_df %>% select(c("scan_id", 
                                                   "euler_mean")), by = c("scan_id"), 
                       all.x = TRUE, all.y = FALSE)
    }
  }
  
  ##New Variables
  cent_df <- cent_df %>%
    rename(age_days_adj = age_days)
  
  ###Age Variables
  cent_df$age_at_scan_years <- cent_df$age_at_scan/365.25
  cent_df$age_years_adj <- cent_df$age_days_adj/365.25
  
  #Assign preterm categories
  cent_df$preterm <- as.factor(cut(cent_df$gestational_age, 
                                   breaks = c(-Inf, 31.9*7, 36.9*7, Inf), 
                                   labels = c("VPM", "LPM", "Term"), include.lowest = TRUE))
  
  #Calculate birth weight percentile
  fenton_table_path <- "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/data_files/Fenton_2013_bw_LMS_upload.csv"
  fenton_lms_table <- read.csv(fenton_table_path)

  cent_df <- merge(cent_df %>%
                     mutate(gestational_age_weeks_complete = floor(gestational_age/7)),
                   fenton_lms_table %>% 
                     rename(gestational_age_weeks_complete = wk,
                            sex = Sex) %>%
                     select(L,M,S,gestational_age_weeks_complete,sex), 
                   by = c("gestational_age_weeks_complete", "sex"), all.x = TRUE, all.y = FALSE) %>%
    mutate(bwz_fen = ((((birth_weight_kg*1000/M)^L)-1)/(L*S)) ) %>%
    mutate(bwp_fen = pnorm(bwz_fen) ) %>%
    select(-c("L", "M", "S"))

  ##Subset to one data point per participant, youngest scan available
  cent_df_distinct <- cent_df %>% 
    group_by(participant_id) %>%
    slice_min(age_days_adj, n = 1, with_ties = FALSE) %>% 
    ungroup()
  
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