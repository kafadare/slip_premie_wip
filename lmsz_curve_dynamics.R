require(dplyr)
require(glue)
require(magrittr)
require(gamlssTools)
# For figures/plotting
require(ggplot2)
require(ggseg)
require(egg)
require(patchwork)

# script to find the plateau point for eTIV (and whatever other phenotype of interest)
# set folder paths
code_path <- "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/"
source(glue("{code_path}R_util_lmsz.R"))

vars_path <-  "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/vars/"

output_folder = "/mnt/isilon/bgdlab_processing/Eren/slip_premie_wip/output/lmsz_median_Th-MPR_25.11/"

csv_save_path <-  paste0(output_folder, "results_csv/")
plot_save_path <- paste0(output_folder, "plots/")

full_vars <- read.csv(paste0(vars_path,"all_vars.csv")) %>% pull(x)
global_vars <- read.csv(paste0(vars_path,"global_vars.csv")) %>% pull(x)
vol_vars <- read.csv(paste0(vars_path,"vol_vars.csv")) %>% pull(x)
sa_vars <- read.csv(paste0(vars_path,"sa_vars.csv")) %>% pull(x)
th_vars <- read.csv(paste0(vars_path,"th_vars.csv")) %>% pull(x)

## Load previously calculated centile fans
fans <- read.csv(glue("{csv_save_path}lmsz_centfans_vpm_uncorrected.csv"))

lmsz_vars <- unique(fans$pheno)

#Scale all phenotypes by 1000
fans <- fans %>% 
  mutate(SLIP_input_fans = SLIP_input_fans*1000,
         LMSz_fans = LMSz_fans*1000)

out_df <- data.frame(
  pheno = character(),
  cent_curve = numeric(),
  cent_deriv_adj_age_days_log_slip = numeric(),
  cent_deriv2_adj_age_days_log_slip = numeric(),
  cent_deriv_adj_age_years_slip = numeric(),
  cent_deriv2_adj_age_years_slip = numeric(),
  adj_age_days_log_peak50_slip = numeric(),
  rank_peak50_slip = numeric(),
  rank_centderiv_slip = numeric(),
  rank_centderiv2_slip = numeric(),
  cent_deriv_adj_age_days_log_lmsz = numeric(),
  cent_deriv2_adj_age_days_log_lmsz = numeric(),
  cent_deriv_adj_age_years_lmsz = numeric(),
  cent_deriv2_adj_age_years_lmsz = numeric(),
  adj_age_days_log_peak50_lmsz = numeric(),
  rank_peak50_lmsz = numeric(),
  rank_centderiv_lmsz = numeric(),
  rank_centderiv2_lmsz = numeric(),
  cent_deriv_adj_age_days_log_slip_female = numeric(),
  cent_deriv2_adj_age_days_log_slip_female = numeric(),
  cent_deriv_adj_age_years_slip_female = numeric(),
  cent_deriv2_adj_age_years_slip_female = numeric(),
  adj_age_days_log_peak50_slip_female = numeric(),
  rank_peak50_slip_female = numeric(),
  rank_centderiv_slip_female = numeric(),
  rank_centderiv2_slip_female = numeric(),
  cent_deriv_adj_age_days_log_lmsz_female = numeric(),
  cent_deriv2_adj_age_days_log_lmsz_female = numeric(),
  cent_deriv_adj_age_years_lmsz_female = numeric(),
  cent_deriv2_adj_age_years_lmsz_female = numeric(),
  adj_age_days_log_peak50_lmsz_female = numeric(),
  rank_peak50_lmsz_female = numeric(),
  rank_centderiv_lmsz_female = numeric(),
  rank_centderiv2_lmsz_female = numeric(),
  cent_deriv_adj_age_days_log_slip_male = numeric(),
  cent_deriv2_adj_age_days_log_slip_male = numeric(),
  cent_deriv_adj_age_years_slip_male = numeric(),
  cent_deriv2_adj_age_years_slip_male = numeric(),
  adj_age_days_log_peak50_slip_male = numeric(),
  rank_peak50_slip_male = numeric(),
  rank_centderiv_slip_male = numeric(),
  rank_centderiv2_slip_male = numeric(),
  cent_deriv_adj_age_days_log_lmsz_male = numeric(),
  cent_deriv2_adj_age_days_log_lmsz_male = numeric(),
  cent_deriv_adj_age_years_lmsz_male = numeric(),
  cent_deriv2_adj_age_years_lmsz_male = numeric(),
  adj_age_days_log_peak50_lmsz_male = numeric(),
  rank_peak50_lmsz_male = numeric(),
  rank_centderiv_lmsz_male = numeric(),
  rank_centderiv2_lmsz_male = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:length(lmsz_vars)) {
  ph <- lmsz_vars[i]
  print(ph)
  out_df[i, "pheno"] <- ph
  
  tmp <- fans %>% filter(pheno == ph)
  
  tmp_50 <- tmp %>% filter(centile_cat == "cent_0.5") %>%
    group_by(adjusted_age_days_log) %>%
    summarise(sex_average_slip = mean(SLIP_input_fans),
              sex_average_lmsz = mean(LMSz_fans)) %>%
    ungroup() %>%
    mutate(adjusted_age_years = exp(adjusted_age_days_log)/365.25)
  
  tmp_50_f <- tmp %>% filter(centile_cat == "cent_0.5" & sex == "Female") %>%
    mutate(adjusted_age_years = exp(adjusted_age_days_log)/365.25)
  
  tmp_50_m <- tmp %>% filter(centile_cat == "cent_0.5" & sex == "Male") %>%
    mutate(adjusted_age_years = exp(adjusted_age_days_log)/365.25)
  
  dx <- diff(tmp_50$adjusted_age_days_log)
  dy_slip <- diff(tmp_50$sex_average_slip)
  dy_lmsz <- diff(tmp_50$sex_average_lmsz)
  
  dy_slip_f <- diff(tmp_50_f$SLIP_input_fans)
  dy_lmsz_f <- diff(tmp_50_f$LMSz_fans)
  
  dy_slip_m <- diff(tmp_50_m$SLIP_input_fans)
  dy_lmsz_m <- diff(tmp_50_m$LMSz_fans)
  
  tmp_50$cent_deriv_slip <- c(99999, (dy_slip / dx))
  tmp_50$cent_deriv2_slip <- c(99999, 99999, diff(dy_slip / dx))
  tmp_50$cent_deriv_lmsz <- c(99999, (dy_lmsz / dx))
  tmp_50$cent_deriv2_lmsz <- c(99999, 99999, diff(dy_lmsz / dx))
  
  tmp_50_f$cent_deriv_slip <- c(99999, (dy_slip_f / dx))
  tmp_50_f$cent_deriv2_slip <- c(99999, 99999, diff(dy_slip_f / dx))
  tmp_50_f$cent_deriv_lmsz <- c(99999, (dy_lmsz_f / dx))
  tmp_50_f$cent_deriv2_lmsz <- c(99999, 99999, diff(dy_lmsz_f / dx))
  
  tmp_50_m$cent_deriv_slip <- c(99999, (dy_slip_m / dx))
  tmp_50_m$cent_deriv2_slip <- c(99999, 99999, diff(dy_slip_m / dx))
  tmp_50_m$cent_deriv_lmsz <- c(99999, (dy_lmsz_m / dx))
  tmp_50_m$cent_deriv2_lmsz <- c(99999, 99999, diff(dy_lmsz_m / dx))
  
  # Sex averaged
  out_df[i, "cent_deriv_adj_age_days_log_slip"] <- tmp_50[which.min(tmp_50$cent_deriv_slip), "adjusted_age_days_log"]
  out_df[i, "cent_deriv2_adj_age_days_log_slip"] <- tmp_50[which.min(tmp_50$cent_deriv2_slip), "adjusted_age_days_log"]
  out_df[i, "cent_deriv_adj_age_years_slip"] <- tmp_50[which.min(tmp_50$cent_deriv_slip), "adjusted_age_years"]
  out_df[i, "cent_deriv2_adj_age_years_slip"] <- tmp_50[which.min(tmp_50$cent_deriv2_slip), "adjusted_age_years"]
  
  out_df[i, "cent_deriv_adj_age_days_log_lmsz"] <- tmp_50[which.min(tmp_50$cent_deriv_lmsz), "adjusted_age_days_log"]
  out_df[i, "cent_deriv2_adj_age_days_log_lmsz"] <- tmp_50[which.min(tmp_50$cent_deriv2_lmsz), "adjusted_age_days_log"]
  out_df[i, "cent_deriv_adj_age_years_lmsz"] <- tmp_50[which.min(tmp_50$cent_deriv_lmsz), "adjusted_age_years"]
  out_df[i, "cent_deriv2_adj_age_years_lmsz"] <- tmp_50[which.min(tmp_50$cent_deriv2_lmsz), "adjusted_age_years"]
  
  age_at_peak_slip <- tmp_50[which.max(tmp_50$sex_average_slip), "adjusted_age_days_log"]
  out_df[i,"adj_age_days_log_peak50_slip"] <- age_at_peak_slip
  age_at_peak_lmsz <- tmp_50[which.max(tmp_50$sex_average_lmsz), "adjusted_age_days_log"]
  out_df[i,"adj_age_days_log_peak50_lmsz"] <- age_at_peak_lmsz
  
  # Female
  out_df[i, "cent_deriv_adj_age_days_log_slip_female"] <- tmp_50_f[which.min(tmp_50_f$cent_deriv_slip), "adjusted_age_days_log"]
  out_df[i, "cent_deriv2_adj_age_days_log_slip_female"] <- tmp_50_f[which.min(tmp_50_f$cent_deriv2_slip), "adjusted_age_days_log"]
  out_df[i, "cent_deriv_adj_age_years_slip_female"] <- tmp_50_f[which.min(tmp_50_f$cent_deriv_slip), "adjusted_age_years"]
  out_df[i, "cent_deriv2_adj_age_years_slip_female"] <- tmp_50_f[which.min(tmp_50_f$cent_deriv2_slip), "adjusted_age_years"]
  
  out_df[i, "cent_deriv_adj_age_days_log_lmsz_female"] <- tmp_50_f[which.min(tmp_50_f$cent_deriv_lmsz), "adjusted_age_days_log"]
  out_df[i, "cent_deriv2_adj_age_days_log_lmsz_female"] <- tmp_50_f[which.min(tmp_50_f$cent_deriv2_lmsz), "adjusted_age_days_log"]
  out_df[i, "cent_deriv_adj_age_years_lmsz_female"] <- tmp_50_f[which.min(tmp_50_f$cent_deriv_lmsz), "adjusted_age_years"]
  out_df[i, "cent_deriv2_adj_age_years_lmsz_female"] <- tmp_50_f[which.min(tmp_50_f$cent_deriv2_lmsz), "adjusted_age_years"]
  
  age_at_peak_slip_female <- tmp_50_f[which.max(tmp_50_f$SLIP_input_fans), "adjusted_age_days_log"]
  out_df[i,"adj_age_days_log_peak50_slip_female"] <- age_at_peak_slip_female
  age_at_peak_lmsz_female <- tmp_50_f[which.max(tmp_50_f$LMSz_fans), "adjusted_age_days_log"]
  out_df[i,"adj_age_days_log_peak50_lmsz_female"] <- age_at_peak_lmsz_female
  
  # Male
  out_df[i, "cent_deriv_adj_age_days_log_slip_male"] <- tmp_50_m[which.min(tmp_50_m$cent_deriv_slip), "adjusted_age_days_log"]
  out_df[i, "cent_deriv2_adj_age_days_log_slip_male"] <- tmp_50_m[which.min(tmp_50_m$cent_deriv2_slip), "adjusted_age_days_log"]
  out_df[i, "cent_deriv_adj_age_years_slip_male"] <- tmp_50_m[which.min(tmp_50_m$cent_deriv_slip), "adjusted_age_years"]
  out_df[i, "cent_deriv2_adj_age_years_slip_male"] <- tmp_50_m[which.min(tmp_50_m$cent_deriv2_slip), "adjusted_age_years"]
  
  out_df[i, "cent_deriv_adj_age_days_log_lmsz_male"] <- tmp_50_m[which.min(tmp_50_m$cent_deriv_lmsz), "adjusted_age_days_log"]
  out_df[i, "cent_deriv2_adj_age_days_log_lmsz_male"] <- tmp_50_m[which.min(tmp_50_m$cent_deriv2_lmsz), "adjusted_age_days_log"]
  out_df[i, "cent_deriv_adj_age_years_lmsz_male"] <- tmp_50_m[which.min(tmp_50_m$cent_deriv_lmsz), "adjusted_age_years"]
  out_df[i, "cent_deriv2_adj_age_years_lmsz_male"] <- tmp_50_m[which.min(tmp_50_m$cent_deriv2_lmsz), "adjusted_age_years"]
  
  age_at_peak_slip_male <- tmp_50_m[which.max(tmp_50_m$SLIP_input_fans), "adjusted_age_days_log"]
  out_df[i,"adj_age_days_log_peak50_slip_male"] <- age_at_peak_slip_male
  age_at_peak_lmsz_male <- tmp_50_m[which.max(tmp_50_m$LMSz_fans), "adjusted_age_days_log"]
  out_df[i,"adj_age_days_log_peak50_lmsz_male"] <- age_at_peak_lmsz_male
}

out_df$cent_curve <- 0.5
out_df <- out_df %>%
  mutate(adj_age_years_peak50_slip = exp(adj_age_days_log_peak50_slip)/365.25,
         adj_age_years_peak50_lmsz = exp(adj_age_days_log_peak50_lmsz)/365.25,
         adj_age_years_peak50_slip_female = exp(adj_age_days_log_peak50_slip_female)/365.25,
         adj_age_years_peak50_lmsz_female = exp(adj_age_days_log_peak50_lmsz_female)/365.25,
         adj_age_years_peak50_slip_male = exp(adj_age_days_log_peak50_slip_male)/365.25,
         adj_age_years_peak50_lmsz_male = exp(adj_age_days_log_peak50_lmsz_male)/365.25)

# add columns for spin testing and plotting
out_df$hemisphere <- ifelse(grepl("lh", out_df$pheno), "L", 
                            ifelse(grepl("rh", out_df$pheno), "R", NA))
out_df$domain <- ifelse(grepl("aparc_GrayVol", out_df$pheno), "vol", 
                        ifelse(grepl("aparc_SurfArea", out_df$pheno), "sa",
                               ifelse(grepl("aparc_ThickAvg", out_df$pheno), "th",
                                      ifelse(grepl("global", out_df$pheno), "global", 
                                             ifelse(grepl("aseg", out_df$pheno), "aseg",NA)))))

out_df$label <- sub("^(aparc_(GrayVol|SurfArea|ThickAvg)_(lh|rh)_)", "", out_df$pheno)

#Assign Ranks
## Sex Average
out_df[which(out_df$domain == "vol"),]$rank_peak50_slip <- out_df %>% filter(domain == "vol") %>% frank(adj_age_days_log_peak50_slip, ties.method = "dense")
out_df[which(out_df$domain == "sa"),]$rank_peak50_slip <- out_df %>% filter(domain == "sa") %>% frank(adj_age_days_log_peak50_slip, ties.method = "dense")
out_df[which(out_df$domain == "th"),]$rank_peak50_slip <- out_df %>% filter(domain == "th") %>% frank(adj_age_days_log_peak50_slip, ties.method = "dense")

out_df[which(out_df$domain == "vol"),]$rank_centderiv_slip <- out_df %>% filter(domain == "vol") %>% frank(cent_deriv_adj_age_days_log_slip, ties.method = "dense")
out_df[which(out_df$domain == "sa"),]$rank_centderiv_slip <- out_df %>% filter(domain == "sa") %>% frank(cent_deriv_adj_age_days_log_slip, ties.method = "dense")
out_df[which(out_df$domain == "th"),]$rank_centderiv <- out_df %>% filter(domain == "th") %>% frank(cent_deriv_adj_age_days_log_slip, ties.method = "dense")

out_df[which(out_df$domain == "vol"),]$rank_centderiv2_slip <- out_df %>% filter(domain == "vol") %>% frank(cent_deriv2_adj_age_days_log_slip, ties.method = "dense")
out_df[which(out_df$domain == "sa"),]$rank_centderiv2_slip <- out_df %>% filter(domain == "sa") %>% frank(cent_deriv2_adj_age_days_log_slip, ties.method = "dense")
out_df[which(out_df$domain == "th"),]$rank_centderiv2_slip <- out_df %>% filter(domain == "th") %>% frank(cent_deriv2_adj_age_days_log_slip, ties.method = "dense")


out_df[which(out_df$domain == "vol"),]$rank_peak50_lmsz <- out_df %>% filter(domain == "vol") %>% frank(adj_age_days_log_peak50_lmsz, ties.method = "dense")
out_df[which(out_df$domain == "sa"),]$rank_peak50_lmsz <- out_df %>% filter(domain == "sa") %>% frank(adj_age_days_log_peak50_lmsz, ties.method = "dense")
out_df[which(out_df$domain == "th"),]$rank_peak50_lmsz <- out_df %>% filter(domain == "th") %>% frank(adj_age_days_log_peak50_lmsz, ties.method = "dense")

out_df[which(out_df$domain == "vol"),]$rank_centderiv_lmsz <- out_df %>% filter(domain == "vol") %>% frank(cent_deriv_adj_age_days_log_lmsz, ties.method = "dense")
out_df[which(out_df$domain == "sa"),]$rank_centderiv_lmsz <- out_df %>% filter(domain == "sa") %>% frank(cent_deriv_adj_age_days_log_lmsz, ties.method = "dense")
out_df[which(out_df$domain == "th"),]$rank_centderiv <- out_df %>% filter(domain == "th") %>% frank(cent_deriv_adj_age_days_log_lmsz, ties.method = "dense")

out_df[which(out_df$domain == "vol"),]$rank_centderiv2_lmsz <- out_df %>% filter(domain == "vol") %>% frank(cent_deriv2_adj_age_days_log_lmsz, ties.method = "dense")
out_df[which(out_df$domain == "sa"),]$rank_centderiv2_lmsz <- out_df %>% filter(domain == "sa") %>% frank(cent_deriv2_adj_age_days_log_lmsz, ties.method = "dense")
out_df[which(out_df$domain == "th"),]$rank_centderiv2_lmsz <- out_df %>% filter(domain == "th") %>% frank(cent_deriv2_adj_age_days_log_lmsz, ties.method = "dense")

## Female
out_df[which(out_df$domain == "vol"),]$rank_peak50_slip_female <- out_df %>% filter(domain == "vol") %>% frank(adj_age_days_log_peak50_slip_female, ties.method = "dense")
out_df[which(out_df$domain == "sa"),]$rank_peak50_slip_female <- out_df %>% filter(domain == "sa") %>% frank(adj_age_days_log_peak50_slip_female, ties.method = "dense")
out_df[which(out_df$domain == "th"),]$rank_peak50_slip_female <- out_df %>% filter(domain == "th") %>% frank(adj_age_days_log_peak50_slip_female, ties.method = "dense")

out_df[which(out_df$domain == "vol"),]$rank_centderiv_slip_female <- out_df %>% filter(domain == "vol") %>% frank(cent_deriv_adj_age_days_log_slip_female, ties.method = "dense")
out_df[which(out_df$domain == "sa"),]$rank_centderiv_slip_female <- out_df %>% filter(domain == "sa") %>% frank(cent_deriv_adj_age_days_log_slip_female, ties.method = "dense")
out_df[which(out_df$domain == "th"),]$rank_centderiv <- out_df %>% filter(domain == "th") %>% frank(cent_deriv_adj_age_days_log_slip_female, ties.method = "dense")

out_df[which(out_df$domain == "vol"),]$rank_centderiv2_slip_female <- out_df %>% filter(domain == "vol") %>% frank(cent_deriv2_adj_age_days_log_slip_female, ties.method = "dense")
out_df[which(out_df$domain == "sa"),]$rank_centderiv2_slip_female <- out_df %>% filter(domain == "sa") %>% frank(cent_deriv2_adj_age_days_log_slip_female, ties.method = "dense")
out_df[which(out_df$domain == "th"),]$rank_centderiv2_slip_female <- out_df %>% filter(domain == "th") %>% frank(cent_deriv2_adj_age_days_log_slip_female, ties.method = "dense")


out_df[which(out_df$domain == "vol"),]$rank_peak50_lmsz_female <- out_df %>% filter(domain == "vol") %>% frank(adj_age_days_log_peak50_lmsz_female, ties.method = "dense")
out_df[which(out_df$domain == "sa"),]$rank_peak50_lmsz_female <- out_df %>% filter(domain == "sa") %>% frank(adj_age_days_log_peak50_lmsz_female, ties.method = "dense")
out_df[which(out_df$domain == "th"),]$rank_peak50_lmsz_female <- out_df %>% filter(domain == "th") %>% frank(adj_age_days_log_peak50_lmsz_female, ties.method = "dense")

out_df[which(out_df$domain == "vol"),]$rank_centderiv_lmsz_female <- out_df %>% filter(domain == "vol") %>% frank(cent_deriv_adj_age_days_log_lmsz_female, ties.method = "dense")
out_df[which(out_df$domain == "sa"),]$rank_centderiv_lmsz_female <- out_df %>% filter(domain == "sa") %>% frank(cent_deriv_adj_age_days_log_lmsz_female, ties.method = "dense")
out_df[which(out_df$domain == "th"),]$rank_centderiv <- out_df %>% filter(domain == "th") %>% frank(cent_deriv_adj_age_days_log_lmsz_female, ties.method = "dense")

out_df[which(out_df$domain == "vol"),]$rank_centderiv2_lmsz_female <- out_df %>% filter(domain == "vol") %>% frank(cent_deriv2_adj_age_days_log_lmsz_female, ties.method = "dense")
out_df[which(out_df$domain == "sa"),]$rank_centderiv2_lmsz_female <- out_df %>% filter(domain == "sa") %>% frank(cent_deriv2_adj_age_days_log_lmsz_female, ties.method = "dense")
out_df[which(out_df$domain == "th"),]$rank_centderiv2_lmsz_female <- out_df %>% filter(domain == "th") %>% frank(cent_deriv2_adj_age_days_log_lmsz_female, ties.method = "dense")

## Male
out_df[which(out_df$domain == "vol"),]$rank_peak50_slip_male <- out_df %>% filter(domain == "vol") %>% frank(adj_age_days_log_peak50_slip_male, ties.method = "dense")
out_df[which(out_df$domain == "sa"),]$rank_peak50_slip_male <- out_df %>% filter(domain == "sa") %>% frank(adj_age_days_log_peak50_slip_male, ties.method = "dense")
out_df[which(out_df$domain == "th"),]$rank_peak50_slip_male <- out_df %>% filter(domain == "th") %>% frank(adj_age_days_log_peak50_slip_male, ties.method = "dense")

out_df[which(out_df$domain == "vol"),]$rank_centderiv_slip_male <- out_df %>% filter(domain == "vol") %>% frank(cent_deriv_adj_age_days_log_slip_male, ties.method = "dense")
out_df[which(out_df$domain == "sa"),]$rank_centderiv_slip_male <- out_df %>% filter(domain == "sa") %>% frank(cent_deriv_adj_age_days_log_slip_male, ties.method = "dense")
out_df[which(out_df$domain == "th"),]$rank_centderiv <- out_df %>% filter(domain == "th") %>% frank(cent_deriv_adj_age_days_log_slip_male, ties.method = "dense")

out_df[which(out_df$domain == "vol"),]$rank_centderiv2_slip_male <- out_df %>% filter(domain == "vol") %>% frank(cent_deriv2_adj_age_days_log_slip_male, ties.method = "dense")
out_df[which(out_df$domain == "sa"),]$rank_centderiv2_slip_male <- out_df %>% filter(domain == "sa") %>% frank(cent_deriv2_adj_age_days_log_slip_male, ties.method = "dense")
out_df[which(out_df$domain == "th"),]$rank_centderiv2_slip_male <- out_df %>% filter(domain == "th") %>% frank(cent_deriv2_adj_age_days_log_slip_male, ties.method = "dense")


out_df[which(out_df$domain == "vol"),]$rank_peak50_lmsz_male <- out_df %>% filter(domain == "vol") %>% frank(adj_age_days_log_peak50_lmsz_male, ties.method = "dense")
out_df[which(out_df$domain == "sa"),]$rank_peak50_lmsz_male <- out_df %>% filter(domain == "sa") %>% frank(adj_age_days_log_peak50_lmsz_male, ties.method = "dense")
out_df[which(out_df$domain == "th"),]$rank_peak50_lmsz_male <- out_df %>% filter(domain == "th") %>% frank(adj_age_days_log_peak50_lmsz_male, ties.method = "dense")

out_df[which(out_df$domain == "vol"),]$rank_centderiv_lmsz_male <- out_df %>% filter(domain == "vol") %>% frank(cent_deriv_adj_age_days_log_lmsz_male, ties.method = "dense")
out_df[which(out_df$domain == "sa"),]$rank_centderiv_lmsz_male <- out_df %>% filter(domain == "sa") %>% frank(cent_deriv_adj_age_days_log_lmsz_male, ties.method = "dense")
out_df[which(out_df$domain == "th"),]$rank_centderiv <- out_df %>% filter(domain == "th") %>% frank(cent_deriv_adj_age_days_log_lmsz_male, ties.method = "dense")

out_df[which(out_df$domain == "vol"),]$rank_centderiv2_lmsz_male <- out_df %>% filter(domain == "vol") %>% frank(cent_deriv2_adj_age_days_log_lmsz_male, ties.method = "dense")
out_df[which(out_df$domain == "sa"),]$rank_centderiv2_lmsz_male <- out_df %>% filter(domain == "sa") %>% frank(cent_deriv2_adj_age_days_log_lmsz_male, ties.method = "dense")
out_df[which(out_df$domain == "th"),]$rank_centderiv2_lmsz_male <- out_df %>% filter(domain == "th") %>% frank(cent_deriv2_adj_age_days_log_lmsz_male, ties.method = "dense")

filename <- glue(("{csv_save_path}slip_lmsz_centile_peaks.csv"))
write.csv(out_df, filename)

x_breaks <- log(c(280, 365.25*1+280, 365.25*2+280, 
                  365.25*6+280, 365.25*21+280))
x_labs <- c("0y","1y","2y","6y","21y")

overlap_plots <- list()

brain_plot_df <- data.frame(
  pheno = full_vars,
  value = rep(NA, length(full_vars)),
  label = sub("(aparc_(GrayVol|SurfArea|ThickAvg)|aseg)_", "", full_vars),
  hemi = rep(NA, length(full_vars))
)

rows_edit <- grepl("aseg(?!_CC)", brain_plot_df$pheno, perl = TRUE)

brain_plot_df$label <- brain_plot_df$label %>% 
  sub("X", "x", .) %>%
  sub("rd_Ventricle", "rd_ventricle", .) %>%
  sub("th_Ventricle", "th_ventricle", .) %>% 
  sub("Thalamus", "Thalamus-Proper", .)

brain_plot_df$label[rows_edit] <- gsub("_", "-", brain_plot_df$label[rows_edit])

brain_plot_df$hemi <- ifelse(grepl("lh", brain_plot_df$pheno), "left",
                             ifelse(grepl("rh", brain_plot_df$pheno), "right", NA))
  
pheno_list <- lmsz_vars
  
  for (pheno_of_int in pheno_list) {
    #df_plot$pheno <- df_plot[[glue("{pheno_of_int}.raw")]]
    #df_plot$pheno[is.na(df_plot[[glue("{pheno_of_int}.lmsz")]])] <- NA
    if(grepl("global", pheno_of_int)){
      ylab_text = ""
      title_text = sub("global_", "", pheno_of_int)
    } else if(grepl("GrayVol", pheno_of_int))  {
      ylab_text = "Volume"
      title_text = sub("aparc_GrayVol_", "", pheno_of_int)
    } else if(grepl("SurfArea", pheno_of_int))  {
      ylab_text = "Surface Area"
      title_text = sub("aparc_SurfArea_", "", pheno_of_int)
    } else if(grepl("ThickAvg", pheno_of_int))  {
      ylab_text = "Thickness"
      title_text = sub("aparc_ThickAvg_", "", pheno_of_int)
    } else if(grepl("aseg", pheno_of_int))  {
      ylab_text = "Aseg/Subcortical"
      title_text = sub("aseg_", "", pheno_of_int)
    } else{
      ylab_text = pheno_of_int
      title_text = pheno_of_int
    }
    
    fans_plot <- fans %>% 
      filter(pheno == pheno_of_int,
             !(centile_cat %in% c("cent_0.25", "cent_0.75"))) %>%
      mutate(adjusted_age_years = exp(adjusted_age_days_log)/365.25)
    
    
    out_df_plot <- out_df %>%
      filter(pheno == pheno_of_int)
    
    if (grepl("Thick",pheno_of_int)) {
      fans_plot <- fans_plot %>% filter(adjusted_age_days_log >= log(2*365.25+280))
    }
    
    # peak at 50th
    i_peak50_lmsz_m = fans_plot[fans_plot$adjusted_age_days_log == out_df_plot$adj_age_days_log_peak50_lmsz_male,] %>% 
      filter(sex == "Male" & centile_cat == "cent_0.5")
    i_peak50_lmsz_f = fans_plot[fans_plot$adjusted_age_days_log == out_df_plot$adj_age_days_log_peak50_lmsz_female,] %>% 
      filter(sex == "Female" & centile_cat == "cent_0.5")
    i_peak50_lmsz_avg = fans_plot[fans_plot$adjusted_age_days_log == out_df_plot$adj_age_days_log_peak50_lmsz,] %>% 
      filter(centile_cat == "cent_0.5") %>%
      mutate(sex_average_lmsz = mean(LMSz_fans))
    
    i_peak50_slip_m = fans_plot[fans_plot$adjusted_age_days_log == out_df_plot$adj_age_days_log_peak50_slip_male,] %>% 
      filter(sex == "Male" & centile_cat == "cent_0.5")
    i_peak50_slip_f = fans_plot[fans_plot$adjusted_age_days_log == out_df_plot$adj_age_days_log_peak50_slip_female,] %>% 
      filter(sex == "Female" & centile_cat == "cent_0.5")
    i_peak50_slip_avg = fans_plot[fans_plot$adjusted_age_days_log == out_df_plot$adj_age_days_log_peak50_slip,] %>% 
      filter(centile_cat == "cent_0.5") %>%
      mutate(sex_average_slip = mean(SLIP_input_fans))
    
    # min d1
    i_d1_lmsz_m = fans_plot[fans_plot$adjusted_age_days_log == out_df_plot$cent_deriv_adj_age_days_log_lmsz_male,] %>% 
      filter(sex == "Male" & centile_cat == "cent_0.5")
    i_d1_lmsz_f = fans_plot[fans_plot$adjusted_age_days_log == out_df_plot$cent_deriv_adj_age_days_log_lmsz_female,] %>% 
      filter(sex == "Female" & centile_cat == "cent_0.5")
    i_d1_lmsz_avg = fans_plot[fans_plot$adjusted_age_days_log == out_df_plot$cent_deriv_adj_age_days_log_lmsz,] %>% 
      filter(centile_cat == "cent_0.5") %>%
      mutate(sex_average_lmsz = mean(LMSz_fans))
    
    i_d1_slip_m = fans_plot[fans_plot$adjusted_age_days_log == out_df_plot$cent_deriv_adj_age_days_log_slip_male,] %>% 
      filter(sex == "Male" & centile_cat == "cent_0.5")
    i_d1_slip_f = fans_plot[fans_plot$adjusted_age_days_log == out_df_plot$cent_deriv_adj_age_days_log_slip_female,] %>% 
      filter(sex == "Female" & centile_cat == "cent_0.5")
    i_d1_slip_avg = fans_plot[fans_plot$adjusted_age_days_log == out_df_plot$cent_deriv_adj_age_days_log_slip,] %>% 
      filter(centile_cat == "cent_0.5") %>%
      mutate(sex_average_slip = mean(SLIP_input_fans))
    
    # min d2
    i_d2_lmsz_m = fans_plot[fans_plot$adjusted_age_days_log == out_df_plot$cent_deriv2_adj_age_days_log_lmsz_male,] %>% 
      filter(sex == "Male" & centile_cat == "cent_0.5")
    i_d2_lmsz_f = fans_plot[fans_plot$adjusted_age_days_log == out_df_plot$cent_deriv2_adj_age_days_log_lmsz_female,] %>% 
      filter(sex == "Female" & centile_cat == "cent_0.5")
    i_d2_lmsz_avg = fans_plot[fans_plot$adjusted_age_days_log == out_df_plot$cent_deriv2_adj_age_days_log_lmsz,] %>% 
      filter(centile_cat == "cent_0.5") %>%
      mutate(sex_average_lmsz = mean(LMSz_fans))
    
    i_d2_slip_m = fans_plot[fans_plot$adjusted_age_days_log == out_df_plot$cent_deriv2_adj_age_days_log_slip_male,] %>% 
      filter(sex == "Male" & centile_cat == "cent_0.5")
    i_d2_slip_f = fans_plot[fans_plot$adjusted_age_days_log == out_df_plot$cent_deriv2_adj_age_days_log_slip_female,] %>% 
      filter(sex == "Female" & centile_cat == "cent_0.5")
    i_d2_slip_avg = fans_plot[fans_plot$adjusted_age_days_log == out_df_plot$cent_deriv2_adj_age_days_log_slip,] %>% 
      filter(centile_cat == "cent_0.5") %>%
      mutate(sex_average_slip = mean(SLIP_input_fans))
      
    p_male <- ggplot(fans_plot %>% filter(sex == "Male")) +
      geom_line(aes(x = adjusted_age_days_log, y = SLIP_input_fans,
                    group = centile_cat,
                    color = "SLIP",
                    linewidth = line_size,
                    linetype = line_type)) +
      geom_line(aes(x = adjusted_age_days_log, y = LMSz_fans,
                    group = centile_cat,
                    color = "LMSz",
                    linewidth = line_size,
                    linetype = line_type)) +
      geom_point(data = i_peak50_slip_m, 
                 aes(adjusted_age_days_log, SLIP_input_fans,  color = "SLIP", shape = "50th centile peak"),
                 size = 2.5) +
      geom_point(data = i_d1_slip_m, 
                 aes(adjusted_age_days_log, SLIP_input_fans,  color = "SLIP", shape = "First Deriv Min"),
                 size = 2.5) +
      geom_point(data = i_d2_slip_m, 
                 aes(adjusted_age_days_log, SLIP_input_fans,  color = "SLIP", shape = "Second Deriv Min"),
                 size = 2.5) +
      geom_point(data = i_peak50_lmsz_m, 
                 aes(adjusted_age_days_log, LMSz_fans,  color = "LMSz", shape = "50th centile peak"),
                 size = 2.5) +
      geom_point(data = i_d1_lmsz_m, 
                 aes(adjusted_age_days_log, LMSz_fans,  color = "LMSz", shape = "First Deriv Min"),
                 size = 2.5) +
      geom_point(data = i_d2_lmsz_m, 
                 aes(adjusted_age_days_log, LMSz_fans,  color = "LMSz", shape = "Second Deriv Min"),
                 size = 2.5) +
      scale_linewidth_identity() +
      scale_linetype_identity() +
      scale_colour_manual(
        values = c("LMSz" = "#566fa3", 
                   "SLIP" = "grey")) +
      scale_shape_manual(
        values = c("50th centile peak" = 8,
                   "First Deriv Min"   = 18,
                   "Second Deriv Min"  = 16)) +
      labs(y = "Male", x = "", colour = "Model", shape = "") +
      scale_x_continuous(breaks = x_breaks, labels = x_labs) +
      guides(colour = "none") +
      theme_classic() +
      theme(plot.title = element_text(size = 12),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 8))
    
    p_female <- ggplot(fans_plot %>% filter(sex == "Female")) +
      geom_line(aes(x = adjusted_age_days_log, y = SLIP_input_fans,
                    group = centile_cat,
                    color = "SLIP",
                    linewidth = line_size,
                    linetype = line_type)) +
      geom_line(aes(x = adjusted_age_days_log, y = LMSz_fans,
                    group = centile_cat,
                    color = "LMSz",
                    linewidth = line_size,
                    linetype = line_type)) +
      geom_point(data = i_peak50_slip_f, 
                 aes(adjusted_age_days_log, SLIP_input_fans,  color = "SLIP", shape = "50th centile peak"),
                 size = 2.5) +
      geom_point(data = i_d1_slip_f, 
                 aes(adjusted_age_days_log, SLIP_input_fans,  color = "SLIP", shape = "First Deriv Min"),
                 size = 2.5) +
      geom_point(data = i_d2_slip_f, 
                 aes(adjusted_age_days_log, SLIP_input_fans,  color = "SLIP", shape = "Second Deriv Min"),
                 size = 2.5) +
      geom_point(data = i_peak50_lmsz_f, 
                 aes(adjusted_age_days_log, LMSz_fans,  color = "LMSz", shape = "50th centile peak"),
                 size = 2.5) +
      geom_point(data = i_d1_lmsz_f, 
                 aes(adjusted_age_days_log, LMSz_fans,  color = "LMSz", shape = "First Deriv Min"),
                 size = 2.5) +
      geom_point(data = i_d2_lmsz_f, 
                 aes(adjusted_age_days_log, LMSz_fans,  color = "LMSz", shape = "Second Deriv Min"),
                 size = 2.5) +
      scale_linewidth_identity() +
      scale_linetype_identity() +
      scale_colour_manual(
        values = c("LMSz" = "#E396B7", 
                   "SLIP" = "grey")) +
      scale_shape_manual(
        values = c("50th centile peak" = 8,
                   "First Deriv Min"   = 18,
                   "Second Deriv Min"  = 16)) +
      labs(y = "Female", x = "", colour = "Model", shape = "") +
      scale_x_continuous(breaks = x_breaks, labels = x_labs) +
      guides(colour = "none", shape = "none") +
      theme_classic() +
      theme(plot.title = element_text(size = 12),
            axis.title = element_text(size = 10),
            axis.text = element_text(size = 8))
  
  # take sex average centile fans
    fans_plot_avg <- fans_plot %>% 
      group_by(model, adjusted_age_days_log, 
               centile_cat, line_size, line_type) %>%
      summarize(across(contains("fans"), ~mean(.x, na.rm = T)))
    
  p_sexavg <- ggplot(fans_plot_avg) +
    geom_line(aes(x = adjusted_age_days_log, y = SLIP_input_fans,
                  group = centile_cat,
                  color = "SLIP",
                  linewidth = line_size,
                  linetype = line_type)) +
    geom_line(aes(x = adjusted_age_days_log, y = LMSz_fans,
                  group = centile_cat,
                  color = "LMSz",
                  linewidth = line_size,
                  linetype = line_type)) +
    geom_point(data = i_peak50_slip_avg, 
               aes(adjusted_age_days_log, sex_average_slip,  color = "SLIP", shape = "50th centile peak"),
               size = 2.5) +
    geom_point(data = i_d1_slip_avg, 
               aes(adjusted_age_days_log, sex_average_slip,  color = "SLIP", shape = "First Deriv Min"),
               size = 2.5) +
    geom_point(data = i_d2_slip_avg, 
               aes(adjusted_age_days_log, sex_average_slip,  color = "SLIP", shape = "Second Deriv Min"),
               size = 2.5) +
    geom_point(data = i_peak50_lmsz_avg, 
               aes(adjusted_age_days_log, sex_average_lmsz,  color = "LMSz", shape = "50th centile peak"),
               size = 2.5) +
    geom_point(data = i_d1_lmsz_avg, 
               aes(adjusted_age_days_log, sex_average_lmsz,  color = "LMSz", shape = "First Deriv Min"),
               size = 2.5) +
    geom_point(data = i_d2_lmsz_avg, 
               aes(adjusted_age_days_log, sex_average_lmsz,  color = "LMSz", shape = "Second Deriv Min"),
               size = 2.5) +
    scale_linewidth_identity() +
    scale_linetype_identity() +
    scale_colour_manual(
      values = c("LMSz" = "#5FA052", 
                 "SLIP" = "grey")) +
    scale_shape_manual(
      values = c("50th centile peak" = 8,
                 "First Deriv Min"   = 18,
                 "Second Deriv Min"  = 16)) +
    labs(y = "Sex Avg", x = "", colour = "Model", shape = "") +
    scale_x_continuous(breaks = x_breaks, labels = x_labs) +
    guides(colour = "none", shape  = "none") +
    theme_classic() +
    theme(plot.title = element_text(size = 12),
          axis.title = element_text(size = 10),
          axis.text = element_text(size = 8))


  # Combine all three plots
  all_plots <- (p_female / p_male / p_sexavg) +
    plot_annotation(title = title_text, subtitle = ylab_text) &
    theme(
      plot.title.position = "plot",
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(
        angle = 0,
        vjust = 0.5,
        hjust = 0.5
      ))

  # Add brain plot
  limits = c(0,1)
  high_color = "#0F766E"
  plot_df <- brain_plot_df %>% select(pheno, value, label)
  pheno_base <- sub("aparc_(GrayVol|SurfArea|ThickAvg)_", "", pheno_of_int)
  plot_df[grepl(pheno_base, plot_df$pheno),"value"] <- 1
  
  if(grepl("aparc", pheno_of_int)) {
    plot_data_dkt <- merge(plot_df, ggseg::dk, by = "label")
    hemi_plot <- plot_data_dkt[which.max(plot_data_dkt$value), "hemi"]
    side_plot <- plot_data_dkt[which.max(plot_data_dkt$value), "side"]
    
    brain_plot <- ggseg(plot_data_dkt, 
                        atlas = dk, 
                        colour = "black",
                        size = .2,
                        hemisphere = hemi_plot,
                        view = side_plot,
                        mapping = aes(fill = value)) +
      scale_fill_gradient2(high = high_color, limits = limits) +
      theme_void()+
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "none")
    
    all_plots <- all_plots +inset_element(brain_plot, left = 0.95, bottom = 0.2, right = 1.5, top = 1, align_to = "plot")
    
    
  } else if (grepl("aseg", pheno_of_int)){ 
    plot_data_aseg <- merge(plot_df, ggseg::aseg, by = "label")
    brain_plot <- ggplot(plot_data_aseg) +
      geom_brain(atlas = aseg,
                 aes(fill = value)) +
      scale_fill_gradient2(high = high_color, limits = limits) +
      theme_void()+
      theme(axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "none")
    
    all_plots <- all_plots +  inset_element(brain_plot, left = 1.02, bottom = 0.2, right = 1.4, top = 1, align_to = "plot")
    
  }
  
overlap_plots[[pheno_of_int]] <- all_plots

}

curves_filename <- glue(("{plot_save_path}slip_lmsz_centile_peaks_plots.rds"))
saveRDS(overlap_plots, curves_filename)
