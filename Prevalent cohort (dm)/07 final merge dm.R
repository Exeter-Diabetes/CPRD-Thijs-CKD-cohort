############################################################################################

#Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "diabetes-jun2024",cprdConf = "C:/Users/tj358/OneDrive - University of Exeter/CPRD/aurum.yaml")
codesets = cprd$codesets()
codes_2024 = codesets$getAllCodeSetVersion(v = "01/06/2024")

analysis_prefix = "ckd"


############################################################################################

## Cohort and patient characteristics
analysis = cprd$analysis("all")
ckd_cohort <- ckd_cohort %>% analysis$cached("diabetes_ckd_cohort")
diabetes_cohort <- diabetes_cohort %>% analysis$cached("diabetes_cohort")
death_causes <- death_causes %>% analysis$cached("death_causes")
townsend_score <- townsend_score %>% analysis$cached("townsend_score")

## Get index date

# get dates at 6 month intervals
dates <- seq(from = as.Date("2019-03-01"),
             to   = as.Date("2024-03-01"),
             by   = "6 months")

date_strings <- format(dates, "%Y-%m-%d")

# create empty dataframe for counts of total population / subset with CKD
counts <- data.frame()

for (d in date_strings) {
  
  index_date <- as.Date(d)
  print(d)
  
  analysis = cprd$analysis(paste0(analysis_prefix, "_", d))
  
  
  ## Biomarkers plus CKD stage
  ckd_stages <- ckd_stages %>% analysis$cached("ckd_stages")
  baseline_biomarkers <- baseline_biomarkers %>% analysis$cached("baseline_biomarkers")
  
  ## Comorbidities
  comorbidities <- comorbidities %>% analysis$cached("comorbidities")
  
  ## Smoking status
  smoking <- smoking %>% analysis$cached("smoking")
  
  ## Medications
  medications <- medications %>% analysis$cached("medications")
  
  
  ############################################################################################
  
  # Define prevalent cohort and add in variables from other tables plus age and diabetes duration at index date and QRISK2 and QDiabetes-HF
  ## Prevalent cohort: registered on 01/02/2020 and with diagnosis at/before then and with linked HES records (and n_patid_hes<=20).
  
  cohort_ids <- ckd_cohort %>%
    filter(first_ckd_date<=index_date & 
             dm_diag_date_all<=index_date & 
             regstartdate<=index_date & 
             gp_end_date>=index_date & 
             (is.na(death_date) | death_date>=index_date)) %>%
    select(patid) %>%
    analysis$cached("cohort_ids", unique_indexes="patid")
  
  # get counts of all patients + patients with CKD at index date
  total_count <- diabetes_cohort %>% 
    filter(diabetes_type == "type 2" & dm_diag_date_all<=index_date & regstartdate<=index_date & gp_end_date>=index_date & (is.na(death_date) | death_date>=index_date)) %>%
    select(patid) %>%
    count() %>% 
    collect() %>% 
    pull(n)
  total_valid_egfr_count <- diabetes_cohort %>% 
    filter(diabetes_type == "type 2" & dm_diag_date_all<=index_date & regstartdate<=index_date & gp_end_date>=index_date & (is.na(death_date) | death_date>=index_date)) %>%
    left_join(baseline_biomarkers, by = "patid") %>%
    filter(!is.na(preegfr)) %>%
    select(patid) %>%
    count() %>% 
    collect() %>% 
    pull(n)
  total_valid_egfr_uacr_count <- diabetes_cohort %>% 
    filter(diabetes_type == "type 2" & dm_diag_date_all<=index_date & regstartdate<=index_date & gp_end_date>=index_date & (is.na(death_date) | death_date>=index_date)) %>%
    left_join(baseline_biomarkers, by = "patid") %>%
    filter(!is.na(preegfr) & (!is.na(preacr) | !is.na(preacr_from_separate))) %>%
    select(patid) %>%
    count() %>% 
    collect() %>% 
    pull(n)
  ckd_count <- cohort_ids %>% 
    count() %>% 
    collect() %>% 
    pull(n)
  percentage <- round(ckd_count / total_count * 100, 1)
  count_at_date <- data.frame(ckd_count = ckd_count, 
                              total_count = total_count, 
                              total_valid_egfr_count = total_valid_egfr_count, 
                              total_valid_egfr_uacr_count = total_valid_egfr_uacr_count, 
                              percentage = percentage, 
                              date = d)
  counts <- rbind(counts, count_at_date)
  rm(count_at_date)
  print(paste0("CKD prevalence: ", percentage, "% (", ckd_count, " out of ", total_count, " at ", d, ")"))
  
  final_merge <- cohort_ids %>%
    left_join(ckd_cohort, by="patid") %>%
    left_join(ckd_stages, by="patid") %>%
    left_join(baseline_biomarkers, by="patid") %>%
    left_join(comorbidities, by="patid") %>%
    left_join(smoking, by="patid") %>%
    left_join(medications, by="patid") %>%
    left_join(townsend_score %>% select(patid, tds_2011), by = "patid") %>%
    left_join(death_causes, by = "patid") %>%
    mutate(index_date_age=datediff(index_date, dob)/365.25,
           index_date_ckd_dur_all=datediff(index_date, first_ckd_date)/365.25,
           index_date = index_date) %>%
    relocate(c(index_date_age, index_date_ckd_dur_all), .before=gender) %>%
    analysis$cached("final_merge", unique_indexes="patid")
  
  ###############################################################
  
  ## Make separate table with additional variables for QRISK2
  
  qscore_vars <- final_merge %>%
    mutate(precholhdl=pretotalcholesterol/prehdl,
           ckd45=!is.na(preckdstage) & (preckdstage=="stage_4" | preckdstage=="stage_5"),
           cvd=pre_index_date_myocardialinfarction==1 | pre_index_date_angina==1 | pre_index_date_stroke==1,
           sex=ifelse(gender==1, "male", ifelse(gender==2, "female", "NA")),
           # dm_duration_cat=0L,
           
           dm_duration_cat=ifelse(dm_dur_all<=1, 0L,
                                  ifelse(dm_dur_all<4, 1L,
                                         ifelse(dm_dur_all<7, 2L,
                                                ifelse(dm_dur_all<11, 3L, 4L)))),
           
           earliest_bp_med=pmin(
             ifelse(is.na(pre_index_date_earliest_ace_inhibitors),as.Date("2050-01-01"),pre_index_date_earliest_ace_inhibitors),
             ifelse(is.na(pre_index_date_earliest_beta_blockers),as.Date("2050-01-01"),pre_index_date_earliest_beta_blockers),
             ifelse(is.na(pre_index_date_earliest_ca_channel_blockers),as.Date("2050-01-01"),pre_index_date_earliest_ca_channel_blockers),
             ifelse(is.na(pre_index_date_earliest_thiazide_diuretics),as.Date("2050-01-01"),pre_index_date_earliest_thiazide_diuretics),
             na.rm=TRUE
           ),
           latest_bp_med=pmax(
             ifelse(is.na(pre_index_date_latest_ace_inhibitors),as.Date("1900-01-01"),pre_index_date_latest_ace_inhibitors),
             ifelse(is.na(pre_index_date_latest_beta_blockers),as.Date("1900-01-01"),pre_index_date_latest_beta_blockers),
             ifelse(is.na(pre_index_date_latest_ca_channel_blockers),as.Date("1900-01-01"),pre_index_date_earliest_ca_channel_blockers),
             ifelse(is.na(pre_index_date_latest_thiazide_diuretics),as.Date("1900-01-01"),pre_index_date_latest_thiazide_diuretics),
             na.rm=TRUE
           ),
           bp_meds=ifelse(earliest_bp_med!=as.Date("2050-01-01") & latest_bp_med!=as.Date("1900-01-01") & datediff(index_date, latest_bp_med)<=28 & earliest_bp_med!=latest_bp_med, 1L, 0L),
           
           type1=0L,
           type2=0L,
           surv_5yr=5L,
           surv_10yr=10L) %>%
    
    select(patid, sex, index_date_age, ethnicity_qrisk2, qrisk2_smoking_cat, dm_duration_cat, bp_meds, type1, type2, cvd, ckd45, pre_index_date_fh_premature_cvd, pre_index_date_af, pre_index_date_rheumatoidarthritis, prehba1c, precholhdl, presbp, prebmi, tds_2011, surv_5yr, surv_10yr) %>%
    
    analysis$cached( "qrisk_vars", indexes=c("patid"))
  
  
  
  ## Calculate 10 year QRISK2 scores
  
  
  ## Remove QRISK2 score for those with biomarker values outside of range:
  ### CholHDL: missing or 1-12
  ### SBP: missing or 70-210
  ### Age: 25-84
  ### Also exclude if BMI<20 as v. different from development cohort
  
  qscores <- qscore_vars %>%
    
    mutate(sex2=ifelse(sex=="male", "male", ifelse(sex=="female", "female", NA))) %>%
    
    calculate_qrisk2(sex=sex2, age=index_date_age, ethrisk=ethnicity_qrisk2, smoking=qrisk2_smoking_cat, type1=type1, type2=type2, fh_cvd=pre_index_date_fh_premature_cvd, renal=ckd45, af=pre_index_date_af, rheumatoid_arth=pre_index_date_rheumatoidarthritis, cholhdl=precholhdl, sbp=presbp, bmi=prebmi, bp_med=bp_meds, town=tds_2011, surv=surv_10yr) %>%
    
    mutate(
      qrisk2_10yr_score=ifelse((is.na(precholhdl) | (precholhdl>=1 & precholhdl<=12)) &
                                 (is.na(presbp) | (presbp>=70 & presbp<=210)) &
                                 index_date_age>=25 & index_date_age<=84 &
                                 (is.na(prebmi) | prebmi>=20), qrisk2_score, NA),
      
      qrisk2_lin_predictor=ifelse((is.na(precholhdl) | (precholhdl>=1 & precholhdl<=12)) &
                                    (is.na(presbp) | (presbp>=70 & presbp<=210)) &
                                    index_date_age>=25 & index_date_age<=84 &
                                    (is.na(prebmi) | prebmi>=20), qrisk2_lin_predictor, NA)) %>%
    
    select(patid, qrisk2_10yr_score, qrisk2_lin_predictor) %>%
    
    analysis$cached("qrisk_scores", indexes=c("patid"))
  
  ## join with main dataset
  
  final_merge <- final_merge %>% left_join(qscores, by = "patid")
  ############################################################################################
  
  # Export to R data object
  ## Convert integer64 datatypes to double
  
  prev_cohort <- collect(final_merge %>% mutate(patid=as.character(patid)))
  
  is.integer64 <- function(x){
    class(x)=="integer64"
  }
  
  prev_cohort <- prev_cohort %>%
    mutate_if(is.integer64, as.integer) %>%
    mutate(index_date = as.Date(d))
  
  # Create a valid name (no dashes)
  df_name <- paste0("prev_", gsub("-", "_", d), "_dm")
  
  # Assign name
  assign(df_name, prev_cohort, envir = .GlobalEnv)
  
  today <- format(Sys.Date(), "%Y%m%d")
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Raw data/")
  save(list = df_name, file=paste0(today, "_prev_ckd_cohort_dm_", d, ".Rda"))
  
  
  rm(medications)
  rm(baseline_biomarkers)
  rm(comorbidities)
  rm(ckd_stages)
  rm(smoking)
  rm(cohort_ids)
  rm(final_merge)
}

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Raw data/")
save(counts, file=paste0(today, "_ckd_counts_dm.Rda"))