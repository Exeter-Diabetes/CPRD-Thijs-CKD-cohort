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


## Get index date

# get dates at 6 month intervals
dates <- seq(from = as.Date("2024-03-01"),
             to   = as.Date("2024-03-01"),
             by   = "6 months")

date_strings <- format(dates, "%Y-%m-%d")


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
  
  cohort_ids %>% count()
  
  final_merge <- cohort_ids %>%
    left_join(ckd_cohort, by="patid") %>%
    left_join(ckd_stages, by="patid") %>%
    left_join(baseline_biomarkers, by="patid") %>%
    left_join(comorbidities, by="patid") %>%
    left_join(smoking, by="patid") %>%
    left_join(medications, by="patid") %>%
    mutate(index_date_age=datediff(index_date, dob)/365.25,
           index_date_ckd_dur_all=datediff(index_date, first_ckd_date)/365.25,
           index_date = index_date) %>%
    relocate(c(index_date_age, index_date_ckd_dur_all), .before=gender) %>%
    analysis$cached("final_merge", unique_indexes="patid")
  
  
  # ############################################################################################
  # 
  # # Add in 5 year QRISK2 score
  # 
  # ## Make separate table with additional variables for QRISK2 
  # 
  # qscore_vars <- final_merge %>%
  #   mutate(precholhdl=pretotalcholesterol/prehdl,
  #          ckd45=!is.na(preckdstage) & (preckdstage=="stage_4" | preckdstage=="stage_5"),
  #          cvd=pre_index_date_myocardialinfarction==1 | pre_index_date_angina==1 | pre_index_date_stroke==1,
  #          sex=ifelse(gender==1, "male", ifelse(gender==2, "female", "NA")),
  #          dm_duration_cat=0L,
  #          
  #          earliest_bp_med=pmin(
  #            ifelse(is.na(pre_index_date_earliest_ace_inhibitors),as.Date("2050-01-01"),pre_index_date_earliest_ace_inhibitors),
  #            ifelse(is.na(pre_index_date_earliest_beta_blockers),as.Date("2050-01-01"),pre_index_date_earliest_beta_blockers),
  #            ifelse(is.na(pre_index_date_earliest_calcium_channel_blockers),as.Date("2050-01-01"),pre_index_date_earliest_calcium_channel_blockers),
  #            ifelse(is.na(pre_index_date_earliest_thiazide_diuretics),as.Date("2050-01-01"),pre_index_date_earliest_thiazide_diuretics),
  #            na.rm=TRUE
  #          ),
  #          latest_bp_med=pmax(
  #            ifelse(is.na(pre_index_date_latest_ace_inhibitors),as.Date("1900-01-01"),pre_index_date_latest_ace_inhibitors),
  #            ifelse(is.na(pre_index_date_latest_beta_blockers),as.Date("1900-01-01"),pre_index_date_latest_beta_blockers),
  #            ifelse(is.na(pre_index_date_latest_calcium_channel_blockers),as.Date("1900-01-01"),pre_index_date_latest_calcium_channel_blockers),
  #            ifelse(is.na(pre_index_date_latest_thiazide_diuretics),as.Date("1900-01-01"),pre_index_date_latest_thiazide_diuretics),
  #            na.rm=TRUE
  #          ),
  #          bp_meds=ifelse(earliest_bp_med!=as.Date("2050-01-01") & latest_bp_med!=as.Date("1900-01-01") & datediff(index_date, latest_bp_med)<=28 & earliest_bp_med!=latest_bp_med, 1L, 0L),
  #          
  #          type1=0L,
  #          type2=0L,
  #          surv_5yr=5L,
  #          surv_10yr=10L) %>%
  #   
  #   select(patid, sex, index_date_age, ethnicity_qrisk2, qrisk2_smoking_cat, dm_duration_cat, bp_meds, type1, type2, cvd, ckd45, pre_index_date_fh_premature_cvd, pre_index_date_af, pre_index_date_rheumatoidarthritis, prehba1c2yrs, precholhdl, presbp, prebmi, tds_2011, surv_5yr, surv_10yr) %>%
  #   
  #   analysis$cached("final_merge_im_q1", unique_indexes="patid")
  # 
  # 
  # 
  # ## Calculate 5 year and 10 year QRISK2 scores
  # ### For some reason it doesn't like collation of sex variable unless remake it
  # 
  # 
  # ## Remove QRISK2 score for those with biomarker values outside of range:
  # ### CholHDL: missing or 1-12
  # ### SBP: missing or 70-210
  # ### Age: 25-84
  # ### Also exclude if BMI<20 as v. different from development cohort
  # 
  # 
  # 
  # qscores <- qscore_vars %>%
  #   
  #   mutate(sex2=ifelse(sex=="male", "male", ifelse(sex=="female", "female", NA))) %>%
  #   
  #   calculate_qrisk2(sex=sex2, age=index_date_age, ethrisk=ethnicity_qrisk2, smoking=qrisk2_smoking_cat, type1=type1, type2=type2, fh_cvd=pre_index_date_fh_premature_cvd, renal=ckd45, af=pre_index_date_af, rheumatoid_arth=pre_index_date_rheumatoidarthritis, cholhdl=precholhdl, sbp=presbp, bmi=prebmi, bp_med=bp_meds, town=tds_2011, surv=surv_5yr) %>%
  #   
  #   rename(qrisk2_5yr_score=qrisk2_score) %>%
  #   
  #   select(-qrisk2_lin_predictor) %>%
  #   
  #   analysis$cached("final_merge_im_q3", unique_indexes="patid")
  # 
  # 
  # 
  # qscores <- qscores %>%
  #   
  #   mutate(sex2=ifelse(sex=="male", "male", ifelse(sex=="female", "female", NA))) %>%
  #   
  #   calculate_qrisk2(sex=sex2, age=index_date_age, ethrisk=ethnicity_qrisk2, smoking=qrisk2_smoking_cat, type1=type1, type2=type2, fh_cvd=pre_index_date_fh_premature_cvd, renal=ckd45, af=pre_index_date_af, rheumatoid_arth=pre_index_date_rheumatoidarthritis, cholhdl=precholhdl, sbp=presbp, bmi=prebmi, bp_med=bp_meds, town=tds_2011, surv=surv_10yr) %>%
  #   
  #   rename(qrisk2_10yr_score=qrisk2_score) %>%
  #   
  #   mutate(across(starts_with("qrisk2"),
  #                 ~ifelse((is.na(precholhdl) | (precholhdl>=1 & precholhdl<=12)) &
  #                           (is.na(presbp) | (presbp>=70 & presbp<=210)) &
  #                           index_date_age>=25 & index_date_age<=84 &
  #                           (is.na(prebmi) | prebmi>=20), .x, NA))) %>%
  #   
  #   select(patid, qdiabeteshf_5yr_score, qdiabeteshf_lin_predictor, qrisk2_5yr_score, qrisk2_10yr_score, qrisk2_lin_predictor) %>%
  #   
  #   analysis$cached("final_merge_im_q4", unique_indexes="patid")
  # 
  # 
  # ## Join with main dataset
  # 
  # final_merge <- final_merge %>%
  #   left_join(qscores, by="patid") %>%
  #   analysis$cached("final_merge_im_2", unique_indexes="patid")
  # 
  # 
  # ############################################################################################
  # 
  # # Add in vars for kidney risk score
  # 
  # ## Make separate table with additional variables
  # 
  # ckdpc_score_vars <- final_merge %>%
  #   
  #   mutate(sex=ifelse(gender==1, "male", ifelse(gender==2, "female", "NA")),
  #          
  #          black_ethnicity=ifelse(!is.na(ethnicity_5cat) & ethnicity_5cat==2, 1L, ifelse(is.na(ethnicity_5cat), NA, 0L)),
  #          
  #          cvd=pre_index_date_myocardialinfarction==1 | pre_index_date_revasc==1 | pre_index_date_heartfailure==1 | pre_index_date_stroke==1,
  #          
  #          oha=ifelse((!is.na(pre_index_date_latest_acarbose) & datediff(index_date, pre_index_date_latest_acarbose)<=183) | 
  #                       (!is.na(pre_index_date_latest_mfn) & datediff(index_date, pre_index_date_latest_mfn)<=183) |
  #                       (!is.na(pre_index_date_latest_dpp4) & datediff(index_date, pre_index_date_latest_dpp4)<=183) |
  #                       (!is.na(pre_index_date_latest_glinide) & datediff(index_date, pre_index_date_latest_glinide)<=183) |
  #                       (!is.na(pre_index_date_latest_glp1) & datediff(index_date, pre_index_date_latest_glp1)<=183) |
  #                       (!is.na(pre_index_date_latest_sglt2) & datediff(index_date, pre_index_date_latest_sglt2)<=183) |
  #                       (!is.na(pre_index_date_latest_su) & datediff(index_date, pre_index_date_latest_su)<=183) |
  #                       (!is.na(pre_index_date_latest_tzd) & datediff(index_date, pre_index_date_latest_tzd)<=183), 1L, 0L),
  #          
  #          # INS=ifelse(!is.na(pre_index_date_latest_insulin) & datediff(index_date, pre_index_date_latest_insulin)<=183, 1L, 0L),
  #          INS = 0L,
  #          
  #          ever_smoker=ifelse(!is.na(smoking_cat) & (smoking_cat=="Ex-smoker" | smoking_cat=="Active smoker"), 1L, ifelse(is.na(smoking_cat), NA, 0L)),
  #          
  #          latest_bp_med=pmax(
  #            ifelse(is.na(pre_index_date_latest_ace_inhibitors), as.Date("1900-01-01"), pre_index_date_latest_ace_inhibitors),
  #            ifelse(is.na(pre_index_date_latest_beta_blockers), as.Date("1900-01-01"), pre_index_date_latest_beta_blockers),
  #            ifelse(is.na(pre_index_date_latest_calcium_channel_blockers), as.Date("1900-01-01"), pre_index_date_latest_calcium_channel_blockers),
  #            ifelse(is.na(pre_index_date_latest_thiazide_diuretics), as.Date("1900-01-01"), pre_index_date_latest_thiazide_diuretics),
  #            na.rm=TRUE
  #          ),
  #          
  #          bp_meds=ifelse(latest_bp_med!=as.Date("1900-01-01") & datediff(index_date, latest_bp_med)<=183, 1L, 0L),
  #          
  #          hypertension=ifelse((!is.na(presbp) & presbp>=140) | (!is.na(predbp) & predbp>=90) | bp_meds==1, 1L,0L),
  #          
  #          uacr=ifelse(!is.na(preacr), preacr, ifelse(!is.na(preacr_from_separate), preacr_from_separate, NA)),
  #          
  #          chd=pre_index_date_myocardialinfarction==1 | pre_index_date_revasc==1,
  #          
  #          current_smoker=ifelse(!is.na(smoking_cat) & smoking_cat=="Active smoker", 1L, ifelse(is.na(smoking_cat), NA, 0L)),
  #          
  #          ex_smoker=ifelse(!is.na(smoking_cat) & smoking_cat=="Ex-smoker", 1L, ifelse(is.na(smoking_cat), NA, 0L))) %>%
  #   
  #   select(patid, index_date_age, sex, black_ethnicity, preegfr, cvd, prehba1c2yrs, INS, oha, ever_smoker, hypertension, prebmi, uacr, presbp, bp_meds, pre_index_date_heartfailure, chd, pre_index_date_af, current_smoker, ex_smoker, preckdstage) %>%
  #   
  #   analysis$cached("final_merge_im_ckd1", unique_indexes="patid")
  # 
  # 
  # 
  # 
  # ## Join with main dataset
  # 
  # final_merge <- final_merge %>%
  #   left_join(ckdpc_score_vars, by="patid") %>%
  #   analysis$cached("final_merge", unique_indexes="patid")
  # 
  # 
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
}
