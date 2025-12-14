############################################################################################

# Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
rm(list=ls())


cprd = CPRDData$new(cprdEnv = "nondiabetes-jun2024",cprdConf = "C:/Users/tj358/OneDrive - University of Exeter/CPRD/aurum.yaml")

today <- format(Sys.Date(), "%Y%m%d")

codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "01/06/2024")

analysis_prefix <- "ckd"

############################################################################################

## Cohort and patient characteristics
analysis = cprd$analysis(analysis_prefix)
ckd_cohort <- ckd_cohort %>% analysis$cached("ckd_cohort")
analysis = cprd$analysis("all_patid")
all_ids <- all_ids %>% analysis$cached("all_ids")

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
    filter(first_ckd_date<=index_date & regstartdate<=index_date & gp_end_date>=index_date & (is.na(death_date) | death_date>=index_date)) %>%
    select(patid) %>%
    analysis$cached("cohort_ids", unique_indexes="patid")
  
  # get counts of all patients + patients with CKD at index date
  total_count <- all_ids %>% 
    filter(regstartdate<=index_date & gp_end_date>=index_date & (is.na(death_date) | death_date>=index_date)) %>%
    select(patid) %>%
    count() %>% 
    collect() %>% 
    pull(n)
  total_valid_egfr_count <- all_ids %>% 
    filter(regstartdate<=index_date & gp_end_date>=index_date & (is.na(death_date) | death_date>=index_date)) %>%
    left_join(baseline_biomarkers, by = "patid") %>%
    filter(!is.na(preegfr)) %>%
    select(patid) %>%
    count() %>% 
    collect() %>% 
    pull(n)
  total_valid_egfr_uacr_count <- all_ids %>% 
    filter(regstartdate<=index_date & gp_end_date>=index_date & (is.na(death_date) | death_date>=index_date)) %>%
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
    mutate(index_date_age=datediff(index_date, dob)/365.25,
           index_date_ckd_dur_all=datediff(index_date, first_ckd_date)/365.25,
           index_date = index_date) %>%
    relocate(c(index_date_age, index_date_ckd_dur_all), .before=gender) %>%
    analysis$cached("final_merge", unique_indexes="patid")
  
  
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
  df_name <- paste0("prev_", gsub("-", "_", d), "_nondm")
  
  # Assign name
  assign(df_name, prev_cohort, envir = .GlobalEnv)
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Raw data/")
  save(list = df_name, file=paste0(today, "_prev_ckd_cohort_nondm_", d, ".Rda"))
  
  rm(medications)
  rm(baseline_biomarkers)
  rm(comorbidities)
  rm(ckd_stages)
  rm(smoking)
  rm(cohort_ids)
  rm(final_merge)
}

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Raw data/")
save(counts, file=paste0(today, "_ckd_counts_nondm.Rda"))
