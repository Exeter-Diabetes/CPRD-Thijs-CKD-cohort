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
start_date = "2024-03-01"
end_date = "2024-03-01"

############################################################################################

## Cohort and patient characteristics
analysis = cprd$analysis(analysis_prefix)
hf_cohort <- hf_cohort %>% analysis$cached("hf_cohort")
analysis = cprd$analysis("all_patid")
all_ids <- all_ids %>% analysis$cached("all_ids")
townsend_score <- townsend_score %>% analysis$cached("townsend_score")
analysis = cprd$analysis("all")
death_causes <- death_causes %>% analysis$cached("death_causes")


## Get index date

# get dates at 6 month intervals
dates <- seq(from = as.Date(start_date),
             to   = as.Date(end_date),
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
  
  hf_ids <- hf_cohort %>%
    filter(first_hf_date<=index_date & regstartdate<=index_date & gp_end_date>=index_date & (is.na(death_date) | death_date>=index_date)) %>%
    select(patid) %>%
    analysis$cached("hf_valid_ids", unique_indexes="patid")
   
  final_merge <- hf_ids %>%
    left_join(hf_cohort, by="patid") %>%
    left_join(ckd_stages, by="patid") %>%
    left_join(baseline_biomarkers, by="patid") %>%
    left_join(comorbidities, by="patid") %>%
    left_join(smoking, by="patid") %>%
    left_join(medications, by="patid") %>%
    left_join(townsend_score %>% select(patid, tds_2011), by = "patid") %>% 
    left_join(death_causes, by = "patid") %>%
    mutate(index_date_age=datediff(index_date, dob)/365.25,
           index_date_hf_dur_all=datediff(index_date, first_hf_date)/365.25,
           index_date = index_date) %>%
    relocate(c(index_date_age, index_date_hf_dur_all), .before=gender) %>%
    analysis$cached("final_merge_hf", unique_indexes="patid")
  
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
  df_name <- paste0("prev_", gsub("-", "_", d), "_hf_nondm")
  
  # Assign name
  assign(df_name, prev_cohort, envir = .GlobalEnv)
  
  setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Raw data/")
  save(list = df_name, file=paste0(today, "_prev_hf_cohort_nondm_", d, ".Rda"))
  
  rm(medications)
  rm(baseline_biomarkers)
  rm(comorbidities)
  rm(ckd_stages)
  rm(smoking)
  rm(cohort_ids)
  rm(final_merge)
}