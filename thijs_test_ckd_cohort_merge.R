############################################################################################

# Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "test-remote-full", cprdConf = "C:/Users/tj358/OneDrive - University of Exeter/CPRD/aurum.yaml")
codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "31/10/2021")
current_list <- codesets$listCodeSets()


minimum_date <- as.Date("2010-01-01")
end_date <- as.Date("2020-10-31")

analysis_prefix <- "thijs_test"
analysis = cprd$analysis(analysis_prefix)

############################################################################################


advanced_ckd <- advanced_ckd %>%
  analysis$cached("advanced_ckd_index_date")

## combine with ethnicity and basic data

analysis = cprd$analysis("all")
ethnicity <- ethnicity %>% analysis$cached("gp_ethnicity")
language <- language %>% analysis$cached("gp_language")

analysis = cprd$analysis(analysis_prefix)
comorbidities <- comorbidities %>% analysis$cached("comorbidities_advanced_ckd", unique_indexes="patid")
medications <- medications %>% analysis$cached("medications_advanced_ckd", unique_indexes = "patid")
baseline_biomarkers <- baseline_biomarkers %>% analysis$cached("biomarkers_advanced_ckd")

#### join up and cache
analysis = cprd$analysis(analysis_prefix)
advanced_ckd <- advanced_ckd %>%
  inner_join((dob %>% select(patid, dob)), by="patid") %>%
  left_join((cprd$tables$patient %>% select(patid, gender, regstartdate, regenddate, pracid, cprd_ddate)), by="patid") %>%
  left_join((cprd$tables$practice %>% select(pracid, lcd, region)), by="pracid") %>%
  ## I have disabled the following lines in the script referring to IMD, HES, and ONS (linkage not available)
  #  left_join((cprd$tables$patientImd2015 %>% select(patid, imd2015_10)), by="patid") %>% # IMD 2015 not present
  #  left_join((cprd$tables$validDateLookup %>% select(patid, ons_death)), by="patid") %>% # ONS linkage not present
  #  left_join((cprd$tables$patidsWithLinkage %>% select(patid, n_patid_hes)), by="patid") %>% # HES linkage not present
  
  mutate(gp_record_end=pmin(if_else(is.na(lcd), end_date, lcd),
                            if_else(is.na(regenddate), end_date, regenddate),
                            end_date, na.rm=TRUE),
         
         #         death_date=pmin(if_else(is.na(cprd_ddate), as.Date("2050-01-01"), cprd_ddate),
         #                         if_else(is.na(ons_death), as.Date("2050-01-01"), ons_death), na.rm=TRUE),
         #         death_date=if_else(death_date==as.Date("2050-01-01"), as.Date(NA), death_date),
         #         
         #         with_hes=ifelse(!is.na(n_patid_hes) & n_patid_hes<=20, 1L, 0L)
  ) %>%
  
  left_join(ethnicity, by = "patid") %>%
  left_join(language,  by = "patid") %>%
  
  select(patid, gender, dob, index_date, pracid, prac_region=region, 
         preckdstage, preckdstagedate, 
         gp_5cat_ethnicity, gp_16cat_ethnicity, gp_qrisk2_ethnicity, 
         language_cat=gp_language_cat,
         #         imd2015_10, 
         regstartdate, gp_record_end, 
         #         death_date, with_hes
  ) %>% 
  
  analysis$cached("ckd_cohort_im", unique_indexes="patid", indexes=c("gender", "dob", "index_date"))

ckd_cohort <- advanced_ckd %>%
  inner_join(comorbidities, by="patid") %>%
  inner_join(medications, by="patid") %>%
  inner_join(baseline_biomarkers, by="patid")
analysis$cached("ckd_cohort", unique_indexes="patid", indexes=c("gender", "dob", "index_date"))