############################################################################################

# Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
rm(list=ls())


cprd = CPRDData$new(cprdEnv = "nondiabetes-jun2024",cprdConf = "C:/Users/tj358/OneDrive - University of Exeter/CPRD/aurum.yaml")


codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "01/06/2024")

analysis_prefix <- "ckd"

############################################################################################

biomarkers <- c("creatinine_blood", "acr", "pcr", "albumin_urine", "creatinine_urine",
                "albumin_blood", "haemoglobin", 
                "dbp", "sbp", "weight", "height", "bmi", "totalcholesterol", "hba1c", "hdl", "potassium")
                


############################################################################################

# Pull out all raw biomarker values and cache

analysis = cprd$analysis("all_patid")

for (i in biomarkers) {
  
  print(i)
  
  raw_tablename <- paste0("raw_", i, "_medcodes")
  
  data <- cprd$tables$observation %>%
    inner_join(codes[[i]], by="medcodeid") %>%
    analysis$cached(raw_tablename, indexes=c("patid", "obsdate", "testvalue", "numunitid"))
  
  assign(raw_tablename, data)
  
}


analysis = cprd$analysis("all_patid")


for (i in biomarkers) {
  
  print(i)
  
  raw_tablename <- paste0("raw_", i, "_medcodes")
  clean_tablename <- paste0("clean_", i, "_medcodes")
  
  if (i=="haemoglobin") {
    message("Converting haemoglobin values to g/L")
    raw_data <- get(raw_tablename) %>%
      mutate(testvalue=ifelse(testvalue<30, testvalue*10, testvalue))
  }
  else {
    raw_data <- get(raw_tablename)
  }
  
  
  if (i=="albumin_urine") {
    data <- raw_data %>%
      filter(numunitid==183)
  }
  else if (i=="creatinine_urine") {
    data <- raw_data %>%
      filter(numunitid==218 | numunitid==285) %>%
      mutate(testvalue=ifelse(numunitid==285, testvalue/1000, testvalue))
  }
  else if (i == "hba1c") {
    
    raw_data <- get(raw_tablename) %>%    
      mutate(testvalue=ifelse(testvalue<=20,((testvalue-2.152)/0.09148),testvalue)) 
    
    data <- raw_data %>%
      
      clean_biomarker_values(testvalue, "hba1c") %>%
      clean_biomarker_units(numunitid, "hba1c") 
    
  } else {
    data <- raw_data %>%
      clean_biomarker_units(testvalue, i) %>%
      #clean_biomarker_values(testvalue, i) %>%
      clean_biomarker_units(numunitid, i)
  }
  data <- data %>%
    group_by(patid,obsdate) %>%
    summarise(testvalue=mean(testvalue, na.rm=TRUE)) %>%
    ungroup() %>%
    
    inner_join(cprd$tables$validDateLookup, by="patid") %>%
    #filter(obsdate>=min_dob & obsdate<=gp_ons_end_date) %>%  #gp_ons_end_date not available on this dataset
    filter(obsdate>=min_dob & obsdate<=gp_end_date) %>%
    
    select(patid, date=obsdate, testvalue) %>%
    
    analysis$cached(clean_tablename, indexes=c("patid", "date", "testvalue"))
  
  assign(clean_tablename, data)
  
}



# egfr
analysis = cprd$analysis("all")
dob <- dob %>% analysis$cached("dob")

analysis = cprd$analysis("all_patid")
clean_egfr_medcodes <- clean_creatinine_blood_medcodes %>%
  
  inner_join((dob %>% select(patid, dob)), by="patid") %>%
  inner_join((cprd$tables$patient %>% select(patid, gender)), by="patid") %>%
  mutate(age_at_creat=(datediff(date, dob))/365.25,
         sex=ifelse(gender==1, "male", ifelse(gender==2, "female", NA))) %>%
  select(-c(dob, gender)) %>%
  
  ckd_epi_2021_egfr(creatinine=testvalue, sex=sex, age_at_creatinine=age_at_creat) %>%
  select(-c(testvalue, sex, age_at_creat)) %>%
  
  rename(testvalue=ckd_epi_2021_egfr) %>%
  filter(!is.na(testvalue)) %>%
  analysis$cached("clean_egfr_medcodes", indexes=c("patid", "date", "testvalue"))

biomarkers <- c("egfr", biomarkers)

# Make ACR from separate urine albumin and urine creatinine measurements on the same day
# Then clean values

clean_acr_from_separate_medcodes <- clean_albumin_urine_medcodes %>%
  inner_join((clean_creatinine_urine_medcodes %>% select(patid, creat_date=date, creat_value=testvalue)), by="patid") %>%
  filter(date==creat_date) %>%
  mutate(new_testvalue=testvalue/creat_value) %>%
  select(patid, date, testvalue=new_testvalue) %>%
 #  clean_biomarker_units_acr(testvalue, i) %>%
  clean_biomarker_values(testvalue, "acr") %>%
  analysis$cached("clean_acr_from_separate_medcodes", indexes=c("patid", "date", "testvalue"))

biomarkers <- setdiff(biomarkers, c("albumin_urine", "creatinine_urine"))
biomarkers <- c("acr_from_separate", biomarkers)

######################################################################################
analysis = cprd$analysis(analysis_prefix)

# 6-monthly dates for 2019-2021 (prevalent cohort), then 3-monthly from 2021 onwards
# (3-monthly required for sequential trial emulation of SGLT2i in non-DM CKD)
dates <- unique(c(
  seq(from = as.Date("2019-03-01"), to = as.Date("2020-09-01"), by = "6 months"),
  seq(from = as.Date("2021-03-01"), to = as.Date("2024-03-01"), by = "3 months")
))

date_strings <- format(dates, "%Y-%m-%d")


for (d in date_strings) {
  print(d)
  index_date = as.Date(d)
  
  for (i in biomarkers) {
    
    clean_tablename <- paste0("clean_", i, "_medcodes")
    index_date_merge_tablename <- paste0(d, "_full_", i, "_merge")
    
    data <- get(clean_tablename) %>%
      mutate(datediff=datediff(date, index_date))
    
    assign(index_date_merge_tablename, data)
    
  }
  
  
  
  ############################################################################################
  
  # Find baseline values
  ## Within period defined above (-2 years to +7 days for all except height)
  ## Then use closest date to index date
  ## May be multiple values; use minimum test result, except for eGFR - use maximum
  ## Can get duplicates where person has identical results on the same day/days equidistant from the index date - choose first row when ordered by datediff
  
  baseline_biomarkers <- cprd$tables$patient %>%
    select(patid)
  
  
  ## For all except height: between 2 years prior and 7 days after index date
  
  biomarkers_no_height <- setdiff(biomarkers, "height")
  
  for (i in biomarkers_no_height) {
    
    print(i)
    
    index_date_merge_tablename <- paste0(d, "_full_", i, "_merge")
    interim_baseline_biomarker_table <- paste0(d, "_biomarkers_im_", i)
    pre_biomarker_variable <- paste0("pre", i)
    pre_biomarker_date_variable <- paste0("pre", i, "date")
    pre_biomarker_datediff_variable <- paste0("pre", i, "datediff")
    
    
    data <- get(index_date_merge_tablename) %>%
      filter(datediff<=7 & datediff>=-730) %>%
      
      group_by(patid) %>%
      
      mutate(min_timediff=min(abs(datediff), na.rm=TRUE)) %>%
      filter(abs(datediff)==min_timediff) %>%
      
      mutate(pre_biomarker=ifelse(i=="egfr", max(testvalue, na.rm=TRUE), min(testvalue, na.rm=TRUE))) %>%
      filter(pre_biomarker==testvalue) %>%
      
      dbplyr::window_order(datediff) %>%
      filter(row_number()==1) %>%
      
      ungroup() %>%
      
      relocate(pre_biomarker, .after=patid) %>%
      relocate(date, .after=pre_biomarker) %>%
      relocate(datediff, .after=date) %>%
      
      rename({{pre_biomarker_variable}}:=pre_biomarker,
             {{pre_biomarker_date_variable}}:=date,
             {{pre_biomarker_datediff_variable}}:=datediff) %>%
      
      select(-c(testvalue, min_timediff))
    
    
    baseline_biomarkers <- baseline_biomarkers %>%
      left_join(data, by="patid") %>%
      analysis$cached(interim_baseline_biomarker_table, unique_indexes="patid")
    
  }
  
  
  ## Height - only keep readings at/post-index date, and find mean
  
  table_name = paste0(d, "_full_height_merge")
  
  baseline_height <- get(table_name) %>%
    filter(datediff>=0) %>%
    group_by(patid) %>%
    summarise(height=mean(testvalue, na.rm=TRUE)) %>%
    ungroup()
  
  baseline_biomarkers <- baseline_biomarkers %>%
    left_join(baseline_height, by="patid") %>%
    analysis$cached(paste0(d, "_baseline_biomarkers"), unique_indexes="patid")
  
}
