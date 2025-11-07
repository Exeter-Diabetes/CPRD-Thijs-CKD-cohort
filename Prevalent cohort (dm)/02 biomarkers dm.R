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

biomarkers <- c("egfr", "acr", "pcr", "acr_from_separate",
                "albumin_blood", "haemoglobin", 
                "dbp", "sbp", "weight", "height", "bmi", "totalcholesterol", "hba1c", "hdl")



############################################################################################

# Pull out all clean biomarker values

analysis = cprd$analysis("all_patid")

for (i in biomarkers) {
  
  print(i)
  
  clean_tablename <- paste0("clean_", i, "_medcodes")
  
  data <- data %>% analysis$cached(clean_tablename)
  
  assign(clean_tablename, data)
  
}

######################################################################################
analysis = cprd$analysis(analysis_prefix)

# get dates at 6 month intervals
dates <- seq(from = as.Date("2019-03-01"),
             to   = as.Date("2024-03-01"),
             by   = "6 months")

date_strings <- format(dates, "%Y-%m-%d")


for (d in date_strings) {
  index_date = as.Date(d)
  print(d)
  
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
