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


# Define medications

meds <- c("ace_inhibitors",
          "arb",
          "statins",
          "finerenone",
          "sglt2",
          "glp1",
          "beta_blockers",
          "ca_channel_blockers",
          "thiazide_diuretics",
          "loop_diuretics",
          "ksparing_diuretics",
          "definite_genital_infection_meds",
          "topical_candidal_meds")


############################################################################################

# Pull out raw script instances and cache with 'all_patid' prefix
## Some of these already exist from previous analyses

analysis = cprd$analysis("all_patid")

for (i in meds) {
  
  print(i)
  
  raw_tablename <- paste0("raw_", i, "_prodcodes")
  
  
  if (i == "sglt2" | i == "glp1") {
    
    empty_variable = paste0(i, "_cat")
    
    data <- cprd$tables$drugIssue %>% inner_join(
      readr::read_tsv(
        here::here(paste0("C:/Users/tj358/OneDrive - University of Exeter/CPRD/Aurum codelists/prodcodes/exeter_prodcodelist_", i, ".txt")),
        col_types = cols(.default=col_character())) %>%
        rename(prodcodeid=ProdCodeId) %>%
        select(prodcodeid) %>%
        mutate(!!sym(empty_variable) := NA), 
      by="prodcodeid", copy = T) %>%
      select(patid, date=issuedate) %>%
      analysis$cached(raw_tablename, indexes=c("patid", "date"))
    
  } else {
    
    data <- cprd$tables$drugIssue %>%
      inner_join(codes[[i]], by="prodcodeid") %>%
      select(patid, date=issuedate) %>%
      analysis$cached(raw_tablename, indexes=c("patid", "date"))
    
  }
  
  
  assign(raw_tablename, data)
  
}


############################################################################################

analysis = cprd$analysis(analysis_prefix)

# get dates at 6 month intervals
dates <- seq(from = as.Date("2019-03-01"),
             to   = as.Date("2024-03-01"),
             by   = "6 months")

date_strings <- format(dates, "%Y-%m-%d")


for (d in date_strings) {
  print(d)
  index_date <- as.Date(d)
 
  
  for (i in meds) {
    
    print(i)
    
      raw_tablename <- paste0("raw_", i, "_prodcodes")
      index_date_merge_tablename <- paste0(d, "_full_", i, "_merge")
      
      data <- get(raw_tablename) %>%
        
        inner_join(cprd$tables$validDateLookup, by="patid") %>%
        # filter(date>=min_dob & date<=gp_ons_end_date) %>%
        filter(date>=min_dob & date<=gp_end_date) %>%
        select(patid, date) %>%
        
        mutate(datediff=datediff(date, index_date)) %>%
        
        analysis$cached(index_date_merge_tablename, indexes="patid")
      
      assign(index_date_merge_tablename, data)
      
    }
  
  
  ############################################################################################
  
  # Find earliest pre-index date, latest pre-index date and first post-index date dates
  
  
  medications <- cprd$tables$patient %>%
    select(patid)
  
  
  for (i in meds) {
    
    print(paste("working out pre- and post- index date code occurrences for", i, " at ", d))
    
    index_date_merge_tablename <- paste0(d, "_full_", i, "_merge")
    interim_medications_table <- paste0(d, "_meds_im_", i)
    pre_index_date_earliest_date_variable <- paste0("pre_index_date_earliest_", i, "")
    pre_index_date_latest_date_variable <- paste0("pre_index_date_latest_", i, "")
    post_index_date_date_variable <- paste0("post_index_date_first_", i, "")
    
    pre_index_date <- get(index_date_merge_tablename) %>%
      filter(date<=index_date) %>%
      group_by(patid) %>%
      summarise({{pre_index_date_earliest_date_variable}}:=min(date, na.rm=TRUE),
                {{pre_index_date_latest_date_variable}}:=max(date, na.rm=TRUE)) %>%
      ungroup()
    
    post_index_date <- get(index_date_merge_tablename) %>%
      filter(date>index_date) %>%
      group_by(patid) %>%
      summarise({{post_index_date_date_variable}}:=min(date, na.rm=TRUE)) %>%
      ungroup()
    
    medications <- medications %>%
      left_join(pre_index_date, by="patid") %>%
      left_join(post_index_date, by="patid") %>%
      analysis$cached(interim_medications_table, unique_indexes="patid")
    
  }
  
  
  # Cache final version
  
  medications <- medications %>% analysis$cached(paste0(d, "_medications"), unique_indexes="patid")
   
}