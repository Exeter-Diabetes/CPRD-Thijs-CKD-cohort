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


# Define medications

meds <- c("ace_inhibitors",
          "arb",
          "statins",
          "finerenone",
          "beta_blockers",
          "ca_channel_blockers",
          "thiazide_diuretics",
          "loop_diuretics",
          "ksparing_diuretics",
          "definite_genital_infection_meds",
          "topical_candidal_meds")


############################################################################################

# Pull out clean script instances

analysis = cprd$analysis("all_patid")

for (i in meds) {
  
  print(i)
  
  raw_tablename <- paste0("raw_", i, "_prodcodes")
  
  data <- data %>% analysis$cached(raw_tablename)
  
  assign(raw_tablename, data)
  
}

raw_oha_prodcodes <- cprd$tables$drugIssue %>%
  inner_join(cprd$tables$ohaLookup, by="prodcodeid") %>%
  analysis$cached("raw_oha_prodcodes", indexes=c("patid", "issuedate"))

raw_insulin_prodcodes <- cprd$tables$drugIssue %>%
  inner_join(codes$insulin, by="prodcodeid") %>%
  analysis$cached("raw_insulin_prodcodes", indexes=c("patid", "issuedate"))

clean_oha_prodcodes <- raw_oha_prodcodes %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  filter(issuedate>=min_dob & issuedate<=gp_end_date) %>%
  rename(date=issuedate) %>%
  analysis$cached("clean_oha_prodcodes", indexes=c("patid", "issuedate"))

clean_glp1_prodcodes <- clean_oha_prodcodes %>%
  filter(drug_class_1=="GLP1" | drug_class_2=="GLP1") %>%
  select(patid, date)

clean_sglt2_prodcodes <- clean_oha_prodcodes %>%
  filter(drug_class_1=="SGLT2" | drug_class_2=="SGLT2") %>%
  select(patid, date)

clean_insulin_prodcodes <- raw_insulin_prodcodes %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  filter(issuedate>=min_dob & issuedate<=gp_end_date) %>%
  rename(date=issuedate) %>%
  analysis$cached("clean_insulin_prodcodes", indexes=c("patid", "issuedate"))

clean_insulin_prodcodes <- clean_insulin_prodcodes %>%
  select(patid, date) %>%
  union_all(clean_oha_prodcodes %>%
              filter(drug_class_1=="INS" | drug_class_2=="INS") %>%
              select(patid, date))


meds <- c(meds, "sglt2", "glp1", "insulin")

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
    
    if (i=="glp1" | i=="sglt2" | i=="insulin") {
      
      clean_tablename <- paste0("clean_", i, "_prodcodes")
      index_date_merge_tablename <- paste0(d, "_full_", i, "_merge")
      
      data <- get(clean_tablename) %>%
        mutate(datediff=datediff(date, index_date)) %>%
        analysis$cached(index_date_merge_tablename, indexes="patid")
      
      assign(index_date_merge_tablename, data)
      
    } else {
    
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
