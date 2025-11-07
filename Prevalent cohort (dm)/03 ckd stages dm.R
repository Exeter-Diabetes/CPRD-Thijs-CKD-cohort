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

# get ckd stages
analysis = cprd$analysis("all_patid")

ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  analysis$cached("ckd_stages_from_algorithm_with_acr",
                  indexes = c("patid"))


######################################################################################
analysis = cprd$analysis(analysis_prefix)

# get dates at 6 month intervals
dates <- seq(from = as.Date("2019-03-01"),
             to   = as.Date("2024-03-01"),
             by   = "6 months")

date_strings <- format(dates, "%Y-%m-%d")

# define CKD stage at each time point
for (d in date_strings) {
  print(d)
  index_date <- as.Date(d)
  
  
  ckd_stage_drug_merge <- cprd$tables$patient %>%
    select(patid) %>%
    left_join(ckd_stages_from_algorithm, by="patid") %>%
    mutate(preckdstage=ifelse(!is.na(stage_5) & datediff(stage_5, index_date)<=7, "stage_5",
                              ifelse(!is.na(stage_4) & datediff(stage_4, index_date)<=7, "stage_4",
                                     ifelse(!is.na(stage_3b) & datediff(stage_3b, index_date)<=7, "stage_3b",
                                            ifelse(!is.na(stage_3a) & datediff(stage_3a, index_date)<=7, "stage_3a",
                                                   ifelse(!is.na(stage_2) & datediff(stage_2, index_date)<=7, "stage_2",
                                                          ifelse(!is.na(stage_1) & datediff(stage_1, index_date)<=7, "stage_1", NA)))))),
           
           preckdstagedate=ifelse(preckdstage=="stage_5", stage_5,
                                  ifelse(preckdstage=="stage_4", stage_4,
                                         ifelse(preckdstage=="stage_3b", stage_3b,
                                                ifelse(preckdstage=="stage_3a", stage_3a,
                                                       ifelse(preckdstage=="stage_2", stage_2,
                                                              ifelse(preckdstage=="stage_1", stage_1, NA)))))),
           
           preckdstagedatediff=datediff(preckdstagedate, index_date)) %>%
    
    select(patid, preckdstage, preckdstagedate, preckdstagedatediff) %>%
    
    analysis$cached(paste0(d, "_ckd_stages"), unique_indexes="patid")
  
}
