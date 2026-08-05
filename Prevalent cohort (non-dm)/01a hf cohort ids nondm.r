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

#Data quality check - should only include acceptable' patients (see CPRD data specification for definition)
cprd$tables$patient %>% count() #45,037,869 - total patient count in download
cprd$tables$patient %>% filter(acceptable ==1) %>% count() #45,037,869
cprd$tables$patient %>% filter(patienttypeid ==3) %>% count() #45,037,869
#All are 'acceptable' and have patienttypeid==3 ('Regular')

############################################################################################

##CPRD recommend excluding 44 practices (as below) that appear likely to have merged into other contributing practices (patient data could be duplicated)
##Define patients to remove later
analysis = cprd$analysis("all_patid")

practice_exclusion_ids <- cprd$tables$patient %>% 
  filter(pracid == "20024" | pracid == "20036" |pracid == "20091" |pracid == "20171" | pracid == "20178" |pracid == "20202" | pracid == "20254" | pracid == "20389" |pracid == "20430" |pracid == "20452" |
           pracid == "20469" | pracid == "20487" | pracid == "20552" | pracid == "20554" | pracid == "20640" | pracid == "20717" | pracid == "20734" | pracid == "20737" | pracid == "20740" | pracid == "20790" |
           pracid == "20803" | pracid == "20822" | pracid == "20868" | pracid == "20912" | pracid == "20996" | pracid == "21001" | pracid == "21015" | pracid == "21078" | pracid == "21112" | pracid == "21118" |
           pracid == "21172" | pracid == "21173" | pracid == "21277" | pracid == "21281" | pracid == "21331" | pracid == "21334" | pracid == "21390" | pracid == "21430" | pracid == "21444" | pracid == "21451" |
           pracid == "21529" | pracid == "21553" | pracid == "21558" | pracid == "21585") %>%
  analysis$cached("practice_exclusion_ids")

practice_exclusion_ids %>% count() #672,504


############################################################################################

##Define patients with gender=3 (indeterminate) to remove later

gender_exclusion_ids <- cprd$tables$patient %>% 
  filter(gender==3) %>%
  analysis$cached("gender_exclusion_ids")

gender_exclusion_ids %>% count() #1767

cprd$tables$patient %>% anti_join(practice_exclusion_ids, by="patid") %>% anti_join(gender_exclusion_ids, by="patid") %>% count() #44,363,638


############################################################################################

analysis = cprd$analysis(analysis_prefix)

# load table with heart failure diagnosis (from comorbidities table)
heartfailure_table <- heartfailure_table %>%
  analysis$cached("2024-03-01_full_heartfailure_merge", indexes = c("patid"))



hf_ids <- heartfailure_table %>%
  mutate(first_hf_date = date) %>%
  group_by(patid) %>%
  dbplyr::window_order(first_hf_date) %>%
  distinct(patid, .keep_all = TRUE) %>%
  ungroup() %>%
  distinct(patid, .keep_all = TRUE) %>%
  select(patid, first_hf_date, source) %>%
  anti_join(practice_exclusion_ids, by="patid") %>% anti_join(gender_exclusion_ids, by="patid") %>%
  analysis$cached("hf_ids", unique_indexes="patid")

hf_ids %>% count() #1349229


############################################################################################

# join with tables

# dob 
analysis = cprd$analysis("all")

dob <- cprd$tables$observation %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  filter(obsdate>=min_dob) %>%
  group_by(patid) %>%
  summarise(earliest_medcode=min(obsdate, na.rm=TRUE)) %>%
  ungroup() %>%
  analysis$cached("earliest_medcode", unique_indexes="patid")

#### Check count
dob %>% count() #44,960,468 - almost everyone in download

#### No-one has missing dob or earliest_medcode so pmin (runs as 'LEAST' in MySQL) works
dob <- dob %>%
  inner_join(cprd$tables$patient, by="patid") %>%
  mutate(dob=as.Date(ifelse(is.na(mob), paste0(yob,"-06-30"), paste0(yob, "-",mob,"-15")))) %>%
  inner_join(cprd$tables$validDateLookup, by = "patid") %>%
  mutate(dob=pmin(dob, earliest_medcode, na.rm=TRUE)) %>%
  mutate(dob=ifelse(regstartdate>=min_dob & regstartdate<dob, regstartdate, dob)) %>%
  select(patid, dob = dob, mob, yob, regstartdate) %>%
  analysis$cached("dob", unique_indexes="patid")

dob <- dob %>% analysis$cached("dob", unique_indexes="patid")

analysis = cprd$analysis("all_patid")
ethnicity <- ethnicity %>% analysis$cached("ethnicity", unique_indexes="patid")

# get list of all ids
all_ids <- dob %>%
  anti_join(practice_exclusion_ids, by="patid") %>% 
  anti_join(gender_exclusion_ids, by="patid") %>%
  left_join((cprd$tables$patient %>% select(patid, gender, regenddate, pracid)), by="patid") %>%
  left_join((cprd$tables$practice %>% select(pracid, lcd, region)), by="pracid") %>%
  left_join((cprd$tables$onsDeath %>% select(patid, reg_date_of_death)), by="patid") %>%
  left_join((cprd$tables$patientImd %>% select(patid, imd_decile)), by="patid") %>%
  left_join((cprd$tables$validDateLookup %>% select(patid, gp_end_date)), by="patid") %>%
  left_join((cprd$tables$patidsWithLinkage %>% mutate(with_hes=1L) %>% select(patid, with_hes, hes_end_date)), by="patid") %>%
  mutate(with_hes=ifelse(is.na(with_hes), 0L, 1L)) %>%
  left_join(ethnicity, by="patid") %>%
  select(patid, gender, dob, pracid, prac_region=region, ethnicity_5cat, ethnicity_16cat, ethnicity_qrisk2, imd_decile, regstartdate, gp_end_date, death_date=reg_date_of_death, with_hes, hes_end_date) %>%
  analysis$cached("all_ids", unique_indexes="patid", indexes=c("gender", "dob"))

all_ids %>% count() 
#44,363,638

# join ids with dob and other data
analysis = cprd$analysis(analysis_prefix)

hf_cohort <- hf_ids %>%
  left_join(dob, by="patid") %>%
  left_join((cprd$tables$patient %>% select(patid, gender, regenddate, pracid)), by="patid") %>%
  left_join((cprd$tables$practice %>% select(pracid, lcd, region)), by="pracid") %>%
  left_join((cprd$tables$onsDeath %>% select(patid, reg_date_of_death)), by="patid") %>%
  left_join((cprd$tables$patientImd %>% select(patid, imd_decile)), by="patid") %>%
  left_join((cprd$tables$validDateLookup %>% select(patid, gp_end_date)), by="patid") %>%
  left_join((cprd$tables$patidsWithLinkage %>% mutate(with_hes=1L) %>% select(patid, with_hes, hes_end_date)), by="patid") %>%
  mutate(with_hes=ifelse(is.na(with_hes), 0L, 1L)) %>%
  left_join(ethnicity, by="patid") %>%
  select(patid, gender, dob, pracid, prac_region=region, ethnicity_5cat, ethnicity_16cat, ethnicity_qrisk2, 
  imd_decile, regstartdate, gp_end_date, death_date=reg_date_of_death, with_hes, hes_end_date, first_hf_date) %>%
  analysis$cached("hf_cohort", unique_indexes="patid", indexes=c("gender", "dob"))
                  
                  
hf_cohort %>% count() # 1,349,229

############################################################################################
