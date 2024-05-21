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
index_date <- as.Date("2020-07-01")
end_date <- as.Date("2020-10-31")

analysis_prefix <- "thijs_test"
analysis = cprd$analysis(analysis_prefix)

############################################################################################

# Get clean creatinine readings and convert to eGFR

## Get raw creatinine readings
raw_creatinine_blood_medcodes <- cprd$tables$observation %>%
  inner_join(codes$creatinine_blood, by="medcodeid") %>%
  filter(obsdate > minimum_date) %>%  # select observations after "minimum date" only
  analysis$cached("raw_creatinine_blood_medcodes", indexes=c("patid", "obsdate", "testvalue", "numunitid"))

## Clean creatinine readings
clean_creatinine_blood_medcodes <- raw_creatinine_blood_medcodes %>%
  group_by(patid, obsdate) %>%
  summarise(testvalue=mean(testvalue, na.rm=TRUE)) %>%
  ungroup() %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  # filter(obsdate>=min_dob & obsdate<=gp_ons_end_date) %>% ## "full" dataset does not have linkage to ONS data, therefore variable gp_ons_end_date substituted with gp_end_date
  filter(obsdate>=min_dob & obsdate<=gp_end_date) %>%
  select(patid, date=obsdate, testvalue) %>%
  analysis$cached("clean_creatinine_blood_medcodes", indexes=c("patid", "date", "testvalue"))

clean_creatinine_blood_medcodes %>% count()
#74,511,626

## Convert to eGFR

### Need DOB from DOB table
dob <- cprd$tables$observation %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  filter(obsdate>=min_dob) %>%
  group_by(patid) %>%
  summarise(earliest_medcode=min(obsdate, na.rm=TRUE)) %>%
  analysis$cached("earliest_medcode", unique_indexes="patid")

#### Check count
dob %>% count()
#### 15,705,720 

#### No-one has missing dob or earliest_medcode so pmin (runs as 'LEAST' in MySQL) works
dob <- dob %>%
  inner_join(cprd$tables$patient, by="patid") %>%
  mutate(dob=as.Date(ifelse(is.na(mob), paste0(yob,"-07-01"), paste0(yob, "-",mob,"-15")))) %>%
  mutate(dob=pmin(dob, earliest_medcode, na.rm=TRUE)) %>%
  select(patid, dob, mob, yob, regstartdate) %>%
  analysis$cached("dob", unique_indexes="patid")

# there are a number of "inplausible" dates, however they do not match any of the rows in clean_egfr_medcodes so they do not pose a problem

# join with dob table to get age/sex; remove negative / zero values of creatinine
clean_egfr_medcodes <- clean_creatinine_blood_medcodes %>%
  inner_join((dob %>% select(patid, dob)), by="patid") %>%
  inner_join((cprd$tables$patient %>% select(patid, gender)), by="patid") %>%
  mutate(age_at_creat=(datediff(date, dob))/365.25,
         sex=ifelse(gender==1, "male", ifelse(gender==2, "female", NA))) %>%
  select(-c(dob, gender)) %>%
  rename(creatinine=testvalue) %>%
  mutate(creatinine=ifelse(creatinine < 0, -creatinine, ifelse(creatinine == 0, NA, creatinine))) %>%
  analysis$cached("clean_egfr_medcodes_interim", indexes=c("patid", "date"))

clean_egfr_medcodes %>% count()
# 6,081,640 - # 614,840 with missing dob (patid does not appear in dob table)

# calculate egfr and save
clean_egfr_medcodes <- clean_egfr_medcodes %>%
  ckd_epi_2021_egfr(creatinine=creatinine, sex=sex, age_at_creatinine=age_at_creat) %>%
  select(-c(sex, age_at_creat)) %>%
  rename(egfr=ckd_epi_2021_egfr) %>%
  filter(!is.na(egfr)) %>%
  analysis$cached("clean_egfr_medcodes", indexes=c("patid", "date", "egfr"))

clean_egfr_medcodes %>% count()
#6,040,559 - lose a 196 readings of people with sex == NA and 40,886 with missing creatinine (1 with both missing) therefore in total 41,081

full_egfr_index_date_merge <- clean_egfr_medcodes %>%
  mutate(datediff=datediff(date, index_date)) %>%
  filter(datediff<=7 & datediff>=-730) %>%
  
  group_by(patid) %>%
  
  mutate(min_timediff=min(abs(datediff), na.rm=TRUE)) %>%
  filter(abs(datediff)==min_timediff) %>%
  
  mutate(pre_egfr=max(egfr, na.rm=TRUE)) %>%
  filter(pre_egfr==egfr) %>%
  
  dbplyr::window_order(datediff) %>%
  filter(row_number()==1) %>%
  
  ungroup() %>%
  
  relocate(pre_egfr, .after=patid) %>%
  relocate(date, .after=pre_egfr) %>%
  relocate(datediff, .after=date) %>%
  
  rename(preegfr=pre_egfr,
         preegfrdate=date,
         preegfrdatediff=datediff) %>%
  
  select(-c(egfr, min_timediff)) %>%
  analysis$cached("full_egfr_index_date_merge", unique_indexes = "patid")


################################################################################################################################

# Convert eGFR to CKD stage

ckd_stages_from_all_egfr <- clean_egfr_medcodes %>%
  mutate(ckd_stage=ifelse(egfr<15, "stage_5",
                          ifelse(egfr<30, "stage_4",
                                 ifelse(egfr<45, "stage_3b",
                                        ifelse(egfr<60, "stage_3a",
                                               ifelse(egfr<90, "stage_2",
                                                      ifelse(egfr>=90, "stage_1", NA)))))))


################################################################################################################################

# Only keep CKD stages if >1 consecutive test with the same stage, and if time between earliest and latest consecutive test with same stage are >=90 days apart

## For each patient:
### A) Define period from current test until next test as having the ckd_stage of current test
### B) Join together consecutive periods with the same ckd_stage
### C) If period contains >1 test, and there is >=90 days between the first and last test in the period, it is 'confirmed'


### A) Define period from current test until next test as having the ckd_stage of current test

#### Add in row labelling within each patient's values + max number of rows for each patient

ckd_stages_from_algorithm <- ckd_stages_from_all_egfr %>%
  group_by(patid) %>%
  dbplyr::window_order(date) %>%
  mutate(patid_row_id=row_number()) %>%
  mutate(patid_total_rows=max(patid_row_id, na.rm=TRUE)) %>%
  ungroup()


#### For rows where there is a next test, use this as end date; for last row, use start date as end date

ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  mutate(next_row=patid_row_id+1) %>%
  left_join(ckd_stages_from_algorithm, by=c("patid","next_row"="patid_row_id")) %>%
  mutate(ckd_start=date.x,
         ckd_end=if_else(is.na(date.y),date.x,date.y),
         ckd_stage=ckd_stage.x) %>%
  select(patid, patid_row_id, ckd_stage, ckd_start, ckd_end)


### B) Join together consecutive periods with the same ckd_stage

ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  group_by(patid, ckd_stage) %>%
  dbplyr::window_order(patid, ckd_stage, patid_row_id) %>%
  mutate(lead_var=lead(ckd_start),
         cummax_var=cummax(ckd_end)) %>%
  mutate(compare=cumsum(lead_var>cummax_var)) %>%
  mutate(indx=ifelse(row_number()==1, 0L, lag(compare))) %>%
  ungroup() %>%
  group_by(patid, ckd_stage ,indx) %>%
  summarise(first_test_date=min(ckd_start,na.rm=TRUE),
            last_test_date=max(ckd_start,na.rm=TRUE),
            end_date=max(ckd_end,na.rm=TRUE),
            test_count=max(patid_row_id, na.rm=TRUE)-min(patid_row_id, na.rm=TRUE)+1) %>%
  ungroup() %>%
  analysis$cached("ckd_stages_from_algorithm_interim_1",indexes=c("patid", "ckd_stage", "test_count", "first_test_date", "last_test_date"))

ckd_stages_from_algorithm %>% count()
#2,959,465

ckd_stages_from_algorithm %>% summarise(total=sum(test_count, na.rm=TRUE))
#total number of tests: 6,040,559 as above


### C) Remove periods with 1 reading, or with multiple readings but <90 days between first and last test, and cache

ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  filter(test_count>1 & datediff(last_test_date, first_test_date)>=90) %>%
  analysis$cached("ckd_stages_from_algorithm_interim_2",indexes=c("patid","ckd_stage","first_test_date"))

ckd_stages_from_algorithm %>% count()
#1,136,411


################################################################################################################################

# Combine with CKD5 medcodes/ICD10/OPCS4 codes

## Get raw CKD5 codes and clean
### All are already in all_patid tables on MySQL from 4_mm_comorbidities script
analysis = cprd$analysis("all_patid")
### Medcodes
raw_ckd5_code_medcodes <- raw_ckd5_medcodes %>% analysis$cached("raw_ckd5_code_medcodes")

### ICD10 codes
raw_ckd5_code_icd10 <- raw_ckd5_icd10 %>% analysis$cached("raw_ckd5_code_icd10")

### OPCS4 codes
raw_ckd5_code_opcs4 <- raw_ckd5_opcs4 %>% analysis$cached("raw_ckd5_code_opcs4")

analysis = cprd$analysis(analysis_prefix)
## Clean, find earliest date per person, and re-cache

earliest_clean_ckd5 <- raw_ckd5_code_medcodes %>%
  select(patid, date=obsdate) %>%
  mutate(source="gp") %>%
  union_all((raw_ckd5_code_icd10 %>% select(patid, date=epistart) %>% mutate(source="hes"))) %>%
  union_all((raw_ckd5_code_opcs4 %>% select(patid, date=evdate) %>% mutate(source="hes"))) %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  #filter(date>=min_dob & ((source=="gp" & date<=gp_ons_end_date) | (source=="hes" & (is.na(gp_ons_death_date) | date<=gp_ons_death_date)))) %>% ## as above - ONS variables substituted
  filter(date>=min_dob & ((source=="gp" & date<=gp_end_date) | (source=="hes" & (is.na(gp_end_date) | date<=gp_end_date)))) %>%
  group_by(patid) %>%
  summarise(first_test_date=min(date, na.rm=TRUE))%>%
  ungroup() %>%
  analysis$cached("earliest_clean_ckd5",indexes=c("patid", "first_test_date"))


## Combine CKD5 and other codes

ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  select(patid, ckd_stage, first_test_date) %>%
  union_all(earliest_clean_ckd5 %>% mutate(ckd_stage="stage_5")) %>%
  analysis$cached("ckd_stages_from_algorithm_interim_3",indexes=c("patid","ckd_stage","first_test_date"))

ckd_stages_from_algorithm %>% count()        
#1,188,059


################################################################################################################################

# Define date of onset for each stage

## For each person, define date of onset of each stage (earliest incident) - assume no returning to less severe stages

ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  group_by(patid, ckd_stage) %>%
  summarise(ckd_stage_start=min(first_test_date, na.rm=TRUE)) %>%
  ungroup()


## Remove where start date of less severe stage is later than start date of more severe stage
### Reshape wide first

ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  pivot_wider(id_cols=patid,
              names_from=ckd_stage,
              values_from=ckd_stage_start) %>%
  mutate(stage_1=ifelse(!is.na(stage_1) & !is.na(stage_2) & stage_1>stage_2, NA, stage_1),
         stage_1=ifelse(!is.na(stage_1) & !is.na(stage_3a) & stage_1>stage_3a, NA, stage_1),
         stage_1=ifelse(!is.na(stage_1) & !is.na(stage_3b) & stage_1>stage_3b, NA, stage_1),
         stage_1=ifelse(!is.na(stage_1) & !is.na(stage_4) & stage_1>stage_4, NA, stage_1),
         stage_1=ifelse(!is.na(stage_1) & !is.na(stage_5) & stage_1>stage_5, NA, stage_1),
         stage_2=ifelse(!is.na(stage_2) & !is.na(stage_3a) & stage_2>stage_3a, NA, stage_2),
         stage_2=ifelse(!is.na(stage_2) & !is.na(stage_3b) & stage_2>stage_3b, NA, stage_2),
         stage_2=ifelse(!is.na(stage_2) & !is.na(stage_4) & stage_2>stage_4, NA, stage_2),
         stage_2=ifelse(!is.na(stage_2) & !is.na(stage_5) & stage_2>stage_5, NA, stage_2),
         stage_3a=ifelse(!is.na(stage_3a) & !is.na(stage_3b) & stage_3a>stage_3b, NA, stage_3a),
         stage_3a=ifelse(!is.na(stage_3a) & !is.na(stage_4) & stage_3a>stage_4, NA, stage_3a),
         stage_3a=ifelse(!is.na(stage_3a) & !is.na(stage_5) & stage_3a>stage_5, NA, stage_3a),
         stage_3b=ifelse(!is.na(stage_3b) & !is.na(stage_4) & stage_3b>stage_4, NA, stage_3b),
         stage_3b=ifelse(!is.na(stage_3b) & !is.na(stage_5) & stage_3b>stage_5, NA, stage_3b),
         stage_4=ifelse(!is.na(stage_4) & !is.na(stage_5) & stage_4>stage_5, NA, stage_4)) %>%
  analysis$cached("ckd_stages_from_algorithm", unique_indexes="patid")

ckd_stages_from_algorithm %>% count()        
#1,095,842

# keep CKD stage at index date

analysis = cprd$analysis(analysis_prefix)


ckd_stages_2020 <- ckd_stages_from_algorithm %>%
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
  
  filter(!preckdstage == "stage_1" & !preckdstage == "stage_5") %>%
  
  analysis$cached("ckd_stages_2020_im", unique_indexes="patid")

## combine with ethnicity and basic data

analysis = cprd$analysis("all")
ethnicity <- ethnicity %>% analysis$cached("gp_ethnicity")
language <- language %>% analysis$cached("gp_language")

analysis = cprd$analysis(analysis_prefix)
comorbidities <- comorbidities %>% analysis$cached("comorbidities_2020_jul", unique_indexes="patid")
#medications <- medications %>% analysis$cached("medications_2020_jul", unique_indexes = "patid")

#### join up and cache
analysis = cprd$analysis(analysis_prefix)
ckd_stages_2020 <- ckd_stages_2020 %>%
  select(patid, preckdstage, preckdstagedate) %>%
  inner_join((dob %>% select(patid, dob)), by="patid") %>%
  left_join((full_egfr_index_date_merge %>% select(patid, preegfr, creatinine, preegfrdate)), by = "patid") %>%
  left_join((cprd$tables$patient %>% select(patid, gender, regstartdate, regenddate, pracid, cprd_ddate)), by="patid") %>%
  left_join((cprd$tables$practice %>% select(pracid, lcd, region)), by="pracid") %>%
  ## I have disabled the lines in the script referring to IMD, HES, and ONS (linkage not available)
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
  
  select(patid, gender, dob, pracid, prac_region=region, 
         preckdstage, preckdstagedate, 
         preegfr, creatinine, preegfrdate, 
         gp_5cat_ethnicity, gp_16cat_ethnicity, gp_qrisk2_ethnicity, 
         language_cat=gp_language_cat,
#         imd2015_10, 
         regstartdate, gp_record_end, 
#         death_date, with_hes
         ) %>% 
  
  analysis$cached("ckd_cohort_temp3", unique_indexes="patid", indexes=c("gender", "dob", "preckdstage"))

ckd_cohort <- ckd_stages_2020 %>%
  inner_join(comorbidities, by="patid") %>%
  #inner_join(medications, by="patid") %>%
  analysis$cached("ckd_cohort_temp4", unique_indexes="patid", indexes=c("gender", "dob", "preckdstage"))

