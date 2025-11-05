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
analysis = cprd$analysis(analysis_prefix)

############################################################################################

# get egfr
analysis = cprd$analysis("all_patid")

clean_egfr_medcodes <- clean_egfr_medcodes %>%
   analysis$cached("clean_egfr_medcodes", indexes=c("patid", "date", "testvalue"))

clean_egfr_medcodes %>% count()
#109,785,372 - lose readings of people with sex == NA or with missing creatinine


################################################################################################################################

# Convert eGFR to CKD stage

ckd_stages_from_all_egfr <- clean_egfr_medcodes %>%
  rename(egfr = testvalue) %>%
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
         ckd_stage=ckd_stage.x,
         egfr=egfr.x) %>%
  select(patid, patid_row_id, ckd_stage, ckd_start, ckd_end, egfr)


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
            maximum_date=max(ckd_end,na.rm=TRUE),
            test_count=max(patid_row_id, na.rm=TRUE)-min(patid_row_id, na.rm=TRUE)+1) %>%
  ungroup() %>%
  analysis$cached("ckd_stages_from_algorithm_interim_1",indexes=c("patid", "ckd_stage", "test_count", "first_test_date", "last_test_date"))

ckd_stages_from_algorithm %>% count()
#39,323,495

ckd_stages_from_algorithm %>% summarise(total=sum(test_count, na.rm=TRUE))
#total number of tests: 109,785,372 as above


### C) Remove periods with 1 reading, or with multiple readings but <90 days between first and last test, and cache

ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  filter(test_count>1 & datediff(last_test_date, first_test_date)>=90) %>%
  analysis$cached("ckd_stages_from_algorithm_interim_2",indexes=c("patid","ckd_stage","first_test_date"))

ckd_stages_from_algorithm %>% count()
#17,886,174


################################################################################################################################

# Combine with CKD5 medcodes/ICD10/OPCS4 codes

## Get raw CKD5 codes and clean
analysis = cprd$analysis("all_patid")
### Medcodes
raw_ckd5_code_medcodes <- raw_ckd5_medcodes %>% analysis$cached("raw_ckd5_code_medcodes")

### ICD10 codes
raw_ckd5_code_icd10 <- raw_ckd5_icd10 %>% analysis$cached("raw_ckd5_code_icd10")

### OPCS4 codes
raw_ckd5_code_opcs4 <- raw_ckd5_opcs4 %>% analysis$cached("raw_ckd5_code_opcs4")

## Clean, find earliest date per person, and re-cache

earliest_clean_ckd5 <- raw_ckd5_code_medcodes %>%
  select(patid, date=obsdate) %>%
  mutate(source="gp") %>%
  union_all((raw_ckd5_code_icd10 %>% select(patid, date=epistart) %>% mutate(source="hes"))) %>%
  union_all((raw_ckd5_code_opcs4 %>% select(patid, date=evdate) %>% mutate(source="hes"))) %>%
  inner_join(cprd$tables$validDateLookup, by="patid") %>%
  #filter(date>=min_dob & ((source=="gp" & date<=gp_ons_maximum_date) | (source=="hes" & (is.na(gp_ons_death_date) | date<=gp_ons_death_date)))) %>% ## as above - ONS variables substituted
  filter(date>=min_dob & ((source=="gp" & date<=gp_end_date) | (source=="hes" & (is.na(gp_end_date) | date<=gp_end_date)))) %>%
  group_by(patid) %>%
  summarise(first_test_date=min(date, na.rm=TRUE)) %>%
  ungroup() %>%
  analysis$cached("earliest_clean_ckd5",indexes=c("patid", "first_test_date"))


## Combine CKD5 and other codes

ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  select(patid, ckd_stage, first_test_date) %>%
  union_all(earliest_clean_ckd5 %>% mutate(ckd_stage="stage_5")) %>%
  analysis$cached("ckd_stages_from_algorithm_interim_3",indexes=c("patid","ckd_stage","first_test_date"))

ckd_stages_from_algorithm %>% count()        
#12,130,677


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
  analysis$cached("ckd_stages_from_algorithm_interim_4", unique_indexes="patid")

ckd_stages_from_algorithm %>% count()        
#8,466,065

################################################################################################################################

# load in acr data
clean_acr_medcodes <- clean_acr_medcodes %>%
  analysis$cached("clean_acr_medcodes", indexes=c("patid", "date", "testvalue"))

clean_acr_from_separate_medcodes <- clean_acr_from_separate_medcodes %>%
  analysis$cached("clean_acr_from_separate_medcodes", indexes=c("patid", "date", "testvalue"))

all_acr <- clean_acr_medcodes %>%
  select(patid, date, testvalue) %>%
  union_all(clean_acr_from_separate_medcodes %>%
              select(patid, date, testvalue))

# select those with acr >=3 mg/mmol
acr_high <- all_acr %>%
  filter(testvalue >= 3)

acr_span <- acr_high %>%
  group_by(patid) %>%
  summarise(
    min_date = min(date, na.rm = TRUE),
    max_date = max(date, na.rm = TRUE),
    n_tests = n()
  )

# confirm 2 readings 3 months apart or longer
confirmed_acr3 <- acr_span %>%
  filter(n_tests >= 2 & datediff(max_date, min_date) >= 90) %>%
  mutate(confirmed_acr3_date = min_date) %>%
  select(patid, confirmed_acr3_date) %>%
  analysis$cached("confirmed_acr3", indexes = c("patid", "confirmed_acr3_date"))

# join with ckd stage
ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  left_join(confirmed_acr3, by = "patid") %>%
  mutate(
    stage_1 = case_when(
      is.na(stage_1) ~ sql("NULL"),
      is.na(confirmed_acr3_date) ~ sql("NULL"),
      confirmed_acr3_date <= stage_1 ~ stage_1,
      confirmed_acr3_date > stage_1  ~ confirmed_acr3_date
    ),
    stage_2 = case_when(
      is.na(stage_2) ~ sql("NULL"),
      is.na(confirmed_acr3_date) ~ sql("NULL"),
      confirmed_acr3_date <= stage_2 ~ stage_2,
      confirmed_acr3_date > stage_2  ~ confirmed_acr3_date
    )
  ) %>%
  analysis$cached("ckd_stages_from_algorithm",
                  indexes = c("patid"))


######################################################################################
analysis = cprd$analysis(analysis_prefix)

# get dates at 6 month intervals
dates <- seq(from = as.Date("2019-03-01"),
             to   = as.Date("2024-03-01"),
             by   = "6 months")

date_strings <- format(dates, "%Y-%m-%d")


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