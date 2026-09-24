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

# create table for ckd stage 5 (by diagnostic codes) and ckd_stages_from_algorithm (by eGFR/ACR algorithm) - to be used for cohort definition
analysis = cprd$analysis("all_patid")
comorbids = "ckd5_code"

for (i in comorbids) {
  if (length(codes[[i]]) > 0) {
    print(paste("making", i, "medcode table"))
    
    raw_tablename <- paste0("raw_", i, "_medcodes")
    
    data <- cprd$tables$observation %>%
      inner_join(codes[[i]], by="medcodeid") %>%
      analysis$cached(raw_tablename, indexes=c("patid", "obsdate"))
    
    assign(raw_tablename, data)
    
  }
  
  if (i!="hypertension" && length(codes[[paste0("icd10_", i)]]) > 0) {
    print(paste("making", i, "ICD10 code table"))
    
    raw_tablename <- paste0("raw_", i, "_icd10")
    
    data <- cprd$tables$hesDiagnosisEpi %>%
      inner_join(codes[[paste0("icd10_",i)]], sql_on="LHS.ICD LIKE CONCAT(icd10,'%')") %>%
      analysis$cached(raw_tablename, indexes=c("patid", "epistart"))
    
    assign(raw_tablename, data)
    
  }
  
  if (length(codes[[paste0("opcs4_", i)]]) > 0) {
    print(paste("making", i, "OPCS4 code table"))
    
    raw_tablename <- paste0("raw_", i, "_opcs4")
    
    data <- cprd$tables$hesProceduresEpi %>%
      inner_join(codes[[paste0("opcs4_",i)]], sql_on="LHS.OPCS LIKE CONCAT(opcs4,'%')") %>%
      analysis$cached(raw_tablename, indexes=c("patid", "evdate"))
    
    assign(raw_tablename, data)
    
  }
}

############################################################################################

# get creatinine_blood biomarker table (for eGFR calculation) and ACR biomarker table (for albuminuria calculation)

biomarkers <- c("creatinine_blood", "acr", "albumin_urine", "creatinine_urine")

for (i in biomarkers) {
  
  print(paste0("Making raw ", i, " biomarker table"))
  
  raw_tablename <- paste0("raw_", i, "_medcodes")
  
  data <- cprd$tables$observation %>%
    inner_join(codes[[i]], by="medcodeid") %>%
    analysis$cached(raw_tablename, indexes=c("patid", "obsdate", "testvalue", "numunitid"))
  
  
  assign(raw_tablename, data)
  
}

for (i in biomarkers) {
  
  print(paste0("Cleaning ", i, " biomarker table"))
  
  raw_tablename <- paste0("raw_", i, "_medcodes")
  clean_tablename <- paste0("clean_", i, "_medcodes")
  
  
    raw_data <- get(raw_tablename)
 
  
  # select valid numunitid
  
  if (i=="albumin_urine") {
    data <- raw_data %>%
      filter(numunitid==183)
  }
  else if (i=="creatinine_urine") {
    data <- raw_data %>%
      filter(numunitid==218 | numunitid==285) %>%
      mutate(testvalue=ifelse(numunitid==285, testvalue/1000, testvalue))
  }
  else {
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
#40,051,589

ckd_stages_from_algorithm %>% summarise(total=sum(test_count, na.rm=TRUE))
#total number of tests: 112,440,442 as above


### C) Remove periods with 1 reading, or with multiple readings but <90 days between first and last test, and cache

ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>%
  filter(test_count>1 & datediff(last_test_date, first_test_date)>=90) %>%
  analysis$cached("ckd_stages_from_algorithm_interim_2",indexes=c("patid","ckd_stage","first_test_date"))

ckd_stages_from_algorithm %>% count()
#17,886,174


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


#################################################################################################################################

# get cohort ids for all ckd stages (1-5) and advanced ckd (stages 4-5)
analysis = cprd$analysis(analysis_prefix)

ckd_ids <- ckd_stages_from_algorithm %>% 
  filter(!(is.na(stage_1) & is.na(stage_2) & is.na(stage_3a) & 
             is.na(stage_3b) & is.na(stage_4) & is.na(stage_5))) %>%
  mutate(
    first_ckd_date = as.Date(
      pmin(
        ifelse(is.na(stage_1), as.Date("2050-01-01"), stage_1),
        ifelse(is.na(stage_2), as.Date("2050-01-01"), stage_2),
        ifelse(is.na(stage_3a), as.Date("2050-01-01"), stage_3a),
        ifelse(is.na(stage_3b), as.Date("2050-01-01"), stage_3b),
        ifelse(is.na(stage_4), as.Date("2050-01-01"), stage_4),
        ifelse(is.na(stage_5), as.Date("2050-01-01"), stage_5),
        na.rm = TRUE
      )
    )
  ) %>%
  mutate(first_ckd_date = ifelse(first_ckd_date == as.Date("2050-01-01"), NA, first_ckd_date)) %>%
  group_by(patid) %>%
  dbplyr::window_order(first_ckd_date) %>%
  distinct(patid, .keep_all = TRUE) %>%
  ungroup() %>%
  distinct(patid, .keep_all = TRUE) %>%
  select(-contains("stage"), -confirmed_acr3_date) %>%
  analysis$cached("ckd_ids_im", unique_indexes="patid")

ckd_ids %>% count() #1452649

ckd_ids %>% anti_join(practice_exclusion_ids, by="patid") %>% anti_join(gender_exclusion_ids, by="patid") %>% count() #2110415


ckd_ids <- ckd_ids %>%
  anti_join(practice_exclusion_ids, by="patid") %>% 
  anti_join(gender_exclusion_ids, by="patid") %>%
  analysis$cached("ckd_ids", unique_indexes="patid")

ckd_ids %>% count() #1452649

## create table for ids with advanced ckd only (ckd stages 4 or 5)
analysis = cprd$analysis("rk")

advanced_ckd_ids <- ckd_stages_from_algorithm %>% 
  filter(!(is.na(stage_4) & is.na(stage_5))) %>%
  mutate(
    index_date = as.Date(
      pmin(
        ifelse(is.na(stage_4), as.Date("2050-01-01"), stage_4),
        ifelse(is.na(stage_5), as.Date("2050-01-01"), stage_5),
        na.rm = TRUE
      )
    )
  ) %>%
  mutate(index_date = ifelse(index_date == as.Date("2050-01-01"), NA, index_date)) %>%
  group_by(patid) %>%
  dbplyr::window_order(index_date) %>%
  distinct(patid, .keep_all = TRUE) %>%
  ungroup() %>%
  distinct(patid, .keep_all = TRUE) %>%
  select(-contains("stage"), -confirmed_acr3_date) %>%
  analysis$cached("advanced_ckd_ids_im", unique_indexes="patid")

advanced_ckd_ids %>% count() 

advanced_ckd_ids %>% anti_join(practice_exclusion_ids, by="patid") %>% anti_join(gender_exclusion_ids, by="patid") %>% count()


advanced_ckd_ids <- advanced_ckd_ids %>%
  anti_join(practice_exclusion_ids, by="patid") %>% 
  anti_join(gender_exclusion_ids, by="patid") %>%
  analysis$cached("advanced_ckd_ids", unique_indexes="patid")

advanced_ckd_ids %>% count()


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

ckd_cohort <- ckd_ids %>%
  left_join(dob, by="patid") %>%
  left_join((cprd$tables$patient %>% select(patid, gender, regenddate, pracid)), by="patid") %>%
  left_join((cprd$tables$practice %>% select(pracid, lcd, region)), by="pracid") %>%
  left_join((cprd$tables$onsDeath %>% select(patid, reg_date_of_death)), by="patid") %>%
  left_join((cprd$tables$patientImd %>% select(patid, imd_decile)), by="patid") %>%
  left_join((cprd$tables$validDateLookup %>% select(patid, gp_end_date)), by="patid") %>%
  left_join((cprd$tables$patidsWithLinkage %>% mutate(with_hes=1L) %>% select(patid, with_hes, hes_end_date)), by="patid") %>%
  mutate(with_hes=ifelse(is.na(with_hes), 0L, 1L)) %>%
  left_join(ethnicity, by="patid") %>%
  select(patid, gender, dob, pracid, prac_region=region, ethnicity_5cat, ethnicity_16cat, ethnicity_qrisk2, imd_decile, regstartdate, gp_end_date, death_date=reg_date_of_death, with_hes, hes_end_date, first_ckd_date) %>%
  analysis$cached("ckd_cohort", unique_indexes="patid", indexes=c("gender", "dob"))
                  
                  
ckd_cohort %>% count() # 1,452,649

# do similar for advanced ckd cohort
analysis = cprd$analysis("rk")

advanced_ckd_cohort <- advanced_ckd_ids %>%
  left_join(dob, by="patid") %>%
  left_join((cprd$tables$patient %>% select(patid, gender, regenddate, pracid)), by="patid") %>%
  left_join((cprd$tables$practice %>% select(pracid, lcd, region)), by="pracid") %>%
  left_join((cprd$tables$onsDeath %>% select(patid, reg_date_of_death)), by="patid") %>%
  left_join((cprd$tables$patientImd %>% select(patid, imd_decile)), by="patid") %>%
  left_join((cprd$tables$validDateLookup %>% select(patid, gp_end_date)), by="patid") %>%
  left_join((cprd$tables$patidsWithLinkage %>% mutate(with_hes=1L) %>% select(patid, with_hes, hes_end_date)), by="patid") %>%
  mutate(with_hes=ifelse(is.na(with_hes), 0L, 1L)) %>%
  left_join(ethnicity, by="patid") %>%
  select(patid, gender, dob, pracid, prac_region=region, ethnicity_5cat, ethnicity_16cat, ethnicity_qrisk2, imd_decile, regstartdate, gp_end_date, death_date=reg_date_of_death, with_hes, hes_end_date, index_date) %>%
  analysis$cached("advanced_ckd_cohort", unique_indexes="patid", indexes=c("gender", "dob"))

  advanced_ckd_cohort %>% count()

############################################################################################
