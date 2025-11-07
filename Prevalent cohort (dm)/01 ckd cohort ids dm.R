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

###############################################################################################

# load ckd stages based on egfr only

analysis = cprd$analysis("all_patid")

ckd_stages_from_algorithm <- ckd_stages_from_algorithm %>% analysis$cached("ckd_stages_from_algorithm")

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
  )  %>% 
  filter(!(is.na(stage_1) & is.na(stage_2) & is.na(stage_3a) & 
             is.na(stage_3b) & is.na(stage_4) & is.na(stage_5))) %>%
  analysis$cached("ckd_stages_from_algorithm_with_acr",
                  indexes = c("patid"))

# create list of ids of people with ckd and combine with type 2 diabetes cohort

analysis = cprd$analysis("diabetes_cohort")
practice_exclusion_ids <- practice_exclusion_ids %>% analysis$cached("practice_exclusion_ids")
gender_exclusion_ids <- gender_exclusion_ids %>% analysis$cached("gender_exclusion_ids")

ckd_ids <- ckd_stages_from_algorithm %>%
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
  select(-contains("stage"), -confirmed_acr3_date) %>% 
  anti_join(practice_exclusion_ids, by="patid") %>% anti_join(gender_exclusion_ids, by="patid") %>% 
  analysis$cached("ckd_ids", unique_indexes="patid")

ckd_ids %>% count() # 913743

# load diabetes cohort
analysis = cprd$analysis("all")

diabetes_cohort <- diabetes_cohort %>% 
  analysis$cached("diabetes_cohort")

diabetes_ckd_cohort <- diabetes_cohort %>%
  filter(diabetes_type == "type 2") %>%
  inner_join(ckd_ids, by = "patid") %>%
  analysis$cached("diabetes_ckd_cohort", unique_indexes = "patid")

diabetes_ckd_cohort  %>% count() # 797526