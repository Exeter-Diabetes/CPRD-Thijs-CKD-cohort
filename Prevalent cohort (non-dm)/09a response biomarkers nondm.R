
# Kidney/AKI outcome ascertainment for the non-DM CKD sequential trial cohort, at each
# 3-monthly trial index date (2022-03-01 to 2024-03-01).
#
# For each index date, finds (all relative to that index date, using post-index-date
# measurements/codes only):
#   - Date of 40%/50% decl in eGFR from baseline (both raw and confirmed by a second
#     qualifying measurement >=28 days later)
#   - Acute kidney injury (AKI), defined as the EARLIEST of:
#       (a) lab-defined AKI per KDIGO creatinine criteria (see below), using creatinine
#           readings up to 2 years post-index date
#       (b) a relevant ICD10 AKI code post-index date (post_index_date_first_aki, from
#           "04 comorbidities nondm.R" - HES-only, see that script for codelist details)
#     aki_stage / aki_stage_max / aki_date_max are LAB-DERIVED only (ICD10 codes carry no
#     KDIGO stage information).


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

# Kidney outcomes

## eGFR

### Add in next eGFR measurement

analysis = cprd$analysis("all_patid")

egfr_long <- egfr_long %>% analysis$cached("clean_egfr_medcodes")

creatinine_long <- creatinine_long %>% analysis$cached("clean_creatinine_blood_medcodes")

# (3-monthly required for sequential trial emulation of SGLT2i in non-DM CKD)
dates <- unique(c(
  seq(from = as.Date("2022-03-01"), to = as.Date("2024-03-01"), by = "3 months")
))


date_strings <- format(dates, "%Y-%m-%d")


for (d in date_strings) {
  print(d)
  index_date <- as.Date(d)

  # Clear stale handles from the previous iteration so a missing per-date table errors
  # loudly instead of silently re-caching the previous date's data under this date's name
  # (baseline_biomarkers/comorbidities_d are read via analysis$cached() before being reassigned below)
  rm(list = intersect(c("baseline_biomarkers", "comorbidities_d"), ls()))

analysis = cprd$analysis(paste0(analysis_prefix, "_", d))

baseline_biomarkers <- baseline_biomarkers %>% analysis$cached("baseline_biomarkers")

# next_egfr <- baseline_biomarkers %>%
#  select(patid, preegfrdate) %>%
#  left_join(egfr_long, by="patid") %>%
#  filter(datediff(date, preegfrdate)>0) %>%
#  group_by(patid) %>%
#  summarise(next_egfr_date=min(date, na.rm=TRUE)) %>%
#  analysis$cached("response_biomarkers_next_egfr", indexes=c("patid"))


### Add in 40% decl in eGFR outcome: both based on single measurement and where confirmed by a second measurement at least 28 days later
### Join drug start dates with all longitudinal eGFR measurements, and only keep later eGFR measurements which are <=60% of the baseline value
### Checked and those with null eGFR do get dropped
egfr_40_dates <- baseline_biomarkers %>%
  select(patid, preegfr, preegfrdate) %>%
  left_join(egfr_long, by="patid") %>%
  filter(datediff(date, preegfrdate)>0 & testvalue<=0.6*preegfr) %>%
  rename(egfr_40_date=date) %>%
  analysis$cached("resp_biomarkers_egfr40_dates", indexes=c("patid", "egfr_40_date"))

egfr40_decl <- egfr_40_dates %>%
  group_by(patid) %>%
  summarise(egfr_40_decl_date=min(egfr_40_date, na.rm=TRUE),
            preegfr=max(preegfr, na.rm=TRUE)) %>%  # preegfr is constant per patid; max() used as first() is unsupported remotely
  ungroup() %>%
  analysis$cached("resp_biomarkers_egfr40_decl", indexes=c("patid", "egfr_40_decl_date"))

egfr40_decl_next_egfr <- egfr40_decl %>%
  left_join(egfr_long, by="patid") %>%
  filter(datediff(date, egfr_40_decl_date)>=28 & testvalue<=0.6*preegfr) %>% 
  group_by(patid) %>%
  summarise(next_egfr=min(date, na.rm=T)) %>% 
  ungroup() %>%
  analysis$cached("resp_biomarkers_egfr40_decl_next_egfr", indexes=c("patid", "next_egfr"))

egfr40_decl_confirmed <- egfr40_decl %>%
  inner_join(egfr40_decl_next_egfr, by="patid") %>%
  mutate(egfr_40_decl_date_confirmed=egfr_40_decl_date) %>%
  select(patid, egfr_40_decl_date_confirmed) %>%
  analysis$cached("resp_biomarkers_egfr40_decl_confirmed", indexes=c("patid"))

print(paste0("Finished eGFR 40% decline outcome ascertainment for index date ", d))
## AKI (acute kidney injury), per KDIGO creatinine-based criteria (up to 2 years after drug start date, but only post-drug start readings)
### Baseline creatinine = precreatinine_blood (already computed above: closest value -2 years to +7 days pre-drug start).
### Patients with no baseline creatinine in this window cannot have AKI assessed - left as NA rather than assumed non-AKI.
### Stage 1: >=26.5 umol/L rise within 48 hours (approximated as <=2 calendar days, since CPRD dates have no time component) OR >=1.5x baseline
### Stage 2: >=2.0x baseline
### Stage 3: >=3.0x baseline OR absolute creatinine >=353.6 umol/L (and rise of >=26.5 umol/L within 48 hours)
### NB: HES codes for AKI to be combined separately from comorbidities script.

aki_creatinine <- creatinine_long %>%
 mutate(datediff=datediff(date, index_date)) %>%
  filter(datediff > 0 & datediff < 730) %>%
  select(patid, date, testvalue) %>%
  analysis$cached("resp_biomarkers_aki_creatinine", indexes=c("patid", "date"))

  print(paste0("Finished creatinine extraction for index date ", d))

### For each post-index-date creatinine reading, find the lowest creatinine value in the
### preceding 48 hours (self-join within patid, since we need the minimum over a date range
### rather than just the adjacent reading)

aki_48h_pairs <- aki_creatinine %>%
  inner_join(
    (aki_creatinine %>% select(patid, earlier_date=date, earlier_testvalue=testvalue)),
    by=c("patid")
  ) %>%
  filter(earlier_date<=date & datediff(date, earlier_date)<=2) %>%
  analysis$cached("resp_biomarkers_aki_48h_pairs", indexes=c("patid", "date"))
  print(paste0("Finished 48h AKI pairwise comparison for index date ", d))

aki_48h_rise <- aki_48h_pairs %>%
  group_by(patid, date) %>%
  summarise(testvalue=max(testvalue, na.rm=TRUE),
            min_earlier_testvalue=min(earlier_testvalue, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(rise_48h=(testvalue-min_earlier_testvalue)>=26.5) %>%
  select(patid, date, rise_48h) %>%
  analysis$cached("resp_biomarkers_aki_48h_rise", indexes=c("patid", "date"))

  print(paste0("Finished kidney/AKI outcome ascertainment for index date ", d))

### Join with baseline creatinine and stage every qualifying postdrug reading

aki_staged <- aki_creatinine %>%
  left_join(aki_48h_rise, by=c("patid", "date")) %>%
  left_join((baseline_biomarkers %>% select(patid, precreatinine_blood)), by=c("patid")) %>%
  filter(!is.na(precreatinine_blood)) %>%
  mutate(rise_48h=coalesce(rise_48h, FALSE),
         creat_ratio=testvalue/precreatinine_blood,
         # Stage 3: >=3x baseline, OR (absolute creatinine >=353.6 AND a qualifying 48h rise)
         aki_stage=case_when(
           creat_ratio>=3 | (testvalue>=353.6 & rise_48h) ~ 3L,
           creat_ratio>=2 ~ 2L,
           creat_ratio>=1.5 | rise_48h ~ 1L,
           TRUE ~ NA_integer_
         )) %>%
  filter(!is.na(aki_stage)) %>%
  analysis$cached("resp_biomarkers_aki_staged", indexes=c("patid", "date"))

  print(paste0("Finished AKI staging for index date ", d))

### Earliest AKI event (any stage) and its stage; highest stage EVER reached and the
### (earliest) date that stage was first reached. Done via separate group_by/join steps
### rather than window functions, since aki_date_max must be the date THAT MATCHES
### aki_stage_max, not simply the latest date on record for the patient.

aki_first_event <- aki_staged %>%
  group_by(patid) %>%
  summarise(aki_date=min(date, na.rm=TRUE)) %>%
  ungroup()

aki_stage_at_first_event <- aki_first_event %>%
  inner_join(aki_staged, by="patid") %>%
  filter(date==aki_date) %>%
  group_by(patid, aki_date) %>%
  summarise(aki_stage=max(aki_stage, na.rm=TRUE)) %>%  # if >1 code on the same day, take the higher stage
  ungroup()

aki_stage_max <- aki_staged %>%
  group_by(patid) %>%
  summarise(aki_stage_max=max(aki_stage, na.rm=TRUE)) %>%
  ungroup()

aki_date_of_max_stage <- aki_stage_max %>%
  inner_join(aki_staged, by="patid") %>%
  filter(aki_stage==aki_stage_max) %>%
  group_by(patid) %>%
  summarise(aki_date_max=min(date, na.rm=TRUE)) %>%  # earliest date the max stage was first reached
  ungroup()

aki_lab_outcome <- aki_stage_at_first_event %>%
  inner_join(aki_stage_max, by="patid") %>%
  left_join(aki_date_of_max_stage, by="patid") %>%
  select(patid, aki_date, aki_stage, aki_stage_max, aki_date_max) %>%
  analysis$cached("resp_biomarkers_aki_lab_outcome", indexes=c("patid"))

  print(paste0("Finished lab-based AKI outcome ascertainment for index date ", d))
### Combine with ICD10-coded AKI (HES only - see 04 comorbidities nondm.R). Table
### "comorbidities" for this index date is cached under the SAME analysis prefix (ckd_{d}),
### and post_index_date_first_aki is already restricted to post-index-date occurrences.
comorbidities_d <- comorbidities_d %>% analysis$cached("comorbidities") %>%
  select(patid, aki_icd10_date=post_index_date_first_aki)

# MySQL has no FULL OUTER JOIN (and dbplyr's full_join() can't translate it remotely), so
# emulate it as left_join + (anti_join of the unmatched RHS rows, with LHS-only columns
# filled as NA) + union_all, matching columns/types between the two halves.
aki_outcome_matched <- aki_lab_outcome %>%
  left_join(comorbidities_d, by="patid")

aki_outcome_icd10_only <- comorbidities_d %>%
  anti_join(aki_lab_outcome, by="patid") %>%
  mutate(aki_date=as.Date(NA), aki_stage=NA_integer_, aki_stage_max=NA_integer_, aki_date_max=as.Date(NA)) %>%
  select(patid, aki_date, aki_stage, aki_stage_max, aki_date_max, aki_icd10_date)

# Remote LEAST()-style min() propagates NULL if either side is NULL on MySQL, so replace
# NA with a distant sentinel date before combining, then revert the sentinel back to NA
# (same pattern used in "10 final merge nondm.R" for combining nullable dates remotely).
aki_outcome <- aki_outcome_matched %>%
  union_all(aki_outcome_icd10_only) %>%
  mutate(
    aki_date=pmin(
      ifelse(is.na(aki_date), as.Date("2050-01-01"), aki_date),
      ifelse(is.na(aki_icd10_date), as.Date("2050-01-01"), aki_icd10_date),
      na.rm=TRUE
    ),
    aki_date=ifelse(aki_date==as.Date("2050-01-01"), as.Date(NA), aki_date)
  ) %>%
  select(patid, aki_date, aki_stage, aki_stage_max, aki_date_max) %>%
  analysis$cached("resp_biomarkers_aki_outcome", indexes=c("patid"))

  print(paste0("Finished AKI outcome ascertainment for index date ", d))
### Add in 50% decl in eGFR outcome: confirmed and not
egfr_50_dates <- baseline_biomarkers %>%
  select(patid, preegfr, preegfrdate) %>%
  left_join(egfr_long, by="patid") %>%
  filter(datediff(date, preegfrdate)>0 & testvalue<=0.5*preegfr) %>%
  rename(egfr_50_date=date) %>%
  analysis$cached("resp_biomarkers_egfr50_dates", indexes=c("patid", "egfr_50_date"))

egfr50_decl <- egfr_50_dates %>%
  group_by(patid) %>%
    summarise(egfr_50_decl_date=min(egfr_50_date, na.rm=TRUE),
            preegfr=max(preegfr, na.rm=TRUE)) %>%  # preegfr is constant per patid; max() used as first() is unsupported remotely
  ungroup() %>%
  analysis$cached("resp_biomarkers_egfr50_decl", indexes=c("patid"))

egfr50_decl_confirmed <- egfr_50_dates %>%
  select(-testvalue) %>%
  left_join(egfr_long, by="patid") %>%
  filter(datediff(date, egfr_50_date)>=28 & testvalue<=0.5*preegfr) %>%
  group_by(patid) %>%
  summarise(egfr_50_decl_date_confirmed=min(egfr_50_date, na.rm=TRUE)) %>%
  ungroup() %>%
  analysis$cached("resp_biomarkers_egfr50_decl_confirmed", indexes=c("patid"))

print(paste0("Finished eGFR 50% decline outcome ascertainment for index date ", d))

############################################################################################
#
## add number of egfr measurements in 12 months following baseline
# egfr_counts_12m <- baseline_biomarkers %>%
#  select(patid, drug_substance, dstartdate, preegfrdate) %>%
#  left_join(egfr_long, by = "patid") %>%
#  # Keep only measurements after the baseline (preegfrdate)
#  filter(datediff(date, preegfrdate) > 0,
#         datediff(date, preegfrdate) <= 365) %>%   # within 12 months
#  group_by(patid, drug_substance, dstartdate) %>%
#  summarise(
#    egfr_count_12m = n(),   # count how many eGFR measurements
#    .groups = "drop"
#  ) %>%
#  analysis$cached("response_biomarkers_egfr_count_12m",
#                  indexes = c("patid", "dstartdate", "drug_substance"))

# preegfr_counts_12m <- baseline_biomarkers %>%
#  select(patid, drug_substance, dstartdate, preegfrdate) %>%
#  left_join(egfr_long, by = "patid") %>%
#  # Keep only measurements after the baseline (preegfrdate)
#  filter(datediff(date, preegfrdate) <= 0,
#         datediff(date, preegfrdate) >= -365) %>%   # within 12 months
#  group_by(patid, drug_substance, dstartdate) %>%
#  summarise(
#    preegfr_count_12m = n(),   # count how many eGFR measurements#
#    .groups = "drop"
#  ) %>%
#  analysis$cached("response_biomarkers_preegfr_count_12m",
#                  indexes = c("patid", "dstartdate", "drug_substance"))


############################################################################################

# Join tables together

# Initialise from the patid list (baseline_biomarkers is 1 row per patient at this index
# date) - the prior version of this script referenced response_biomarkers before it was
# ever assigned, which only works once the final cached table already exists.
response_biomarkers <- baseline_biomarkers %>%
  select(patid) %>%
#  left_join(next_egfr, by=c("patid")) %>%
  left_join(egfr40_decl, by=c("patid")) %>%
  left_join(egfr40_decl_next_egfr, by=c("patid")) %>%
  left_join(egfr40_decl_confirmed, by=c("patid")) %>%
    analysis$cached("response_biomarkers_im", indexes=c("patid"))

response_biomarkers <- response_biomarkers %>%
  left_join(egfr50_decl, by=c("patid")) %>%
  left_join(egfr50_decl_confirmed, by=c("patid")) %>%
  left_join(aki_outcome, by=c("patid")) %>%
  select(-contains("preegfr")) %>%
    analysis$cached("response_biomarkers", indexes=c("patid"))

}
