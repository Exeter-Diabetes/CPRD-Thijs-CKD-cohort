# Kidney-protective medication sorting for the non-DM CKD cohort
# Drug classes tracked: RASi (ACE inhibitors + ARBs combined), SGLT2i, GLP1-RA
#
# Approach mirrors 01_mm_drug_sorting_and_combos.R from the MASTERMIND cohort,
# but simplified to class level only (no substance tracking needed) and restricted
# to kidney-protective medications relevant in a non-diabetic CKD population.
#
# A new drug episode starts if the gap from the previous prescription exceeds
# 183 days (6 months), or if there is no prior prescription.
#
# Raw prodcode tables (patid, date) are created by 05 medications nondm.R under
# the 'all_patid' MySQL prefix. analysis$cached is idempotent - if those tables
# already exist, no recomputation occurs.
#
# Output tables:
#   ckd_kp_meds_long  : one row per patid/date/drug_class with start/stop flags
#   ckd_kp_meds       : one row per patid/date (wide), with:
#                         on_RASi, on_SGLT2, on_GLP1     (1 = currently on that class)
#                         dstart_RASi, dstop_RASi, ...   (episode start/stop flags)
#                         kp_combo                        (e.g. "RASi", "RASi_SGLT2")
#                         n_kp_classes                    (number of classes currently on)

############################################################################################

# Setup
library(tidyverse)
library(aurum)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "nondiabetes-jun2024", cprdConf = "C:/Users/tj358/OneDrive - University of Exeter/CPRD/aurum.yaml")

analysis = cprd$analysis("ckd")


############################################################################################

# Get handles to raw prodcode tables (created by 05 medications nondm.R)

analysis_all_patid <- cprd$analysis("all_patid")

raw_ace_inhibitors_prodcodes <- raw_ace_inhibitors_prodcodes %>% analysis_all_patid$cached("raw_ace_inhibitors_prodcodes")
# raw_ace_inhibitors_prodcodes %>% count()

raw_arb_prodcodes <- raw_arb_prodcodes %>% analysis_all_patid$cached("raw_arb_prodcodes")
# raw_arb_prodcodes %>% count()

raw_sglt2_prodcodes <- raw_sglt2_prodcodes %>% analysis_all_patid$cached("raw_sglt2_prodcodes")
# raw_sglt2_prodcodes %>% count()

raw_glp1_prodcodes <- raw_glp1_prodcodes %>% analysis_all_patid$cached("raw_glp1_prodcodes")
# raw_glp1_prodcodes %>% count()


############################################################################################

# Build long table: one row per patid / prescription date / drug class
# ACE inhibitors and ARBs are merged into a single 'RASi' class

kp_meds_long <- raw_ace_inhibitors_prodcodes %>%
  mutate(drug_class = "RASi") %>%
  union_all(
    raw_arb_prodcodes %>%
      mutate(drug_class = "RASi")
  ) %>%
  union_all(
    raw_sglt2_prodcodes %>%
      mutate(drug_class = "SGLT2")
  ) %>%
  union_all(
    raw_glp1_prodcodes %>%
      mutate(drug_class = "GLP1")
  ) %>%
  # Apply date validity filters
  inner_join(cprd$tables$validDateLookup, by = "patid") %>%
  filter(date >= min_dob & date <= gp_end_date) %>%
  select(patid, date, drug_class) %>%
  # Remove same-day same-class duplicates (e.g. ACEi + ARB on same day both become RASi)
  distinct() %>%
  analysis$cached("kp_meds_long_interim_1", indexes = c("patid", "date"))

# kp_meds_long %>% count()


############################################################################################

# Define episode start and stop dates per drug class
# dstart = 1 if no prior prescription, or gap from previous > 183 days
# dstop  = 1 if no next prescription,  or gap to next     > 183 days

kp_meds_long <- kp_meds_long %>%
  group_by(patid, drug_class) %>%
  dbplyr::window_order(date) %>%
  mutate(
    dnextuse = datediff(lead(date), date),
    dprevuse = datediff(date, lag(date)),
    dstart   = dprevuse > 183 | is.na(dprevuse),
    dstop    = dnextuse > 183 | is.na(dnextuse)
  ) %>%
  ungroup() %>%
  analysis$cached("kp_meds_long_interim_2", indexes = c("patid", "date"))


# Number of classes starting / stopping on each prescription date
kp_meds_long <- kp_meds_long %>%
  group_by(patid, date) %>%
  mutate(
    numstart = sum(dstart, na.rm = TRUE),
    numstop  = sum(dstop,  na.rm = TRUE)
  ) %>%
  ungroup() %>%
  analysis$cached("kp_meds_long", indexes = c("patid", "date"))

# kp_meds_long %>% count()


############################################################################################

# Reshape wide: one row per patid / date, with dstart/dstop columns per class

kp_meds_wide <- kp_meds_long %>%
  pivot_wider(
    id_cols     = c(patid, date, numstart, numstop),
    names_from  = drug_class,
    values_from = c(dstart, dstop),
    values_fill = list(dstart = FALSE, dstop = FALSE)
  ) %>%
  analysis$cached("kp_meds_wide_interim", indexes = c("patid", "date"))

# kp_meds_wide %>% count()


############################################################################################

# Determine whether patient is currently on each drug class at each prescription date
# Uses cumulative sum of episode starts minus stops (same logic as MASTERMIND)

kp_meds <- kp_meds_wide %>%
  group_by(patid) %>%
  dbplyr::window_order(date) %>%
  mutate(
    on_RASi  = cumsum(dstart_RASi)  > cumsum(dstop_RASi)  | dstart_RASi  == 1 | dstop_RASi  == 1,
    on_SGLT2 = cumsum(dstart_SGLT2) > cumsum(dstop_SGLT2) | dstart_SGLT2 == 1 | dstop_SGLT2 == 1,
    on_GLP1  = cumsum(dstart_GLP1)  > cumsum(dstop_GLP1)  | dstart_GLP1  == 1 | dstop_GLP1  == 1,

    # Cumulative drug count (drug stopped on its dstop date still counts as 'on')
    cu_numstart  = cumsum(numstart),
    cu_numstop   = cumsum(numstop),
    numdrugs     = cu_numstart - cu_numstop + numstop
  ) %>%
  ungroup() %>%
  analysis$cached("kp_meds_interim_2", indexes = c("patid", "date"))


# Build combination string and direct drug count; drop intermediate columns
# Empty string (not NA) used in paste0 to avoid MySQL CONCAT NULL-propagation

kp_meds <- kp_meds %>%
  select(-c(starts_with("dstart"), starts_with("dstop"),
            cu_numstart, cu_numstop)) %>%
  mutate(
    kp_combo = paste0(
      ifelse(on_RASi  == 1, "RASi_",  ""),
      ifelse(on_SGLT2 == 1, "SGLT2_", ""),
      ifelse(on_GLP1  == 1, "GLP1_",  "")
    ),
    # Strip trailing underscore
    kp_combo = ifelse(
      str_sub(kp_combo, -1, -1) == "_",
      str_sub(kp_combo, 1, -2),
      kp_combo
    ),
    # Empty string (no active class on this date) → NA
    kp_combo     = ifelse(kp_combo == "", NA_character_, kp_combo),
    n_kp_classes = as.integer(on_RASi) + as.integer(on_SGLT2) + as.integer(on_GLP1)
  ) %>%
  analysis$cached("kp_meds", indexes = c("patid", "date"))

# kp_meds %>% count()


############################################################################################

# Sense-check: n_kp_classes should match numdrugs on every row
# kp_meds %>% filter(numdrugs != n_kp_classes | is.na(numdrugs) | is.na(n_kp_classes)) %>% count()
# expected: 0

# Distribution of combinations
# kp_meds %>% count(kp_combo) %>% collect() %>% arrange(desc(n))
