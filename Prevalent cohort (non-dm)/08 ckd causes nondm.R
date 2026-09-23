# Creates a patient-level table flagging non-standard CKD causes (one row per patid)
# Useful for excluding patients with specific CKD aetiology from analyses
# (e.g., excluding polycystic kidney disease or glomerulonephritis when predicting
#  long-term kidney failure risk)
#
# Causes covered: aav, adpkd, alport, antigbm, fabry, fsgs, gn_nos, igan,
#                 mcd, membranous, mpgn, other_pkd, sle
#
# Raw GP medcode and HES ICD-10 tables are shared with 04 comorbidities nondm.R:
# if those already exist (analysis$cached is idempotent), no recomputation occurs.
#
# Output table: ckd_ckd_causes (one row per patid)
#   - ckd_cause_{cause}          : 1 if patient has any ever-recorded code, 0 otherwise
#   - ckd_cause_earliest_{cause} : date of earliest recorded code (any source)
#   - any_pkd_cause              : 1 if adpkd or other_pkd
#   - any_gn_cause               : 1 if any glomerulonephritis / vasculitis cause
#   - any_nonstandard_ckd_cause  : 1 if any of the 13 causes above

############################################################################################

# Setup
library(tidyverse)
library(aurum)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "nondiabetes-jun2024", cprdConf = "C:/Users/tj358/OneDrive - University of Exeter/CPRD/aurum.yaml")

codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "01/06/2024")

analysis_prefix <- "ckd"

codepath <- "C:/Users/tj358/OneDrive - University of Exeter/CPRD/Aurum codelists/medcodes/ckd causes/"


############################################################################################

ckd_causes <- c(
  "aav", "adpkd", "alport", "antigbm", "fabry", "fsgs",
  "gn_nos", "igan", "mcd", "membranous", "mpgn", "other_pkd", "sle"
)


# Read codelists from local files
cause_codes <- list()

for (i in ckd_causes) {

  cause_codes[[i]] <- readr::read_tsv(
    paste0(codepath, "exeter_medcodelist_", i, ".tsv"),
    col_types = cols(.default = col_character())
  ) %>%
    rename(medcodeid = MedCodeId) %>%
    select(medcodeid)

  cause_codes[[paste0("icd10_", i)]] <- readr::read_tsv(
    paste0(codepath, "exeter_icd10_", i, ".txt"),
    col_types = cols(.default = col_character())
  ) %>%
    select(ICD10)

}


############################################################################################

# Pull raw code instances and cache under 'all_patid' prefix
# These may already exist from 04 comorbidities nondm.R - analysis$cached is idempotent

analysis <- cprd$analysis("all_patid")

for (i in ckd_causes) {

  # GP medcodes
  if (nrow(cause_codes[[i]]) > 0) {

    print(paste("making", i, "medcode table"))
    raw_tablename <- paste0("raw_", i, "_medcodes")

    data <- cprd$tables$observation %>%
      inner_join(cause_codes[[i]], by = "medcodeid", copy = TRUE) %>%
      analysis$cached(raw_tablename, indexes = c("patid", "obsdate"))

    assign(raw_tablename, data)
  }

  # HES ICD-10 codes
  if (nrow(cause_codes[[paste0("icd10_", i)]]) > 0) {

    print(paste("making", i, "ICD10 table"))
    raw_tablename <- paste0("raw_", i, "_icd10")

    data <- cprd$tables$hesDiagnosisEpi %>%
      inner_join(cause_codes[[paste0("icd10_", i)]], sql_on = "LHS.ICD LIKE CONCAT(ICD10,'%')", copy = TRUE) %>%
      analysis$cached(raw_tablename, indexes = c("patid", "epistart"))

    assign(raw_tablename, data)
  }

}


############################################################################################

# Create patient-level summary table: ever/never flag and earliest date per cause
# All codes are included (not restricted to pre-index date) because these are used
# for patient-level exclusions, not time-varying covariates

analysis <- cprd$analysis(analysis_prefix)

ckd_causes_table <- cprd$tables$patient %>%
  select(patid)


for (i in ckd_causes) {

  print(paste("Processing CKD cause:", i))

  interim_table    <- paste0("ckd_causes_im_", i)
  ever_flag_var    <- paste0("ckd_cause_", i)
  earliest_date_var <- paste0("ckd_cause_earliest_", i)

  medcode_tablename <- paste0("raw_", i, "_medcodes")
  icd10_tablename   <- paste0("raw_", i, "_icd10")


  # Combine GP and HES codes, applying date validity filters

  if (exists(medcode_tablename)) {

    gp_codes <- get(medcode_tablename) %>%
      inner_join(cprd$tables$validDateLookup, by = "patid") %>%
      filter(obsdate >= min_dob & obsdate <= gp_end_date) %>%
      select(patid, date = obsdate)

    all_codes <- gp_codes

    if (exists(icd10_tablename)) {

      hes_codes <- get(icd10_tablename) %>%
        inner_join(cprd$tables$validDateLookup, by = "patid") %>%
        filter(epistart >= min_dob & epistart <= as.Date("2023-03-31")) %>%
        select(patid, date = epistart)

      all_codes <- all_codes %>% union_all(hes_codes)

    }

  } else if (exists(icd10_tablename)) {

    all_codes <- get(icd10_tablename) %>%
      inner_join(cprd$tables$validDateLookup, by = "patid") %>%
      filter(epistart >= min_dob & epistart <= as.Date("2023-03-31")) %>%
      select(patid, date = epistart)

  } else {

    # No tables available for this cause - skip
    next

  }


  # Aggregate to patient level: earliest ever code date
  cause_summary <- all_codes %>%
    group_by(patid) %>%
    summarise(
      !!earliest_date_var := min(date, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    ungroup()

  earliest_date_var_sym <- sym(earliest_date_var)

  # Left join into main table and add binary ever/never flag
  ckd_causes_table <- ckd_causes_table %>%
    left_join(cause_summary, by = "patid") %>%
    mutate(
      !!ever_flag_var := as.integer(!is.na(!!earliest_date_var_sym))
    ) %>%
    analysis$cached(interim_table, indexes = "patid")

}


############################################################################################

# Add summary flags

ckd_causes_table <- ckd_causes_table %>%
  mutate(

    # Any PKD cause (autosomal dominant or other)
    any_pkd_cause = as.integer(
      ckd_cause_adpkd    == 1 |
      ckd_cause_other_pkd == 1
    ),

    # Any GN / vasculitis / immune-mediated cause, EXCLUDING IgAN
    # IgAN is treated separately because patients with IgAN may be included in most
    # analyses, unlike the other immune-mediated GN causes
    any_gn_cause = as.integer(
      ckd_cause_aav        == 1 |
      ckd_cause_alport     == 1 |
      ckd_cause_antigbm    == 1 |
      ckd_cause_fsgs       == 1 |
      ckd_cause_gn_nos     == 1 |
      ckd_cause_mcd        == 1 |
      ckd_cause_membranous == 1 |
      ckd_cause_mpgn       == 1 |
      ckd_cause_sle        == 1
    ),

    # any_nonstandard_ckd_cause: any of the 13 causes (most restrictive exclusion criterion)
    any_nonstandard_ckd_cause = as.integer(
      ckd_cause_aav        == 1 |
      ckd_cause_adpkd      == 1 |
      ckd_cause_alport     == 1 |
      ckd_cause_antigbm    == 1 |
      ckd_cause_fabry      == 1 |
      ckd_cause_fsgs       == 1 |
      ckd_cause_gn_nos     == 1 |
      ckd_cause_igan       == 1 |
      ckd_cause_mcd        == 1 |
      ckd_cause_membranous == 1 |
      ckd_cause_mpgn       == 1 |
      ckd_cause_other_pkd  == 1 |
      ckd_cause_sle        == 1
    ),

    # any_nonstandard_ckd_cause_excl_igan: any non-standard cause excluding IgAN
    # Use this flag when IgAN patients should be included in an analysis
    any_nonstandard_ckd_cause_excl_igan = as.integer(
      ckd_cause_aav        == 1 |
      ckd_cause_adpkd      == 1 |
      ckd_cause_alport     == 1 |
      ckd_cause_antigbm    == 1 |
      ckd_cause_fabry      == 1 |
      ckd_cause_fsgs       == 1 |
      ckd_cause_gn_nos     == 1 |
      ckd_cause_mcd        == 1 |
      ckd_cause_membranous == 1 |
      ckd_cause_mpgn       == 1 |
      ckd_cause_other_pkd  == 1 |
      ckd_cause_sle        == 1
    )

  ) %>%
  analysis$cached("ckd_causes", unique_indexes = "patid")

ckd_causes_table %>% count()
