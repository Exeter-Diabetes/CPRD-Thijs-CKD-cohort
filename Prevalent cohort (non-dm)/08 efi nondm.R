# Calculates electronic frailty index (https://pubmed.ncbi.nlm.nih.gov/26944937/) at each index date
# Polypharmacy assumed to be 1 for all patients
# Adapted from 08 efi dm.R for the non-DM CKD cohort

############################################################################################

# Setup
library(tidyverse)
library(aurum)
library(EHRBiomarkr)
rm(list=ls())

cprd = CPRDData$new(cprdEnv = "nondiabetes-jun2024", cprdConf = "C:/Users/tj358/OneDrive - University of Exeter/CPRD/aurum.yaml")
codesets = cprd$codesets()
codes_2024 = codesets$getAllCodeSetVersion(v = "01/06/2024")

analysis_prefix = "ckd"


############################################################################################

# Create a vector of all codelist names that start with "efi_"
efi_deficits <- grep("^efi_", names(codes_2024), value = TRUE)

# List of deficit names to shorten (to avoid 64-character MySQL table name limit)
short_deficit_map <- list(
  "efi_mobility_and_transfer_problems" = "efi_mobility_transfer"
)

# Function to shorten table name where necessary; otherwise use the original name
get_short_deficit <- function(deficit) {
  if (deficit %in% names(short_deficit_map)) {
    return(short_deficit_map[[deficit]])
  } else {
    return(deficit)
  }
}


############################################################################################

# Pull out all raw code instances and cache with 'all_patid' prefix

analysis <- cprd$analysis("all_patid")

for (deficit in efi_deficits) {

  # If the codelist is not empty
  if (length(codes_2024[[deficit]]) > 0) {

    print(paste("making", deficit, "medcode table"))

    # Shorten deficit name for table naming
    deficit_short <- get_short_deficit(deficit)

    raw_tablename <- paste0("raw_", deficit_short, "_medcodes")

    data <- cprd$tables$observation %>%
      inner_join(codes_2024[[deficit]], by = "medcodeid") %>%
      analysis$cached(
        raw_tablename,
        indexes = c("patid", "obsdate")
      )

    assign(raw_tablename, data)
  }
}


############################################################################################

# Clean medcodes then merge with index dates
# Remove medcodes before DOB or after lcd/deregistration

# 6-monthly dates for 2019-2021 (prevalent cohort), then 3-monthly from 2021 onwards
# (3-monthly required for sequential trial emulation of SGLT2i in non-DM CKD)
dates <- unique(c(
  seq(from = as.Date("2019-03-01"), to = as.Date("2020-09-01"), by = "6 months"),
  seq(from = as.Date("2021-03-01"), to = as.Date("2024-03-01"), by = "3 months")
))

date_strings <- format(dates, "%Y-%m-%d")


for (d in date_strings) {

  analysis <- cprd$analysis(paste0(analysis_prefix, "_", d))

  index_date <- as.Date(d)
  print(d)


  # Clean deficit data and combine with index dates
  for (deficit in efi_deficits) {

    print(paste("merging index date with", deficit, "code occurrences"))

    deficit_short <- get_short_deficit(deficit)

    medcode_tablename     <- paste0("raw_", deficit_short, "_medcodes")
    index_date_m_tablename <- paste0("full_", deficit_short, "_m") # use _m instead of _merge as column name too long for MySQL otherwise

    if (exists(medcode_tablename)) {

      medcodes <- get(medcode_tablename) %>%
        select(patid, date = obsdate, code = medcodeid)

      medcodes_clean <- medcodes %>%
        inner_join(cprd$tables$validDateLookup, by = "patid") %>%
        filter(date >= min_dob & date <= gp_end_date) %>%
        select(patid, date, code)

      rm(medcodes)

      data <- medcodes_clean %>%
        mutate(datediff = datediff(date, index_date)) %>%
        analysis$cached(index_date_m_tablename, index = "patid")

      rm(medcodes_clean)

      assign(index_date_m_tablename, data)

      rm(data)
    }
  }


  ############################################################################################

  # Find if there has been any pre-index-date occurrence of each deficit

  # Initialise efi table from all patients
  efi <- cprd$tables$patient %>%
    select(patid)


  for (deficit in efi_deficits) {

    print(paste("Working out pre_index_date code occurrences for", deficit))

    deficit_short <- get_short_deficit(deficit)

    index_date_m_tablename              <- paste0("full_", deficit_short, "_m")
    interim_efi_table                   <- paste0("efi_im_", deficit_short)
    pre_index_date_indicator            <- paste0("pre_index_date_", deficit_short)
    pre_index_date_earliest_date_variable <- paste0("pre_index_date_earliest_", deficit_short)

    # Get earliest date of pre-index-date occurrence
    pre_index_date <- get(index_date_m_tablename) %>%
      filter(date <= index_date) %>%
      group_by(patid) %>%
      summarise(
        !!pre_index_date_earliest_date_variable := min(date, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      ungroup()

    pre_index_date_earliest_date_variable_sym <- sym(pre_index_date_earliest_date_variable)

    # Merge into efi table and create binary indicator (1 if ever recorded on/before index date)
    efi <- efi %>%
      left_join(pre_index_date, by = "patid") %>%
      mutate(
        !!pre_index_date_indicator := !is.na(!!pre_index_date_earliest_date_variable_sym)
      ) %>%
      analysis$cached(
        interim_efi_table,
        indexes = c("patid")
      )
  }


  # Drop all 'pre_index_date_earliest_*' columns for now
  efi <- efi %>%
    select(-matches("^pre_index_date_earliest_"))


  # Polypharmacy: assumed to be 1 for all patients (assume everyone is on >= 5 medications)
  # NB: unlike the DM cohort, efi_diabetes is NOT forced to 1 here - it is computed from
  #     actual codes above, and should be ~0 for all patients in this non-DM cohort
  efi <- efi %>%
    mutate(pre_index_date_efi_polypharmacy = 1L) %>%
    analysis$cached(
      "efi_deficits",
      indexes = c("patid")
    )


  ############################################################################################

  # Sum deficits and calculate EFI score (sum / 36 deficits)

  efi_deficits_short <- lapply(efi_deficits, get_short_deficit)
  deficit_vars       <- paste0("pre_index_date_", efi_deficits_short)

  sql_adding_expr <- paste(deficit_vars, collapse = " + ")

  efi <- efi %>%
    mutate(
      efi_n_deficits         = dbplyr::sql(sql_adding_expr),
      pre_index_date_efi_score = efi_n_deficits / 36,
      pre_index_date_efi_cat = case_when(
        pre_index_date_efi_score <  0.12                                    ~ "fit",
        pre_index_date_efi_score >= 0.12 & pre_index_date_efi_score < 0.24 ~ "mild",
        pre_index_date_efi_score >= 0.24 & pre_index_date_efi_score < 0.36 ~ "moderate",
        pre_index_date_efi_score >= 0.36                                    ~ "severe"
      )
    ) %>%
    analysis$cached(
      "efi",
      indexes = c("patid")
    )
}
