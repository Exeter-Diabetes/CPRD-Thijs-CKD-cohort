
# Calculates electronic frailty index (https://pubmed.ncbi.nlm.nih.gov/26944937/) at index date
# Polypharmacy assumed to be 1

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

# Create a vector of all codelist names that start with "efi_"
efi_deficits <- grep("^efi_", names(codes_2024), value = TRUE)

# List of deficit names to shorten
short_deficit_map <- list(
  "efi_mobility_and_transfer_problems" = "efi_mobility_transfer"
)


# Function to shorten table name to avoid 64-character limit;
# otherwise use the original name
get_short_deficit <- function(deficit) {
  if (deficit %in% names(short_deficit_map)) {
    return(short_deficit_map[[deficit]])
  } else {
    return(deficit)
  }
}


# Pull out all raw code instances and cache with 'all_patid' prefix

analysis <- cprd$analysis("all_patid")

for (deficit in efi_deficits) {
  
  # If the codelist is not empty
  if (length(codes_2024[[deficit]]) > 0) {
    
    print(paste("making", deficit, "medcode table"))
    
    # Shorten deficit name
    deficit <- get_short_deficit(deficit)
    
    # Name intermediate table (e.g., "raw_efi_anaemia_haematinic_deficiency_medcodes")
    raw_tablename <- paste0("raw_", deficit, "_medcodes")
    
    # Get all relevant observation rows for this deficit
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



# get dates at 6 month intervals
dates <- seq(from = as.Date("2013-01-01"),
             to   = as.Date("2013-01-01"),
             by   = "6 months")

date_strings <- format(dates, "%Y-%m-%d")


for (d in date_strings) {
  
  analysis <- cprd$analysis(paste0(analysis_prefix, "_", d))
  
  
  index_date <- as.Date(d)
  print(d)
  
  
  # Clean deficit data and combine with index dates
  for (deficit in efi_deficits) {
    
    print(paste("merging index date with", deficit, "code occurrences"))
    
    # Shorten deficit name
    deficit <- get_short_deficit(deficit)
    
    # Define table names
    medcode_tablename <- paste0("raw_", deficit, "_medcodes")
    index_date_m_tablename <- paste0("full_", deficit, "_m") # use _m instead of _merge as column name too lon for MySQL otherwise
    
    # If the table exists
    if (exists(medcode_tablename)) {
      
      # Select relevant columns from the raw medcode table
      medcodes <- get(medcode_tablename) %>%
        select(patid, date = obsdate, code = medcodeid)
    }
    
    # Clean data by joining valid date info
    medcodes_clean <- medcodes %>%
      inner_join(cprd$tables$validDateLookup, by = "patid") %>%
      filter(date >= min_dob & date <= gp_end_date) %>%
      select(patid, date, code)
    
    rm(medcodes)
    
    # Merge with index info
    data <- medcodes_clean %>%
      mutate(datediff=datediff(date, index_date)) %>%
      analysis$cached(index_date_m_tablename, index="patid")
    
    rm(medcodes_clean)
    
    assign(index_date_m_tablename, data)
    
    rm(data)
  }
  
  ############################################################################################
  
  # Find if there has been any instance of pre_index_date occurrence of deficit
  
  # Initialise efi table
  efi <- cprd$tables$patient %>%
    select(patid)
  
  
  for (deficit in efi_deficits) {
    
    print(paste("Working out pre_index_date code occurrences for", deficit))
    
    # Shorten deficit name
    deficit <- get_short_deficit(deficit)
    
    # Define table and column names
    index_date_m_tablename <- paste0("full_", deficit, "_m")
    interim_efi_table <- paste0("efi_im_", deficit)
    pre_index_date_indicator <- paste0("pre_index_date_", deficit)
    pre_index_date_earliest_date_variable <- paste0("pre_index_date_earliest_", deficit)
    
    # Get earliest date of pre_index_date occurrence
    pre_index_date <- get(index_date_m_tablename) %>%
      filter(date <= index_date) %>%
      group_by(patid) %>%
      summarise(
        !!pre_index_date_earliest_date_variable := min(date, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      ungroup()
    
    # Convert the earliest date variable to a symbol
    pre_index_date_earliest_date_variable  <- sym(pre_index_date_earliest_date_variable)
    
    # Merge pre_index_date data into the eFI table and create a boolean indicator (1 if 
    # there was ever an occurrence of this deficit on or before index date)
    efi <- efi %>%
      left_join(pre_index_date, by = c("patid")) %>%
      mutate(
        !!pre_index_date_indicator := !is.na(!!pre_index_date_earliest_date_variable )
      ) %>%
      analysis$cached(
        interim_efi_table,
        indexes = c("patid")
      )
  }
  
  
  # Drop all 'pre_index_date_earliest_*' columns from the final 'efi' table for now
  efi <- efi %>%
    select(-matches("^pre_index_date_earliest_"))
  
  
  #   Fill all pre_index_date_diabetes occurrences with 1
  # - Fill all pre_index_date_polypharmacy with 1 (assume every patient is on 5 or more medications for now).
  
  efi <- efi %>%
    mutate(pre_index_date_efi_diabetes = 1L,
           pre_index_date_efi_polypharmacy = 1L) %>%
    analysis$cached(
      "efi_deficits",
      indexes = c("patid")
    )
  
  
  # Sum deficits and calculate score
  
  efi_deficits_short <- lapply(efi_deficits, get_short_deficit)
  deficit_vars <- paste0("pre_index_date_", efi_deficits_short)
  
  sql_adding_expr <- paste(deficit_vars, collapse = " + ")
  
  efi <- efi %>%
    mutate(efi_n_deficits=dbplyr::sql(sql_adding_expr),
           pre_index_date_efi_score = efi_n_deficits / 36,
           pre_index_date_efi_cat = case_when(
             pre_index_date_efi_score < 0.12 ~ "fit",
             pre_index_date_efi_score >= 0.12 & pre_index_date_efi_score < 0.24 ~ "mild",
             pre_index_date_efi_score >= 0.24 & pre_index_date_efi_score < 0.36 ~ "moderate",
             pre_index_date_efi_score >=0.36 ~ "severe" )) %>%
    analysis$cached(
      "efi",
      indexes = c("patid")
    )
}
