
# Produces all_patid_townsend_deprivation_score table

# See: https://github.com/drkgyoung/Exeter_Diabetes_codelists/blob/main/readme.md#townsend-deprivation-scores


################################################################################################################################

##################### SETUP ####################################################################################################

library(aurum)
library(tidyverse)
library(readxl)

cprd = CPRDData$new(cprdEnv = "nondiabetes-jun2024",cprdConf = "C:/Users/tj358/OneDrive - University of Exeter/CPRD/aurum.yaml")

analysis = cprd$analysis("all_patid")


################################################################################################################################

##################### IMPORT IMD/TDS/LSOA LOOKUPS ##############################################################################

setwd("C:/Users/tj358/OneDrive - University of Exeter/CPRD/2024/Scripts/CPRD-Thijs-CKD-cohort/townsend/")

imd_lsoa <- read_excel("File_1_-_IMD2019_Index_of_Multiple_Deprivation.xlsx", sheet="IMD2019") %>%
  select(lsoa_2011='LSOA code (2011)',
         imd_decile='Index of Multiple Deprivation (IMD) Decile')

townsend_lsoa <- read_csv("Scores- 2011 UK LSOA.csv") %>%
  select(lsoa_2011=GEO_CODE, tds_2011=TDS)


################################################################################################################################

##################### FIND MEDIAN TDS BY IMD DECILE ############################################################################

imd_townsend <- imd_lsoa %>%
  inner_join(townsend_lsoa, by="lsoa_2011") %>%
  select(-lsoa_2011) %>%
  group_by(imd_decile) %>%
  summarise(tds_2011=round(median(tds_2011), 3)) %>%
  ungroup()


################################################################################################################################

##################### MAKE MYSQL TABLE OF TOWNSEND SCORES ######################################################################


townsend_score <- cprd$tables$patientImd %>%
  select(patid, imd_decile) %>%
  inner_join(imd_townsend, by="imd_decile", copy=TRUE) %>%
  analysis$cached("townsend_score", unique_index="patid")
