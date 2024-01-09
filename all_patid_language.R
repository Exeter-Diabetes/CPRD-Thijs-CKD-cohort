############################################################################################

# Setup

library(aurum)
library(tidyverse)

cprd = CPRDData$new(cprdEnv = "test-remote-full", cprdConf = "C:/Users/tj358/OneDrive - University of Exeter/CPRD/aurum.yaml")
codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "31/10/2021")

analysis = cprd$analysis("all")


############################################################################################

# Find language codes

raw_gp_language <- cprd$tables$observation %>% 
  inner_join(codes$language, by="medcodeid") %>%
  select(patid, obsdate, medcodeid) %>%
  analysis$cached("raw_gp_language", indexes=c("patid", "obsdate", "medcodeid"))


gp_language <- raw_gp_language %>% 
  inner_join(codes$language, by="medcodeid") %>%
  
  group_by(patid, language_cat) %>%
  summarise(lang_code_count=n(),                                                    # code count per person per language category
            latest_date_per_cat=max(obsdate, na.rm=TRUE)) %>%                       
  ungroup() %>%
  
  group_by(patid) %>%
  filter(lang_code_count==max(lang_code_count, na.rm=TRUE)) %>%                       # only keep categories with most counts
  filter(n()==1 | latest_date_per_cat==max(latest_date_per_cat, na.rm=TRUE)) %>%    # keep if 1 row per person, or if on latest date
  filter(n()==1) %>%                                                                # keep if 1 row per person
  ungroup() %>% 
  
  select(patid, gp_language_cat=language_cat) %>%
  analysis$cached("gp_language",unique_indexes="patid", indexes="gp_language_cat")

