####################################################################################
# Purpose: Read in Archbold demographic data from eleven species, clean and process data to have consistent format for vital rate analyses
# Author: Joshua Fowler
# Date: Aug 16, 2023
####################################################################################

##### Set up #####

# library(renv) # track package versions
# # renv::init()
# # renv::snapshot()
# renv::restore()
 

library(tidyverse)
library(readxl)
library(lubridate)


quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}
####################################################################################
###### Loading in each data for each species #######################################
####################################################################################

## Filepath for stored demographic data
filepath <- c("/Users/joshuacfowler/Dropbox/UofMiami/Demographic Data")

# Reading in data for Eryngium cuneifolium (ERYCUN)
# for ERYCUN, there is demographic data from rosemary balds at Archbold as well as at Royce Ranch (one dataset from a firelane transect and another from .
ERYCUN_abs_raw <- read_excel(path = paste0(filepath,"/UM demographic models/Ec demog MASTERSHEET 2022.xlsx"), sheet = "data-ecabs", col_types = "text") %>% 
  mutate(Site_tag = "archbold")

ERYCUN_apfi_raw <- read_excel(path = paste0(filepath,"/UM demographic models/Ec demog MASTERSHEET 2022.xlsx"), sheet = "data-ecapfi", col_types = "text") %>% 
  mutate(Site_tag = "royce_ranch")

ERYCUN_apsc_raw <- read_excel(path = paste0(filepath,"/UM demographic models/Ec demog MASTERSHEET 2022.xlsx"), sheet = "data-ecapsc", col_types = "text") %>% 
  mutate(Site_tag = "royce_ranch")


####################################################################################
###### Cleaning and merging together ERYCUN ########################################
####################################################################################
# Eryngium cuneifolium data was collected starting in 1994, and is some of the longest-term monitoring at Archbold, initiated by Eric Menges and Pedro Quintana-Ascencio

#### Loading in the data from the main Archbold sites
# removing columns related to XY coordinates only applicable in some plots as well as outdated info on burn history in the sites. We will merge in the most updated burn history later
# removing redundant demographic information columns (several contain the change in flowering stem number or rosette diameter between years)
# removing columns containing ground cover information about other species in the plots, cool information, but seems to mostly have been collected before 1994
# removing columns containing record of sand accretion on the metal id tag %>% This might be worth including as a covariate, but I don't know if it really tells us much and some of it is after the plant is already dead.
# Making a unique id for each plant in each patch in each bald
# Then pivoting the data to long format so that the demographic measures are in combined columns for each row
# We will also remove the assigned stage information from the dataset since we will work with the raw size measurements
# found one plant that had a typo in the spreadsheet, where the 2022 census was missing survival info and columns for size and reproduction where shifted one column to the left. so recoding this
# And recode the survival column. Archbold uses an idiosyncratic code to represent whether plants are seedling, alive, dead, previously dead. (0 = dead; 1 = alive; 2 = not found; 3 = new adult; 4 = not yet born; 5 = seedling; 6 = loose tag/pulled; 7 = putative seedling; 9 = previously dead, flag pulled)
# There are a few instances where tags were not found during a census but relocated in the following year or the second year after. Recoding these to be alive during those censuses, but we will obviously be missing size information during those years.
# There is one plant which has outlier size. I think this is likely just a decimal point typo based on size of the plant in preceding years, so moving that data point
# most of the time when there are no reproductive structures, the cell are left blank. particularly for new seedlings,  with rosette measurements, I think this is safe to assume is 0 stems. There are often times cases where this data is actually missing however

ERYCUN_abs <- ERYCUN_abs_raw %>%
  dplyr::select(-gps18,-ltreb,-pull,-X,-Y,-x.cor,-y.cor,-xy.cor.notes, 
                -contains("byr2"), -tsf, -breg, -sdno95, -stat, pull98, -cohort, -oldX, -oldtag, -contains("burn2"), -hobo0710, -rx0710, -starts_with("rx"), -tsf2018) %>% 
  dplyr::select(-contains("cs9"), -contains("cr9"),-contains("cr0"), 
                -contains("agr0"), -contains("agr9"), -contains("annsur9"), -contains("annsur0"), -contains("annsur1"),
                -contains("hstg9"), -contains("age9")) %>% 
  dplyr::select(-ag, -ca, -cev, -cp, -cs, -ec, -hc, -ld, -lc, -pc, -pr, -pb, -sab, -sel, -qi, -qu, -licania, -ceratiol, -groundco, -perlitte, -sppshrub, -distshru, -oakht02, -oakdis02, -quad, -otherspp) %>% 
  dplyr::select(-starts_with("sa"), -contains("pull"), -master) %>% 
  dplyr::select(-starts_with("stg")) %>% 
  mutate(row_id = row_number()) %>% 
  mutate(plant_id = paste(bald,patch,plant,TP, row_id, sep = "_")) %>% 
  mutate(across(everything(), as.character)) %>% 
  rename(first_year = yr1) %>% 
  pivot_longer(cols = !c(plant_id, Site_tag, bald, patch, plant, TP, row_id, first_year), names_to = c("measurement", "census_year"), names_sep = "(?<=[A-Za-z])(?=[0-9])") %>% 
  mutate(measurement = case_match(measurement, 
                                 "s" ~ "ARCHBOLD_surv",
                                 "stg" ~ "assigned_stage",
                                 "rdm" ~ "ros_diameter",
                                 "stm" ~ "flw_stem", "fl" ~ "flw_stem",
                                 "ht" ~ "max_stem_height",
                                 "hstm" ~ "herb_count",
                                 "h" ~ "flw_head",
                                 "sa" ~ "sand_accretion",
                                 "comm" ~ "comment",
                                 .default = as.character(measurement)),
         census_year = case_match(as.numeric(census_year),
                                  22 ~ 2022,
                                  21 ~ 2021,
                                  20 ~ 2020,
                                  19 ~ 2019,
                                  18 ~ 2018,
                                  17 ~ 2017,
                                  16 ~ 2016,
                                  15 ~ 2015,
                                  14 ~ 2014,
                                  13 ~ 2013,
                                  12 ~ 2012, 2012 ~ 2012,
                                  11 ~ 2011, 2011 ~ 2011,
                                  10 ~ 2010, 2010 ~ 2010,
                                  09 ~ 2009,
                                  08 ~ 2008,
                                  07 ~ 2007,
                                  06 ~ 2006,
                                  05 ~ 2005,
                                  04 ~ 2004,
                                  03 ~ 2003,
                                  02 ~ 2002,
                                  01 ~ 2001,
                                  00 ~ 2000,
                                  99 ~ 1999,
                                  98 ~ 1998,
                                  97 ~ 1997,
                                  96 ~ 1996,
                                  95 ~ 1995,
                                  94 ~ 1994,
                                  93 ~ 1993,
                                  92 ~ 1992,
                                  91 ~ 1991,
                                  90 ~ 1990,
                                  89 ~ 1989,
                                  88 ~ 1988,
                                  .default = as.numeric(census_year))) %>% 
  pivot_wider(id_cols = c(plant_id, Site_tag, bald, patch, plant, TP, row_id, first_year, census_year), names_from = measurement, values_from = value) %>% 
  mutate(ARCHBOLD_surv = case_when(plant_id == "57_2_1107_BA_6503" & census_year == 2022 ~ "1", TRUE ~ ARCHBOLD_surv),
         ros_diameter = case_when(plant_id == "57_2_1107_BA_6503" & census_year == 2022  ~ "4.7", TRUE ~ ros_diameter),
         flw_stem = case_when(plant_id == "57_2_1107_BA_6503" & census_year == 2022  ~ "1", TRUE ~ flw_stem),
         flw_head = case_when(plant_id == "57_2_1107_BA_6503" & census_year == 2022  ~ "7", TRUE ~ flw_head)) %>% 
  mutate(across(c(plant_id, Site_tag, bald, patch, plant, TP, row_id, comment), as.character)) %>% 
  mutate(across(c(first_year, census_year, ARCHBOLD_surv, flw_stem, flw_head, herb_count), function(x) suppressWarnings(as.integer(x)))) %>% 
  mutate(across(c(ros_diameter, max_stem_height), function(x) suppressWarnings(as.numeric(x)))) %>% 
  mutate(ARCHBOLD_surv = case_when(ARCHBOLD_surv == 20 ~ 2, TRUE ~ ARCHBOLD_surv),
         surv = case_match(as.numeric(ARCHBOLD_surv),
                           0 ~ 0,
                           1 ~ 1,
                           2 ~ 0, # 2 means tag not found, but need to adjust because they recorded 2 sometimes for multiple years
                           3 ~ 1, # 3 is new adult
                           4 ~ NA,
                           5 ~ 1, # 5 is new seedling
                           6 ~ 0, # 6 stands for tag not found, so this is assigning those plants as dead. I will go through and check that plants don't reappear later as alive.
                           7 ~ NA, # There's only one of these, which is "putative seedling", so I think maybe plant idea was uncertain on small individual. but it doesn't have any other measurements anyways.
                           8 ~ 0,
                           9 ~ NA,
                           12 ~ NA, # This is for 1 site, which had change in layout, and so plants were outside of new census, True NA's not necessarily dead))
                           .default = as.numeric(ARCHBOLD_surv))) %>% 
  group_by(plant_id) %>%
  mutate(surv = case_when(lead(ARCHBOLD_surv) == 9 & ARCHBOLD_surv == 2 ~ NA,
                          lead(ARCHBOLD_surv) == 2 & ARCHBOLD_surv == 2 ~ NA,
                          dplyr::lead(ARCHBOLD_surv) == 2 & ARCHBOLD_surv == 0 ~ NA,
                          dplyr::lead(ARCHBOLD_surv) == 0 & ARCHBOLD_surv == 2 ~ NA,
                          dplyr::lead(ARCHBOLD_surv) == 0 & ARCHBOLD_surv == 0 ~ NA,
                          dplyr::lead(ARCHBOLD_surv) == 6 & ARCHBOLD_surv == 6 ~ NA,
                          dplyr::lead(ARCHBOLD_surv) == 0 & ARCHBOLD_surv == 6 ~ NA,
                          dplyr::lead(ARCHBOLD_surv) == 2 & ARCHBOLD_surv == 6 ~ NA,
                          
                          dplyr::lag(ARCHBOLD_surv) == 1 & ARCHBOLD_surv == 2 ~ 1,
                          dplyr::lag(ARCHBOLD_surv) == 1 & ARCHBOLD_surv == 0 ~ 1,
                          dplyr::lag(ARCHBOLD_surv, n = 2) == 1 & dplyr::lag(ARCHBOLD_surv) == 2 & ARCHBOLD_surv == 2 ~ 1,
                          TRUE ~ surv)) %>% 
  mutate(birth_year = case_when(ARCHBOLD_surv == 5 ~ census_year,
                                ARCHBOLD_surv == 3 & ros_diameter<3 & flw_stem == 0 ~ census_year,
                                TRUE ~ NA)) %>% 
  fill(birth_year, .direction = "updown") %>% 
  mutate(ros_diameter = case_when(plant_id == "54_9_155_NA_4531" & census_year == 2018 ~ ros_diameter*.1, TRUE ~ ros_diameter)) %>% 
  mutate(flw_stem = case_when(first_year == census_year & !is.na(ros_diameter) & is.na(flw_stem) ~ 0, TRUE ~ flw_stem)) %>% 
  mutate(flw_stem = case_when(!is.na(ros_diameter) & is.na(flw_stem) & is.na(flw_head) ~ 0,
                              !is.na(ros_diameter) & is.na(max_stem_height) & is.na(flw_head) ~ 0, TRUE ~ flw_stem)) %>% 
  dplyr::select(plant_id,Site_tag,bald, patch, plant, TP, row_id, first_year, birth_year, census_year, ARCHBOLD_surv, surv, ros_diameter, max_stem_height, flw_stem, flw_head, herb_count, comment)


#### Loading in the Royce Ranch fire lane data 
# Removing some extraneous location info
# removing redundant demographic information columns (columns with change in size and in number of stems)
# removing measurement of sand accretion
# creating unique id for each plant
# Then pivoting to combine the columns for each vital rate measurement
# And recode the survival column. Archbold uses an idiosyncratic code to represent whether plants are seedling, alive, dead, previously dead. (0 = dead; 1 = alive; 2 = not found; 3 = new adult; 4 = not yet born; 5 = seedling; 6 = loose tag/pulled; 7 = putative seedling; 9 = previously dead, flag pulled)
# For this dataset, there are several instances, particulary recruits in the 90s, where new seedlings have a size when they appear but NA survival, so I am correcting this
# And correcting cases where tag was not found  in multiple years. Assuming plant is dead, unless it is re-found in the following year or the next

ERYCUN_apfi <- ERYCUN_apfi_raw %>%
  dplyr::select(-pull) %>% 
  dplyr::select(-contains("rgr"), -contains("chst")) %>% 
  dplyr::select(-starts_with("sa")) %>% 
  mutate(row_id = row_number(), patch = "Firelane") %>% 
  mutate(plant_id = paste(patch, quad, plant,TP, row_id, sep = "_")) %>% 
  mutate(across(everything(), as.character)) %>% 
  rename(first_year = yr1) %>% 
  pivot_longer(cols = !c(plant_id, Site_tag, patch, quad, plant, TP, row_id, first_year), names_to = c("measurement", "census_year"), names_sep = "(?<=[A-Za-z])(?=[0-9])") %>% 
  mutate(measurement = case_match(measurement, 
                                "s" ~ "ARCHBOLD_surv", "a" ~ "ARCHBOLD_surv",
                                "stg" ~ "assigned_stage",
                                "rdm" ~ "ros_diameter", "rd" ~ "ros_diameter", "bdm" ~ "ros_diameter",
                                "stm" ~ "flw_stem", "fl" ~ "flw_stem", "st" ~ "flw_stem",
                                "ht" ~ "max_stem_height",
                                "hstm" ~ "herb_count",
                                "h" ~ "flw_head", "hd" ~ "flw_head", "he" ~ "flw_head",
                                "sa" ~ "sand_accretion",
                                "comm" ~ "comment",
                                .default = as.character(measurement)),
       census_year = case_match(as.numeric(census_year),
                                22 ~ 2022,
                                21 ~ 2021,
                                20 ~ 2020,
                                19 ~ 2019,
                                18 ~ 2018,
                                17 ~ 2017,
                                16 ~ 2016,
                                15 ~ 2015,
                                14 ~ 2014,
                                13 ~ 2013, 2013 ~ 2013,
                                12 ~ 2012, 2012 ~ 2012,
                                11 ~ 2011, 2011 ~ 2011,
                                10 ~ 2010, 2010 ~ 2010,
                                09 ~ 2009,
                                08 ~ 2008,
                                07 ~ 2007,
                                06 ~ 2006,
                                05 ~ 2005,
                                04 ~ 2004,
                                03 ~ 2003,
                                02 ~ 2002,
                                01 ~ 2001,
                                00 ~ 2000,
                                99 ~ 1999,
                                98 ~ 1998,
                                97 ~ 1997,
                                96 ~ 1996,
                                95 ~ 1995,
                                94 ~ 1994,
                                93 ~ 1993,
                                92 ~ 1992,
                                91 ~ 1991,
                                90 ~ 1990,
                                89 ~ 1989,
                                88 ~ 1988,
                                .default = as.numeric(census_year))) %>% 
  pivot_wider(id_cols = c(plant_id, Site_tag, patch, quad, plant, TP, row_id, first_year, census_year), names_from = measurement, values_from = value) %>% 
  mutate(across(c(plant_id, Site_tag, patch, plant, TP, row_id, comment), as.character)) %>% 
  mutate(across(c(first_year, census_year, ARCHBOLD_surv, flw_stem, flw_head, herb_count), as.integer)) %>% 
  mutate(across(c(ros_diameter, max_stem_height), as.numeric)) %>% 
  mutate(surv = case_match(ARCHBOLD_surv,
                           0 ~ 0,
                           1 ~ 1,
                           2 ~ 0, # 2 means tag not found, but need to adjust because they recorded 2 sometimes for multiple years
                           3 ~ 1, # 3 is new adult
                           4 ~ NA,
                           5 ~ 1, # 5 is new seedling
                           6 ~ 0, # 6 stands for tag not found, so this is assigning those plants as dead. I will go through and check that plants don't reappear later as alive.
                           7 ~ 1, # There's only one of these, which is "putative seedling", so I think maybe plant idea was uncertain on small individual. but it doesn't have any other measurements anyways.
                           8 ~ 0,
                           9 ~ NA,
                           12 ~ NA, # This is for 1 site, which had change in layout, and so plants were outside of new census))
                           .default = as.numeric(ARCHBOLD_surv))) %>% 
  mutate(surv = case_when(is.na(surv) & !is.na(ros_diameter) & first_year == census_year ~ 1,
                          TRUE ~ surv)) %>% 
  mutate(surv = case_when(lead(ARCHBOLD_surv) == 9 & ARCHBOLD_surv == 2 ~ NA,
                          lead(ARCHBOLD_surv) == 2 & ARCHBOLD_surv == 2 ~ NA,
                          lead(ARCHBOLD_surv) == 9 & ARCHBOLD_surv == 8 ~ NA,
                          lead(ARCHBOLD_surv) == 8 & ARCHBOLD_surv == 8 ~ NA,
                          lead(ARCHBOLD_surv) == 8 & ARCHBOLD_surv == 0 ~ NA,
                          lead(ARCHBOLD_surv) == 0 & ARCHBOLD_surv == 8 ~ NA,
                          lead(ARCHBOLD_surv) == 0 & ARCHBOLD_surv == 0 ~ NA,
                          lead(ARCHBOLD_surv) == 1 & ARCHBOLD_surv == 9 ~ 0,
                          TRUE ~ surv)) %>% 
  mutate(birth_year = case_when(ARCHBOLD_surv == 5 ~ census_year,
                                ARCHBOLD_surv == 3 & ros_diameter<3 & flw_stem == 0 ~ census_year,
                                TRUE ~ NA)) %>% 
  fill(birth_year, .direction = "updown") %>% 
  dplyr::select(plant_id,Site_tag, patch, quad, plant, TP, row_id, first_year,  birth_year, census_year, ARCHBOLD_surv, surv, ros_diameter, max_stem_height, flw_stem, flw_head, herb_count, comment)

  
#### Loading in the Royce Ranch Scrub data
# removing some extraneous meta columns
# removing redundant demographic information columns (several contain the change in flowering stem number or rosette diameter between years)
# removing measurement of sand accretion
# removing comments and record of burn/mowing
# Removing the assigned stage information from the dataset since we will work with the raw size measurements
# creating a unique id for each plant
# removing the column for the first census year for this species because in this dataset, the column not complete and I'll just calculate it later
# creating unique id for each plant
# Then pivoting to combine the columns for each vital rate measurement
# And recode the survival column. Archbold uses an idiosyncratic code to represent whether plants are seedling, alive, dead, previously dead. (0 = dead; 1 = alive; 2 = not found; 3 = new adult; 4 = not yet born; 5 = seedling; 6 = loose tag/pulled; 7 = putative seedling; 9 = previously dead, flag pulled)
# recoding a two size measurements, which are marked as 999.0 to be NA.
# Then there are a few cases where plants where missing but found alive in later years. Correcting survival to reflect this
# And calculating the first observation year for each individual for this dataset
ERYCUN_apsc <- ERYCUN_apsc_raw %>%
  dplyr::select(-microhabitat, -pull, -pull_temp, -pull98, -oldtag, -bigq, -s18, -subquad) %>% 
  dplyr::select(-contains("annsur"), -contains("rgr"), -contains("chst")) %>% 
  dplyr::select(-starts_with("sa")) %>% 
  dplyr::select(-starts_with("comm"), -contains("mow"), -contains("burn")) %>% 
  dplyr::select(-starts_with("stg")) %>% 
  mutate(row_id = row_number(), patch = "Scrub") %>% 
  mutate(plant_id = paste(patch, quad, plant,TP, row_id, sep = "_")) %>% 
  mutate(across(everything(), as.character)) %>% 
  pivot_longer(cols = !c(plant_id, Site_tag, patch, quad, plant, TP, row_id), names_to = c("measurement", "census_year"), names_sep = "(?<=[A-Za-z])(?=[0-9])") %>% 
  mutate(measurement = case_match(measurement, 
                                  "s" ~ "ARCHBOLD_surv", "a" ~ "ARCHBOLD_surv",
                                  "stg" ~ "assigned_stage",
                                  "rdm" ~ "ros_diameter", "rd" ~ "ros_diameter", "bdm" ~ "ros_diameter",
                                  "stm" ~ "flw_stem", "fl" ~ "flw_stem", "st" ~ "flw_stem",
                                  "ht" ~ "max_stem_height",
                                  "hstm" ~ "herb_count",
                                  "h" ~ "flw_head", "hd" ~ "flw_head", "he" ~ "flw_head", "hea" ~ "flw_head",
                                  "sa" ~ "sand_accretion",
                                  "comm" ~ "comment",
                                  .default = as.character(measurement)),
         census_year = case_match(as.numeric(census_year),
                                  22 ~ 2022,
                                  21 ~ 2021,
                                  20 ~ 2020,
                                  19 ~ 2019,
                                  18 ~ 2018,
                                  17 ~ 2017,
                                  16 ~ 2016,
                                  15 ~ 2015,
                                  14 ~ 2014,
                                  13 ~ 2013, 2013 ~ 2013,
                                  12 ~ 2012, 2012 ~ 2012,
                                  11 ~ 2011, 2011 ~ 2011,
                                  10 ~ 2010, 2010 ~ 2010,
                                  09 ~ 2009,
                                  08 ~ 2008,
                                  07 ~ 2007,
                                  06 ~ 2006,
                                  05 ~ 2005,
                                  04 ~ 2004,
                                  03 ~ 2003,
                                  02 ~ 2002,
                                  01 ~ 2001,
                                  00 ~ 2000,
                                  99 ~ 1999,
                                  98 ~ 1998,
                                  97 ~ 1997,
                                  96 ~ 1996,
                                  95 ~ 1995,
                                  94 ~ 1994,
                                  93 ~ 1993,
                                  92 ~ 1992,
                                  91 ~ 1991,
                                  90 ~ 1990,
                                  89 ~ 1989,
                                  88 ~ 1988,
                                  .default = as.numeric(census_year))) %>% 
  pivot_wider(id_cols = c(plant_id, Site_tag, patch, quad, plant, TP, row_id, census_year), names_from = measurement, values_from = value) %>% 
  mutate(across(c(plant_id, Site_tag, patch, plant, TP, row_id), as.character)) %>% 
  mutate(across(c( census_year, ARCHBOLD_surv, flw_stem, flw_head, herb_count), function(x) suppressWarnings(as.integer(x)))) %>%  
  mutate(across(c(ros_diameter, max_stem_height), function(x) suppressWarnings(as.numeric(x)))) %>% 
  mutate(surv = case_match(ARCHBOLD_surv,
                           0 ~ 0,
                           1 ~ 1,
                           2 ~ 0, # 2 means tag not found, but need to adjust because they recorded 2 sometimes for multiple years
                           3 ~ 1, # 3 is new adult
                           4 ~ NA,
                           5 ~ 1, # 5 is new seedling
                           6 ~ 0, # 6 stands for tag not found, so this is assigning those plants as dead. I will go through and check that plants don't reappear later as alive.
                           7 ~ 1, # There's only one of these, which is "putative seedling", so I think maybe plant idea was uncertain on small individual. but it doesn't have any other measurements anyways.
                           8 ~ 0, 88 ~ 0,
                           9 ~ NA,
                           12 ~ NA, # This is for 1 site, which had change in layout, and so plants were outside of new census))
                           .default = as.numeric(ARCHBOLD_surv))) %>% 
  mutate(ros_diameter = case_when(ros_diameter == 999.0 ~ NA, TRUE ~ ros_diameter)) %>% 
  mutate(surv = case_when(lead(ARCHBOLD_surv) == 9 & ARCHBOLD_surv == 2 ~ NA,
                          lead(ARCHBOLD_surv) == 2 & ARCHBOLD_surv == 2 ~ NA,
                          lead(ARCHBOLD_surv) == 9 & ARCHBOLD_surv == 8 ~ NA,
                          lead(ARCHBOLD_surv) == 8 & ARCHBOLD_surv == 8 ~ NA,
                          lead(ARCHBOLD_surv) == 8 & ARCHBOLD_surv == 2~ NA,
                          lead(ARCHBOLD_surv) == 8 & ARCHBOLD_surv == 0 ~ NA,
                          lead(ARCHBOLD_surv) == 0 & ARCHBOLD_surv == 8 ~ NA,
                          lead(ARCHBOLD_surv) == 0 & ARCHBOLD_surv == 0 ~ NA,
                          lead(ARCHBOLD_surv) == 1 & ARCHBOLD_surv == 9 ~ 0,
                          lag(ARCHBOLD_surv) == 1 & ARCHBOLD_surv == 8 ~ 1,
                          TRUE ~ surv))  %>% 
  group_by(plant_id) %>% 
  mutate(census_temp = case_when(!is.na(surv) ~ census_year),first_year = min(census_temp, na.rm = T)) %>% 
  mutate(birth_year = case_when(ARCHBOLD_surv == 5 ~ census_year,
                                ARCHBOLD_surv == 3 & ros_diameter<3 & flw_stem == 0 ~ census_year,
                                TRUE ~ NA)) %>% 
  fill(birth_year, .direction = "updown") %>% 
  ungroup()  %>% 
  dplyr::select(plant_id,Site_tag, patch,  quad, plant, TP, row_id, first_year, birth_year, census_year, ARCHBOLD_surv, surv, ros_diameter, max_stem_height, flw_stem, flw_head, herb_count)


# Still need to decide on the best way/if to combine the spatial hierarchy across sites. 

# combining data from the three sites
# Then creating columns for lagged census year size and reproduction 

ERYCUN <- ERYCUN_abs %>% 
  full_join(ERYCUN_apfi) %>% 
  full_join(ERYCUN_apsc) %>% 
  # filter(!is.na(surv)) %>% 
  mutate(flw_status = case_when(flw_stem>0 ~ 1, flw_stem==0 ~ 0,
                                surv == 0 ~ NA)) %>% 
  group_by(plant_id) %>% 
  arrange(plant_id, census_year) %>% 
  mutate(year.t1 = census_year,
         surv.t1 = surv,
         ros_diameter.t1 = ros_diameter,
         max_stem_height.t1 = max_stem_height,
         flw_status.t1 = flw_status,
         flw_stem.t1 = flw_stem,
         flw_head.t1 = flw_head,
         herb_count.t1 = herb_count,
         age.t1 = census_year-birth_year) %>%  
  mutate(year.t = dplyr::lag(year.t1, n = 1, default = NA),
         ros_diameter.t = dplyr::lag(ros_diameter.t1, n = 1, default = NA),
         max_stem_height.t = dplyr::lag(max_stem_height.t1, n = 1, default = NA),
         flw_status.t = dplyr::lag(flw_status.t1, n = 1, default = NA),
         flw_stem.t = dplyr::lag(flw_stem.t1, n = 1, default = NA),
         flw_head.t = dplyr::lag(flw_head.t1, n = 1, default = NA),
         herb_count.t = dplyr::lag(herb_count.t1, n = 1, default = NA),
         age.t = dplyr::lag(age.t1, n = 1, default = NA)) %>% 
  filter(!is.na(surv.t1)) %>% 
  select(plant_id,Site_tag, bald, patch, quad, plant, TP, row_id, first_year, birth_year,
         year.t1, ARCHBOLD_surv, surv.t1, ros_diameter.t1, max_stem_height.t1, flw_status.t1, flw_stem.t1, flw_head.t1, herb_count.t1, age.t1,
         year.t, ros_diameter.t, max_stem_height.t, flw_status.t, flw_stem.t, flw_head.t, herb_count.t, age.t)



####################################################################################
###### Processing the soil nutrient analysis data ##################################
####################################################################################
nutrients_raw <- read_csv(file = "/Users/joshuacfowler/Dropbox/UofMiami/Citrus Sampling Logistics/IFAS soil nutrient data/R8090.csv", skip = 14, skip_empty_rows = TRUE) %>% 
    rename(ID = `ID#`, Lab_Number = `Lab Number`)
  
nutrients_id_key <- read_xlsx(path = "/Users/joshuacfowler/Dropbox/UofMiami/Citrus Sampling Logistics/soil_analysis_id.xlsx", sheet = "soil nutrient analysis id") %>% 
  mutate(bald = case_when(project == "2024 rosemary bald sampling" ~ str_remove(Bald_id, ".+ ")),
         maehr_site = case_when(project == "citrus sampling" ~ str_remove(Bald_id, ".+ ")))


nutrients <- nutrients_raw %>% 
  filter(if_any(everything(), ~ !is.na(.))) %>% 
  filter(ID %in% 1:154) %>% 
  select(-...14, -...15) %>% 
  mutate(across(ID:pH, ~as.numeric(.))) %>% 
  left_join(nutrients_id_key, by = c("ID" = "Soil_analysis_number"))

# ggplot(data = nutrients)+
#   geom_histogram(aes(x = pH, fill = project))
# 
# ggplot(data = nutrients)+
#   geom_histogram(aes(x = OrgMat, fill = project))
# 
# ggplot(data = nutrients)+
#   geom_histogram(aes(x = Fe, fill = project))

####################################################################################
###### Merging the datasets with the environmental covariates ######################
####################################################################################
ARCHBOLD_fire_to2022 <- read_xlsx(path = "~/Dropbox/UofMiami/Balds2009_FireIntensityArea_Through2022.xlsx", sheet = "Balds2009_FireIntensityArea", guess_max = 1048576)# guess_max makes the function look deeper in the columns to assign type

fire_summary <- ARCHBOLD_fire_to2022 %>% 
  group_by(Bald_U, Bald_) %>% 
  filter(INTENSITY != 0) %>% 
  summarize(last_fire = max(as.numeric(Year), na.rm = TRUE),
            time_since_fire = 2023-max(as.numeric(Year), na.rm = TRUE),
            fire_list = toString(unique(as.numeric(Year))),
            fire_frequency = length(unique(as.numeric(Year)))) %>% 
  ungroup() %>% 
  add_row(Bald_U = "1S", Bald_ = 1, last_fire = NA, time_since_fire = NA, fire_list = NA, fire_frequency = 0) %>% 
  mutate(bald = case_when(as.character(Bald_) %in% c("1", "45","95") ~ as.character(Bald_),
                          TRUE ~ as.character(Bald_U)),
         bald = case_when(Bald_U == "1S" ~ "1S", 
                          Bald_U == "95W" ~ "95W", Bald_U == "95N" ~ "95N", TRUE ~ bald),
         bald_simple = as.character(parse_number(bald)))
  

# Now getting the relative elevation from the old fire history file

elev.df <- read_xlsx(path = "~/Dropbox/UofMiami/Experiment Set up/firehistory_thru2018.xlsx", sheet = "Rx_Freq", guess_max = 1048576) %>% # guess_max makes the function look deeper in the columns to assign type
  rename(rel_elev = rel.eve) %>% 
  dplyr::select(bald, rel_elev) %>% 
  mutate(
    bald = case_when(bald == "01S" ~ "1S",
                     bald == "01N" ~ "1",
                     bald == "02" ~ "2",
                     bald == "05E" ~ "5E",
                     bald == "07N" ~ "7N",
                     bald == "35N" ~ "35",
                     bald == "65E" ~ "65",
                     bald == "85N" ~ "85",
                     bald == "45N" ~ "45",
                     bald == "70N"~ "70",
                     bald == "72N" ~ "72",
                     bald == "95" ~ "95S",
                     TRUE ~ bald),
    bald_simple = as.character(parse_number(bald)))



# Pulling out fire history relative to the observation year for each datapoint
# needing to change a few of the names of balds to make them line up with names in elev and fire datasets. this also should be double checked because some you have to choose N/S
ERYCUN_covariates <- ERYCUN %>% 
  left_join(elev.df) %>% 
  left_join(fire_summary, by = join_by(bald)) %>% 
    mutate(rel_elev = case_when(bald == 95 ~ 0.68, TRUE ~ rel_elev)) %>% 
  filter(Site_tag == "archbold")






#
last_fire_fx <- function(x,y){ max(as.numeric(x)[as.numeric(x)<=y])}
fire_frequency_fx <- function(x,y){length(as.numeric(x)[as.numeric(x)<=y])/(y - 1950)}



#
ERYCUN_covariates$last_fire_actual <- unlist(Map(f = last_fire_fx, x = strsplit(unlist(ERYCUN_covariates$fire_list), ", "), y= ERYCUN_covariates$year.t1))
ERYCUN_covariates$time_since_fire_actual <- ERYCUN_covariates$year.t1 - ERYCUN_covariates$last_fire_actual
ERYCUN_covariates$fire_frequency_actual <- unlist(Map(f = fire_frequency_fx, x = strsplit(unlist(ERYCUN_covariates$fire_list), ", "), y= ERYCUN_covariates$year.t1))

ERYCUN_covariates <- ERYCUN_covariates %>% 
  mutate(fire_frequency_actual = case_when(bald == "1S" ~ 0/(year.t1-1950),
                                           TRUE ~ fire_frequency_actual)) %>% 
  dplyr::select(-last_fire, -fire_frequency, -time_since_fire)

write.csv(ERYCUN_covariates, "ERYCUN_covariates.csv")



#########

# write_csv(ERYCUN_covariates, paste0("Documents/R_projects/seed_predation_at_Archbold","/cleaned_data", "/ERYCUN_covariates.csv"))
# write_csv(HYPCUM_covariates, paste0(filepath,"/cleaned_data", "/HYPCUM_covariates.csv"))
# write_csv(POLBAS_covariates, paste0(filepath,"/cleaned_data", "/POLBAS_covariates.csv"))
# write_csv(BALANG_covariates, paste0(filepath,"/cleaned_data", "/BALANG_covariates.csv"))
# write_csv(CHAFAS_covariates, paste0(filepath,"/cleaned_data", "/CHAFAS_covariates.csv"))
# write_csv(PARCHA_covariates, paste0(filepath,"/cleaned_data", "/PARCHA_covariates.csv"))
# write_csv(CHAFLO, paste0(filepath,"/cleaned_data", "/CHAFLO.csv"))
# write_csv(LIAOHL, paste0(filepath,"/cleaned_data", "/LIAOHL.csv"))


