### This script calculates recoveries for EMSL

## Created: 2023-12-19 by AMP for ionic strength experiment

# load packages
require(pacman)
pacman::p_load(tidyverse, # keep things tidy
               janitor, # useful for simplifying column names
               googlesheets4, # read_sheet 
               googledrive) # drive_ functions

#double check your wd. should be ../tempest_ionic_strength
#if not you need to do new relative file pathing

getwd()

## Set Github filepath for NPOC raw data files:

directory = "./Data/Processed Data/DOC"

# 1. Read in DOC data to calculate theoretical DOC concentrations in extracts -------

doc_actual <- read_csv(paste0(directory, "/ISTMP_NPOC_TDN_L1.csv")) %>% 
  select(sample_name, doc_mg_l) %>% #only need these columns
  mutate(sample_name = stringr::str_c("TMP_",sample_name , "_EXP")) #make the sample names match what's on the SPE vials.

duplicates <- doc_actual %>% subset(duplicated(sample_name))

View(duplicates)

readme <- readxl::read_excel("./SPE/TEMPEST_EXP_SPE_volumes_metadata_all.xlsx") %>%
  rename(sample_name = Sample_ID)

duplicates <- readme %>% subset(duplicated(sample_name))

View(duplicates)

#2. Read in Recoveries DOC data ---------

doc_recoveries <- read_csv(paste0(directory, "/Recoveries_NPOC_L1.csv")) %>% 
  select(sample_name, doc_mg_l, sample_vol, total_vol) %>% #only need these columns
rename(spe_doc_mg_l = doc_mg_l)
  
duplicates <- doc_recoveries %>% subset(duplicated(sample_name))

View(duplicates)


# 3. Do some math things to get theoretical DOC and compute the recovery ---------

#The dilution factor = (Volume used to redissolved dried down DOC * Methanol Elution Volume) / (Volume of methanol used to dry down * the volume of sample that was SPEd)
# The recovery = (The recovered DOC concentration * Dilution factor * 100) / Ambient NPOC concentration

doc_rec_cal <- doc_actual %>%
  inner_join(readme, by = "sample_name") %>%
  inner_join(doc_recoveries, by = "sample_name") %>%
  mutate(dil_fac = (total_vol * Actual_Vol_MeOH_mL) / (sample_vol * Actual_Vol_SPEd_mL),
         doc_recovery_per = (spe_doc_mg_l * dil_fac * 100) / doc_mg_l)

#these look insane but...
write_csv(doc_rec_cal, "../tempest_ionic_strength/Data/SPE_recoveries_calculated.csv")




