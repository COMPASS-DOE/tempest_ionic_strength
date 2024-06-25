## Created: 2024-06-24 by AMP for ionic strength experiment

# 1. Setup ---------------------------------------------------------------------

# load packages
require(pacman)
pacman::p_load(tidyverse, # keep things tidy
               janitor, # useful for simplifying column names
               googlesheets4, # read_sheet 
               googledrive) # drive_ functions

#double check your wd. should be ../tempest_ionic_strength
#if not you need to do new relative file pathing

getwd()

# 2. Load in L1 Data ---------------------------------------------------------------------

doc_tdn_ASW <- read_csv("../tempest_ionic_strength/Data/Processed Data/DOC/DOC_L1/TEMPEST_ASW_NPOC_TDN_L1.csv") %>%
  group_by(Treatment, Wash) %>%
  dplyr::summarise(doc_mg_l_mean = mean(doc_mg_l, na.rm=TRUE), doc_mg_l_sd = sd(doc_mg_l, na.rm=TRUE)) %>%
  mutate(doc_mg_l_mean = round(doc_mg_l_mean,2),
         doc_mg_l_sd = round(doc_mg_l_sd, 2))
