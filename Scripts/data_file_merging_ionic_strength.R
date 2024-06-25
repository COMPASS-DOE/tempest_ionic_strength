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
  group_by(Exp_Type, Treatment, Wash) %>%
  mutate(Wash = as.numeric(Wash)) %>%
  dplyr::summarise(doc_mg_l_mean = mean(doc_mg_l, na.rm=TRUE), doc_mg_l_sd = sd(doc_mg_l, na.rm=TRUE)) %>%
  mutate(doc_mg_l_mean = round(doc_mg_l_mean,2),
         doc_mg_l_sd = round(doc_mg_l_sd, 2))

EEMs_ASW <- read_csv("../tempest_ionic_strength/Data/Processed Data/EEMs/EEMs_L1/20231114_TEMPEST_Corrected_RSU_Data.csv") %>%
  filter(grepl("EXP", Sample_ID)) %>%
  mutate(Treatment = stringr::str_extract(Sample_ID, "\\d+(?=\\.)"),
         Wash = stringr::str_extract(Sample_ID, "(?<=\\.)(\\d+)(?=_EXP)"),
         Exp_Type = "ASW") %>%
  mutate(Treatment = case_when(Treatment == "01" ~ "0.1",
                               TRUE ~ Treatment)) %>%
  mutate(Treatment = as.numeric(Treatment)) %>%
  mutate(Wash = as.numeric(Wash)) %>%
  select(Exp_Type, Treatment, Wash, everything(), -Sample_ID) %>%
  arrange(Treatment, Wash)

ABS_ASW <- read_csv("../tempest_ionic_strength/Data/Processed Data/EEMs/EEMs_L1/20231114_TEMPEST_SpectralIndices.csv") %>%
  filter(grepl("EXP", Sample_ID)) %>%
  mutate(Treatment = stringr::str_extract(Sample_ID, "\\d+(?=\\.)"),
         Wash = stringr::str_extract(Sample_ID, "(?<=\\.)(\\d+)(?=_EXP)"),
         Exp_Type = "ASW") %>%
  mutate(Treatment = case_when(Treatment == "01" ~ "0.1",
                               TRUE ~ Treatment)) %>%
  mutate(Treatment = as.numeric(Treatment)) %>%
  mutate(Wash = as.numeric(Wash)) %>%
  select(Exp_Type, Treatment, Wash, everything(), -Sample_ID) %>%
  arrange(Treatment, Wash)

icp_ASW <- readxl::read_excel("../tempest_ionic_strength/Data/Processed Data/ICP/ICP_L1/Salinity_ICP_results_formatted.xlsx") %>%
  rename(Wash = Rinse) %>%
  mutate(Exp_Type= "ASW") %>%
  mutate(Treatment = as.numeric(Treatment)) %>%
  mutate(Wash = as.numeric(Wash)) %>%
  select(Exp_Type, Treatment, Wash, everything()) %>%
  arrange(Treatment, Wash)

ph_cond_ASW <- readxl::read_excel("../tempest_ionic_strength/Data/Processed Data/cond_ph/Salinity_Cond_pH_results.xlsx") %>%
  select(`Exp-Type`:`pH-stdev`) %>%
  rename(Exp_Type = `Exp-Type`,
         Cond_stdev = `Cond-stdev`,
         pH_stdev = `pH-stdev`) %>%
  mutate(pH = as.numeric(pH)) %>%
  mutate(Treatment = as.numeric(Treatment)) %>%
  mutate(Wash = as.numeric(Wash)) %>%
  select(Exp_Type, Treatment, Wash, everything()) %>%
  arrange(Treatment, Wash)

asw_all_chem <-   ph_cond_ASW %>%
  full_join(doc_tdn_ASW, by= c("Exp_Type","Treatment","Wash")) %>%
  full_join(icp_ASW, by= c("Exp_Type","Treatment","Wash")) %>%
  full_join(EEMs_ASW, by= c("Exp_Type","Treatment","Wash")) %>%
  full_join(ABS_ASW, by= c("Exp_Type","Treatment","Wash")) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

write_csv(asw_all_chem, "../tempest_ionic_strength/Data/Processed Data/TEMPEST_Ionic_Strength_ASW_chem_L2.csv")
