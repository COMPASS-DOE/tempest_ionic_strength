# Merging all column experiment and batch experiment data for TMP porewater manuscript analysis 

# 1. Setup ---------------------------------------------------------------------

# load packages
require(pacman)
pacman::p_load(tidyverse, # keep things tidy
               janitor, # useful for simplifying column names
               readxl # read excel 
               ) 

#double check your wd. should be ../tempest_ionic_strength
#if not you need to do new relative file pathing

getwd()

#clean out your environment 
rm(list = ls())

#2. Load in L2 Data ---------------------------------------------------------------------

exp1_asw <- read_csv("~/GitHub/tempest_ionic_strength/Data/Exp 1 Ionic Strength ASW Data/TEMPEST_Ionic_Strength_ASW_chem_L2.csv") %>%
  rename(conductivity_uSpercm = Cond_uScm,
         conductivity_uSpercm_stdev = Cond_stdev,
         doc_mgperL = doc_mg_l_mean,
         doc_mgperL_stdev = doc_mg_l_sd,
         Na_ppm = Na,
         Na_ppm_stdev = Na_stdev,
         K_ppm = K,
         K_ppm_stdev = K_stdev,
         Ca_ppm = Ca,
         Ca_ppm_stdev = Ca_stdev,
         Mg_ppm = Mg,
         Mg_ppm_stdev = Mg_stdev,
         Al_ppm = Al,
         Al_ppm_stdev = Al_stdev,
         Fe_ppm = Fe,
         Fe_ppm_stdev = Fe_stdev,
         S_ppm = S,
         S_ppm_stdev = S_stdev
  )

exp2_nacl <- readxl::read_excel("~/GitHub/tempest_ionic_strength/Data/Exp 2 Ionic Strength NaCl Data/TEMPEST_Ionic_Strength_NaCl_chem_L2.xlsx")

exp3_colex <- readxl::read_excel("~/GitHub/tempest_ionic_strength/Data/Exp 3 Column Experiment ColEx Data/ColEx_data_allreps_alltreatments_formatted.xlsx") 

#3. Merge! ---------------------------------------------------------------------

common_cols <- intersect(names(exp1_asw), names(exp2_nacl))

exp_1n2 <- exp1_asw %>%
  full_join(exp2_nacl, by = common_cols) %>%
  select(all_of(common_cols)) %>%
  rename(Timepoint = Wash) %>%
  mutate(Treatment = as.character(Treatment),
         Column = "Average of Reps A, B, C" )

common_cols_2  <- intersect(names(exp_1n2), names(exp3_colex))

exp_all <- exp_1n2 %>%
  full_join(exp3_colex) %>%
  select(all_of(common_cols_2),  Column) %>%
  relocate(Exp_Type, Treatment, Timepoint, Column) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

write_csv(exp_all, "~/GitHub/tempest_ionic_strength/BHM/Data/TEMPEST_lab_experiments_all_L2.csv")

  