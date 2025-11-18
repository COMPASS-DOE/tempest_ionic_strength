library(tidyverse)
library(readxl)

# averages from Ed's spreadsheet: `Column_data.xlsx` and `Control column redo.xlsx`

colex_data_avgs_fromed <- read_excel("~/GitHub/tempest_ionic_strength/Data/Exp 3 Column Experiment ColEx Data/ColEx_data_avgs_formatted.xlsx")

#Bring in CDOM:
colex_eems <- read_excel("~/GitHub/tempest_ionic_strength/Data/Exp 3 Column Experiment ColEx Data/CDOM/24-10-30_TMP_ColEx_CDOM/7_Spectral indices/TMP_ColEx_SpectralIndices_combined.xlsx", sheet= "fluor_indices_DOC")

colex_abs <- read_excel("~/GitHub/tempest_ionic_strength/Data/Exp 3 Column Experiment ColEx Data/CDOM/24-10-30_TMP_ColEx_CDOM/7_Spectral indices/TMP_ColEx_SpectralIndices_combined.xlsx", sheet= "abs_indices")

colex_cdom <- colex_eems %>% 
  full_join(colex_abs, by=c("index", "analysis_date", "ID", "tmt", "wash", "soil_column")) %>%
  rename(Treatment = tmt,
         Timepoint = wash) %>%
  mutate(Treatment = case_when(Treatment == "ASW" ~ "SW" ,
                               TRUE ~ Treatment)) %>%
  filter(Treatment %in% c("SW", "FW")) %>%
  select(Treatment, Timepoint, pB:a254)

colex_cdom_averages <- colex_cdom %>%
  group_by(Treatment, Timepoint) %>%
  reframe(across(everything(),
      list(avg = mean, sd = sd),
      .names = "{.col}_{.fn}"))

colex_control_cdom <- read_excel("~/GitHub/tempest_ionic_strength/Data/Exp 3 Column Experiment ColEx Data/2025_control_col_data/CDOM/7_Spectral indices/SpectralIndices_combined_sample_only_DOC_norm.xlsx") %>%
  rename(Treatment = Exp_Type,
         Timepoint = wash) %>%
  select(Treatment, Timepoint, SUVA254:BIX) %>%
  filter(Timepoint != "CondSol") %>%
  filter(Timepoint != "Cmix") %>%
  mutate(Timepoint = as.numeric(Timepoint))

colex_cdom_control_averages <- colex_control_cdom %>%
  group_by(Treatment, Timepoint) %>%
  reframe(across(everything(),
                 list(avg = mean, sd = sd),
                 .names = "{.col}_{.fn}"))


colex_avgs_cdom_all <- colex_cdom_averages %>%
  full_join(colex_cdom_control_averages) %>%
  ungroup()

#merge them all together: 

colex_avgs_all <- colex_data_avgs_fromed %>%
  full_join(colex_avgs_cdom_all, by=c("Treatment", "Timepoint")) %>%
  mutate(across(where(is.numeric), ~ round(., 2)))

write_csv(colex_avgs_all, "../tempest_ionic_strength/Data/Exp 3 Column Experiment ColEx Data/TEMPEST_Col_Ex_all_chem_averaged_L2.csv")


