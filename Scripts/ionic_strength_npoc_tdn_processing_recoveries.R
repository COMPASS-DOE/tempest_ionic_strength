## This script imports raw data for NPOC measured using a Shimadzu TOC-L
## at PNNL MCRL and exports clean, Level 1 QC'ed data for recovery calculations
## Raw Data are read in from GitHub
## 
## Created: 2022-01-15 by Peter Regier for EXCHANGE
## Updated: 2022-06-26 by Allison Myers-Pigg for TEMPEST
## Updated: 2023-12-19 by AMP for ionic strength experiment
## 
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

## Set Github filepath for NPOC raw data files:

directory = "../tempest_ionic_strength/Data/DOC"

recovery_directory = "../tempest_ionic_strength/SPE/Metadata"

# 2. Functions -----------------------------------------------------------------

## Create a function to read in data
read_data <- function(data){
  # First, scrape date from filename
  rundate <- str_extract(data, "[0-9]{8}")
  # Second, read in data
  read_delim(file = data, skip = 10, delim = "\t") %>% 
    rename(sample_name = `Sample Name`, 
           npoc_raw = `Result(NPOC)`, 
           run_datetime = `Date / Time`) %>% 
    select(sample_name, npoc_raw, run_datetime) %>% 
    mutate(rundate = rundate)
}

read_mes <- function(readme){
  # First, scrape date from filename
  rundate <- str_extract(readme, "[0-9]{8}")
  # Second, read in Read Me
  readxl::read_excel(path = readme, sheet = 1) %>% 
    rename(sample_name = `Sample Name`,
           sample_vol = `Sample wt`,
           total_vol = `Total vol`) %>% 
    select(sample_name, Action) %>% 
    mutate(rundate = rundate)
}

recovery_dilutions <- function(recovery){
  # read in recovery info
  readxl::read_excel(path = recovery, sheet = 1) %>% 
    rename(short_name = Short_ID,
           sample_vol = mL_MeOH,
           total_vol = water_added_mL) %>% 
    mutate(sample_name = stringr::str_c("TMP_",short_name, "_EXP"))   %>%
  select(sample_name, Action, sample_vol, total_vol) 
}

# 3. Import data ---------------------------------------------------------------

## Create a list of files to download
files <- list.files(path = directory, pattern = "Summary.+RECOVERY", full.names = TRUE) 
files
ReadMes <- list.files(path = directory, pattern = "Readme.+RECOVERY", full.names = TRUE) 
ReadMes

recoveries <- list.files(path=recovery_directory, pattern = "recovery", full.names = TRUE )
recoveries

recoveries_raw <- recoveries %>% 
  map_df(recovery_dilutions) %>% 
  filter(grepl("\\.[a-zA-Z]\\.", sample_name)) %>% # filter to samples only
  filter(!grepl("omit", Action)) %>%
  select(-Action) %>%
  bind_rows() 

recoveries_raw_action <- recoveries %>% 
  map_df(recovery_dilutions) %>% 
  filter(grepl("\\.[a-zA-Z]\\.", sample_name)) %>% # filter to samples only
  filter(grepl("omit", Action)) %>%
  bind_rows() 

npoc_raw <- files %>% 
  map_df(read_data) %>% 
  filter(grepl("\\.[a-zA-Z]\\.", sample_name)) %>% # filter to samples only
  bind_rows() 

blanks_raw <- files %>% 
  map_df(read_data) %>% 
  filter(grepl("^Blank", sample_name)) %>% # filter to blanks only
  bind_rows() 

readmes_dilution_action <- ReadMes %>% 
  map_df(read_mes) %>% 
  filter(grepl("\\.[a-zA-Z]\\.", sample_name)) %>% # filter to samples only
  filter(grepl("ilution correction", Action)) %>%
  bind_rows() 

readmes_all <- ReadMes %>% 
  map_df(read_mes) %>% 
  filter(grepl("\\.[a-zA-Z]\\.", sample_name)) %>% # filter to samples only
  bind_rows() 

curvepts <- files %>% 
  map_df(read_data) %>% 
  filter(grepl("STD_0-50ppmNPOC", sample_name)) %>% # filter to the right curve only
  rename(standard_high_C = npoc_raw) %>%
  select(rundate,standard_high_C) %>%
  #this part of the code would matter if we actually ran different curve ranges between the two runs then applied the other curve to the other dataset (which is what was done). 
  #It doesn't matter now since functionally the same for the same concentration ranges. 
  pivot_longer(cols = c(standard_high_C)) %>%
  na.omit() %>%
  group_by(rundate) %>%
  distinct()%>%
  pivot_wider(names_from= name, values_from = value)%>%
  bind_rows()

#per the read me the second run used the first run's curve...
second_run_curve <- as.data.frame(rbind(c("20240201",50))) %>%
  rename(rundate = V1,
         standard_high_C = V2) %>%
  mutate(standard_high_C = as.double(standard_high_C))

curvepts <- curvepts %>% full_join(second_run_curve, by= c("rundate","standard_high_C"))

View(curvepts)
# 4. Calculate blanks and add to data ------------------------------------------

blanks <- blanks_raw %>% 
  filter(!run_datetime %in% NA) %>% 
  mutate(npoc_raw = ifelse(npoc_raw > 0, npoc_raw, NA)) %>%
  group_by(rundate) %>% 
  summarize(npoc_blank= round(mean(npoc_raw[!is.na(npoc_raw)]), 2),
            npoc_blank_SD= round(sd(npoc_raw[!is.na(npoc_raw)]), 2)) %>% #add SD columns
  select(rundate, npoc_blank, npoc_blank_SD)

View(blanks) # Check out the blank data 

### for recoveries, you need to blank correct

npoc_blk_corrected <- npoc_raw %>%
  inner_join(blanks, by = "rundate") %>% 
  mutate(npoc_corrected = npoc_raw - npoc_blank)

# 5. Flag sketch data -----------------------------------------------------------

npoc_flagged <- npoc_blk_corrected %>% 
  filter(grepl("\\.[a-zA-Z]\\.", sample_name)) %>% # filter to samples only
  inner_join(curvepts, by= "rundate") %>%
  mutate(
         #most curves only to 50, those samples were not above it. making 100 for the August and September, which used 0-100
         npoc_flag = case_when(npoc_raw > standard_high_C ~ "value above cal curve",
                               npoc_blank > 0.15*npoc_raw ~ "blank is ≥ 15% of sample value" # flagging if blank concentration is > 20% of the sample concentration
                               )
         # npoc_raw = case_when(npoc_flag == "incorrect sample naming, cannot resolve" ~ NA,
         #                      TRUE ~ npoc_raw),
         # tdn_raw = case_when(tdn_flag == "incorrect sample naming, cannot resolve" ~ NA,
         #                     TRUE ~ tdn_raw)
  )

# 6. Dilution Corrections ------------------------------------------------------
# total_volume = 7/.2
# 
 dilutions = 
   readmes_dilution_action %>% 
   left_join(recoveries_raw,  by="sample_name") %>%
#   mutate(Dilution =  total_vol/sample_vol) %>% 
   dplyr::select(rundate, sample_name, total_vol, sample_vol) %>% 
#   mutate(Dilution = case_when(is.na(Dilution) ~ total_volume,
#          TRUE ~ Dilution)) %>%
  force()

#samples_to_dilution_corrected = 
#  npoc_flagged %>%
 # left_join(dilutions, by = c("sample_name", "rundate")) %>% 
 # filter(grepl("ilution correction", Action)) %>%
 # filter(!Action %in% "Omit") %>% 
 # mutate(doc_mg_l= npoc_corrected * Dilution , # True concentration x dilution correction
  #       doc_mg_l = as.numeric(doc_mg_l), doc_mg_l = round(doc_mg_l, 2)
#         ) %>%
 # mutate(doc_mg_l = case_when(Dilution > 30 & npoc_flag == "blank is ≥ 15% of sample value" ~ NA,
  #                            TRUE ~ doc_mg_l), # removing values if high blanks and high dilution ratios, potentially large source of error. 
 #        npoc_flag = case_when(is.na(doc_mg_l) ~ "omitted for high dilution and blank values",
  #                             TRUE ~ npoc_flag)) # removing values if high blanks and high dilution ratios, potentially large source of error. 

all_samples_dilution_corrected =
  npoc_flagged %>%
  inner_join(readmes_all, by = c("sample_name", "rundate")) %>% 
  inner_join(dilutions, by = c("sample_name", "rundate")) %>%
  mutate(doc_mg_l = npoc_corrected) %>%
 # filter(!grepl("ilution correction", Action)) %>% 
  filter(!Action %in% "Omit") %>%
 # bind_rows(samples_to_dilution_corrected) %>%
  dplyr::select(sample_name, rundate, doc_mg_l, npoc_flag, sample_vol, total_vol)%>%
  mutate(doc_mg_l = if_else(doc_mg_l < 0, "NA", as.character(doc_mg_l)),
         doc_mg_l = as.numeric(doc_mg_l), doc_mg_l = round(doc_mg_l, 2))

#Identify if any duplicates were run, this should return an empty data frame if not:#

duplicates <- all_samples_dilution_corrected %>% subset(duplicated(sample_name))

View(duplicates)


## 7. Clean data ----------------------------------------------------------------


## Flagging data
npoc_flags <- all_samples_dilution_corrected %>% 
  ## add flags 
  # Below blank 
  mutate(npoc_flag = case_when(
                               doc_mg_l == 'NaN' ~ "value below blank",
                               grepl("value above cal curve",npoc_flag) ~ "value above cal curve",
                               TRUE ~ npoc_flag), 
         doc_mg_l = case_when(doc_mg_l == "NaN" ~ NA,
                     TRUE ~ doc_mg_l))

npoc_wmeta <- npoc_flags %>%
  mutate(sample_name = stringr::str_replace(sample_name,"O","0"), #fix weird sample names 
    Treatment = stringr::str_extract(sample_name, "\\d+(?=\\.)"),
         Wash = stringr::str_extract(sample_name, "(?<=\\.[a-zA-Z].)\\d+")
         ) %>%
  mutate(Treatment = case_when(Treatment == "01" ~ "0.1",
                               TRUE ~ Treatment),
         sample_vol = case_when(is.na(sample_vol) ~ 0.2,
                                TRUE ~sample_vol),
         total_vol = case_when(is.na(total_vol) ~ 7,
                                TRUE ~ total_vol),) 

View(npoc_wmeta)

# 8. Write data ----------------------------------------------------------------

#not sure the blank is >25% is staying to the end of this data frame 
write_csv(npoc_wmeta, "../tempest_ionic_strength/Data/Processed Data/DOC/Recoveries_NPOC_L1.csv")

treatment_order <- c('0','0.1','1','5', '25', '100')


npoc_wmeta %>%
  group_by(Treatment, Wash) %>%
  dplyr::summarise(mean_doc_mg_l = mean(doc_mg_l, na.rm=TRUE), sd_doc_mg_l = sd(doc_mg_l, na.rm=TRUE)) %>%
ggplot()+
  geom_pointrange(aes(x=Wash, y=mean_doc_mg_l, ymin = mean_doc_mg_l- sd_doc_mg_l, ymax = mean_doc_mg_l + sd_doc_mg_l, color= factor(Treatment, levels= treatment_order))) +
  geom_path(aes(x=Wash, y=mean_doc_mg_l, color= factor(Treatment, levels= treatment_order), group=factor(Treatment, levels= treatment_order) )) +
  theme_classic() +
  labs(x = "Wash", y = "DOC mgC/L", color = "% ASW")

