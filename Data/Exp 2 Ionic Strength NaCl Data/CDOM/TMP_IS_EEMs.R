# Processing EEMs data from Horiba Aqualog

## Things to be done before running the code
# 1. Gather your raw data: export your data from Aqualog (3 files for each sample)
#    1) Abs Spectra Graphs (axis in log)
#    2) Waterfall Plot Blank
#    3) Waterfall Plot Sample
# 2. Complete 'sample metadata sheet.xlsx'
# 3. Create a separate folder with the metadata file and all raw data files (.dat)
# 4. Make sure you have the right working directory - check with getwd()

# Load libraries
library(stringr)
library(readxl) # Read in excel files
library(fewsdom) # Pre-process eems - Katie's package
library(eemR) # Pre-process eems
library(staRdom) # PARAFAC

# If you don't have the fewsdom package installed, use devtools:
#library(devtools)
#install_github("katiewampler/fewsdom")


## 1. Set up a project path
getwd()
wd <- "/Users/kimj704/Github/tempest_ionic_strength/Data/Exp 2 Ionic Strength NaCl Data/CDOM" # replace this with your working directory
setwd(wd)

prjpath <- wd # where the raw data files are


## 2. Read in data

# Read in metadata
metadata <- read_xlsx(file.path(prjpath, "sample_metadata_TMP_IS.xlsx"))
metadata


run_eems(prjpath = prjpath, meta_name = "sample_metadata_TMP_IS.xlsx",
         get_doc = F, rayleigh_mask = c(30, 60, 25, 25))#, replace_blank = "blank fil1 1s") # replace blank = identifier of the process blank
# if getting error, check if index column has blanks



## 3. Combine spectral indices data and clean up sample names
read_fluor_data <- function(data){
  read_xlsx(path = data, sheet = "fluor_indices_DOC")
}

read_abs_data <- function(data){
  read_xlsx(path = data, sheet = "abs_indices")
}


filesall <- list.files(path = prjpath, pattern = "SpectralIndices_CDOM", full.names = TRUE, recursive = TRUE) 

# abs + DOC_normalized fluor data
fluor_ind_DOC <- filesall %>% 
  map_df(read_fluor_data) %>% 
  bind_rows() 


data_sample <- filesall %>% 
  map_df(read_abs_data) %>% 
  bind_rows() %>%
  left_join(fluor_ind_DOC) %>%
  filter(!grepl("BLK|pre|pst|Pre|Pos", sample))


data_sample2 <- data_sample %>%
  mutate(sample = sub("\\.2s_1_2s|\\.1s_1_1s|\\.4s_1_4s| 2s_1_2s| 1s_1_1s", "", sample)) %>%
  mutate(sample = sub("\\.IS\\.|IS", "_IS_", sample)) %>%
  slice(mixedorder(sample))


write_xlsx(data_sample2, paste0(prjpath, "/combined_spectralindices_TMP_IS.xlsx"))






