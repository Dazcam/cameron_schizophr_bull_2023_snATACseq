#--------------------------------------------------------------------------------------
#
#     Prep public dataset peaks for sLDSR
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Prep Ziffra et al (2021) peak files for sLDSR

##  Load Packages  --------------------------------------------------------------------
library(readxl)
library(tidyverse)

##  Define global variables  -----------------------------------------------------------
cat('\nDefining variables ... \n')
ZIFFRA_DIR <- "~/Desktop/fetal_brain_snATACseq_V3_010323/resources/public_datasets/ziffra_2021/"
PEAK_DIR <-  "~/Desktop/fetal_brain_snATACseq_V3_010323/results/05PEAKS/ZIFFRA_PEAKS/"
dir.create(PEAK_DIR)

##  Load public data  -----------------------------------------------------------------
# peak list = all peaks key; mac2 all peaks per cell type; specific; specific peaks per cell type
cat('\nLoading public data ... \n')
options(scipen = 999) # required to prevent peak coords. being abbr. in sci' notation
peak_list <- read_excel(paste0(ZIFFRA_DIR, 'Ziffra_2021_supp_tables_2_13.xlsx'), sheet = 'ST2 AllPrimaryPeaks') %>%
  dplyr::select(seqnames, start, end, peak_name) 
macs2_peak_list <- read_excel(paste0(ZIFFRA_DIR, 'Ziffra_2021_supp_tables_2_13.xlsx'), sheet = 'ST3 MACSpeaks_byCelltype') 
#specific_peak_list <- read_excel(paste0(ZIFFRA_DIR, 'Ziffra_2021_supp_tables_2_13.xlsx'), sheet = 'ST4 Specificpeaks_byCelltype') 


## Set up pairwise comparison df  -----------------------------------------------------
# Create df for pairwise tests
cat('\nCreating pairwise df ... \n')
pairwise_df <- data.frame(bray  = c("FC.ExN", "FC.ExN", "FC.ExN", "FC.InN", "FC.RG", 
                                    "FC.MG", "FC.undef", "LGE.InN", "MGE.InN", "MGE.InN",
                                    "CGE.InN", "GE.RG", "GE.Proj"),
                          ziffra = c("earlyEN", "dlEN", "ulEN", "IN_MGE", "RG",
                                     "Microglia", "EndoMural", "IN_MGE", "IN_MGE", "MGE",
                                     "IN_CGE", "RG", "Insula")
)

## Munge public data into list of peaks per cell-type  --------------------------------
## Join key and cell type lists
cat('\nCreat Ziffra peaks into bedr format ... \n')
ziffra_peaks_all <- macs2_peak_list %>% 
  inner_join(peak_list) %>%
  relocate(seqnames, start, end, peak_name) %>%
  rename(chr = seqnames)  %>%
  rename_with(~str_remove(., '_MACSpeaks')) %>%
  select(chr, start, end, peak_name, dlEN, earlyEN, ulEN) %>%
  gather(cell_type, val, -peak_name, -chr, -start, -end) %>%
  filter(val == 1) %>%
  group_split(cell_type) %>%
  purrr::set_names(purrr::map_chr(., ~.x$cell_type[1])) %>%
  purrr::imap(~write_tsv(.x[,1:4], paste0(PEAK_DIR,  .y, '_ziffra_macs2peaks.bed'), 
                                   col_names = FALSE)) 

# ziffra_peaks_specific <- specific_peak_list %>% 
#   inner_join(peak_list) %>%
#   relocate(seqnames, start, end, peak_name) %>%
#   rename_with(~str_remove(., '_Specificpeaks')) %>%
#   gather(cell_type, val, -peak_name, -seqnames, -start, -end) %>%
#   filter(val == 1) %>%
#   group_split(cell_type) %>%
#   purrr::set_names(purrr::map_chr(., ~.x$cell_type[1])) %>%
#   purrr::imap(~write_tsv(.x[,1:4], paste0(PEAK_DIR,  .y, '_ziffra_specificpeaks_100UP_100DOWN.bed'), 
#                          col_names = FALSE))


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------