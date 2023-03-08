## Load local packages for Signac processing packages  --------------------------------
library(Signac)
library(Seurat)
library(tidyverse)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(data.table)
library(cowplot)
library(GenomicRanges)
library(future)
library(tictoc)
library(SeuratWrappers) # For RunFastMNN
set.seed(1234)

# Set multiprocessing for dataset merging
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 15000 * 1024^2)

# Set global variables
SAMPLE_DIR <- '~/Desktop/fetal_brain_snATACseq_V3_010323/results/01CELLRANGER/'
SIGNAC_DIR <- '~/Desktop/fetal_brain_snATACseq_V3_010323/results/03SIGNAC/'
SAMPLES <- c("14510_WGE_ATAC", "14611_WGE_ATAC", "14993_WGE_ATAC")
SAMPLE_IDs <- c("GE_510", "GE_611", "GE_993")

GE_MARKER_GENES <-  c('SLC17A7', # ExN
                      'LHX6', "PLS3", "NXPH1", "SFTA3", "SOX6", # MGE
                      'PROX1', "NFIB", "PDZRN3", "AP1S2", "CALB2", "SCGN",
                      "PCDH9", "KLHL35", "ANKS1B", # CGE
                      'SIX3', 'TSHZ1', "ZNF503", "SERTAD4", "ISL1", # LGE
                      'GAD1', 'GAD2', 'SLC32A1', # InN
                      'GLI3', 'TNC', # RG
                      'C3', 'SPI1', # MG
                      'EOMES', # IP
                      'OLIG1', 'OLIG2', # OPC
                      'ITM2A') # Endothelial

cat('\nAvailable Vars loaded from signac_local_env.R: \n\nSAMPLES:', SAMPLES, 
    '\n\nSAMPLE_IDs:', SAMPLE_IDs, 
    '\n\nSAMPLE_DIR:', SAMPLE_DIR, 
    '\n\nSIGNAC_DIR:', SIGNAC_DIR, 
    '\n\nGE_MARKER_GENES:', GE_MARKER_GENES, '\n\n')


