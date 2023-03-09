#--------------------------------------------------------------------------------------
#
#    ArchR - Load default env variables and files
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------

#  Run the analysis up until cluster QC cell removal 

## Initialise R library  --------------------------------------------------------------
.libPaths( c( "/scratch/c.c1477909/R/library", .libPaths() ) )

##  Load Packages  --------------------------------------------------------------------
library(ArchR)
library(pheatmap)
library(tidyverse)
library(rmarkdown)
library(BSgenome.Hsapiens.UCSC.hg38) 
library(ComplexHeatmap)
library(clustree)
library(cowplot)
library(rmarkdown)
library(argparser)
library(plyr)
library(gtools)

## Parse region / set region variable -------------------------------------------------
cat('\nParsing args ... \n')
p <- arg_parser("\nRead brain region and output directory for snATACseq QC ... \n")
p <- add_argument(p, "region", help = "No brain region specified")
p <- add_argument(p, "data_dir", help = "No input data directory specified")
p <- add_argument(p, "archR_out_dir", help = "No ArchR output directory specified")
p <- add_argument(p, "markdown_file", help = "No markdown file path specified")
p <- add_argument(p, "report_dir", help = "No report output directory specified")
p <- add_argument(p, "report_file", help = "No report filename specified")
args <- parse_args(p)
print(args)

##  Define global variables  -----------------------------------------------------------
cat('\nDefining variables ... \n')
REGION <- args$region
DATA_DIR <- args$data_dir
OUT_DIR <- args$archR_out_dir
MARKDOWN_FILE <- args$markdown_file
REPORT_DIR <- args$report_dir
REPORT_FILE <- args$report_file
FRAGS_THRESH <- 3000
TSS_THRESH <- 4 
MAX_CLUSTERS <- 6
VAR_FEATURES <- 25000
N_START <- 10
BATCH_VARS <- c('Sample', 'Age', 'Prep', 'Region', 'Study')

addArchRThreads(threads = 24) # Set Hawk to 32 cores so 0.75 of total
addArchRGenome("hg38")

# Create ArchR output directory
cat('\nCreate output directory for Arch R project  ... \n')
dir.create(OUT_DIR, recursive = TRUE) # Required ArchR doesn't create this for you

# Loop to extract sample IDs 
cat('\nSet sample IDs  ... \n')
  
SAMPLES <- c("14510_WGE_ATAC", "14611_WGE_ATAC", "14993_WGE_ATAC")
  
SAMPLE_IDs <- SAMPLES %>% 
  str_remove("14") %>% 
  str_remove("_ATAC") %>%
  str_replace("510_WGE", "WGE_510") %>%
  str_replace("611_WGE", "WGE_611") %>%
  str_replace("993_WGE", "WGE_993") 


# Fix this not needed on hawk
LEVELS <- SAMPLE_IDs # For stacked barplots

  
MARKER_GENES <-  c('GAD1', 'GAD2', 'SLC32A1', 'GLI3', 'SLC17A7',
                   'TNC', 'PROX1', 'SCGN', 'LHX6', 'NXPH1',
                   'MEIS2','ZFHX3', 'SPI1', 'LHX8', 'ISL1', 'GBX2')
  

# Not needed anymore???
#cat(paste0('\nLoading Seurat object for ', REGION, ' ... \n'))
#seurat.obj <- readRDS("../resources/R_objects/seurat.wge.final.rds")
#seurat.obj$cellIDs <- gsub('GE-', '', seurat.obj$cellIDs)



#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

