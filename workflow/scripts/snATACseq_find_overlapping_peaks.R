#--------------------------------------------------------------------------------------
#
#    Find overlapping snATACseq peaks
#
#--------------------------------------------------------------------------------------

##  Resources  ------------------------------------------------------------------------

# https://www.biostars.org/p/430015/
# https://www.biostars.org/p/453725/

##  Load Packages  --------------------------------------------------------------------
library(tidyverse)
library(ChIPpeakAnno)
library(cowplot)
library(VennDiagram)
library(readxl)
library(rmarkdown)

##  Set variables  --------------------------------------------------------------------
PUBLIC_DIR <- '~/Desktop/fetal_brain_snATACseq_V3_010323/resources/public_datasets/'
ZIFFRA_DIR <- paste0(PUBLIC_DIR, 'ziffra_2021/')
MSCOFF_DIR <- paste0(PUBLIC_DIR, 'markenscoff_2020/')
PEAKS_DIR <- '~/Desktop/fetal_brain_snATACseq_V3_010323/results/05PEAKS/'
SCRIPT_DIR <- '~/Desktop/fetal_brain_snATACseq_V3_010323/workflow/scripts/'
RESULTS_DIR <- '~/Desktop/fetal_brain_snATACseq_V3_010323/results/'

# Load functions
source('~/Desktop/fetal_brain_snATACseq_V3_010323/workflow/scripts/snATACseq_functions.R')

# Load hg19 peaks
for (REGION in c('LGE', 'MGE', 'CGE')) {
  
  PEAKS <- read_tsv(paste0(PEAKS_DIR, REGION, '.hg19.ext250bp.bed'), 
                    col_names = c('chr', 'start', 'end')) %>%
    toGRanges(format = "BED", header = FALSE) 
  assign(paste0(tolower(REGION), '_hg19_peaks'), PEAKS)
  
  
}

# GE overlaps - Markenscoff 2020  -----------------------------------------------------
# Note that this was done on hg19!!!!
mscoff_peaks <- readxl::read_excel(paste0(MSCOFF_DIR, 'markenscoff_2020_supp_table_2.xlsx'), 
                                   sheet = 'S2B') %>%
  select(`OCR coordinates (hg19)`, overlaps_cge_peak, overlaps_lge_peak, overlaps_mge_peak) %>%
  dplyr::rename('peaks' = 'OCR coordinates (hg19)') %>%
  separate(col = peaks, sep = ':', c("chr", "region")) %>%
  separate(col = region, sep = '-', c("start", "end")) 

for (REGION in c('lge', 'cge', 'mge')) {
  
  cat(REGION)
  
  PEAKS <- mscoff_peaks %>%
    filter(get(paste0('overlaps_', REGION, '_peak')) == 1) %>%
    select(chr, start, end) %>%
    toGRanges(format = "BED", header = FALSE) 
  
  assign(paste0(REGION, '_mscoff_peaks'), PEAKS)
  
}


# Find peak overlaps
for (REGION in c('lge', 'cge', 'mge')) {
  
  PEAKS <- get(paste0(REGION, '_hg19_peaks'))
  MSCOFF_PEAKS <- get(paste0(REGION, '_mscoff_peaks'))
  
  OVERLAPS <- findOverlapsOfPeaks(PEAKS, MSCOFF_PEAKS, minoverlap = 100)
  VENN <- makeVennDiagram(OVERLAPS, minoverlap = 100)
  VENN2 <- pairwise_venn(VENN, toupper(REGION), 
                         paste0(toupper(REGION), '-Mscoff')) 
  
  UNIQUE_PEAKS <- peaks <- OVERLAPS$peaklist$PEAKS
  MSCOFF_UNIQ_PEAKS <- OVERLAPS$peaklist$MSCOFF_PEAKS
  
  assign(paste0(REGION, '_overlaps'), OVERLAPS)
  assign(paste0(REGION, '_venn1'), VENN)
  assign(paste0(REGION, '_venn2'), VENN2)
  assign(paste0(REGION, '_uniq_peaks'), UNIQUE_PEAKS)
  assign(paste0(REGION, '_mscoff_uniq_peaks'), MSCOFF_UNIQ_PEAKS)
  
  
}

all_venn_plot <- plot_grid(mge_venn2, lge_venn2, cge_venn2)


# Get unique peaks
mge_ziffra_hg38_peaks <- mge_hg38_overlaps$peaklist$ziffra_mge_peaks
mge_hg38_peaks <- mge_hg38_overlaps$peaklist$mge_hg38_peaks


# Run motif enrichment
mge_uniq_motif_plot <- run_motif_analysis(mge_uniq_peaks, BSgenome.Hsapiens.UCSC.hg19, 30)
lge_uniq_motif_plot <- run_motif_analysis(lge_uniq_peaks, BSgenome.Hsapiens.UCSC.hg19, 30)
cge_uniq_motif_plot <- run_motif_analysis(cge_uniq_peaks, BSgenome.Hsapiens.UCSC.hg19, 30)
cge_mscoff_uniq_motif_plot <- run_motif_analysis(cge_mscoff_uniq_peaks, BSgenome.Hsapiens.UCSC.hg19, 30)
lge_mscoff_uniq_motif_plot <- run_motif_analysis(lge_mscoff_uniq_peaks, BSgenome.Hsapiens.UCSC.hg19, 30)
mge_mscoff_uniq_motif_plot <- run_motif_analysis(mge_mscoff_uniq_peaks, BSgenome.Hsapiens.UCSC.hg19, 30)


##  Create markdown html  -------------------------------------------------------------
render(paste0(SCRIPT_DIR, 'snATACseq_peak_overlap_and_motif_enrich.Rmd'),
       output_file = paste0('snATACseq_peak_overlap_and_motif_enrich.html'),
       output_dir = paste0(RESULTS_DIR, '03MARKDOWN'))

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------