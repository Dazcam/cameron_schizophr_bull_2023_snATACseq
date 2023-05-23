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
for (REGION in c('LGE', 'MGE', 'CGE', 'progenitor')) {
  
  PEAKS_HG19 <- read_tsv(paste0(PEAKS_DIR, REGION, '.hg19.ext250bp.bed'), 
                    col_names = c('chr', 'start', 'end')) %>%
    toGRanges(format = "BED", header = FALSE) 
  PEAKS_HG38 <- read_tsv(paste0(PEAKS_DIR, REGION, '.hg38.ext250bp.bed'), 
                         col_names = c('chr', 'start', 'end')) %>%
    toGRanges(format = "BED", header = FALSE) 
  
  assign(paste0(tolower(REGION), '_hg19_peaks'), PEAKS_HG19)
  assign(paste0(tolower(REGION), '_hg38_peaks'), PEAKS_HG38)
  
  
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


# Find peak overlaps - 2-way
for (REGION in c('lge', 'cge', 'mge')) {
  
  PEAKS <- get(paste0(REGION, '_hg19_peaks'))
  MSCOFF_PEAKS <- get(paste0(REGION, '_mscoff_peaks'))
  
  OVERLAPS <- findOverlapsOfPeaks(PEAKS, MSCOFF_PEAKS, minoverlap = 100)
  VENN <- makeVennDiagram(OVERLAPS, minoverlap = 100)
  VENN2 <- pairwise_venn(VENN, toupper(REGION), 
                         paste0(toupper(REGION), '-Mscoff')) 
  
  UNIQUE_PEAKS <- peaks <- OVERLAPS$peaklist$PEAKS
  MSCOFF_UNIQ_PEAKS <- OVERLAPS$peaklist$MSCOFF_PEAKS
  
  assign(paste0(REGION, '_2way_overlaps'), OVERLAPS)
  assign(paste0(REGION, '_2way_venn1'), VENN)
  assign(paste0(REGION, '_2way_venn2'), VENN2)
  assign(paste0(REGION, '_2way_uniq_peaks'), UNIQUE_PEAKS)
  assign(paste0(REGION, '_2way_mscoff_uniq_peaks'), MSCOFF_UNIQ_PEAKS)
  
  
}


# Find GE InNs peak overlaps hg38 - 3-way
cat('\nFinding overlaps for GE InNs ... \n\n')
overlaps_ge_ins <- findOverlapsOfPeaks(cge_hg38_peaks, lge_hg38_peaks, mge_hg38_peaks, minoverlap = 100)
venn_ge_ins <- makeVennDiagram(overlaps_ge_ins, minoverlap = 100)
venn2_ge_ins <- threeway_venn(venn_ge_ins, 'CGE-InN','LGE-InN','MGE-InN') 
all_venn_plot <- plot_grid(venn2_ge_ins, mge_2way_venn2, lge_2way_venn2, cge_2way_venn2, 
                           scale = c(1, 0.8, 0.8, 0.8))
all_venn_plot


# Run motif enrichment
mge_uniq_motif_plot <- run_motif_analysis(mge_2way_uniq_peaks, BSgenome.Hsapiens.UCSC.hg19, 30, 'mge')
lge_uniq_motif_plot <- run_motif_analysis(lge_2way_uniq_peaks, BSgenome.Hsapiens.UCSC.hg19, 30, 'lge')
cge_uniq_motif_plot <- run_motif_analysis(cge_2way_uniq_peaks, BSgenome.Hsapiens.UCSC.hg19, 30, 'cge')
#cge_mscoff_uniq_motif_plot <- run_motif_analysis(cge_2way_mscoff_uniq_peaks, BSgenome.Hsapiens.UCSC.hg19, 30)
#lge_mscoff_uniq_motif_plot <- run_motif_analysis(lge_2way_mscoff_uniq_peaks, BSgenome.Hsapiens.UCSC.hg19, 30)
#mge_mscoff_uniq_motif_plot <- run_motif_analysis(mge_2way_mscoff_uniq_peaks, BSgenome.Hsapiens.UCSC.hg19, 30)


##  Create markdown html  -------------------------------------------------------------
render(paste0(SCRIPT_DIR, 'snATACseq_peak_overlap_and_motif_enrich.Rmd'),
       output_file = paste0('snATACseq_peak_overlap_and_motif_enrich.html'),
       output_dir = paste0(RESULTS_DIR, '03MARKDOWN'))

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------