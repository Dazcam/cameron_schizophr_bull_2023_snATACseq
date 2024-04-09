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
library(Seurat) # For no legend

##  Set variables  --------------------------------------------------------------------
ROOT_DIR <- '~/Desktop/fetal_brain_snATACseq_V3_010323/'
RESULTS_DIR <- paste0(ROOT_DIR, 'results/')
PUBLIC_DIR <- paste0(ROOT_DIR, 'resources/public_datasets/')
MSCOFF_DIR <- paste0(PUBLIC_DIR, 'markenscoff_2020/')
PEAKS_DIR <- paste0(RESULTS_DIR, '05PEAKS/')
SCRIPT_DIR <- paste0(ROOT_DIR, 'workflow/scripts/')

TABLE_DIR <- paste0(RESULTS, '07TABLES/')

RUN_PEAK_MOTIFS <- FALSE
RUN_UNIQUE_PEAK_MOTIFS <- FALSE
RUN_MRKDN <- FALSE

# Load functions
source(paste0(SCRIPT_DIR, 'snATACseq_functions.R'))

# Load snATACseq peaks - hg19 and hg38
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
# Markenscoff data aligned to hg19
mscoff_peaks <- readxl::read_excel(paste0(MSCOFF_DIR, 'markenscoff_2020_supp_table_2.xlsx'), 
                                   sheet = 'S2B') %>%
  dplyr::select(`OCR coordinates (hg19)`, overlaps_cge_peak, overlaps_lge_peak, overlaps_mge_peak) %>%
  dplyr::rename('peaks' = 'OCR coordinates (hg19)') %>%
  separate(col = peaks, sep = ':', c("chr", "region")) %>%
  separate(col = region, sep = '-', c("start", "end")) 

for (REGION in c('lge', 'cge', 'mge')) {
  
  cat(REGION)
  
  PEAKS <- mscoff_peaks %>%
    filter(get(paste0('overlaps_', REGION, '_peak')) == 1) %>%
    dplyr::select(chr, start, end) %>%
    toGRanges(format = "BED", header = FALSE) 
  
  assign(paste0(REGION, '_mscoff_peaks'), PEAKS)
  
}


## Find peak overlaps - 2-way  -------
for (REGION in c('lge', 'cge', 'mge')) {
  
  SHI_PEAKS <- get(paste0(REGION, '_hg19_peaks'))
  MSCOFF_PEAKS <- get(paste0(REGION, '_mscoff_peaks'))
  
  OVERLAPS <- findOverlapsOfPeaks(SHI_PEAKS, MSCOFF_PEAKS, minoverlap = 100)
  VENN <- makeVennDiagram(OVERLAPS, minoverlap = 100)
  VENN2 <- pairwise_venn(VENN, paste0(toupper(REGION), '-N'), 
                         paste0(toupper(REGION), '-Bulk')) 
  
  UNIQUE_PEAKS <- OVERLAPS$peaklist$SHI_PEAKS
  MSCOFF_UNIQ_PEAKS <- OVERLAPS$peaklist$MSCOFF_PEAKS
  
  assign(paste0(REGION, '_2way_overlaps'), OVERLAPS)
  assign(paste0(REGION, '_2way_venn1'), VENN)
  assign(paste0(REGION, '_2way_venn2'), VENN2)
  assign(paste0(REGION, '_2way_uniq_peaks'), UNIQUE_PEAKS)
  assign(paste0(REGION, '_2way_mscoff_uniq_peaks'), MSCOFF_UNIQ_PEAKS)
  
  
}


## Create bed for unique peaks  -----------
dir.create(paste0(PEAKS_DIR, 'UNIQUE_PEAKS/'))
for (region in c('lge', 'cge', 'mge')) {
  
  bed <- valr::gr_to_bed(get(paste0(region, '_2way_overlaps'))$uniquePeaks) %>%
    mutate(names = names(get(paste0(region, '_2way_overlaps'))$uniquePeaks)) %>%
    filter(str_detect(names, 'SHI_PEAKS')) %>%
    dplyr::select(!names) %>%
    write_tsv(paste0(PEAKS_DIR, 'UNIQUE_PEAKS/', toupper(region), '.vs.MSCOFF_unique_peaks.hg19.ext250bp.bed'))
    message('Region peak cnt: ', nrow(bed))
}

all.peaks <- OVERLAPS$all.peaks
gr1.renamed <- all.peaks$gr1
gr2.renamed <- all.peaks$gr2
peakNames <- melt(ol$peaklist[['gr1///gr2']]$peakNames, value.name="merged.peak.id")
gr1.sub <- gr1.renamed[peakNames[grepl("^gr1", peakNames[, 3]), 3]]
gr2.sub <- gr2.renamed[peakNames[grepl("^gr2", peakNames[, 3]), 3]]

# Find GE InNs peak overlaps hg38 - 3-way
cat('\nFinding overlaps for GE InNs ... \n\n')
overlaps_ge_ins <- findOverlapsOfPeaks(cge_hg38_peaks, lge_hg38_peaks, mge_hg38_peaks, minoverlap = 100)
venn_ge_ins <- makeVennDiagram(overlaps_ge_ins, minoverlap = 100)
venn2_ge_ins <- threeway_venn(venn_ge_ins, 'CGE-N','LGE-N','MGE-N') 
all_venn_plot <- plot_grid(venn2_ge_ins, mge_2way_venn2, lge_2way_venn2, cge_2way_venn2, 
                           scale = c(1, 0.8, 0.8, 0.8))
all_venn_plot

# Run motif enrichment - takes about 40mins
if (RUN_UNIQUE_PEAK_MOTIFS) {
  
  
  mge_uniq_motif_plot <- run_motif_analysis(mge_2way_uniq_peaks, 'hg19', 30, 'mge_uniq')
  lge_uniq_motif_plot <- run_motif_analysis(lge_2way_uniq_peaks, 'hg19', 30, 'lge_uniq')
  cge_uniq_motif_plot <- run_motif_analysis(cge_2way_uniq_peaks, 'hg19', 30, 'cge_uniq')
  #cge_mscoff_uniq_motif_plot <- run_motif_analysis(cge_2way_mscoff_uniq_peaks, BSgenome.Hsapiens.UCSC.hg19, 30)
  #lge_mscoff_uniq_motif_plot <- run_motif_analysis(lge_2way_mscoff_uniq_peaks, BSgenome.Hsapiens.UCSC.hg19, 30)
  #mge_mscoff_uniq_motif_plot <- run_motif_analysis(mge_2way_mscoff_uniq_peaks, BSgenome.Hsapiens.UCSC.hg19, 30)

}

# Run motif enrichment - takes about 40mins
if (RUN_PEAK_MOTIFS) {
  
  
  mge_motif_plot <- run_motif_analysis(mge_hg38_peaks, 'hg38', 30, 'mge')
  lge_motif_plot <- run_motif_analysis(lge_hg38_peaks, 'hg38', 30, 'lge')
  cge_motif_plot <- run_motif_analysis(cge_hg38_peaks, 'hg38', 30, 'cge')
  progenitor_motif_plot <- run_motif_analysis(progenitor_hg38_peaks, 'hg38', 30, 'progenitor')
  
}

# Save motif dfs
motifs_list <- list(cge_motifs_df, lge_motifs_df, mge_motifs_df, progenitor_motifs_df,
                    cge_uniq_motifs_df, lge_uniq_motifs_df, mge_uniq_motifs_df)
names(motifs_list) <- c('cge_motifs_df', 'lge_motifs_df', 'mge_motifs_df', 'progenitor_motifs_df',
                        'cge_uniq_motifs_df', 'lge_uniq_motifs_df', 'mge_uniq_motifs_df')
openxlsx::write.xlsx(motifs_list, paste0(TABLE_DIR, 'snATACseq_motif_tables.xlsx'))

## Create final plots for motif enrichment  ----------
for (REGION in c('LGE', 'MGE', 'CGE', 'Progenitor')) {
  # https://stackoverflow.com/questions/76055759
  
  motifs_factor <- get(paste0(tolower(REGION), '_motifs_df')) %>%
    arrange(desc(enrich_log2), desc(negLog10Padj)) %>%
    na.omit() %>%
    top_n(n = 30, wt = enrich_log2) %>%
    arrange(enrich_log2) %>%
    pull(motifs)

  motifs_df <- get(paste0(tolower(REGION), '_motifs_df'))  %>%
    arrange(desc(enrich_log2), desc(negLog10Padj)) %>%
    na.omit() %>%
    top_n(n = 30, wt = enrich_log2) %>%
    pivot_longer(!motifs, names_to = "name", values_to = "value") %>%
    arrange(desc(value)) %>%
    mutate(motifs = fct_relevel(motifs, motifs_factor))
  
  df_split <- split(motifs_df, motifs_df$name)
  pal_drug <- c("negLog10Padj" = "Blues", "enrich_log2" = "Greens")
  
  if (REGION == 'Progenitor') {
    
    TITLE <- REGION
    
  } else {
    
    TITLE <- REGION %>%
      toupper() %>%
      str_replace(., 'GE', 'GE-N') 
    
  }

  motif_plot <- ggplot(motifs_df, aes(y = motifs, x = name, label = round(value, 3))) +
    purrr::imap(df_split, function(x, y) {
      if (y == 'enrich_log2') lim <- (2.5)  else lim <- (320)
      if (y == 'enrich_log2') legend_order <- 1 else legend_order <- 0
      if (y == 'enrich_log2') legend_title <- expression(log[2](Enrichment)) else legend_title <- expression(-log[10](P[adj]))
      list(
        geom_tile(data = x, aes(fill = value), colour = "black", linewidth = 0.3),
        scale_fill_distiller(
          palette = pal_drug[[y]], name = legend_title,
          direction = 1, limits = c(0, lim), guide = 'colourbar'
        ),
        #guides(fill = guide_colourbar(order = legend_order)),
        ggnewscale::new_scale_fill()
      )
    }) +
    #geom_text(color = "black") +
    theme_bw() +
    theme(
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      axis.ticks = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 13),
      axis.title = element_blank(),
      plot.title = element_text(size = 15, hjust = 0.5, face = 'bold'),
      legend.position = "right",
      legend.text = element_text(size = 11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    ) +
    scale_x_discrete(expand = c(0,0), position = "top") +
    ggtitle(TITLE) +
    coord_fixed(ratio = 0.7) 
  
  
  assign(paste0(tolower(REGION), "_top_motifs_plot"), motif_plot, .GlobalEnv)

}

legend <- get_legend(cge_top_motifs_plot)

plot_grid(cge_top_motifs_plot + NoLegend(),
          lge_top_motifs_plot + NoLegend(),
          mge_top_motifs_plot + NoLegend(),
          progenitor_top_motifs_plot + NoLegend(), legend,
          nrow = 1, scale = c(0.9, 0.9, 0.9, 0.9, 1))



##  Create markdown html  -------------------------------------------------------------
if (RUN_MRKDN == TRUE) {

    render(paste0(SCRIPT_DIR, 'snATACseq_peak_overlap_and_motif_enrich.Rmd'),
         output_file = paste0('snATACseq_peak_overlap_and_motif_enrich.html'),
         output_dir = paste0(RESULTS_DIR, '03MARKDOWN'))

}

# save.image(file='~/Desktop/motifs_env.RData')
load('~/Desktop/motifs_env.RData')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------