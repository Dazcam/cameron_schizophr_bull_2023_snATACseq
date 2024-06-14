#--------------------------------------------------------------------------------------
#
#    ArchR - Plots for paper
#
#--------------------------------------------------------------------------------------

##  Plots for snATACseq figure 4 and S14

## Load packages  ---------------------------------------------------------------------
if (!require("Require")) install.packages("Require")
Require::Require(c("tidyverse", "readxl", "ArchR", "Signac", "Seurat", 
                   "SeuratWrappers", "cowplot", "scCustomize", "rmarkdown", "SeuratDisk",
                   "ArchRtoSignac", "pheatmap", "ComplexHeatmap", "plyr", "gtools",
                   "BSgenome.Hsapiens.UCSC.hg38")) 
library(ArchR)
library(Seurat)

## Load env variables  ----------------------------------------------------------------
SCRIPT_DIR <- '~/Desktop/fetal_brain_snATACseq_V3_010323/workflow/scripts/'
RESULTS_DIR <- '~/Desktop/fetal_brain_snATACseq_V3_010323/results/'
ARCHR_DIR <- paste0(RESULTS_DIR, '02ARCHR/')
FRAGS_DIR <- paste0(RESULTS_DIR, '04FRAGMENT_FILES/')
PEAKS_DIR <- paste0(RESULTS_DIR, '05PEAKS/')
PLOTS_DIR <- paste0(RESULTS_DIR, '07PLOTS/')

# ATAC gene lists
MGE_GENES <- c('LHX6', 'NKX2-1', 'SFTA3', 'MAF', 
               'SLC32A1', 'GAD1', 'GAD2')
LGE_GENES <- c('MEIS2', 'FOXP1', 'FOXP2', 'ISL1', 'SERTAD4', 'ZNF503', 
               'SIX3', 'ZFHX3', 'PBX1', 'ZNF521', 
               'PDYN', 'EBF1', 'GRIA2', 'CNTNAP2', 'ZFHX3')
CGE_GENES <- c('SCGN','NR2F1', 'NFIB', 'PAX6', 'CALB2')
PROGENITOR_GENES <- c('HES1', 'TNC', 'GLI3', 'PHLDA1')
ALL_GENES <- c('SCGN','NR2F1', 'NFIB', 'PAX6',
               'FOXP1', 'FOXP2', 'ISL1', 'PBX1',  
               'LHX6', 'NKX2-1', 'SLC32A1', 'GAD2',
                'HES1', 'TNC', 'GLI3', 'PHLDA1')

ALL_GENES2 <- c("LHX6", "NKX2-1", "SFTA3", "MAF", "SLC32A1", "GAD1", "GAD2", 
               "SCGN", "NR2F1", "NFIB", "PAX6", "FOXP1", "FOXP2", "ISL1", "PBX1", 
               "HES1", "TNC", "GLI3", "PHLDA1", "MEIS2", "SERTAD4", "ZNF503", 
               "SIX3", "ZFHX3", "ZNF521", "PDYN", "EBF1", "GRIA2", "CNTNAP2", 
               "CALB2")


## Load functions  --------------------------------------------------------------------
source(paste0(SCRIPT_DIR, 'snATACseq_functions.R'))
source(paste0(SCRIPT_DIR, 'snATACseq_find_overlapping_peaks.R')) # For venns

# Load data ---------------------------------------------------------------------------
archR <- loadArchRProject(paste0(ARCHR_DIR, 'GE')) 
archR_50 <- loadArchRProject(paste0(ARCHR_DIR, 'GE_pred_id_50')) 

# Fig S14 ---------------------------------------------------------------------------
test_510 <- create_archR_tss_frag_plot(archR_50, 'GE_510')
test_611 <- create_archR_tss_frag_plot(archR_50, 'GE_611')
test_993 <- create_archR_tss_frag_plot(archR_50, 'GE_993')

frag_size_plot <- create_frag_size_plot(archR_50)
tss_plot <- create_tss_plot(archR_50)

# Figure 4  ---------------------------------------------------------------------------
## Clusters UMAP - 4A -----
archR_50$Clusters_pred <- str_replace_all(archR_50$Clusters_pred, 'GE', 'GE-N')
archR_50_clusters <- plotEmbedding(ArchRProj = archR_50, 
                                  colorBy = "cellColData", 
                                  name = "Clusters_pred", 
                                  embedding = "UMAP", 
                                  pal = c(`CGE-N` = '#6098ab',
                                          `MGE-N` = '#f18e2a', 
                                          `LGE-N` = '#e1575a', 
                                          Progenitor = '#58a14e'),
                                 # labelMeans = FALSE, # This removes Cluster labels altogether
                                  labelAsFactors = FALSE,
                                  labelSize = 4) +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.8),
        plot.title = element_text(hjust = 0.5, face = 'bold', size = 14),
        axis.title.x = element_text(colour = "#000000", size = 13, vjust = -0.5),
        axis.title.y = element_text(colour = "#000000", size = 13),
        axis.text.x  = element_text(colour = "#000000", size = 13),
        axis.text.y  = element_text(colour = "#000000", size = 13),
#        legend.key.size = unit(3, "line"),
        legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(size = 12),
        legend.title = element_blank()) + 
  ggtitle('') + xlab('UMAP_1') + ylab('UMAP_2') +
  NoLegend()
  #guides(color = guide_legend(override.aes = list(size = 3, shape = 15))) # Squares for 


## Marker gene UMAPs - 4B  -----
for (REGION in c('mge', 'lge', 'cge', 'progenitor')) {
  
  GET_MARKERS <- function(x) { # Instead of nested ifelse
    switch(x,
           "mge" = 'LHX6',
           "lge" = 'ZNF503',
           "cge" = 'NR2F1',
           "progenitor" = 'HES1',
           stop("Invalid SUFFIX value")
    )
    
  }
  
  MARKERS <- GET_MARKERS(REGION)
  
  PLOT <- plotEmbedding(
    ArchRProj = archR_50, 
    colorBy = "GeneScoreMatrix", 
    name = MARKERS, 
    embedding = 'UMAP',
    imputeWeights = getImputeWeights(archR_50),
    pal = c('#f7fcfd', '#e0ecf4', '#bfd3e6', '#9ebcda', '#8c96c6', '#8c6bb1', '#88419d', '#810f7c', '#4d004b'),
    plotAs = 'points',  # See streaking issue https://github.com/GreenleafLab/ArchR/issues/1731
  ) + ggtitle(MARKERS) 
  
  PLOT <- PLOT + 
    theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", size = 0.8),
          plot.title = element_text(hjust = 0.5, face = 'bold', size = 16),
          axis.title.x = element_text(colour = "#000000", size = 13, vjust = -0.5),
          axis.title.y = element_text(colour = "#000000", size = 13),
          axis.text.x  = element_text(colour = "#000000", size = 13),
          axis.text.y  = element_text(colour = "#000000", size = 13),
          legend.position = "right",
          legend.direction = "vertical",
          legend.text = element_text(size = 12), 
          legend.title = element_blank()) + 
    ggtitle(MARKERS) + xlab('UMAP_1') + ylab('UMAP_2')
  
  assign(paste0(REGION, '_umaps_1_gene'), PLOT, .GlobalEnv)
  
}

# UMAPs by marker gene
umaps_1_gene <- plot_grid(lge_umaps_1_gene, mge_umaps_1_gene, cge_umaps_1_gene, progenitor_umaps_1_gene)

## Peak summary plot - 4C -----
peak_call_summary <- metadata(archR_50@peakSet)$PeakCallSummary %>%
  mutate(Group = str_replace(Group, 'E\\(', 'E-N \\(')) %>%
  mutate(Group = str_replace(Group, 'or\\(', 'or \\('))

peak_call_summary_plot <- ggplot(peak_call_summary, 
                                 aes(fill=Var1, y=Freq, x=Group)) + 
  geom_bar(position="stack", stat="identity") +
  viridis::scale_fill_viridis(discrete = T) +

  theme_bw() + 
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
       # panel.grid.major = element_blank(), 
       # panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", size = 0.8),
        plot.title = element_text(hjust = 0.5, face = 'bold', size = 16),
        axis.title.x = element_text(colour = "#000000", size = 13, vjust = -0.5),
        axis.title.y = element_text(colour = "#000000", size = 13),
        axis.text.x  = element_text(colour = "#000000", size = 13, angle = 45, hjust=1),
        axis.text.y  = element_text(colour = "#000000", size = 13),
        legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(size = 12),
        legend.title = element_blank()) +
  xlab("") +
  ylab(expression("No. of Peaks"~(x10^{"3"}))) 


## PLOTS  -----
# Venns loaded from snATACseq_find_overlapping_peaks.R
top_plot <- plot_grid(archR_50_clusters, umaps_1_gene, peak_call_summary_plot, venn2_ge_ins,
          labels = 'AUTO', label_size = 24, scale = c(1, 0.9, 0.9, 0.9))

two_way_venns <- plot_grid(mge_2way_venn2, lge_2way_venn2, cge_2way_venn2, nrow = 1,
                           labels = c('E', '', ''), label_size = 28, 
                           scale = c(0.9, 0.9, 0.9))

Fig_4 <- plot_grid(top_plot, two_way_venns, ncol = 1, rel_heights = c(2,1))

tiff(paste0(FIG_DIR, "Fig_4.tiff"), height = 30, width = 30, units='cm', 
     compression = "lzw", res = 300)
Fig_4 
dev.off()

Fig_S14 <- plot_grid(test_510, test_611, test_993, frag_size_plot,
                     ncol = 2, labels = 'AUTO', scale = c(0.9, 0.9 ,0.9 ,0.8),  
                     label_size = 28)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------