#--------------------------------------------------------------------------------------
#
#    ArchR - Local testing
#
#--------------------------------------------------------------------------------------

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


# Change cell type names
getCellColData(archR_50) %>%
  as_tibble() %>%
  #select(Clusters_pred) %>%
  mutate(across(Clusters_pred, str_replace, 'GE', 'GE-InN')) %>%
  select(Clusters_pred) 

# QC ---------------------------------------------------------------------------
test_510 <- create_archR_tss_frag_plot(archR_50, 'GE_510')
test_611 <- create_archR_tss_frag_plot(archR_50, 'GE_611')
test_993 <- create_archR_tss_frag_plot(archR_50, 'GE_993')

frag_size_plot <- create_frag_size_plot(archR_50)
tss_plot <- create_tss_plot(archR_50)


# qc_boxplot <- getCellColData(archR_50) %>% 
#   as_tibble() %>%
#   mutate(Sample = str_replace(Sample, 'GE_', '')) %>%
#   mutate(Sample = factor(Sample, levels = c('993', '611', '510'))) %>%
#   dplyr::rename('N Fragments' = nFrags,
#                 'Reads In Peaks' = ReadsInPeaks,
#                 'TSS Enrichment' = TSSEnrichment,
#                 'Reads Prom.' = ReadsInPromoter,
#                 'Promoter ratio' = PromoterRatio) %>%
#   dplyr::select(Sample, 'N Fragments', 'Reads In Peaks', FRIP, 'TSS Enrichment', 'Reads Prom.', 'Promoter ratio') %>%
#   tidyr::pivot_longer(cols = c('N Fragments', 'Reads In Peaks', FRIP, 'TSS Enrichment', 'Reads Prom.', 'Promoter ratio')) %>%
#   mutate(name = factor(name, levels = c('N Fragments', 'Reads In Peaks', 'FRIP', 
#                                           'TSS Enrichment', 'Reads Prom.', 'Promoter ratio'))) %>%
#   ggplot(aes(factor(Sample), value, fill = Sample)) +
#   geom_boxplot(width = 0.5) +
#   facet_grid(~name, scales = 'free') +
#   theme_bw() +
#   theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
#         panel.border = element_rect(colour = "black", size = 0.8),
#         plot.title = element_text(hjust = 0.5, face = 'bold', size = 16),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.x  = element_text(colour = "#000000", size = 13, angle = 45, hjust = 1),
#         axis.text.y  = element_text(colour = "#000000", size = 13),
#         legend.position = "none",
#         panel.spacing = unit(1.2, "lines"),
#         strip.text.x = element_text(size = 12),
#         strip.background = element_blank()) +
#   scale_fill_manual(values = c("#FF6347", "#EFC000FF", "#0073C2FF")) +
#   scale_y_continuous(breaks = scales::pretty_breaks(n = 3), 
#                      labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
#   coord_flip() 
# 
# a_v <- plot_grid(test_510, test_611, test_993, ncol = 1, labels = 'AUTO', label_size = 28)
# b_h <- plot_grid(frag_size_plot, tss_plot, ncol = 2, labels = c('D', 'E'), label_size = 28, scale = c(0.9, 0.9))
# b_c <- plot_grid(b_h, qc_boxplot, ncol = 1, scale = c(0.9, 0.9, 0.9), labels = c('', 'F'), label_size = 28)
# final_plot <- plot_grid(a_v, b_c, ncol = 2, rel_widths = c(1, 2))


# qc_group_plot <- plot_grid(test_510, test_611, test_993, frag_size_plot, tss_plot, qc_boxplot,
#                            ncol = 3, labels = 'AUTO', label_size = 28, 
#                            scale = c(0.9, 0.9, 0.9, 0.9, 0.9, 0.9),
#                            align = 'tblr')
# 
# tiff(paste0("~/Desktop/test.tiff"), height = 30, width = 60, units='cm', 
#      compression = "lzw", res = 300)
# qc_group_plot
# dev.off()


plot_grid(test_510, test_611, test_993, frag_size_plot,
          ncol = 2, labels = 'AUTO', scale = c(0.9, 0.9 ,0.9 ,0.8),  label_size = 28)
  

# Figure 4  ---------------------------------------------------------------------------
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
        panel.border = element_rect(colour = "black", size = 0.8),
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

# Peak summary plot
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


top_plot <- plot_grid(archR_50_clusters, umaps_1_gene, peak_call_summary_plot, venn2_ge_ins,
          labels = 'AUTO', label_size = 24, scale = c(1, 0.9, 0.9, 0.9))

two_way_venns <- plot_grid(mge_2way_venn2, lge_2way_venn2, cge_2way_venn2, nrow = 1,
                           labels = c('E', '', ''), label_size = 28, 
                           scale = c(0.9, 0.9, 0.9))

plot_grid(top_plot, two_way_venns, ncol = 1, rel_heights = c(2,1))

# UMAPs by marker gene ----------------------------------------------------------------
for (REGION in c('mge', 'lge', 'cge', 'progenitor', 'all')) {
  
  MARKER_GENES <- get(paste0(toupper(REGION), '_GENES'))
  
  plot_UMAPs_by_marker_genes(archR_50, 'UMAP', MARKER_GENES)
  assign(paste0(REGION, '_umaps'), all_genes_UMAP, .GlobalEnv)
  
}

mge_umaps
lge_umaps
cge_umaps
progenitor_umaps


## Motif Footprinting  -  Chptr 14  ---------------------------------------------------
#  cat(paste0('\nRunning footprinting ... \n'))
# Loop to pull out cisbp codes for motifs of interest
cat(paste0('\nAdding motif annotations  ... \n'))
archR_50 <- addMotifAnnotations(ArchRProj = archR_50, motifSet = "cisbp", name = "Motif", force = TRUE)

cat(paste0('\nGenerating motif footprints for GE ... \n'))
motifPositions <- getPositions(archR_50)
motifPositions

cat(paste0('\nAdding marker peaks for each cell type ... \n'))
# Add marker peaks - returns summarisedExperiment with 6 assays
markersPeaks <- getMarkerFeatures(
  ArchRProj = archR_50, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters_pred",
  bias = c("TSSEnrichment", "log10(nFrags)"), # Correction fof diffs in data quality
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

cat(paste0('\nAssessing motif enrichment ... \n'))
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = archR_50,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
enrichMotifs
colnames(enrichMotifs) <- c("CGE-InN", "LGE-InN", "MGE-InN", "Progenitor")
as_tibble(unlist(enrichMotifs@assays@data), rownames = 'gene') %>%
  separate_wider_delim('gene', delim = ".", names = c(NA, "gene")) %>%
  pivot_longer(cols = c("CGE", "LGE", "MGE", "Progenitor")) %>%
  mutate(value = as.numeric(value)) %>%
 
  
  group_by(name) %>%
  arrange(value, .by_group = TRUE) %>%
  slice_max(value, n = 10)

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 10, transpose = TRUE)
ComplexHeatmap::draw(heatmapEM, 
                     heatmap_legend_side = "bot", 
                     annotation_legend_side = "bot")

rownames(enrichMotifs)

MOTIFS_recode_all <- vector()

for (MOTIF in rownames(enrichMotifs)) {
  
  MOTIFS_recode <- grep(MOTIF, rownames(enrichMotifs), value = TRUE)
  print(MOTIFS_recode)
  MOTIFS_recode_all <<- c(MOTIFS_recode_all, MOTIFS_recode)
  
}


markerMotifs <- unlist(lapply(MOTIFS_recode_all, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs

# Compute footprints
cat(paste0('Computing footprints ... \n'))
seFoot <- getFootprints(
  ArchRProj = archR_50,
  positions = motifPositions[markerMotifs],
  groupBy = "Clusters_pred"
)



# Plot footprint - plot = FALSE required to get grob object
cat(paste0('Plotting ... \n'))
footprint_grob <-  plotFootprints(
  seFoot = seFoot,
  ArchRProj = archR_50,
  normMethod = "Subtract",
  plotName = "Footprints_subtract_bias",
  addDOC = FALSE,
  smoothWindow = 5,
  plot = TRUE
  
)


for (MOTIF in unique(MOTIFS_recode_all)) {
  
print(ggplotify::as.ggplot(footprint_grob[[MOTIF]]))

}
                     