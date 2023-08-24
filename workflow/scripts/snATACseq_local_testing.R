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


## Load env variables  ----------------------------------------------------------------
SCRIPT_DIR <- '~/Desktop/fetal_brain_snATACseq_V3_010323/workflow/scripts/'
RESULTS_DIR <- '~/Desktop/fetal_brain_snATACseq_V3_010323/results/'
ARCHR_DIR <- paste0(RESULTS_DIR, '02ARCHR/')
FRAGS_DIR <- paste0(RESULTS_DIR, '04FRAGMENT_FILES/')
PEAKS_DIR <- paste0(RESULTS_DIR, '05PEAKS/')

REGION <- 'GE'
SET_FILTER <- NULL

MAX_CLUSTERS <- 6
VAR_FEATURES <- 25000
N_START <- 10

RUN_HARMONY <- TRUE
COVARIATES <- c('Sample')
RUN_INTEGRATION <- TRUE

MARKER_GENES <-  c('GAD1', 'GAD2', 'SLC32A1', 'GLI3', 'SLC17A7',
                   'TNC', 'PROX1', 'SCGN', 'LHX6', 'NXPH1',
                   'MEIS2','ZFHX3', 'SPI1', 'LHX8', 'ISL1', 'GBX2') 

SAMPLES <- c("14510_WGE_ATAC", "14611_WGE_ATAC", "14993_WGE_ATAC")
  
SAMPLE_IDs <- SAMPLES %>% 
  str_remove("14|_ATAC") %>% 
  str_replace("510_WGE", "GE_510") %>%
  str_replace("611_WGE", "GE_611") %>%
  str_replace("993_WGE", "GE_993") 
  

## Load functions  --------------------------------------------------------------------
source(paste0(SCRIPT_DIR, 'snATACseq_functions.R'))
source(paste0(SCRIPT_DIR, 'snATACseq_plots_LDSR_barplots.R'))

## Load and prep archR project  -------------------------------------------------------
archR <- loadArchRProject(paste0(ARCHR_DIR, 'GE'))  

## Get counts  -----------------------------------------------------------------------
cell_cnts <- as_tibble(table(archR$Sample), .name_repair = make.names) %>%
  dplyr::rename("Sample" = "X","cell_cnts_postQC" = "n")



##  Run Cluster analysis  -------------------------------------------------------------
# run_cluster_analysis(ARCHR_ID = archR) - done on hawk
create_archR_group_plot(archR, 'Clusters', 'UMAP')
group_plot_Clusters_raw <- group_plot_Clusters 
plot_UMAPs_by_marker_genes(archR, 'UMAP', MARKER_GENES)


##  Run Harmony -----------------------------------------------------------------------
if(RUN_HARMONY) {
  
  run_cluster_analysis(ARCHR_ID = archR, 
                       HARMONY_NAME = 'Harmony',
                       COVARIATE = COVARIATES)

  create_archR_group_plot(archR.3, 'Clusters_Harmony', 'UMAPHarmony')
  plot_UMAPs_by_marker_genes(archR.3, 'UMAPHarmony', MARKER_GENES)

}


##  Run peak calling   ----------------------------------------------------------------
# Need a peak matrix before running archR to Signac
ENV_PATH <- "/Users/darren/opt/miniconda3/envs/p2_macs"
MACS_PATH <- "/Users/darren/opt/miniconda3/envs/p2_macs/bin/macs2"
reticulate::use_condaenv(ENV_PATH) 
run_peak_calling()

## ArchR to Signac  -------------------------------------------------------------------
# To save hassle match archR$Sample IDs to the Cellranger sample out files
seurat_atac <- run_archR2Signac(archR_pks_Clusters, FRAGS_DIR)

# Run integration  --------------------------------------------------------------------
signac_snRNAseq_integration(seurat_atac, 'Shi', 'seurat_atac', 'IterativeLSI')

# Assess integration using integration prediction score  ------------------------------
# Subset based on integration prediction score - Going with prediction.score.max >= 0.50
seurat_atac_50 <- subset(seurat_int_seurat_atac, prediction.score.max >= 0.50)
seurat_atac_50 <- subset(x = seurat_atac_50, 
                         subset = predicted.id %in% c('LGE', 'MGE', 'CGE', 'progenitor'))
                                            
for (SCORE in c('50')) {
  
  SEURAT_OBJ <- get(paste0('seurat_atac_', SCORE))
  
  UMAP <- DimPlot(
    object = SEURAT_OBJ,
    group.by = 'predicted.id',
    label = TRUE,
    repel = TRUE) + NoLegend() + 
    ggtitle(paste0('scATAC-seq_', SCORE))
  
  SEURAT_META <- SEURAT_OBJ[[]]
  
  SEURAT_ID_BY_SAMP <- SEURAT_META %>%
    group_by(Sample) %>%
    count(predicted.id) %>%
    spread(Sample, n)
  
  PRED_ID_SCORE <- table(SEURAT_META$predicted.id) 
  SAMPLE_CNT <- table(SEURAT_META$Sample) 
  
  assign(paste0('umap', SCORE), UMAP, .GlobalEnv)
  assign(paste0('pred_id_score_', SCORE), PRED_ID_SCORE, .GlobalEnv)
  assign(paste0('sample_cnt_', SCORE), SAMPLE_CNT, .GlobalEnv)
  assign(paste0('seurat_meta_', SCORE), SEURAT_META, .GlobalEnv)
  assign(paste0('seurat_id_by_sample_', SCORE), SEURAT_ID_BY_SAMP, .GlobalEnv)
  
} 

# Plot
umap_archR <- DimPlot(
  object = seurat_atac,
  group.by = 'Clusters',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('archR')
umap_pred_id <- DimPlot(
  object = seurat_int_seurat_atac,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
umap_pred_score <- FeaturePlot(
  object = seurat_int_seurat_atac,
  features = 'prediction.score.max') + ggtitle('pred-score')

umap_pred_score <- plot_grid(umap_archR, umap_pred_id, umap_pred_score, 
          umap50, umap60, umap75)


# Check prediction score distributions
hist(seurat_int_seurat_atac$prediction.score.max)
seurat_meta <- seurat_int_seurat_atac[[]]
seurat_meta %>%
  dplyr::select(Clusters, prediction.score.max) %>%
  group_by(Clusters) %>%
  summarise(mean_score = mean(prediction.score.max), 
            median_score = median(prediction.score.max))

##  Subset ArchR object
# Note that ArchR has a # between sample and cell ID and Signac uses _
for (SCORE in c('50')) {
  
  cat('\nSubsetting ArchR object based on Signac int pred cut off:', SCORE, '\n')
  SEURAT_OBJ <- get(paste0('seurat_atac_', SCORE))

  archR_cells2keep <- as_tibble(Cells(SEURAT_OBJ)) %>%
    separate(value, into = c("a", "b"), sep = "(?<=[0-9])_") %>%
    unite(value, a:b, sep = '#') %>%
    pull(value)

  ARCHR_OBJ <- archR_pks_Clusters[archR_cells2keep, ]
  ARCHR_OBJ$predicted.id <- SEURAT_OBJ@meta.data %>%
    pull(predicted.id)

  # Re-cluster - run_cluster_analysis returns archR.2
  run_cluster_analysis(ARCHR_OBJ)
  create_archR_group_plot(archR.2, 'Clusters', 'UMAP')
  plot_UMAPs_by_marker_genes(archR.2, 'UMAP', MARKER_GENES)
  
  # Compare UMAPs
  UMAP_clusters <- plotEmbedding(ArchRProj = archR.2, colorBy = "cellColData", 
                                    name = 'Clusters', 
                                    embedding = 'UMAP',
                                    plotAs = 'points') + ggtitle('Clusters')
  
  UMAP_clusters_pred <- plotEmbedding(ArchRProj = archR.2, colorBy = "cellColData", 
                                         name = 'predicted.id', 
                                         embedding = 'UMAP',
                                         plotAs = 'points') + ggtitle('Clusters_predicted')
  
  UMAP_compare <- plot_grid(UMAP_clusters, UMAP_clusters_pred)
  
  assign(paste0('UMAP_compare_', SCORE), UMAP_compare, .GlobalEnv)
  assign(paste0('group_plot_', SCORE), group_plot_Clusters, .GlobalEnv)
  assign(paste0('all_genes_UMAP_Clusters_', SCORE), all_genes_UMAP, .GlobalEnv)
  assign(paste0('archR_', SCORE), archR.2, .GlobalEnv)
  
}

# Run peak calling only archR subsets
for (SCORE in c('50')) {
  
  cat('\nRunning peaks on ArchR object: ArchR_', SCORE, '\n', sep = '')
  ARCHR_OBJECT <- get(paste0('archR_', SCORE))
  
  old_names <- gtools::mixedsort(unique(ARCHR_OBJECT$Clusters))
  
  if (SCORE == 50) {
    
    new_names <- c('LGE', 'LGE', 'Progenitor', 'Progenitor', 'MGE',
                   'CGE', 'MGE', 'MGE', 'MGE', 'MGE', 'MGE')
    
    } else if (SCORE == 60) {
    
    new_names <- c('LGE', 'LGE', 'Progenitor', 'MGE', 'CGE',
                   'Undef', 'MGE', 'MGE', 'MGE')
    
    } else {
      
      new_names <- c('Progenitor', 'LGE', 'LGE', 'MGE',
                     'CGE', 'MGE', 'MGE', 'MGE')

  }
  
  run_peak_calling(ARCHR_OBJECT, 'Clusters', 'Clusters_pred',
                   new_names, old_names)
  
  assign(paste0('archR_', SCORE), archR_pks_Clusters_pred, .GlobalEnv)
  assign(paste0('peakCallParams_summary_df_', SCORE), peakCallParams_summary_df_Clusters_pred, .GlobalEnv)
  assign(paste0('peak_call_summary_plot_', SCORE), peak_call_summary_plot_Clusters_pred, .GlobalEnv)

  
}

# Need to save using dropCells here. See https://github.com/GreenleafLab/ArchR/issues/483
# devtools::install_github("immunogenomics/presto") issue was caused by not having presto
saveArchRProject(archR_50, paste0(ARCHR_DIR, 'GE_pred_id_50'), dropCells = TRUE)  

# Create bed files
for (CELL_TYPE in c('CGE', 'MGE', 'LGE', 'progenitor', 'union')) {
  
  if (CELL_TYPE == 'union') {
    
    PEAKS <- getPeakSet(archR_50)
    
  } else {
    
    # Load reproducible peak set for each cell-type
    PEAKS <- readRDS(paste0(ARCHR_DIR, 'GE_pred_id_50/PeakCalls/', 
                            CELL_TYPE, '-reproduciblePeaks.gr.rds'))
    
  }
  
  # Convert and write bed file 
  PEAKS_DF <- tibble(data.frame(seqnames=seqnames(PEAKS),
                         starts=start(PEAKS)-1,
                         ends=end(PEAKS),
                         names=c(rep(".", length(PEAKS))),
                         scores=c(rep(".", length(PEAKS))),
                         strands=strand(PEAKS))) %>%
    write_tsv(file = paste0(PEAKS_DIR, CELL_TYPE, '.hg38.ext250bp.bed'))
  
  # Assign Granges object 
  assign(paste0(CELL_TYPE, '_peaks'), PEAKS)
  
}

# Create union peak set for Neuronal cells only - see here: https://github.com/GreenleafLab/ArchR/discussions/2007
# Note that this will replace the union peak set stored in ArchR object if saved
get_neuronal_union_peaks <- FALSE
union_peaks <- getPeakSet(archR_50)
if (get_neuronal_union_peaks) {
  
  # Create and load separate ArchR project to make sure main project peak / union files are not overwritten
  saveArchRProject(archR_50, paste0(ARCHR_DIR, 'GE_pred_id_50_Nrn_union'), dropCells = TRUE)  
  archR_50_N_union <- loadArchRProject(paste0(ARCHR_DIR, 'GE_pred_id_50_Nrn_union')) 
  
  # Set progenitor to NA to consider only Neuronal cell types for union
  archR_50_N_union$neurons_only <- str_replace(archR_50_N_union$Clusters_pred, "Progenitor", NA_character_)
  
  # Add coverages
  archR_50_N_union <- addGroupCoverages(ArchRProj = archR_50_N_union, groupBy = 'neurons_only', force = TRUE)
  
  # Call peaks
  archR_50_N_union <- addReproduciblePeakSet(
    ArchRProj = archR_50_N_union, 
    groupBy = 'neurons_only', 
    pathToMacs2 = MACS_PATH,
    cutOff = 0.05, 
    extendSummits = 250)
  
  union_peaks_Ns <- getPeakSet(archR_50_N_union)
  
  # Convert and write bed file 
  union_peaks_Ns_df <- tibble(data.frame(seqnames = seqnames(union_peaks_Ns),
                         starts = start(union_peaks_Ns)-1,
                         ends = end(union_peaks_Ns),
                         names = c(rep(".", length(union_peaks_Ns))),
                         scores = c(rep(".", length(union_peaks_Ns))),
                         strands = strand(union_peaks_Ns))) %>%
    write_tsv(file = paste0(PEAKS_DIR, 'union_neurons.hg38.ext250bp.bed'))

  
}

## Motif enrichment   -----------------------------------------------------------------
archR_50 <- loadArchRProject(paste0(ARCHR_DIR, 'GE_pred_id_50'))  
archR_50 <- addMotifAnnotations(ArchRProj = archR_50, motifSet = "cisbp", name = "Motif")
markersPeaks <- getMarkerFeatures(
  ArchRProj = archR_50, 
  useMatrix = "PeakMatrix", 
  groupBy = 'predicted.id',
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon")

enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = archR_50,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

colnames(enrichMotifs) <- c(' CGE-N', ' LGE-N', ' MGE-N', ' Progenitor') # Add space to make plot look nice
rownames(enrichMotifs) <- rownames(enrichMotifs) %>% str_extract(".*(?=\\_)") 
rownames(enrichMotifs) <- paste0(" ", rownames(enrichMotifs)) # Add space to make plot look nice

  
heatmap_mat <- plotEnrichHeatmap(enrichMotifs, n = 10, transpose = TRUE, returnMatrix = T)
heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 10, transpose = TRUE)
heatmapEM <- Heatmap(heatmap_mat, name = "mat", col = paletteContinuous(set = "comet", n = 100),
        show_column_dend = FALSE, show_row_dend = FALSE, rect_gp = gpar(col = "grey", lwd = 1),
        heatmap_legend_param = list(title = NULL, legend_direction = "horizontal", 
                                    labels_gp = gpar(fontsize = 16),
                                    grid_height = unit(1, "cm"), 
                                    grid_width = unit(5, "mm")), 
        column_names_side = c("top"), row_names_gp = gpar(fontsize = 18))
ComplexHeatmap::draw(heatmapEM, 
                     heatmap_legend_side = "bot", 
                     annotation_legend_side = "bot")
         


## Peak coaccessibility  --------------------------------------------------------------
archR_50 <- run_peak_coaccesibility(archR_50, PEAKS_DIR)

track_plot <- plotBrowserTrack(
  ArchRProj = archR_50, 
  groupBy = "predicted.id", 
  geneSymbol = 'DLX1', 
  upstream = 50000,
  downstream = 50000,
  loops = getCoAccessibility(archR_50)
)

grid::grid.newpage()
grid::grid.draw(track_plot$DLX1)


plot_UMAPs_by_marker_genes(archR_50, 'UMAP', MARKER_GENES)

# Save 
saveArchRProject(archR_50, paste0(ARCHR_DIR, 'GE_pred_id_50'), dropCells = TRUE)  

save.image(file = '~/Desktop/fetal_brain_snATACseq_V3_010323/resources/R_obj/atac_workspace_50.RData')
load('~/Desktop/fetal_brain_snATACseq_V3_010323/resources/R_obj/atac_workspace_50.RData')

##  Create markdown html  -------------------------------------------------------------
render(paste0(SCRIPT_DIR, 'snATACseq_local_testing.Rmd'),
       output_file = paste0('snATACseq_local_testing.html'),
       output_dir = paste0(RESULTS_DIR, '03MARKDOWN'))
       

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

