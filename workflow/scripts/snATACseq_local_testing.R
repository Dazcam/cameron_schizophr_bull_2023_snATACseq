#--------------------------------------------------------------------------------------
#
#    ArchR - Local testing
#
#--------------------------------------------------------------------------------------

## Load packages  ---------------------------------------------------------------------
library(ArchR)
library(Seurat)
library(Signac)
library(pheatmap)
library(tidyverse)
library(rmarkdown)
library(ComplexHeatmap)
library(cowplot)
library(rmarkdown)
library(plyr)
library(gtools)
library(BSgenome.Hsapiens.UCSC.hg38) # For peak calling
library(readxl)
library(ArchRtoSignac)
#library(tictoc) # This causes an error with addReproduciblePeaks function
library(SeuratWrappers) # For RunFastMNN used during integration


## Load env variables  ----------------------------------------------------------------
SCRIPT_DIR <- '~/Desktop/fetal_brain_snATACseq_V3_010323/workflow/scripts/'
RESULTS_DIR <- '~/Desktop/fetal_brain_snATACseq_V3_010323/results/'
FRAGS_DIR <- paste0(RESULTS_DIR, '04FRAGMENT_FILES/')
ARCHR_DIR <- paste0(RESULTS_DIR, '02ARCHR/')
REGION <- 'GE'
SET_FILTER <- NULL

MAX_CLUSTERS <- 6
VAR_FEATURES <- 25000
N_START <- 10

RUN_HARMONY <- TRUE
COVARIATES <- c('Sample')
RUN_INTEGRATION <- TRUE

MARKER_GENES <-  c('SLC17A7', # ExN
                   'TLE3', 'LHX2', # UL ExN
                   'CRYM', 'FEZF2', # DL ExN
                   'GAD1', 'GAD2', 'SLC32A1', # InN
                   'GLI3', 'TNC', # RG
                   'C3', 'SPI1', # MG
                   'EOMES', # IP
                   'OLIG1', 'OLIG2', # OPC
                   'ITM2A') # Endothelial

SAMPLES <- c("14510_WGE_ATAC", "14611_WGE_ATAC", "14993_WGE_ATAC")
  
SAMPLE_IDs <- SAMPLES %>% 
  str_remove("14|_ATAC") %>% 
  str_replace("510_WGE", "GE_510") %>%
  str_replace("611_WGE", "GE_611") %>%
  str_replace("993_WGE", "GE_993") 
  

## Load functions  --------------------------------------------------------------------
source(paste0(SCRIPT_DIR, 'snATACseq_functions.R'))

## Load and prep archR project  -------------------------------------------------------
archR <- loadArchRProject(paste0(ARCHR_DIR, 'GE'))  

## Get counts  -----------------------------------------------------------------------
cell_cnts <- as_tibble(table(archR$Sample), .name_repair = make.names) %>%
  dplyr::rename("Sample" = "X","cell_cnts_postQC" = "n")

##  Run Cluster analysis  -------------------------------------------------------------
# run_cluster_analysis(ARCHR_ID = archR) - done on hawk
create_archR_group_plot(archR, 'Clusters', 'UMAP')
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
# Subset based on integration prediction score
seurat_atac_75 <- subset(seurat_int_seurat_atac, prediction.score.max >= 0.75)
seurat_atac_60 <- subset(seurat_int_seurat_atac, prediction.score.max >= 0.60)
seurat_atac_50 <- subset(seurat_int_seurat_atac, prediction.score.max >= 0.50)

for (SCORE in c('50', '60', '75')) {
  
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

plotEmbedding(ArchRProj = archR_50, colorBy = "cellColData", 
              name = CLUSTERS_ID, 
              embedding = UMAP_ID,
              plotAs = 'points') + ggtitle('Clusters')


##  Subset ArchR object
# Note that ArchR has a # between sample and cell ID and Signac uses _
for (SCORE in c('50', '60', '75')) {
  
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
for (SCORE in c('50', '60', '75')) {
  
  cat('\nRunning peaks on ArchR object: ArchR_', SCORE, '\n', sep = '')
  ARCHR_OBJECT <- get(paste0('archR_', SCORE))
  
  old_names <- gtools::mixedsort(unique(ARCHR_OBJECT$Clusters))
  
  if (SCORE == 50) {
    
    new_names <- c('LGE', 'LGE', 'Progenitor', 'Progenitor', 'MGE',
                   'MGE', 'MGE', 'MGE', 'MGE', 'CGE', 'Undef')
    
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

save.image(file = '~/Desktop/fetal_brain_snATACseq_V3_010323/resources/R_obj/atac_workspace.RData')
load('~/Desktop/fetal_brain_snATACseq_V3_010323/resources/R_obj/atac_workspace.RData')

##  Create markdown html  -------------------------------------------------------------
render(paste0(SCRIPT_DIR, 'snATACseq_local_testing.Rmd'),
       output_file = paste0('snATACseq_local_testing.html'),
       output_dir = paste0(RESULTS_DIR, '03MARKDOWN'))
       

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

       