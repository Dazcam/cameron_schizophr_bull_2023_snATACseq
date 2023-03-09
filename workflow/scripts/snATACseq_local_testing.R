#--------------------------------------------------------------------------------------
#
#    ArchR - Local testing
#
#--------------------------------------------------------------------------------------

## Load packages  ---------------------------------------------------------------------
library(ArchR)
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

## Load env variables  ----------------------------------------------------------------
SCRIPT_DIR <- '~/Desktop/fetal_brain_snATACseq_091222/workflow/scripts/'
BATCH_VARS <- c('Sample', 'Age', 'Prep', 'Region', 'Study')
REGION <- 'FC'
SET_FILTER <- 'ziffra_only_default'

MAX_CLUSTERS <- 6
VAR_FEATURES <- 25000
N_START <- 10

RUN_HARMONY <- TRUE
COVARIATES <- c('Sample', 'Prep')
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

if (REGION == "FC") {
  
  SAMPLES <- c("14510_PFC_ATAC", "14611_PFC_ATAC", "14993_PFC_ATAC_AGGR",
               "GW17_Cortex", "Cortex_GW18_Frozen", "Cortex_GW21_Fresh",
               "PFC_GW20_SC", "PFC_Twin34_Frozen", "Twin31_PFC")
  
  SAMPLE_IDs <- SAMPLES %>% 
    str_remove("14") %>% 
    str_remove("win|W") %>% 
    str_remove("_Frozen|_Fresh|_SC|_ATAC_AGGR|_ATAC") %>% 
    str_replace("Cortex", "CTX") %>%
    str_replace("510_PFC", "PFC_510") %>%
    str_replace("611_PFC", "PFC_611") %>%
    str_replace("993_PFC", "PFC_993") %>%
    str_replace("G17_CTX", "CTX_G17") %>%
    str_replace("T31_PFC", "PFC_T31")
  
  
} else {
  
  SAMPLES <- c("14510_WGE_ATAC", "14611_WGE_ATAC", "14993_WGE_ATAC", 
               "MGE_GW20_SC", "MGE_Twin34_Frozen")
  
  SAMPLE_IDs <- SAMPLES %>% 
    str_remove("win") %>% 
    str_remove("14") %>% 
    str_remove("_Frozen|_SC|_ATAC") %>%
    str_replace("MGE_GW20", "MGE_G20") %>%
    str_replace("510_WGE", "WGE_510") %>%
    str_replace("611_WGE", "WGE_611") %>%
    str_replace("993_WGE", "WGE_993") 
  
} 

## Load functions  --------------------------------------------------------------------
source(paste0(SCRIPT_DIR, 'snATACseq_functions.R'))

## Load and prep archR project  -------------------------------------------------------
archR <- loadArchRProject(path ='/Users/darren/Desktop/temp/snATACseq/FC')  

## Subset by cell names  --------------------------------------------------------------
archR$log10nFrags <- log10(archR$nFrags)

# Counts before trim
cnts_preTrim_df <- as_tibble(table(archR$Sample), .name_repair = make.names) %>%
  dplyr::rename("Sample" = "X","cnts_preTrim" = "n")

# Set hard trim levels for tss and frags
all_pass <- set_filter_params(SET_FILTER)
cellsPass <- archR$cellNames[all_pass]
archR.trim <- archR[cellsPass, ]

cnts_postTrim_df <- as_tibble(table(archR.trim$Sample), .name_repair = make.names) %>%
  dplyr::rename("Sample" = "X",
         "cnts_postTrim" = "n") %>%
  left_join(cnts_preTrim_df) %>%
  relocate(Sample, cnts_preTrim, cnts_postTrim) %>%
  janitor::adorn_totals("row")

##  Run Cluster analysis  -------------------------------------------------------------
run_cluster_analysis(ARCHR_ID = archR.trim)

create_archR_group_plot(archR.2, 'Clusters', 'UMAP')
create_umap_batch_plots(archR.2, 'Clusters', 'UMAP')
plot_UMAPs_by_marker_genes(archR.2, 'UMAP', MARKER_GENES)


##  Run Harmony -----------------------------------------------------------------------
if(RUN_HARMONY) {
  
  COVARIATES <- c('Sample', 'Prep')
  run_cluster_analysis(ARCHR_ID = archR.2, 
                       HARMONY_NAME = 'Harmony',
                       COVARIATE = COVARIATES)

  create_archR_group_plot(archR.3, 'Clusters_Harmony', 'UMAPHarmony')
  create_umap_batch_plots(archR.3, 'Clusters_Harmony', 'UMAPHarmony')
  plot_UMAPs_by_marker_genes(archR.3, 'UMAPHarmony', MARKER_GENES)

}

##  Run unconstrained integration   ----------------------------------------------------
if(RUN_INTEGRATION) {
  
  SEURAT_DIR <- '~/Desktop/fetal_brain_snRNAseq_110122/resources/R_objects/'
  seurat.pfc <- readRDS(paste0(SEURAT_DIR, 'seurat.pfc.final.rds'))
  run_unconstrained_intergration(archR.2, 'Clusters', 'UMAP', 'IterativeLSI', seurat.pfc)
  
  if(RUN_HARMONY) {
    
    run_unconstrained_intergration(archR.3, 'Clusters_Harmony', 
                                   'UMAPHarmony', 'Harmony', seurat.pfc)
    
  }
  
}

##  Run peak calling   -----------------------------------------------------------------
RUN_PEAK_CALLING <- TRUE
RUN_ZIFFRA_COMPARE <- TRUE
ENV_PATH <- "/Users/darren/opt/miniconda3/envs/p2_macs"
MACS_PATH <- "/Users/darren/opt/miniconda3/envs/p2_macs/bin/macs2"
reticulate::use_condaenv(ENV_PATH) 
new_cellIDs_frags_3.5 <- c("FC-RG", "FC-RG", "FC-InN", "FC-InN", "FC-InN", 
                           "FC-InN", "FC-InN", "FC-InN", "FC-OPC", "FC-IP", 
                           "FC-ExN-UL", "FC-ExN-UL", "FC-ExN-UL", "FC-ExN-UL", "FC-ExN-UL", 
                           "FC-ExN-DL", "FC-MG", "FC-ExN-UL", "FC-ExN-UL", "FC-ExN-UL",
                           "FC-ExN-DL", "FC-ExN-DL", "FC-Undef", "FC-Undef", "FC-Undef")

run_peak_calling(archR.2, 'Clusters',
                 new_cellIDs_frags_3.5, 
                 gtools::mixedsort(unique(archR.2$Clusters)),
                 MACS_PATH, 0.05, 250)

create_archR_group_plot(archR.4, 'Clusters_broad', 'UMAP')
create_umap_batch_plots(archR.4, 'Clusters_broad', 'UMAP')

## ArchR to Signac  -------------------------------------------------------------------
seurat_atac <- readRDS(paste0(R_DIR, 'seurat_atac.rds'))

# Run integration
signac_snRNAseq_integration(seurat_atac, 'Shi', 'seurat_atac', 'IterativeLSI')

# Subset based on integration prediction score
seurat_atac_75 <- subset(seurat_int_seurat_atac, prediction.score.max >= 0.75)
seurat_atac_60 <- subset(seurat_int_seurat_atac, prediction.score.max >= 0.60)
seurat_atac_50 <- subset(seurat_int_seurat_atac, prediction.score.max >= 0.50)

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

for (SCORE in c('50', '60', '75')) {
  
  SEURAT_OBJ <- get(paste0('seurat_atac_', SCORE))
  
  UMAP <- DimPlot(
    object = SEURAT_OBJ,
    group.by = 'predicted.id',
    label = TRUE,
    repel = TRUE) + NoLegend() + 
    ggtitle(paste0('scATAC-seq_', SCORE))
  
  SEURAT_META <- SEURAT_OBJ[[]]
  
  PRED_ID_SCORE <- table(SEURAT_META$predicted.id) 
  SAMPLE_CNT <- table(SEURAT_META$Sample) 
  
  CELL_IDs <- rownames(SEURAT_OBJ)
  
  assign(paste0('umap', SCORE), UMAP, .GlobalEnv)
  assign(paste0('pred_id_score_', SCORE), PRED_ID_SCORE, .GlobalEnv)
  assign(paste0('sample_cnt_', SCORE), SAMPLE_CNT, .GlobalEnv)
  assign(paste0('cell_ID_', SCORE), CELL_IDs, .GlobalEnv) 
  
} 

plot_grid(umap_archR, umap_pred_id, umap_pred_score, 
          umap50, umap60, umap75)

# Check prediction score distributions
hist(seurat_int_seurat_atac$prediction.score.max)
seurat_meta <- seurat_int_seurat_atac[[]]
seurat_meta %>%
  dplyr::select(Clusters, prediction.score.max) %>%
  group_by(Clusters) %>%
  summarise(mean_score = mean(prediction.score.max), 
            median_score= median(prediction.score.max))

##  Peak overlap  ---------------------------------------------------------------------
library(ChIPpeakAnno)
library(VennDiagram)
ZIFFRA_DIR <- "~/Desktop/fetal_brain_snATACseq_070222/resources/public_datasets/ziffra_2021/"
peak_list <- read_excel(paste0(ZIFFRA_DIR, 'Ziffra_2021_supp_tables_2_13.xlsx'), sheet = 'ST2 AllPrimaryPeaks') %>%
  dplyr::select(seqnames, start, end, peak_name) 
macs2_peak_list <- read_excel(paste0(ZIFFRA_DIR, 'Ziffra_2021_supp_tables_2_13.xlsx'), sheet = 'ST3 MACSpeaks_byCelltype') 

cat('\nRemove MGE peaks from Ziffra data ... \n')
ziffra_peaks_all <- macs2_peak_list %>% 
  inner_join(peak_list) %>%
  relocate(seqnames, start, end, peak_name) %>%
  rename(chr = seqnames)  %>%
  rename_with(~str_remove(., '_MACSpeaks')) %>%
  gather(cell_type, val, -peak_name, -chr, -start, -end) %>%
  filter(val == 1) %>%
  filter(!cell_type == 'MGE') %>%
  select(chr, start, end, peak_name) %>%
  unique()

PEAKS_UNION <- getPeakSet(archR.4)

ziffra_gr <- toGRanges(ziffra_peaks_all, format="BED", header=FALSE) 

FC_overlaps <- findOverlapsOfPeaks(PEAKS_UNION, ziffra_gr, minoverlap = 100)
FC_venn <- makeVennDiagram(FC_overlaps, minoverlap = 100)
FC_venn2 <- pairwise_venn(FC_venn, 'Cameron', 'Ziffra') 

ziffra_uniq_peaks <- FC_overlaps$peaklist$ziffra_gr
cameron_uniq_peaks <- FC_overlaps$peaklist$PEAKS_UNION


# Motif enrichment - Human species number is 9606
library(monaLisa)
library(JASPAR2020)
library(TFBSTools)

pwms <- getMatrixSet(JASPAR2020, opts = list(matrixtype = "PWM", species = 9606))

# Get sequence for peakset
cameron_uniq_peaks_seq <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, cameron_uniq_peaks)
ziffra_uniq_peaks_seq <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38, ziffra_uniq_peaks)

cameron_motif_obj <- calcBinnedMotifEnrR(seqs = cameron_uniq_peaks_seq,
                                         pwmL = pwms,
                                         background = "genome",
                                         genome = BSgenome.Hsapiens.UCSC.hg38,
                                         genome.regions = NULL, # sample from full genome
                                         genome.oversample = 2, 
                                         BPPARAM = BiocParallel::SerialParam(RNGseed = 42),
                                         verbose = TRUE)

ziffra_motif_obj <- calcBinnedMotifEnrR(seqs = ziffra_uniq_peaks_seq,
                                         pwmL = pwms,
                                         background = "genome",
                                         genome = BSgenome.Hsapiens.UCSC.hg38,
                                         genome.regions = NULL, # sample from full genome
                                         genome.oversample = 2, 
                                         BPPARAM = BiocParallel::SerialParam(RNGseed = 42),
                                         verbose = TRUE)

# Pull out most significant
cameron_motif_obj_sig <- assay(cameron_motif_obj, "negLog10Padj")[, 1] > 20.0
cameron_motif_obj_sig[is.na(cameron_motif_obj_sig)] <- FALSE

pval <- as_tibble(assay(cameron_motif_obj, "negLog10Padj")) %>% rownames_to_column() %>% rename(negLog10Padj = `1`)
enrich <- as_tibble(assay(cameron_motif_obj, "log2enr")) %>% rownames_to_column() %>% rename(enrich_log2 = `1`)
motifs <- as_tibble(rowData(cameron_motif_obj)$motif.name) %>% rownames_to_column() %>% rename(motifs = value)
cameron_motif_df <- pval %>% 
  left_join(enrich) %>%
  left_join(motifs) %>% 
  select(!rowname) %>%
  arrange(desc(negLog10Padj))

ziffra_motif_obj_sig <- assay(ziffra_motif_obj, "negLog10Padj")[, 1] > 318.0
ziffra_motif_obj_sig[is.na(ziffra_motif_obj_sig)] <- FALSE

pval <- as_tibble(assay(ziffra_motif_obj, "negLog10Padj")) %>% rownames_to_column() %>% rename(negLog10Padj = `1`)
enrich <- as_tibble(assay(ziffra_motif_obj, "log2enr")) %>% rownames_to_column() %>% rename(enrich_log2 = `1`)
motifs <- as_tibble(rowData(ziffra_motif_obj)$motif.name) %>% rownames_to_column() %>% rename(motifs = value)
ziffra_motif_df <- pval %>% 
  left_join(enrich) %>%
  left_join(motifs) %>%
  select(!rowname) %>% 
  arrange(desc(negLog10Padj))

cameron_motif_plot <- plotMotifHeatmaps(x = cameron_motif_obj[cameron_motif_obj_sig], 
                                        which.plots = c("log2enr", "negLog10Padj"), 
                                        width = 1.8, maxEnr = 2, maxSig = 100,
                                        show_seqlogo = TRUE, ) 

ziffra_motif_plot <- plotMotifHeatmaps(x = ziffra_motif_obj[ziffra_motif_obj_sig], 
                                        which.plots = c("log2enr", "negLog10Padj"), 
                                        width = 1.8, maxEnr = 2, maxSig = 350,
                                        show_seqlogo = TRUE) 

unlist(plot_list) <- list(cameron_motif_plot, ziffra_motif_plot)
plot_grid(plotlist = unlist(plot_list), ncol = 3)

save.image(file = 'Desktop/temp/snATAC_local_env_frags_3.5.RData')
load(file = 'Desktop/temp/snATAC_local_env_frags_3.5.RData')

# PEAKS_UNION_DF <- data.frame(seqnames=seqnames(PEAKS_UNION),
#                              starts=start(PEAKS_UNION)-1,
#                              ends=end(PEAKS_UNION),
#                              names=c(rep(".", length(PEAKS_UNION))),
#                              scores=c(rep(".", length(PEAKS_UNION))),
#                              strands=strand(PEAKS_UNION))
# 
# # Can I work out correlations?
meta_covariates <- getCellColData(archR.2, c('Sample', 'Age', 'Prep', 'Region', 'Study', 'nFrags')) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "cells")
lsi_10 <- getReducedDims(archR.2)[,1:10] %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "cells") %>%
  left_join(meta_covariates)

res <- round(cor(lsi_10 %>% select(!cells)), 2)

corrplot::corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = res, col = col, symm = TRUE)

##  Create markdown html  -------------------------------------------------------------
render(paste0(SCRIPT_DIR, 'snATACseq_local_testing.Rmd'),
       output_file = paste0('snATACseq_local_testing_ziffra_default', 
                            paste0(COVARIATES, collapse = '_'), 
                            '.html'),
       output_dir = '~/Desktop/temp/snATACseq/local_testing/')

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

       