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
library(ArchRtoSignac)

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
source(paste0(SCRIPT_DIR, 'snATACseq_signac_local_functions.R'))

## Load and prep archR project  -------------------------------------------------------
archR <- loadArchRProject(paste0(ARCHR_DIR, 'GE'))  

## Filter samples  --------------------------------------------------------------------
if (is.null(SET_FILTER)) {
  
  print('\nNo filtering specification set. Skipping filtering ...\n')
  
} else {
  
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

}

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


##  Run peak calling   -----------------------------------------------------------------
# Note running peak calling on default clusters first so new_cellIDs = old_cellIDs
# Need a peak matrix before running archR to Signac
ENV_PATH <- "/Users/darren/opt/miniconda3/envs/p2_macs"
MACS_PATH <- "/Users/darren/opt/miniconda3/envs/p2_macs/bin/macs2"
reticulate::use_condaenv(ENV_PATH) 
new_cellIDs <- gtools::mixedsort(unique(archR$Clusters))
run_peak_calling()

peakCallParams_summary_df_default <- peakCallParams_summary_df
peak_call_summary_plot_default <- peak_call_summary_plot

## ArchR to Signac  -------------------------------------------------------------------
# To save hassle match archR$Sample IDs to the Cellranger sample out files
run_archR2Signac(archR.4, FRAGS_DIR)

cat('\nRun ArchR2Signac ... \n')
seurat_atac <- ArchR2Signac(
  ArchRProject = archR.4,
  refversion = "hg38",
  fragments_dir = FRAGS_DIR,
  pm = pkm, 
  fragments_fromcellranger = "Yes", 
  fragments_file_extension = NULL,
  annotation = annotations 
)

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

       