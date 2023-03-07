#--------------------------------------------------------------------------------------
#
#    ArchR - Signac local testing - GE
#
#--------------------------------------------------------------------------------------

## Load packages  ---------------------------------------------------------------------
SCRIPTS_DIR <- '~/Desktop/fetal_brain_snATACseq_V3_010323/workflow/scripts/'
SAMPLE_DIR <- '~/Desktop/fetal_brain_snATACseq_V3_010323/results/01CELLRANGER/'
SIGNAC_DIR <- '~/Desktop/fetal_brain_snATACseq_V3_010323/results/02SIGNAC/'
source(paste0(SCRIPTS_DIR, 'snATACseq_signac_local_env.R'))
source(paste0(SCRIPTS_DIR, 'snATACseq_signac_local_functions.R'))

## Run Signac on default settings -----------------------------------------------------
SAMPLES <- "14611_WGE_ATAC"
SAMPLE_IDs <- "GE_611"
default_signac_processing('14611_WGE_ATAC', 
                          SAMPLE_DIR,
                          GE_MARKER_GENES)

DimPlot(object = seurat_obj_GE_611 , label = TRUE) + 
  NoLegend() + 
  ggtitle(SAMPLE_IDs)


## Integration with RNA object
# Load the pre-processed scRNA-seq data for seurat.objs


source(paste0(SHI_DIR, 'snRNAseq_GE_prep_shi_data_for_Seurat.R'))
seurat.shi <- CreateSeuratObject(counts = shi_data, meta.data = shi_meta)
seurat.shi <- NormalizeData(seurat.shi)
seurat.shi <- FindVariableFeatures(seurat.shi, selection.method = "vst", nfeatures = 2000)
seurat.shi <- RunFastMNN(object.list = SplitObject(seurat.shi, split.by = "pcw"))
seurat.shi <- RunUMAP(seurat.shi, reduction = "mnn", dims = 1:10)
seurat.shi <- FindNeighbors(seurat.shi, reduction = "mnn", dims = 1:10)
seurat.shi <- FindClusters(seurat.shi, resolution = 0.5)

DefaultAssay(seurat_obj_GE_611) <- "RNA"
DefaultAssay(seurat.shi) <- "mnn.reconstructed"

transfer.anchors <- FindTransferAnchors(
  reference = seurat.shi,
  query = seurat_obj_GE_611,
  reduction = 'cca'
)

predicted.labels.shi <- TransferData(
  anchorset = transfer.anchors,
  refdata = seurat.shi$ClusterID,
  weight.reduction = seurat_obj_GE_611[['lsi']],
  dims = 2:30
)

seurat_obj_GE_611 <- AddMetaData(object = seurat_obj_GE_611, metadata = predicted.labels.shi)

plot1 <- DimPlot(
  object = seurat.shi,
  group.by = 'ClusterID',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = seurat_obj_GE_611,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')

plot1 + plot2





DefaultAssay(seurat_obj_GE_611) <- 'RNA'
DefaultAssay(seurat_obj_GE_611) <- 'peaks'
dir.create(paste0(SIGNAC_DIR, 'peaks'), recursive = T)

# Peaks by cluster
peaks_by_clust <- CallPeaks(
  object = seurat_obj_GE_611,
  group.by = "seurat_clusters",
  macs2.path = '/Users/darren/opt/miniconda3/bin/macs3',
  outdir = paste0(SIGNAC_DIR, 'peaks')
)

peaks_all <- CallPeaks(
  object = seurat_obj_GE_611,
  macs2.path = '/Users/darren/opt/miniconda3/bin/macs3',
  outdir = paste0(SIGNAC_DIR, 'peaks_all')
)




UMAP_list <- list()

for (SAMPLE in SAMPLES) {
  
  SEURAT_OBJ <- get(paste0('seurat_obj_', SAMPLE))

  UMAP_PLOT <- DimPlot(object = SEURAT_OBJ, label = TRUE) + 
    NoLegend() + 
    ggtitle(SAMPLE)
  
  UMAP_list[[SAMPLE]] <- UMAP_PLOT
  QC_list[[SAMPLE]] <- get(paste0('qc_plot_', SAMPLE)) +
    ggtitle(SAMPLE)
  
}

## Merge datasets
for (SAMPLE in SAMPLES) {
  
  SAMPLE_DIR <- paste0('~/Desktop/temp/snATACseq/signac/', SAMPLE, '/')
  
  # read in peak sets
  PEAKS <- read.table(
    file = paste0(SAMPLE_DIR, "peaks.bed"),
    col.names = c("chr", "start", "end"))
  
  # convert to genomic ranges
  RANGES <- makeGRangesFromDataFrame(PEAKS)

  assign(paste0('gr_', SAMPLE), RANGES, .GlobalEnv)
  
}



# Create a unified set of peaks to quantify in each dataset
combined.cameron <- reduce(x = c(gr_14510_PFC_ATAC, gr_14611_PFC_ATAC, gr_14993_PFC_ATAC_AGGR))
combined.ziffra <- reduce(x = c(gr_Cortex_GW18_Frozen, gr_Cortex_GW21_Fresh, gr_GW17_Cortex,
                                gr_PFC_GW20_SC, gr_PFC_Twin34_Frozen, gr_Twin31_PFC))
combined.all <- reduce(x = c(gr_14510_PFC_ATAC, gr_14611_PFC_ATAC, gr_14993_PFC_ATAC_AGGR,
                             gr_Cortex_GW18_Frozen, gr_Cortex_GW21_Fresh, gr_GW17_Cortex,
                             gr_PFC_GW20_SC, gr_PFC_Twin34_Frozen, gr_Twin31_PFC))

# Filter out bad peaks based on length
for (NAME in c('cameron', 'ziffra', 'all')) {
  
  PEAKS <- get(paste0('combined.', NAME))
  
  peakwidths <- width(PEAKS)
  PEAKS <- PEAKS[peakwidths  < 10000 & peakwidths > 20]
  PEAKS
  
  assign(paste0('combined.', NAME), PEAKS, .GlobalEnv)

}

for (SAMPLE in SAMPLES) {
  
  SAMPLE_DIR <- paste0('~/Desktop/temp/snATACseq/signac/', SAMPLE, '/')

  # load metadata
  md.500 <- read.table(
    file = paste0(SAMPLE_DIR, "singlecell.csv"),
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  

}

