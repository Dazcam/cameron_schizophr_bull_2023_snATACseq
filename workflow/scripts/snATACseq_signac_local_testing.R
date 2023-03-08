#--------------------------------------------------------------------------------------
#
#    ArchR - Signac local testing - GE - draft
#
#--------------------------------------------------------------------------------------

## Load packages  ---------------------------------------------------------------------
SCRIPTS_DIR <- '~/Desktop/fetal_brain_snATACseq_V3_010323/workflow/scripts/'
source(paste0(SCRIPTS_DIR, 'snATACseq_signac_local_env.R'))
source(paste0(SCRIPTS_DIR, 'snATACseq_signac_local_functions.R'))

## Run Signac on default settings -----------------------------------------------------
for (SAMPLE in SAMPLES) {
  
  signac_default_processing(SAMPLES, 
                            SAMPLE_DIR,
                            GE_MARKER_GENES)
  
}

## Get stats
for (SAMPLE in SAMPLE_IDs) {
  
  SEURAT_PRE <- get(paste0('seurat_obj_', SAMPLE))
  SEURAT_POST <- get(paste0('seurat_sub_', SAMPLE))
  
  DefaultAssay(SEURAT_POST) <- "peaks"
  counts_df <- data.frame("Sample" = SAMPLE,
                          "Nuclei_pre_filter" = ncol(SEURAT_PRE),
                          "Nuclei_post_filter" = ncol(SEURAT_POST),
                          "Peaks_pre_filter" = nrow(SEURAT_PRE),
                          "Peaks_post_filter" = nrow(SEURAT_POST),
                          "Median_TSS_Enrich_pre" = median(SEURAT_PRE$TSS.enrichment),
                          "Median_TSS_Enrich_post" = median(SEURAT_POST$TSS.enrichment))
  
  if (exists('counts_df_all')) {
    
    counts_df_all <<- rbind(counts_df_all, counts_df)
    
  } else {
    
    counts_df_all <<- counts_df
    
  }
  
}

# UMAPs
cluster_umaps <- signac_UMAP_per_sample()


## Run Integration with RNA object on individual samples  -----------------------------
for (SAMPLE in SAMPLE_IDs) {

signac_snRNAseq_integration(get(paste0('seurat_sub_', SAMPLE)), 'Shi', SAMPLE)
signac_snRNAseq_integration(get(paste0('seurat_sub_', SAMPLE)), 'Cam', SAMPLE)
  
}


## Call peaks  ------------------------------------------------------------------------
for (SAMPLE in SAMPLE_IDs) {
  
  SEURAT_INT <- get(paste0('seurat_obj_int_', SAMPLE))
  SEURAT_INT$ClusterID <- paste0("C", SEURAT_INT$seurat_clusters)

  
  DefaultAssay(SEURAT_INT) <- 'peaks'
  dir.create(paste0(SIGNAC_DIR, SAMPLE, '/peaks_by_cluster'), recursive = T)
  dir.create(paste0(SIGNAC_DIR, SAMPLE, '/peaks_union'), recursive = T)
  
  # Peaks by cluster
  peaks_by_clust <- CallPeaks(
    object = SEURAT_INT,
    group.by = "ClusterID",
    macs2.path = '/Users/darren/opt/miniconda3/bin/macs3',
    outdir = paste0(SIGNAC_DIR, SAMPLE, '/peaks_by_cluster')
  )

  # Peaks union
  peaks_union <- CallPeaks(
    object = SEURAT_INT,
    macs2.path = '/Users/darren/opt/miniconda3/bin/macs3',
    outdir = paste0(SIGNAC_DIR, SAMPLE, '/peaks_union')
  )

  assign(paste0('peaks_by_clust_', SAMPLE), peaks_by_clust, .GlobalEnv)
  assign(paste0('peaks_union_', SAMPLE), peaks_union, .GlobalEnv)

}


## Merge datasets  --------------------------------------------------------------------
# Only ran this using defaults to get it working
# Vingette uses cellranger peaks file. I'm plugging in the Granges generated above.
# for (SAMPLE in 1:length(SAMPLES)) {
#   
#   # Read in peak sets
#   PEAKS <- read.table(
#     file = paste0(SAMPLE_DIR, SAMPLES[SAMPLE], "/outs/peaks.bed"),
#     col.names = c("chr", "start", "end"))
#   
#   # convert to genomic ranges
#   RANGES <- makeGRangesFromDataFrame(PEAKS)
# 
#   assign(paste0('gr_', SAMPLE_IDs[SAMPLE]), RANGES, .GlobalEnv)
#   
# }

# Create a unified set of peaks to quantify in each dataset
# combined.peaks <- reduce(x = c(gr_GE_510, gr_GE_611, gr_GE_993))
combined.peaks <- reduce(x = c(peaks_by_clust_GE_510, 
                               peaks_by_clust_GE_611, 
                               peaks_by_clust_GE_993))

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

for (SAMPLE in 1:length(SAMPLES)) {
  
  # load metadata
  MD_OBJ <- read.table(
    file = paste0(SAMPLE_DIR, SAMPLES[SAMPLE], "/outs/singlecell.csv"),
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ] # remove the first row
  
  # Filter
  MD_FILT <- MD_OBJ[MD_OBJ$passed_filters > 500, ]
  
  assign(paste0('metadata_', SAMPLE_IDs[SAMPLE]), MD_OBJ, .GlobalEnv)
  assign(paste0('metadata_filt_', SAMPLE_IDs[SAMPLE]), MD_FILT, .GlobalEnv)

}

for (SAMPLE in 1:length(SAMPLES)) {
  
  # create fragment objects
  FRAGS <- CreateFragmentObject(
    path = paste0(SAMPLE_DIR, SAMPLES[SAMPLE], "/outs/fragments.tsv.gz"),
    cells = rownames(get(paste0('metadata_filt_', SAMPLE_IDs[SAMPLE]))))
  
  FRAGS_COUNTS <- FeatureMatrix(
    fragments = FRAGS,
    features = combined.peaks,
    cells = rownames(get(paste0('metadata_filt_', SAMPLE_IDs[SAMPLE]))))

  
  assign(paste0('frags_', SAMPLE_IDs[SAMPLE]), FRAGS, .GlobalEnv)
  assign(paste0('frags_counts', SAMPLE_IDs[SAMPLE]), FRAGS_COUNTS, .GlobalEnv)
  
}

for (SAMPLE in 1:length(SAMPLES)) {
  
  ASSAY <- CreateChromatinAssay(get(paste0('frags_counts', SAMPLE_IDs[SAMPLE])),
                                fragments = get(paste0('frags_', SAMPLE_IDs[SAMPLE])))
  SEURAT_OBJ <- CreateSeuratObject(ASSAY, assay = "ATAC", 
                                meta.data = get(paste0('metadata_filt_', SAMPLE_IDs[SAMPLE])))
  SEURAT_OBJ$dataset <- SAMPLE_IDs[SAMPLE]
  
  assign(paste0('seurat_pre_mrg_', SAMPLE_IDs[SAMPLE]), SEURAT_OBJ, .GlobalEnv)
  
}


# Merge all datasets, adding a cell ID to make sure cell names are unique
seurat.combined <- merge(
  x = seurat_pre_mrg_GE_510,
  y = list(seurat_pre_mrg_GE_611, seurat_pre_mrg_GE_993),
  add.cell.ids = c("GE_510", "GE_611", "GE_993")
)
seurat.combined[["ATAC"]]

seurat.combined <- RunTFIDF(seurat.combined)
seurat.combined <- FindTopFeatures(seurat.combined, min.cutoff = 20)
seurat.combined <- RunSVD(seurat.combined)
seurat.combined <- RunUMAP(seurat.combined, dims = 2:50, reduction = 'lsi')

DimPlot(seurat.combined, group.by = 'dataset', pt.size = 0.1)
Annotation(seurat.combined) <- annotations
cat('\nGenerate gene activity matrix ... \n\n')
gene.activities <- GeneActivity(seurat.combined)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
seurat.combined[['RNA']] <- CreateAssayObject(counts = gene.activities)
seurat.combined <- NormalizeData(
  object = seurat.combined,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(seurat.combined$nCount_RNA)
)

DefaultAssay(seurat.combined) <- 'RNA'

# add the gene information to the object
signac_snRNAseq_integration(seurat.combined, 'Shi', 'combined')
signac_snRNAseq_integration(seurat.combined, 'Cam', 'combined')

R_DIR <- '~/Desktop/fetal_brain_snATACseq_V3_010323/resources/R_obj/'
dir.create(paste0(R_DIR), recursive = T)
save.image(file=paste0(R_DIR, 'snATACseq_signac_local_testing.RData'))
load(file = paste0(R_DIR, 'snATACseq_signac_local_testing.RData'))
#To_do:
#Fix peak_names - done still does not save peak files to dir 
# (no need can generate this myself anyway)

# Number of cells is far higher after merging than we have after individual QC
# Not adding the filtered and QCed objects properly to the merging process
# Move merge to a function
# Need to work out how to merge QC'd objects
# Also no doublet removal in vingette - authors recommend Amulet

