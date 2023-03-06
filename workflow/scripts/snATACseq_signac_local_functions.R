default_signac_processing <- function( 
  
  SAMPLES = NULL,
  SAMPLE_DIR = NULL) {
    
    #' Run basic Signac processing for set of samples  
    #' 
    #' @description Run basic signac processing for set of samples 
    #' 
    #' @param SAMPLES A vector of samples names
    #' @param SAMPLE_DIR Dase directory where cell ranger output for each 
    #' sample is stored
  
  for (SAMPLE in SAMPLES) {

  cat('\nCreating seurat object for:', SAMPLE, '\n\n')
  SAMPLE_DIR <- paste0(SAMPLE_DIR, SAMPLE, '/')
  
  # Load data
  counts <- Read10X_h5(filename = paste0(SAMPLE_DIR, 'filtered_peak_bc_matrix.h5'))
  metadata <- read.csv(
    file = paste0(SAMPLE_DIR, 'singlecell.csv'),
    header = TRUE,
    row.names = 1
  )

  chrom_assay <- CreateChromatinAssay(
    counts = counts,
    sep = c(":", "-"),
    genome = 'hg38',
    fragments = paste0(SAMPLE_DIR, 'fragments.tsv.gz'),
    min.cells = 10,
    min.features = 200
  )
  
  seurat.obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata
  )
  
  seurat.obj[['peaks']]
  
  
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  
  # change to UCSC style since the data was mapped to hg19
  seqlevelsStyle(annotations) <- 'UCSC'
  
  # add the gene information to the object
  Annotation(seurat.obj) <- annotations
  
  # QCs
  # compute nucleosome signal score per cell
  seurat.obj <- NucleosomeSignal(object = seurat.obj)
  
  # compute TSS enrichment score per cell
  seurat.obj <- TSSEnrichment(object = seurat.obj, fast = FALSE)
  
  # add blacklist ratio and fraction of reads in peaks
  seurat.obj$pct_reads_in_peaks <- seurat.obj$peak_region_fragments / seurat.obj$passed_filters * 100
  seurat.obj$blacklist_ratio <- seurat.obj$blacklist_region_fragments / seurat.obj$peak_region_fragments
  
  
  seurat.obj$high.tss <- ifelse(seurat.obj$TSS.enrichment > 2, 'High', 'Low')
  TSSPlot(seurat.obj, group.by = 'high.tss') + NoLegend()
  
  seurat.obj$nucleosome_group <- ifelse(seurat.obj$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  FragmentHistogram(object = seurat.obj, group.by = 'nucleosome_group')
  
  QC_PLOT <- VlnPlot(
    object = seurat.obj,
    features = c('pct_reads_in_peaks', 'peak_region_fragments',
                 'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
    pt.size = 0.1,
    ncol = 5
  )
  
  seurat.obj <- subset(
    x = seurat.obj,
    subset = peak_region_fragments > 3000 &
      peak_region_fragments < 20000 &
      pct_reads_in_peaks > 15 &
      blacklist_ratio < 0.05 &
      nucleosome_signal < 4 &
      TSS.enrichment > 2
  )
  seurat.obj
  
  seurat.obj <- RunTFIDF(seurat.obj)
  seurat.obj <- FindTopFeatures(seurat.obj, min.cutoff = 'q0')
  seurat.obj <- RunSVD(seurat.obj)
  
  DepthCor(seurat.obj)
  
  seurat.obj <- RunUMAP(object = seurat.obj, reduction = 'lsi', dims = 2:30)
  seurat.obj <- FindNeighbors(object = seurat.obj, reduction = 'lsi', dims = 2:30)
  seurat.obj <- FindClusters(object = seurat.obj, verbose = FALSE, algorithm = 3)
  DimPlot(object = seurat.obj, label = TRUE) + NoLegend()
  
  gene.activities <- GeneActivity(seurat.obj)
  
  # add the gene activity matrix to the Seurat object as a new assay and normalize it
  seurat.obj[['RNA']] <- CreateAssayObject(counts = gene.activities)
  seurat.obj <- NormalizeData(
    object = seurat.obj,
    assay = 'RNA',
    normalization.method = 'LogNormalize',
    scale.factor = median(seurat.obj$nCount_RNA)
  )
  
  DefaultAssay(seurat.obj) <- 'RNA'
  
  GENE_UMAP <- FeaturePlot(
    object = seurat.obj,
    features = FC_MARKER_GENES,
    pt.size = 0.1,
    max.cutoff = 'q95',
    ncol = 3
  )
  
  
  assign(paste0('qc_plot_', SAMPLE), QC_PLOT, .GlobalEnv)
  assign(paste0('gene_UMAP_', SAMPLE), GENE_UMAP, .GlobalEnv)
  assign(paste0('seurat_obj_', SAMPLE), seurat.obj)

  }
  
} 

