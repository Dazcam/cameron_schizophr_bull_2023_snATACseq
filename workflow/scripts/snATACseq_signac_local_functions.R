toc_min <- function(tic,toc,msg="") {
  mins <- round((((toc-tic)/60)),2)
  outmsg <- paste0(mins, " minutes elapsed")
}

toc_hr <- function(tic,toc,msg="") {
  mins <- round((((toc-tic)/3600)),2)
  outmsg <- paste0(mins, " hours elapsed")
}


signac_default_processing <- function( 
    
  SAMPLES = NULL,
  SAMPLE_DIR = NULL,
  MARKER_GENES = NULL) {
  
  #' Run basic Signac processing for set of samples 
  #' (loaded in snATACseq_signac_local_env.R)
  #' 
  #' @description Run basic Signac processing for set of samples 
  #' 
  #' @param SAMPLES A vector of samples names
  #' @param SAMPLE_DIR Base cell ranger output directory
  #' @param MARKER_GENES A vector of (HGNC ID) marker genes 
  
  for (SAMPLE in 1:length(SAMPLES)) {
    
    tic()
    
    cat('\nCreating seurat object for:', SAMPLES[SAMPLE])
    cat('\nLoading Cell Ranger data ... \n\n')
    counts <- Read10X_h5(filename = paste0(SAMPLE_DIR, SAMPLES[SAMPLE],
                                           '/outs/filtered_peak_bc_matrix.h5'))
    metadata <- read.csv(
      file = paste0(SAMPLE_DIR, SAMPLES[SAMPLE],
                    '/outs/singlecell.csv'),
      header = TRUE,
      row.names = 1
    )
    
    chrom_assay <- CreateChromatinAssay(
      counts = counts,
      sep = c(":", "-"),
      genome = 'hg38',
      fragments = paste0(SAMPLE_DIR, SAMPLES[SAMPLE],
                         '/outs/fragments.tsv.gz'),
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
    
    cat('\nRunning QCs ... \n\n')
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
    
    cat('\nSubsetting ... \n\n')
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
    
    cat('\nRunning main Signac processes ... \n\n')
    seurat.obj <- RunTFIDF(seurat.obj)
    seurat.obj <- FindTopFeatures(seurat.obj, min.cutoff = 'q0')
    seurat.obj <- RunSVD(seurat.obj)
    seurat.obj <- RunUMAP(object = seurat.obj, reduction = 'lsi', dims = 2:30)
    seurat.obj <- FindNeighbors(object = seurat.obj, reduction = 'lsi', dims = 2:30)
    seurat.obj <- FindClusters(object = seurat.obj, verbose = FALSE, algorithm = 3)
    DimPlot(object = seurat.obj, label = TRUE) + NoLegend()
    
    cat('\nGenerate gene activity matrix ... \n\n')
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
    
    cat('\nCreate gene marker feature plot ... \n\n')
    GENE_UMAP <- FeaturePlot(
      object = seurat.obj,
      features = MARKER_GENES,
      pt.size = 0.1,
      max.cutoff = 'q95',
      ncol = 3
    )
    
    
    assign(paste0('qc_plot_', SAMPLE_IDs[SAMPLE]), QC_PLOT, .GlobalEnv)
    assign(paste0('gene_UMAP_', SAMPLE_IDs[SAMPLE]), GENE_UMAP, .GlobalEnv)
    assign(paste0('seurat_obj_', SAMPLE_IDs[SAMPLE]), seurat.obj, .GlobalEnv)
    assign(paste0('depth_cor_', SAMPLE_IDs[SAMPLE]), DepthCor(seurat.obj), .GlobalEnv)
    
    toc(func.toc=toc_min)
    
  }
  
} 

signac_snRNAseq_integration <- function(
    
  SEURAT_OBJ = NULL,
  PUBLIC_DATA = NULL) {
  
  #' Run Signac snRNAseq integration on GE data
  #' 
  #' @description Run Signac snRNAseq integration on GE data 
  #' 
  #' @param SEURAT_OBJ A Seurat object for chromatin data 
  #' @param PUBLIC_DATA A GE public dataset to use for integration. Choices
  #' Shi or Cameron.
  
  SHI_DIR <- '~/Desktop/fetal_brain_snRNAseq_GE_270922/workflow/scripts/'
  CAMERON_DIR <- '~/Desktop/fetal_brain_snRNAseq_110122/resources/R_objects/'
  DefaultAssay(SEURAT_OBJ) <- "RNA"
  
  if (PUBLIC_DATA == 'Shi') {
    
    cat('\nGenerating Seurat object for Shi data ...\n\n')
    source(paste0(SHI_DIR, 'snRNAseq_GE_prep_shi_data_for_Seurat.R'))
    seurat.pub <- CreateSeuratObject(counts = shi_data, meta.data = shi_meta)
    seurat.pub <- NormalizeData(seurat.pub)
    seurat.pub <- FindVariableFeatures(seurat.pub, selection.method = "vst", nfeatures = 2000)
    seurat.pub <- RunFastMNN(object.list = SplitObject(seurat.pub, split.by = "pcw"))
    seurat.pub <- RunUMAP(seurat.pub, reduction = "mnn", dims = 1:10)
    seurat.pub <- FindNeighbors(seurat.pub, reduction = "mnn", dims = 1:10)
    seurat.pub <- FindClusters(seurat.pub, resolution = 0.5)
    DefaultAssay(seurat.pub) <- "mnn.reconstructed"
  
} else {
  
  cat('\nLoading Seurat object for Cameron data ...\n\n')
  seurat.pub <- readRDS(paste0(CAMERON_DIR, 'seurat.wge.final.rds'))
  seurat.pub$ClusterID <- seurat.pub$cellIDs
  
  }
  
  cat('\nRunning CCA ...\n\n')  
  transfer.anchors <- FindTransferAnchors(
    reference = seurat.pub,
    query = SEURAT_OBJ,
    reduction = 'cca'
  )
  
  cat('\nTranferring labels  ...\n\n') 
  predicted.labels <- TransferData(
    anchorset = transfer.anchors,
    refdata = seurat.pub$ClusterID,
    weight.reduction = SEURAT_OBJ[['lsi']],
    dims = 2:30
  )
  
  SEURAT_OBJ <- AddMetaData(object = SEURAT_OBJ, metadata = predicted.labels)
  
  if (PUBLIC_DATA == 'Shi') {
    
    SEURAT_OBJ$predicted.labels.shi <<- SEURAT_OBJ$predicted.labels
    
  } else {
    
    SEURAT_OBJ$predicted.labels.cam <<- SEURAT_OBJ$predicted.cam
    
  }
  
  plot1 <- DimPlot(
    object = seurat.pub,
    group.by = 'ClusterID',
    label = TRUE,
    repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')
  
  plot2 <- DimPlot(
    object = SEURAT_OBJ,
    group.by = 'predicted.id',
    label = TRUE,
    repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')
  
  PLOT <- plot1 + plot2
  
  if (PUBLIC_DATA == 'Shi') {
    
    assign('rna_int_plot_shi', PLOT, .GlobalEnv)
    cat('\nReturned rna_int_plot_shi.')
    
  } else {
    
    assign('rna_int_plot_cam', PLOT, .GlobalEnv)
    cat('\nReturned rna_int_plot_cam.')
    
  }

}
